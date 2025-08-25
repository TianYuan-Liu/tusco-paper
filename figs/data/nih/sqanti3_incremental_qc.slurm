#!/bin/bash

# Re-exec with bash if invoked via 'sh'
if [ -z "${BASH_VERSION:-}" ]; then
	exec /usr/bin/env bash "$0" "$@"
fi

#SBATCH --job-name=sqanti3_intersections_all
#SBATCH --output=sqanti3_intersections_all_%A_%a.out
#SBATCH --error=sqanti3_intersections_all_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=24G
#SBATCH --qos=short
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tianyuan.liu@csic.es
# NOTE: Set the array range to cover all combos under
# $REPO_ROOT/figs/data/nih/intersection_all_comb. For 5 K and 5 B
# samples (all 2..5 combos per group), there are 52 combos total.
# You can override this at submission with: sbatch --array=0-51
#SBATCH --array=0-51

#############################################
# 1) Load environment
#############################################
# Load modules if available (HPC); ignore on machines without 'module'
if command -v module >/dev/null 2>&1; then
	module load anaconda || true
fi

# Initialize conda for non-interactive shells and activate env if present
if command -v conda >/dev/null 2>&1; then
	eval "$(conda shell.bash hook)" 2>/dev/null || true
	conda activate SQANTI3.env || true
fi

#############################################
# 2) Paths and references (repo-local)
#############################################
# Resolve repo root (this script lives in figs/data/nih)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Reference data in repository
ANNOTATION="$REPO_ROOT/figs/data/reference/mm39.ncbiRefSeq_SIRV.gtf"
GENOME="$REPO_ROOT/figs/data/reference/mouse/mm39_SIRV.fa"

# CAGE peaks removed - previously optional

# Directory containing all new combinations (generated beforehand)
INTERSECTION_DIR="$REPO_ROOT/figs/data/nih/intersection_all_comb"

# Directory containing single-sample GTFs
SINGLE_SAMPLE_DIR="$REPO_ROOT/figs/data/nih/single_sample"

#############################################
# 3) Build combo list from INTERSECTION_DIR
#############################################
# Gather all combination directories (portable to macOS; BSD find lacks -printf,
# and macOS bash 3.2 lacks mapfile). Build the array via a while-read loop.
ALL_COMBOS=()
if [ -d "$INTERSECTION_DIR" ]; then
	while IFS= read -r dirpath; do
		ALL_COMBOS+=( "$(basename "$dirpath")" )
	done < <(find "$INTERSECTION_DIR" -mindepth 1 -maxdepth 1 -type d -print 2>/dev/null | sort)
else
    echo "WARN: INTERSECTION_DIR not found: $INTERSECTION_DIR (skipping multi-sample)"
fi

#############################################
# 3a) Resolve SQANTI3 executable once
#############################################
SQANTI3_BIN="/Users/tianyuan/Desktop/GitHub/SQANTI3/sqanti3_qc.py"
if [ ! -f "$SQANTI3_BIN" ]; then
	if command -v sqanti3_qc.py >/dev/null 2>&1; then
		SQANTI3_BIN="sqanti3_qc.py"
	else
		echo "ERROR: SQANTI3 not found at /Users/tianyuan/Desktop/GitHub/SQANTI3 or in PATH"
		exit 1
	fi
fi

#############################################
# 3b) Local-run mode: process ALL combos (always, no SLURM)
#############################################
if true; then
    #############################################
    # 3b) Single-sample processing (idempotent)
    #############################################
	if [ -d "$SINGLE_SAMPLE_DIR" ]; then
		SINGLE_SAMPLE_GTFS=()
		while IFS= read -r gtf; do
			SINGLE_SAMPLE_GTFS+=( "$gtf" )
		done < <(find "$SINGLE_SAMPLE_DIR" -mindepth 1 -maxdepth 1 -type f \( -name "*.gtf" -o -name "*.GTF" \) -print 2>/dev/null | sort)

		if [ ${#SINGLE_SAMPLE_GTFS[@]} -eq 0 ]; then
			echo "[Single-sample] No .gtf files found under $SINGLE_SAMPLE_DIR"
		else
			echo "[Single-sample] Processing ${#SINGLE_SAMPLE_GTFS[@]} samples under $SINGLE_SAMPLE_DIR"
			for SAMPLE_GTF in "${SINGLE_SAMPLE_GTFS[@]}"; do
				if [ ! -s "$SAMPLE_GTF" ]; then
					echo "[Single-sample] WARN: Skipping empty or missing GTF: $SAMPLE_GTF"
					continue
				fi

				SAMPLE_BASE="$(basename "$SAMPLE_GTF")"
				SAMPLE_NAME="${SAMPLE_BASE%.gtf}"
				SAMPLE_NAME="${SAMPLE_NAME%.GTF}"
				SAMPLE_OUTDIR="$SINGLE_SAMPLE_DIR/$SAMPLE_NAME"
				mkdir -p "$SAMPLE_OUTDIR"

				SAMPLE_CLASS_TXT="$SAMPLE_OUTDIR/${SAMPLE_NAME}_classification.txt"
				if [ -s "$SAMPLE_CLASS_TXT" ]; then
					echo "[Single-sample] Skip $SAMPLE_NAME (already finished: $SAMPLE_CLASS_TXT)"
					continue
				fi

				SQANTI_ARGS=(
				  --isoforms "$SAMPLE_GTF"
				  --refGTF "$ANNOTATION"
				  --refFasta "$GENOME"
				  -o "$SAMPLE_NAME"
				  -d "$SAMPLE_OUTDIR"
				  --skipORF
				  --report skip
				)

				"$SQANTI3_BIN" "${SQANTI_ARGS[@]}"
				if [ $? -ne 0 ]; then
					echo "[Single-sample] ERROR: SQANTI3 QC failed on ${SAMPLE_NAME}"
					exit 1
				fi

				echo "[Single-sample] Completed for ${SAMPLE_NAME}"
			done
		fi
	else
		echo "[Single-sample] WARN: SINGLE_SAMPLE_DIR not found: $SINGLE_SAMPLE_DIR (skipping single-sample)"
	fi

	echo "[Local mode] Processing all ${#ALL_COMBOS[@]} combos under $INTERSECTION_DIR"
	for combo_name in "${ALL_COMBOS[@]}"; do
		COMBO_DIR="$INTERSECTION_DIR/$combo_name"
		mkdir -p "$COMBO_DIR"

		CAND_GTF="$COMBO_DIR/${combo_name}_common.gtf"
		if [ ! -s "$CAND_GTF" ]; then
			echo "ERROR: Candidate GTF not found or empty: $CAND_GTF"
			exit 1
		fi

		COMBO_CLASS_TXT="$COMBO_DIR/${combo_name}_classification.txt"
		if [ -s "$COMBO_CLASS_TXT" ]; then
			echo "[SQANTI3 QC] Skip ${combo_name} (already finished: $COMBO_CLASS_TXT)"
			continue
		fi

		SQANTI_ARGS=(
		  --isoforms "$CAND_GTF"
		  --refGTF "$ANNOTATION"
		  --refFasta "$GENOME"
		  -o "$combo_name"
		  -d "$COMBO_DIR"
		  --skipORF
		  --report skip
		)

		"$SQANTI3_BIN" "${SQANTI_ARGS[@]}"
		if [ $? -ne 0 ]; then
			echo "[SQANTI3 QC] ERROR: SQANTI3 QC failed on ${combo_name}"
			exit 1
		fi

		echo "[SQANTI3 QC] Completed for ${combo_name}"
	done

	echo "[Local mode] Completed all combos"
	exit 0
fi