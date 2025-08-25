#!/usr/bin/env python3
import os
import itertools
from pathlib import Path

# Import the selection utilities from the local script
import sys
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import select_transcripts as st  # type: ignore

BASE = SCRIPT_DIR
SINGLE_SAMPLE_DIR = BASE / "single_sample"
OUT_BASE = BASE / "intersection_all_comb"
ANNOTATION = Path("/Users/tianyuan/Desktop/github_dev/tusco-paper/figs/data/reference/mm39.ncbiRefSeq_SIRV.gtf")

GROUPS = {
    "K": sorted([p for p in SINGLE_SAMPLE_DIR.glob("K*.isoforms.gtf")]),
    "B": sorted([p for p in SINGLE_SAMPLE_DIR.glob("B*.isoforms.gtf")]),
}

assert ANNOTATION.exists(), f"Annotation not found: {ANNOTATION}"

# Parse annotation once
annotation_transcripts = st.parse_gtf_file(str(ANNOTATION))
overlap_index = st.build_span_overlap_index(annotation_transcripts)
annotation_index = st.build_annotation_index(annotation_transcripts)

# Pre-parse all sample GTFs once
parsed_by_file = {}
for group, files in GROUPS.items():
    for fpath in files:
        parsed_by_file[str(fpath)] = st.parse_gtf_file(str(fpath))

OUT_BASE.mkdir(parents=True, exist_ok=True)

num_written = 0
for group, files in GROUPS.items():
    if len(files) < 2:
        continue
    for r in range(2, len(files) + 1):
        for combo in itertools.combinations(files, r):
            combo_names = [p.stem for p in combo]  # e.g., K31.isoforms
            combo_names = [name.replace(".isoforms", "") if name.endswith(".isoforms") else name for name in combo_names]
            combo_name = "_".join(combo_names)
            out_dir = OUT_BASE / combo_name
            out_dir.mkdir(parents=True, exist_ok=True)
            out_gtf = out_dir / f"{combo_name}_common.gtf"
            if out_gtf.exists() and out_gtf.stat().st_size > 0:
                continue

            # Build transcripts_by_file for this combo with pre-filtering against annotation overlap
            transcripts_by_file = {}
            for fpath in combo:
                tx_dict = parsed_by_file[str(fpath)]
                # Pre-filter: keep only transcripts overlapping annotation on same strand
                filtered = {tx_id: tx for tx_id, tx in tx_dict.items() if st.transcript_overlaps_annotation(tx, overlap_index)}
                transcripts_by_file[str(fpath)] = filtered

            # Skip if empty after filtering
            total_remaining = sum(len(d) for d in transcripts_by_file.values())
            if total_remaining == 0:
                # write empty file to indicate processed
                out_gtf.write_text("")
                continue

            groups = st.group_transcripts(transcripts_by_file)
            selected = st.select_transcripts_by_mode(groups, mode="intersect", annotation_index=annotation_index)
            if not selected:
                out_gtf.write_text("")
                continue

            st.write_gtf(selected, str(out_gtf))
            num_written += 1

print(f"Wrote {num_written} combination GTF(s) under {OUT_BASE}")
