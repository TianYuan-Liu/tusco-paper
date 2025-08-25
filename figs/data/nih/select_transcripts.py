#!/usr/bin/env python3
"""
select_common_transcript.py

This script selects transcript representatives from one or more GTF files by
grouping transcripts with identical structure and then choosing a single
representative per structure using an annotation-aware rule.

Modes
- intersect: keep structures present in ≥2 input files.
- merge: union of all structures across inputs (deduplicated).

Selection Rule (applies when duplicates exist within the same structure key)
1) If an annotation GTF is provided:
   For each candidate transcript and each annotation transcript on the same
   chromosome and strand, compute:
       score_binary = 0
       if |TSS_candidate – TSS_annotation| > 50 → score_binary += 1
       if |TTS_candidate – TTS_annotation| > 50 → score_binary += 1
   For that candidate, take the annotation transcript that minimizes
   score_binary (ties allowed). If there is a tie, compute
   D = |TSS_candidate – TSS_annotation| + |TTS_candidate – TTS_annotation|
   and take the smallest D. Candidates are ranked by (score_binary, D).
2) If no annotation is provided (or no same-chrom/strand annotations exist),
   pick the transcript with the longest length (end – start).

Usage:
    ./select_common_transcript.py -f file1.gtf file2.gtf -o out.gtf \
        --mode intersect|merge [--annotation anno.gtf]

Author: Tianyuan Liu
Date: 2025-02-09
"""

import argparse
import sys
from typing import Dict, List, Optional, Tuple
from bisect import bisect_left, bisect_right


def parse_gtf_attributes(attr_str):
    """Parse a GTF attribute string into a dictionary."""
    attributes = {}
    # Split by semicolon; ignore empty entries.
    for attr in attr_str.strip().split(";"):
        attr = attr.strip()
        if not attr:
            continue
        # Expecting key "value"
        parts = attr.split(" ", 1)
        if len(parts) != 2:
            continue
        key, value = parts
        # Remove any surrounding quotes.
        value = value.strip().strip('"')
        attributes[key] = value
    return attributes


def parse_gtf_file(filename):
    """
    Parse a GTF file and return a dictionary mapping transcript_id to a transcript object.
    The transcript object is a dict with:
      - chrom, strand
      - exons: a list of (start, end) tuples
      - tx_start: minimum exon start
      - tx_end: maximum exon end
      - length: tx_end - tx_start
      - junctions: for transcripts with >=2 exons, a tuple of internal junction pairs:
                   ((exon1_end, exon2_start), (exon2_end, exon3_start), ...)
      - attributes: the attribute dictionary from the first exon (used for transcript_id, gene_id, etc.)
    """
    transcripts = {}
    with open(filename) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature, start_str, end_str, score, strand, frame, attr_str = fields
            # We only care about exon lines.
            if feature.lower() != "exon":
                continue
            try:
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                sys.exit("Error converting start/end to integer in line:\n" + line)
            attrs = parse_gtf_attributes(attr_str)
            if "transcript_id" not in attrs:
                continue
            tx_id = attrs["transcript_id"]
            if tx_id not in transcripts:
                transcripts[tx_id] = {
                    "chrom": chrom,
                    "strand": strand,
                    "exons": [],
                    "attributes": attrs
                }
            transcripts[tx_id]["exons"].append((start, end))
    # Process each transcript
    for tx_id, tx in transcripts.items():
        # Sort exons by start coordinate
        tx["exons"].sort(key=lambda x: x[0])
        tx["tx_start"] = min(start for start, end in tx["exons"])
        tx["tx_end"] = max(end for start, end in tx["exons"])
        tx["length"] = tx["tx_end"] - tx["tx_start"]
        if len(tx["exons"]) >= 2:
            # Compute internal junctions: (exon[i].end, exon[i+1].start)
            junctions = tuple((tx["exons"][i][1], tx["exons"][i + 1][0]) for i in range(len(tx["exons"]) - 1))
            tx["junctions"] = junctions
        else:
            tx["junctions"] = None
    return transcripts


def group_transcripts(transcripts_by_file):
    """
    Given a dict mapping filename -> {transcript_id: transcript_object},
    group transcripts across files by a “common” key.

    For multi-exon transcripts (junctions not None):
         Key = (chrom, strand, junctions)
    For single-exon transcripts:
         Key = (chrom, strand, tx_start, tx_end)

    Returns a dictionary mapping key -> dict with:
       - "files": a set of filenames in which at least one transcript with that key was found
       - "transcripts": a list of transcript objects from all files with that key
    """
    groups = {}
    for fname, tx_dict in transcripts_by_file.items():
        for tx_id, tx in tx_dict.items():
            if tx["junctions"] is not None:
                key = (tx["chrom"], tx["strand"], tx["junctions"])
            else:
                key = (tx["chrom"], tx["strand"], tx["tx_start"], tx["tx_end"])
            if key not in groups:
                groups[key] = {"files": set(), "transcripts": []}
            groups[key]["files"].add(fname)
            groups[key]["transcripts"].append(tx)
    return groups


def compute_tss_tts(transcript: Dict) -> Tuple[int, int]:
    """Return (TSS, TTS) positions given transcript strand and tx bounds.

    For '+' strand: TSS = tx_start, TTS = tx_end
    For '-' strand: TSS = tx_end,   TTS = tx_start
    """
    strand = transcript["strand"]
    tx_start = transcript.get("tx_start", min(s for s, _ in transcript["exons"]))
    tx_end = transcript.get("tx_end", max(e for _, e in transcript["exons"]))
    if strand == "+":
        return tx_start, tx_end
    else:
        return tx_end, tx_start


def build_annotation_index(annotation_transcripts_by_txid: Optional[Dict[str, Dict]]):
    """Index annotation transcripts by (chrom, strand) for fast TSS/TTS lookup.

    Returns dict: (chrom, strand) -> list of (tss, tts)
    """
    if not annotation_transcripts_by_txid:
        return {}
    index: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
    for tx in annotation_transcripts_by_txid.values():
        key = (tx["chrom"], tx["strand"])
        tss, tts = compute_tss_tts(tx)
        index.setdefault(key, []).append((tss, tts))
    return index


def build_span_overlap_index(annotation_transcripts_by_txid: Optional[Dict[str, Dict]]) -> Dict[Tuple[str, str], Dict[str, List[int]]]:
    """Build an index of annotation transcript spans for overlap queries.

    Returns dict: (chrom, strand) -> { "starts": [sorted starts], "ends": [sorted ends] }
    """
    index: Dict[Tuple[str, str], Dict[str, List[int]]] = {}
    if not annotation_transcripts_by_txid:
        return index
    for tx in annotation_transcripts_by_txid.values():
        chrom = tx["chrom"]
        strand = tx["strand"]
        tx_start = tx.get("tx_start", min(s for s, _ in tx["exons"]))
        tx_end = tx.get("tx_end", max(e for _, e in tx["exons"]))
        bucket = index.setdefault((chrom, strand), {"starts": [], "ends": []})
        bucket["starts"].append(tx_start)
        bucket["ends"].append(tx_end)
    for bucket in index.values():
        bucket["starts"].sort()
        bucket["ends"].sort()
    return index


def interval_overlaps_any(query_start: int, query_end: int, starts_sorted: List[int], ends_sorted: List[int]) -> bool:
    """Return True if [query_start, query_end] overlaps any interval in the index."""
    if not starts_sorted:
        return False
    # Count of intervals with start <= query_end
    num_start_leq_end = bisect_right(starts_sorted, query_end)
    # Minus intervals that end before query_start
    num_end_lt_start = bisect_left(ends_sorted, query_start)
    return (num_start_leq_end - num_end_lt_start) > 0


def transcript_overlaps_annotation(transcript: Dict, overlap_index: Dict[Tuple[str, str], Dict[str, List[int]]]) -> bool:
    """Check if transcript overlaps any annotated transcript on same chrom/strand.

    Overlap means any overlap across the transcript span (tx_start..tx_end). This
    also implicitly covers exon-overlap cases.
    """
    key = (transcript["chrom"], transcript["strand"])
    bucket = overlap_index.get(key)
    if not bucket:
        return False
    tx_start = transcript.get("tx_start", min(s for s, _ in transcript["exons"]))
    tx_end = transcript.get("tx_end", max(e for _, e in transcript["exons"]))
    return interval_overlaps_any(tx_start, tx_end, bucket["starts"], bucket["ends"])


def best_annotation_metrics_for_candidate(candidate: Dict, annotation_index: Dict) -> Optional[Tuple[int, int]]:
    """For a candidate transcript, compute (min_score_binary, min_D) vs annotations.

    Restrict comparisons to same (chrom, strand). If no annotations on that key,
    return None.
    """
    key = (candidate["chrom"], candidate["strand"])
    ann_list = annotation_index.get(key)
    if not ann_list:
        return None
    cand_tss, cand_tts = compute_tss_tts(candidate)
    best_score = None
    best_D = None
    THRESHOLD = 50
    for ann_tss, ann_tts in ann_list:
        score = 0
        if abs(cand_tss - ann_tss) > THRESHOLD:
            score += 1
        if abs(cand_tts - ann_tts) > THRESHOLD:
            score += 1
        D = abs(cand_tss - ann_tss) + abs(cand_tts - ann_tts)
        if best_score is None or (score < best_score) or (score == best_score and D < best_D):
            best_score = score
            best_D = D
    return best_score, best_D


def select_representative(transcripts: List[Dict], annotation_index: Dict) -> Dict:
    """Select one representative from a list using the Selection Rule.

    If annotation is available on same chrom/strand, choose by (score_binary, D).
    Otherwise choose the longest.
    """
    # Try annotation-guided selection
    scored: List[Tuple[Tuple[int, int], Dict]] = []
    for tx in transcripts:
        metrics = best_annotation_metrics_for_candidate(tx, annotation_index)
        if metrics is not None:
            scored.append((metrics, tx))
    if scored:
        scored.sort(key=lambda item: (item[0][0], item[0][1], -item[1]["length"]))
        return scored[0][1]
    # Fallback: longest length
    return max(transcripts, key=lambda t: t["length"]) if transcripts else None


def select_transcripts_by_mode(groups: Dict, mode: str, annotation_index: Dict, min_files_for_intersect: int = 2) -> List[Dict]:
    """Select transcripts across groups according to mode and selection rule."""
    selected: List[Dict] = []
    for key, group in groups.items():
        if mode == "intersect":
            if len(group["files"]) < min_files_for_intersect:
                continue
        elif mode == "merge":
            pass  # keep all groups
        else:
            raise ValueError(f"Unsupported mode: {mode}")
        rep = select_representative(group["transcripts"], annotation_index)
        if rep is not None:
            selected.append(rep)
    return selected


def write_gtf(transcripts, output_filename):
    """
    Write selected transcripts to a GTF file with both transcript and exon features.
    For each transcript, emit one 'transcript' line spanning tx_start..tx_end, followed by
    one 'exon' line per exon.
    """
    with open(output_filename, "w") as outfh:
        for tx in transcripts:
            chrom = tx["chrom"]
            source = "selected"
            strand = tx["strand"]
            score = "."
            frame = "."
            # Attributes
            tx_id = tx["attributes"].get("transcript_id", "NA")
            gene_id = tx["attributes"].get("gene_id", "NA")
            attr_str = f'gene_id "{gene_id}"; transcript_id "{tx_id}";'

            # Transcript feature line
            tx_start = tx.get("tx_start", min(start for start, _ in tx["exons"]))
            tx_end = tx.get("tx_end", max(end for _, end in tx["exons"]))
            tx_fields = [
                chrom,
                source,
                "transcript",
                str(tx_start),
                str(tx_end),
                score,
                strand,
                frame,
                attr_str,
            ]
            outfh.write("\t".join(tx_fields) + "\n")

            # Exon feature lines
            for exon in tx["exons"]:
                start, end = exon
                exon_fields = [chrom, source, "exon", str(start), str(end), score, strand, frame, attr_str]
                outfh.write("\t".join(exon_fields) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Select transcript representatives across multiple GTF files. "
            "If --annotation is provided, first remove any transcript that does not "
            "overlap an annotated transcript on the same strand (by span or exons). "
            "Mode 'intersect' keeps structures present in ≥2 files; "
            "mode 'merge' deduplicates the union across files."
        )
    )
    parser.add_argument("-f", "--files", nargs="+", required=True, help="List of GTF files (same genome build)")
    parser.add_argument("-o", "--output", required=True, help="Output GTF file for selected transcripts")
    parser.add_argument("--mode", choices=["intersect", "merge"], default="intersect", help="Selection mode")
    parser.add_argument(
        "--annotation",
        required=False,
        help=(
            "Optional annotation GTF. Used to pre-filter inputs to only transcripts "
            "overlapping annotation (same strand), and to guide duplicate selection "
            "via TSS/TTS proximity."
        ),
    )
    args = parser.parse_args()

    input_files = args.files
    transcripts_by_file: Dict[str, Dict[str, Dict]] = {}
    for fname in input_files:
        transcripts_by_file[fname] = parse_gtf_file(fname)

    # Parse annotation once (if provided), build indices, and pre-filter inputs
    annotation_transcripts: Optional[Dict[str, Dict]] = None
    annotation_index: Dict = {}
    overlap_index: Dict = {}
    if args.annotation:
        annotation_transcripts = parse_gtf_file(args.annotation)
        overlap_index = build_span_overlap_index(annotation_transcripts)
        # Pre-filter: keep only transcripts that overlap annotation (same strand)
        if overlap_index:
            for fname, tx_dict in list(transcripts_by_file.items()):
                filtered = {tx_id: tx for tx_id, tx in tx_dict.items() if transcript_overlaps_annotation(tx, overlap_index)}
                transcripts_by_file[fname] = filtered

    # If nothing remains after filtering, exit early with a clear message
    total_remaining = sum(len(d) for d in transcripts_by_file.values())
    if total_remaining == 0:
        sys.exit("No transcripts remain after filtering for annotation overlap.")

    groups = group_transcripts(transcripts_by_file)

    # Build annotation index for TSS/TTS-guided selection if provided
    if annotation_transcripts:
        annotation_index = build_annotation_index(annotation_transcripts)

    selected_transcripts = select_transcripts_by_mode(groups, args.mode, annotation_index)

    if not selected_transcripts:
        if args.mode == "intersect":
            sys.exit("No intersecting transcript structures found (need ≥2 files per structure).")
        else:
            sys.exit("No transcripts found to merge.")

    write_gtf(selected_transcripts, args.output)
    print(f"Selected {len(selected_transcripts)} transcript(s) in {args.mode} mode. Output written to {args.output}")


if __name__ == "__main__":
    main()