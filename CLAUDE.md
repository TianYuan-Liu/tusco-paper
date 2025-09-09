# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains code and data for the TUSCO paper (Transcriptome Universal Single-isoform COntrol). TUSCO is an endogenous benchmarking framework for long-read RNA-seq that uses a curated set of widely expressed single-isoform genes (human and mouse) to evaluate transcript reconstruction precision and sensitivity without external spike-ins.

## Repository Structure

### Figure Organization
- All figures are organized under `figs/fig-<n>/` (main figures) or `figs/fig-s<n>/` (supplementary figures)
- Each figure directory contains:
  - `code/`: R scripts and analysis code to generate the figure
  - `plot/`: Final vector outputs (PDF preferred)
  - `tsv/`: Data files and intermediate results

### Data Organization
- Core data is stored under `figs/data/` with subdirectories:
  - `tusco/`: TUSCO gene sets and GTF files (`tusco_human.gtf`, `tusco_mouse.gtf`, TSV files)
  - `lrgasp/`: LRGASP benchmark data
  - `spike-ins/`: SIRV spike-in data
  - `reference/`: Reference annotations
  - `RIN/`: RNA integrity number analysis data
  - `nih/`: NIH dataset analysis

## Key Development Commands

### Running R Scripts
Most figure generation scripts are executable R scripts with shebangs:
```bash
# Execute individual figure scripts
./figs/fig-3/code/figure3a-human.R
./figs/fig-3/code/figure3c.R

# Or using Rscript
Rscript figs/fig-3/code/figure3a-human.R
```

### TUSCO Analysis Pipeline
Key shell scripts for data processing:
- `figs/fig-5/analysis/tusco_replication_sqanti_pipeline.sh`
- `figs/fig-5/analysis/rerun_all_intersection_sqanti.sh`

## Code Architecture and Conventions

### R Script Structure
All R scripts follow consistent patterns:
1. **Library loading**: Use `suppressPackageStartupMessages()` for clean output
2. **Path resolution**: Use `resolve_path()` helper function for portable file paths
3. **Environment setup**: Set OpenMP variables to avoid SHM issues
4. **Parameterized output**: Accept `out_dir`, `width`, `height` parameters

### TUSCO Classification System
- **TP (True Positive)**: SQANTI3 `reference_match` reproducing TUSCO exon-intron chain with TSS/TTS within ±50 bp
- **PTP (Partial True Positive)**: SQANTI3 `full-splice_match` (non-RM) or `incomplete-splice_match`
- **FP (False Positive)**: SQANTI3 `NIC`, `NNC`, `Genic Intron`, `Genic Genomic`, or `Fusion` within TUSCO genes
- **FN (False Negative)**: TUSCO genes lacking TP or PTP calls

### Metrics Calculations
- `Sn = 100 * TP_genes / G` (Sensitivity, G = # TUSCO genes)
- `nrPre = 100 * TP / N` (Non-redundant Precision)
- `rPre = 100 * (TP + PTP) / N` (Redundant Precision)
- `FDR = 100 * (N - TP) / N` (False Discovery Rate)
- `PDR = 100 * detected_genes(TP or PTP) / G` (Partial Detection Rate)

### Plotting Standards
- Use `theme_classic(base_family = "Helvetica", base_size = 7)`
- Export with Cairo PDF device with explicit figure sizes
- Output to `plots/` directory by default
- Use relative paths only - never hardcode absolute paths

### Path Management
Scripts use a `resolve_path()` function to handle both absolute and relative paths:
```r
resolve_path <- function(candidates, is_dir = FALSE) {
  for (p in candidates) {
    if (!is_dir && file.exists(p)) return(p)
    if (is_dir && dir.exists(p)) return(p)
  }
  return(candidates[[1]])
}
```

## Development Guidelines

### File Creation Rules
- **DO**: Use relative paths rooted at repository base
- **DO**: Accept `out_dir`, `width`, `height` as parameters
- **DO**: Output vector formats (PDF preferred) to `plots/` directory
- **DON'T**: Hardcode absolute filesystem paths
- **DON'T**: Write intermediate files to `plots/` (use `tmp/` if needed)

### Code Style
- Keep R scripts modular and parameterized
- Use consistent library loading patterns
- Map SQANTI3 classifications to TUSCO labels before computing metrics
- For replicate studies, compute intersection sets by identical junction chains before labeling

## Common Tasks

### Generating Figures
1. Navigate to appropriate figure directory: `cd figs/fig-<n>/code/`
2. Execute R script: `./script_name.R` or `Rscript script_name.R`
3. Check output in `../plot/` directory

### Adding New Figures
1. Create directory structure: `figs/fig-<n>/code/` and `figs/fig-<n>/plot/`
2. Follow existing R script patterns for portability
3. Use relative paths for data access under `figs/data/`
4. Export PDFs with explicit dimensions to `plot/` directory

### Working with TUSCO Data
- Human TUSCO genes: `figs/data/tusco/tusco_human.tsv`
- Mouse TUSCO genes: `figs/data/tusco/tusco_mouse.tsv`
- GTF files available for both species with tissue-specific variants