# Repository Guidelines

## Project Structure & Module Organization
- `figs/`: Figure sources and assets. Each figure lives in a subfolder (for example, `fig-1`, `fig-3`, `fig-s3`). Shared inputs are under `figs/data/` (reference annotations, expression matrices, intermediates) with stable paths used across figures.
- `manuscript/`: Draft materials (`TUSCO.docx`).
- Root: Repository metadata (`README.md`, `.gitattributes`, `.gitignore`) and occasional artifacts (for example, `Rplots.pdf`).

## Build, Test, and Development Commands
- This repo has no monolithic build; run scripts per figure or dataset.
- Python example: `python figs/data/nih/generate_intersections_all_comb.py` (generates intersection/union artifacts for NIH datasets).
- R usage: Place ad‑hoc analysis scripts in the relevant figure folder and run with `Rscript path/to/script.R` when present. Add a brief `README.md` in each new `fig-*` folder with exact commands and expected outputs.

## Coding Style & Naming Conventions
- Python: PEP 8, 4‑space indents, descriptive names, prefer `snake_case` for files and functions. Add docstrings and minimal type hints where helpful.
- R: Tidy, consistent style; comment section headers; prefer explicit arguments.
- Folders: Keep existing pattern `fig-<n>` and `fig-s<n>`; avoid renaming shared data paths in `figs/data/` to maintain reproducibility.

## Testing Guidelines
- Prefer small, deterministic checks over heavy test suites. Include shape/row counts, header presence, and ID uniqueness assertions in scripts.
- If feasible, add a tiny sample input alongside a script (for example, `samples/`) and document expected outputs.
- Save derived assets inside each figure folder under `output/` to keep provenance clear.

## Commit & Pull Request Guidelines
- Commits: Imperative mood, scoped prefixes help triage: `[fig-3]`, `[data]`, `[docs]`, `[infra]` (for LFS/ignore updates).
- PRs: Include what changed, affected figures, reproduction commands, and paths to new outputs (screenshots or PDFs if applicable). Link issues and note data sources.
- Large files: This project uses Git LFS. Ensure LFS is installed (`git lfs install`) and do not commit >10 MB binaries outside tracked patterns.

## Security & Configuration Tips
- No secrets should be committed. External credentials (if any) must be sourced locally.
- Pin environments locally (for example, Python ≥3.10, recent R) and record version notes in per‑figure READMEs to aid reproducibility.

