# tusco-paper

Project structure (reorganized for clarity and reproducibility):

- `data/`
  - `raw/`: External inputs and large downloads (reference, lrgasp, expression, nih, spike-ins)
  - `processed/`: Derived, versioned outputs (e.g., `processed/tusco/{hsa,mmu}`)
- `figs/`
  - `figure-0N/` and `supp-fig-0N/`: Each with `code/`, `plots/`, `tables/`
- `src/`: Reusable Python code (`tusco_selector`, `tusco_novel_simulator`)
- `R/`: Shared R helpers (e.g., `R/paths.R`)
- `envs/`: Conda environments (e.g., `envs/tusco_selector.yml`)
- `config/`: Central configuration (e.g., `config/project.yml`)
- `workflows/`: Pipelines and SLURM jobs (optional)

Notes
- All scripts now read from `data/raw` and write reusable results to `data/processed`.
- Run Python modules from repo root with `export PYTHONPATH=src`.
