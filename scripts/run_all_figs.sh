#!/usr/bin/env bash
set -euo pipefail

# Minimal runner to execute all figure R scripts and report success/failure.
# Uses repo-rooted data layout (data/raw, data/processed). Does not install R packages.
# For heavy pipelines, set LOCAL_ONLY=1 where supported to avoid unnecessary work.

ROOT_DIR=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT_DIR"

# Mitigate OpenMP SHM issues in restricted environments
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=FALSE
export OMP_WAIT_POLICY=PASSIVE
export KMP_INIT_AT_FORK=0
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

SCRIPTS=(
  figs/figure-01/code/figure1d.R
  figs/figure-01/code/figure1c_figure-s1.R
  figs/figure-03/code/figure3a-human.R
  figs/figure-03/code/figure3a-mouse.R
  figs/figure-03/code/figure3b-human.R
  figs/figure-03/code/figure3b-mouse.R
  figs/figure-03/code/figure3c.R
  figs/figure-03/code/figure3d-human.R
  figs/figure-03/code/figure3d-mouse.R
  figs/figure-03/code/table-s1.R
  figs/figure-04/code/fig-4b.R
  figs/figure-04/code/fig-s5.R
  figs/figure-05/code/fig-5b-5c-s6.R
  figs/supp-fig-02/code/fig-s2.R
  figs/supp-fig-03/code/FigureS3_correlation.R
  figs/supp-fig-04/code/fig-s4.R
)

STATUS_KEYS=()
STATUS_VALS=()

for s in "${SCRIPTS[@]}"; do
  if [[ ! -f "$s" ]]; then
    STATUS_KEYS+=("$s"); STATUS_VALS+=("SKIP (missing script)")
    echo "[skip] $s (not found)"
    continue
  fi
  echo "[run] $s"
  # Special env flags for specific scripts
  if [[ "$s" == *"fig-5b-5c-s6.R"* ]]; then
    if LOCAL_ONLY=1 Rscript "$s" >/dev/null 2>"${s%.R}.err"; then
      STATUS_KEYS+=("$s"); STATUS_VALS+=("OK")
    else
      STATUS_KEYS+=("$s"); STATUS_VALS+=("FAIL")
    fi
  else
    if Rscript "$s" >/dev/null 2>"${s%.R}.err"; then
      STATUS_KEYS+=("$s"); STATUS_VALS+=("OK")
    else
      STATUS_KEYS+=("$s"); STATUS_VALS+=("FAIL")
    fi
  fi
done

echo
echo "=== Run Summary ==="
pad() { printf "%-60s" "$1"; }
for i in "${!STATUS_KEYS[@]}"; do
  s=${STATUS_KEYS[$i]}
  st=${STATUS_VALS[$i]}
  pad "$s"; echo " $st"
done

echo
echo "Logs (stderr) for failures are at figs/.../*.err next to scripts."
