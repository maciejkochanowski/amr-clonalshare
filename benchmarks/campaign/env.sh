# Sourced by every batch script of the campaign. The campaign directory holds
# source/ (a checkout of the release commit), env/ (a virtual environment
# built from requirements-lock.txt), outputs/ and tmp/; nothing else is read.
set -euo pipefail
ROOT=${AMR_CAMPAIGN_ROOT:?export AMR_CAMPAIGN_ROOT, the campaign directory}
PY=${AMR_CAMPAIGN_ENV:-$ROOT/env}/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONPATH="$ROOT/source/src:$ROOT/source"
export PYTHONDONTWRITEBYTECODE=1
export TMPDIR="$ROOT/tmp/${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$TMPDIR" "$ROOT/outputs"
cd "$ROOT/source"
