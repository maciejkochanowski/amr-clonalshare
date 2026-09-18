#!/bin/bash
#SBATCH --account=plgsoa27eupahw-cpu
#SBATCH --partition=plgrid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --job-name=mic-release-empirical
# Release 1.0.0 reanalysis of the S. suis agents; inputs from prepare_empirical.py.
set -euo pipefail
ROOT=${MIC_CAMPAIGN_ROOT:?campaign root}
ENV=/net/storage/pr3/plgrid/plggzchsdbbpiwet/amr_clonalshare_final_20260915/env-offline
export PYTHONPATH="$ROOT/source/src"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export TMPDIR="$ROOT/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR" "$ROOT/outputs/empirical"
cd "$ROOT"
"$ENV/bin/python" source/benchmarks/mic_inference/run_empirical.py --manifest empirical_inputs/manifest.json \
  --agent-index "$SLURM_ARRAY_TASK_ID" --output outputs/empirical --workers "$SLURM_CPUS_PER_TASK"
