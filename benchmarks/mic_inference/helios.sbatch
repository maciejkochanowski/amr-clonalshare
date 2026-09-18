#!/bin/bash
#SBATCH --account=plgsoa27eupahw-cpu
#SBATCH --partition=plgrid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --job-name=clonalshare-mic-ci
set -euo pipefail
ROOT=/net/storage/pr3/plgrid/plggzchsdbbpiwet/amr_clonalshare_mic_inference_20260915
ENV=/net/storage/pr3/plgrid/plggzchsdbbpiwet/amr_clonalshare_final_20260915/env-offline
PHASE=${1:?pilot or confirmation}
TASKS=${2:?task table}
read -r CELL START REPS BOOT < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$TASKS")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONPATH="$ROOT/source/src:$ROOT/source"
# User scratch quota is exhausted; use only this campaign's project storage.
export TMPDIR="$ROOT/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"
cd "$ROOT/source"
"$ENV/bin/python" benchmarks/mic_inference/calibrate.py --cell "$CELL" --start "$START" \
  --replicates "$REPS" --bootstrap "$BOOT" --phase "$PHASE" \
  --output "$ROOT/outputs/$PHASE/cell_${CELL}_${START}.csv"
