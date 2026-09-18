#!/bin/bash
#SBATCH --account=plgsoa27eupahw-cpu
#SBATCH --partition=plgrid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --job-name=mic-null-calibration
set -euo pipefail
ROOT=${MIC_CAMPAIGN_ROOT:-/net/storage/pr3/plgrid/plggzchsdbbpiwet/amr_clonalshare_mic_null_20260916}
ENV=/net/storage/pr3/plgrid/plggzchsdbbpiwet/amr_clonalshare_final_20260915/env-offline
PHASE=${1:?phase}
TASKS=${2:?task_table}
read -r CELL START REPS BOOT BOUNDS < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$TASKS")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONPATH="$ROOT/source/src:$ROOT/source"
export PYTHONDONTWRITEBYTECODE=1
export TMPDIR="$ROOT/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR" "$ROOT/outputs/$PHASE"
cd "$ROOT/source"
PIDS=()
for WORKER in 0 1 2 3; do
  FIRST=$((START+WORKER*REPS))
  "$ENV/bin/python" benchmarks/mic_inference/calibrate_null.py --cell "$CELL" --start "$FIRST" \
    --replicates "$REPS" --bootstrap "$BOOT" --root-seed "${MIC_ROOT_SEED:-20260916073}" --phase "$PHASE" --bounds-every "$BOUNDS" \
    --output "$ROOT/outputs/$PHASE/cell_${CELL}_${FIRST}.csv" &
  PIDS+=("$!")
done
FAILED=0
for PID in "${PIDS[@]}"; do
  if ! wait "$PID"; then FAILED=1; fi
done
exit "$FAILED"
