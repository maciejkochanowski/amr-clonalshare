#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --job-name=mic-null
# One row of tasks_null_confirmation.tsv per array task, its replicates split
# over four worker processes:
#   sbatch --array=0-583 mic_null.sbatch
source "${AMR_CAMPAIGN_ROOT:?export AMR_CAMPAIGN_ROOT}/source/benchmarks/campaign/env.sh"
TASKS=benchmarks/mic_inference/tasks_null_confirmation.tsv
read -r CELL START REPS BOOT BOUNDS < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$TASKS")
OUT="$ROOT/outputs/mic_null"
mkdir -p "$OUT"
PIDS=()
for WORKER in 0 1 2 3; do
  FIRST=$((START+WORKER*REPS))
  "$PY" benchmarks/mic_inference/calibrate_null.py --cell "$CELL" --start "$FIRST" \
    --replicates "$REPS" --bootstrap "$BOOT" --bounds-every "$BOUNDS" \
    --output "$OUT/cell_${CELL}_${FIRST}.csv" &
  PIDS+=("$!")
done
FAILED=0
for PID in "${PIDS[@]}"; do
  if ! wait "$PID"; then FAILED=1; fi
done
exit "$FAILED"
