#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --job-name=mic-extension
# One row of tasks_extension.tsv per array task, its replicates split over
# four worker processes; see calibrate_extension.py for the designs:
#   sbatch --array=0-254 mic_extension.sbatch
source "${AMR_CAMPAIGN_ROOT:?export AMR_CAMPAIGN_ROOT}/source/benchmarks/campaign/env.sh"
TASKS=benchmarks/mic_inference/tasks_extension.tsv
read -r CELL START REPS BOOT BOUNDS < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$TASKS")
OUT="$ROOT/outputs/mic_extension"
mkdir -p "$OUT"
PIDS=()
for WORKER in 0 1 2 3; do
  FIRST=$((START+WORKER*REPS))
  "$PY" benchmarks/mic_inference/calibrate_extension.py --cell "$CELL" --start "$FIRST" \
    --replicates "$REPS" --bootstrap "$BOOT" --bounds-every "$BOUNDS" \
    --output "$OUT/cell_${CELL}_${FIRST}.csv" &
  PIDS+=("$!")
done
FAILED=0
for PID in "${PIDS[@]}"; do
  if ! wait "$PID"; then FAILED=1; fi
done
exit "$FAILED"
