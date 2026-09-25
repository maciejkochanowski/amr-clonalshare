#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=04:00:00
#SBATCH --job-name=mic-exact
# One row of tasks_confirmation_exact.tsv per array task:
#   sbatch --array=0-41 mic_exact.sbatch
source "${AMR_CAMPAIGN_ROOT:?export AMR_CAMPAIGN_ROOT}/source/benchmarks/campaign/env.sh"
TASKS=benchmarks/mic_inference/tasks_confirmation_exact.tsv
read -r CELL START REPS BOOT < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$TASKS")
mkdir -p "$ROOT/outputs/mic_exact"
"$PY" benchmarks/mic_inference/calibrate.py --cell "$CELL" --start "$START" \
  --replicates "$REPS" --bootstrap "$BOOT" \
  --output "$ROOT/outputs/mic_exact/cell_${CELL}_${START}.csv"
