#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=360G
#SBATCH --time=06:00:00
#SBATCH --job-name=amr-queue
# One array task of a queue file on one node (run_queue.py):
#   sbatch --array=0-59 queue.sbatch benchmarks/campaign/commands/queue_main.txt
source "${AMR_CAMPAIGN_ROOT:?export AMR_CAMPAIGN_ROOT}/source/benchmarks/campaign/env.sh"
export ROOT PY
"$PY" benchmarks/campaign/run_queue.py "${1:?queue file}" "$SLURM_ARRAY_TASK_ID" --cores "$SLURM_CPUS_PER_TASK"
