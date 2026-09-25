#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=amr-lines
# Runs line SLURM_ARRAY_TASK_ID+1 of a command file from the source directory;
# cores, memory and time are given on the sbatch command line:
#   sbatch --array=0-134 --cpus-per-task=1 --mem=4G --time=12:00:00 \
#       lines.sbatch benchmarks/campaign/commands/estimator_grid.txt
source "${AMR_CAMPAIGN_ROOT:?export AMR_CAMPAIGN_ROOT}/source/benchmarks/campaign/env.sh"
COMMANDS=${1:?command file}
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$COMMANDS")
[ -n "$LINE" ] || { echo "no line $((SLURM_ARRAY_TASK_ID+1)) in $COMMANDS" >&2; exit 2; }
echo "$LINE"
eval "$LINE"
