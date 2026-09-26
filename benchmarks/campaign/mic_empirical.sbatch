#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --job-name=mic-empirical
# One S. suis agent per array task, from the inputs prepare_empirical.py wrote:
#   sbatch --array=0-15 mic_empirical.sbatch INPUTS BOOTSTRAP OUTPUT [COVARIATE]
# INPUTS and OUTPUT are directories under $AMR_CAMPAIGN_ROOT/outputs; the
# optional fourth argument names the input column whose levels enter the
# model as fixed effects.
source "${AMR_CAMPAIGN_ROOT:?export AMR_CAMPAIGN_ROOT}/source/benchmarks/campaign/env.sh"
INPUTS=${1:?inputs directory}
BOOT=${2:?simulated datasets per tested rho}
OUTDIR=${3:?output directory}
COVARIATE=${4:-}
"$PY" benchmarks/mic_inference/run_empirical.py --manifest "$ROOT/outputs/$INPUTS/manifest.json" \
  --agent-index "$SLURM_ARRAY_TASK_ID" --output "$ROOT/outputs/$OUTDIR" \
  --workers "$SLURM_CPUS_PER_TASK" --bootstrap "$BOOT" \
  ${COVARIATE:+--covariate-column "$COVARIATE"}
