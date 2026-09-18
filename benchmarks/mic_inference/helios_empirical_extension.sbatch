#!/bin/bash
#SBATCH --account=plgsoa27eupahw-cpu
#SBATCH --partition=plgrid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --job-name=mic-empirical-ext
# S. suis reanalyses beyond the main run: a single-country, single-laboratory
# subset, and more simulated datasets per tested value of rho.
#   sbatch --array=0-15 helios_empirical_extension.sbatch empirical_inputs_uk 199 empirical_uk
#   sbatch --array=2,5,9 helios_empirical_extension.sbatch empirical_inputs 999 empirical_b999
#   sbatch --array=0-15 helios_empirical_extension.sbatch empirical_inputs 199 empirical_adjusted panel
# A fourth argument names the input column whose levels enter the model as
# fixed effects (the testing laboratory, held in the panel column).
set -euo pipefail
ROOT=${MIC_CAMPAIGN_ROOT:?campaign root}
ENV=/net/storage/pr3/plgrid/plggzchsdbbpiwet/amr_clonalshare_final_20260915/env-offline
INPUTS=${1:?manifest directory under the campaign root}
BOOT=${2:?simulated datasets per tested rho}
OUTDIR=${3:?output directory under the campaign root}
COVARIATE=${4:-}
export PYTHONPATH="$ROOT/source/src"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export TMPDIR="$ROOT/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR" "$ROOT/outputs/$OUTDIR"
cd "$ROOT"
"$ENV/bin/python" source/benchmarks/mic_inference/run_empirical.py --manifest "$INPUTS/manifest.json" \
  --agent-index "$SLURM_ARRAY_TASK_ID" --output "outputs/$OUTDIR" --workers "$SLURM_CPUS_PER_TASK" --bootstrap "$BOOT" \
  ${COVARIATE:+--covariate-column "$COVARIATE"}
