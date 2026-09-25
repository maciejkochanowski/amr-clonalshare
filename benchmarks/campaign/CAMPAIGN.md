# The validation campaign

Every result under `benchmarks/results_*` and every validation file shipped
inside the package (`validation_grid.json`, `population_probit_validation.json`,
`general_probit_protocol.json`) is produced by the commands below, in this
order, from one commit of this repository on one software stack. Each step
writes a receipt (`run_logged.py`) naming the command, the commit, a digest of
every Python file of the package, the interpreter and library versions, the
Slurm job and the sha256 of every file it wrote. The batch scripts take the
account and partition from the command line and read nothing outside the
campaign directory.

## 0. Environment

```bash
export AMR_CAMPAIGN_ROOT=/path/to/campaign        # holds source/, env/, outputs/
git clone <repository> "$AMR_CAMPAIGN_ROOT/source"
git -C "$AMR_CAMPAIGN_ROOT/source" checkout <commit>
uv python install 3.11.15
uv venv --python 3.11.15 "$AMR_CAMPAIGN_ROOT/env"
VIRTUAL_ENV="$AMR_CAMPAIGN_ROOT/env" uv pip install -r "$AMR_CAMPAIGN_ROOT/source/requirements-lock.txt"
VIRTUAL_ENV="$AMR_CAMPAIGN_ROOT/env" uv pip install --no-deps -e "$AMR_CAMPAIGN_ROOT/source"
cd "$AMR_CAMPAIGN_ROOT/source"
S="sbatch -A <account> -p <partition>"
C=benchmarks/campaign
```

## 1. Gate calibrations and grids

```bash
$S --array=0-4   --cpus-per-task=1  --mem=4G  --time=08:00:00 $C/lines.sbatch $C/commands/serial_runs.txt
$S --array=0-1   --cpus-per-task=24 --mem=24G --time=08:00:00 $C/lines.sbatch $C/commands/parallel_runs.txt
$S --array=0-134 --cpus-per-task=1  --mem=4G  --time=24:00:00 $C/lines.sbatch $C/commands/estimator_grid.txt
$S --array=0-58  --cpus-per-task=1  --mem=4G  --time=24:00:00 $C/lines.sbatch $C/commands/realised_calibration.txt
$S --array=0-23  --cpus-per-task=16 --mem=16G --time=08:00:00 $C/lines.sbatch $C/commands/repeated_looks.txt
```

The constants the package reads off these results (`SUPPORT_THRESHOLD`,
`FEW_LINEAGES`, `KURTOSIS_LIMIT`, `DEFAULT_MIN_SHARED_SUPPORT`,
`CENSORED_GROUP_LIMIT`, `SINGLE_CUT_PREVALENCE`) are checked against them
before step 2; a change of a constant is a new commit, and step 1 is repeated
on it.

## 2. MIC interval and population model

```bash
$S --array=0-41   $C/mic_exact.sbatch
$S --array=0-583  $C/mic_null.sbatch
$S --array=0-254  $C/mic_extension.sbatch
$AMR_CAMPAIGN_ROOT/env/bin/python benchmarks/mic_inference/prepare_empirical.py \
    --output "$AMR_CAMPAIGN_ROOT/outputs/empirical_inputs"
$AMR_CAMPAIGN_ROOT/env/bin/python benchmarks/mic_inference/prepare_empirical.py \
    --country "United Kingdom" --seed-base 20260917200 \
    --output "$AMR_CAMPAIGN_ROOT/outputs/empirical_inputs_uk"
$S --array=0-15  $C/mic_empirical.sbatch empirical_inputs    199 mic_empirical country
$S --array=0-15  $C/mic_empirical.sbatch empirical_inputs    199 empirical_unadjusted
$S --array=0-15  $C/mic_empirical.sbatch empirical_inputs_uk 199 empirical_uk
$S --array=2,5,9 $C/mic_empirical.sbatch empirical_inputs    499 empirical_b499 country
$S --array=2,5,9 $C/mic_empirical.sbatch empirical_inputs    999 empirical_b999 country

$S --array=0-179 --cpus-per-task=48 --mem=96G --time=24:00:00 $C/lines.sbatch $C/commands/population_calibration_packed.txt
$S --array=0-324 --cpus-per-task=48 --mem=180G --time=48:00:00 $C/lines.sbatch $C/commands/population_validation_general_packed.txt
$S --array=0-6   --cpus-per-task=48 --mem=96G --time=24:00:00 $C/lines.sbatch $C/commands/population_benefit_packed.txt
$S --array=0     --cpus-per-task=1 --mem=4G  --time=12:00:00 $C/lines.sbatch $C/commands/population_anchor.txt
```

Each line of a `_packed` command file runs 48 replicate chunks at once, one per
core, and fails if any of them fails; `make_population_commands.py` writes them.

When the calibration has finished:

```bash
O=$AMR_CAMPAIGN_ROOT/outputs/population
$PY -m benchmarks.population_model.aggregate calibrate $O/lr_calibration $O/calibration.json
python -c "import json;print(json.load(open('$O/calibration.json'))['critical_value'])" > $O/critical_value.txt
$S --array=0-129 --cpus-per-task=48 --mem=96G --time=24:00:00 $C/lines.sbatch $C/commands/population_validation_lr_packed.txt
```

## 3. Summaries

Each summary runs under `run_logged.py` with the directory it reads as an
output, so its receipt carries the sha256 of every file the array wrote.

```bash
R="$PY benchmarks/campaign/run_logged.py"
O=$AMR_CAMPAIGN_ROOT/outputs
$R --receipt $O/estimator_grid/RUN_RECEIPT.json --outputs $O/estimator_grid -- \
    $PY benchmarks/estimator_benchmark.py --aggregate --out $O/estimator_grid
$R --receipt $O/realised_calibration/RUN_RECEIPT.json --outputs $O/realised_calibration -- \
    $PY benchmarks/realised_calibration.py --aggregate --out $O/realised_calibration
$R --receipt $O/mic_exact/RUN_RECEIPT.json --outputs $O/mic_exact -- \
    $PY benchmarks/mic_inference/summarize.py $O/mic_exact $O/exact_confirmation_summary.csv
$R --receipt $O/mic_null/RUN_RECEIPT.json --outputs $O/mic_null -- \
    $PY -m benchmarks.mic_inference.summarize_null $O/mic_null $O/null_confirmation_summary.csv \
    --tasks benchmarks/mic_inference/tasks_null_confirmation.tsv
$R --receipt $O/mic_extension/RUN_RECEIPT.json --outputs $O/mic_extension -- \
    $PY -m benchmarks.mic_inference.summarize_null $O/mic_extension $O/null_extension_summary.csv \
    --tasks benchmarks/mic_inference/tasks_extension.tsv
for F in validation robustness; do for M in lr general; do
  $R --receipt $O/population/${M}_$F/RUN_RECEIPT.json --outputs $O/population/${M}_$F -- \
      $PY -m benchmarks.population_model.aggregate cover $O/population/${M}_$F $O/population/${F}_$M.csv
done; done
$R --receipt $O/population/benefit/RUN_RECEIPT.json --outputs $O/population/benefit -- \
    $PY -m benchmarks.population_model.aggregate benefit $O/population/benefit $O/population/benefit.csv
$PY benchmarks/mic_inference/account_jobs.py $O/accounting.json <pytest log> <name>=<job> ...
```

## 4. On one workstation, from the same commit and stack

```bash
python benchmarks/estimator_grid_summary.py <estimator_grid>/estimator_benchmark.json \
    src/amr_clonalshare/validation_grid.json
python -m benchmarks.population_model.aggregate evidence <population results> src/amr_clonalshare
python -m benchmarks.mic_inference.write_extension_grid benchmarks/results_mic_release/extension/extension_design_grid.json
python scripts/update_golden.py
amr-clonalshare --config examples/ssuis/config.yaml --results-dir examples/ssuis/expected --overwrite
amr-clonalshare --config examples/salmonella_poultry/config.yaml --results-dir examples/salmonella_poultry/expected_serovar --overwrite
amr-clonalshare --config examples/salmonella_poultry/config_cluster.yaml --results-dir examples/salmonella_poultry/expected_cluster --overwrite
python benchmarks/empirical/ssuis_analysis.py
python benchmarks/empirical/salmonella_analysis.py
OMP_NUM_THREADS=1 python benchmarks/seed_stability.py 40 18 benchmarks/results_seed_stability/seeds.json
OMP_NUM_THREADS=1 python benchmarks/profile_run.py examples/ssuis/config.yaml benchmarks/results_profile/profile.json
```

each under `run_logged.py` as above.
