# Tutorial: one run on a public collection

The steps below reproduce the *Streptococcus suis* analyses of the release on a workstation. Commands are run from the root of the repository; the data, their provenance and the cut-offs are described in `examples/ssuis/DATA_PROVENANCE.md`.

## Step 1. Install the package and obtain the example data

```
python -m pip install amr-clonalshare==1.0.0
git clone --branch v1.0.0 https://github.com/maciejkochanowski/amr-clonalshare
cd amr-clonalshare
```

## Step 2. Check the input before any estimate is produced

The check writes `input_qc.md`, which lists the analysed identifiers, lineage sizes, lineage coverage, exclusions and outcome counts per agent.

```
amr-clonalshare --config examples/ssuis/config.yaml --check-input --results-dir out/ssuis_qc
```

## Step 3. Run the configured analysis

The example configuration names the testing laboratory as the panel column, the country of isolation as a fixed effect of the MIC model and as the stratifying column, and switches the calibrated MIC interval off (`censored: {calibrated_interval: false}`) to keep the run short, about fifteen minutes on one core with the three strata. Deleting that line adds the calibrated interval for all sixteen agents: 7 to 37 minutes per agent on 16 cores in the records of this release, a median of 23, and 5.4 hours for the whole panel. The run writes the files listed in the [results manual](05-results.md).

```
amr-clonalshare --config examples/ssuis/config.yaml --results-dir out/ssuis
```

Open `out/ssuis/report.html`. Table 1 gives the lineage share of every binary endpoint with its limits, e-value and reading; Table 5 the share of the dilution scale per agent; Tables 5a to 5c the panel of each laboratory, the fixed-effect coefficients with their reference level and the status of every calibrated interval; the heat map the readings per lineage and dilution well; and Table 5d the same shares within each country.

## Step 4. Compare two lineage definitions on matched records

The comparison tool reads one row per isolate with a 0/1 outcome and two lineage columns. A short script prepares this table for ceftiofur:

```python
import pandas as pd
d = "examples/ssuis/data/"
meta = pd.read_csv(d + "metadata.csv", dtype=str)
calls = pd.read_csv(d + "calls_long.csv", dtype=str)
cef = calls[calls.antibiotic == "ceftiofur"].merge(meta, on="genome_id")
cef["positive"] = cef.call.map({"non-susceptible": 1, "susceptible": 0}).astype("Int64")  # missing stays empty
cef[["genome_id", "positive", "baps_cluster", "mlst"]].to_csv("aligned.csv", index=False)
```

```
amr-clonalshare-compare --input aligned.csv --id-column genome_id --outcome positive \
  --lineage-a baps_cluster --lineage-b mlst --seed 20260915 --permutations 199 \
  --bootstraps 399 --output out/ceftiofur_comparison
```

`arms.csv` lists the four analyses: population clusters on all 677 isolates and on the 458 shared records, and sequence types on the shared records and on their own; the sequence-type score is flagged because its lineage coverage is 85.81%.

## Step 5. Compute the calibrated MIC interval for one agent

The interval tool reads a table with the columns `lo`, `hi` (log2 bounds) and `lineage`, and optionally a panel and covariate column per reading; the inputs of the release analysis are stored in the repository.

```
amr-clonalshare-mic benchmarks/results_mic_release/mic_empirical_inputs/ceftiofur.csv \
  --method bootstrap --seed 20260916102 --panel-column panel --covariate-column country \
  --panel-edges '{"LGC Fordham": [-5,-4,-3,-2,-1,0,1,2,3], "OUCRU Ho Chi Minh City": [-4,-3,-2]}' \
  --workers 16 --output out/ceftiofur_rho.json
```

The output records the estimate and the interval together with every tested value of ρ, its bootstrap sample, the fixed-effect coefficients, the input digest and the software versions; the computation takes about ten minutes on 16 cores. `--bootstrap` (default 199) sets the simulated datasets per tested value, `--alpha` (0.05) the level and `--tolerance` (0.002) the bisection step in ρ; `--method null-test --rho R` runs the test at one candidate value instead of inverting it, and `--method exact` applies Wald's pivot when every reading is exact. The comparison tool takes `--folds` (5) and `--repeats` (20) for its cross-validated scores beside `--permutations` and `--bootstraps`.

## Step 6. Read the results

Start from `results.csv`: check the status and reason of every row, the analysed records and the interval kind, and only then compare estimates between agents, strata or lineage definitions (see the [reporting checklist](13-reporting-checklist.md)).
