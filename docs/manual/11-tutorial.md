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

The example configuration reads the hierBAPS population cluster as the lineage, the year of collection as the intake of a programme of looks, the testing laboratory as the panel column and the country of isolation as a covariate of the MIC reading and as the stratifying column. The laboratory and the country define the strata within which the MIC ordering is read. The run writes the files listed in the [results manual](05-results.md).

```
amr-clonalshare --config examples/ssuis/config.yaml --results-dir out/ssuis
```

Open `out/ssuis/report.html`. Table 1 gives the lineage share of every binary endpoint with its interval, p-value, q-value, e-value, selection and reading; Table 3 the bounds every call places on the share of the latent ordering, with their one-sided 95% lower limit; Table 4b the running e-value over the yearly intakes; Table 5 the lineage share of the MIC ordering per agent with its interval, p-value, q-value, selection across agents, the bounds on the share of the latent MIC ordering and the share of readings on an end well; Table 5a the panel of each laboratory; the heat map the readings per lineage and dilution well; and Table 5d the same shares within each country.

## Step 4. Compare two lineage definitions on matched records

The comparison tool reads one row per isolate with a 0/1 outcome and two lineage columns. A short script prepares this table for ceftiofur:

```python
import pandas as pd
d = "examples/ssuis/data/"
meta = pd.read_csv(d + "metadata.csv", dtype=str)
calls = pd.read_csv(d + "calls_long.csv", dtype=str)
cef = calls[calls.antibiotic == "ceftiofur"].merge(meta, on="genome_id")
cef["positive"] = cef.call.map({"NWT": 1, "WT": 0}).astype("Int64")  # missing stays empty
cef[["genome_id", "positive", "baps_cluster", "mlst"]].to_csv("aligned.csv", index=False)
```

```
amr-clonalshare-compare --input aligned.csv --id-column genome_id --outcome positive \
  --lineage-a baps_cluster --lineage-b mlst --seed 20260915 --permutations 199 \
  --bootstraps 399 --output out/ceftiofur_comparison
```

`arms.csv` lists the six analyses: population clusters on all 677 isolates, on the 458 shared records and on the 393 isolates in a repeated lineage under both definitions, and sequence types on those 393, on the shared records and on their own. `comparison.json` splits the difference between the first and the last into five terms, with a paired 95% interval for every term (`paired_intervals`) for these records and labels, and `report.md` sets them out in one table. The sequence types reach a support of 85.81% of the typed isolates, but only 43 of the 108 types repeat. The tool also takes `--folds` (5) and `--repeats` (20) for its cross-validated scores beside `--permutations` and `--bootstraps`.

## Step 5. Read the results

Start from `results.csv`: check the status and reason of every row, the analysed records and the interval kind, and only then compare estimates between agents, strata or lineage definitions (see the [reporting checklist](12-reporting-checklist.md)).
