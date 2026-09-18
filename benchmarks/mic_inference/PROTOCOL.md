# MIC interval calibration: protocol and results

This folder holds the simulation drivers, task tables and Slurm scripts behind the
calibration of the MIC variance-ratio interval. Results, receipts and job accounting
are in `benchmarks/results_mic_release/`; `account_jobs.py` reads the Slurm records of
the campaigns and writes `accounting.json`, the CPU-hours the article quotes.

## Target and model
rho = tau^2/(tau^2+sigma^2) in a Gaussian random-intercept model with independent
lineages, independent homoscedastic residuals and group sizes unrelated to the
effects. Readings are exact log2 values or intervals from a declared panel.

## Methods
- Exact readings: inversion of Wald's F pivot (`mic_inference.exact_gaussian_interval`).
- Interval readings: likelihood-ratio test at each candidate rho, with the mean and
  total scale refitted under that rho and 199 datasets simulated from the refitted
  model with the same sizes and panel (`mic_null_bootstrap`). p = (1 + #{T* >= T})/(B + 1);
  a failed simulated fit counts as T* = infinity. The interval inverts this test on a
  grid with bisection. Simulated fits may end on the bounded nuisance box; the
  observed-data fit must be interior.

## Design
`calibrate.py` with `tasks_confirmation_exact.tsv` (root seed 20260915071): 24 Gaussian
exact-reading designs with 5,000 datasets each, 6 stress designs and 4 singleton-only
controls with 1,000 each.
`calibrate_null.py` with `tasks_null_confirmation.tsv` (root seed 20260916301): 48
Gaussian panel designs with 1,000 datasets each, 16 designs with other mean, scale and
rho with 500 each, and 12 stress designs with 200 each. The full interval is also
computed for every 500th dataset.
`prepare_null_campaign.py` writes the task table without reading any result.

## Decision rule
For panel readings the primary measure is the share of datasets in which the test at
the true rho keeps that value. Datasets without a result are counted separately. A
design passes when the upper Wilson 95% limit of its coverage reaches 0.95. Stress
designs describe sensitivity to non-normal residuals and are not part of the
calibration claim.

## Pilot and amendment
A pilot showed that taking the critical value from the unrestricted fit undercovered
(909 of 1,000 in the 8-lineage, half-singleton, wide-panel design at rho = 0.5), so the
test refits the nuisance parameters at every candidate. Before the final run, the
S. suis reanalysis returned the full range for erythromycin and tylosin because 14 and
18 of 199 simulated fits at rho = 0.995 ended on the nuisance box. Box maxima are now
accepted in the simulated fits. This can only turn an infinite simulated statistic
into a finite one. The full panel study was then run once on the final code with the
seed above, and only that run is reported.

## Extension after the main run
Three features of recorded MIC collections lie outside the main designs: panels on
which most readings fall on an end well, a two-component residual distribution
within lineages, and a laboratory effect partly aligned with lineage.
`calibrate_extension.py` with `tasks_extension.tsv` (root seed 20260917401) adds 45
designs with 500 datasets each, run after the main results were known and reported
separately: 24 heavy-censoring designs (narrow panel, mean shifted so that 60% or 75%
of readings fall on an end well, 8 or 30 lineages, uneven or half-singleton sizes,
rho 0.1, 0.5 or 0.9), 6 designs with a two-component residual (30% of readings
shifted by 2.5 residual standard deviations, standardised to unit total variance,
wide or narrow panel), and 18 laboratory designs (two laboratories differing by a
fixed shift equal to the total standard deviation, 85% of each lineage's isolates in one
of them, analysed pooled without a laboratory term, restricted to one laboratory, or
pooled with the laboratory as a fixed effect of the model). The same decision rule
applies. The loop, the outputs and the receipts are those of `calibrate_null.py`.
The six fixed-effect designs (cells 42 to 47) were added after the pooled designs had
shown that an analysis without the term is centred on the between-lineage share of the
mixed collection; they were run with the implementation of the fixed effect, whose
code path without a covariate is unchanged (the likelihood and its gradient agree to
machine precision with the earlier implementation), so the earlier cells were not
rerun. Three further designs (cells 48 to 50; narrow panel, 30 uneven lineages, rho
0.1, 0.5 or 0.9) add a second covariate to the adjusted laboratory design: a year of
isolation drawn independently of lineage and laboratory from three values with a
drift of 0.3 total standard deviations per year, and the analysis adjusts for both
the laboratory and the year. They were run after the several-covariate
implementation was in place, in a root of their own (`amr_clonalshare_ext4_20260917`,
same root seed), and are summarised separately
(`null_extension_summary_country_year.csv`). In these three designs the full interval
of replicate 0 had not finished after 100 minutes with two covariates, so the chunks
holding replicates 0 to 24 were rerun without the full interval (`--bounds-every 0`,
one worker per chunk, `helios_extension.sbatch` otherwise unchanged); the three designs
therefore have the test at the true rho for all 500 datasets and no full interval. The receipts record the source hashes of each run. Two tasks of the heavy-censoring
designs exceeded the two-hour limit of the array and their unfinished chunks were rerun
in smaller pieces; for design 20 (75% on an end well, 30 uneven lineages, rho 0.9) the
full interval of replicate 0 did not finish within a further six hours and was left out,
so that design has the test at the true rho for all 500 datasets and no full interval.

The S. suis reanalysis was repeated with the panel inferred within each testing
laboratory (`prepare_empirical.py`), with the country of isolation as a fixed effect
(`--covariate-column country`; the reported analysis), without it, for the isolates of
one laboratory only (`--country "United Kingdom"`), and, with the fixed effect, with 499
and 999 simulated datasets per candidate for three agents
(`helios_empirical_extension.sbatch`).

## Changes after the runs that do not affect results
Worker processes limit BLAS threads to one, the pipeline counts only the CPUs the
process may use, and some docstrings changed. Receipts therefore list earlier hashes
for a few source files than the released files.

## Reproduce
Adapt the paths in `helios.sbatch`, `helios_null.sbatch`, `helios_empirical.sbatch`,
`helios_extension.sbatch` and `helios_empirical_extension.sbatch`, then run
`summarize.py`, `summarize_null.py --tasks tasks_null_confirmation.tsv` and, on the
extension outputs, `summarize_null.py` with the first 210 rows of `tasks_extension.tsv`
for designs 0 to 41, rows 211 to 240 for designs 42 to 47 and rows 241 to 255 for
designs 48 to 50 (the three runs have different source hashes, which the summariser
refuses to mix; the summaries are
`benchmarks/results_mic_release/extension/null_extension_summary_main.csv`,
`null_extension_summary_adjusted.csv` and `null_extension_summary_country_year.csv`),
which check every result file against its receipt before summarising.
