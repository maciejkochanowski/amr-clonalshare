# MIC interval calibration: protocol

This folder holds the simulation drivers and task tables behind the calibration
of the MIC variance-ratio interval. The batch scripts and the order in which
they run are in `benchmarks/campaign/`; results, receipts and job accounting
are in `benchmarks/results_mic_release/`, and `account_jobs.py` reads the Slurm
records of the campaign and writes `accounting.json`, the CPU-hours the article
quotes.

## Target and model
rho = tau^2/(tau^2+sigma^2) in a Gaussian random-intercept model with independent
lineages, independent homoscedastic residuals and group sizes unrelated to the
effects. Readings are exact log2 values or intervals from a declared panel.

## Methods
- Exact readings: inversion of Wald's F pivot (`mic_inference.exact_gaussian_interval`).
- Interval readings: likelihood-ratio test at each candidate rho, with the mean,
  total scale and any fixed effects refitted under that rho and 199 datasets
  simulated from the refitted model with the same sizes, panels and covariate
  levels (`mic_null_bootstrap`). The statistic of every simulated dataset is
  computed by the same function as the observed one: both likelihoods are
  maximised over the same bounded nuisance box, from the same starting rules.
  A simulated dataset in which a covariate level fell wholly beyond the panel
  is scored like any other: the coefficient of that level stays within the
  box, where the observed data would be refused as unidentified.
  p = (1 + #{T* >= T})/(B + 1); a simulated dataset whose statistic cannot be
  computed counts as T* = infinity, so it can retain a candidate but never
  reject it. The interval inverts this test on a grid with bisection.

## Design
`calibrate.py` with `tasks_confirmation_exact.tsv` (root seed 20260915071): 24
Gaussian exact-reading designs with 5,000 datasets each, 6 stress designs and 12
singleton-only controls with 1,000 each.
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

## Additional designs
Four features of recorded MIC collections lie outside the main designs and are
studied in separate designs with the same decision rule, reported apart from the
main claim. `calibrate_extension.py` with `tasks_extension.tsv` (root seed
20260917401) has 51 designs with 500 datasets each: 24 heavy-censoring designs
(narrow panel, mean shifted so that 60% or 75% of readings fall on an end well, 8 or
30 lineages, uneven or half-singleton sizes, rho 0.1, 0.5 or 0.9); 6 designs with a
two-component residual (30% of readings shifted by 2.5 residual standard
deviations, standardised to unit total variance, wide or narrow panel); 18
laboratory designs (two laboratories differing by a fixed shift of one total
standard deviation, 85% of each lineage's isolates in one of them, analysed pooled
without a laboratory term, restricted to one laboratory, or pooled with the
laboratory as a fixed effect); and 3 designs that add a year of isolation, drawn
independently of lineage and laboratory from three values with a drift of 0.3 total
standard deviations per year, adjusted for both the laboratory and the year. The
full interval is computed for the first dataset of each design, except in design 20
(75% on an end well, 30 uneven lineages, rho 0.9) and in the three designs with two
covariates, where a single interval takes several hours; the fifth column of the
task table records this.

The *S. suis* reanalysis reads the panel within each testing laboratory
(`prepare_empirical.py`) with the country of isolation as a fixed effect (the reported
analysis), without it, for the isolates of one laboratory only (`--country "United
Kingdom"`), and, with the fixed effect, with 499 and 999 simulated datasets per
candidate for three agents (`run_empirical.py`).

## Reproduce
`benchmarks/campaign/CAMPAIGN.md` lists every command. `summarize.py` and
`summarize_null.py --tasks <task table>` check every result file against its
receipt before summarising.
