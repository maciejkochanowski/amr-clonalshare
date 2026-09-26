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
  reject it. A simulated dataset that the analysis would refuse outright
  (every reading in one interval, or one lineage left with readings) is
  replaced by a fresh draw, up to 20 B times, since the observed dataset is one
  the analysis reports on; the number replaced is recorded. The interval
  inverts this test on a grid with bisection.
- Check of the Gaussian model (`_mic_shape`): at the estimated rho, the deviance
  of the readings of every panel and covariate level against the fitted normal
  law, and the Pearson statistic of every lineage's readings against their
  prediction given the rest of the lineage, both calibrated by the simulated
  datasets of the test at the estimate; rejected when twice the smaller
  p-value is at most 0.05.

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
task table records this. A design whose first dataset is refused has no full
interval.

The *S. suis* reanalysis reads the panel within each testing laboratory
(`prepare_empirical.py`) with the country of isolation as a fixed effect (the reported
analysis), without it, for the isolates of one laboratory only (`--country "United
Kingdom"`), and, with the fixed effect, with 499 and 999 simulated datasets per
candidate for three agents (`run_empirical.py`).

## Model check, several panels and mixed readings
These designs are reported apart from the calibration claim.

- Check of the Gaussian model (`shape_check.py`), 400 datasets per design, the
  datasets of the designs above (same generators, seeds and replicate numbers)
  so that each joins its coverage row: every interval-reading design of
  `calibrate_null.py` except the all-singleton controls, every design of
  `calibrate_extension.py`, and the plasmode designs below. The rejection rate
  in the Gaussian designs is the size of the check; the check is kept if no
  Gaussian design has a lower Wilson 95% limit of its rejection rate above
  0.05. The rejection rates under t4, contaminated-normal and two-component
  residuals, and in the mixture designs, are its power. Coverage of the test is
  also reported among the datasets the check passes.
- Several panels (`calibrate_panels.py --mode multipanel`, root seed 20261001),
  12 designs with 1,000 datasets each: 30 lineages of uneven size read by
  three laboratories (a wide panel with cut points -4 to 4, a narrow panel from
  -2 to 1, one cut point at 0), laboratory membership aligned with lineage
  (85%) or random, the laboratories shifted by -0.5, 0 and +0.5 total standard
  deviations or not, rho 0.1, 0.5 or 0.9. Every dataset is analysed with each
  laboratory on its own panel (the shift, where present, as a fixed effect) and
  on the records of the wide-panel laboratory alone. The full interval is
  computed for every 50th dataset. Decision rule as for the main designs.
- Plasmode (`calibrate_panels.py --mode plasmode`, root seed 20261002), 500
  datasets per *S. suis* agent: readings simulated from the Gaussian model
  fitted to the agent (country as a fixed effect), on its recorded lineages,
  laboratories, panels and countries, tested at the fitted rho. Decision rule
  as for the main designs.
- Wild-type / non-wild-type mixtures (`calibrate_mixture.py`, root seed
  20261003), 8 designs with 500 datasets each, on the lineage sizes of the
  *S. suis* collection: a latent log2 MIC equal to the wild-type mode (-5.5,
  below the panel, or -1.0, inside it) plus 5 dilutions for non-wild-type
  isolates, a lineage effect (SD 0.3) and a residual (SD 0.8); non-wild-type
  membership from a liability threshold model with prevalence 0.3 or 0.6 and
  liability ICC 0.3 or 0.7; one panel with cut points -5 to 3. Each dataset
  records the calibrated test at the lineage share of latent MIC variance, the
  check of the Gaussian model, and the population interval on the
  non-wild-type calls at the liability ICC, with its mixing check. The
  designs are reported and do not decide. The study is void if the lineage
  share of latent variance, computed by quadrature, differs by more than 0.002
  from a Monte Carlo evaluation of the same formula on 2,000,000 lineages.

## Reproduce
`benchmarks/campaign/CAMPAIGN.md` lists every command. `summarize.py` and
`summarize_null.py --tasks <task table>` check every result file against its
receipt before summarising.
