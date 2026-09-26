# Population liability ICC: protocol

This folder holds the calibration and validation of the optional population
model (`population_probit` in a configuration, `population_probit_icc` and
`general_probit_icc` in Python). The rules below were fixed before any result
was computed; `design.py` holds every cell and every seed.

## Target and model

The target is the liability intraclass correlation rho of a threshold model:
an isolate of lineage g has liability u_g + e, with u_g ~ N(0, rho) and
e ~ N(0, 1 - rho) independent, and is positive when the liability exceeds the
(1 - p) quantile of a standard normal, so that p is the marginal prevalence.
Under this model rho does not depend on p, which is the reason to report it
beside the observed-scale lineage share, whose value does.

`design.draw` generates data from that statement, by thresholding simulated
liabilities; it does not use the package's likelihood.

## Methods

- Fixed cut-off (`fixed_cutoff`, the default when the model is enabled): the set of
  rho whose profile likelihood-ratio statistic does not exceed a fixed
  critical value c.
- General method (`interval_method: general`): the fitted-nuisance bootstrap
  hull with 4,999 reference draws and its fixed internal threshold of 0.04;
  a one-sided upper limit when only one lineage has repeated isolates.

## Families

- **Calibration** (216 cells, 2,000 datasets each, seed 20260926101): the
  Gaussian model with 10, 30 or 100 lineages of 5 or 20 isolates or of unequal
  sizes (1, 2, 5, 10, 20, 50 repeated), p = 0.05, 0.15, 0.30, 0.50 and rho = 0,
  0.05, 0.2, 0.4, 0.6, 0.8. The model is symmetric in p and 1 - p.
- **Validation** (144 cells, 2,000 datasets each for the fixed cut-off and 500
  for the general method, seed 20260926202): the Gaussian model off the
  calibration grid, with 15, 50 or 200 lineages of 10 isolates or of skewed
  sizes (a few large lineages and a tail of singletons), p = 0.02, 0.10, 0.25,
  0.85 and rho = 0.1 to 0.9; and the 30 lineage sizes of the *S. suis* example
  at p = 0.04, 0.24, 0.54, 0.85 and rho = 0 to 0.9.
- **Robustness** (12 cells, same replication, seed 20260926303): one assumption
  broken at a time, 50 lineages, p = 0.10 or 0.25, rho = 0.3 or 0.7: lineage
  effects from a t distribution on four degrees of freedom, from a two-point
  (carrier) law with 20% of lineages carrying the effect, both scaled by their
  theoretical variance, and lineage sizes that grow with the lineage effect.
- **Benefit** (12 cells, 500 datasets each, seed 20260926404): 30 lineages of 20
  isolates, rho = 0.3 or 0.6, p = 0.02, 0.05, 0.10, 0.20, 0.35, 0.50. Each
  dataset is read by the profile estimate and by the observed-scale
  lineage-membership share of `clonal_share` (5 folds, 5 repeats, 49
  permutations).
- **Anchor** (300 datasets from the validation cells, seed 20260926505):
  `anchor.py` recomputes the log-likelihood at the fitted point and at the
  generating rho by adaptive quadrature written from the model statement.

- **Mixing check** (validation and robustness cells, 500 datasets each, the
  datasets of the fixed-cut-off arrays): every dataset records the
  fixed-cut-off interval and the package's check of the Gaussian law of
  lineage effects against a nonparametric one (`_population_mixing`, 99
  simulated datasets). The rejection rate in the Gaussian validation cells is
  its size, the rates in the robustness cells its power, and coverage is
  reported among the datasets the check passes.

## Rules

1. The critical value c is the largest, over the calibration cells, of the
   95th percentile (the ceil(0.95 n)-th order statistic) of the
   likelihood-ratio statistic at the generating rho. A dataset whose set is
   [0, 1] because the outcome is constant or fewer than two lineages repeat
   counts as a statistic of zero; a failed fit counts as infinite.
2. A dataset is covered when its reported set contains the generating rho;
   [0, 1] covers; a failed fit does not. A cell passes when the upper Wilson 95%
   limit of its coverage reaches 0.95, the rule of the MIC interval.
3. Void conditions: the study is void if the anchor differs from the package
   by more than 1e-6 in any log-likelihood, or if any calibration cell has more
   than 0.5% failed fits.
4. Decision, fixed in advance. The population model is described in the
   article if (a) every Gaussian validation cell passes for the fixed cut-off,
   and (b) in each benefit row with p from 0.05 to 0.50 the mean profile
   estimate lies within 0.05 of rho and varies less across p than the mean
   observed-scale share. If (a) fails, the population model is removed from the
   release. The general method is kept only if every Gaussian validation cell
   passes for it as well. Robustness cells are reported and do not decide. The mixing check is kept
   if no Gaussian validation cell has a lower Wilson 95% limit of its
   rejection rate above 0.05.

## Reproduce

The commands are in `benchmarks/campaign/CAMPAIGN.md`. `run_cell.py` runs
replicates of one cell, `aggregate.py` applies the rules above and writes
`src/amr_clonalshare/population_probit_validation.json`, and `anchor.py` runs
the anchor.
