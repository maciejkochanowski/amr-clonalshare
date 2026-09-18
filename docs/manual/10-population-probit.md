# 10. Optional population liability ICC

The optional grouped-binomial probit profile estimates the Gaussian
random-intercept **population liability ICC**. It is a separate estimand
from the classical finite-collection lineage-membership share. Its results
do not change the count of admitted membership-share cells.

Enable it explicitly in an otherwise valid configuration:

```yaml
population_probit:
  enabled: true
```

The default is disabled. To run the optional method without classical
attribution, also set `attribution: {enabled: false}`. Other analysis sections
retain their own settings. No cutoff override is accepted in configuration.
The record stores results under `metadata_diagnostics.population_probit_icc`;
both reports add a distinct population-model section.

The grouped API takes successes and tested counts per group:

```python
from amr_clonalshare import population_probit_icc

result = population_probit_icc([0, 1, 3, 4], [4, 4, 4, 4])
record = result.as_dict()
```

Each vector must be nonempty, one-dimensional and integer-valued, with
`0 <= counts <= sizes`. The implementation size cap is 1,000,000, which is
not an empirical-validation or practical-runtime guarantee. Invalid grouped
inputs raise `ValueError`. Pipeline results retain each agent's finite,
tested-and-typed subset and both exclusion counts; an absent result is never
converted into a negative outcome. Nonbinary retained calls are reported as
invalid input.

| Status | Meaning |
|---|---|
| `ok` | Numerical fit and finite profile diagnostics completed; model validity is not established by this status. |
| `constant_outcome_unidentified` | No ICC point estimate; explicitly uninformative [0,1] set. |
| `insufficient_repeated_groups` | Fewer than two groups contain repeated observations; no point estimate and [0,1]. This is a reporting policy. |
| `numerical_failure` | No point estimate or interval; the numerical reason is retained. |
| `no_retained_calls` | No finite result with a recorded lineage for this agent. |
| `invalid_binary_input` or `invalid_grouped_input` | The input violates the binary/grouped-count contract or implementation bounds. |

The minimum of two repeated groups is separate from the classical 0.90
support policy. The model assumes independent Gaussian group effects,
conditional binomial sampling, group sizes fixed independently of the effects,
and ignorable selection for the modeled population. Marginal prevalence is
jointly fitted and can differ from the raw observed proportion when group
sizes differ. Neither ICC nor its interval identifies transmission,
horizontal transfer, selection or intervention effects.

The fixed likelihood-ratio cutoff is **6.505145233310685**, from 36,000
calibration draws (seed 2026091124). The independent validation used 68,000
fresh draws (seed 2026091125), across 68 cells, with zero numerical failures.
All-draw cell coverage ranged from 96.2% to 100.0%; this includes uninformative
sets. Overcoverage can cost precision. The code, cutoff and validation-summary
hashes, all 68 tested designs and available/informative denominators are
packaged in `population_probit_validation.json`; each result identifies that
record and the underlying evidence.

The tested designs used group counts 10, 20, 30 and 100, model prevalences
0.02, 0.085, 0.25 and 0.5, and balanced or specified unequal group sizes up to
1,000. Some designs were outside calibration, but had been seen during V2
numerical development; the final V3 draws were held out. These finite-grid
results are **not a universal 95% coverage guarantee**. Range inclusion
does not validate every intermediate design. Result records flag group-count,
group-size and fitted-prevalence extrapolation and retain the assumptions.
A finite optimization/topology scan also cannot prove global correctness for
every possible dataset. Large legal groups can be expensive because the
quadrature and its cache grow with group size.
Many distinct group sizes can also be expensive: the frozen implementation
retains 32 normal-quadrature rules, so a fit with more distinct sizes may
recompute them repeatedly. The accepted numerical source is preserved byte
for byte; performance changes require their own documented verification.

## General method

An additional explicit choice, `interval_method: general`, computes for the
actual retained group sizes. See
[General population-model confidence limits](../GENERAL_POPULATION_ICC.md) for
its separate geometry branches, approximate 95% coverage target, resource
planning and resumable-computation contract. A finite-panel empirical
confirmation of the coverage target is recorded with the method; the approximation
and model assumptions remain explicit. This option does not change the default `fixed_v3` choice or turn the
finite-grid evidence above into a general coverage guarantee.
