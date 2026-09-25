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

The minimum of two repeated groups is separate from the classical 0.80
support policy. The model assumes independent Gaussian group effects,
conditional binomial sampling, group sizes fixed independently of the effects,
and ignorable selection for the modeled population. Marginal prevalence is
jointly fitted and can differ from the raw observed proportion when group
sizes differ. Neither ICC nor its interval identifies transmission,
horizontal transfer, selection or intervention effects.

The fixed likelihood-ratio cut-off is **5.443134202054537**. It is the largest,
over 216 calibration designs with 2,000 datasets each (432,000 datasets, seed
20260926101), of the 95th percentile of the statistic at the generating value;
the rule, the designs and every seed were fixed before any result was computed
(`benchmarks/population_model/PROTOCOL.md`). The validation used fresh
datasets from 144 designs off the calibration grid (15, 50 or 200 lineages of
equal or skewed sizes, and the 30 lineage sizes of the *S. suis* example;
prevalence 0.02 to 0.85), 2,000 datasets each (288,000 datasets, seed
20260926202), with no numerical failure. Coverage ranged from 95.9% to 100.0%
and every design met the prespecified rule, an upper Wilson 95% limit of at
least 0.95; an uninformative [0,1] set counts as covering. The median width
of the interval was 0.45. The code, cut-off and validation-summary hashes, all
144 designs and their denominators are packaged in
`population_probit_validation.json`; each result names that record.

These finite-grid results are **not a universal 95% coverage guarantee**, and
range inclusion does not validate every intermediate design. Result records
flag group counts, group sizes and fitted prevalences outside the tested
designs, and retain the assumptions. When an assumption fails, coverage can
fall far below the nominal level. With lineage effects from a t distribution
on four degrees of freedom it stayed between 90.8% and 95.1%; with lineage
sizes that grow with the lineage effect it fell to between 61.9% and 85.8%;
and with a two-point (carrier) law of lineage effects it fell to between 1.4%
and 86.9%, lowest at a latent ICC of 0.7. Large legal groups can be expensive
because the quadrature and its cache grow with group size, and a fit with more
than 32 distinct group sizes recomputes normal-quadrature rules repeatedly.

Under the model the latent ICC does not depend on prevalence, whereas the
observed-scale lineage share does. In 500 datasets of 30 lineages of 20
isolates with a latent ICC of 0.3, the mean profile estimate stayed between
0.28 and 0.29 at prevalences from 0.05 to 0.50, while the mean observed-scale
share rose from 0.05 to 0.19; at a latent ICC of 0.6 the estimate stayed
between 0.55 and 0.58 and the share rose from 0.15 to 0.40. At a prevalence
of 0.02 the estimate fell to 0.26 and 0.51. The two quantities answer
different questions and are reported side by side.

## General method

An additional explicit choice, `interval_method: general`, computes for the
actual retained group sizes. See
[General population-model confidence limits](../GENERAL_POPULATION_ICC.md) for
its separate geometry branches, approximate 95% coverage target, resource
planning and resumable-computation contract. A finite-panel empirical
validation of the coverage target is recorded with the method; the approximation
and model assumptions remain explicit. This option does not change the default
`fixed_cutoff` choice or turn the finite-grid evidence above into a general
coverage guarantee.
