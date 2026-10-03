# amr-clonalshare run report

Every number below is read from `clonal_share_result.json`, written by the same run; none is recomputed here.

- **Run:** C61ADB50
- **Issued:** 2026-10-03T05:58Z
- **Software:** amr-clonalshare 1.0.0
- **Isolates:** 96
- **Lineage:** lineage
- **Seed:** 42
- **Configuration:** sha256:e55bce1fa118
- **Record digest:** sha256:30dce49612c9a27fd11f59b9b7ed5343ac80b1535b5ba2ecd6dcfed1dc9412b3

**1 antimicrobials; 7 analysis outcomes.** 7 computed, 0 full ranges, 0 unavailable. Each result retains its own target and data subset.

## 1. Measurement summary

**Quantity measured.** The share of the variation in the recorded binary outcome across this collection that lineage membership accounts for, read from a lineage label and an interpreted result and scored on isolates the estimator did not see.

- **Collection:** 96 isolates in 12 lineages, typed by `lineage`
- **Antimicrobials read:** 1
- **Membership shares estimable:** 1 of 1

**Reporting conditions.** 7 computed, 0 full ranges, 0 unavailable. Each result retains its own target and data subset. An unavailable value failed a declared reporting condition of its estimator; the recorded reason distinguishes input limitations from numerical failure.

**Table 0.** Available analyses and their data subsets. A full range is a completed, uninformative result; an unavailable value is not zero.

| Agent | Analysis | N | Status | Data scope |
|---|---|---:|---|---|
| demo_agent | finite_collection_prevalence | 96 | computed | recorded finite collection |
| demo_agent | collection_membership | 96 | computed | observed and typed subset |
| demo_agent | call_latent_bounds | 96 | computed | observed and typed subset |
| demo_agent | call_latent_lower_limit | 96 | computed | observed and typed subset |
| demo_agent | lineage_evidence | 96 | computed | observed and typed subset |
| demo_agent | decomposition_composition | 96 | computed | observed and typed subset |
| demo_agent | decomposition_within_lineage | 96 | computed | observed and typed subset |

**Table 0a.** Prevalence in the recorded collection: bounds allow every missing outcome to be either negative or positive. These are identification bounds, not confidence intervals. Observed prevalence uses observed outcomes only.

| Agent | Recorded | Observed | Missing | Observed prevalence | Collection bounds |
|---|---:|---:|---:|---:|---|
| demo_agent | 96 | 96 | 0 | 40.6 % | 0.406 to 0.406 |

**Comparison scope.** Decomposition components describe the observed, labelled subsets. Equal typing coverage or a nonsignificant missingness test does not establish representativeness. Full observed prevalence and subset prevalence are recorded separately; missing data prevent automatic collection-wide generalization.

**Table 1.** Every trait, ordered by share, with the 95 % interval, the permutation p-value and the Benjamini-Yekutieli q-value across traits; selected is the Benjamini-Yekutieli selection at a false-discovery level of 0.05, which holds whatever the dependence between traits. Control is the share on shuffled lineage labels. A p-value written with ≤ is the smallest the permutations can give. The last column gives the conclusion of the results table of the form: established, lineage structure established, when the interval lies above zero and the trait is selected, or the lower confidence limit is above zero; not established; whole range, when the interval spans it; or not estimable. The conclusion says whether this collection shows lineage structure in the measurement, with the values behind it. It is not a statement about transmission or a resistance mechanism.

| Trait | Share | 95 % interval | p | q | Control | e-value | Selected | Conclusion |
|---|---:|---:|---:|---:|---:|---:|---|---|
| demo_agent | -0.050 | 0.000 to 0.405 | 1.000 | 1.000 | 0.04 | 0.6 | no | not established |

**Table 1b.** Lineage structure behind each share: the lineages with at least two isolates on which the share is scored, the isolates of singleton lineages set aside, the effective number of scored lineages (inverse of the sum of squared lineage shares) and the share of isolates in lineages with at least two members (support). Few effective lineages mean the share rests on a few lineages.

| Trait | Lineages scored | Singletons set aside | Effective | Support |
|---|---:|---:|---:|---:|
| demo_agent | 12 | 0 | 12.0 | 100.0 % |

> **Not evaluated on this run, and why**
>
> Reading at the recorded dilution: no dilutions were supplied; shares are read from binary calls.
>

> **Interpretation (generated from Table 1 by fixed rules)**
>
> For demo_agent, no lineage effect is distinguishable from none in this collection at the chosen resolution. This does not establish independence from lineage or identify the reason for a prevalence change.
>
> These readings describe association in the sampled collection. The analysis does not identify transmission, horizontal transfer, selection, or the effect of an intervention.
>

## 2. Admissibility of the input

- **Lineage column:** `lineage`
- **Phenotype interpretation:** clinical_sir
- **Positive outcome:** R only in fictional format example
- **Applied coding:** R (resistant)
- **Intermediate policy:** susceptible
- **Interpretation source:** fictional teaching fixture with no AST interpretation
- **AST standard/version:** not supplied / not supplied
- **Susceptibility calls:** 1 antimicrobials on 96 isolates
- **Resampling:** 2 folds, 2 repeats, 20 bootstrap draws, 19 permutations per antimicrobial

Antimicrobials with fewer than 20 isolates of the rarer outcome: 0 of 1. Each method uses its own reporting conditions; sparse outcomes may yield wide intervals or an unavailable result.

Lineage groups in the input: 12, of which 0 hold a single isolate. Support, the share of isolates in lineages of at least two, is 100.0 %; the share is scored on those isolates and the singletons are set aside. At least two lineages hold two or more isolates.

*Figure 1 (drawn on the page).* Isolates per lineage, largest first. A lineage of one isolate, drawn in orange, cannot be predicted out of sample and is set aside; support is the share of isolates in the other lineages, 100.0 % here. The largest lineage holds 8.3 % of the isolates, which sets how much one lineage can weigh in the share.

**Table 2.** Conditions the estimator requires before any result is reported.

| Condition | Observed | Required | Verdict |
|---|---:|---:|---|
| Lineages with at least two isolates, fewest over the antimicrobials read | 12 | ≥ 2 | accepted |
| Support, lowest over the antimicrobials read | 100.0 % | reported | singletons set aside |

Each row uses only isolates with both a readable result for that agent and a recorded lineage. Support is the share of them in lineages of at least two isolates; the share is scored on those isolates and the singletons are set aside. The rarer-outcome count is a warning threshold, not a reporting condition of the estimator.

**Table 2a.** Per-agent input feasibility.

| Agent | Retained | Support | Input failures | Rarer outcome |
|---|---|---|---|---|
| demo_agent | 96 | 100.0 % | none at input level | 39 |

## 3. The lineage share of the call, trait by trait

The interval is for the share of the lineages in this collection: the lineages and their sizes are held fixed and the isolates of each lineage are drawn again from the smoothed distribution of its calls. For **1 of the 1 traits** shown the interval reaches zero, so no lineage effect is distinguishable from none for that trait.

The control column is the same estimator run on shuffled lineage labels; it should sit near zero, and a share is read against it rather than against zero.

*Figure 2 (drawn on the page).* Lineage share of the call by trait, point estimate with 95 % interval for the represented lineages, 1 largest of 1. Traits whose interval reaches zero are drawn in grey and marked ‡.

The last column asks what the call says about the ordering it was cut from: the lineage share of the latent MIC, or of any continuous value the call thresholds, the rank intraclass correlation of that value. A call does not identify that share; it bounds it. The bracket gives the smallest share any ordering consistent with the calls allows, the lower bound the calls establish, and the largest; an upper end marked ≤ is a certified bound rather than a share attained. The figure after it is the one-sided 95 % lower confidence limit of the lower bound, from splitting the isolates of every lineage at random, choosing a direction on one half and testing it on the other: a share of the latent ordering that the lineages of this collection are shown to account for, with no model for the latent values beyond their agreeing with the readings. A MIC reading on a panel that contains the cut-off refines the call, and refining the readings can only narrow the bounds they allow.

**Table 3.** The share of every trait with its interval, and the bounds the call places on the share of the latent ordering.

| Trait | Share | 95 % interval, represented lineages | Latent ordering: bounds; 95 % lower limit |
|---|---:|---:|---:|
| demo_agent | -0.050 | 0.000 to 0.405 | 0.000 to 0.430; ≥ 0.000 |

## 4. Evidence that survives re-reading

A p-value is a statement about one look at the data. A surveillance panel is looked at again every year, and a p-value recomputed each time loses its error control. The e-value is evidence on a scale made for that: this run's e-value is a statement about this collection, and a programme that adds an intake each year multiplies the e-value of each new intake, scored against the lineage rates learned from the earlier ones, into a running product whose error control holds at whatever intake it is read (sequential_e_process). The e-BH procedure controls the false-discovery rate across the panel whatever the dependence between traits. Larger is stronger; 1 is no evidence.

**Table 4.** e-value per trait and the e-BH selection at level 0.05. No trait reached the e-BH threshold on this run.

| Trait | e-value | natural log | Selected |
|---|---:|---:|---|
| demo_agent | 0.6 | -0.43 | no |

0 of 1 traits are selected. Note: e-value in the betting sense of Vovk and Wang, not the BLAST expectation value and not the E-value of VanderWeele and Ding.

*Figure 3 (drawn on the page).* Evidence per trait on the natural-log scale, 1 largest of 1. The dashed rule is 1/α, the evidence one trait alone needs at level 0.05; no trait reached the e-BH threshold on this run. Traits the e-BH procedure selected are drawn in blue; the scale is logarithmic, so equal steps are equal factors of evidence.

## 5. Reading at the recorded resolution

No recorded dilutions were supplied on this run; the shares above are read from binary calls only.

## 6. A change of lineages or a change within them

Contrast period: earlier minus later: of 1 traits, 0 have a composition component (a change in lineage mix) and 0 a within-lineage component (a change in rate) selected by the Benjamini-Yekutieli step-up within each component family; the intervals and the step-up rest on nominal percentile-bootstrap tail probabilities (their calibration is reported with the package).

With at least 20 finite draw(s), the best attainable BY-adjusted value over 1 agents is 0.095, above the target 0.05; increase surveillance.n_boot before interpreting an absence of discoveries.

**Table 6.** The prevalence difference of each trait, earlier minus later, split into a change in lineage composition and a change in rate within lineages. The two components sum to the difference. The limits are bootstrap percentiles; p is the bootstrap p-value of each component and q its Benjamini-Yekutieli q-value within the component family; the last column names the components the false-discovery procedure selected.

| Trait | Difference | Composition | 95 % interval | p | q | Within lineage | 95 % interval | p | q | Selected |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| demo_agent | -0.021 | 0.000 | -0.033 to 0.068 | 0.667 | 0.667 | -0.021 | -0.236 to 0.177 | 1.000 | 1.000 | neither |

**Table 6b.** The two collections behind each difference: the prevalence and the isolates of each, the lineages both hold and those seen in one only, the share of isolates in lineages both hold (shared support, the mean of the two collections; the within-lineage component needs at least 0.90) and the turnover share, the part of the difference carried by lineages seen in one collection only.

| Trait | Prevalence, earlier | n | Prevalence, later | n | Lineages shared; only first; only second | Shared support | Turnover share |
|---|---:|---:|---:|---:|---:|---:|---:|
| demo_agent | 39.6 % | 48 | 41.7 % | 48 | 12; 0; 0 | 100.0 % | 0.000 |

## 7. Provenance and terms

- **Software:** amr-clonalshare 1.0.0, record schema 1.0
- **Seed:** 42
- **Configuration:** `sha256:e55bce1fa118`
- **Record:** `sha256:30dce49612c9a27fd11f59b9b7ed5343ac80b1535b5ba2ecd6dcfed1dc9412b3`

> **Terms used in this report**
>
> **Share.** How much of the variation of the call between isolates is associated with their recorded lineage. Near 1: a strong association on the measured scale. Near 0: little association on that scale. This does not identify transmission or its mechanism.
>
> **Bounds on the latent ordering.** The smallest and the largest lineage share of the ordering of the latent values (the MIC behind a reading or a call) that the readings allow; the smallest is the lower bound the readings establish. They are identification bounds, not a confidence interval; the one-sided lower confidence limit beside them adds the sampling of isolates.
>
> **Interval for the represented lineages.** The 95 % interval of a share holds the lineages of this collection and their sizes fixed and draws the isolates of every lineage again from the smoothed distribution of its readings, scoring every draw as the data were scored and studentizing it by its own standard error (the bootstrap-t interval). It does not describe a fresh draw of lineages from the species.
>
> **Control.** The same calculation on shuffled lineage labels, which should return roughly zero; a share is read against it.
>
> **Support.** The share of isolates in lineages with at least two members. A singleton lineage cannot be predicted out of sample, so the share is scored on the other isolates and the singletons are set aside; the share then speaks for singleton lineages only if they are drawn from the same population of lineages as the repeated ones.
>
> **Reporting condition.** A condition fixed before the run. If it does not hold, the software withholds the estimate and names the condition that failed.
>
> **e-value.** A fixed-look e-value and the running evidence over intakes have different uses. Repeated-look validity requires the sequential construction and its stated null-model assumptions; a fixed-look value alone does not justify repeated inspection. Larger values are stronger evidence against the specified null.
>

## Notes on reading this report

- ‡ The 95 % interval includes zero. No lineage effect is distinguishable from none for this trait, and the point estimate must not be read on its own.
- § Derived arithmetically from a quantity recorded elsewhere in the record rather than estimated in this run.

amr-clonalshare 1.0.0 · record schema 1.0 · seed 42 · configuration sha256:e55bce1fa118
Cite this run as: “amr-clonalshare 1.0.0, run C61ADB50, record sha256:30dce49612c9.”
This report supersedes any earlier report bearing the same run identifier. It is regenerated from the record and holds no value that the record does not. The symbols carry the same wording in every run of this software; no symbol against a value means only that none of the listed conditions fired.
