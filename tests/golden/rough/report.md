# amr-clonalshare run report

Every number below is read from `clonal_share_result.json`, written by the same run; none is recomputed here.

- **Run:** DDB501EB
- **Issued:** 2026-10-03T05:58Z
- **Software:** amr-clonalshare 1.0.0
- **Isolates:** 25
- **Lineage:** lineage
- **Seed:** 42
- **Configuration:** sha256:32c7c068053a
- **Record digest:** sha256:0b9410c5f5fbc920175f00073c2d84b3033dd27fac63f2e14419cd2e719d41f0

**3 antimicrobials; 13 analysis outcomes.** 7 computed, 1 full ranges, 5 unavailable. Each result retains its own target and data subset.

## 1. Measurement summary

**Quantity measured.** The share of the variation in the recorded binary outcome across this collection that lineage membership accounts for, read from a lineage label and an interpreted result and scored on isolates the estimator did not see.

- **Collection:** 25 isolates in 6 lineages, typed by `lineage`
- **Antimicrobials read:** 3
- **Membership shares estimable:** 1 of 2

**Reporting conditions.** 7 computed, 1 full ranges, 5 unavailable. Each result retains its own target and data subset. An unavailable value failed a declared reporting condition of its estimator; the recorded reason distinguishes input limitations from numerical failure.

**Table 0.** Available analyses and their data subsets. A full range is a completed, uninformative result; an unavailable value is not zero.

| Agent | Analysis | N | Status | Data scope |
|---|---|---:|---|---|
| agent_a | finite_collection_prevalence | 25 | computed | recorded finite collection |
| agent_b | finite_collection_prevalence | 25 | computed | recorded finite collection |
| agent_c | finite_collection_prevalence | 25 | full_range | recorded finite collection |
| agent_a | collection_membership | 22 | computed | observed and typed subset |
| agent_a | call_latent_bounds | 22 | computed | observed and typed subset |
| agent_a | call_latent_lower_limit | 22 | computed | observed and typed subset |
| agent_b | collection_membership | 22 | unavailable | observed and typed subset |
| agent_b | call_latent_bounds | 22 | unavailable | observed and typed subset |
| agent_b | call_latent_lower_limit | 22 | unavailable | observed and typed subset |
| agent_a | lineage_evidence | 22 | computed | observed and typed subset |
| agent_b | lineage_evidence | 22 | computed | observed and typed subset |
| agent_c | collection_membership | 0 | unavailable | observed and typed subset |
| agent_c | lineage_evidence | 0 | unavailable | observed and typed subset |

**Table 0b.** Reasons and next checks for results requiring interpretation. These checks do not change the recorded status or justify changing reporting conditions.

| Agent | Analysis / status | Reason | Next check |
|---|---|---|---|
| agent_c | finite_collection_prevalence / full_range | The completed interval or bounds span the entire admissible range [0, 1] | Review missing-outcome counts; the recorded frame permits every prevalence from 0 to 1. |
| agent_b | collection_membership / unavailable | Constant retained outcome; no variation to attribute | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| agent_b | call_latent_bounds / unavailable | Constant retained outcome; no variation to attribute | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| agent_b | call_latent_lower_limit / unavailable | Constant retained outcome; no variation to attribute | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| agent_c | collection_membership / unavailable | No readable call available for this requested analysis | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| agent_c | lineage_evidence / unavailable | No readable call available for this requested analysis | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |

**Table 0a.** Prevalence in the recorded collection: bounds allow every missing outcome to be either negative or positive. These are identification bounds, not confidence intervals. Observed prevalence uses observed outcomes only.

| Agent | Recorded | Observed | Missing | Observed prevalence | Collection bounds |
|---|---:|---:|---:|---:|---|
| agent_a | 25 | 25 | 0 | 36.0 % | 0.360 to 0.360 |
| agent_b | 25 | 23 | 2 | 0.0 % | 0.000 to 0.080 |
| agent_c | 25 | 0 | 25 | not computed | 0.000 to 1.000 |

**Table 1.** Every trait, ordered by share, with the 95 % interval, the permutation p-value and the Benjamini-Yekutieli q-value across traits; selected is the Benjamini-Yekutieli selection at a false-discovery level of 0.05, which holds whatever the dependence between traits. Control is the share on shuffled lineage labels. A p-value written with ≤ is the smallest the permutations can give. The last column gives the conclusion of the results table of the form: established, lineage structure established, when the interval lies above zero and the trait is selected, or the lower confidence limit is above zero; not established; whole range, when the interval spans it; or not estimable. The conclusion says whether this collection shows lineage structure in the measurement, with the values behind it. It is not a statement about transmission or a resistance mechanism.

| Trait | Share | 95 % interval | p | q | Control | e-value | Selected | Conclusion |
|---|---:|---:|---:|---:|---:|---:|---|---|
| agent_a | 0.007 | 0.000 to 0.302 | 1.000 | 1.000 | 0.12 | 0.1 | no | not established |

**Table 1b.** Lineage structure behind each share: the lineages with at least two isolates on which the share is scored, the isolates of singleton lineages set aside, the effective number of scored lineages (inverse of the sum of squared lineage shares) and the share of isolates in lineages with at least two members (support). Few effective lineages mean the share rests on a few lineages.

| Trait | Lineages scored | Singletons set aside | Effective | Support |
|---|---:|---:|---:|---:|
| agent_a | 4 | 2 | 4.0 | 90.9 % |

> **Not evaluated on this run, and why**
>
> `agent_b`: not scored; every isolate carries the same call, so there is no variance to attribute.
>
> Decomposition of a prevalence difference into lineage composition and within-lineage rate: needs two collections; this record holds one.
>
> Reading at the recorded dilution: no dilutions were supplied; shares are read from binary calls.
>

> **Interpretation (generated from Table 1 by fixed rules)**
>
> For agent_a, no lineage effect is distinguishable from none in this collection at the chosen resolution. This does not establish independence from lineage or identify the reason for a prevalence change.
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
- **Susceptibility calls:** 3 antimicrobials on 25 isolates
- **Resampling:** 2 folds, 2 repeats, 20 bootstrap draws, 19 permutations per antimicrobial

Antimicrobials with fewer than 20 isolates of the rarer outcome: 2 of 3. Each method uses its own reporting conditions; sparse outcomes may yield wide intervals or an unavailable result.

Lineage groups in the input: 6, of which 2 hold a single isolate. Support, the share of isolates in lineages of at least two, is 90.9 %; the share is scored on those isolates and the singletons are set aside. At least two lineages hold two or more isolates.

*Figure 1 (drawn on the page).* Isolates per lineage, largest first. A lineage of one isolate, drawn in orange, cannot be predicted out of sample and is set aside; support is the share of isolates in the other lineages, 90.9 % here. The largest lineage holds 22.7 % of the isolates, which sets how much one lineage can weigh in the share.

**Table 2.** Conditions the estimator requires before any result is reported.

| Condition | Observed | Required | Verdict |
|---|---:|---:|---|
| Lineages with at least two isolates, fewest over the antimicrobials read | 4 | ≥ 2 | accepted |
| Support, lowest over the antimicrobials read | 90.9 % | reported | singletons set aside |

Each row uses only isolates with both a readable result for that agent and a recorded lineage. Support is the share of them in lineages of at least two isolates; the share is scored on those isolates and the singletons are set aside. The rarer-outcome count is a warning threshold, not a reporting condition of the estimator.

**Table 2a.** Per-agent input feasibility.

| Agent | Retained | Support | Input failures | Rarer outcome |
|---|---|---|---|---|
| agent_a | 22 | 90.9 % | none at input level | 9 (below 20) |
| agent_b | 22 | 90.9 % | constant trait | 0 (below 20) |
| agent_c | 0 | not defined | fewer than 2 tested and typed isolates; fewer than 2 lineages with at least 2 isolates | 0 (below 20) |

## 3. The lineage share of the call, trait by trait

The interval is for the share of the lineages in this collection: the lineages and their sizes are held fixed and the isolates of each lineage are drawn again from the smoothed distribution of its calls. For **1 of the 1 traits** shown the interval reaches zero, so no lineage effect is distinguishable from none for that trait.

The control column is the same estimator run on shuffled lineage labels; it should sit near zero, and a share is read against it rather than against zero.

*Figure 2 (drawn on the page).* Lineage share of the call by trait, point estimate with 95 % interval for the represented lineages, 1 largest of 1. Traits whose interval reaches zero are drawn in grey and marked ‡.

The last column asks what the call says about the ordering it was cut from: the lineage share of the latent MIC, or of any continuous value the call thresholds, the rank intraclass correlation of that value. A call does not identify that share; it bounds it. The bracket gives the smallest share any ordering consistent with the calls allows, the lower bound the calls establish, and the largest; an upper end marked ≤ is a certified bound rather than a share attained. The figure after it is the one-sided 95 % lower confidence limit of the lower bound, from splitting the isolates of every lineage at random, choosing a direction on one half and testing it on the other: a share of the latent ordering that the lineages of this collection are shown to account for, with no model for the latent values beyond their agreeing with the readings. A MIC reading on a panel that contains the cut-off refines the call, and refining the readings can only narrow the bounds they allow.

**Table 3.** The share of every trait with its interval, and the bounds the call places on the share of the latent ordering.

| Trait | Share | 95 % interval, represented lineages | Latent ordering: bounds; 95 % lower limit |
|---|---:|---:|---:|
| agent_a | 0.007 | 0.000 to 0.302 | 0.000 to 0.254; ≥ 0.000 |

## 4. Evidence that survives re-reading

A p-value is a statement about one look at the data. A surveillance panel is looked at again every year, and a p-value recomputed each time loses its error control. The e-value is evidence on a scale made for that: this run's e-value is a statement about this collection, and a programme that adds an intake each year multiplies the e-value of each new intake, scored against the lineage rates learned from the earlier ones, into a running product whose error control holds at whatever intake it is read (sequential_e_process). The e-BH procedure controls the false-discovery rate across the panel whatever the dependence between traits. Larger is stronger; 1 is no evidence.

**Table 4.** e-value per trait and the e-BH selection at level 0.05. No trait reached the e-BH threshold on this run.

| Trait | e-value | natural log | Selected |
|---|---:|---:|---|
| agent_b | 1.0 | 0.00 | no |
| agent_a | 0.1 | -2.42 | no |

0 of 2 traits are selected. Note: e-value in the betting sense of Vovk and Wang, not the BLAST expectation value and not the E-value of VanderWeele and Ding.

*Figure 3 (drawn on the page).* Evidence per trait on the natural-log scale, 2 largest of 2. The dashed rule is 1/α, the evidence one trait alone needs at level 0.05; no trait reached the e-BH threshold on this run. Traits the e-BH procedure selected are drawn in blue; the scale is logarithmic, so equal steps are equal factors of evidence.

## 5. Reading at the recorded resolution

No recorded dilutions were supplied on this run; the shares above are read from binary calls only.

## 6. A change of lineages or a change within them

A prevalence difference between two collections can be split into a change in lineage composition and a change in rate within lineages. That decomposition needs two collections and was not run here: this record holds one.

## 7. Provenance and terms

- **Software:** amr-clonalshare 1.0.0, record schema 1.0
- **Seed:** 42
- **Configuration:** `sha256:32c7c068053a`
- **Record:** `sha256:0b9410c5f5fbc920175f00073c2d84b3033dd27fac63f2e14419cd2e719d41f0`

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

amr-clonalshare 1.0.0 · record schema 1.0 · seed 42 · configuration sha256:32c7c068053a
Cite this run as: “amr-clonalshare 1.0.0, run DDB501EB, record sha256:0b9410c5f5fb.”
This report supersedes any earlier report bearing the same run identifier. It is regenerated from the record and holds no value that the record does not. The symbols carry the same wording in every run of this software; no symbol against a value means only that none of the listed conditions fired.
