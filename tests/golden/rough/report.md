# amr-clonalshare run report

Every number below is read from `clonal_share_result.json`, written by the same run, or from the release's validation grid where the text says so; none is recomputed here.

- **Run:** 53B7F822
- **Issued:** 2026-09-26T01:26Z
- **Software:** amr-clonalshare 1.0.0
- **Isolates:** 25
- **Lineage:** lineage
- **Seed:** 42
- **Configuration:** sha256:3856283ad807
- **Record digest:** sha256:23a9875e767bafa735447ac8de12eeaf9ae58cdce655efae01efb06ca7d19d02

**3 antimicrobials; 12 analysis outcomes.** 6 computed, 1 full ranges, 5 unavailable or incomplete. Each result retains its own target and data subset.

## 1. Measurement summary

**Quantity measured.** The share of the variation in the recorded binary outcome across this collection that lineage membership accounts for, read from a lineage label and an interpreted result and scored on isolates the estimator did not see.

- **Cohort:** 25 isolates in 6 lineages, typed by `lineage`
- **Antimicrobials read:** 3
- **Membership shares estimable:** 1 of 2

**Gates.** 6 computed, 1 full ranges, 5 unavailable or incomplete. Each result retains its own target and data subset. A withheld value failed a declared reporting condition of its estimator; the recorded reason distinguishes input limitations from numerical failure. Passing a gate does not verify the model.

**Table 0.** Available analyses and their data subsets. A full range is a completed, uninformative result; an unavailable value is not zero.

| Agent | Analysis | N | Status | Data scope |
|---|---|---:|---|---|
| agent_a | finite_collection_prevalence | 25 | computed | recorded finite collection |
| agent_b | finite_collection_prevalence | 25 | computed | recorded finite collection |
| agent_c | finite_collection_prevalence | 25 | full_range | recorded finite collection |
| agent_a | collection_membership | 22 | computed | observed and typed subset |
| agent_b | collection_membership | 22 | unavailable | observed and typed subset |
| agent_a | realised_component | 22 | computed | observed and typed subset |
| agent_b | realised_component | 22 | unavailable | observed and typed subset |
| agent_a | lineage_evidence | 22 | computed | observed and typed subset |
| agent_b | lineage_evidence | 22 | computed | observed and typed subset |
| agent_c | collection_membership | 0 | unavailable | observed and typed subset |
| agent_c | realised_component | 0 | unavailable | observed and typed subset |
| agent_c | lineage_evidence | 0 | unavailable | observed and typed subset |

**Table 0b.** Reasons and next checks for results requiring interpretation. These checks do not change the recorded status or justify changing reporting gates.

| Agent | Analysis / status | Reason | Next check |
|---|---|---|---|
| agent_c | finite_collection_prevalence / full_range | The completed interval or bounds span the entire admissible range [0, 1] | Review missing-outcome counts; the recorded frame permits every prevalence from 0 to 1. |
| agent_b | collection_membership / unavailable | Constant retained outcome; no variation to attribute | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| agent_b | realised_component / unavailable | every lineage is constant, so the within-lineage scale is zero and no share is identified | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| agent_c | collection_membership / unavailable | No readable call available for this requested analysis | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| agent_c | realised_component / unavailable | No readable call available for this requested analysis | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| agent_c | lineage_evidence / unavailable | No readable call available for this requested analysis | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |

**Table 0a.** Prevalence in the recorded collection: bounds allow every missing outcome to be either negative or positive. These are identification bounds, not confidence intervals. Observed prevalence uses observed outcomes only.

| Agent | Recorded | Observed | Missing | Observed prevalence | Collection bounds |
|---|---:|---:|---:|---:|---|
| agent_a | 25 | 25 | 0 | 36.0 % | 0.360 to 0.360 |
| agent_b | 25 | 23 | 2 | 0.0 % | 0.000 to 0.080 |
| agent_c | 25 | 0 | 25 | not computed | 0.000 to 1.000 |

**Table 1.** Every trait, ordered by share. The reading in the last column is fixed by two conditions the record holds: whether the interval excludes zero and whether the Benjamini-Yekutieli step-up on the permutation p-values selected the trait at a false-discovery level of 0.05, which holds whatever the dependence between traits.

| Trait | Share | 95 % interval | Control | e-value |  | Reading |
|---|---:|---:|---:|---:|---|---|
| agent_a | 0.038 | -0.624 to 0.258 | -0.47 | 0.1 | ‡ | no detectable lineage effect |

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
| Lineages with at least two isolates, fewest over the antimicrobials read | 6 | ≥ 2 | accepted |
| Support, lowest over the antimicrobials read | 90.9 % | reported | singletons set aside |

Each row uses only isolates with both a readable result for that agent and a recorded lineage. Support is the share of them in lineages of at least two isolates; the share is scored on those isolates and the singletons are set aside. The rarer-outcome count is a warning threshold, not an estimator gate.

**Table 2a.** Per-agent input feasibility.

| Agent | Retained | Support | Input failures | Rarer outcome |
|---|---|---|---|---|
| agent_a | 22 | 90.9 % | none at input level | 9 (below 20) |
| agent_b | 22 | 90.9 % | constant trait | 0 (below 20) |
| agent_c | 0 | not defined | fewer than 2 tested and typed isolates; fewer than 2 lineages with at least 2 isolates | 0 (below 20) |

## 3. The clonal share, trait by trait

The intervals have a nominal 95 % level; their measured coverage depends on the cohort design and model assumptions. For **1 of the 1 traits** shown the interval includes zero, so no lineage effect is distinguishable from none for that trait.

The control column is the same estimator run on shuffled lineage labels; it should sit near zero, and a share is read against it rather than against zero.

*Figure 2 (drawn on the page).* Clonal share by trait, point estimate with 95 % interval, 1 largest of 1. Traits whose interval crosses zero are drawn in grey and marked ‡. The thin line beneath each interval is the species interval of Table 3, for a fresh draw of lineages; it is floored at zero by construction. The intervals are drawn as computed. The quantity lies between 0 and 1, so the part of an interval below zero carries no information: reading each interval as its overlap with that range leaves the coverage unchanged, and a lower limit at or below zero means the same thing either way.

On the release's validation grid, the interval this estimator prints for a binary trait contained the truth in 0.970 of 28,000 runs over 70 simulated cohorts (by cohort, 0.825 to 1.000). That figure belongs to the release, not to this run; it is what the phrase "95 % interval" was measured to mean.

Two intervals answer two questions. The interval above is for the share the lineages in this collection carry. The second interval below is for the share a fresh draw of lineages from the species would show, stated on the scale of the realised share of Table 3b below; it adds the sampling of the lineages themselves, on 5 degrees of freedom, is widened to the envelope of the first taken on that scale (with a lower end no smaller than zero, since a species share is not negative), and is the one to quote when the figure is read as a property of the species rather than of this collection. With few lineages it is markedly wider; with many the two nearly coincide.

This collection holds fewer than ten lineages. On the release's validation grid the interval for lineage membership contained the truth in only 0.82 to 0.96 of runs at 5 lineages. A resample of so few lineages has too few distinct outcomes, and the species interval, wider in the normal-effect cells, answers a different question and does not repair uncertainty about this collection.

The last column is an exploratory, model-equivalent restatement under a Gaussian probit threshold model. It transforms the collection's observed-scale share and interval at the estimated prevalence, treated as fixed. It does not estimate the actual latent variance of these particular lineages or a generic binomial mixed-model intraclass correlation. Monotonicity preserves coverage only for a matching target at a known, fixed prevalence; it does not guarantee coverage after estimating prevalence. Read the observed-scale interval as primary. Values at or below zero are displayed at the zero boundary.

**Table 3.** The two intervals for every trait: for lineage membership in this collection and for the species, with the share on the latent scale.

| Trait | Share | 95 % interval, these lineages | Species share (95 % interval) | Latent restatement (exploratory interval) |
|---|---:|---:|---:|---:|
| agent_a | 0.038 | -0.624 to 0.258 | 0.051 (0.000 to 1.000) | 0.062 (0.000 to 0.399) |

**The realised share of the same call.** The variance-component ratio of the lineages in hand, with an interval that is exact under a Gaussian within-lineage law. Its gate reads the excess kurtosis of the within-lineage residuals and withholds the interval above 0.99. Rare binary calls often fail this check; prevalence alone does not determine the residual kurtosis. The interval is approximate on binary data even when the gate opens; the estimate is printed either way, the interval only where the gate opened, for 1 of 2 traits here.

**Table 3b.** The realised share per trait and the verdict of its kurtosis gate.

| Trait | Realised share | Gaussian-model 95 % interval | Residual excess kurtosis | Gate |
|---|---:|---:|---:|---|
| agent_a | 0.000 | 0.000 to 0.200 | -2.02 | open |
| agent_b | not computed | withheld | not computed | closed |

## 4. Evidence that survives re-reading

A p-value is a statement about one look at the data. A surveillance panel is looked at again every year, and a p-value recomputed each time loses its guarantee. The e-value is evidence on a scale made for that: this run's e-value is a statement about this cohort, and a programme that adds an intake each year multiplies the e-value of each new intake, scored against the lineage rates learned from the earlier ones, into a running product whose guarantee holds at whatever intake it is read (sequential_e_process). The e-BH procedure controls the false-discovery rate across the panel whatever the dependence between traits. Larger is stronger; 1 is no evidence.

**Table 4.** e-value per trait and the e-BH selection at level 0.05. No trait reached the e-BH threshold on this run.

| Trait | e-value | natural log | Selected |
|---|---:|---:|---|
| agent_b | 1.0 | -0.00 | no |
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
- **Configuration:** `sha256:3856283ad807`
- **Record:** `sha256:23a9875e767bafa735447ac8de12eeaf9ae58cdce655efae01efb06ca7d19d02`

> **Terms used in this report**
>
> **Share.** How much of the variation of the call between isolates is associated with their recorded lineage. Near 1: a strong association on the measured scale. Near 0: little association on that scale. This does not identify transmission or its mechanism.
>
> **Collection-bootstrap interval.** The lineage-membership interval targets the sampled collection. Separate Gaussian-model component and species intervals retain their own targets and assumptions; coverage measured on a finite simulation grid is not a universal guarantee.
>
> **Control.** The same calculation on shuffled lineage labels, which should return roughly zero; a share is read against it.
>
> **Support.** The share of isolates in lineages with at least two members. A singleton lineage cannot be predicted out of sample, so the share is scored on the other isolates and the singletons are set aside; the share then speaks for singleton lineages only if they are drawn from the same population of lineages as the repeated ones.
>
> **Gate.** A condition fixed before the run. If it does not hold, the software withholds the estimate and names the condition that failed.
>
> **e-value.** A fixed-look e-value and a sequential e-process have different uses. Repeated-look validity requires the sequential construction and its stated null-model assumptions; a fixed-look value alone does not justify repeated inspection. Larger values are stronger evidence against the specified null.
>

## Notes on reading this report

- ‡ The 95 % interval includes zero. No lineage effect is distinguishable from none for this trait, and the point estimate must not be read on its own.
- § Derived arithmetically from a quantity recorded elsewhere in the record rather than estimated in this run.

amr-clonalshare 1.0.0 · record schema 1.0 · seed 42 · configuration sha256:3856283ad807
Cite this run as: “amr-clonalshare 1.0.0, run 53B7F822, record sha256:23a9875e767b.”
This report supersedes any earlier report bearing the same run identifier. It is regenerated from the record and holds no value that the record does not. The symbols carry the same wording in every run of this software; no symbol against a value means only that none of the listed conditions fired.
Classical collection-bootstrap endpoints are printed without clipping; the separately labeled species interval is floored at zero as part of its construction.
