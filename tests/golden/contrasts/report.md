# amr-clonalshare run report

Every number below is read from `clonal_share_result.json`, written by the same run, or from the release's validation grid where the text says so; none is recomputed here.

- **Run:** E1E7B29E
- **Issued:** 2026-09-25T15:42Z
- **Software:** amr-clonalshare 1.0.0
- **Isolates:** 96
- **Lineage:** lineage
- **Seed:** 42
- **Configuration:** sha256:f12b49aad646
- **Record digest:** sha256:9c2eb99ec5e20b8e4bc9b9d6bb87b5de54f83e15556afbbebc2e3ad0e4cbf0a9

**1 antimicrobials; 6 analysis outcomes.** 5 computed, 0 full ranges, 1 unavailable or incomplete. Each result retains its own target and data subset.

## 1. Measurement summary

**Quantity measured.** The share of the variation in the recorded binary outcome across this collection that lineage membership accounts for, read from a lineage label and an interpreted result and scored on isolates the estimator did not see.

- **Cohort:** 96 isolates in 12 lineages, typed by `lineage`
- **Antimicrobials read:** 1
- **Membership shares estimable:** 1 of 1

**Gates.** 5 computed, 0 full ranges, 1 unavailable or incomplete. Each result retains its own target and data subset. A withheld value failed a declared reporting condition of its estimator; the recorded reason distinguishes input limitations from numerical failure. Passing a gate does not verify the model.

**Table 0.** Available analyses and their data subsets. A full range is a completed, uninformative result; an unavailable value is not zero.

| Agent | Analysis | N | Status | Data scope |
|---|---|---:|---|---|
| demo_agent | finite_collection_prevalence | 96 | computed | recorded finite collection |
| demo_agent | collection_membership | 96 | computed | observed and typed subset |
| demo_agent | realised_component | 96 | unavailable | observed and typed subset |
| demo_agent | lineage_evidence | 96 | computed | observed and typed subset |
| demo_agent | decomposition_composition | 96 | computed | observed and typed subset |
| demo_agent | decomposition_within_lineage | 96 | computed | observed and typed subset |

**Table 0b.** Reasons and next checks for results requiring interpretation. These checks do not change the recorded status or justify changing reporting gates.

| Agent | Analysis / status | Reason | Next check |
|---|---|---|---|
| demo_agent | realised_component / unavailable | the observed mean-square ratio 0.267 falls below the 0.025 point of the central F on 11 and 84 degrees of freedom, so the inverted confidence set for the noncentrality is empty and the interval closes on zero | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |

**Table 0a.** Prevalence in the recorded collection: bounds allow every missing outcome to be either negative or positive. These are identification bounds, not confidence intervals. Observed prevalence uses observed outcomes only.

| Agent | Recorded | Observed | Missing | Observed prevalence | Collection bounds |
|---|---:|---:|---:|---:|---|
| demo_agent | 96 | 96 | 0 | 40.6 % | 0.406 to 0.406 |

**Comparison scope.** Decomposition components describe the observed, labelled subsets. Equal typing coverage or a nonsignificant missingness test does not establish representativeness. Full observed prevalence and subset prevalence are recorded separately; missing data prevent automatic collection-wide generalization.

**Table 1.** Every trait, ordered by share. The reading in the last column is fixed by two conditions the record holds: whether the interval excludes zero and whether the e-BH procedure selected the trait at level 0.05.

| Trait | Share | 95 % interval | Control | e-value |  | Reading |
|---|---:|---:|---:|---:|---|---|
| demo_agent | -0.034 | -0.339 to 0.047 | -0.27 | 0.6 | ‡ | no detectable lineage effect |

**Table 1b.** Lineage structure behind each share: lineages with at least two isolates, the effective number of lineages (inverse of the sum of squared lineage shares) and the share of isolates in lineages with at least two members (support). A high support with few effective lineages means the share rests on a few lineages.

| Trait | Lineages | Repeated | Effective | Support |
|---|---:|---:|---:|---:|
| demo_agent | 12 | 12 | 12.0 | 100.0 % |

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

Lineage groups in the input: 12, of which 0 hold a single isolate. Support 100.0 % against the 80.0 % the estimator needs: accepted.

*Figure 1 (drawn on the page).* Isolates per lineage, largest first. A lineage of one isolate, drawn in orange, cannot be predicted out of sample and counts against support; support is the share of isolates in the other lineages, 100.0 % here against the 80.0 % the estimator requires. The largest lineage holds 8.3 % of the isolates, which sets how much one lineage can weigh in the share.

**Table 2.** Conditions the estimator requires before any result is reported.

| Condition | Observed | Required | Verdict |
|---|---:|---:|---|
| Lineage support, lowest over the antimicrobials read | 100.0 % | ≥ 80.0 % | accepted |
| Lineage groups used | 12 | reported | accepted |

Each row uses only isolates with both a readable result for that agent and a recorded lineage. The last column is an algebraic support scenario: one genuinely new tested isolate in each of that many distinct singleton lineages. It does not guarantee interval precision, coverage, or overall estimability; it cannot repair a constant trait or a missing lineage contrast. Duplicating existing rows adds no evidence. The rarer-outcome count is a warning threshold, not an additional estimator gate.

**Table 2a.** Per-agent input feasibility and support scenario.

| Agent | Retained | Support | Input failures | Rarer outcome | New isolates for support |
|---|---|---|---|---|---|
| demo_agent | 96 | 100.0 % | none at input level | 39 | 0 |

## 3. The clonal share, trait by trait

The intervals have a nominal 95 % level; their measured coverage depends on the cohort design and model assumptions. For **1 of the 1 traits** shown the interval includes zero, so no lineage effect is distinguishable from none for that trait.

The control column is the same estimator run on shuffled lineage labels; it should sit near zero, and a share is read against it rather than against zero.

*Figure 2 (drawn on the page).* Clonal share by trait, point estimate with 95 % interval, 1 largest of 1. Traits whose interval crosses zero are drawn in grey and marked ‡. The thin line beneath each interval is the species interval of Table 3, for a fresh draw of lineages; it is floored at zero by construction. The intervals are drawn as computed. The quantity lies between 0 and 1, so the part of an interval below zero carries no information: reading each interval as its overlap with that range leaves the coverage unchanged, and a lower limit at or below zero means the same thing either way.

On the release's validation grid, the interval this estimator prints for a binary trait contained the truth in 0.970 of 28,000 runs over 70 simulated cohorts (by cohort, 0.825 to 1.000). That figure belongs to the release, not to this run; it is what the phrase "95 % interval" was measured to mean.

Two intervals answer two questions. The interval above is for the share the lineages in this collection carry. The second interval below is for the share a fresh draw of lineages from the species would show, stated on the scale of the realised share of Table 3b below; it adds the sampling of the lineages themselves, on 11 degrees of freedom, is widened to the envelope of the first taken on that scale (with a lower end no smaller than zero, since a species share is not negative), and is the one to quote when the figure is read as a property of the species rather than of this collection. With few lineages it is markedly wider; with many the two nearly coincide.

The last column is an exploratory, model-equivalent restatement under a Gaussian probit threshold model. It transforms the collection's observed-scale share and interval at the estimated prevalence, treated as fixed. It does not estimate the actual latent variance of these particular lineages or a generic binomial mixed-model intraclass correlation. Monotonicity preserves coverage only for a matching target at a known, fixed prevalence; it does not guarantee coverage after estimating prevalence. Read the observed-scale interval as primary. Values at or below zero are displayed at the zero boundary.

**Table 3.** The two intervals for every trait: for lineage membership in this collection and for the species, with the share on the latent scale.

| Trait | Share | 95 % interval, these lineages | Species share (95 % interval) | Latent restatement (exploratory interval) |
|---|---:|---:|---:|---:|
| demo_agent | -0.034 | -0.339 to 0.047 | -0.037 (0.000 to 0.269) | 0.000 (0.000 to 0.075) |

**The realised share of the same call.** The variance-component ratio of the lineages in hand, with an interval that is exact under a Gaussian within-lineage law. Its gate reads the excess kurtosis of the within-lineage residuals and withholds the interval above 0.99. Rare binary calls often fail this check; prevalence alone does not determine the residual kurtosis. The interval is approximate on binary data even when the gate opens; the estimate is printed either way, the interval only where the gate opened, which it did for no trait here.

**Table 3b.** The realised share per trait and the verdict of its kurtosis gate.

| Trait | Realised share | Gaussian-model 95 % interval | Residual excess kurtosis | Gate |
|---|---:|---:|---:|---|
| demo_agent | 0.000 | withheld | -1.76 | closed |

## 4. Evidence that survives re-reading

A p-value is a statement about one look at the data. A surveillance panel is looked at again every year, and a p-value recomputed each time loses its guarantee. The e-value is evidence on a scale made for that: this run's e-value is a statement about this cohort, and a programme that adds an intake each year multiplies the e-value of each new intake, scored against the lineage rates learned from the earlier ones, into a running product whose guarantee holds at whatever intake it is read (sequential_e_process). The e-BH procedure controls the false-discovery rate across the panel whatever the dependence between traits. Larger is stronger; 1 is no evidence.

**Table 4.** e-value per trait and the e-BH selection at level 0.05. No trait reached the e-BH threshold on this run.

| Trait | e-value | natural log | Selected |
|---|---:|---:|---|
| demo_agent | 0.6 | -0.43 | no |

0 of 1 traits are selected. Note: e-value in the betting sense of Vovk and Wang, not the BLAST expectation value and not the E-value of VanderWeele and Ding.

*Figure 3 (drawn on the page).* Evidence per trait on the natural-log scale, 1 largest of 1. The dashed rule is 1/α, the evidence one trait alone needs at level 0.05; no trait reached the e-BH threshold on this run. Traits the e-BH procedure selected are drawn in blue; the scale is logarithmic, so equal steps are equal factors of evidence.

## 5. Reading at the recorded resolution

No recorded dilutions were supplied on this run; the shares above are read from binary calls only.

## 6. A change of lineages or a change within them

Contrast period: earlier vs later: of 1 traits, 0 have a detected composition component (a change in lineage mix), and 0 have a detected within-lineage component (a change in rate), after false-discovery control within each component family.

With at least 20 finite draw(s), the best attainable BY-adjusted value over 1 agents is 0.095, above the target 0.05; increase surveillance.n_boot before interpreting an absence of discoveries.

**Table 6.** The prevalence difference of each trait split into a change in lineage composition and a change in rate within lineages. The two components sum to the difference. The limits are bootstrap percentiles and are not adjusted for multiplicity; the last column names the components the false-discovery procedure selected.

| Trait | Difference | Composition | 95 % interval | Within lineage | 95 % interval | Selected |
|---|---:|---:|---:|---:|---:|---|
| demo_agent | -0.021 | 0.000 | -0.033 to 0.068 | -0.021 | -0.236 to 0.177 | neither |

## 7. Provenance and terms

- **Software:** amr-clonalshare 1.0.0, record schema 1.0
- **Seed:** 42
- **Configuration:** `sha256:f12b49aad646`
- **Record:** `sha256:9c2eb99ec5e20b8e4bc9b9d6bb87b5de54f83e15556afbbebc2e3ad0e4cbf0a9`

> **Terms used in this report**
>
> **Share.** How much of the variation of the call between isolates is associated with their recorded lineage. Near 1: a strong association on the measured scale. Near 0: little association on that scale. This does not identify transmission or its mechanism.
>
> **Collection-bootstrap interval.** The lineage-membership interval targets the sampled collection. Separate Gaussian-model component and species intervals retain their own targets and assumptions; coverage measured on a finite simulation grid is not a universal guarantee.
>
> **Control.** The same calculation on shuffled lineage labels, which should return roughly zero; a share is read against it.
>
> **Classical support policy.** The share of isolates in repeated lineages; below 80 % the classical share is not admitted for interpretation; diagnostic values may remain in the record. This is separate from the optional population model's repeated-group policy.
>
> **Gate.** A condition fixed before the run. If it does not hold, the software withholds the estimate and names the condition that failed.
>
> **e-value.** A fixed-look e-value and a sequential e-process have different uses. Repeated-look validity requires the sequential construction and its stated null-model assumptions; a fixed-look value alone does not justify repeated inspection. Larger values are stronger evidence against the specified null.
>

## Notes on reading this report

- ‡ The 95 % interval includes zero. No lineage effect is distinguishable from none for this trait, and the point estimate must not be read on its own.
- † Support below 80 %. Too few isolates sit in lineages large enough to inform the estimate.
- § Derived arithmetically from a quantity recorded elsewhere in the record rather than estimated in this run.

amr-clonalshare 1.0.0 · record schema 1.0 · seed 42 · configuration sha256:f12b49aad646
Cite this run as: “amr-clonalshare 1.0.0, run E1E7B29E, record sha256:9c2eb99ec5e2.”
This report supersedes any earlier report bearing the same run identifier. It is regenerated from the record and holds no value that the record does not. The symbols carry the same wording in every run of this software; no symbol against a value means only that none of the listed conditions fired.
Classical collection-bootstrap endpoints are printed without clipping; the separately labeled species interval is floored at zero as part of its construction.
