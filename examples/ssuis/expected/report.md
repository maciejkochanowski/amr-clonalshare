# amr-clonalshare run report

Every number below is read from `clonal_share_result.json`, written by the same run, or from the release's validation grid where the text says so; none is recomputed here.

- **Run:** 000284F4
- **Issued:** 2026-09-25T15:49Z
- **Software:** amr-clonalshare 1.0.0
- **Isolates:** 677
- **Lineage:** baps_cluster
- **Seed:** 42
- **Configuration:** sha256:bfa04c09a84c
- **Record digest:** sha256:c8a0b5ce7cb1d8e25511f7bf9aa632cc656804ceded82c06c69b466d96f57365

**16 antimicrobials; 94 analysis outcomes.** 85 computed, 0 full ranges, 9 unavailable or incomplete. Each result retains its own target and data subset.

## 1. Measurement summary

**Quantity measured.** The share of the variation in the recorded binary outcome across this collection that lineage membership accounts for, read from a lineage label and an interpreted result and scored on isolates the estimator did not see.

- **Cohort:** 677 isolates in 30 lineages, typed by `baps_cluster`
- **Antimicrobials read:** 16
- **Membership shares estimable:** 13 of 13

**Gates.** 85 computed, 0 full ranges, 9 unavailable or incomplete. Each result retains its own target and data subset. A withheld value failed a declared reporting condition of its estimator; the recorded reason distinguishes input limitations from numerical failure. Passing a gate does not verify the model.

**Table 0.** Available analyses and their data subsets. A full range is a completed, uninformative result; an unavailable value is not zero.

| Agent | Analysis | N | Status | Data scope |
|---|---|---:|---|---|
| amoxicillin | finite_collection_prevalence | 677 | computed | recorded finite collection |
| cefquinome | finite_collection_prevalence | 677 | computed | recorded finite collection |
| ceftiofur | finite_collection_prevalence | 677 | computed | recorded finite collection |
| doxycycline | finite_collection_prevalence | 677 | computed | recorded finite collection |
| erythromycin | finite_collection_prevalence | 677 | computed | recorded finite collection |
| lincomycin | finite_collection_prevalence | 677 | computed | recorded finite collection |
| penicillin | finite_collection_prevalence | 677 | computed | recorded finite collection |
| spectinomycin | finite_collection_prevalence | 677 | computed | recorded finite collection |
| tetracycline | finite_collection_prevalence | 677 | computed | recorded finite collection |
| tiamulin | finite_collection_prevalence | 677 | computed | recorded finite collection |
| tilmicosin | finite_collection_prevalence | 677 | computed | recorded finite collection |
| trimethoprim | finite_collection_prevalence | 677 | computed | recorded finite collection |
| tylosin | finite_collection_prevalence | 677 | computed | recorded finite collection |
| amoxicillin | collection_membership | 677 | computed | observed and typed subset |
| cefquinome | collection_membership | 677 | computed | observed and typed subset |
| ceftiofur | collection_membership | 677 | computed | observed and typed subset |
| doxycycline | collection_membership | 677 | computed | observed and typed subset |
| erythromycin | collection_membership | 677 | computed | observed and typed subset |
| lincomycin | collection_membership | 677 | computed | observed and typed subset |
| penicillin | collection_membership | 677 | computed | observed and typed subset |
| spectinomycin | collection_membership | 677 | computed | observed and typed subset |
| tetracycline | collection_membership | 677 | computed | observed and typed subset |
| tiamulin | collection_membership | 677 | computed | observed and typed subset |
| tilmicosin | collection_membership | 677 | computed | observed and typed subset |
| trimethoprim | collection_membership | 677 | computed | observed and typed subset |
| tylosin | collection_membership | 677 | computed | observed and typed subset |
| amoxicillin | realised_component | 677 | unavailable | observed and typed subset |
| cefquinome | realised_component | 677 | unavailable | observed and typed subset |
| ceftiofur | realised_component | 677 | unavailable | observed and typed subset |
| doxycycline | realised_component | 677 | unavailable | observed and typed subset |
| erythromycin | realised_component | 677 | computed | observed and typed subset |
| lincomycin | realised_component | 677 | computed | observed and typed subset |
| penicillin | realised_component | 677 | unavailable | observed and typed subset |
| spectinomycin | realised_component | 677 | unavailable | observed and typed subset |
| tetracycline | realised_component | 677 | unavailable | observed and typed subset |
| tiamulin | realised_component | 677 | unavailable | observed and typed subset |
| tilmicosin | realised_component | 677 | computed | observed and typed subset |
| trimethoprim | realised_component | 677 | computed | observed and typed subset |
| tylosin | realised_component | 677 | computed | observed and typed subset |
| amoxicillin | population_probit | 677 | computed | population model conditional on retained grouped data |
| cefquinome | population_probit | 677 | computed | population model conditional on retained grouped data |
| ceftiofur | population_probit | 677 | computed | population model conditional on retained grouped data |
| doxycycline | population_probit | 677 | computed | population model conditional on retained grouped data |
| erythromycin | population_probit | 677 | computed | population model conditional on retained grouped data |
| lincomycin | population_probit | 677 | computed | population model conditional on retained grouped data |
| penicillin | population_probit | 677 | computed | population model conditional on retained grouped data |
| spectinomycin | population_probit | 677 | computed | population model conditional on retained grouped data |
| tetracycline | population_probit | 677 | computed | population model conditional on retained grouped data |
| tiamulin | population_probit | 677 | computed | population model conditional on retained grouped data |
| tilmicosin | population_probit | 677 | computed | population model conditional on retained grouped data |
| trimethoprim | population_probit | 677 | computed | population model conditional on retained grouped data |
| tylosin | population_probit | 677 | computed | population model conditional on retained grouped data |
| amoxicillin | censored_mic | 677 | computation_incomplete | population model fitted to readable and typed MIC subset |
| cefquinome | censored_mic | 677 | computed | population model fitted to readable and typed MIC subset |
| ceftiofur | censored_mic | 677 | computed | population model fitted to readable and typed MIC subset |
| doxycycline | censored_mic | 677 | computed | population model fitted to readable and typed MIC subset |
| enrofloxacin | censored_mic | 677 | computed | population model fitted to readable and typed MIC subset |
| erythromycin | censored_mic | 677 | computed | population model fitted to readable and typed MIC subset |
| florfenicol | censored_mic | 677 | computed | population model fitted to readable and typed MIC subset |
| lincomycin | censored_mic | 677 | computed | population model fitted to readable and typed MIC subset |
| marbofloxacin | censored_mic | 677 | computed | population model fitted to readable and typed MIC subset |
| penicillin | censored_mic | 677 | computed | population model fitted to readable and typed MIC subset |
| spectinomycin | censored_mic | 677 | computed | population model fitted to readable and typed MIC subset |
| tetracycline | censored_mic | 677 | computed | population model fitted to readable and typed MIC subset |
| tiamulin | censored_mic | 677 | computed | population model fitted to readable and typed MIC subset |
| tilmicosin | censored_mic | 677 | computed | population model fitted to readable and typed MIC subset |
| trimethoprim | censored_mic | 677 | computed | population model fitted to readable and typed MIC subset |
| tylosin | censored_mic | 677 | computed | population model fitted to readable and typed MIC subset |
| amoxicillin | lineage_evidence | 677 | computed | observed and typed subset |
| cefquinome | lineage_evidence | 677 | computed | observed and typed subset |
| ceftiofur | lineage_evidence | 677 | computed | observed and typed subset |
| doxycycline | lineage_evidence | 677 | computed | observed and typed subset |
| erythromycin | lineage_evidence | 677 | computed | observed and typed subset |
| lincomycin | lineage_evidence | 677 | computed | observed and typed subset |
| penicillin | lineage_evidence | 677 | computed | observed and typed subset |
| spectinomycin | lineage_evidence | 677 | computed | observed and typed subset |
| tetracycline | lineage_evidence | 677 | computed | observed and typed subset |
| tiamulin | lineage_evidence | 677 | computed | observed and typed subset |
| tilmicosin | lineage_evidence | 677 | computed | observed and typed subset |
| trimethoprim | lineage_evidence | 677 | computed | observed and typed subset |
| tylosin | lineage_evidence | 677 | computed | observed and typed subset |
| amoxicillin | sequential_lineage_evidence | 651 | computed | readable typed outcomes in prespecified ordered batches |
| cefquinome | sequential_lineage_evidence | 651 | computed | readable typed outcomes in prespecified ordered batches |
| ceftiofur | sequential_lineage_evidence | 651 | computed | readable typed outcomes in prespecified ordered batches |
| doxycycline | sequential_lineage_evidence | 651 | computed | readable typed outcomes in prespecified ordered batches |
| erythromycin | sequential_lineage_evidence | 651 | computed | readable typed outcomes in prespecified ordered batches |
| lincomycin | sequential_lineage_evidence | 651 | computed | readable typed outcomes in prespecified ordered batches |
| penicillin | sequential_lineage_evidence | 651 | computed | readable typed outcomes in prespecified ordered batches |
| spectinomycin | sequential_lineage_evidence | 651 | computed | readable typed outcomes in prespecified ordered batches |
| tetracycline | sequential_lineage_evidence | 651 | computed | readable typed outcomes in prespecified ordered batches |
| tiamulin | sequential_lineage_evidence | 651 | computed | readable typed outcomes in prespecified ordered batches |
| tilmicosin | sequential_lineage_evidence | 651 | computed | readable typed outcomes in prespecified ordered batches |
| trimethoprim | sequential_lineage_evidence | 651 | computed | readable typed outcomes in prespecified ordered batches |
| tylosin | sequential_lineage_evidence | 651 | computed | readable typed outcomes in prespecified ordered batches |

**Table 0b.** Reasons and next checks for results requiring interpretation. These checks do not change the recorded status or justify changing reporting gates.

| Agent | Analysis / status | Reason | Next check |
|---|---|---|---|
| amoxicillin | realised_component / unavailable | within-lineage residuals have excess kurtosis 14.18 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| cefquinome | realised_component / unavailable | within-lineage residuals have excess kurtosis 21.55 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| ceftiofur | realised_component / unavailable | within-lineage residuals have excess kurtosis 2.41 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| doxycycline | realised_component / unavailable | within-lineage residuals have excess kurtosis 1.23 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| penicillin | realised_component / unavailable | within-lineage residuals have excess kurtosis 3.18 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| spectinomycin | realised_component / unavailable | within-lineage residuals have excess kurtosis 4.92 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| tetracycline | realised_component / unavailable | within-lineage residuals have excess kurtosis 1.39 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| tiamulin | realised_component / unavailable | within-lineage residuals have excess kurtosis 2.06 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| amoxicillin | censored_mic / computation_incomplete | moment iteration limit reached; point retained for diagnostics only; 11 of 30 lineages lie wholly beyond the panel, holding 19% of isolates, below the reporting threshold of 50% | Inspect the recorded numerical or resource failure and method diagnostics before retrying in a new output directory. |
| cefquinome | censored_mic / computed | ok | Read the interval kind, assumptions and retained cohort before comparing results. |
| ceftiofur | censored_mic / computed | ok | Read the interval kind, assumptions and retained cohort before comparing results. |
| doxycycline | censored_mic / computed | ok | Read the interval kind, assumptions and retained cohort before comparing results. |
| enrofloxacin | censored_mic / computed | ok | Read the interval kind, assumptions and retained cohort before comparing results. |
| erythromycin | censored_mic / computed | 4 of 30 lineages lie wholly beyond the panel, holding 2% of isolates, below the reporting threshold of 50% | Read the interval kind, assumptions and retained cohort before comparing results. |
| florfenicol | censored_mic / computed | ok | Read the interval kind, assumptions and retained cohort before comparing results. |
| lincomycin | censored_mic / computed | 7 of 30 lineages lie wholly beyond the panel, holding 5% of isolates, below the reporting threshold of 50% | Read the interval kind, assumptions and retained cohort before comparing results. |
| marbofloxacin | censored_mic / computed | ok | Read the interval kind, assumptions and retained cohort before comparing results. |
| penicillin | censored_mic / computed | 4 of 30 lineages lie wholly beyond the panel, holding 8% of isolates, below the reporting threshold of 50% | Read the interval kind, assumptions and retained cohort before comparing results. |
| spectinomycin | censored_mic / computed | 1 of 30 lineages lie wholly beyond the panel, holding 0% of isolates, below the reporting threshold of 50% | Read the interval kind, assumptions and retained cohort before comparing results. |
| tetracycline | censored_mic / computed | ok | Read the interval kind, assumptions and retained cohort before comparing results. |
| tiamulin | censored_mic / computed | 3 of 30 lineages lie wholly beyond the panel, holding 1% of isolates, below the reporting threshold of 50% | Read the interval kind, assumptions and retained cohort before comparing results. |
| tilmicosin | censored_mic / computed | ok | Read the interval kind, assumptions and retained cohort before comparing results. |
| trimethoprim | censored_mic / computed | ok | Read the interval kind, assumptions and retained cohort before comparing results. |
| tylosin | censored_mic / computed | 6 of 30 lineages lie wholly beyond the panel, holding 4% of isolates, below the reporting threshold of 50% | Read the interval kind, assumptions and retained cohort before comparing results. |

**Table 0a.** Prevalence in the recorded collection: bounds allow every missing outcome to be either negative or positive. These are identification bounds, not confidence intervals. Observed prevalence uses observed outcomes only.

| Agent | Recorded | Observed | Missing | Observed prevalence | Collection bounds |
|---|---:|---:|---:|---:|---|
| amoxicillin | 677 | 677 | 0 | 5.5 % | 0.055 to 0.055 |
| cefquinome | 677 | 677 | 0 | 3.8 % | 0.038 to 0.038 |
| ceftiofur | 677 | 677 | 0 | 24.1 % | 0.241 to 0.241 |
| doxycycline | 677 | 677 | 0 | 84.3 % | 0.843 to 0.843 |
| erythromycin | 677 | 677 | 0 | 54.1 % | 0.541 to 0.541 |
| lincomycin | 677 | 677 | 0 | 62.6 % | 0.626 to 0.626 |
| penicillin | 677 | 677 | 0 | 23.0 % | 0.230 to 0.230 |
| spectinomycin | 677 | 677 | 0 | 11.5 % | 0.115 to 0.115 |
| tetracycline | 677 | 677 | 0 | 84.8 % | 0.848 to 0.848 |
| tiamulin | 677 | 677 | 0 | 19.4 % | 0.194 to 0.194 |
| tilmicosin | 677 | 677 | 0 | 53.8 % | 0.538 to 0.538 |
| trimethoprim | 677 | 677 | 0 | 28.1 % | 0.281 to 0.281 |
| tylosin | 677 | 677 | 0 | 54.1 % | 0.541 to 0.541 |

**Table 1.** Every trait, ordered by share. The reading in the last column is fixed by two conditions the record holds: whether the interval excludes zero and whether the e-BH procedure selected the trait at level 0.05.

| Trait | Share | 95 % interval | Control | e-value |  | Reading |
|---|---:|---:|---:|---:|---|---|
| penicillin | 0.512 | 0.244 to 0.693 | -0.06 | 4.7 × 10¹⁹ |  | evidence of a lineage effect |
| ceftiofur | 0.501 | 0.247 to 0.683 | -0.06 | 1.4 × 10¹⁸ |  | evidence of a lineage effect |
| tiamulin | 0.406 | 0.231 to 0.540 | -0.06 | 2.4 × 10¹⁴ |  | evidence of a lineage effect |
| trimethoprim | 0.304 | 0.166 to 0.393 | -0.06 | 1.1 × 10¹⁰ |  | evidence of a lineage effect |
| amoxicillin | 0.241 | -0.032 to 0.446 | -0.06 | 1.6 × 10⁵ | ‡ | lineage effect selected by e-BH; its size is not resolved |
| lincomycin | 0.236 | 0.119 to 0.354 | -0.06 | 2.0 × 10⁸ |  | evidence of a lineage effect |
| cefquinome | 0.188 | -0.028 to 0.401 | -0.06 | 1.9 × 10³ | ‡ | lineage effect selected by e-BH; its size is not resolved |
| erythromycin | 0.167 | 0.071 to 0.279 | -0.06 | 7.5 × 10⁴ |  | evidence of a lineage effect |
| tylosin | 0.158 | 0.062 to 0.256 | -0.06 | 5.9 × 10⁴ |  | evidence of a lineage effect |
| tilmicosin | 0.152 | 0.058 to 0.252 | -0.06 | 3.6 × 10⁴ |  | evidence of a lineage effect |
| spectinomycin | 0.129 | -0.079 to 0.327 | -0.06 | 2.3 × 10³ | ‡ | lineage effect selected by e-BH; its size is not resolved |
| tetracycline | 0.037 | -0.038 to 0.076 | -0.06 | 2.0 | ‡ | no detectable lineage effect |
| doxycycline | 0.030 | -0.046 to 0.066 | -0.06 | 2.7 | ‡ | no detectable lineage effect |

**Table 1b.** Lineage structure behind each share: lineages with at least two isolates, the effective number of lineages (inverse of the sum of squared lineage shares) and the share of isolates in lineages with at least two members (support). A high support with few effective lineages means the share rests on a few lineages.

| Trait | Lineages | Repeated | Effective | Support |
|---|---:|---:|---:|---:|
| penicillin | 30 | 27 | 10.2 | 99.6 % |
| ceftiofur | 30 | 27 | 10.2 | 99.6 % |
| tiamulin | 30 | 27 | 10.2 | 99.6 % |
| trimethoprim | 30 | 27 | 10.2 | 99.6 % |
| amoxicillin | 30 | 27 | 10.2 | 99.6 % |
| lincomycin | 30 | 27 | 10.2 | 99.6 % |
| cefquinome | 30 | 27 | 10.2 | 99.6 % |
| erythromycin | 30 | 27 | 10.2 | 99.6 % |
| tylosin | 30 | 27 | 10.2 | 99.6 % |
| tilmicosin | 30 | 27 | 10.2 | 99.6 % |
| spectinomycin | 30 | 27 | 10.2 | 99.6 % |
| tetracycline | 30 | 27 | 10.2 | 99.6 % |
| doxycycline | 30 | 27 | 10.2 | 99.6 % |

> **Not evaluated on this run, and why**
>
> Decomposition of a prevalence difference into lineage composition and within-lineage rate: needs two collections; this record holds one.
>

> **Interpretation (generated from Table 1 by fixed rules)**
>
> For penicillin and ceftiofur, lineage membership is associated with the recorded outcome: the estimated share is at or above one half in this collection at the chosen typing resolution.
>
> For tiamulin, trimethoprim, lincomycin, erythromycin, tylosin and tilmicosin, a lineage association is detected and the estimated share is below one half. Most variation is not explained by the lineage labels at this typing resolution.
>
> For amoxicillin, cefquinome and spectinomycin, the e-value selects a lineage effect while the interval for its size still reaches zero: there is evidence that positive outcomes are not spread evenly across the lineages, and this collection is too small, or too uneven across its lineages, to say how much of it the lineages carry. Read the e-value as the finding and the share as not yet resolved.
>
> For tetracycline and doxycycline, no lineage effect is distinguishable from none in this collection at the chosen resolution. This does not establish independence from lineage or identify the reason for a prevalence change.
>
> These readings describe association in the sampled collection. The analysis does not identify transmission, horizontal transfer, selection, or the effect of an intervention.
>

## 2. Admissibility of the input

- **Lineage column:** `baps_cluster`
- **Phenotype interpretation:** wt_nwt
- **Positive outcome:** NWT (non-wild-type)
- **Applied coding:** NWT (non-wild-type)
- **Intermediate policy:** not applicable
- **Interpretation source:** MIC against the cut-off of each agent in DATA_PROVENANCE.md section 2: EUCAST epidemiological cut-offs for tetracycline, doxycycline and erythromycin, collection-specific cut-offs for the other ten
- **AST standard/version:** not supplied / not supplied
- **Susceptibility calls:** 13 antimicrobials on 677 isolates
- **Recorded dilutions:** 10832 rows read, 10832 joined; 677 of 677 isolates carry a recorded dilution
- **Resampling:** 5 folds, 20 repeats, 400 bootstrap draws, 200 permutations per antimicrobial

Antimicrobials with fewer than 20 isolates of the rarer outcome: 0 of 13. Each method uses its own reporting conditions; sparse outcomes may yield wide intervals or an unavailable result.

Lineage groups in the input: 30, of which 3 hold a single isolate. Support 99.6 % against the 80.0 % the estimator needs: accepted.

*Figure 1 (drawn on the page).* Isolates per lineage, largest first. A lineage of one isolate, drawn in orange, cannot be predicted out of sample and counts against support; support is the share of isolates in the other lineages, 99.6 % here against the 80.0 % the estimator requires. The largest lineage holds 23.8 % of the isolates, which sets how much one lineage can weigh in the share.

**Table 2.** Conditions the estimator requires before any result is reported.

| Condition | Observed | Required | Verdict |
|---|---:|---:|---|
| Lineage support, lowest over the antimicrobials read | 99.6 % | ≥ 80.0 % | accepted |
| Lineage groups used | 30 | reported | accepted |

Each row uses only isolates with both a readable result for that agent and a recorded lineage. The last column is an algebraic support scenario: one genuinely new tested isolate in each of that many distinct singleton lineages. It does not guarantee interval precision, coverage, or overall estimability; it cannot repair a constant trait or a missing lineage contrast. Duplicating existing rows adds no evidence. The rarer-outcome count is a warning threshold, not an additional estimator gate.

**Table 2a.** Per-agent input feasibility and support scenario.

| Agent | Retained | Support | Input failures | Rarer outcome | New isolates for support |
|---|---|---|---|---|---|
| amoxicillin | 677 | 99.6 % | none at input level | 37 | 0 |
| cefquinome | 677 | 99.6 % | none at input level | 26 | 0 |
| ceftiofur | 677 | 99.6 % | none at input level | 163 | 0 |
| doxycycline | 677 | 99.6 % | none at input level | 106 | 0 |
| erythromycin | 677 | 99.6 % | none at input level | 311 | 0 |
| lincomycin | 677 | 99.6 % | none at input level | 253 | 0 |
| penicillin | 677 | 99.6 % | none at input level | 156 | 0 |
| spectinomycin | 677 | 99.6 % | none at input level | 78 | 0 |
| tetracycline | 677 | 99.6 % | none at input level | 103 | 0 |
| tiamulin | 677 | 99.6 % | none at input level | 131 | 0 |
| tilmicosin | 677 | 99.6 % | none at input level | 313 | 0 |
| trimethoprim | 677 | 99.6 % | none at input level | 190 | 0 |
| tylosin | 677 | 99.6 % | none at input level | 311 | 0 |

## 3. The clonal share, trait by trait

The intervals have a nominal 95 % level; their measured coverage depends on the cohort design and model assumptions. For **5 of the 13 traits** shown the interval includes zero, so no lineage effect is distinguishable from none for that trait.

The control column is the same estimator run on shuffled lineage labels; it should sit near zero, and a share is read against it rather than against zero.

*Figure 2 (drawn on the page).* Clonal share by trait, point estimate with 95 % interval, 13 largest of 13. Traits whose interval crosses zero are drawn in grey and marked ‡. The thin line beneath each interval is the species interval of Table 3, for a fresh draw of lineages; it is floored at zero by construction. The intervals are drawn as computed. The quantity lies between 0 and 1, so the part of an interval below zero carries no information: reading each interval as its overlap with that range leaves the coverage unchanged, and a lower limit at or below zero means the same thing either way.

On the release's validation grid, the interval this estimator prints for a binary trait contained the truth in 0.970 of 28,000 runs over 70 simulated cohorts (by cohort, 0.825 to 1.000). That figure belongs to the release, not to this run; it is what the phrase "95 % interval" was measured to mean.

Two intervals answer two questions. The interval above is for the share the lineages in this collection carry. The second interval below is for the share a fresh draw of lineages from the species would show, stated on the scale of the realised share of Table 3b below; it adds the sampling of the lineages themselves, on 29 degrees of freedom, is widened to the envelope of the first taken on that scale (with a lower end no smaller than zero, since a species share is not negative), and is the one to quote when the figure is read as a property of the species rather than of this collection. With few lineages it is markedly wider; with many the two nearly coincide.

For 2 traits one lineage carries more than half of the between-lineage variation (cefquinome 69.2 %; spectinomycin 58.2 %). The species interval describes lineage effects drawn from one law, and a collection in which one lineage carries the positive calls is not that: on the validation grid a carrier law of this kind took the species interval below its level with few lineages while the interval for the lineages in hand held. For them read the first interval as the statement about this collection and the species interval with that reservation.

The last column is an exploratory, model-equivalent restatement under a Gaussian probit threshold model. It transforms the collection's observed-scale share and interval at the estimated prevalence, treated as fixed. It does not estimate the actual latent variance of these particular lineages or a generic binomial mixed-model intraclass correlation. Monotonicity preserves coverage only for a matching target at a known, fixed prevalence; it does not guarantee coverage after estimating prevalence. Read the observed-scale interval as primary. Values at or below zero are displayed at the zero boundary.

**Table 3.** The two intervals for every trait: for lineage membership in this collection and for the species, with the share on the latent scale.

| Trait | Share | 95 % interval, these lineages | Species share (95 % interval) | Latent restatement (exploratory interval) |
|---|---:|---:|---:|---:|
| penicillin | 0.512 | 0.244 to 0.693 | 0.538 (0.264 to 0.714) | 0.750 (0.415 to 0.900) |
| ceftiofur | 0.501 | 0.247 to 0.683 | 0.527 (0.266 to 0.705) | 0.737 (0.415 to 0.892) |
| tiamulin | 0.406 | 0.231 to 0.540 | 0.431 (0.250 to 0.592) | 0.646 (0.411 to 0.787) |
| trimethoprim | 0.304 | 0.166 to 0.393 | 0.326 (0.181 to 0.488) | 0.485 (0.279 to 0.603) |
| amoxicillin | 0.241 | -0.032 to 0.446 | 0.261 (0.000 to 0.472) | 0.547 (0.000 to 0.776) |
| lincomycin | 0.236 | 0.119 to 0.354 | 0.255 (0.130 to 0.397) | 0.371 (0.191 to 0.535) |
| cefquinome | 0.188 | -0.028 to 0.401 | 0.204 (0.000 to 0.426) | 0.499 (0.000 to 0.757) |
| erythromycin | 0.167 | 0.071 to 0.279 | 0.182 (0.078 to 0.305) | 0.260 (0.111 to 0.425) |
| tylosin | 0.158 | 0.062 to 0.256 | 0.172 (0.068 to 0.291) | 0.246 (0.097 to 0.392) |
| tilmicosin | 0.152 | 0.058 to 0.252 | 0.166 (0.063 to 0.284) | 0.237 (0.090 to 0.386) |
| spectinomycin | 0.129 | -0.079 to 0.327 | 0.141 (0.000 to 0.350) | 0.287 (0.000 to 0.593) |
| tetracycline | 0.037 | -0.038 to 0.076 | 0.040 (0.000 to 0.116) | 0.082 (0.000 to 0.163) |
| doxycycline | 0.030 | -0.046 to 0.066 | 0.033 (0.000 to 0.110) | 0.066 (0.000 to 0.141) |

**The realised share of the same call.** The variance-component ratio of the lineages in hand, with an interval that is exact under a Gaussian within-lineage law. Its gate reads the excess kurtosis of the within-lineage residuals and withholds the interval above 0.99. Rare binary calls often fail this check; prevalence alone does not determine the residual kurtosis. The interval is approximate on binary data even when the gate opens; the estimate is printed either way, the interval only where the gate opened, for 5 of 13 traits here.

**Table 3b.** The realised share per trait and the verdict of its kurtosis gate.

| Trait | Realised share | Gaussian-model 95 % interval | Residual excess kurtosis | Gate |
|---|---:|---:|---:|---|
| penicillin | 0.554 | withheld | 3.18 | closed |
| ceftiofur | 0.553 | withheld | 2.41 | closed |
| tiamulin | 0.439 | withheld | 2.06 | closed |
| trimethoprim | 0.331 | 0.273 to 0.388 | 0.51 | open |
| amoxicillin | 0.312 | withheld | 14.18 | closed |
| cefquinome | 0.237 | withheld | 21.55 | closed |
| lincomycin | 0.229 | 0.173 to 0.288 | -1.11 | open |
| spectinomycin | 0.191 | withheld | 4.92 | closed |
| erythromycin | 0.166 | 0.114 to 0.223 | -1.38 | open |
| tylosin | 0.160 | 0.108 to 0.217 | -1.38 | open |
| tilmicosin | 0.155 | 0.103 to 0.211 | -1.41 | open |
| doxycycline | 0.034 | withheld | 1.23 | closed |
| tetracycline | 0.027 | withheld | 1.39 | closed |

## 4. Evidence that survives re-reading

A p-value is a statement about one look at the data. A surveillance panel is looked at again every year, and a p-value recomputed each time loses its guarantee. The e-value is evidence on a scale made for that: this run's e-value is a statement about this cohort, and a programme that adds an intake each year multiplies the e-value of each new intake, scored against the lineage rates learned from the earlier ones, into a running product whose guarantee holds at whatever intake it is read (sequential_e_process). The e-BH procedure controls the false-discovery rate across the panel whatever the dependence between traits. Larger is stronger; 1 is no evidence.

**Table 4.** e-value per trait and the e-BH selection at level 0.05. The selection threshold on this run is 23.6; a trait at or above it is selected.

| Trait | e-value | natural log | Selected |
|---|---:|---:|---|
| penicillin | 4.7 × 10¹⁹ | 45.29 | yes |
| ceftiofur | 1.4 × 10¹⁸ | 41.76 | yes |
| tiamulin | 2.4 × 10¹⁴ | 33.12 | yes |
| trimethoprim | 1.1 × 10¹⁰ | 23.13 | yes |
| lincomycin | 2.0 × 10⁸ | 19.09 | yes |
| amoxicillin | 1.6 × 10⁵ | 12.01 | yes |
| erythromycin | 7.5 × 10⁴ | 11.23 | yes |
| tylosin | 5.9 × 10⁴ | 10.98 | yes |
| tilmicosin | 3.6 × 10⁴ | 10.49 | yes |
| spectinomycin | 2.3 × 10³ | 7.74 | yes |
| cefquinome | 1.9 × 10³ | 7.56 | yes |
| doxycycline | 2.7 | 0.98 | no |
| tetracycline | 2.0 | 0.70 | no |

11 of 13 traits are selected. Note: e-value in the betting sense of Vovk and Wang, not the BLAST expectation value and not the E-value of VanderWeele and Ding.

*Figure 3 (drawn on the page).* Evidence per trait on the natural-log scale, 13 largest of 13. The dashed rule is 1/α, the evidence one trait alone needs at level 0.05; the solid rule is the e-BH selection threshold on this run, 23.6, which rises with the number of traits read together, and a trait at or beyond it is selected. Traits the e-BH procedure selected are drawn in blue; the scale is logarithmic, so equal steps are equal factors of evidence.

This cohort also carries intakes, 25 of them read in the order of collection_year, from 1 to 177 isolates each, with 26 isolates carrying no intake and set aside. Each intake is scored against the lineage rates learned from the intakes before it, and the running product is the sequential e-value: it may be read after any intake without spending the guarantee. 2 of 13 traits are selected on the product, against 11 at this single look, which is the price of a guarantee that holds at whatever intake the programme is read at.

**Table 4b.** The sequential e-value per trait after the last intake, and the e-BH selection taken on the log scale, where a product over intakes does not overflow.

| Trait | natural log of the product | Selected |
|---|---:|---|
| ceftiofur | 61.56 | yes |
| penicillin | 53.99 | yes |
| tiamulin | -9.05 | no |
| cefquinome | -10.02 | no |
| amoxicillin | -15.70 | no |
| trimethoprim | -35.37 | no |
| lincomycin | -93.49 | no |
| tilmicosin | -96.40 | no |
| doxycycline | -97.03 | no |
| tetracycline | -100.17 | no |
| tylosin | -101.10 | no |
| erythromycin | -105.09 | no |
| spectinomycin | -133.07 | no |

*Figure 4 (drawn on the page).* The running product of e-values over the 25 intakes, one line per trait, on the natural-log scale. Each intake is scored against the lineage rates learned from the intakes before it, so the first intakes carry no evidence and a line may fall as well as rise. The dashed rule is 1/α; the solid rule is the e-BH threshold on the product, 130.0. Traits selected on the product are drawn in blue and named at the right. Stopping-time panel control requires validity in the joint panel filtration; repeated rejection unions are not controlled. The frame is cut at −25: a line below it gives little or no evidence against lineage independence. It does not establish absence of a lineage effect.

## 5. Reading at the recorded resolution

Where a dilution was recorded, the share is also read from the dilution itself rather than from the call derived from it. A call and a dilution are both intervals on one concentration scale, so the dilution carries more of the reading; a reading at an end well is censored and enters as an interval. This table reads the share of that scale carried by lineage, beside the share of the binary call above. These are different quantities on independently retained readable and typed subsets; their denominators may differ.

**Table 5.** Share of the dilution scale carried by lineage, per agent, from the interval-censored likelihood. The share and its approximate F interval come from the moment iteration; that interval undercovered in simulation. The calibrated column gives the maximum-likelihood share with the interval validated by simulation for the Gaussian population model, subject to its sampling and coarsening assumptions. The realised interval is for the lineages this cohort holds; here it is an approximation built on the fitted mean-square ratio and not the exact interval of uncensored data.

| Agent | Share | Approximate F interval | Calibrated share and 95 % interval | Realised interval | Censored readings | Estimable |
|---|---:|---:|---:|---:|---:|---|
| amoxicillin | 0.557 | 0.426 to 0.704 | not computed | 0.489 to 0.616 | 89.1 % | no: moment iteration limit reached; point retained for diagnostics only; 11 of 30 lineages lie wholly beyond the panel, holding 19% of isolates, below the reporting threshold of 50% |
| cefquinome | 0.450 | 0.330 to 0.605 | not computed | 0.395 to 0.501 | 6.9 % | yes |
| ceftiofur | 0.672 | 0.558 to 0.791 | not computed | 0.634 to 0.705 | 1.8 % | yes |
| doxycycline | 0.070 | 0.028 to 0.150 | not computed | 0.031 to 0.120 | 2.5 % | yes |
| enrofloxacin | 0.221 | 0.137 to 0.357 | not computed | 0.164 to 0.280 | 0.7 % | yes |
| erythromycin | 0.228 | 0.143 to 0.367 | not computed | 0.171 to 0.288 | 56.4 % | yes |
| florfenicol | 0.240 | 0.151 to 0.382 | not computed | 0.181 to 0.300 | 3.8 % | yes |
| lincomycin | 0.357 | 0.247 to 0.513 | not computed | 0.297 to 0.414 | 55.8 % | yes |
| marbofloxacin | 0.036 | 0.006 to 0.097 | not computed | 0.006 to 0.080 | 5.5 % | yes |
| penicillin | 0.677 | 0.560 to 0.796 | not computed | 0.630 to 0.717 | 73.0 % | yes |
| spectinomycin | 0.461 | 0.341 to 0.615 | not computed | 0.407 to 0.511 | 8.7 % | yes |
| tetracycline | 0.077 | 0.033 to 0.161 | not computed | 0.037 to 0.128 | 10.6 % | yes |
| tiamulin | 0.698 | 0.588 to 0.809 | not computed | 0.663 to 0.728 | 10.0 % | yes |
| tilmicosin | 0.172 | 0.101 to 0.295 | not computed | 0.119 to 0.229 | 2.4 % | yes |
| trimethoprim | 0.437 | 0.318 to 0.592 | not computed | 0.381 to 0.488 | 11.5 % | yes |
| tylosin | 0.270 | 0.175 to 0.416 | not computed | 0.211 to 0.329 | 52.4 % | yes |

*Figure 5 (drawn on the page).* Share of the dilution scale carried by lineage, per agent, 16 largest of 16: the point estimate with its 95 % interval (thin) and the realised interval for these lineages (thick); an agent whose reading is refused is drawn in grey. No calibrated interval was computed on this run, so every interval drawn is the approximate F interval beside the moment estimate; that interval undercovered in simulation. The hollow diamond is the share of the binary call from Table 1 for the same agent, a different quantity whose retained isolates may differ, placed here so that the two readings can be seen side by side. The figure in parentheses at the right is the share of readings at an end well, which enter as censored.

The calibrated interval reads each laboratory on the wells that laboratory tested. The table below gives, per agent and laboratory, the panel the readings were taken to come from and the cut points the model uses; an end-well reading is censored at the panel edge of its own laboratory, not at the widest edge in the collection.

**Table 5a.** Panel geometry per agent and testing laboratory. Cut points are the boundaries between adjacent wells on which the likelihood is written; the admissible modes say which readings of the end wells the panel supports.

| Agent | Laboratory | Wells | Range | Cut points | Doubling | Admissible modes |
|---|---|---:|---:|---:|---|---|
| amoxicillin | LGC Fordham | 8 | 0.0312 to 4.0000 | none | yes | binary and interval |
| amoxicillin | OUCRU Ho Chi Minh City | 2 | 0.0156 to 0.0312 | none | yes | binary |
| cefquinome | LGC Fordham | 11 | 0.0020 to 2.0000 | none | yes | binary, interval and point |
| cefquinome | OUCRU Ho Chi Minh City | 3 | 0.0078 to 0.0312 | none | yes | binary and interval |
| ceftiofur | LGC Fordham | 10 | 0.0312 to 16.0000 | none | yes | binary, interval and point |
| ceftiofur | OUCRU Ho Chi Minh City | 4 | 0.0625 to 0.5000 | none | yes | binary, interval and point |
| doxycycline | LGC Fordham | 12 | 0.0312 to 64.0000 | none | yes | binary, interval and point |
| doxycycline | OUCRU Ho Chi Minh City | 9 | 0.0312 to 8.0000 | none | yes | binary and interval |
| enrofloxacin | LGC Fordham | 11 | 0.0078 to 8.0000 | none | yes | binary, interval and point |
| enrofloxacin | OUCRU Ho Chi Minh City | 5 | 0.2500 to 4.0000 | none | yes | binary, interval and point |
| erythromycin | LGC Fordham | 12 | 0.0156 to 32.0000 | none | yes | binary and interval |
| erythromycin | OUCRU Ho Chi Minh City | 12 | 0.0156 to 32.0000 | none | yes | binary and interval |
| florfenicol | LGC Fordham | 4 | 0.5000 to 4.0000 | none | yes | binary, interval and point |
| florfenicol | OUCRU Ho Chi Minh City | 5 | 0.2500 to 4.0000 | none | yes | binary, interval and point |
| lincomycin | LGC Fordham | 12 | 0.0625 to 128.0000 | none | yes | binary and interval |
| lincomycin | OUCRU Ho Chi Minh City | 12 | 0.0625 to 128.0000 | none | yes | binary and interval |
| marbofloxacin | LGC Fordham | 11 | 0.0156 to 16.0000 | none | yes | binary, interval and point |
| marbofloxacin | OUCRU Ho Chi Minh City | 3 | 0.2500 to 1.0000 | none | yes | binary and interval |
| penicillin | LGC Fordham | 10 | 0.0312 to 16.0000 | none | yes | binary and interval |
| penicillin | OUCRU Ho Chi Minh City | 5 | 0.0312 to 0.5000 | none | yes | binary and interval |
| spectinomycin | LGC Fordham | 10 | 1.0000 to 512.0000 | none | yes | binary and interval |
| spectinomycin | OUCRU Ho Chi Minh City | 8 | 2.0000 to 256.0000 | none | yes | binary and interval |
| tetracycline | LGC Fordham | 12 | 0.0625 to 128.0000 | none | yes | binary and interval |
| tetracycline | OUCRU Ho Chi Minh City | 9 | 0.5000 to 128.0000 | none | yes | binary and interval |
| tiamulin | LGC Fordham | 10 | 0.1250 to 64.0000 | none | yes | binary and interval |
| tiamulin | OUCRU Ho Chi Minh City | 7 | 0.0312 to 2.0000 | none | yes | binary, interval and point |
| tilmicosin | LGC Fordham | 11 | 0.2500 to 256.0000 | none | yes | binary, interval and point |
| tilmicosin | OUCRU Ho Chi Minh City | 12 | 0.5000 to 1024.0000 | none | yes | binary and interval |
| trimethoprim | LGC Fordham | 12 | 0.0156 to 32.0000 | none | yes | binary and interval |
| trimethoprim | OUCRU Ho Chi Minh City | 8 | 0.0625 to 8.0000 | none | yes | binary and interval |
| tylosin | LGC Fordham | 12 | 0.1250 to 256.0000 | none | yes | binary and interval |
| tylosin | OUCRU Ho Chi Minh City | 12 | 0.1250 to 256.0000 | none | yes | binary and interval |

**Table 5c.** Status of the calibrated interval per agent. The status names the outcome of the interval inversion; the reason is given whenever no interval is reported.

| Agent | Status | Method | Draws | Failed refits | Reason |
|---|---|---|---:|---:|---|
| amoxicillin | not computed |  |  |  | calibrated_interval is off in the configuration |
| cefquinome | not computed |  |  |  | calibrated_interval is off in the configuration |
| ceftiofur | not computed |  |  |  | calibrated_interval is off in the configuration |
| doxycycline | not computed |  |  |  | calibrated_interval is off in the configuration |
| enrofloxacin | not computed |  |  |  | calibrated_interval is off in the configuration |
| erythromycin | not computed |  |  |  | calibrated_interval is off in the configuration |
| florfenicol | not computed |  |  |  | calibrated_interval is off in the configuration |
| lincomycin | not computed |  |  |  | calibrated_interval is off in the configuration |
| marbofloxacin | not computed |  |  |  | calibrated_interval is off in the configuration |
| penicillin | not computed |  |  |  | calibrated_interval is off in the configuration |
| spectinomycin | not computed |  |  |  | calibrated_interval is off in the configuration |
| tetracycline | not computed |  |  |  | calibrated_interval is off in the configuration |
| tiamulin | not computed |  |  |  | calibrated_interval is off in the configuration |
| tilmicosin | not computed |  |  |  | calibrated_interval is off in the configuration |
| trimethoprim | not computed |  |  |  | calibrated_interval is off in the configuration |
| tylosin | not computed |  |  |  | calibrated_interval is off in the configuration |

*Figure 6 (drawn on the page).* Readings per lineage and dilution interval for 6 of 16 agents, the 14 largest of 30 lineages: each cell counts the isolates of one lineage whose reading fell in one interval of the panel, darker for more. Open intervals at the panel edges are the censored readings. A lineage whose readings sit in one or two adjacent cells contributes to the lineage share; readings spread along a row do not.

## 6. The same analyses within each stratum

The run was repeated on the isolates of each level of `isolation_country` (Canada (n = 205), United Kingdom (n = 423) and Vietnam (n = 49); 0 isolates without a level were left out). Each stratum is an analysis of its own, on its own lineages and its own panel, and the shares below are not adjusted for one another. A share that holds within every stratum is a property of the lineages; a share seen only in the pooled run and in no stratum is carried by the differences between strata.

**Table 5d.** Clonal share per stratum and agent: the share of the binary call with its interval, and the calibrated share of the dilution scale where a MIC table was supplied. The same rows are written to strata_results.csv.

| Stratum | Agent | n | Call share | 95 % interval | Calibrated MIC share (95 % interval) |
|---|---|---:|---:|---:|---:|
| Canada | amoxicillin | 205 | 0.259 | -0.059 to 0.533 | not computed |
| Canada | cefquinome | 205 | 0.181 | -0.060 to 0.486 | not computed |
| Canada | ceftiofur | 205 | 0.488 | 0.235 to 0.652 | not computed |
| Canada | doxycycline | 205 | -0.001 | -0.170 to 0.111 | not computed |
| Canada | enrofloxacin |  | not computed |  | not computed |
| Canada | erythromycin | 205 | 0.202 | 0.011 to 0.328 | not computed |
| Canada | florfenicol |  | not computed |  | not computed |
| Canada | lincomycin | 205 | 0.234 | 0.048 to 0.338 | not computed |
| Canada | marbofloxacin |  | not computed |  | not computed |
| Canada | penicillin | 205 | 0.549 | 0.213 to 0.745 | not computed |
| Canada | spectinomycin | 205 | 0.018 | -0.204 to 0.165 | not computed |
| Canada | tetracycline | 205 | 0.067 | 0.007 to 0.153 | not computed |
| Canada | tiamulin | 205 | 0.422 | 0.143 to 0.594 | not computed |
| Canada | tilmicosin | 205 | 0.174 | -0.034 to 0.307 | not computed |
| Canada | trimethoprim | 205 | 0.373 | 0.157 to 0.509 | not computed |
| Canada | tylosin | 205 | 0.176 | -0.021 to 0.308 | not computed |
| United Kingdom | amoxicillin | 423 | 0.271 | 0.036 to 0.499 | not computed |
| United Kingdom | cefquinome | 423 | 0.329 | 0.040 to 0.474 | not computed |
| United Kingdom | ceftiofur | 423 | 0.512 | 0.199 to 0.713 | not computed |
| United Kingdom | doxycycline | 423 | 0.105 | -0.039 to 0.263 | not computed |
| United Kingdom | enrofloxacin |  | not computed |  | not computed |
| United Kingdom | erythromycin | 423 | 0.229 | 0.081 to 0.349 | not computed |
| United Kingdom | florfenicol |  | not computed |  | not computed |
| United Kingdom | lincomycin | 423 | 0.354 | 0.185 to 0.577 | not computed |
| United Kingdom | marbofloxacin |  | not computed |  | not computed |
| United Kingdom | penicillin | 423 | 0.481 | 0.181 to 0.683 | not computed |
| United Kingdom | spectinomycin | 423 | 0.403 | -0.004 to 0.704 | not computed |
| United Kingdom | tetracycline | 423 | 0.105 | -0.039 to 0.263 | not computed |
| United Kingdom | tiamulin | 423 | 0.381 | 0.118 to 0.576 | not computed |
| United Kingdom | tilmicosin | 423 | 0.233 | 0.087 to 0.360 | not computed |
| United Kingdom | trimethoprim | 423 | 0.336 | 0.174 to 0.429 | not computed |
| United Kingdom | tylosin | 423 | 0.251 | 0.120 to 0.377 | not computed |
| Vietnam | amoxicillin | 49 | not computed |  | not computed |
| Vietnam | cefquinome | 49 | not computed |  | not computed |
| Vietnam | ceftiofur | 49 | not computed |  | not computed |
| Vietnam | doxycycline | 49 | not computed |  | not computed |
| Vietnam | enrofloxacin |  | not computed |  | not computed |
| Vietnam | erythromycin | 49 | not computed |  | not computed |
| Vietnam | florfenicol |  | not computed |  | not computed |
| Vietnam | lincomycin | 49 | not computed |  | not computed |
| Vietnam | marbofloxacin |  | not computed |  | not computed |
| Vietnam | penicillin | 49 | not computed |  | not computed |
| Vietnam | spectinomycin | 49 | not computed |  | not computed |
| Vietnam | tetracycline | 49 | not computed |  | not computed |
| Vietnam | tiamulin | 49 | not computed |  | not computed |
| Vietnam | tilmicosin | 49 | not computed |  | not computed |
| Vietnam | trimethoprim | 49 | not computed |  | not computed |
| Vietnam | tylosin | 49 | not computed |  | not computed |

## 7. A change of lineages or a change within them

A prevalence difference between two collections can be split into a change in lineage composition and a change in rate within lineages. That decomposition needs two collections and was not run here: this record holds one.

## 8. Population liability ICC

This optional result targets the Gaussian random-intercept population liability ICC, with marginal prevalence jointly fitted under a conditional binomial model. Independent Gaussian group effects and group sizes fixed independently of those effects are assumptions. Selection into the modeled population must be ignorable. The model does not identify transmission or intervention effects.

**Table P1.** Population-model profile sets; finite-grid calibration is not a universal 95% coverage guarantee. These results are separate from classical lineage-membership estimates and admissions.

| Agent | Population liability ICC | Model profile interval | Fitted prevalence | Status |
|---|---:|---:|---:|---|
| amoxicillin | 0.584 | 0.303 to 0.843 | 0.105 | ok |
| cefquinome | 0.522 | 0.222 to 0.836 | 0.063 | ok |
| ceftiofur | 0.755 | 0.554 to 0.904 | 0.401 | ok |
| doxycycline | 0.124 | 0.012 to 0.380 | 0.853 | ok |
| erythromycin | 0.355 | 0.167 to 0.620 | 0.591 | ok |
| lincomycin | 0.459 | 0.253 to 0.708 | 0.713 | ok |
| penicillin | 0.722 | 0.517 to 0.884 | 0.391 | ok |
| spectinomycin | 0.327 | 0.139 to 0.587 | 0.171 | ok |
| tetracycline | 0.113 | 0.009 to 0.358 | 0.861 | ok |
| tiamulin | 0.655 | 0.432 to 0.849 | 0.294 | ok |
| tilmicosin | 0.320 | 0.146 to 0.578 | 0.576 | ok |
| trimethoprim | 0.521 | 0.305 to 0.751 | 0.359 | ok |
| tylosin | 0.331 | 0.153 to 0.589 | 0.581 | ok |

- **Fixed LR cutoff:** 5.443134202054537
- **Calibration:** population-model-calibration-20260926101
- **Calibration SHA-256:** 2c2206be26ea5f5ca88bcb07964ec64fded065676291df9a1e2a68b4a0e911e4
- **Independent validation:** population-model-validation-20260926202
- **Validation status:** accepted_finite_grid_model_dependent
- **Validation summary SHA-256:** e026f78eb7559ec1de3b914377d8c052eebc1aa45cb64042881a04c6d634faf0

Constant outcomes or too few repeated groups retain an explicitly uninformative [0,1] set with no point estimate. Numerical failure has no estimate or interval. The repeated-group policy is separate from the support gate of the clonal share. Counts, exclusions, numerical checks and the complete model-domain diagnostics are in the result record.

doxycycline: extrapolation warning — fitted prevalence outside the tested finite-grid range.

tetracycline: extrapolation warning — fitted prevalence outside the tested finite-grid range.

Tested finite-grid group-count range: [15, 200]; maximum tested group size: 161; model prevalence range: [0.02, 0.85]. Range inclusion does not establish validation for every design within these bounds.

## 9. Provenance and terms

- **Software:** amr-clonalshare 1.0.0, record schema 1.0
- **Seed:** 42
- **Configuration:** `sha256:bfa04c09a84c`
- **Record:** `sha256:c8a0b5ce7cb1d8e25511f7bf9aa632cc656804ceded82c06c69b466d96f57365`

> **Terms used in this report**
>
> **Share.** How much of the variation of the call between isolates is associated with their recorded lineage. Near 1: a strong association on the measured scale. Near 0: little association on that scale. This does not identify transmission or its mechanism.
>
> **Collection-bootstrap interval.** The lineage-membership interval targets the sampled collection. Separate Gaussian-model component and species intervals retain their own targets and assumptions; coverage measured on a finite simulation grid is not a universal guarantee.
>
> **Population-model interval.** The grouped-probit profile targets a Gaussian population liability ICC, using its recorded fixed cutoff and finite-grid calibration/validation status. It does not target the sampled collection's membership share.
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

amr-clonalshare 1.0.0 · record schema 1.0 · seed 42 · configuration sha256:bfa04c09a84c
Cite this run as: “amr-clonalshare 1.0.0, run 000284F4, record sha256:c8a0b5ce7cb1.”
This report supersedes any earlier report bearing the same run identifier. It is regenerated from the record and holds no value that the record does not. The symbols carry the same wording in every run of this software; no symbol against a value means only that none of the listed conditions fired.
Classical collection-bootstrap endpoints are printed without clipping; the separately labeled species interval is floored at zero as part of its construction.
Population-model profile endpoints are computed within [0,1]. A numerical failure has missing endpoints, and an uninformative full set is labeled separately. Its population target and finite-grid evidence are independent of the classical collection-bootstrap interval.
