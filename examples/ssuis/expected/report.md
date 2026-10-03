# amr-clonalshare run report

Every number below is read from `clonal_share_result.json`, written by the same run; none is recomputed here.

- **Run:** 7DEE408A
- **Issued:** 2026-10-03T06:13Z
- **Software:** amr-clonalshare 1.0.0
- **Isolates:** 677
- **Lineage:** baps_cluster
- **Seed:** 42
- **Configuration:** sha256:7675bb625cb2
- **Record digest:** sha256:3ee3d16084cc3520502340537a20203de1e6d99fcbd49914274e3e1ee1920b46

**16 antimicrobials; 126 analysis outcomes.** 124 computed, 2 full ranges, 0 unavailable. Each result retains its own target and data subset.

## 1. Measurement summary

**Quantity measured.** The share of the variation in the recorded binary outcome across this collection that lineage membership accounts for, read from a lineage label and an interpreted result and scored on isolates the estimator did not see.

- **Collection:** 677 isolates in 30 lineages, typed by `baps_cluster`
- **Antimicrobials read:** 16
- **Membership shares estimable:** 13 of 13

**Strata of the calls.** The analyses of every call were read within the 3 levels of `isolation_country`: lineage labels are exchanged only within a level, the interval draws each isolate's level with its call, and the bounds a call places on the latent ordering are the bounds within levels, so a shift between levels is not read as a lineage difference. The e-values take no strata.

**Reporting conditions.** 124 computed, 2 full ranges, 0 unavailable. Each result retains its own target and data subset. An unavailable value failed a declared reporting condition of its estimator; the recorded reason distinguishes input limitations from numerical failure.

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
| amoxicillin | call_latent_bounds | 677 | computed | observed and typed subset |
| amoxicillin | call_latent_lower_limit | 677 | computed | observed and typed subset |
| cefquinome | collection_membership | 677 | computed | observed and typed subset |
| cefquinome | call_latent_bounds | 677 | computed | observed and typed subset |
| cefquinome | call_latent_lower_limit | 677 | computed | observed and typed subset |
| ceftiofur | collection_membership | 677 | computed | observed and typed subset |
| ceftiofur | call_latent_bounds | 677 | computed | observed and typed subset |
| ceftiofur | call_latent_lower_limit | 677 | computed | observed and typed subset |
| doxycycline | collection_membership | 677 | computed | observed and typed subset |
| doxycycline | call_latent_bounds | 677 | full_range | observed and typed subset |
| doxycycline | call_latent_lower_limit | 677 | computed | observed and typed subset |
| erythromycin | collection_membership | 677 | computed | observed and typed subset |
| erythromycin | call_latent_bounds | 677 | computed | observed and typed subset |
| erythromycin | call_latent_lower_limit | 677 | computed | observed and typed subset |
| lincomycin | collection_membership | 677 | computed | observed and typed subset |
| lincomycin | call_latent_bounds | 677 | computed | observed and typed subset |
| lincomycin | call_latent_lower_limit | 677 | computed | observed and typed subset |
| penicillin | collection_membership | 677 | computed | observed and typed subset |
| penicillin | call_latent_bounds | 677 | computed | observed and typed subset |
| penicillin | call_latent_lower_limit | 677 | computed | observed and typed subset |
| spectinomycin | collection_membership | 677 | computed | observed and typed subset |
| spectinomycin | call_latent_bounds | 677 | computed | observed and typed subset |
| spectinomycin | call_latent_lower_limit | 677 | computed | observed and typed subset |
| tetracycline | collection_membership | 677 | computed | observed and typed subset |
| tetracycline | call_latent_bounds | 677 | full_range | observed and typed subset |
| tetracycline | call_latent_lower_limit | 677 | computed | observed and typed subset |
| tiamulin | collection_membership | 677 | computed | observed and typed subset |
| tiamulin | call_latent_bounds | 677 | computed | observed and typed subset |
| tiamulin | call_latent_lower_limit | 677 | computed | observed and typed subset |
| tilmicosin | collection_membership | 677 | computed | observed and typed subset |
| tilmicosin | call_latent_bounds | 677 | computed | observed and typed subset |
| tilmicosin | call_latent_lower_limit | 677 | computed | observed and typed subset |
| trimethoprim | collection_membership | 677 | computed | observed and typed subset |
| trimethoprim | call_latent_bounds | 677 | computed | observed and typed subset |
| trimethoprim | call_latent_lower_limit | 677 | computed | observed and typed subset |
| tylosin | collection_membership | 677 | computed | observed and typed subset |
| tylosin | call_latent_bounds | 677 | computed | observed and typed subset |
| tylosin | call_latent_lower_limit | 677 | computed | observed and typed subset |
| amoxicillin | mic_order | 677 | computed | readable and typed MIC subset |
| amoxicillin | mic_latent_bounds | 677 | computed | readable and typed MIC subset |
| amoxicillin | mic_latent_lower_limit | 677 | computed | readable and typed MIC subset |
| cefquinome | mic_order | 677 | computed | readable and typed MIC subset |
| cefquinome | mic_latent_bounds | 677 | computed | readable and typed MIC subset |
| cefquinome | mic_latent_lower_limit | 677 | computed | readable and typed MIC subset |
| ceftiofur | mic_order | 677 | computed | readable and typed MIC subset |
| ceftiofur | mic_latent_bounds | 677 | computed | readable and typed MIC subset |
| ceftiofur | mic_latent_lower_limit | 677 | computed | readable and typed MIC subset |
| doxycycline | mic_order | 677 | computed | readable and typed MIC subset |
| doxycycline | mic_latent_bounds | 677 | computed | readable and typed MIC subset |
| doxycycline | mic_latent_lower_limit | 677 | computed | readable and typed MIC subset |
| enrofloxacin | mic_order | 677 | computed | readable and typed MIC subset |
| enrofloxacin | mic_latent_bounds | 677 | computed | readable and typed MIC subset |
| enrofloxacin | mic_latent_lower_limit | 677 | computed | readable and typed MIC subset |
| erythromycin | mic_order | 677 | computed | readable and typed MIC subset |
| erythromycin | mic_latent_bounds | 677 | computed | readable and typed MIC subset |
| erythromycin | mic_latent_lower_limit | 677 | computed | readable and typed MIC subset |
| florfenicol | mic_order | 677 | computed | readable and typed MIC subset |
| florfenicol | mic_latent_bounds | 677 | computed | readable and typed MIC subset |
| florfenicol | mic_latent_lower_limit | 677 | computed | readable and typed MIC subset |
| lincomycin | mic_order | 677 | computed | readable and typed MIC subset |
| lincomycin | mic_latent_bounds | 677 | computed | readable and typed MIC subset |
| lincomycin | mic_latent_lower_limit | 677 | computed | readable and typed MIC subset |
| marbofloxacin | mic_order | 677 | computed | readable and typed MIC subset |
| marbofloxacin | mic_latent_bounds | 677 | computed | readable and typed MIC subset |
| marbofloxacin | mic_latent_lower_limit | 677 | computed | readable and typed MIC subset |
| penicillin | mic_order | 677 | computed | readable and typed MIC subset |
| penicillin | mic_latent_bounds | 677 | computed | readable and typed MIC subset |
| penicillin | mic_latent_lower_limit | 677 | computed | readable and typed MIC subset |
| spectinomycin | mic_order | 677 | computed | readable and typed MIC subset |
| spectinomycin | mic_latent_bounds | 677 | computed | readable and typed MIC subset |
| spectinomycin | mic_latent_lower_limit | 677 | computed | readable and typed MIC subset |
| tetracycline | mic_order | 677 | computed | readable and typed MIC subset |
| tetracycline | mic_latent_bounds | 677 | computed | readable and typed MIC subset |
| tetracycline | mic_latent_lower_limit | 677 | computed | readable and typed MIC subset |
| tiamulin | mic_order | 677 | computed | readable and typed MIC subset |
| tiamulin | mic_latent_bounds | 677 | computed | readable and typed MIC subset |
| tiamulin | mic_latent_lower_limit | 677 | computed | readable and typed MIC subset |
| tilmicosin | mic_order | 677 | computed | readable and typed MIC subset |
| tilmicosin | mic_latent_bounds | 677 | computed | readable and typed MIC subset |
| tilmicosin | mic_latent_lower_limit | 677 | computed | readable and typed MIC subset |
| trimethoprim | mic_order | 677 | computed | readable and typed MIC subset |
| trimethoprim | mic_latent_bounds | 677 | computed | readable and typed MIC subset |
| trimethoprim | mic_latent_lower_limit | 677 | computed | readable and typed MIC subset |
| tylosin | mic_order | 677 | computed | readable and typed MIC subset |
| tylosin | mic_latent_bounds | 677 | computed | readable and typed MIC subset |
| tylosin | mic_latent_lower_limit | 677 | computed | readable and typed MIC subset |
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

**Table 0b.** Reasons and next checks for results requiring interpretation. These checks do not change the recorded status or justify changing reporting conditions.

| Agent | Analysis / status | Reason | Next check |
|---|---|---|---|
| amoxicillin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.897 | Read the interval kind, assumptions and retained collection before comparing results. |
| cefquinome | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.900 | Read the interval kind, assumptions and retained collection before comparing results. |
| ceftiofur | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.869 | Read the interval kind, assumptions and retained collection before comparing results. |
| doxycycline | call_latent_bounds / full_range | The largest share attained by an arrangement of the lineages inside the readings is 0.742 | Report the full range as uninformative; inspect repeated-lineage support and outcome variation. |
| erythromycin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.645 | Read the interval kind, assumptions and retained collection before comparing results. |
| lincomycin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.758 | Read the interval kind, assumptions and retained collection before comparing results. |
| penicillin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.863 | Read the interval kind, assumptions and retained collection before comparing results. |
| spectinomycin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.839 | Read the interval kind, assumptions and retained collection before comparing results. |
| tetracycline | call_latent_bounds / full_range | The largest share attained by an arrangement of the lineages inside the readings is 0.746 | Report the full range as uninformative; inspect repeated-lineage support and outcome variation. |
| tiamulin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.833 | Read the interval kind, assumptions and retained collection before comparing results. |
| tilmicosin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.639 | Read the interval kind, assumptions and retained collection before comparing results. |
| trimethoprim | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.762 | Read the interval kind, assumptions and retained collection before comparing results. |
| tylosin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.644 | Read the interval kind, assumptions and retained collection before comparing results. |
| amoxicillin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.860 | Read the interval kind, assumptions and retained collection before comparing results. |
| cefquinome | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.596 | Read the interval kind, assumptions and retained collection before comparing results. |
| ceftiofur | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.704 | Read the interval kind, assumptions and retained collection before comparing results. |
| doxycycline | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.435 | Read the interval kind, assumptions and retained collection before comparing results. |
| enrofloxacin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.577 | Read the interval kind, assumptions and retained collection before comparing results. |
| erythromycin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.382 | Read the interval kind, assumptions and retained collection before comparing results. |
| florfenicol | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.812 | Read the interval kind, assumptions and retained collection before comparing results. |
| lincomycin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.596 | Read the interval kind, assumptions and retained collection before comparing results. |
| marbofloxacin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.506 | Read the interval kind, assumptions and retained collection before comparing results. |
| penicillin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.817 | Read the interval kind, assumptions and retained collection before comparing results. |
| spectinomycin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.524 | Read the interval kind, assumptions and retained collection before comparing results. |
| tetracycline | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.438 | Read the interval kind, assumptions and retained collection before comparing results. |
| tiamulin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.593 | Read the interval kind, assumptions and retained collection before comparing results. |
| tilmicosin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.589 | Read the interval kind, assumptions and retained collection before comparing results. |
| trimethoprim | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.365 | Read the interval kind, assumptions and retained collection before comparing results. |
| tylosin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.508 | Read the interval kind, assumptions and retained collection before comparing results. |

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

**Table 1.** Every trait, ordered by share, with the 95 % interval, the permutation p-value and the Benjamini-Yekutieli q-value across traits; selected is the Benjamini-Yekutieli selection at a false-discovery level of 0.05, which holds whatever the dependence between traits. Control is the share on shuffled lineage labels. A p-value written with ≤ is the smallest the permutations can give. The last column gives the conclusion of the results table of the form: established, lineage structure established, when the interval lies above zero and the trait is selected, or the lower confidence limit is above zero; not established; whole range, when the interval spans it; or not estimable. The conclusion says whether this collection shows lineage structure in the measurement, with the values behind it. It is not a statement about transmission or a resistance mechanism.

| Trait | Share | 95 % interval | p | q | Control | e-value | Selected | Conclusion |
|---|---:|---:|---:|---:|---:|---:|---|---|
| penicillin | 0.510 | 0.448 to 0.586 | ≤ 0.001 | 0.003 | 0.00 | 4.7 × 10¹⁹ | yes | established |
| ceftiofur | 0.498 | 0.440 to 0.568 | ≤ 0.001 | 0.003 | 0.00 | 1.4 × 10¹⁸ | yes | established |
| tiamulin | 0.387 | 0.309 to 0.474 | ≤ 0.001 | 0.003 | 0.00 | 2.4 × 10¹⁴ | yes | established |
| trimethoprim | 0.283 | 0.219 to 0.362 | ≤ 0.001 | 0.003 | 0.00 | 1.1 × 10¹⁰ | yes | established |
| amoxicillin | 0.263 | 0.149 to 0.398 | ≤ 0.001 | 0.003 | 0.00 | 1.6 × 10⁵ | yes | established |
| lincomycin | 0.224 | 0.169 to 0.281 | ≤ 0.001 | 0.003 | 0.00 | 2.0 × 10⁸ | yes | established |
| cefquinome | 0.211 | 0.048 to 0.386 | ≤ 0.001 | 0.003 | 0.00 | 1.9 × 10³ | yes | established |
| spectinomycin | 0.181 | 0.115 to 0.258 | ≤ 0.001 | 0.003 | 0.00 | 2.3 × 10³ | yes | established |
| erythromycin | 0.153 | 0.103 to 0.204 | ≤ 0.001 | 0.003 | 0.00 | 7.5 × 10⁴ | yes | established |
| tylosin | 0.151 | 0.100 to 0.201 | ≤ 0.001 | 0.003 | 0.00 | 5.9 × 10⁴ | yes | established |
| tilmicosin | 0.149 | 0.098 to 0.198 | ≤ 0.001 | 0.003 | 0.00 | 3.6 × 10⁴ | yes | established |
| doxycycline | 0.035 | 0.001 to 0.078 | ≤ 0.001 | 0.003 | 0.00 | 2.7 | yes | established |
| tetracycline | 0.035 | 0.005 to 0.073 | 0.006 | 0.019 | 0.00 | 2.0 | yes | established |

**Table 1b.** Lineage structure behind each share: the lineages with at least two isolates on which the share is scored, the isolates of singleton lineages set aside, the effective number of scored lineages (inverse of the sum of squared lineage shares) and the share of isolates in lineages with at least two members (support). Few effective lineages mean the share rests on a few lineages.

| Trait | Lineages scored | Singletons set aside | Effective | Support |
|---|---:|---:|---:|---:|
| penicillin | 27 | 3 | 10.1 | 99.6 % |
| ceftiofur | 27 | 3 | 10.1 | 99.6 % |
| tiamulin | 27 | 3 | 10.1 | 99.6 % |
| trimethoprim | 27 | 3 | 10.1 | 99.6 % |
| amoxicillin | 27 | 3 | 10.1 | 99.6 % |
| lincomycin | 27 | 3 | 10.1 | 99.6 % |
| cefquinome | 27 | 3 | 10.1 | 99.6 % |
| spectinomycin | 27 | 3 | 10.1 | 99.6 % |
| erythromycin | 27 | 3 | 10.1 | 99.6 % |
| tylosin | 27 | 3 | 10.1 | 99.6 % |
| tilmicosin | 27 | 3 | 10.1 | 99.6 % |
| doxycycline | 27 | 3 | 10.1 | 99.6 % |
| tetracycline | 27 | 3 | 10.1 | 99.6 % |

> **Not evaluated on this run, and why**
>
> Decomposition of a prevalence difference into lineage composition and within-lineage rate: needs two collections; this record holds one.
>

> **Interpretation (generated from Table 1 by fixed rules)**
>
> For penicillin, lineage membership is associated with the recorded outcome: the estimated share is at or above one half in this collection at the chosen typing resolution.
>
> For ceftiofur, tiamulin, trimethoprim, amoxicillin, lincomycin, cefquinome, spectinomycin, erythromycin, tylosin, tilmicosin, doxycycline and tetracycline, a lineage association is detected and the estimated share is below one half. Most variation is not explained by the lineage labels at this typing resolution.
>
> These readings describe association in the sampled collection. The analysis does not identify transmission, horizontal transfer, selection, or the effect of an intervention.
>

## 2. Admissibility of the input

- **Lineage column:** `baps_cluster`
- **Phenotype interpretation:** binary
- **Positive outcome:** MIC above the cut-off of its agent: non-wild-type for tetracycline (EUCAST ECOFF), doxycycline and erythromycin (tentative EUCAST ECOFFs); above a cut-off derived from this collection for the other ten
- **Applied coding:** 1 (positive)
- **Intermediate policy:** not applicable
- **Interpretation source:** MIC against the cut-off of each agent in DATA_PROVENANCE.md section 2: EUCAST epidemiological cut-offs for tetracycline, doxycycline and erythromycin, collection-specific cut-offs for the other ten
- **AST standard/version:** not supplied / not supplied
- **Susceptibility calls:** 13 antimicrobials on 677 isolates
- **Recorded dilutions:** 10832 rows read, 10832 joined; 677 of 677 isolates carry a recorded dilution
- **Resampling:** 5 folds, 20 repeats, 999 bootstrap draws, 999 permutations per antimicrobial

Antimicrobials with fewer than 20 isolates of the rarer outcome: 0 of 13. Each method uses its own reporting conditions; sparse outcomes may yield wide intervals or an unavailable result.

Lineage groups in the input: 30, of which 3 hold a single isolate. Support, the share of isolates in lineages of at least two, is 99.6 %; the share is scored on those isolates and the singletons are set aside. At least two lineages hold two or more isolates.

*Figure 1 (drawn on the page).* Isolates per lineage, largest first. A lineage of one isolate, drawn in orange, cannot be predicted out of sample and is set aside; support is the share of isolates in the other lineages, 99.6 % here. The largest lineage holds 23.8 % of the isolates, which sets how much one lineage can weigh in the share.

**Table 2.** Conditions the estimator requires before any result is reported.

| Condition | Observed | Required | Verdict |
|---|---:|---:|---|
| Lineages with at least two isolates, fewest over the antimicrobials read | 27 | ≥ 2 | accepted |
| Support, lowest over the antimicrobials read | 99.6 % | reported | singletons set aside |

Each row uses only isolates with both a readable result for that agent and a recorded lineage. Support is the share of them in lineages of at least two isolates; the share is scored on those isolates and the singletons are set aside. The rarer-outcome count is a warning threshold, not a reporting condition of the estimator.

**Table 2a.** Per-agent input feasibility.

| Agent | Retained | Support | Input failures | Rarer outcome |
|---|---|---|---|---|
| amoxicillin | 677 | 99.6 % | none at input level | 37 |
| cefquinome | 677 | 99.6 % | none at input level | 26 |
| ceftiofur | 677 | 99.6 % | none at input level | 163 |
| doxycycline | 677 | 99.6 % | none at input level | 106 |
| erythromycin | 677 | 99.6 % | none at input level | 311 |
| lincomycin | 677 | 99.6 % | none at input level | 253 |
| penicillin | 677 | 99.6 % | none at input level | 156 |
| spectinomycin | 677 | 99.6 % | none at input level | 78 |
| tetracycline | 677 | 99.6 % | none at input level | 103 |
| tiamulin | 677 | 99.6 % | none at input level | 131 |
| tilmicosin | 677 | 99.6 % | none at input level | 313 |
| trimethoprim | 677 | 99.6 % | none at input level | 190 |
| tylosin | 677 | 99.6 % | none at input level | 311 |

## 3. The lineage share of the call, trait by trait

The interval is for the share of the lineages in this collection: the lineages and their sizes are held fixed and the isolates of each lineage are drawn again from the smoothed distribution of its calls. For **0 of the 13 traits** shown the interval reaches zero, so no lineage effect is distinguishable from none for that trait.

The control column is the same estimator run on shuffled lineage labels; it should sit near zero, and a share is read against it rather than against zero.

*Figure 2 (drawn on the page).* Lineage share of the call by trait, point estimate with 95 % interval for the represented lineages, 13 largest of 13. Traits whose interval reaches zero are drawn in grey and marked ‡.

The last column asks what the call says about the ordering it was cut from: the lineage share of the latent MIC, or of any continuous value the call thresholds, the rank intraclass correlation of that value. A call does not identify that share; it bounds it. The bracket gives the smallest share any ordering consistent with the calls allows, the lower bound the calls establish, and the largest; an upper end marked ≤ is a certified bound rather than a share attained. The figure after it is the one-sided 95 % lower confidence limit of the lower bound, from splitting the isolates of every lineage at random, choosing a direction on one half and testing it on the other: a share of the latent ordering that the lineages of this collection are shown to account for, with no model for the latent values beyond their agreeing with the readings. A MIC reading on a panel that contains the cut-off refines the call (Table 5), and refining the readings can only narrow the bounds they allow.

**Table 3.** The share of every trait with its interval, and the bounds the call places on the share of the latent ordering.

| Trait | Share | 95 % interval, represented lineages | Latent ordering: bounds; 95 % lower limit |
|---|---:|---:|---:|
| penicillin | 0.510 | 0.448 to 0.586 | 0.120 to ≤ 1.000; ≥ 0.035 |
| ceftiofur | 0.498 | 0.440 to 0.568 | 0.104 to ≤ 1.000; ≥ 0.020 |
| tiamulin | 0.387 | 0.309 to 0.474 | 0.039 to ≤ 1.000; ≥ 0.000 |
| trimethoprim | 0.283 | 0.219 to 0.362 | 0.020 to ≤ 1.000; ≥ 0.000 |
| amoxicillin | 0.263 | 0.149 to 0.398 | 0.014 to ≤ 1.000; ≥ 0.000 |
| lincomycin | 0.224 | 0.169 to 0.281 | 0.000 to ≤ 1.000; ≥ 0.000 |
| cefquinome | 0.211 | 0.048 to 0.386 | 0.005 to ≤ 1.000; ≥ 0.000 |
| spectinomycin | 0.181 | 0.115 to 0.258 | 0.029 to ≤ 1.000; ≥ 0.000 |
| erythromycin | 0.153 | 0.103 to 0.204 | 0.000 to ≤ 1.000; ≥ 0.000 |
| tylosin | 0.151 | 0.100 to 0.201 | 0.000 to ≤ 1.000; ≥ 0.000 |
| tilmicosin | 0.149 | 0.098 to 0.198 | 0.000 to ≤ 1.000; ≥ 0.000 |
| doxycycline | 0.035 | 0.001 to 0.078 | 0.000 to ≤ 1.000; ≥ 0.000 |
| tetracycline | 0.035 | 0.005 to 0.073 | 0.000 to ≤ 1.000; ≥ 0.000 |

## 4. Evidence that survives re-reading

A p-value is a statement about one look at the data. A surveillance panel is looked at again every year, and a p-value recomputed each time loses its error control. The e-value is evidence on a scale made for that: this run's e-value is a statement about this collection, and a programme that adds an intake each year multiplies the e-value of each new intake, scored against the lineage rates learned from the earlier ones, into a running product whose error control holds at whatever intake it is read (sequential_e_process). The e-BH procedure controls the false-discovery rate across the panel whatever the dependence between traits. Larger is stronger; 1 is no evidence.

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

This collection is also divided into intakes, 25 of them read in the order of collection_year, from 1 to 177 isolates each, with 26 isolates carrying no intake and set aside. Each intake is scored against the lineage rates learned from the intakes before it, and the running product is the sequential e-value: it may be read after any intake without spending its error control. 2 of 13 traits are selected on the product, against 11 at this single look, which is the price of error control that holds at whatever intake the programme is read at.

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

Where a dilution was recorded, the lineage question is also read from the dilution itself rather than from the call derived from it. Every reading is placed by its position among the readings of its own stratum, a level of `testing_laboratory × isolation_country`: the share of those readings below it plus half the share at the same reading, with a censored reading placed by the nonparametric maximum-likelihood distribution of the readings. The share of these positions that lineage accounts for is estimated as the share of a call is, and its control shuffles the lineage labels only within strata. No distribution is assumed for the readings, a change of the MIC scale such as mg/L to log2 leaves the share unchanged, and a constant offset or a different panel between strata drops out. The share of the call and this share are read on independently retained readable and typed subsets; their denominators may differ.

The readings also bound the lineage share of the MIC ordering they were read from, the ordering of the latent MICs within a level of `testing_laboratory × isolation_country`: the rank intraclass correlation of the MIC itself. A reading says only that the MIC lies in its range, so the readings do not identify that share; they bound it. The lower bound is the smallest share that any arrangement of the latent MICs within their readings allows, the lower bound the readings establish, computed with a certified error; the upper end is the largest share, attained by an arrangement, or, marked ≤, a certified bound on it. Every share between them is allowed. The figure after the bounds is the one-sided 95 % lower confidence limit of the lower bound, and so of the share: the isolates of every lineage are split at random, one half chooses the lineage scores on which the lower bound rests and the other half tests them, and then the other way round; the tests of 25 such splits are averaged, and their average spread sets the limit. A finer panel never lowers the lower bound, and a call, one cut of the panel, never raises it. Resolution is the share of the ordering the panel resolves, one less the expected squared width of a reading on the scale of positions; it is the factor by which the share of the midpoint arrangement falls short of the plug-in share of the scores.

The share of the latent ordering is pooled over strata: a lineage placed high in one stratum and low in another does not count as ordered. For 3 agents the readings of some stratum on its own establish a larger lower bound than the pooled one, so lineage effects differ between strata: cefquinome (pooled 0.198; LGC Fordham / Canada 0.376); lincomycin (pooled 0.059; LGC Fordham / United Kingdom 0.215); penicillin (pooled 0.220; LGC Fordham / Canada 0.389).

**Table 5.** Lineage share of the MIC ordering per agent, with the 95 % interval for the lineages of this collection and the permutation p-value and the Benjamini-Yekutieli q-value across agents; selected is the Benjamini-Yekutieli selection at level 0.05. Control is the share on shuffled lineage labels. The bounds are those the readings place on the share of the latent MIC ordering, with the 95 % lower confidence limit of the lower bound. Readings on an end well are tied at the panel's edge and carry less of the ordering; their share is given beside the estimate. A p-value written with ≤ is the smallest the permutations can give. The last column gives the conclusion: established, lineage structure established, when the interval lies above zero and the agent is selected, or the lower confidence limit is above zero; not established; whole range, when the interval spans it; or not estimable. The conclusion says whether this collection shows lineage structure in the measurement, with the values behind it. It is not a statement about transmission or a resistance mechanism.

| Agent | Readings | Share | 95 % interval, represented lineages | p | q | Control | Selected | Latent ordering: bounds; 95 % lower limit | Resolution | End wells | Strata | Conclusion |
|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---|
| amoxicillin | 677 | 0.311 | 0.212 to 0.417 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.038 to ≤ 1.000; ≥ 0.001 | 0.272 | 89.1 % | 3 | established |
| cefquinome | 677 | 0.380 | 0.319 to 0.442 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.198 to ≤ 0.692; ≥ 0.079 | 0.889 | 6.9 % | 3 | established |
| ceftiofur | 677 | 0.485 | 0.424 to 0.551 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.255 to ≤ 0.845; ≥ 0.128 | 0.826 | 1.8 % | 3 | established |
| doxycycline | 677 | 0.139 | 0.089 to 0.196 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.031 to ≤ 0.594; ≥ 0.000 | 0.843 | 2.5 % | 3 | established |
| enrofloxacin | 677 | 0.214 | 0.157 to 0.281 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.030 to ≤ 0.856; ≥ 0.000 | 0.777 | 0.7 % | 3 | established |
| erythromycin | 677 | 0.106 | 0.063 to 0.155 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.011 to ≤ 0.517; ≥ 0.000 | 0.858 | 56.4 % | 3 | established |
| florfenicol | 677 | 0.228 | 0.168 to 0.292 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.001 to ≤ 1.000; ≥ 0.000 | 0.505 | 3.8 % | 3 | established |
| lincomycin | 677 | 0.239 | 0.172 to 0.307 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.059 to ≤ 0.798; ≥ 0.004 | 0.776 | 55.8 % | 3 | established |
| marbofloxacin | 677 | 0.158 | 0.105 to 0.220 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.011 to ≤ 0.747; ≥ 0.000 | 0.769 | 5.5 % | 3 | established |
| penicillin | 677 | 0.544 | 0.483 to 0.602 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.220 to ≤ 1.000; ≥ 0.115 | 0.585 | 73.0 % | 3 | established |
| spectinomycin | 677 | 0.192 | 0.133 to 0.256 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.055 to ≤ 0.713; ≥ 0.000 | 0.794 | 8.7 % | 3 | established |
| tetracycline | 677 | 0.107 | 0.060 to 0.160 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.010 to ≤ 0.612; ≥ 0.000 | 0.825 | 10.6 % | 3 | established |
| tiamulin | 677 | 0.462 | 0.392 to 0.543 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.351 to ≤ 0.637; ≥ 0.215 | 0.943 | 10.0 % | 3 | established |
| tilmicosin | 677 | 0.158 | 0.110 to 0.213 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.006 to ≤ 0.869; ≥ 0.000 | 0.762 | 2.4 % | 3 | established |
| trimethoprim | 677 | 0.247 | 0.188 to 0.313 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.181 to ≤ 0.386; ≥ 0.064 | 0.947 | 11.5 % | 3 | established |
| tylosin | 677 | 0.130 | 0.082 to 0.181 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.005 to ≤ 0.731; ≥ 0.000 | 0.787 | 52.4 % | 3 | established |

Resolution of the selection: with 999 permutations one agent alone cannot be selected among 16 at alpha = 0.05; at least 1081 permutations would let a single very small p-value pass the step-up.

**Table 5b.** Lineage structure behind each share of the MIC ordering: the lineages with at least two readings on which the share is scored, the readings of singleton lineages set aside, the effective number of scored lineages (inverse of the sum of squared lineage shares) and the share of readings in lineages with at least two members (support). Few effective lineages mean the share rests on a few lineages.

| Agent | Lineages scored | Singletons set aside | Effective | Support |
|---|---:|---:|---:|---:|
| amoxicillin | 27 | 3 | 10.1 | 99.6 % |
| cefquinome | 27 | 3 | 10.1 | 99.6 % |
| ceftiofur | 27 | 3 | 10.1 | 99.6 % |
| doxycycline | 27 | 3 | 10.1 | 99.6 % |
| enrofloxacin | 27 | 3 | 10.1 | 99.6 % |
| erythromycin | 27 | 3 | 10.1 | 99.6 % |
| florfenicol | 27 | 3 | 10.1 | 99.6 % |
| lincomycin | 27 | 3 | 10.1 | 99.6 % |
| marbofloxacin | 27 | 3 | 10.1 | 99.6 % |
| penicillin | 27 | 3 | 10.1 | 99.6 % |
| spectinomycin | 27 | 3 | 10.1 | 99.6 % |
| tetracycline | 27 | 3 | 10.1 | 99.6 % |
| tiamulin | 27 | 3 | 10.1 | 99.6 % |
| tilmicosin | 27 | 3 | 10.1 | 99.6 % |
| trimethoprim | 27 | 3 | 10.1 | 99.6 % |
| tylosin | 27 | 3 | 10.1 | 99.6 % |

*Figure 5 (drawn on the page).* Lineage share of the MIC ordering per agent, 16 largest of 16: the point estimate with its 95 % interval; agents selected across the panel are drawn in blue. The hollow diamond is the share of the binary call from Table 1 for the same agent, a different quantity whose retained isolates may differ, placed here so that the two readings can be seen side by side. The figure in parentheses at the right is the share of readings on an end well.

Each laboratory is read on the wells that laboratory tested: an end-well reading is censored at the panel edge of its own laboratory, not at the widest edge in the collection. The table below gives, per agent and laboratory, the panel the readings were taken to come from.

**Table 5a.** Panel geometry per agent and testing laboratory: the wells taken as tested, their range, whether they form a doubling series, and the shares of readings on the lowest and the highest well.

| Agent | Laboratory | Wells | Range | Doubling | On lowest well | On highest well |
|---|---|---:|---:|---|---:|---:|
| amoxicillin | LGC Fordham | 8 | 0.0312 to 4.0000 | yes | 87.7 % | 0.5 % |
| amoxicillin | OUCRU Ho Chi Minh City | 2 | 0.0156 to 0.0312 | yes | 98.0 % | 2.0 % |
| cefquinome | LGC Fordham | 11 | 0.0020 to 2.0000 | yes | 0.8 % | 0.8 % |
| cefquinome | OUCRU Ho Chi Minh City | 3 | 0.0078 to 0.0312 | yes | 4.1 % | 71.4 % |
| ceftiofur | LGC Fordham | 10 | 0.0312 to 16.0000 | yes | 0.5 % | 1.0 % |
| ceftiofur | OUCRU Ho Chi Minh City | 4 | 0.0625 to 0.5000 | yes | 4.1 % | 2.0 % |
| doxycycline | LGC Fordham | 12 | 0.0312 to 64.0000 | yes | 0.2 % | 0.2 % |
| doxycycline | OUCRU Ho Chi Minh City | 9 | 0.0312 to 8.0000 | yes | 8.2 % | 22.4 % |
| enrofloxacin | LGC Fordham | 11 | 0.0078 to 8.0000 | yes | 0.2 % | 0.3 % |
| enrofloxacin | OUCRU Ho Chi Minh City | 5 | 0.2500 to 4.0000 | yes | 2.0 % | 2.0 % |
| erythromycin | LGC Fordham | 12 | 0.0156 to 32.0000 | yes | 8.9 % | 46.3 % |
| erythromycin | OUCRU Ho Chi Minh City | 12 | 0.0156 to 32.0000 | yes | 57.1 % | 14.3 % |
| florfenicol | LGC Fordham | 4 | 0.5000 to 4.0000 | yes | 0.6 % | 3.2 % |
| florfenicol | OUCRU Ho Chi Minh City | 5 | 0.2500 to 4.0000 | yes | 2.0 % | 2.0 % |
| lincomycin | LGC Fordham | 12 | 0.0625 to 128.0000 | yes | 0.8 % | 57.5 % |
| lincomycin | OUCRU Ho Chi Minh City | 12 | 0.0625 to 128.0000 | yes | 2.0 % | 22.4 % |
| marbofloxacin | LGC Fordham | 11 | 0.0156 to 16.0000 | yes | 4.3 % | 0.3 % |
| marbofloxacin | OUCRU Ho Chi Minh City | 3 | 0.2500 to 1.0000 | yes | 6.1 % | 10.2 % |
| penicillin | LGC Fordham | 10 | 0.0312 to 16.0000 | yes | 70.7 % | 0.2 % |
| penicillin | OUCRU Ho Chi Minh City | 5 | 0.0312 to 0.5000 | yes | 95.9 % | 4.1 % |
| spectinomycin | LGC Fordham | 10 | 1.0000 to 512.0000 | yes | 0.2 % | 8.3 % |
| spectinomycin | OUCRU Ho Chi Minh City | 8 | 2.0000 to 256.0000 | yes | 2.0 % | 10.2 % |
| tetracycline | LGC Fordham | 12 | 0.0625 to 128.0000 | yes | 0.2 % | 8.6 % |
| tetracycline | OUCRU Ho Chi Minh City | 9 | 0.5000 to 128.0000 | yes | 10.2 % | 24.5 % |
| tiamulin | LGC Fordham | 10 | 0.1250 to 64.0000 | yes | 2.1 % | 8.3 % |
| tiamulin | OUCRU Ho Chi Minh City | 7 | 0.0312 to 2.0000 | yes | 4.1 % | 2.0 % |
| tilmicosin | LGC Fordham | 11 | 0.2500 to 256.0000 | yes | 0.2 % | 0.2 % |
| tilmicosin | OUCRU Ho Chi Minh City | 12 | 0.5000 to 1024.0000 | yes | 4.1 % | 24.5 % |
| trimethoprim | LGC Fordham | 12 | 0.0156 to 32.0000 | yes | 5.4 % | 5.3 % |
| trimethoprim | OUCRU Ho Chi Minh City | 8 | 0.0625 to 8.0000 | yes | 20.4 % | 2.0 % |
| tylosin | LGC Fordham | 12 | 0.1250 to 256.0000 | yes | 0.2 % | 54.1 % |
| tylosin | OUCRU Ho Chi Minh City | 12 | 0.1250 to 256.0000 | yes | 10.2 % | 18.4 % |

*Figure 6 (drawn on the page).* Readings per lineage and dilution interval for 6 of 16 agents, the 14 largest of 30 lineages: each cell counts the isolates of one lineage whose reading fell in one interval of the panel, darker for more. Open intervals at the panel edges are the censored readings. A lineage whose readings sit in one or two adjacent cells contributes to the lineage share; readings spread along a row do not.

## 6. The same analyses within each stratum

The run was repeated on the isolates of each level of `isolation_country` (Canada (n = 205), United Kingdom (n = 423) and Vietnam (n = 49); 0 isolates without a level were left out). Each stratum is an analysis of its own, on its own lineages and its own panel, and the shares below are not adjusted for one another. A share that holds within every stratum is a property of the lineages; a share seen only in the pooled run and in no stratum is accounted for by the differences between strata.

**Table 5d.** Lineage share per stratum and agent: the share of the binary call with its interval and p-value, and the lineage share of the MIC ordering with its interval and p-value where a MIC table was supplied; a share not computed carries its reason. strata_results.csv holds these and the other analyses of every stratum.

| Stratum | Agent | n | Call share | 95 % interval | p | MIC ordering share (95 % interval) | p |
|---|---|---:|---:|---:|---:|---:|---:|
| Canada | amoxicillin | 205 | 0.284 | 0.165 to 0.459 | ≤ 0.001 | 0.439 (0.289 to 0.601) | ≤ 0.001 |
| Canada | cefquinome | 205 | 0.223 | 0.055 to 0.442 | ≤ 0.001 | 0.528 (0.393 to 0.670) | ≤ 0.001 |
| Canada | ceftiofur | 205 | 0.519 | 0.391 to 0.652 | ≤ 0.001 | 0.546 (0.418 to 0.687) | ≤ 0.001 |
| Canada | doxycycline | 205 | -0.011 | 0.000 to 0.165 | 0.485 | 0.218 (0.099 to 0.354) | ≤ 0.001 |
| Canada | enrofloxacin |  | not computed |  |  | 0.195 (0.081 to 0.356) | ≤ 0.001 |
| Canada | erythromycin | 205 | 0.206 | 0.089 to 0.328 | ≤ 0.001 | 0.068 (0.000 to 0.181) | 0.025 |
| Canada | florfenicol |  | not computed |  |  | 0.149 (0.007 to 0.366) | ≤ 0.001 |
| Canada | lincomycin | 205 | 0.206 | 0.061 to 0.370 | ≤ 0.001 | 0.190 (0.048 to 0.343) | ≤ 0.001 |
| Canada | marbofloxacin |  | not computed |  |  | 0.184 (0.072 to 0.341) | ≤ 0.001 |
| Canada | penicillin | 205 | 0.585 | 0.476 to 0.702 | ≤ 0.001 | 0.617 (0.505 to 0.737) | ≤ 0.001 |
| Canada | spectinomycin | 205 | 0.091 | 0.000 to 0.239 | 0.004 | 0.158 (0.044 to 0.297) | ≤ 0.001 |
| Canada | tetracycline | 205 | 0.009 | 0.000 to 0.208 | 0.701 | 0.159 (0.047 to 0.300) | ≤ 0.001 |
| Canada | tiamulin | 205 | 0.454 | 0.307 to 0.610 | ≤ 0.001 | 0.436 (0.285 to 0.619) | ≤ 0.001 |
| Canada | tilmicosin | 205 | 0.191 | 0.078 to 0.302 | ≤ 0.001 | 0.204 (0.070 to 0.344) | ≤ 0.001 |
| Canada | trimethoprim | 205 | 0.360 | 0.221 to 0.513 | ≤ 0.001 | 0.341 (0.225 to 0.467) | ≤ 0.001 |
| Canada | tylosin | 205 | 0.193 | 0.078 to 0.306 | ≤ 0.001 | 0.169 (0.035 to 0.304) | ≤ 0.001 |
| United Kingdom | amoxicillin | 423 | 0.279 | 0.000 to 0.685 | ≤ 0.001 | 0.256 (0.000 to 0.569) | ≤ 0.001 |
| United Kingdom | cefquinome | 423 | 0.336 | 0.000 to 0.898 | ≤ 0.001 | 0.350 (0.273 to 0.432) | ≤ 0.001 |
| United Kingdom | ceftiofur | 423 | 0.546 | 0.454 to 0.634 | ≤ 0.001 | 0.514 (0.431 to 0.605) | ≤ 0.001 |
| United Kingdom | doxycycline | 423 | 0.118 | 0.067 to 0.178 | ≤ 0.001 | 0.154 (0.086 to 0.224) | ≤ 0.001 |
| United Kingdom | enrofloxacin |  | not computed |  |  | 0.339 (0.244 to 0.441) | ≤ 0.001 |
| United Kingdom | erythromycin | 423 | 0.216 | 0.147 to 0.288 | ≤ 0.001 | 0.162 (0.098 to 0.235) | ≤ 0.001 |
| United Kingdom | florfenicol |  | not computed |  |  | 0.339 (0.244 to 0.458) | ≤ 0.001 |
| United Kingdom | lincomycin | 423 | 0.343 | 0.273 to 0.420 | ≤ 0.001 | 0.390 (0.302 to 0.475) | ≤ 0.001 |
| United Kingdom | marbofloxacin |  | not computed |  |  | 0.224 (0.149 to 0.318) | ≤ 0.001 |
| United Kingdom | penicillin | 423 | 0.503 | 0.398 to 0.612 | ≤ 0.001 | 0.558 (0.456 to 0.661) | ≤ 0.001 |
| United Kingdom | spectinomycin | 423 | 0.424 | 0.285 to 0.604 | ≤ 0.001 | 0.293 (0.206 to 0.389) | ≤ 0.001 |
| United Kingdom | tetracycline | 423 | 0.118 | 0.067 to 0.178 | ≤ 0.001 | 0.160 (0.091 to 0.246) | ≤ 0.001 |
| United Kingdom | tiamulin | 423 | 0.403 | 0.313 to 0.499 | ≤ 0.001 | 0.554 (0.467 to 0.651) | ≤ 0.001 |
| United Kingdom | tilmicosin | 423 | 0.220 | 0.154 to 0.287 | ≤ 0.001 | 0.221 (0.148 to 0.302) | ≤ 0.001 |
| United Kingdom | trimethoprim | 423 | 0.331 | 0.246 to 0.423 | ≤ 0.001 | 0.325 (0.250 to 0.400) | ≤ 0.001 |
| United Kingdom | tylosin | 423 | 0.234 | 0.165 to 0.305 | ≤ 0.001 | 0.188 (0.117 to 0.267) | ≤ 0.001 |
| Vietnam | amoxicillin | 49 | not computed |  |  | not estimable |  |
| Vietnam | cefquinome | 49 | not computed |  |  | not estimable |  |
| Vietnam | ceftiofur | 49 | not computed |  |  | not estimable |  |
| Vietnam | doxycycline | 49 | not computed |  |  | not estimable |  |
| Vietnam | enrofloxacin |  | not computed |  |  | not estimable |  |
| Vietnam | erythromycin | 49 | not computed |  |  | not estimable |  |
| Vietnam | florfenicol |  | not computed |  |  | not estimable |  |
| Vietnam | lincomycin | 49 | not computed |  |  | not estimable |  |
| Vietnam | marbofloxacin |  | not computed |  |  | not estimable |  |
| Vietnam | penicillin | 49 | not computed |  |  | not estimable |  |
| Vietnam | spectinomycin | 49 | not computed |  |  | not estimable |  |
| Vietnam | tetracycline | 49 | not computed |  |  | not estimable |  |
| Vietnam | tiamulin | 49 | not computed |  |  | not estimable |  |
| Vietnam | tilmicosin | 49 | not computed |  |  | not estimable |  |
| Vietnam | trimethoprim | 49 | not computed |  |  | not estimable |  |
| Vietnam | tylosin | 49 | not computed |  |  | not estimable |  |

## 7. A change of lineages or a change within them

A prevalence difference between two collections can be split into a change in lineage composition and a change in rate within lineages. That decomposition needs two collections and was not run here: this record holds one.

## 8. Provenance and terms

- **Software:** amr-clonalshare 1.0.0, record schema 1.0
- **Seed:** 42
- **Configuration:** `sha256:7675bb625cb2`
- **Record:** `sha256:3ee3d16084cc3520502340537a20203de1e6d99fcbd49914274e3e1ee1920b46`

> **Terms used in this report**
>
> **Share.** How much of the variation of the call between isolates is associated with their recorded lineage. Near 1: a strong association on the measured scale. Near 0: little association on that scale. This does not identify transmission or its mechanism.
>
> **MIC ordering share.** The share read from the dilution: how much of the ordering of the readings within a laboratory (or another declared stratum) is associated with lineage. It assumes no distribution for the readings.
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

amr-clonalshare 1.0.0 · record schema 1.0 · seed 42 · configuration sha256:7675bb625cb2
Cite this run as: “amr-clonalshare 1.0.0, run 7DEE408A, record sha256:3ee3d16084cc.”
This report supersedes any earlier report bearing the same run identifier. It is regenerated from the record and holds no value that the record does not. The symbols carry the same wording in every run of this software; no symbol against a value means only that none of the listed conditions fired.
