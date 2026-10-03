# amr-clonalshare run report

Every number below is read from `clonal_share_result.json`, written by the same run; none is recomputed here.

- **Run:** F3CC72CB
- **Issued:** 2026-10-03T06:12Z
- **Software:** amr-clonalshare 1.0.0
- **Isolates:** 677
- **Lineage:** mlst
- **Seed:** 42
- **Configuration:** sha256:a9a983a1ba3b
- **Record digest:** sha256:ae9f745ec7a07ea3a49aef8a9bf81adbc5133df5e69d4fb11cdec5a827c7d34b

**16 antimicrobials; 139 analysis outcomes.** 112 computed, 1 full ranges, 26 unavailable. Each result retains its own target and data subset.

## 1. Measurement summary

**Quantity measured.** The share of the variation in the recorded binary outcome across this collection that lineage membership accounts for, read from a lineage label and an interpreted result and scored on isolates the estimator did not see.

- **Collection:** 677 isolates in 108 lineages, typed by `mlst`
- **Antimicrobials read:** 16
- **Membership shares estimable:** 13 of 13

**Strata of the calls.** The analyses of every call were read within the 3 levels of `isolation_country`: lineage labels are exchanged only within a level, the interval draws each isolate's level with its call, and the bounds a call places on the latent ordering are the bounds within levels, so a shift between levels is not read as a lineage difference. The e-values take no strata.

**Reporting conditions.** 112 computed, 1 full ranges, 26 unavailable. Each result retains its own target and data subset. An unavailable value failed a declared reporting condition of its estimator; the recorded reason distinguishes input limitations from numerical failure.

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
| amoxicillin | collection_membership | 458 | computed | observed and typed subset |
| amoxicillin | call_latent_bounds | 458 | computed | observed and typed subset |
| amoxicillin | call_latent_lower_limit | 458 | computed | observed and typed subset |
| cefquinome | collection_membership | 458 | full_range | observed and typed subset |
| cefquinome | call_latent_bounds | 458 | computed | observed and typed subset |
| cefquinome | call_latent_lower_limit | 458 | computed | observed and typed subset |
| ceftiofur | collection_membership | 458 | computed | observed and typed subset |
| ceftiofur | call_latent_bounds | 458 | computed | observed and typed subset |
| ceftiofur | call_latent_lower_limit | 458 | computed | observed and typed subset |
| doxycycline | collection_membership | 458 | computed | observed and typed subset |
| doxycycline | call_latent_bounds | 458 | computed | observed and typed subset |
| doxycycline | call_latent_lower_limit | 458 | computed | observed and typed subset |
| erythromycin | collection_membership | 458 | computed | observed and typed subset |
| erythromycin | call_latent_bounds | 458 | computed | observed and typed subset |
| erythromycin | call_latent_lower_limit | 458 | computed | observed and typed subset |
| lincomycin | collection_membership | 458 | computed | observed and typed subset |
| lincomycin | call_latent_bounds | 458 | computed | observed and typed subset |
| lincomycin | call_latent_lower_limit | 458 | computed | observed and typed subset |
| penicillin | collection_membership | 458 | computed | observed and typed subset |
| penicillin | call_latent_bounds | 458 | computed | observed and typed subset |
| penicillin | call_latent_lower_limit | 458 | computed | observed and typed subset |
| spectinomycin | collection_membership | 458 | computed | observed and typed subset |
| spectinomycin | call_latent_bounds | 458 | computed | observed and typed subset |
| spectinomycin | call_latent_lower_limit | 458 | computed | observed and typed subset |
| tetracycline | collection_membership | 458 | computed | observed and typed subset |
| tetracycline | call_latent_bounds | 458 | computed | observed and typed subset |
| tetracycline | call_latent_lower_limit | 458 | computed | observed and typed subset |
| tiamulin | collection_membership | 458 | computed | observed and typed subset |
| tiamulin | call_latent_bounds | 458 | computed | observed and typed subset |
| tiamulin | call_latent_lower_limit | 458 | computed | observed and typed subset |
| tilmicosin | collection_membership | 458 | computed | observed and typed subset |
| tilmicosin | call_latent_bounds | 458 | computed | observed and typed subset |
| tilmicosin | call_latent_lower_limit | 458 | computed | observed and typed subset |
| trimethoprim | collection_membership | 458 | computed | observed and typed subset |
| trimethoprim | call_latent_bounds | 458 | computed | observed and typed subset |
| trimethoprim | call_latent_lower_limit | 458 | computed | observed and typed subset |
| tylosin | collection_membership | 458 | computed | observed and typed subset |
| tylosin | call_latent_bounds | 458 | computed | observed and typed subset |
| tylosin | call_latent_lower_limit | 458 | computed | observed and typed subset |
| amoxicillin | mic_order | 458 | computed | readable and typed MIC subset |
| amoxicillin | mic_latent_bounds | 458 | computed | readable and typed MIC subset |
| amoxicillin | mic_latent_lower_limit | 458 | computed | readable and typed MIC subset |
| cefquinome | mic_order | 458 | computed | readable and typed MIC subset |
| cefquinome | mic_latent_bounds | 458 | computed | readable and typed MIC subset |
| cefquinome | mic_latent_lower_limit | 458 | computed | readable and typed MIC subset |
| ceftiofur | mic_order | 458 | computed | readable and typed MIC subset |
| ceftiofur | mic_latent_bounds | 458 | computed | readable and typed MIC subset |
| ceftiofur | mic_latent_lower_limit | 458 | computed | readable and typed MIC subset |
| doxycycline | mic_order | 458 | computed | readable and typed MIC subset |
| doxycycline | mic_latent_bounds | 458 | computed | readable and typed MIC subset |
| doxycycline | mic_latent_lower_limit | 458 | computed | readable and typed MIC subset |
| enrofloxacin | mic_order | 458 | computed | readable and typed MIC subset |
| enrofloxacin | mic_latent_bounds | 458 | computed | readable and typed MIC subset |
| enrofloxacin | mic_latent_lower_limit | 458 | computed | readable and typed MIC subset |
| erythromycin | mic_order | 458 | computed | readable and typed MIC subset |
| erythromycin | mic_latent_bounds | 458 | computed | readable and typed MIC subset |
| erythromycin | mic_latent_lower_limit | 458 | computed | readable and typed MIC subset |
| florfenicol | mic_order | 458 | computed | readable and typed MIC subset |
| florfenicol | mic_latent_bounds | 458 | computed | readable and typed MIC subset |
| florfenicol | mic_latent_lower_limit | 458 | computed | readable and typed MIC subset |
| lincomycin | mic_order | 458 | computed | readable and typed MIC subset |
| lincomycin | mic_latent_bounds | 458 | computed | readable and typed MIC subset |
| lincomycin | mic_latent_lower_limit | 458 | computed | readable and typed MIC subset |
| marbofloxacin | mic_order | 458 | computed | readable and typed MIC subset |
| marbofloxacin | mic_latent_bounds | 458 | computed | readable and typed MIC subset |
| marbofloxacin | mic_latent_lower_limit | 458 | computed | readable and typed MIC subset |
| penicillin | mic_order | 458 | computed | readable and typed MIC subset |
| penicillin | mic_latent_bounds | 458 | computed | readable and typed MIC subset |
| penicillin | mic_latent_lower_limit | 458 | computed | readable and typed MIC subset |
| spectinomycin | mic_order | 458 | computed | readable and typed MIC subset |
| spectinomycin | mic_latent_bounds | 458 | computed | readable and typed MIC subset |
| spectinomycin | mic_latent_lower_limit | 458 | computed | readable and typed MIC subset |
| tetracycline | mic_order | 458 | computed | readable and typed MIC subset |
| tetracycline | mic_latent_bounds | 458 | computed | readable and typed MIC subset |
| tetracycline | mic_latent_lower_limit | 458 | computed | readable and typed MIC subset |
| tiamulin | mic_order | 458 | computed | readable and typed MIC subset |
| tiamulin | mic_latent_bounds | 458 | computed | readable and typed MIC subset |
| tiamulin | mic_latent_lower_limit | 458 | computed | readable and typed MIC subset |
| tilmicosin | mic_order | 458 | computed | readable and typed MIC subset |
| tilmicosin | mic_latent_bounds | 458 | computed | readable and typed MIC subset |
| tilmicosin | mic_latent_lower_limit | 458 | computed | readable and typed MIC subset |
| trimethoprim | mic_order | 458 | computed | readable and typed MIC subset |
| trimethoprim | mic_latent_bounds | 458 | computed | readable and typed MIC subset |
| trimethoprim | mic_latent_lower_limit | 458 | computed | readable and typed MIC subset |
| tylosin | mic_order | 458 | computed | readable and typed MIC subset |
| tylosin | mic_latent_bounds | 458 | computed | readable and typed MIC subset |
| tylosin | mic_latent_lower_limit | 458 | computed | readable and typed MIC subset |
| amoxicillin | lineage_evidence | 458 | computed | observed and typed subset |
| cefquinome | lineage_evidence | 458 | computed | observed and typed subset |
| ceftiofur | lineage_evidence | 458 | computed | observed and typed subset |
| doxycycline | lineage_evidence | 458 | computed | observed and typed subset |
| erythromycin | lineage_evidence | 458 | computed | observed and typed subset |
| lincomycin | lineage_evidence | 458 | computed | observed and typed subset |
| penicillin | lineage_evidence | 458 | computed | observed and typed subset |
| spectinomycin | lineage_evidence | 458 | computed | observed and typed subset |
| tetracycline | lineage_evidence | 458 | computed | observed and typed subset |
| tiamulin | lineage_evidence | 458 | computed | observed and typed subset |
| tilmicosin | lineage_evidence | 458 | computed | observed and typed subset |
| trimethoprim | lineage_evidence | 458 | computed | observed and typed subset |
| tylosin | lineage_evidence | 458 | computed | observed and typed subset |
| amoxicillin | decomposition_composition | 409 | unavailable | observed and typed subset |
| amoxicillin | decomposition_within_lineage | 409 | unavailable | observed and typed subset |
| cefquinome | decomposition_composition | 409 | unavailable | observed and typed subset |
| cefquinome | decomposition_within_lineage | 409 | unavailable | observed and typed subset |
| ceftiofur | decomposition_composition | 409 | unavailable | observed and typed subset |
| ceftiofur | decomposition_within_lineage | 409 | unavailable | observed and typed subset |
| doxycycline | decomposition_composition | 409 | unavailable | observed and typed subset |
| doxycycline | decomposition_within_lineage | 409 | unavailable | observed and typed subset |
| erythromycin | decomposition_composition | 409 | unavailable | observed and typed subset |
| erythromycin | decomposition_within_lineage | 409 | unavailable | observed and typed subset |
| lincomycin | decomposition_composition | 409 | unavailable | observed and typed subset |
| lincomycin | decomposition_within_lineage | 409 | unavailable | observed and typed subset |
| penicillin | decomposition_composition | 409 | unavailable | observed and typed subset |
| penicillin | decomposition_within_lineage | 409 | unavailable | observed and typed subset |
| spectinomycin | decomposition_composition | 409 | unavailable | observed and typed subset |
| spectinomycin | decomposition_within_lineage | 409 | unavailable | observed and typed subset |
| tetracycline | decomposition_composition | 409 | unavailable | observed and typed subset |
| tetracycline | decomposition_within_lineage | 409 | unavailable | observed and typed subset |
| tiamulin | decomposition_composition | 409 | unavailable | observed and typed subset |
| tiamulin | decomposition_within_lineage | 409 | unavailable | observed and typed subset |
| tilmicosin | decomposition_composition | 409 | unavailable | observed and typed subset |
| tilmicosin | decomposition_within_lineage | 409 | unavailable | observed and typed subset |
| trimethoprim | decomposition_composition | 409 | unavailable | observed and typed subset |
| trimethoprim | decomposition_within_lineage | 409 | unavailable | observed and typed subset |
| tylosin | decomposition_composition | 409 | unavailable | observed and typed subset |
| tylosin | decomposition_within_lineage | 409 | unavailable | observed and typed subset |

**Table 0b.** Reasons and next checks for results requiring interpretation. These checks do not change the recorded status or justify changing reporting conditions.

| Agent | Analysis / status | Reason | Next check |
|---|---|---|---|
| amoxicillin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.833 | Read the interval kind, assumptions and retained collection before comparing results. |
| cefquinome | collection_membership / full_range | The completed interval or bounds span the entire admissible range [0, 1] | Report the full range as uninformative; inspect repeated-lineage support and outcome variation. |
| cefquinome | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.840 | Read the interval kind, assumptions and retained collection before comparing results. |
| ceftiofur | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.839 | Read the interval kind, assumptions and retained collection before comparing results. |
| doxycycline | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.767 | Read the interval kind, assumptions and retained collection before comparing results. |
| erythromycin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.718 | Read the interval kind, assumptions and retained collection before comparing results. |
| lincomycin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.714 | Read the interval kind, assumptions and retained collection before comparing results. |
| penicillin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.839 | Read the interval kind, assumptions and retained collection before comparing results. |
| spectinomycin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.834 | Read the interval kind, assumptions and retained collection before comparing results. |
| tetracycline | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.767 | Read the interval kind, assumptions and retained collection before comparing results. |
| tiamulin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.837 | Read the interval kind, assumptions and retained collection before comparing results. |
| tilmicosin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.705 | Read the interval kind, assumptions and retained collection before comparing results. |
| trimethoprim | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.810 | Read the interval kind, assumptions and retained collection before comparing results. |
| tylosin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.713 | Read the interval kind, assumptions and retained collection before comparing results. |
| amoxicillin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.825 | Read the interval kind, assumptions and retained collection before comparing results. |
| cefquinome | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.648 | Read the interval kind, assumptions and retained collection before comparing results. |
| ceftiofur | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.706 | Read the interval kind, assumptions and retained collection before comparing results. |
| doxycycline | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.494 | Read the interval kind, assumptions and retained collection before comparing results. |
| enrofloxacin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.647 | Read the interval kind, assumptions and retained collection before comparing results. |
| erythromycin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.494 | Read the interval kind, assumptions and retained collection before comparing results. |
| florfenicol | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.783 | Read the interval kind, assumptions and retained collection before comparing results. |
| lincomycin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.626 | Read the interval kind, assumptions and retained collection before comparing results. |
| marbofloxacin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.554 | Read the interval kind, assumptions and retained collection before comparing results. |
| penicillin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.833 | Read the interval kind, assumptions and retained collection before comparing results. |
| spectinomycin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.645 | Read the interval kind, assumptions and retained collection before comparing results. |
| tetracycline | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.514 | Read the interval kind, assumptions and retained collection before comparing results. |
| tiamulin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.609 | Read the interval kind, assumptions and retained collection before comparing results. |
| tilmicosin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.683 | Read the interval kind, assumptions and retained collection before comparing results. |
| trimethoprim | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.531 | Read the interval kind, assumptions and retained collection before comparing results. |
| tylosin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.590 | Read the interval kind, assumptions and retained collection before comparing results. |
| amoxicillin | decomposition_composition / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| amoxicillin | decomposition_within_lineage / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| cefquinome | decomposition_composition / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| cefquinome | decomposition_within_lineage / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| ceftiofur | decomposition_composition / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| ceftiofur | decomposition_within_lineage / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| doxycycline | decomposition_composition / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| doxycycline | decomposition_within_lineage / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| erythromycin | decomposition_composition / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| erythromycin | decomposition_within_lineage / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| lincomycin | decomposition_composition / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| lincomycin | decomposition_within_lineage / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| penicillin | decomposition_composition / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| penicillin | decomposition_within_lineage / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| spectinomycin | decomposition_composition / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| spectinomycin | decomposition_within_lineage / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| tetracycline | decomposition_composition / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| tetracycline | decomposition_within_lineage / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| tiamulin | decomposition_composition / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| tiamulin | decomposition_within_lineage / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| tilmicosin | decomposition_composition / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| tilmicosin | decomposition_within_lineage / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| trimethoprim | decomposition_composition / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| trimethoprim | decomposition_within_lineage / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| tylosin | decomposition_composition / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |
| tylosin | decomposition_within_lineage / unavailable | ['shared lineage support below the threshold'] | Check both recorded collections, readable outcomes and shared lineage support; retain the reported subset scope. |

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

**Comparison scope.** Decomposition components describe the observed, labelled subsets. Equal typing coverage or a nonsignificant missingness test does not establish representativeness. Full observed prevalence and subset prevalence are recorded separately; missing data prevent automatic collection-wide generalization.

**Table 1.** Every trait, ordered by share, with the 95 % interval, the permutation p-value and the Benjamini-Yekutieli q-value across traits; selected is the Benjamini-Yekutieli selection at a false-discovery level of 0.05, which holds whatever the dependence between traits. Control is the share on shuffled lineage labels. A p-value written with ≤ is the smallest the permutations can give. The last column gives the conclusion of the results table of the form: established, lineage structure established, when the interval lies above zero and the trait is selected, or the lower confidence limit is above zero; not established; whole range, when the interval spans it; or not estimable. The conclusion says whether this collection shows lineage structure in the measurement, with the values behind it. It is not a statement about transmission or a resistance mechanism.

| Trait | Share | 95 % interval | p | q | Control | e-value | Selected | Conclusion |
|---|---:|---:|---:|---:|---:|---:|---|---|
| ceftiofur | 0.593 | 0.379 to 0.871 | ≤ 0.001 | 0.004 | 0.00 | 5.6 × 10⁶ | yes | established |
| penicillin | 0.535 | 0.317 to 0.822 | ≤ 0.001 | 0.004 | 0.00 | 5.7 × 10⁴ | yes | established |
| trimethoprim | 0.406 | 0.274 to 0.547 | ≤ 0.001 | 0.004 | 0.00 | 1.3 × 10⁶ | yes | established |
| cefquinome | 0.392 | 0.000 to 1.000 | ≤ 0.001 | 0.004 | 0.00 | 58.9 | yes | whole range |
| tiamulin | 0.385 | 0.111 to 0.770 | ≤ 0.001 | 0.004 | 0.00 | 266.1 | yes | established |
| erythromycin | 0.287 | 0.200 to 0.371 | ≤ 0.001 | 0.004 | 0.00 | 1.2 × 10⁴ | yes | established |
| lincomycin | 0.277 | 0.196 to 0.358 | ≤ 0.001 | 0.004 | 0.00 | 1.5 × 10³ | yes | established |
| tylosin | 0.276 | 0.191 to 0.356 | ≤ 0.001 | 0.004 | 0.00 | 7.4 × 10³ | yes | established |
| tilmicosin | 0.257 | 0.169 to 0.346 | ≤ 0.001 | 0.004 | 0.00 | 2.3 × 10³ | yes | established |
| doxycycline | 0.254 | 0.151 to 0.364 | ≤ 0.001 | 0.004 | 0.00 | 1.5 × 10³ | yes | established |
| tetracycline | 0.254 | 0.151 to 0.364 | ≤ 0.001 | 0.004 | 0.00 | 1.5 × 10³ | yes | established |
| amoxicillin | 0.110 | 0.000 to 0.580 | 0.010 | 0.034 | 0.00 | 3.4 | yes | not established |
| spectinomycin | 0.072 | 0.000 to 0.235 | 0.016 | 0.051 | 0.00 | 0.8 | no | not established |

**Table 1b.** Lineage structure behind each share: the lineages with at least two isolates on which the share is scored, the isolates of singleton lineages set aside, the effective number of scored lineages (inverse of the sum of squared lineage shares) and the share of isolates in lineages with at least two members (support). Few effective lineages mean the share rests on a few lineages.

| Trait | Lineages scored | Singletons set aside | Effective | Support |
|---|---:|---:|---:|---:|
| ceftiofur | 43 | 65 | 6.0 | 85.8 % |
| penicillin | 43 | 65 | 6.0 | 85.8 % |
| trimethoprim | 43 | 65 | 6.0 | 85.8 % |
| cefquinome | 43 | 65 | 6.0 | 85.8 % |
| tiamulin | 43 | 65 | 6.0 | 85.8 % |
| erythromycin | 43 | 65 | 6.0 | 85.8 % |
| lincomycin | 43 | 65 | 6.0 | 85.8 % |
| tylosin | 43 | 65 | 6.0 | 85.8 % |
| tilmicosin | 43 | 65 | 6.0 | 85.8 % |
| doxycycline | 43 | 65 | 6.0 | 85.8 % |
| tetracycline | 43 | 65 | 6.0 | 85.8 % |
| amoxicillin | 43 | 65 | 6.0 | 85.8 % |
| spectinomycin | 43 | 65 | 6.0 | 85.8 % |

> **Interpretation (generated from Table 1 by fixed rules)**
>
> For ceftiofur and penicillin, lineage membership is associated with the recorded outcome: the estimated share is at or above one half in this collection at the chosen typing resolution.
>
> For trimethoprim, tiamulin, erythromycin, lincomycin, tylosin, tilmicosin, doxycycline and tetracycline, a lineage association is detected and the estimated share is below one half. Most variation is not explained by the lineage labels at this typing resolution.
>
> For cefquinome and amoxicillin, the panel selection finds a lineage effect while the interval for its size still reaches zero: there is evidence that positive outcomes are not spread evenly across the lineages, and this collection is too small, or too uneven across its lineages, to say how much of it the lineages account for. Read the selection as the finding and the share as not yet resolved.
>
> For spectinomycin, no lineage effect is distinguishable from none in this collection at the chosen resolution. This does not establish independence from lineage or identify the reason for a prevalence change.
>
> These readings describe association in the sampled collection. The analysis does not identify transmission, horizontal transfer, selection, or the effect of an intervention.
>

## 2. Admissibility of the input

- **Lineage column:** `mlst`
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

Lineage groups in the input: 108, of which 65 hold a single isolate. Support, the share of isolates in lineages of at least two, is 85.8 %; the share is scored on those isolates and the singletons are set aside. At least two lineages hold two or more isolates.

*Figure 1 (drawn on the page).* Isolates per lineage, largest first, the 40 largest of 108. A lineage of one isolate, drawn in orange, cannot be predicted out of sample and is set aside; support is the share of isolates in the other lineages, 85.8 % here. The largest lineage holds 31.4 % of the isolates, which sets how much one lineage can weigh in the share.

**Table 2.** Conditions the estimator requires before any result is reported.

| Condition | Observed | Required | Verdict |
|---|---:|---:|---|
| Lineages with at least two isolates, fewest over the antimicrobials read | 43 | ≥ 2 | accepted |
| Support, lowest over the antimicrobials read | 85.8 % | reported | singletons set aside |

Each row uses only isolates with both a readable result for that agent and a recorded lineage. Support is the share of them in lineages of at least two isolates; the share is scored on those isolates and the singletons are set aside. The rarer-outcome count is a warning threshold, not a reporting condition of the estimator.

**Table 2a.** Per-agent input feasibility.

| Agent | Retained | Support | Input failures | Rarer outcome |
|---|---|---|---|---|
| amoxicillin | 458 | 85.8 % | none at input level | 8 (below 20) |
| cefquinome | 458 | 85.8 % | none at input level | 9 (below 20) |
| ceftiofur | 458 | 85.8 % | none at input level | 51 |
| doxycycline | 458 | 85.8 % | none at input level | 80 |
| erythromycin | 458 | 85.8 % | none at input level | 221 |
| lincomycin | 458 | 85.8 % | none at input level | 203 |
| penicillin | 458 | 85.8 % | none at input level | 46 |
| spectinomycin | 458 | 85.8 % | none at input level | 34 |
| tetracycline | 458 | 85.8 % | none at input level | 80 |
| tiamulin | 458 | 85.8 % | none at input level | 36 |
| tilmicosin | 458 | 85.8 % | none at input level | 221 |
| trimethoprim | 458 | 85.8 % | none at input level | 77 |
| tylosin | 458 | 85.8 % | none at input level | 220 |

## 3. The lineage share of the call, trait by trait

The interval is for the share of the lineages in this collection: the lineages and their sizes are held fixed and the isolates of each lineage are drawn again from the smoothed distribution of its calls. For **3 of the 13 traits** shown the interval reaches zero, so no lineage effect is distinguishable from none for that trait.

The control column is the same estimator run on shuffled lineage labels; it should sit near zero, and a share is read against it rather than against zero.

*Figure 2 (drawn on the page).* Lineage share of the call by trait, point estimate with 95 % interval for the represented lineages, 13 largest of 13. Traits whose interval reaches zero are drawn in grey and marked ‡.

The last column asks what the call says about the ordering it was cut from: the lineage share of the latent MIC, or of any continuous value the call thresholds, the rank intraclass correlation of that value. A call does not identify that share; it bounds it. The bracket gives the smallest share any ordering consistent with the calls allows, the lower bound the calls establish, and the largest; an upper end marked ≤ is a certified bound rather than a share attained. The figure after it is the one-sided 95 % lower confidence limit of the lower bound, from splitting the isolates of every lineage at random, choosing a direction on one half and testing it on the other: a share of the latent ordering that the lineages of this collection are shown to account for, with no model for the latent values beyond their agreeing with the readings. A MIC reading on a panel that contains the cut-off refines the call (Table 5), and refining the readings can only narrow the bounds they allow.

**Table 3.** The share of every trait with its interval, and the bounds the call places on the share of the latent ordering.

| Trait | Share | 95 % interval, represented lineages | Latent ordering: bounds; 95 % lower limit |
|---|---:|---:|---:|
| ceftiofur | 0.593 | 0.379 to 0.871 | 0.088 to ≤ 1.000; ≥ 0.000 |
| penicillin | 0.535 | 0.317 to 0.822 | 0.077 to ≤ 1.000; ≥ 0.000 |
| trimethoprim | 0.406 | 0.274 to 0.547 | 0.066 to ≤ 1.000; ≥ 0.000 |
| cefquinome | 0.392 | 0.000 to 1.000 | 0.020 to ≤ 1.000; ≥ 0.000 |
| tiamulin | 0.385 | 0.111 to 0.770 | 0.039 to ≤ 1.000; ≥ 0.000 |
| erythromycin | 0.287 | 0.200 to 0.371 | 0.029 to ≤ 1.000; ≥ 0.000 |
| lincomycin | 0.277 | 0.196 to 0.358 | 0.020 to ≤ 1.000; ≥ 0.000 |
| tylosin | 0.276 | 0.191 to 0.356 | 0.024 to ≤ 1.000; ≥ 0.000 |
| tilmicosin | 0.257 | 0.169 to 0.346 | 0.017 to ≤ 1.000; ≥ 0.000 |
| doxycycline | 0.254 | 0.151 to 0.364 | 0.048 to ≤ 1.000; ≥ 0.000 |
| tetracycline | 0.254 | 0.151 to 0.364 | 0.048 to ≤ 1.000; ≥ 0.000 |
| amoxicillin | 0.110 | 0.000 to 0.580 | 0.002 to ≤ 1.000; ≥ 0.000 |
| spectinomycin | 0.072 | 0.000 to 0.235 | 0.011 to ≤ 1.000; ≥ 0.000 |

## 4. Evidence that survives re-reading

A p-value is a statement about one look at the data. A surveillance panel is looked at again every year, and a p-value recomputed each time loses its error control. The e-value is evidence on a scale made for that: this run's e-value is a statement about this collection, and a programme that adds an intake each year multiplies the e-value of each new intake, scored against the lineage rates learned from the earlier ones, into a running product whose error control holds at whatever intake it is read (sequential_e_process). The e-BH procedure controls the false-discovery rate across the panel whatever the dependence between traits. Larger is stronger; 1 is no evidence.

**Table 4.** e-value per trait and the e-BH selection at level 0.05. The selection threshold on this run is 23.6; a trait at or above it is selected.

| Trait | e-value | natural log | Selected |
|---|---:|---:|---|
| ceftiofur | 5.6 × 10⁶ | 15.54 | yes |
| trimethoprim | 1.3 × 10⁶ | 14.06 | yes |
| penicillin | 5.7 × 10⁴ | 10.96 | yes |
| erythromycin | 1.2 × 10⁴ | 9.39 | yes |
| tylosin | 7.4 × 10³ | 8.91 | yes |
| tilmicosin | 2.3 × 10³ | 7.74 | yes |
| lincomycin | 1.5 × 10³ | 7.31 | yes |
| doxycycline | 1.5 × 10³ | 7.31 | yes |
| tetracycline | 1.5 × 10³ | 7.31 | yes |
| tiamulin | 266.1 | 5.58 | yes |
| cefquinome | 58.9 | 4.08 | yes |
| amoxicillin | 3.4 | 1.21 | no |
| spectinomycin | 0.8 | -0.22 | no |

11 of 13 traits are selected. Note: e-value in the betting sense of Vovk and Wang, not the BLAST expectation value and not the E-value of VanderWeele and Ding.

*Figure 3 (drawn on the page).* Evidence per trait on the natural-log scale, 13 largest of 13. The dashed rule is 1/α, the evidence one trait alone needs at level 0.05; the solid rule is the e-BH selection threshold on this run, 23.6, which rises with the number of traits read together, and a trait at or beyond it is selected. Traits the e-BH procedure selected are drawn in blue; the scale is logarithmic, so equal steps are equal factors of evidence.

## 5. Reading at the recorded resolution

Where a dilution was recorded, the lineage question is also read from the dilution itself rather than from the call derived from it. Every reading is placed by its position among the readings of its own stratum, a level of `testing_laboratory × isolation_country`: the share of those readings below it plus half the share at the same reading, with a censored reading placed by the nonparametric maximum-likelihood distribution of the readings. The share of these positions that lineage accounts for is estimated as the share of a call is, and its control shuffles the lineage labels only within strata. No distribution is assumed for the readings, a change of the MIC scale such as mg/L to log2 leaves the share unchanged, and a constant offset or a different panel between strata drops out. The share of the call and this share are read on independently retained readable and typed subsets; their denominators may differ.

The readings also bound the lineage share of the MIC ordering they were read from, the ordering of the latent MICs within a level of `testing_laboratory × isolation_country`: the rank intraclass correlation of the MIC itself. A reading says only that the MIC lies in its range, so the readings do not identify that share; they bound it. The lower bound is the smallest share that any arrangement of the latent MICs within their readings allows, the lower bound the readings establish, computed with a certified error; the upper end is the largest share, attained by an arrangement, or, marked ≤, a certified bound on it. Every share between them is allowed. The figure after the bounds is the one-sided 95 % lower confidence limit of the lower bound, and so of the share: the isolates of every lineage are split at random, one half chooses the lineage scores on which the lower bound rests and the other half tests them, and then the other way round; the tests of 25 such splits are averaged, and their average spread sets the limit. A finer panel never lowers the lower bound, and a call, one cut of the panel, never raises it. Resolution is the share of the ordering the panel resolves, one less the expected squared width of a reading on the scale of positions; it is the factor by which the share of the midpoint arrangement falls short of the plug-in share of the scores.

The share of the latent ordering is pooled over strata: a lineage placed high in one stratum and low in another does not count as ordered. For 5 agents the readings of some stratum on its own establish a larger lower bound than the pooled one, so lineage effects differ between strata: cefquinome (pooled 0.190; LGC Fordham / Canada 0.394); ceftiofur (pooled 0.155; LGC Fordham / Canada 0.258); lincomycin (pooled 0.168; LGC Fordham / United Kingdom 0.399); penicillin (pooled 0.090; LGC Fordham / Canada 0.210); trimethoprim (pooled 0.272; LGC Fordham / Canada 0.477).

**Table 5.** Lineage share of the MIC ordering per agent, with the 95 % interval for the lineages of this collection and the permutation p-value and the Benjamini-Yekutieli q-value across agents; selected is the Benjamini-Yekutieli selection at level 0.05. Control is the share on shuffled lineage labels. The bounds are those the readings place on the share of the latent MIC ordering, with the 95 % lower confidence limit of the lower bound. Readings on an end well are tied at the panel's edge and carry less of the ordering; their share is given beside the estimate. A p-value written with ≤ is the smallest the permutations can give. The last column gives the conclusion: established, lineage structure established, when the interval lies above zero and the agent is selected, or the lower confidence limit is above zero; not established; whole range, when the interval spans it; or not estimable. The conclusion says whether this collection shows lineage structure in the measurement, with the values behind it. It is not a statement about transmission or a resistance mechanism.

| Agent | Readings | Share | 95 % interval, represented lineages | p | q | Control | Selected | Latent ordering: bounds; 95 % lower limit | Resolution | End wells | Strata | Conclusion |
|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---|
| amoxicillin | 458 | 0.258 | 0.000 to 0.676 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.020 to ≤ 1.000; ≥ 0.000 | 0.093 | 96.7 % | 3 | not established |
| cefquinome | 458 | 0.377 | 0.272 to 0.510 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.190 to ≤ 0.788; ≥ 0.021 | 0.827 | 9.8 % | 3 | established |
| ceftiofur | 458 | 0.384 | 0.264 to 0.543 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.155 to ≤ 0.981; ≥ 0.015 | 0.665 | 2.0 % | 3 | established |
| doxycycline | 458 | 0.213 | 0.133 to 0.313 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.105 to ≤ 0.608; ≥ 0.000 | 0.859 | 3.7 % | 3 | established |
| enrofloxacin | 458 | 0.230 | 0.143 to 0.334 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.049 to ≤ 0.885; ≥ 0.000 | 0.759 | 0.9 % | 3 | established |
| erythromycin | 458 | 0.208 | 0.131 to 0.284 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.083 to ≤ 0.600; ≥ 0.000 | 0.870 | 57.2 % | 3 | established |
| florfenicol | 458 | 0.275 | 0.170 to 0.400 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.028 to ≤ 1.000; ≥ 0.000 | 0.525 | 4.6 % | 3 | established |
| lincomycin | 458 | 0.331 | 0.237 to 0.429 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.168 to ≤ 0.816; ≥ 0.027 | 0.816 | 52.0 % | 3 | established |
| marbofloxacin | 458 | 0.191 | 0.097 to 0.307 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.037 to ≤ 0.753; ≥ 0.000 | 0.744 | 9.2 % | 3 | established |
| penicillin | 458 | 0.437 | 0.237 to 0.656 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.090 to ≤ 1.000; ≥ 0.000 | 0.223 | 88.6 % | 3 | established |
| spectinomycin | 458 | 0.211 | 0.117 to 0.323 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.061 to ≤ 0.939; ≥ 0.000 | 0.661 | 5.2 % | 3 | established |
| tetracycline | 458 | 0.199 | 0.119 to 0.290 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.092 to ≤ 0.673; ≥ 0.000 | 0.836 | 9.0 % | 3 | established |
| tiamulin | 458 | 0.399 | 0.284 to 0.537 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.271 to ≤ 0.719; ≥ 0.077 | 0.870 | 5.0 % | 3 | established |
| tilmicosin | 458 | 0.279 | 0.191 to 0.371 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.058 to ≤ 0.981; ≥ 0.000 | 0.758 | 51.1 % | 3 | established |
| trimethoprim | 458 | 0.373 | 0.277 to 0.491 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.272 to ≤ 0.570; ≥ 0.062 | 0.882 | 9.2 % | 3 | established |
| tylosin | 458 | 0.238 | 0.160 to 0.316 | ≤ 0.001 | 0.003 | 0.00 | yes | 0.063 to ≤ 0.805; ≥ 0.000 | 0.803 | 74.9 % | 3 | established |

Resolution of the selection: with 999 permutations one agent alone cannot be selected among 16 at alpha = 0.05; at least 1081 permutations would let a single very small p-value pass the step-up.

**Table 5b.** Lineage structure behind each share of the MIC ordering: the lineages with at least two readings on which the share is scored, the readings of singleton lineages set aside, the effective number of scored lineages (inverse of the sum of squared lineage shares) and the share of readings in lineages with at least two members (support). Few effective lineages mean the share rests on a few lineages.

| Agent | Lineages scored | Singletons set aside | Effective | Support |
|---|---:|---:|---:|---:|
| amoxicillin | 43 | 65 | 6.0 | 85.8 % |
| cefquinome | 43 | 65 | 6.0 | 85.8 % |
| ceftiofur | 43 | 65 | 6.0 | 85.8 % |
| doxycycline | 43 | 65 | 6.0 | 85.8 % |
| enrofloxacin | 43 | 65 | 6.0 | 85.8 % |
| erythromycin | 43 | 65 | 6.0 | 85.8 % |
| florfenicol | 43 | 65 | 6.0 | 85.8 % |
| lincomycin | 43 | 65 | 6.0 | 85.8 % |
| marbofloxacin | 43 | 65 | 6.0 | 85.8 % |
| penicillin | 43 | 65 | 6.0 | 85.8 % |
| spectinomycin | 43 | 65 | 6.0 | 85.8 % |
| tetracycline | 43 | 65 | 6.0 | 85.8 % |
| tiamulin | 43 | 65 | 6.0 | 85.8 % |
| tilmicosin | 43 | 65 | 6.0 | 85.8 % |
| trimethoprim | 43 | 65 | 6.0 | 85.8 % |
| tylosin | 43 | 65 | 6.0 | 85.8 % |

*Figure 4 (drawn on the page).* Lineage share of the MIC ordering per agent, 16 largest of 16: the point estimate with its 95 % interval; agents selected across the panel are drawn in blue. The hollow diamond is the share of the binary call from Table 1 for the same agent, a different quantity whose retained isolates may differ, placed here so that the two readings can be seen side by side. The figure in parentheses at the right is the share of readings on an end well.

Each laboratory is read on the wells that laboratory tested: an end-well reading is censored at the panel edge of its own laboratory, not at the widest edge in the collection. The table below gives, per agent and laboratory, the panel the readings were taken to come from.

**Table 5a.** Panel geometry per agent and testing laboratory: the wells taken as tested, their range, whether they form a doubling series, and the shares of readings on the lowest and the highest well.

| Agent | Laboratory | Wells | Range | Doubling | On lowest well | On highest well |
|---|---|---:|---:|---|---:|---:|
| amoxicillin | LGC Fordham | 8 | 0.0312 to 4.0000 | yes | 96.1 % | 0.2 % |
| amoxicillin | OUCRU Ho Chi Minh City | 2 | 0.0156 to 0.0312 | yes | 98.0 % | 2.0 % |
| cefquinome | LGC Fordham | 11 | 0.0020 to 2.0000 | yes | 1.2 % | 0.7 % |
| cefquinome | OUCRU Ho Chi Minh City | 3 | 0.0078 to 0.0312 | yes | 4.1 % | 71.4 % |
| ceftiofur | LGC Fordham | 10 | 0.0312 to 16.0000 | yes | 0.7 % | 0.7 % |
| ceftiofur | OUCRU Ho Chi Minh City | 4 | 0.0625 to 0.5000 | yes | 4.1 % | 2.0 % |
| doxycycline | LGC Fordham | 12 | 0.0312 to 64.0000 | yes | 0.2 % | 0.2 % |
| doxycycline | OUCRU Ho Chi Minh City | 9 | 0.0312 to 8.0000 | yes | 8.2 % | 22.4 % |
| enrofloxacin | LGC Fordham | 11 | 0.0078 to 8.0000 | yes | 0.2 % | 0.2 % |
| enrofloxacin | OUCRU Ho Chi Minh City | 5 | 0.2500 to 4.0000 | yes | 2.0 % | 2.0 % |
| erythromycin | LGC Fordham | 12 | 0.0156 to 32.0000 | yes | 10.3 % | 45.2 % |
| erythromycin | OUCRU Ho Chi Minh City | 12 | 0.0156 to 32.0000 | yes | 57.1 % | 14.3 % |
| florfenicol | LGC Fordham | 4 | 0.5000 to 4.0000 | yes | 0.7 % | 3.9 % |
| florfenicol | OUCRU Ho Chi Minh City | 5 | 0.2500 to 4.0000 | yes | 2.0 % | 2.0 % |
| lincomycin | LGC Fordham | 12 | 0.0625 to 128.0000 | yes | 1.2 % | 54.0 % |
| lincomycin | OUCRU Ho Chi Minh City | 12 | 0.0625 to 128.0000 | yes | 2.0 % | 22.4 % |
| marbofloxacin | LGC Fordham | 8 | 0.0156 to 2.0000 | yes | 6.6 % | 1.7 % |
| marbofloxacin | OUCRU Ho Chi Minh City | 3 | 0.2500 to 1.0000 | yes | 6.1 % | 10.2 % |
| penicillin | LGC Fordham | 7 | 0.0312 to 2.0000 | yes | 85.8 % | 1.5 % |
| penicillin | OUCRU Ho Chi Minh City | 5 | 0.0312 to 0.5000 | yes | 95.9 % | 4.1 % |
| spectinomycin | LGC Fordham | 10 | 1.0000 to 512.0000 | yes | 0.2 % | 4.2 % |
| spectinomycin | OUCRU Ho Chi Minh City | 8 | 2.0000 to 256.0000 | yes | 2.0 % | 10.2 % |
| tetracycline | LGC Fordham | 12 | 0.0625 to 128.0000 | yes | 0.2 % | 5.6 % |
| tetracycline | OUCRU Ho Chi Minh City | 9 | 0.5000 to 128.0000 | yes | 10.2 % | 24.5 % |
| tiamulin | LGC Fordham | 10 | 0.1250 to 64.0000 | yes | 2.0 % | 2.9 % |
| tiamulin | OUCRU Ho Chi Minh City | 7 | 0.0312 to 2.0000 | yes | 4.1 % | 2.0 % |
| tilmicosin | LGC Fordham | 10 | 0.2500 to 128.0000 | yes | 0.2 % | 53.5 % |
| tilmicosin | OUCRU Ho Chi Minh City | 12 | 0.5000 to 1024.0000 | yes | 4.1 % | 24.5 % |
| trimethoprim | LGC Fordham | 12 | 0.0156 to 32.0000 | yes | 4.9 % | 2.7 % |
| trimethoprim | OUCRU Ho Chi Minh City | 8 | 0.0625 to 8.0000 | yes | 20.4 % | 2.0 % |
| tylosin | LGC Fordham | 10 | 0.5000 to 256.0000 | yes | 27.9 % | 52.6 % |
| tylosin | OUCRU Ho Chi Minh City | 12 | 0.1250 to 256.0000 | yes | 10.2 % | 18.4 % |

*Figure 5 (drawn on the page).* Readings per lineage and dilution interval for 6 of 16 agents, the 14 largest of 108 lineages: each cell counts the isolates of one lineage whose reading fell in one interval of the panel, darker for more. Open intervals at the panel edges are the censored readings. A lineage whose readings sit in one or two adjacent cells contributes to the lineage share; readings spread along a row do not.

## 6. A change of lineages or a change within them

Contrast isolation_country: Canada minus United Kingdom: of 13 traits, 0 have a composition component (a change in lineage mix) and 0 a within-lineage component (a change in rate) selected by the Benjamini-Yekutieli step-up within each component family; the intervals and the step-up rest on nominal percentile-bootstrap tail probabilities (their calibration is reported with the package).

Both components were refused for every trait: the lineage labels are missing unevenly between the two collections and in association with the trait, so the labelled isolates are not the collections and neither component is counted above.

**Table 6.** The prevalence difference of each trait, Canada minus United Kingdom, split into a change in lineage composition and a change in rate within lineages. The two components sum to the difference. The limits are bootstrap percentiles; p is the bootstrap p-value of each component and q its Benjamini-Yekutieli q-value within the component family; the last column names the components the false-discovery procedure selected.

| Trait | Difference | Composition | 95 % interval | p | q | Within lineage | 95 % interval | p | q | Selected |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| amoxicillin | 0.039 | 0.023 | -0.012 to 0.063 | 0.206 | 0.851 | 0.016 | 0.000 to 0.053 | 0.752 | 1.000 | neither |
| cefquinome | 0.036 | 0.036 | 0.001 to 0.075 | 0.044 | 0.260 | 0.000 | 0.000 to 0.000 | 1.000 | 1.000 | neither |
| ceftiofur | 0.037 | 0.032 | -0.034 to 0.104 | 0.373 | 1.000 | 0.005 | 0.000 to 0.016 | 0.956 | 1.000 | neither |
| doxycycline | 0.174 | 0.049 | -0.004 to 0.131 | 0.068 | 0.312 | 0.125 | 0.044 to 0.171 | 0.007 | 0.048 | neither |
| erythromycin | 0.283 | 0.145 | 0.070 to 0.247 | 0.001 | 0.014 | 0.138 | 0.037 to 0.214 | 0.005 | 0.048 | neither |
| lincomycin | 0.250 | 0.098 | 0.024 to 0.199 | 0.016 | 0.110 | 0.151 | 0.048 to 0.224 | 0.003 | 0.041 | neither |
| penicillin | 0.024 | 0.022 | -0.039 to 0.090 | 0.510 | 1.000 | 0.002 | -0.008 to 0.014 | 1.000 | 1.000 | neither |
| spectinomycin | 0.079 | 0.058 | 0.016 to 0.106 | 0.007 | 0.058 | 0.021 | -0.002 to 0.048 | 0.075 | 0.443 | neither |
| tetracycline | 0.174 | 0.049 | -0.004 to 0.131 | 0.068 | 0.312 | 0.125 | 0.044 to 0.171 | 0.007 | 0.048 | neither |
| tiamulin | -0.015 | -0.027 | -0.078 to 0.024 | 0.296 | 1.000 | 0.012 | -0.006 to 0.036 | 0.239 | 1.000 | neither |
| tilmicosin | 0.291 | 0.145 | 0.073 to 0.247 | 0.001 | 0.014 | 0.146 | 0.042 to 0.217 | 0.003 | 0.041 | neither |
| trimethoprim | -0.117 | -0.146 | -0.206 to -0.072 | 0.001 | 0.014 | 0.028 | -0.031 to 0.084 | 0.424 | 1.000 | neither |
| tylosin | 0.286 | 0.142 | 0.068 to 0.243 | 0.002 | 0.021 | 0.145 | 0.040 to 0.219 | 0.003 | 0.041 | neither |

**Table 6b.** The two collections behind each difference: the prevalence and the isolates of each, the lineages both hold and those seen in one only, the share of isolates in lineages both hold (shared support, the mean of the two collections; the within-lineage component needs at least 0.90) and the turnover share, the part of the difference carried by lineages seen in one collection only.

| Trait | Prevalence, Canada | n | Prevalence, United Kingdom | n | Lineages shared; only first; only second | Shared support | Turnover share |
|---|---:|---:|---:|---:|---:|---:|---:|
| amoxicillin | 4.7 % | 129 | 0.7 % | 280 | 8; 26; 74 | 53.4 % | 0.856 |
| cefquinome | 4.7 % | 129 | 1.1 % | 280 | 8; 26; 74 | 53.4 % | 1.000 |
| ceftiofur | 14.7 % | 129 | 11.1 % | 280 | 8; 26; 74 | 53.4 % | 0.970 |
| doxycycline | 93.8 % | 129 | 76.4 % | 280 | 8; 26; 74 | 53.4 % | 0.657 |
| erythromycin | 73.6 % | 129 | 45.4 % | 280 | 8; 26; 74 | 53.4 % | 0.621 |
| lincomycin | 76.7 % | 129 | 51.8 % | 280 | 8; 26; 74 | 53.4 % | 0.632 |
| penicillin | 12.4 % | 129 | 10.0 % | 280 | 8; 26; 74 | 53.4 % | 0.949 |
| spectinomycin | 10.1 % | 129 | 2.1 % | 280 | 8; 26; 74 | 53.4 % | 0.624 |
| tetracycline | 93.8 % | 129 | 76.4 % | 280 | 8; 26; 74 | 53.4 % | 0.657 |
| tiamulin | 7.8 % | 129 | 9.3 % | 280 | 8; 26; 74 | 53.4 % | 0.822 |
| tilmicosin | 74.4 % | 129 | 45.4 % | 280 | 8; 26; 74 | 53.4 % | 0.617 |
| trimethoprim | 10.1 % | 129 | 21.8 % | 280 | 8; 26; 74 | 53.4 % | 0.731 |
| tylosin | 73.6 % | 129 | 45.0 % | 280 | 8; 26; 74 | 53.4 % | 0.620 |

## 7. Provenance and terms

- **Software:** amr-clonalshare 1.0.0, record schema 1.0
- **Seed:** 42
- **Configuration:** `sha256:a9a983a1ba3b`
- **Record:** `sha256:ae9f745ec7a07ea3a49aef8a9bf81adbc5133df5e69d4fb11cdec5a827c7d34b`

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

amr-clonalshare 1.0.0 · record schema 1.0 · seed 42 · configuration sha256:a9a983a1ba3b
Cite this run as: “amr-clonalshare 1.0.0, run F3CC72CB, record sha256:ae9f745ec7a0.”
This report supersedes any earlier report bearing the same run identifier. It is regenerated from the record and holds no value that the record does not. The symbols carry the same wording in every run of this software; no symbol against a value means only that none of the listed conditions fired.
