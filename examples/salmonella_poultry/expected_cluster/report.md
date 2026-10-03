# amr-clonalshare run report

Every number below is read from `clonal_share_result.json`, written by the same run; none is recomputed here.

- **Run:** CE314965
- **Issued:** 2026-10-03T06:12Z
- **Software:** amr-clonalshare 1.0.0
- **Isolates:** 7049
- **Lineage:** pds_cluster
- **Seed:** 42
- **Configuration:** sha256:4890efe08a73
- **Record digest:** sha256:d36d68a5934104b98ea8d69f4974452ec886cd8ac788a9980b2759e9fd7cb241

**22 antimicrobials; 110 analysis outcomes.** 95 computed, 2 full ranges, 13 unavailable. Each result retains its own target and data subset.

## 1. Measurement summary

**Quantity measured.** The share of the variation in the recorded binary outcome across this collection that lineage membership accounts for, read from a lineage label and an interpreted result and scored on isolates the estimator did not see.

- **Collection:** 7049 isolates in 1 lineages, typed by `pds_cluster`
- **Antimicrobials read:** 22
- **Membership shares estimable:** 18 of 22

**Reporting conditions.** 95 computed, 2 full ranges, 13 unavailable. Each result retains its own target and data subset. An unavailable value failed a declared reporting condition of its estimator; the recorded reason distinguishes input limitations from numerical failure.

**Table 0.** Available analyses and their data subsets. A full range is a completed, uninformative result; an unavailable value is not zero.

| Agent | Analysis | N | Status | Data scope |
|---|---|---:|---|---|
| amikacin | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| amoxicillin-clavulanic acid | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| ampicillin | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| azithromycin | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| cefoxitin | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| ceftazidime | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| ceftiofur | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| ceftriaxone | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| chloramphenicol | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| ciprofloxacin | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| colistin | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| gentamicin | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| imipenem | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| kanamycin | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| meropenem | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| nalidixic acid | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| spectinomycin | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| streptomycin | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| sulfamethoxazole | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| sulfisoxazole | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| tetracycline | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| trimethoprim-sulfamethoxazole | finite_collection_prevalence | 7049 | computed | recorded finite collection |
| amikacin | collection_membership | 2401 | full_range | observed and typed subset |
| amikacin | call_latent_bounds | 2401 | full_range | observed and typed subset |
| amikacin | call_latent_lower_limit | 2401 | computed | observed and typed subset |
| amoxicillin-clavulanic acid | collection_membership | 6917 | computed | observed and typed subset |
| amoxicillin-clavulanic acid | call_latent_bounds | 6917 | computed | observed and typed subset |
| amoxicillin-clavulanic acid | call_latent_lower_limit | 6917 | computed | observed and typed subset |
| ampicillin | collection_membership | 7049 | computed | observed and typed subset |
| ampicillin | call_latent_bounds | 7049 | computed | observed and typed subset |
| ampicillin | call_latent_lower_limit | 7049 | computed | observed and typed subset |
| azithromycin | collection_membership | 4779 | computed | observed and typed subset |
| azithromycin | call_latent_bounds | 4779 | computed | observed and typed subset |
| azithromycin | call_latent_lower_limit | 4779 | computed | observed and typed subset |
| cefoxitin | collection_membership | 6917 | computed | observed and typed subset |
| cefoxitin | call_latent_bounds | 6917 | computed | observed and typed subset |
| cefoxitin | call_latent_lower_limit | 6917 | computed | observed and typed subset |
| ceftazidime | collection_membership | 132 | computed | observed and typed subset |
| ceftazidime | call_latent_bounds | 132 | computed | observed and typed subset |
| ceftazidime | call_latent_lower_limit | 132 | computed | observed and typed subset |
| ceftiofur | collection_membership | 3616 | computed | observed and typed subset |
| ceftiofur | call_latent_bounds | 3616 | computed | observed and typed subset |
| ceftiofur | call_latent_lower_limit | 3616 | computed | observed and typed subset |
| ceftriaxone | collection_membership | 6917 | computed | observed and typed subset |
| ceftriaxone | call_latent_bounds | 6917 | computed | observed and typed subset |
| ceftriaxone | call_latent_lower_limit | 6917 | computed | observed and typed subset |
| chloramphenicol | collection_membership | 7049 | computed | observed and typed subset |
| chloramphenicol | call_latent_bounds | 7049 | computed | observed and typed subset |
| chloramphenicol | call_latent_lower_limit | 7049 | computed | observed and typed subset |
| ciprofloxacin | collection_membership | 6935 | computed | observed and typed subset |
| ciprofloxacin | call_latent_bounds | 6935 | computed | observed and typed subset |
| ciprofloxacin | call_latent_lower_limit | 6935 | computed | observed and typed subset |
| colistin | collection_membership | 525 | unavailable | observed and typed subset |
| colistin | call_latent_bounds | 525 | unavailable | observed and typed subset |
| colistin | call_latent_lower_limit | 525 | unavailable | observed and typed subset |
| gentamicin | collection_membership | 7047 | computed | observed and typed subset |
| gentamicin | call_latent_bounds | 7047 | computed | observed and typed subset |
| gentamicin | call_latent_lower_limit | 7047 | computed | observed and typed subset |
| imipenem | collection_membership | 132 | unavailable | observed and typed subset |
| imipenem | call_latent_bounds | 132 | unavailable | observed and typed subset |
| imipenem | call_latent_lower_limit | 132 | unavailable | observed and typed subset |
| kanamycin | collection_membership | 3206 | computed | observed and typed subset |
| kanamycin | call_latent_bounds | 3206 | computed | observed and typed subset |
| kanamycin | call_latent_lower_limit | 3206 | computed | observed and typed subset |
| meropenem | collection_membership | 3300 | unavailable | observed and typed subset |
| meropenem | call_latent_bounds | 3300 | unavailable | observed and typed subset |
| meropenem | call_latent_lower_limit | 3300 | unavailable | observed and typed subset |
| nalidixic acid | collection_membership | 6915 | computed | observed and typed subset |
| nalidixic acid | call_latent_bounds | 6915 | computed | observed and typed subset |
| nalidixic acid | call_latent_lower_limit | 6915 | computed | observed and typed subset |
| spectinomycin | collection_membership | 1 | unavailable | observed and typed subset |
| spectinomycin | call_latent_bounds | 1 | unavailable | observed and typed subset |
| spectinomycin | call_latent_lower_limit | 1 | unavailable | observed and typed subset |
| streptomycin | collection_membership | 6520 | computed | observed and typed subset |
| streptomycin | call_latent_bounds | 6520 | computed | observed and typed subset |
| streptomycin | call_latent_lower_limit | 6520 | computed | observed and typed subset |
| sulfamethoxazole | collection_membership | 253 | computed | observed and typed subset |
| sulfamethoxazole | call_latent_bounds | 253 | computed | observed and typed subset |
| sulfamethoxazole | call_latent_lower_limit | 253 | computed | observed and typed subset |
| sulfisoxazole | collection_membership | 6662 | computed | observed and typed subset |
| sulfisoxazole | call_latent_bounds | 6662 | computed | observed and typed subset |
| sulfisoxazole | call_latent_lower_limit | 6662 | computed | observed and typed subset |
| tetracycline | collection_membership | 7047 | computed | observed and typed subset |
| tetracycline | call_latent_bounds | 7047 | computed | observed and typed subset |
| tetracycline | call_latent_lower_limit | 7047 | computed | observed and typed subset |
| trimethoprim-sulfamethoxazole | collection_membership | 6915 | computed | observed and typed subset |
| trimethoprim-sulfamethoxazole | call_latent_bounds | 6915 | computed | observed and typed subset |
| trimethoprim-sulfamethoxazole | call_latent_lower_limit | 6915 | computed | observed and typed subset |
| amikacin | lineage_evidence | 2401 | computed | observed and typed subset |
| amoxicillin-clavulanic acid | lineage_evidence | 6917 | computed | observed and typed subset |
| ampicillin | lineage_evidence | 7049 | computed | observed and typed subset |
| azithromycin | lineage_evidence | 4779 | computed | observed and typed subset |
| cefoxitin | lineage_evidence | 6917 | computed | observed and typed subset |
| ceftazidime | lineage_evidence | 132 | computed | observed and typed subset |
| ceftiofur | lineage_evidence | 3616 | computed | observed and typed subset |
| ceftriaxone | lineage_evidence | 6917 | computed | observed and typed subset |
| chloramphenicol | lineage_evidence | 7049 | computed | observed and typed subset |
| ciprofloxacin | lineage_evidence | 6935 | computed | observed and typed subset |
| colistin | lineage_evidence | 525 | computed | observed and typed subset |
| gentamicin | lineage_evidence | 7047 | computed | observed and typed subset |
| imipenem | lineage_evidence | 132 | computed | observed and typed subset |
| kanamycin | lineage_evidence | 3206 | computed | observed and typed subset |
| meropenem | lineage_evidence | 3300 | computed | observed and typed subset |
| nalidixic acid | lineage_evidence | 6915 | computed | observed and typed subset |
| spectinomycin | lineage_evidence | 1 | unavailable | observed and typed subset |
| streptomycin | lineage_evidence | 6520 | computed | observed and typed subset |
| sulfamethoxazole | lineage_evidence | 253 | computed | observed and typed subset |
| sulfisoxazole | lineage_evidence | 6662 | computed | observed and typed subset |
| tetracycline | lineage_evidence | 7047 | computed | observed and typed subset |
| trimethoprim-sulfamethoxazole | lineage_evidence | 6915 | computed | observed and typed subset |

**Table 0b.** Reasons and next checks for results requiring interpretation. These checks do not change the recorded status or justify changing reporting conditions.

| Agent | Analysis / status | Reason | Next check |
|---|---|---|---|
| amikacin | collection_membership / full_range | The completed interval or bounds span the entire admissible range [0, 1] | Report the full range as uninformative; inspect repeated-lineage support and outcome variation. |
| amikacin | call_latent_bounds / full_range | The largest share attained by an arrangement of the lineages inside the readings is 0.997 | Report the full range as uninformative; inspect repeated-lineage support and outcome variation. |
| amoxicillin-clavulanic acid | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.955 | Read the interval kind, assumptions and retained collection before comparing results. |
| ampicillin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.903 | Read the interval kind, assumptions and retained collection before comparing results. |
| azithromycin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.997 | Read the interval kind, assumptions and retained collection before comparing results. |
| cefoxitin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.984 | Read the interval kind, assumptions and retained collection before comparing results. |
| ceftazidime | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.915 | Read the interval kind, assumptions and retained collection before comparing results. |
| ceftiofur | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.969 | Read the interval kind, assumptions and retained collection before comparing results. |
| ceftriaxone | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.975 | Read the interval kind, assumptions and retained collection before comparing results. |
| chloramphenicol | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.996 | Read the interval kind, assumptions and retained collection before comparing results. |
| ciprofloxacin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.998 | Read the interval kind, assumptions and retained collection before comparing results. |
| colistin | collection_membership / unavailable | Constant retained outcome; no variation to attribute | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| colistin | call_latent_bounds / unavailable | The largest share attained by an arrangement of the lineages inside the readings is 0.986 | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| colistin | call_latent_lower_limit / unavailable | Constant retained outcome; no variation to attribute | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| gentamicin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.961 | Read the interval kind, assumptions and retained collection before comparing results. |
| imipenem | collection_membership / unavailable | Constant retained outcome; no variation to attribute | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| imipenem | call_latent_bounds / unavailable | The largest share attained by an arrangement of the lineages inside the readings is 0.964 | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| imipenem | call_latent_lower_limit / unavailable | Constant retained outcome; no variation to attribute | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| kanamycin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.971 | Read the interval kind, assumptions and retained collection before comparing results. |
| meropenem | collection_membership / unavailable | Constant retained outcome; no variation to attribute | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| meropenem | call_latent_bounds / unavailable | The largest share attained by an arrangement of the lineages inside the readings is 0.994 | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| meropenem | call_latent_lower_limit / unavailable | Constant retained outcome; no variation to attribute | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| nalidixic acid | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.998 | Read the interval kind, assumptions and retained collection before comparing results. |
| spectinomycin | collection_membership / unavailable | Only one retained lineage; no between-lineage comparison | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| spectinomycin | call_latent_bounds / unavailable | Only one retained lineage; no between-lineage comparison | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| spectinomycin | call_latent_lower_limit / unavailable | Only one retained lineage; no between-lineage comparison | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| streptomycin | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.944 | Read the interval kind, assumptions and retained collection before comparing results. |
| sulfamethoxazole | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.971 | Read the interval kind, assumptions and retained collection before comparing results. |
| sulfisoxazole | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.971 | Read the interval kind, assumptions and retained collection before comparing results. |
| tetracycline | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.972 | Read the interval kind, assumptions and retained collection before comparing results. |
| trimethoprim-sulfamethoxazole | call_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.998 | Read the interval kind, assumptions and retained collection before comparing results. |
| spectinomycin | lineage_evidence / unavailable | Only one retained lineage; no between-lineage comparison | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |

**Table 0a.** Prevalence in the recorded collection: bounds allow every missing outcome to be either negative or positive. These are identification bounds, not confidence intervals. Observed prevalence uses observed outcomes only.

| Agent | Recorded | Observed | Missing | Observed prevalence | Collection bounds |
|---|---:|---:|---:|---:|---|
| amikacin | 7049 | 2401 | 4648 | 0.0 % | 0.000 to 0.660 |
| amoxicillin-clavulanic acid | 7049 | 6917 | 132 | 16.8 % | 0.165 to 0.184 |
| ampicillin | 7049 | 7049 | 0 | 26.4 % | 0.264 to 0.264 |
| azithromycin | 7049 | 4779 | 2270 | 0.5 % | 0.003 to 0.325 |
| cefoxitin | 7049 | 6917 | 132 | 9.4 % | 0.092 to 0.111 |
| ceftazidime | 7049 | 132 | 6917 | 28.0 % | 0.005 to 0.987 |
| ceftiofur | 7049 | 3616 | 3433 | 14.6 % | 0.075 to 0.562 |
| ceftriaxone | 7049 | 6917 | 132 | 13.0 % | 0.127 to 0.146 |
| chloramphenicol | 7049 | 7049 | 0 | 6.0 % | 0.060 to 0.060 |
| ciprofloxacin | 7049 | 6935 | 114 | 8.7 % | 0.086 to 0.102 |
| colistin | 7049 | 525 | 6524 | 100.0 % | 0.074 to 1.000 |
| gentamicin | 7049 | 7047 | 2 | 14.3 % | 0.143 to 0.143 |
| imipenem | 7049 | 132 | 6917 | 0.0 % | 0.000 to 0.981 |
| kanamycin | 7049 | 3206 | 3843 | 13.3 % | 0.061 to 0.606 |
| meropenem | 7049 | 3300 | 3749 | 0.0 % | 0.000 to 0.532 |
| nalidixic acid | 7049 | 6915 | 134 | 8.4 % | 0.083 to 0.102 |
| spectinomycin | 7049 | 1 | 7048 | 0.0 % | 0.000 to 1.000 |
| streptomycin | 7049 | 6520 | 529 | 40.8 % | 0.377 to 0.452 |
| sulfamethoxazole | 7049 | 253 | 6796 | 21.3 % | 0.008 to 0.972 |
| sulfisoxazole | 7049 | 6662 | 387 | 28.9 % | 0.273 to 0.328 |
| tetracycline | 7049 | 7047 | 2 | 52.6 % | 0.526 to 0.526 |
| trimethoprim-sulfamethoxazole | 7049 | 6915 | 134 | 3.0 % | 0.030 to 0.049 |

**Table 1.** Every trait, ordered by share, with the 95 % interval, the permutation p-value and the Benjamini-Yekutieli q-value across traits; selected is the Benjamini-Yekutieli selection at a false-discovery level of 0.05, which holds whatever the dependence between traits. Control is the share on shuffled lineage labels. A p-value written with ≤ is the smallest the permutations can give. The last column gives the conclusion of the results table of the form: established, lineage structure established, when the interval lies above zero and the trait is selected, or the lower confidence limit is above zero; not established; whole range, when the interval spans it; or not estimable. The conclusion says whether this collection shows lineage structure in the measurement, with the values behind it. It is not a statement about transmission or a resistance mechanism.

| Trait | Share | 95 % interval | p | q | Control | e-value | Selected | Conclusion |
|---|---:|---:|---:|---:|---:|---:|---|---|
| nalidixic acid | 0.943 | 0.910 to 0.977 | ≤ 0.007 | 0.025 | 0.00 | 2.5 × 10¹⁷⁹ | yes | established |
| ciprofloxacin | 0.912 | 0.874 to 0.953 | ≤ 0.007 | 0.025 | 0.00 | 1.6 × 10¹⁷⁴ | yes | established |
| tetracycline | 0.686 | 0.667 to 0.706 | ≤ 0.007 | 0.025 | 0.00 | 1.1 × 10²⁵⁷ | yes | established |
| sulfisoxazole | 0.646 | 0.627 to 0.668 | ≤ 0.007 | 0.025 | 0.00 | 1.2 × 10²⁰⁴ | yes | established |
| streptomycin | 0.574 | 0.554 to 0.599 | ≤ 0.007 | 0.025 | 0.00 | 7.2 × 10¹⁹¹ | yes | established |
| sulfamethoxazole | 0.529 | 0.388 to 0.694 | ≤ 0.007 | 0.025 | 0.00 | 2.4 × 10⁵ | yes | established |
| chloramphenicol | 0.463 | 0.407 to 0.515 | ≤ 0.007 | 0.025 | 0.00 | 1.3 × 10⁸³ | yes | established |
| kanamycin | 0.408 | 0.365 to 0.455 | ≤ 0.007 | 0.025 | 0.00 | 1.9 × 10⁵¹ | yes | established |
| ceftiofur | 0.403 | 0.371 to 0.442 | ≤ 0.007 | 0.025 | 0.00 | 1.4 × 10⁵⁴ | yes | established |
| ceftazidime | 0.399 | 0.213 to 0.651 | ≤ 0.007 | 0.025 | 0.00 | 95.8 | yes | established |
| cefoxitin | 0.398 | 0.366 to 0.430 | ≤ 0.007 | 0.025 | 0.00 | 2.6 × 10⁸¹ | yes | established |
| ceftriaxone | 0.392 | 0.365 to 0.417 | ≤ 0.007 | 0.025 | 0.00 | 7.5 × 10¹⁰⁰ | yes | established |
| gentamicin | 0.371 | 0.341 to 0.399 | ≤ 0.007 | 0.025 | 0.00 | 9.9 × 10⁹⁸ | yes | established |
| amoxicillin-clavulanic acid | 0.357 | 0.334 to 0.378 | ≤ 0.007 | 0.025 | 0.00 | 2.7 × 10⁹⁷ | yes | established |
| ampicillin | 0.337 | 0.318 to 0.356 | ≤ 0.007 | 0.025 | 0.00 | 1.5 × 10¹¹¹ | yes | established |
| trimethoprim-sulfamethoxazole | 0.324 | 0.261 to 0.382 | ≤ 0.007 | 0.025 | 0.00 | 3.5 × 10⁵¹ | yes | established |
| azithromycin | 0.128 | 0.000 to 0.408 | ≤ 0.007 | 0.025 | 0.00 | 1.0 × 10⁵ | yes | not established |
| amikacin | 0.000 | 0.000 to 1.000 | 0.629 | 1.000 | 0.00 | 0.6 | no | whole range |

**Table 1b.** Lineage structure behind each share: the lineages with at least two isolates on which the share is scored, the isolates of singleton lineages set aside, the effective number of scored lineages (inverse of the sum of squared lineage shares) and the share of isolates in lineages with at least two members (support). Few effective lineages mean the share rests on a few lineages.

| Trait | Lineages scored | Singletons set aside | Effective | Support |
|---|---:|---:|---:|---:|
| nalidixic acid | 307 | 189 | 34.0 | 97.3 % |
| ciprofloxacin | 325 | 206 | 35.5 | 97.0 % |
| tetracycline | 327 | 207 | 34.2 | 97.1 % |
| sulfisoxazole | 297 | 180 | 33.1 | 97.3 % |
| streptomycin | 314 | 199 | 33.6 | 96.9 % |
| sulfamethoxazole | 33 | 28 | 15.3 | 88.9 % |
| chloramphenicol | 327 | 207 | 34.2 | 97.1 % |
| kanamycin | 186 | 124 | 23.2 | 96.1 % |
| ceftiofur | 200 | 136 | 24.9 | 96.2 % |
| ceftazidime | 23 | 24 | 7.4 | 81.8 % |
| cefoxitin | 307 | 189 | 34.0 | 97.3 % |
| ceftriaxone | 307 | 189 | 34.0 | 97.3 % |
| gentamicin | 327 | 207 | 34.2 | 97.1 % |
| amoxicillin-clavulanic acid | 307 | 189 | 34.0 | 97.3 % |
| ampicillin | 327 | 207 | 34.2 | 97.1 % |
| trimethoprim-sulfamethoxazole | 307 | 189 | 34.0 | 97.3 % |
| azithromycin | 242 | 166 | 26.7 | 96.5 % |
| amikacin | 163 | 111 | 23.7 | 95.4 % |

> **Not evaluated on this run, and why**
>
> `colistin`: not scored; every isolate carries the same call, so there is no variance to attribute.
>
> `imipenem`: not scored; every isolate carries the same call, so there is no variance to attribute.
>
> `meropenem`: not scored; every isolate carries the same call, so there is no variance to attribute.
>
> `spectinomycin`: not scored; the isolates with a result sit in a single lineage, so there is no contrast between lineages.
>
> Decomposition of a prevalence difference into lineage composition and within-lineage rate: needs two collections; this record holds one.
>
> Reading at the recorded dilution: no dilutions were supplied; shares are read from binary calls.
>

> **Interpretation (generated from Table 1 by fixed rules)**
>
> For nalidixic acid, ciprofloxacin, tetracycline, sulfisoxazole, streptomycin and sulfamethoxazole, lineage membership is associated with the recorded outcome: the estimated share is at or above one half in this collection at the chosen typing resolution.
>
> For chloramphenicol, kanamycin, ceftiofur, ceftazidime, cefoxitin, ceftriaxone, gentamicin, amoxicillin-clavulanic acid, ampicillin and trimethoprim-sulfamethoxazole, a lineage association is detected and the estimated share is below one half. Most variation is not explained by the lineage labels at this typing resolution.
>
> For azithromycin, the panel selection finds a lineage effect while the interval for its size still reaches zero: there is evidence that positive outcomes are not spread evenly across the lineages, and this collection is too small, or too uneven across its lineages, to say how much of it the lineages account for. Read the selection as the finding and the share as not yet resolved.
>
> For amikacin, no lineage effect is distinguishable from none in this collection at the chosen resolution. This does not establish independence from lineage or identify the reason for a prevalence change.
>
> These readings describe association in the sampled collection. The analysis does not identify transmission, horizontal transfer, selection, or the effect of an intervention.
>

## 2. Admissibility of the input

- **Lineage column:** `pds_cluster`
- **Phenotype interpretation:** binary
- **Positive outcome:** non-susceptible (I or R) as recorded
- **Applied coding:** 1 (positive)
- **Intermediate policy:** not applicable
- **Interpretation source:** NCBI Pathogen Detection AST table, release PDG000000002.4210
- **AST standard/version:** not supplied / not supplied
- **Susceptibility calls:** 22 antimicrobials on 7049 isolates
- **Resampling:** 5 folds, 10 repeats, 300 bootstrap draws, 150 permutations per antimicrobial

Antimicrobials with fewer than 20 isolates of the rarer outcome: 1 of 22. Each method uses its own reporting conditions; sparse outcomes may yield wide intervals or an unavailable result.

Lineage groups in the input: 534, of which 207 hold a single isolate. Support, the share of isolates in lineages of at least two, is 97.1 %; the share is scored on those isolates and the singletons are set aside. At least two lineages hold two or more isolates.

*Figure 1 (drawn on the page).* Isolates per lineage, largest first, the 40 largest of 534. A lineage of one isolate, drawn in orange, cannot be predicted out of sample and is set aside; support is the share of isolates in the other lineages, 97.1 % here. The largest lineage holds 7.5 % of the isolates, which sets how much one lineage can weigh in the share.

**Table 2.** Conditions the estimator requires before any result is reported.

| Condition | Observed | Required | Verdict |
|---|---:|---:|---|
| Lineages with at least two isolates, fewest over the antimicrobials read | 0 | ≥ 2 | refused |
| Support, lowest over the antimicrobials read | 0.0 % | reported | singletons set aside |

Each row uses only isolates with both a readable result for that agent and a recorded lineage. Support is the share of them in lineages of at least two isolates; the share is scored on those isolates and the singletons are set aside. The rarer-outcome count is a warning threshold, not a reporting condition of the estimator.

**Table 2a.** Per-agent input feasibility.

| Agent | Retained | Support | Input failures | Rarer outcome |
|---|---|---|---|---|
| amikacin | 2401 | 95.4 % | none at input level | 1 (below 20) |
| amoxicillin-clavulanic acid | 6917 | 97.3 % | none at input level | 1163 |
| ampicillin | 7049 | 97.1 % | none at input level | 1864 |
| azithromycin | 4779 | 96.5 % | none at input level | 24 |
| cefoxitin | 6917 | 97.3 % | none at input level | 649 |
| ceftazidime | 132 | 81.8 % | none at input level | 37 |
| ceftiofur | 3616 | 96.2 % | none at input level | 529 |
| ceftriaxone | 6917 | 97.3 % | none at input level | 897 |
| chloramphenicol | 7049 | 97.1 % | none at input level | 425 |
| ciprofloxacin | 6935 | 97.0 % | none at input level | 606 |
| colistin | 525 | 93.5 % | constant trait | 0 (below 20) |
| gentamicin | 7047 | 97.1 % | none at input level | 1008 |
| imipenem | 132 | 81.8 % | constant trait | 0 (below 20) |
| kanamycin | 3206 | 96.1 % | none at input level | 427 |
| meropenem | 3300 | 96.5 % | constant trait | 0 (below 20) |
| nalidixic acid | 6915 | 97.3 % | none at input level | 584 |
| spectinomycin | 1 | not defined | fewer than 2 tested and typed isolates; fewer than 2 lineages with at least 2 isolates; constant trait | 0 (below 20) |
| streptomycin | 6520 | 96.9 % | none at input level | 2659 |
| sulfamethoxazole | 253 | 88.9 % | none at input level | 54 |
| sulfisoxazole | 6662 | 97.3 % | none at input level | 1927 |
| tetracycline | 7047 | 97.1 % | none at input level | 3338 |
| trimethoprim-sulfamethoxazole | 6915 | 97.3 % | none at input level | 208 |

## 3. The lineage share of the call, trait by trait

The interval is for the share of the lineages in this collection: the lineages and their sizes are held fixed and the isolates of each lineage are drawn again from the smoothed distribution of its calls. For **2 of the 18 traits** shown the interval reaches zero, so no lineage effect is distinguishable from none for that trait.

The control column is the same estimator run on shuffled lineage labels; it should sit near zero, and a share is read against it rather than against zero.

*Figure 2 (drawn on the page).* Lineage share of the call by trait, point estimate with 95 % interval for the represented lineages, 18 largest of 18. Traits whose interval reaches zero are drawn in grey and marked ‡.

The last column asks what the call says about the ordering it was cut from: the lineage share of the latent MIC, or of any continuous value the call thresholds, the rank intraclass correlation of that value. A call does not identify that share; it bounds it. The bracket gives the smallest share any ordering consistent with the calls allows, the lower bound the calls establish, and the largest; an upper end marked ≤ is a certified bound rather than a share attained. The figure after it is the one-sided 95 % lower confidence limit of the lower bound, from splitting the isolates of every lineage at random, choosing a direction on one half and testing it on the other: a share of the latent ordering that the lineages of this collection are shown to account for, with no model for the latent values beyond their agreeing with the readings. A MIC reading on a panel that contains the cut-off refines the call, and refining the readings can only narrow the bounds they allow.

**Table 3.** The share of every trait with its interval, and the bounds the call places on the share of the latent ordering.

| Trait | Share | 95 % interval, represented lineages | Latent ordering: bounds; 95 % lower limit |
|---|---:|---:|---:|
| nalidixic acid | 0.943 | 0.910 to 0.977 | 0.220 to ≤ 1.000; ≥ 0.211 |
| ciprofloxacin | 0.912 | 0.874 to 0.953 | 0.211 to ≤ 1.000; ≥ 0.192 |
| tetracycline | 0.686 | 0.667 to 0.706 | 0.293 to ≤ 1.000; ≥ 0.240 |
| sulfisoxazole | 0.646 | 0.627 to 0.668 | 0.248 to ≤ 1.000; ≥ 0.205 |
| streptomycin | 0.574 | 0.554 to 0.599 | 0.190 to ≤ 1.000; ≥ 0.138 |
| sulfamethoxazole | 0.529 | 0.388 to 0.694 | 0.152 to ≤ 1.000; ≥ 0.044 |
| chloramphenicol | 0.463 | 0.407 to 0.515 | 0.011 to ≤ 1.000; ≥ 0.000 |
| kanamycin | 0.408 | 0.365 to 0.455 | 0.045 to ≤ 1.000; ≥ 0.018 |
| ceftiofur | 0.403 | 0.371 to 0.442 | 0.060 to ≤ 1.000; ≥ 0.029 |
| ceftazidime | 0.399 | 0.213 to 0.651 | 0.085 to ≤ 1.000; ≥ 0.000 |
| cefoxitin | 0.398 | 0.366 to 0.430 | 0.051 to ≤ 1.000; ≥ 0.025 |
| ceftriaxone | 0.392 | 0.365 to 0.417 | 0.045 to ≤ 1.000; ≥ 0.024 |
| gentamicin | 0.371 | 0.341 to 0.399 | 0.023 to ≤ 1.000; ≥ 0.002 |
| amoxicillin-clavulanic acid | 0.357 | 0.334 to 0.378 | 0.048 to ≤ 1.000; ≥ 0.026 |
| ampicillin | 0.337 | 0.318 to 0.356 | 0.029 to ≤ 1.000; ≥ 0.010 |
| trimethoprim-sulfamethoxazole | 0.324 | 0.261 to 0.382 | 0.002 to ≤ 1.000; ≥ 0.000 |
| azithromycin | 0.128 | 0.000 to 0.408 | 0.001 to ≤ 1.000; ≥ 0.000 |
| amikacin | 0.000 | 0.000 to 1.000 | 0.000 to ≤ 1.000; ≥ 0.000 |

## 4. Evidence that survives re-reading

A p-value is a statement about one look at the data. A surveillance panel is looked at again every year, and a p-value recomputed each time loses its error control. The e-value is evidence on a scale made for that: this run's e-value is a statement about this collection, and a programme that adds an intake each year multiplies the e-value of each new intake, scored against the lineage rates learned from the earlier ones, into a running product whose error control holds at whatever intake it is read (sequential_e_process). The e-BH procedure controls the false-discovery rate across the panel whatever the dependence between traits. Larger is stronger; 1 is no evidence.

**Table 4.** e-value per trait and the e-BH selection at level 0.05. The selection threshold on this run is 25.9; a trait at or above it is selected.

| Trait | e-value | natural log | Selected |
|---|---:|---:|---|
| tetracycline | 1.1 × 10²⁵⁷ | 591.82 | yes |
| sulfisoxazole | 1.2 × 10²⁰⁴ | 469.92 | yes |
| streptomycin | 7.2 × 10¹⁹¹ | 441.77 | yes |
| nalidixic acid | 2.5 × 10¹⁷⁹ | 413.07 | yes |
| ciprofloxacin | 1.6 × 10¹⁷⁴ | 401.13 | yes |
| ampicillin | 1.5 × 10¹¹¹ | 256.01 | yes |
| ceftriaxone | 7.5 × 10¹⁰⁰ | 232.27 | yes |
| gentamicin | 9.9 × 10⁹⁸ | 227.95 | yes |
| amoxicillin-clavulanic acid | 2.7 × 10⁹⁷ | 224.33 | yes |
| chloramphenicol | 1.3 × 10⁸³ | 191.36 | yes |
| cefoxitin | 2.6 × 10⁸¹ | 187.45 | yes |
| ceftiofur | 1.4 × 10⁵⁴ | 124.66 | yes |
| trimethoprim-sulfamethoxazole | 3.5 × 10⁵¹ | 118.68 | yes |
| kanamycin | 1.9 × 10⁵¹ | 118.07 | yes |
| sulfamethoxazole | 2.4 × 10⁵ | 12.39 | yes |
| azithromycin | 1.0 × 10⁵ | 11.55 | yes |
| ceftazidime | 95.8 | 4.56 | yes |
| imipenem | 1.0 | 0.00 | no |
| colistin | 1.0 | 0.00 | no |
| meropenem | 1.0 | 0.00 | no |
| amikacin | 0.6 | -0.47 | no |
| spectinomycin | not computed | not computed | no |

17 of 22 traits are selected. Note: e-value in the betting sense of Vovk and Wang, not the BLAST expectation value and not the E-value of VanderWeele and Ding.

*Figure 3 (drawn on the page).* Evidence per trait on the natural-log scale, 18 largest of 21. The dashed rule is 1/α, the evidence one trait alone needs at level 0.05; the solid rule is the e-BH selection threshold on this run, 25.9, which rises with the number of traits read together, and a trait at or beyond it is selected. Traits the e-BH procedure selected are drawn in blue; the scale is logarithmic, so equal steps are equal factors of evidence.

## 5. Reading at the recorded resolution

No recorded dilutions were supplied on this run; the shares above are read from binary calls only.

## 6. A change of lineages or a change within them

A prevalence difference between two collections can be split into a change in lineage composition and a change in rate within lineages. That decomposition needs two collections and was not run here: this record holds one.

## 7. Provenance and terms

- **Software:** amr-clonalshare 1.0.0, record schema 1.0
- **Seed:** 42
- **Configuration:** `sha256:4890efe08a73`
- **Record:** `sha256:d36d68a5934104b98ea8d69f4974452ec886cd8ac788a9980b2759e9fd7cb241`

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

amr-clonalshare 1.0.0 · record schema 1.0 · seed 42 · configuration sha256:4890efe08a73
Cite this run as: “amr-clonalshare 1.0.0, run CE314965, record sha256:d36d68a59341.”
This report supersedes any earlier report bearing the same run identifier. It is regenerated from the record and holds no value that the record does not. The symbols carry the same wording in every run of this software; no symbol against a value means only that none of the listed conditions fired.
