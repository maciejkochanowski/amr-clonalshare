# amr-clonalshare run report

Every number below is read from `clonal_share_result.json`, written by the same run, or from the release's validation grid where the text says so; none is recomputed here.

- **Run:** 2FDB68DA
- **Issued:** 2026-09-25T15:23Z
- **Software:** amr-clonalshare 1.0.0
- **Isolates:** 7049
- **Lineage:** pds_cluster
- **Seed:** 42
- **Configuration:** sha256:958ab2b3ed48
- **Record digest:** sha256:61ce0b04b4d3c4c78cabd30def8c5bd5a3399f85bd862ae278e8a17831b4cff5

**22 antimicrobials; 88 analysis outcomes.** 63 computed, 0 full ranges, 25 unavailable or incomplete. Each result retains its own target and data subset.

## 1. Measurement summary

**Quantity measured.** The share of the variation in the recorded binary outcome across this collection that lineage membership accounts for, read from a lineage label and an interpreted result and scored on isolates the estimator did not see.

- **Cohort:** 7049 isolates in 534 lineages, typed by `pds_cluster`
- **Antimicrobials read:** 22
- **Membership shares estimable:** 18 of 22

**Gates.** 63 computed, 0 full ranges, 25 unavailable or incomplete. Each result retains its own target and data subset. A withheld value failed a declared reporting condition of its estimator; the recorded reason distinguishes input limitations from numerical failure. Passing a gate does not verify the model.

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
| amikacin | collection_membership | 2401 | computed | observed and typed subset |
| amoxicillin-clavulanic acid | collection_membership | 6917 | computed | observed and typed subset |
| ampicillin | collection_membership | 7049 | computed | observed and typed subset |
| azithromycin | collection_membership | 4779 | computed | observed and typed subset |
| cefoxitin | collection_membership | 6917 | computed | observed and typed subset |
| ceftazidime | collection_membership | 132 | computed | observed and typed subset |
| ceftiofur | collection_membership | 3616 | computed | observed and typed subset |
| ceftriaxone | collection_membership | 6917 | computed | observed and typed subset |
| chloramphenicol | collection_membership | 7049 | computed | observed and typed subset |
| ciprofloxacin | collection_membership | 6935 | computed | observed and typed subset |
| colistin | collection_membership | 525 | unavailable | observed and typed subset |
| gentamicin | collection_membership | 7047 | computed | observed and typed subset |
| imipenem | collection_membership | 132 | unavailable | observed and typed subset |
| kanamycin | collection_membership | 3206 | computed | observed and typed subset |
| meropenem | collection_membership | 3300 | unavailable | observed and typed subset |
| nalidixic acid | collection_membership | 6915 | computed | observed and typed subset |
| spectinomycin | collection_membership | 1 | unavailable | observed and typed subset |
| streptomycin | collection_membership | 6520 | computed | observed and typed subset |
| sulfamethoxazole | collection_membership | 253 | computed | observed and typed subset |
| sulfisoxazole | collection_membership | 6662 | computed | observed and typed subset |
| tetracycline | collection_membership | 7047 | computed | observed and typed subset |
| trimethoprim-sulfamethoxazole | collection_membership | 6915 | computed | observed and typed subset |
| amikacin | realised_component | 2401 | unavailable | observed and typed subset |
| amoxicillin-clavulanic acid | realised_component | 6917 | unavailable | observed and typed subset |
| ampicillin | realised_component | 7049 | computed | observed and typed subset |
| azithromycin | realised_component | 4779 | unavailable | observed and typed subset |
| cefoxitin | realised_component | 6917 | unavailable | observed and typed subset |
| ceftazidime | realised_component | 132 | computed | observed and typed subset |
| ceftiofur | realised_component | 3616 | unavailable | observed and typed subset |
| ceftriaxone | realised_component | 6917 | unavailable | observed and typed subset |
| chloramphenicol | realised_component | 7049 | unavailable | observed and typed subset |
| ciprofloxacin | realised_component | 6935 | unavailable | observed and typed subset |
| colistin | realised_component | 525 | unavailable | observed and typed subset |
| gentamicin | realised_component | 7047 | unavailable | observed and typed subset |
| imipenem | realised_component | 132 | unavailable | observed and typed subset |
| kanamycin | realised_component | 3206 | unavailable | observed and typed subset |
| meropenem | realised_component | 3300 | unavailable | observed and typed subset |
| nalidixic acid | realised_component | 6915 | unavailable | observed and typed subset |
| spectinomycin | realised_component | 1 | unavailable | observed and typed subset |
| streptomycin | realised_component | 6520 | unavailable | observed and typed subset |
| sulfamethoxazole | realised_component | 253 | unavailable | observed and typed subset |
| sulfisoxazole | realised_component | 6662 | unavailable | observed and typed subset |
| tetracycline | realised_component | 7047 | unavailable | observed and typed subset |
| trimethoprim-sulfamethoxazole | realised_component | 6915 | unavailable | observed and typed subset |
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

**Table 0b.** Reasons and next checks for results requiring interpretation. These checks do not change the recorded status or justify changing reporting gates.

| Agent | Analysis / status | Reason | Next check |
|---|---|---|---|
| colistin | collection_membership / unavailable | Constant retained outcome; no variation to attribute | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| imipenem | collection_membership / unavailable | Constant retained outcome; no variation to attribute | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| meropenem | collection_membership / unavailable | Constant retained outcome; no variation to attribute | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| spectinomycin | collection_membership / unavailable | Only one retained lineage; no between-lineage comparison | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| amikacin | realised_component / unavailable | the observed mean-square ratio 0.061 falls below the 0.025 point of the central F on 273 and 2127 degrees of freedom, so the inverted confidence set for the noncentrality is empty and the interval closes on zero | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| amoxicillin-clavulanic acid | realised_component / unavailable | within-lineage residuals have excess kurtosis 2.02 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| azithromycin | realised_component / unavailable | within-lineage residuals have excess kurtosis 217.42 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| cefoxitin | realised_component / unavailable | within-lineage residuals have excess kurtosis 7.68 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| ceftiofur | realised_component / unavailable | within-lineage residuals have excess kurtosis 3.49 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| ceftriaxone | realised_component / unavailable | within-lineage residuals have excess kurtosis 3.70 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| chloramphenicol | realised_component / unavailable | within-lineage residuals have excess kurtosis 11.52 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| ciprofloxacin | realised_component / unavailable | within-lineage residuals have excess kurtosis 98.84 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| colistin | realised_component / unavailable | every lineage is constant, so the within-lineage scale is zero and no share is identified | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| gentamicin | realised_component / unavailable | within-lineage residuals have excess kurtosis 3.35 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| imipenem | realised_component / unavailable | every lineage is constant, so the within-lineage scale is zero and no share is identified | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| kanamycin | realised_component / unavailable | within-lineage residuals have excess kurtosis 4.12 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| meropenem | realised_component / unavailable | every lineage is constant, so the within-lineage scale is zero and no share is identified | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| nalidixic acid | realised_component / unavailable | within-lineage residuals have excess kurtosis 170.19 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| spectinomycin | realised_component / unavailable | 1 populated lineage(s); at least 2 are needed for one between-lineage degree of freedom | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| streptomycin | realised_component / unavailable | within-lineage residuals have excess kurtosis 2.31 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| sulfamethoxazole | realised_component / unavailable | within-lineage residuals have excess kurtosis 3.86 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| sulfisoxazole | realised_component / unavailable | within-lineage residuals have excess kurtosis 4.51 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| tetracycline | realised_component / unavailable | within-lineage residuals have excess kurtosis 4.03 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
| trimethoprim-sulfamethoxazole | realised_component / unavailable | within-lineage residuals have excess kurtosis 16.09 against a limit of 0.99, above which the exact interval was measured to lose its level | Check readable calls, recorded labels, outcome variation and method-specific support in input QC and diagnostics. |
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

**Table 1.** Every trait, ordered by share. The reading in the last column is fixed by two conditions the record holds: whether the interval excludes zero and whether the e-BH procedure selected the trait at level 0.05.

| Trait | Share | 95 % interval | Control | e-value |  | Reading |
|---|---:|---:|---:|---:|---|---|
| nalidixic acid | 0.936 | 0.237 to 0.977 | -0.07 | 2.5 × 10¹⁷⁹ |  | evidence of a lineage effect |
| ciprofloxacin | 0.878 | 0.397 to 0.950 | -0.07 | 1.6 × 10¹⁷⁴ |  | evidence of a lineage effect |
| tetracycline | 0.667 | 0.577 to 0.740 | -0.07 | 1.1 × 10²⁵⁷ |  | evidence of a lineage effect |
| sulfisoxazole | 0.635 | 0.411 to 0.744 | -0.07 | 1.2 × 10²⁰⁴ |  | evidence of a lineage effect |
| streptomycin | 0.564 | 0.468 to 0.631 | -0.08 | 7.2 × 10¹⁹¹ |  | evidence of a lineage effect |
| sulfamethoxazole | 0.518 | 0.114 to 0.737 | -0.20 | 2.4 × 10⁵ |  | evidence of a lineage effect |
| chloramphenicol | 0.451 | 0.128 to 0.494 | -0.07 | 1.3 × 10⁸³ |  | evidence of a lineage effect |
| kanamycin | 0.432 | 0.282 to 0.566 | -0.10 | 1.9 × 10⁵¹ |  | evidence of a lineage effect |
| ceftazidime | 0.411 | -0.254 to 0.530 | -0.28 | 95.8 | ‡ | lineage effect selected by e-BH; its size is not resolved |
| ceftiofur | 0.404 | 0.250 to 0.579 | -0.09 | 1.4 × 10⁵⁴ |  | evidence of a lineage effect |
| ceftriaxone | 0.394 | 0.278 to 0.527 | -0.07 | 7.5 × 10¹⁰⁰ |  | evidence of a lineage effect |
| cefoxitin | 0.385 | 0.226 to 0.553 | -0.07 | 2.6 × 10⁸¹ |  | evidence of a lineage effect |
| gentamicin | 0.372 | 0.267 to 0.432 | -0.07 | 9.9 × 10⁹⁸ |  | evidence of a lineage effect |
| amoxicillin-clavulanic acid | 0.364 | 0.254 to 0.489 | -0.07 | 2.7 × 10⁹⁷ |  | evidence of a lineage effect |
| ampicillin | 0.349 | 0.280 to 0.448 | -0.07 | 1.5 × 10¹¹¹ |  | evidence of a lineage effect |
| trimethoprim-sulfamethoxazole | 0.341 | -0.018 to 0.369 | -0.07 | 3.5 × 10⁵¹ | ‡ | lineage effect selected by e-BH; its size is not resolved |
| azithromycin | 0.148 | -0.003 to 0.471 | -0.08 | 1.0 × 10⁵ | ‡ | lineage effect selected by e-BH; its size is not resolved |
| amikacin | 0.074 | 0.072 to 0.075 | -0.09 | 0.6 |  | interval excludes zero; not selected by e-BH |

**Table 1b.** Lineage structure behind each share: lineages with at least two isolates, the effective number of lineages (inverse of the sum of squared lineage shares) and the share of isolates in lineages with at least two members (support). A high support with few effective lineages means the share rests on a few lineages.

| Trait | Lineages | Repeated | Effective | Support |
|---|---:|---:|---:|---:|
| nalidixic acid | 496 | 307 | 35.9 | 97.3 % |
| ciprofloxacin | 531 | 325 | 37.7 | 97.0 % |
| tetracycline | 534 | 327 | 36.3 | 97.1 % |
| sulfisoxazole | 477 | 297 | 35.0 | 97.3 % |
| streptomycin | 513 | 314 | 35.7 | 96.9 % |
| sulfamethoxazole | 61 | 33 | 19.2 | 88.9 % |
| chloramphenicol | 534 | 327 | 36.3 | 97.1 % |
| kanamycin | 310 | 186 | 25.1 | 96.1 % |
| ceftazidime | 47 | 23 | 10.8 | 81.8 % |
| ceftiofur | 336 | 200 | 26.9 | 96.2 % |
| ceftriaxone | 496 | 307 | 36.0 | 97.3 % |
| cefoxitin | 496 | 307 | 36.0 | 97.3 % |
| gentamicin | 534 | 327 | 36.3 | 97.1 % |
| amoxicillin-clavulanic acid | 496 | 307 | 36.0 | 97.3 % |
| ampicillin | 534 | 327 | 36.3 | 97.1 % |
| trimethoprim-sulfamethoxazole | 496 | 307 | 35.9 | 97.3 % |
| azithromycin | 408 | 242 | 28.7 | 96.5 % |
| amikacin | 274 | 163 | 26.1 | 95.4 % |

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
> For chloramphenicol, kanamycin, ceftiofur, ceftriaxone, cefoxitin, gentamicin, amoxicillin-clavulanic acid and ampicillin, a lineage association is detected and the estimated share is below one half. Most variation is not explained by the lineage labels at this typing resolution.
>
> For ceftazidime, trimethoprim-sulfamethoxazole and azithromycin, the e-value selects a lineage effect while the interval for its size still reaches zero: there is evidence that positive outcomes are not spread evenly across the lineages, and this collection is too small, or too uneven across its lineages, to say how much of it the lineages carry. Read the e-value as the finding and the share as not yet resolved.
>
> For amikacin, the interval for the share excludes zero while the e-value does not select a lineage effect at the panel's false-discovery level: the two readings disagree, and the e-value, which carries the multiplicity correction, is the one to read; treat the share as suggestive.
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

Lineage groups in the input: 534, of which 207 hold a single isolate. Support 97.1 % against the 80.0 % the estimator needs: accepted.

*Figure 1 (drawn on the page).* Isolates per lineage, largest first, the 40 largest of 534. A lineage of one isolate, drawn in orange, cannot be predicted out of sample and counts against support; support is the share of isolates in the other lineages, 97.1 % here against the 80.0 % the estimator requires. The largest lineage holds 7.5 % of the isolates, which sets how much one lineage can weigh in the share.

**Table 2.** Conditions the estimator requires before any result is reported.

| Condition | Observed | Required | Verdict |
|---|---:|---:|---|
| Lineage support, lowest over the antimicrobials read | 81.8 % | ≥ 80.0 % | accepted |
| Lineage groups used | 534 | reported | accepted |

Each row uses only isolates with both a readable result for that agent and a recorded lineage. The last column is an algebraic support scenario: one genuinely new tested isolate in each of that many distinct singleton lineages. It does not guarantee interval precision, coverage, or overall estimability; it cannot repair a constant trait or a missing lineage contrast. Duplicating existing rows adds no evidence. The rarer-outcome count is a warning threshold, not an additional estimator gate.

**Table 2a.** Per-agent input feasibility and support scenario.

| Agent | Retained | Support | Input failures | Rarer outcome | New isolates for support |
|---|---|---|---|---|---|
| amikacin | 2401 | 95.4 % | none at input level | 1 (below 20) | 0 |
| amoxicillin-clavulanic acid | 6917 | 97.3 % | none at input level | 1163 | 0 |
| ampicillin | 7049 | 97.1 % | none at input level | 1864 | 0 |
| azithromycin | 4779 | 96.5 % | none at input level | 24 | 0 |
| cefoxitin | 6917 | 97.3 % | none at input level | 649 | 0 |
| ceftazidime | 132 | 81.8 % | none at input level | 37 | 0 |
| ceftiofur | 3616 | 96.2 % | none at input level | 529 | 0 |
| ceftriaxone | 6917 | 97.3 % | none at input level | 897 | 0 |
| chloramphenicol | 7049 | 97.1 % | none at input level | 425 | 0 |
| ciprofloxacin | 6935 | 97.0 % | none at input level | 606 | 0 |
| colistin | 525 | 93.5 % | constant trait | 0 (below 20) | 0 |
| gentamicin | 7047 | 97.1 % | none at input level | 1008 | 0 |
| imipenem | 132 | 81.8 % | constant trait | 0 (below 20) | 0 |
| kanamycin | 3206 | 96.1 % | none at input level | 427 | 0 |
| meropenem | 3300 | 96.5 % | constant trait | 0 (below 20) | 0 |
| nalidixic acid | 6915 | 97.3 % | none at input level | 584 | 0 |
| spectinomycin | 1 | not defined | fewer than 2 tested and typed isolates; fewer than 2 lineages; constant trait | 0 (below 20) | 1 |
| streptomycin | 6520 | 96.9 % | none at input level | 2659 | 0 |
| sulfamethoxazole | 253 | 88.9 % | none at input level | 54 | 0 |
| sulfisoxazole | 6662 | 97.3 % | none at input level | 1927 | 0 |
| tetracycline | 7047 | 97.1 % | none at input level | 3338 | 0 |
| trimethoprim-sulfamethoxazole | 6915 | 97.3 % | none at input level | 208 | 0 |

## 3. The clonal share, trait by trait

The intervals have a nominal 95 % level; their measured coverage depends on the cohort design and model assumptions. For **3 of the 18 traits** shown the interval includes zero, so no lineage effect is distinguishable from none for that trait.

The control column is the same estimator run on shuffled lineage labels; it should sit near zero, and a share is read against it rather than against zero.

*Figure 2 (drawn on the page).* Clonal share by trait, point estimate with 95 % interval, 18 largest of 18. Traits whose interval crosses zero are drawn in grey and marked ‡. The thin line beneath each interval is the species interval of Table 3, for a fresh draw of lineages; it is floored at zero by construction. The intervals are drawn as computed. The quantity lies between 0 and 1, so the part of an interval below zero carries no information: reading each interval as its overlap with that range leaves the coverage unchanged, and a lower limit at or below zero means the same thing either way.

On the release's validation grid, the interval this estimator prints for a binary trait contained the truth in 0.970 of 28,000 runs over 70 simulated cohorts (by cohort, 0.825 to 1.000). That figure belongs to the release, not to this run; it is what the phrase "95 % interval" was measured to mean.

Two intervals answer two questions. The interval above is for the share the lineages in this collection carry. The second interval below is for the share a fresh draw of lineages from the species would show, stated on the scale of the realised share of Table 3b below; it adds the sampling of the lineages themselves, on 533 degrees of freedom, is widened to the envelope of the first taken on that scale (with a lower end no smaller than zero, since a species share is not negative), and is the one to quote when the figure is read as a property of the species rather than of this collection. With few lineages it is markedly wider; with many the two nearly coincide.

For 6 traits one lineage carries more than half of the between-lineage variation (nalidixic acid 81.7 %; ciprofloxacin 70.9 %; chloramphenicol 66.7 %; ceftazidime 51.9 %; trimethoprim-sulfamethoxazole 79.2 %; amikacin 94.9 %). The species interval describes lineage effects drawn from one law, and a collection in which one lineage carries the positive calls is not that: on the validation grid a carrier law of this kind took the species interval below its level with few lineages while the interval for the lineages in hand held. For them read the first interval as the statement about this collection and the species interval with that reservation.

The last column is an exploratory, model-equivalent restatement under a Gaussian probit threshold model. It transforms the collection's observed-scale share and interval at the estimated prevalence, treated as fixed. It does not estimate the actual latent variance of these particular lineages or a generic binomial mixed-model intraclass correlation. Monotonicity preserves coverage only for a matching target at a known, fixed prevalence; it does not guarantee coverage after estimating prevalence. Read the observed-scale interval as primary. Values at or below zero are displayed at the zero boundary.

**Table 3.** The two intervals for every trait: for lineage membership in this collection and for the species, with the share on the latent scale.

| Trait | Share | 95 % interval, these lineages | Species share (95 % interval) | Latent restatement (exploratory interval) |
|---|---:|---:|---:|---:|
| nalidixic acid | 0.936 | 0.237 to 0.977 | 0.938 (0.243 to 0.978) | 0.997 (0.500 to 1.000) |
| ciprofloxacin | 0.878 | 0.397 to 0.950 | 0.881 (0.403 to 0.951) | 0.988 (0.697 to 0.998) |
| tetracycline | 0.667 | 0.577 to 0.740 | 0.673 (0.584 to 0.745) | 0.866 (0.787 to 0.918) |
| sulfisoxazole | 0.635 | 0.411 to 0.744 | 0.642 (0.418 to 0.749) | 0.852 (0.623 to 0.926) |
| streptomycin | 0.564 | 0.468 to 0.631 | 0.571 (0.475 to 0.638) | 0.777 (0.674 to 0.839) |
| sulfamethoxazole | 0.518 | 0.114 to 0.737 | 0.532 (0.120 to 0.747) | 0.761 (0.212 to 0.928) |
| chloramphenicol | 0.451 | 0.128 to 0.494 | 0.458 (0.131 to 0.509) | 0.775 (0.346 to 0.810) |
| kanamycin | 0.432 | 0.282 to 0.566 | 0.442 (0.290 to 0.575) | 0.704 (0.520 to 0.828) |
| ceftazidime | 0.411 | -0.254 to 0.530 | 0.434 (0.000 to 0.629) | 0.625 (0.000 to 0.758) |
| ceftiofur | 0.404 | 0.250 to 0.579 | 0.413 (0.257 to 0.588) | 0.665 (0.464 to 0.834) |
| ceftriaxone | 0.394 | 0.278 to 0.527 | 0.401 (0.284 to 0.534) | 0.664 (0.517 to 0.797) |
| cefoxitin | 0.385 | 0.226 to 0.553 | 0.391 (0.231 to 0.560) | 0.679 (0.472 to 0.835) |
| gentamicin | 0.372 | 0.267 to 0.432 | 0.379 (0.273 to 0.439) | 0.630 (0.492 to 0.699) |
| amoxicillin-clavulanic acid | 0.364 | 0.254 to 0.489 | 0.370 (0.259 to 0.496) | 0.606 (0.457 to 0.746) |
| ampicillin | 0.349 | 0.280 to 0.448 | 0.355 (0.285 to 0.455) | 0.550 (0.456 to 0.673) |
| trimethoprim-sulfamethoxazole | 0.341 | -0.018 to 0.369 | 0.348 (0.000 to 0.401) | 0.716 (0.000 to 0.743) |
| azithromycin | 0.148 | -0.003 to 0.471 | 0.152 (0.000 to 0.480) | 0.599 (0.000 to 0.884) |
| amikacin | 0.074 | 0.072 to 0.075 | 0.077 (0.055 to 0.097) | 0.601 (0.598 to 0.604) |

**The realised share of the same call.** The variance-component ratio of the lineages in hand, with an interval that is exact under a Gaussian within-lineage law. Its gate reads the excess kurtosis of the within-lineage residuals and withholds the interval above 0.99. Rare binary calls often fail this check; prevalence alone does not determine the residual kurtosis. The interval is approximate on binary data even when the gate opens; the estimate is printed either way, the interval only where the gate opened, for 2 of 22 traits here.

**Table 3b.** The realised share per trait and the verdict of its kurtosis gate.

| Trait | Realised share | Gaussian-model 95 % interval | Residual excess kurtosis | Gate |
|---|---:|---:|---:|---|
| nalidixic acid | 0.942 | withheld | 170.19 | closed |
| ciprofloxacin | 0.915 | withheld | 98.84 | closed |
| tetracycline | 0.689 | withheld | 4.03 | closed |
| sulfisoxazole | 0.646 | withheld | 4.51 | closed |
| streptomycin | 0.574 | withheld | 2.31 | closed |
| sulfamethoxazole | 0.508 | withheld | 3.86 | closed |
| chloramphenicol | 0.463 | withheld | 11.52 | closed |
| cefoxitin | 0.406 | withheld | 7.68 | closed |
| ceftiofur | 0.399 | withheld | 3.49 | closed |
| ceftriaxone | 0.385 | withheld | 3.70 | closed |
| kanamycin | 0.378 | withheld | 4.12 | closed |
| gentamicin | 0.359 | withheld | 3.35 | closed |
| amoxicillin-clavulanic acid | 0.346 | withheld | 2.02 | closed |
| ampicillin | 0.322 | 0.304 to 0.340 | 0.04 | open |
| ceftazidime | 0.312 | 0.117 to 0.472 | 0.92 | open |
| trimethoprim-sulfamethoxazole | 0.298 | withheld | 16.09 | closed |
| azithromycin | 0.185 | withheld | 217.42 | closed |
| amikacin | 0.000 | withheld | 2254.54 | closed |
| colistin | not computed | withheld | not computed | closed |
| imipenem | not computed | withheld | not computed | closed |
| meropenem | not computed | withheld | not computed | closed |
| spectinomycin | not computed | withheld | not computed | closed |

## 4. Evidence that survives re-reading

A p-value is a statement about one look at the data. A surveillance panel is looked at again every year, and a p-value recomputed each time loses its guarantee. The e-value is evidence on a scale made for that: this run's e-value is a statement about this cohort, and a programme that adds an intake each year multiplies the e-value of each new intake, scored against the lineage rates learned from the earlier ones, into a running product whose guarantee holds at whatever intake it is read (sequential_e_process). The e-BH procedure controls the false-discovery rate across the panel whatever the dependence between traits. Larger is stronger; 1 is no evidence.

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
| imipenem | 1.0 | -0.00 | no |
| colistin | 1.0 | -0.00 | no |
| meropenem | 1.0 | -0.00 | no |
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
- **Configuration:** `sha256:958ab2b3ed48`
- **Record:** `sha256:61ce0b04b4d3c4c78cabd30def8c5bd5a3399f85bd862ae278e8a17831b4cff5`

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

amr-clonalshare 1.0.0 · record schema 1.0 · seed 42 · configuration sha256:958ab2b3ed48
Cite this run as: “amr-clonalshare 1.0.0, run 2FDB68DA, record sha256:61ce0b04b4d3.”
This report supersedes any earlier report bearing the same run identifier. It is regenerated from the record and holds no value that the record does not. The symbols carry the same wording in every run of this software; no symbol against a value means only that none of the listed conditions fired.
Classical collection-bootstrap endpoints are printed without clipping; the separately labeled species interval is floored at zero as part of its construction.
