# amr-clonalshare run report

Every number below is read from `clonal_share_result.json`, written by the same run; none is recomputed here.

- **Run:** ADDB37A7
- **Issued:** 2026-09-11T10:23Z
- **Software:** amr-clonalshare 1.0.0
- **Isolates:** 7049
- **Lineage:** pds_cluster
- **Seed:** 42
- **Configuration:** sha256:72686f22ceea
- **Record digest:** sha256:6fe3209471d9bd8edf5ecaa60d4684dec4062942f9bc417a6024c01380742c1e

**16 of 22 antimicrobials estimable.** The estimator withheld 2 of 22; the condition that failed is named per antimicrobial. 4 of 22 could not be scored, for a reason named below.

## 1. Measurement summary

**Quantity measured.** The share of the variation in resistance across this collection that lineage membership accounts for, read from a lineage label and a susceptibility result and scored on isolates the estimator did not see.

- **Cohort:** 7049 isolates in 534 lineages, typed by `pds_cluster`
- **Antimicrobials read:** 22
- **Estimable:** 16 of 22

**Gates.** The estimator withheld 2 of 22; the condition that failed is named per antimicrobial. 4 of 22 could not be scored, for a reason named below. A withheld value is a completed reading: it says the collection cannot identify the quantity at this typing resolution. A value that clears the gate is not thereby correct.

**Table 1.** Every trait, ordered by share. The reading in the last column is fixed by two conditions the record holds: whether the interval excludes zero and whether the e-BH procedure selected the trait at level 0.05.

| Trait | Share | 95 % interval | Control | e-value |  | Reading |
|---|---:|---:|---:|---:|---|---|
| nalidixic acid | 0.937 | 0.233 to 0.975 | -0.07 | 2.7 × 10¹⁷⁶ |  | evidence of a lineage effect |
| ciprofloxacin | 0.878 | 0.392 to 0.950 | -0.07 | 8.8 × 10¹⁶⁶ |  | evidence of a lineage effect |
| tetracycline | 0.668 | 0.559 to 0.735 | -0.07 | 3.5 × 10²⁶¹ |  | evidence of a lineage effect |
| sulfisoxazole | 0.635 | 0.393 to 0.761 | -0.07 | 2.4 × 10²⁰⁸ |  | evidence of a lineage effect |
| streptomycin | 0.563 | 0.474 to 0.642 | -0.08 | 3.2 × 10¹⁸⁷ |  | evidence of a lineage effect |
| sulfamethoxazole | 0.515 | 0.149 to 0.739 | -0.22 | 2.8 × 10⁵ | † | refused |
| chloramphenicol | 0.451 | 0.130 to 0.499 | -0.07 | 7.6 × 10⁸⁵ |  | evidence of a lineage effect |
| kanamycin | 0.430 | 0.265 to 0.606 | -0.09 | 3.3 × 10⁴⁹ |  | evidence of a lineage effect |
| ceftiofur | 0.404 | 0.261 to 0.570 | -0.09 | 9.0 × 10⁵⁵ |  | evidence of a lineage effect |
| ceftazidime | 0.403 | -0.379 to 0.537 | -0.27 | 84.3 | ‡† | refused |
| ceftriaxone | 0.395 | 0.300 to 0.535 | -0.07 | 6.0 × 10¹⁰² |  | evidence of a lineage effect |
| cefoxitin | 0.386 | 0.240 to 0.556 | -0.07 | 3.6 × 10⁸¹ |  | evidence of a lineage effect |
| gentamicin | 0.373 | 0.246 to 0.432 | -0.07 | 1.3 × 10¹⁰⁹ |  | evidence of a lineage effect |
| amoxicillin-clavulanic acid | 0.367 | 0.259 to 0.487 | -0.07 | 2.9 × 10¹⁰⁴ |  | evidence of a lineage effect |
| ampicillin | 0.347 | 0.275 to 0.431 | -0.07 | 8.5 × 10¹¹² |  | evidence of a lineage effect |
| trimethoprim-sulfamethoxazole | 0.343 | -0.033 to 0.371 | -0.07 | 9.7 × 10⁵⁶ | ‡ | lineage effect selected by e-BH; its size is not resolved |
| azithromycin | 0.145 | 0.041 to 0.474 | -0.08 | 2.9 × 10⁴ |  | evidence of a lineage effect |
| amikacin | 0.074 | 0.072 to 0.075 | -0.09 | 0.6 |  | interval excludes zero; not selected by e-BH |

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
> For nalidixic acid, ciprofloxacin, tetracycline, sulfisoxazole and streptomycin, resistance travels with lineage: the share is at or above one half, so it is carried by particular clones and will move when those clones move. Measures that act on the clone, such as movement control and clone-directed surveillance, act on these agents.
>
> For chloramphenicol, kanamycin, ceftiofur, ceftriaxone, cefoxitin, gentamicin, amoxicillin-clavulanic acid, ampicillin and azithromycin, a lineage effect is detected but the share is below one half: most of the variation in resistance sits within lineages, so the resistance moves largely independently of the clone, and selection pressure, dosing and choice of agent are the levers that reach most of it.
>
> For trimethoprim-sulfamethoxazole, the e-value selects a lineage effect while the interval for its size still reaches zero: there is evidence that resistance is not spread evenly across the lineages, and this collection is too small, or too uneven across its lineages, to say how much of it the lineages carry. Read the e-value as the finding and the share as not yet resolved.
>
> For amikacin, the interval for the share excludes zero while the e-value does not select a lineage effect at the panel's false-discovery level: the two readings disagree, and the e-value, which carries the multiplicity correction, is the one to read; treat the share as suggestive.
>
> For sulfamethoxazole and ceftazidime, the estimator refused: the collection cannot identify the quantity. No reading follows.
>
> These sentences follow from the table by fixed rules and carry no judgement beyond it.
>

## 2. Admissibility of the input

- **Lineage column:** `pds_cluster`
- **Susceptibility calls:** 22 antimicrobials on 7049 isolates
- **Resampling:** 5 folds, 10 repeats, 300 bootstrap draws, 150 permutations per antimicrobial

Antimicrobials with fewer than 20 isolates of the rarer outcome: 1 of 22. Below that count the estimate is still computed and its interval says how little it rests on.

Lineage groups in the input: 534, of which 207 hold a single isolate. Support 97.1 % against the 90.0 % the estimator needs: accepted.

*Figure 1 (drawn on the page).* Isolates per lineage, largest first, the 40 largest of 534. A lineage of one isolate, drawn in orange, cannot be predicted out of sample and counts against support; support is the share of isolates in the other lineages, 97.1 % here against the 90.0 % the estimator requires. The largest lineage holds 7.5 % of the isolates, which sets how much one lineage can weigh in the share.

**Table 2.** Conditions the estimator requires before any result is reported.

| Condition | Observed | Required | Verdict |
|---|---:|---:|---|
| Lineage support, lowest over the antimicrobials read | 81.8 % | ≥ 90.0 % | refused |
| Lineage groups used | 534 | reported | accepted |

## 3. The clonal share, trait by trait

Read the intervals as frequencies. If cohorts like this one were drawn again and again, the interval printed for a trait would cover that trait's true share on about 95 draws in 100. For **2 of the 18 traits** shown the interval includes zero, so no lineage effect is distinguishable from none for that trait.

The control column is the same estimator run on shuffled lineage labels; it should sit near zero, and a share is read against it rather than against zero.

*Figure 2 (drawn on the page).* Clonal share by trait, point estimate with 95 % interval, 18 largest of 18. Traits whose interval crosses zero are drawn in grey and marked ‡. The thin line beneath each interval is the species interval of Table 3, for a fresh draw of lineages; it is floored at zero by construction. The intervals are drawn as computed. The quantity lies between 0 and 1, so the part of an interval below zero carries no information: reading each interval as its overlap with that range leaves the coverage unchanged, and a lower limit at or below zero means the same thing either way.

On the release's validation grid, the interval this estimator prints for a binary trait contained the truth in 0.970 of 28,000 runs over 70 simulated cohorts (by cohort, 0.825 to 1.000). That figure belongs to the release, not to this run; it is what the phrase "95 % interval" was measured to mean.

Two intervals answer two questions. The interval above is for the share the lineages in this collection carry. The second interval below is for the share a fresh draw of lineages from the species would show, stated on the scale of the realised share of Table 3b below; it adds the sampling of the lineages themselves, on 533 degrees of freedom, is widened to the envelope of the first taken on that scale (with a lower end no smaller than zero, since a species share is not negative), and is the one to quote when the figure is read as a property of the species rather than of this collection. With few lineages it is markedly wider; with many the two nearly coincide.

For 6 traits one lineage carries more than half of the between-lineage variation (nalidixic acid 81.7 %; ciprofloxacin 70.9 %; chloramphenicol 66.7 %; ceftazidime 51.9 %; trimethoprim-sulfamethoxazole 79.2 %; amikacin 94.9 %). The species interval describes lineage effects drawn from one law, and a collection in which one lineage carries the resistance is not that: on the validation grid a carrier law of this kind took the species interval below its level with few lineages while the interval for the lineages in hand held. For them read the first interval as the statement about this collection and the species interval with that reservation.

The last column restates the share on the latent scale of a threshold model, the scale on which a mixed model with a binomial link reports an intraclass correlation. The two scales are one-to-one at a given prevalence and the map is monotone, so the interval keeps its coverage; the latent figure is larger because a call discards the part of the liability that does not cross the threshold. A share at or below zero has no latent counterpart and is printed as zero.

**Table 3.** The two intervals for every trait: for lineage membership in this collection and for the species, with the share on the latent scale.

| Trait | Share | 95 % interval, these lineages | 95 % interval, species | Latent scale (95 % interval) |
|---|---:|---:|---:|---:|
| nalidixic acid | 0.937 | 0.233 to 0.975 | 0.238 to 0.976 | 0.997 (0.493 to 1.000) |
| ciprofloxacin | 0.878 | 0.392 to 0.950 | 0.398 to 0.951 | 0.988 (0.692 to 0.998) |
| tetracycline | 0.668 | 0.559 to 0.735 | 0.566 to 0.740 | 0.867 (0.770 to 0.915) |
| sulfisoxazole | 0.635 | 0.393 to 0.761 | 0.400 to 0.766 | 0.851 (0.601 to 0.935) |
| streptomycin | 0.563 | 0.474 to 0.642 | 0.481 to 0.648 | 0.776 (0.681 to 0.848) |
| sulfamethoxazole | 0.515 | 0.149 to 0.739 | 0.156 to 0.749 | 0.757 (0.270 to 0.929) |
| chloramphenicol | 0.451 | 0.130 to 0.499 | 0.134 to 0.511 | 0.775 (0.351 to 0.815) |
| kanamycin | 0.430 | 0.265 to 0.606 | 0.273 to 0.616 | 0.702 (0.495 to 0.859) |
| ceftiofur | 0.404 | 0.261 to 0.570 | 0.268 to 0.579 | 0.665 (0.480 to 0.827) |
| ceftazidime | 0.403 | -0.379 to 0.537 | 0.000 to 0.627 | 0.616 (0.000 to 0.764) |
| ceftriaxone | 0.395 | 0.300 to 0.535 | 0.306 to 0.542 | 0.665 (0.547 to 0.804) |
| cefoxitin | 0.386 | 0.240 to 0.556 | 0.245 to 0.563 | 0.681 (0.492 to 0.837) |
| gentamicin | 0.373 | 0.246 to 0.432 | 0.251 to 0.439 | 0.631 (0.460 to 0.698) |
| amoxicillin-clavulanic acid | 0.367 | 0.259 to 0.487 | 0.264 to 0.494 | 0.610 (0.464 to 0.745) |
| ampicillin | 0.347 | 0.275 to 0.431 | 0.280 to 0.438 | 0.549 (0.448 to 0.653) |
| trimethoprim-sulfamethoxazole | 0.343 | -0.033 to 0.371 | 0.000 to 0.407 | 0.718 (0.000 to 0.744) |
| azithromycin | 0.145 | 0.041 to 0.474 | 0.042 to 0.483 | 0.595 (0.338 to 0.885) |
| amikacin | 0.074 | 0.072 to 0.075 | 0.055 to 0.097 | 0.601 (0.598 to 0.604) |

**The realised share of the same call.** The variance-component ratio of the lineages in hand, with an interval that is exact under a Gaussian within-lineage law. Its gate reads the excess kurtosis of the within-lineage residuals and withholds the interval above 0.99, which on a binary call closes wherever the prevalence is far from one half; the estimate is printed either way, the interval only where the gate opened, for 2 of 22 traits here.

**Table 3b.** The realised share per trait and the verdict of its kurtosis gate.

| Trait | Realised share | Exact 95 % interval | Residual excess kurtosis | Gate |
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
| tetracycline | 3.5 × 10²⁶¹ | 602.21 | yes |
| sulfisoxazole | 2.4 × 10²⁰⁸ | 479.80 | yes |
| streptomycin | 3.2 × 10¹⁸⁷ | 431.74 | yes |
| nalidixic acid | 2.7 × 10¹⁷⁶ | 406.26 | yes |
| ciprofloxacin | 8.8 × 10¹⁶⁶ | 384.40 | yes |
| ampicillin | 8.5 × 10¹¹² | 260.03 | yes |
| gentamicin | 1.3 × 10¹⁰⁹ | 251.23 | yes |
| amoxicillin-clavulanic acid | 2.9 × 10¹⁰⁴ | 240.52 | yes |
| ceftriaxone | 6.0 × 10¹⁰² | 236.65 | yes |
| chloramphenicol | 7.6 × 10⁸⁵ | 197.75 | yes |
| cefoxitin | 3.6 × 10⁸¹ | 187.80 | yes |
| trimethoprim-sulfamethoxazole | 9.7 × 10⁵⁶ | 131.22 | yes |
| ceftiofur | 9.0 × 10⁵⁵ | 128.84 | yes |
| kanamycin | 3.3 × 10⁴⁹ | 114.02 | yes |
| sulfamethoxazole | 2.8 × 10⁵ | 12.55 | yes |
| azithromycin | 2.9 × 10⁴ | 10.29 | yes |
| ceftazidime | 84.3 | 4.43 | yes |
| colistin | 1.0 | 0.00 | no |
| imipenem | 1.0 | 0.00 | no |
| meropenem | 1.0 | 0.00 | no |
| amikacin | 0.6 | -0.47 | no |
| spectinomycin | not computed | not computed | no |

17 of 22 traits are selected. Note: e-value in the betting sense of Vovk and Wang, not the BLAST expectation value and not the E-value of VanderWeele and Ding.

*Figure 3 (drawn on the page).* Evidence per trait on the natural-log scale, 18 largest of 21. The dashed rule is 1/α, the evidence one trait alone needs at level 0.05; the solid rule is the e-BH selection threshold on this run, 25.9, which rises with the number of traits read together, and a trait at or beyond it is selected. Traits the e-BH procedure selected are drawn in blue; the scale is logarithmic, so equal steps are equal factors of evidence.

## 5. Reading at the recorded resolution

No recorded dilutions were supplied on this run; the shares above are read from binary calls only.

## 6. How resistance is carried across lineages

Of the **19** traits with a carriage reading, **17** depart from carriage in proportion to lineage size far enough to matter for a prevalence reading: 17 in fewer lineages than chance gives. The remaining 2 are consistent with proportional carriage.

The widest gap between the per-isolate and the per-lineage prevalence is -0.183, on `tetracycline`. Where the two differ, a small number of large lineages carry most of the resistance, and a change in prevalence may be a change in which lineages were sampled rather than a change in rate.

Two prevalences are reported for every trait, and they are two quantities rather than a biased and an unbiased one: per isolate is the clinical burden of the collection, per lineage is the share of its diversity that carries the trait. They separate whenever sampling across lineages is uneven. Beside them the effective number of carrying lineages, the reciprocal of a Herfindahl index over the carriers, is compared with the number that carriage in proportion to lineage size would give, by permutation. A trait departs from proportional carriage where that permutation p-value falls below 0.05; the direction says whether the carriers sit in fewer lineages than chance gives, in more, or in as many but in uneven shares across them.

**Table 6.** Prevalence on two scales and the concentration of carriage, per trait. The per-lineage interval is a two-stage cluster bootstrap over lineages and then isolates; the null interval for the effective number is the permutation reference for carriage in proportion to lineage size.

| Trait | Per isolate | Per lineage (95 % interval) | Carrying lineages | Effective number (null 95 %) | Carriage |
|---|---:|---:|---:|---:|---|
| colistin | 100.0 % | 100.0 % (100.0 % to 100.0 %) | 86 of 86 | 13.8 (13.8 to 13.8) | proportional |
| tetracycline | 52.6 % | 34.3 % (30.4 % to 38.1 %) | 236 of 534 | 17.3 (34.8 to 37.5) | concentrated in fewer lineages |
| streptomycin | 40.8 % | 27.9 % (24.2 % to 31.5 %) | 208 of 513 | 18.8 (33.7 to 37.4) | concentrated in fewer lineages |
| sulfisoxazole | 28.9 % | 23.5 % (19.9 % to 27.1 %) | 165 of 477 | 8.4 (32.4 to 36.9) | concentrated in fewer lineages |
| ceftazidime | 28.0 % | 10.6 % (3.9 % to 18.6 %) | 10 of 47 | 2.0 (5.4 to 15.7) | concentrated in fewer lineages |
| ampicillin | 26.4 % | 18.8 % (15.7 % to 21.9 %) | 172 of 534 | 14.4 (33.4 to 38.5) | concentrated in fewer lineages |
| sulfamethoxazole | 21.3 % | 17.0 % (8.6 % to 26.2 %) | 15 of 61 | 4.4 (10.8 to 21.1) | concentrated in fewer lineages |
| amoxicillin-clavulanic acid | 16.8 % | 15.2 % (12.4 % to 18.2 %) | 130 of 496 | 12.5 (32.2 to 38.6) | concentrated in fewer lineages |
| ceftiofur | 14.6 % | 16.4 % (12.6 % to 20.2 %) | 87 of 336 | 6.9 (22.3 to 29.7) | concentrated in fewer lineages |
| gentamicin | 14.3 % | 9.4 % (7.3 % to 11.7 %) | 120 of 534 | 8.3 (31.7 to 39.1) | concentrated in fewer lineages |
| kanamycin | 13.3 % | 7.0 % (4.6 % to 9.6 %) | 47 of 310 | 6.5 (20.4 to 28.0) | concentrated in fewer lineages |
| ceftriaxone | 13.0 % | 13.6 % (10.8 % to 16.5 %) | 116 of 496 | 7.7 (31.3 to 38.6) | concentrated in fewer lineages |
| cefoxitin | 9.4 % | 13.4 % (10.7 % to 16.3 %) | 117 of 496 | 8.7 (30.2 to 39.1) | concentrated in fewer lineages |
| ciprofloxacin | 8.7 % | 8.1 % (5.9 % to 10.4 %) | 62 of 531 | 1.9 (31.4 to 40.7) | concentrated in fewer lineages |
| nalidixic acid | 8.4 % | 2.2 % (1.1 % to 3.6 %) | 24 of 496 | 1.4 (30.0 to 38.7) | concentrated in fewer lineages |
| chloramphenicol | 6.0 % | 6.0 % (4.1 % to 7.9 %) | 61 of 534 | 2.0 (29.2 to 39.8) | concentrated in fewer lineages |
| trimethoprim-sulfamethoxazole | 3.0 % | 1.4 % (0.6 % to 2.4 %) | 17 of 496 | 1.3 (24.7 to 38.4) | concentrated in fewer lineages |
| azithromycin | 0.5 % | 1.1 % (0.3 % to 2.1 %) | 8 of 408 | 2.7 (8.2 to 20.6) | concentrated in fewer lineages |
| amikacin | 0.0 % | 0.0 % (0.0 % to 0.0 %) | 1 of 274 | 1.0 (1.0 to 1.0) | proportional |

*Figure 4 (drawn on the page).* Left: prevalence per isolate (filled) and per lineage (hollow, with its 95 % interval) for each trait, 18 largest of 19 by prevalence per isolate; a long connector means that a few large lineages carry most of the resistance, or that many small ones do. Right: the effective number of carrying lineages (point) against the permutation interval for carriage in proportion to lineage size (grey band); a point below the band is carriage concentrated in fewer lineages than chance gives, above it carriage spread over more, and a point inside it with a departure is carriage in as many lineages as chance gives but in uneven shares. Traits that do not depart from proportional carriage at level 0.05 are drawn in grey.

A prevalence difference between two collections can be split into a change in lineage composition and a change in rate within lineages. That decomposition needs two collections and was not run here: this record holds one.

## 7. Provenance and terms

- **Software:** amr-clonalshare 1.0.0, record schema 1.0
- **Seed:** 42
- **Configuration:** `sha256:72686f22ceea`
- **Record:** `sha256:6fe3209471d9bd8edf5ecaa60d4684dec4062942f9bc417a6024c01380742c1e`

> **Terms used in this report**
>
> **Share.** How much of the difference in resistance between isolates is explained by which lineage they belong to. Near 1: resistance travels with lineages. Near 0: it moves between them.
>
> **Interval.** The range the true value is expected to lie in. It is an interval for the share the lineages in this collection carry, not for a fresh draw of lineages; the fraction of the time such intervals contain that truth was measured on the release's validation grid.
>
> **Control.** The same calculation on shuffled lineage labels, which should return roughly zero; a share is read against it.
>
> **Support.** The share of isolates sitting in lineages large enough to inform the estimate; below 90 % the software withholds the value.
>
> **Gate.** A condition fixed before the run. If it does not hold, the software withholds the estimate and names the condition that failed.
>
> **e-value.** Evidence on a scale that stays honest when a panel is looked at again every year; larger is stronger, 1 is none.
>

## Notes on reading this report

- ‡ The 95 % interval includes zero. No lineage effect is distinguishable from none for this trait, and the point estimate must not be read on its own.
- † Support below 90 %. Too few isolates sit in lineages large enough to inform the estimate.
- § Derived arithmetically from a quantity recorded elsewhere in the record rather than estimated in this run.

amr-clonalshare 1.0.0 · record schema 1.0 · seed 42 · configuration sha256:72686f22ceea
Cite this run as: “amr-clonalshare 1.0.0, run ADDB37A7, record sha256:6fe3209471d9.”
This report supersedes any earlier report bearing the same run identifier. It is regenerated from the record and holds no value that the record does not. The symbols carry the same wording in every run of this software; no symbol against a value means only that none of the listed conditions fired. The interval for the lineages in hand is printed as computed, without clipping to the natural bounds of the quantity; clipping would leave its coverage unchanged and is left to the reader so that widths stay comparable across runs. The species interval is floored at zero, which is part of its construction.
