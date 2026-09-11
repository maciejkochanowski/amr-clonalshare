# amr-clonalshare run report

Every number below is read from `clonal_share_result.json`, written by the same run; none is recomputed here.

- **Run:** ACCE2632
- **Issued:** 2026-09-11T10:22Z
- **Software:** amr-clonalshare 1.0.0
- **Isolates:** 7049
- **Lineage:** serovar
- **Seed:** 42
- **Configuration:** sha256:046fc826d97d
- **Record digest:** sha256:0356b07181b3d760a3d1139e3da2f445b7adb34f0568f1dfbb9790d077faaed3

**18 of 22 antimicrobials estimable.** 4 of 22 could not be scored, for a reason named below.

## 1. Measurement summary

**Quantity measured.** The share of the variation in resistance across this collection that lineage membership accounts for, read from a lineage label and a susceptibility result and scored on isolates the estimator did not see.

- **Cohort:** 7049 isolates in 85 lineages, typed by `serovar`
- **Antimicrobials read:** 22
- **Estimable:** 18 of 22

**Gates.** 4 of 22 could not be scored, for a reason named below. A withheld value is a completed reading: it says the collection cannot identify the quantity at this typing resolution. A value that clears the gate is not thereby correct.

**Table 1.** Every trait, ordered by share. The reading in the last column is fixed by two conditions the record holds: whether the interval excludes zero and whether the e-BH procedure selected the trait at level 0.05.

| Trait | Share | 95 % interval | Control | e-value |  | Reading |
|---|---:|---:|---:|---:|---|---|
| nalidixic acid | 0.765 | 0.017 to 0.851 | -0.01 | 1.4 × 10¹⁴⁴ |  | evidence of a lineage effect |
| ciprofloxacin | 0.666 | 0.065 to 0.792 | -0.01 | 2.8 × 10¹²⁷ |  | evidence of a lineage effect |
| sulfisoxazole | 0.478 | 0.116 to 0.635 | -0.01 | 3.6 × 10¹⁶² |  | evidence of a lineage effect |
| streptomycin | 0.365 | 0.120 to 0.476 | -0.01 | 1.2 × 10¹³⁰ |  | evidence of a lineage effect |
| chloramphenicol | 0.346 | 0.000 to 0.403 | -0.01 | 3.4 × 10⁷⁵ |  | evidence of a lineage effect |
| tetracycline | 0.335 | 0.176 to 0.532 | -0.01 | 5.3 × 10¹³⁷ |  | evidence of a lineage effect |
| trimethoprim-sulfamethoxazole | 0.276 | -0.003 to 0.300 | -0.01 | 9.8 × 10⁴⁶ | ‡ | lineage effect selected by e-BH; its size is not resolved |
| ceftazidime | 0.274 | -0.153 to 0.378 | -0.16 | 6.0 | ‡ | no detectable lineage effect |
| ceftriaxone | 0.234 | 0.052 to 0.331 | -0.01 | 3.5 × 10⁶⁸ |  | evidence of a lineage effect |
| ceftiofur | 0.224 | 0.058 to 0.347 | -0.02 | 4.7 × 10³⁹ |  | evidence of a lineage effect |
| cefoxitin | 0.209 | 0.029 to 0.356 | -0.01 | 4.5 × 10⁵⁵ |  | evidence of a lineage effect |
| gentamicin | 0.202 | 0.091 to 0.282 | -0.01 | 9.5 × 10⁶⁵ |  | evidence of a lineage effect |
| sulfamethoxazole | 0.198 | -0.118 to 0.375 | -0.10 | 221.3 | ‡ | lineage effect selected by e-BH; its size is not resolved |
| amoxicillin-clavulanic acid | 0.194 | 0.067 to 0.284 | -0.01 | 5.2 × 10⁶¹ |  | evidence of a lineage effect |
| ampicillin | 0.193 | 0.104 to 0.252 | -0.01 | 1.3 × 10⁶⁹ |  | evidence of a lineage effect |
| kanamycin | 0.117 | 0.068 to 0.180 | -0.02 | 1.1 × 10¹⁷ |  | evidence of a lineage effect |
| amikacin | 0.024 | 0.022 to 0.024 | -0.03 | 0.6 |  | interval excludes zero; not selected by e-BH |
| azithromycin | 0.021 | -0.024 to 0.031 | -0.02 | 351.0 | ‡ | lineage effect selected by e-BH; its size is not resolved |

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
> For nalidixic acid and ciprofloxacin, resistance travels with lineage: the share is at or above one half, so it is carried by particular clones and will move when those clones move. Measures that act on the clone, such as movement control and clone-directed surveillance, act on these agents.
>
> For sulfisoxazole, streptomycin, chloramphenicol, tetracycline, ceftriaxone, ceftiofur, cefoxitin, gentamicin, amoxicillin-clavulanic acid, ampicillin and kanamycin, a lineage effect is detected but the share is below one half: most of the variation in resistance sits within lineages, so the resistance moves largely independently of the clone, and selection pressure, dosing and choice of agent are the levers that reach most of it.
>
> For trimethoprim-sulfamethoxazole, sulfamethoxazole and azithromycin, the e-value selects a lineage effect while the interval for its size still reaches zero: there is evidence that resistance is not spread evenly across the lineages, and this collection is too small, or too uneven across its lineages, to say how much of it the lineages carry. Read the e-value as the finding and the share as not yet resolved.
>
> For amikacin, the interval for the share excludes zero while the e-value does not select a lineage effect at the panel's false-discovery level: the two readings disagree, and the e-value, which carries the multiplicity correction, is the one to read; treat the share as suggestive.
>
> For ceftazidime, no lineage effect is distinguishable from none in this collection. Its resistance is not tied to any lineage the data can see; this run gives no basis for clone-directed action on it, and a change in its prevalence is more plausibly a change in selection than in which clones are present.
>
> These sentences follow from the table by fixed rules and carry no judgement beyond it.
>

## 2. Admissibility of the input

- **Lineage column:** `serovar`
- **Susceptibility calls:** 22 antimicrobials on 7049 isolates
- **Resampling:** 5 folds, 10 repeats, 300 bootstrap draws, 150 permutations per antimicrobial

Antimicrobials with fewer than 20 isolates of the rarer outcome: 1 of 22. Below that count the estimate is still computed and its interval says how little it rests on.

Lineage groups in the input: 85, of which 20 hold a single isolate. Support 99.7 % against the 90.0 % the estimator needs: accepted.

*Figure 1 (drawn on the page).* Isolates per lineage, largest first, the 40 largest of 85. A lineage of one isolate, drawn in orange, cannot be predicted out of sample and counts against support; support is the share of isolates in the other lineages, 99.7 % here against the 90.0 % the estimator requires. The largest lineage holds 19.9 % of the isolates, which sets how much one lineage can weigh in the share.

**Table 2.** Conditions the estimator requires before any result is reported.

| Condition | Observed | Required | Verdict |
|---|---:|---:|---|
| Lineage support, lowest over the antimicrobials read | 96.2 % | ≥ 90.0 % | accepted |
| Lineage groups used | 85 | reported | accepted |

## 3. The clonal share, trait by trait

Read the intervals as frequencies. If cohorts like this one were drawn again and again, the interval printed for a trait would cover that trait's true share on about 95 draws in 100. For **4 of the 18 traits** shown the interval includes zero, so no lineage effect is distinguishable from none for that trait.

The control column is the same estimator run on shuffled lineage labels; it should sit near zero, and a share is read against it rather than against zero.

*Figure 2 (drawn on the page).* Clonal share by trait, point estimate with 95 % interval, 18 largest of 18. Traits whose interval crosses zero are drawn in grey and marked ‡. The thin line beneath each interval is the species interval of Table 3, for a fresh draw of lineages; it is floored at zero by construction. The intervals are drawn as computed. The quantity lies between 0 and 1, so the part of an interval below zero carries no information: reading each interval as its overlap with that range leaves the coverage unchanged, and a lower limit at or below zero means the same thing either way.

On the release's validation grid, the interval this estimator prints for a binary trait contained the truth in 0.970 of 28,000 runs over 70 simulated cohorts (by cohort, 0.825 to 1.000). That figure belongs to the release, not to this run; it is what the phrase "95 % interval" was measured to mean.

Two intervals answer two questions. The interval above is for the share the lineages in this collection carry. The second interval below is for the share a fresh draw of lineages from the species would show, stated on the scale of the realised share of Table 3b below; it adds the sampling of the lineages themselves, on 84 degrees of freedom, is widened to the envelope of the first taken on that scale (with a lower end no smaller than zero, since a species share is not negative), and is the one to quote when the figure is read as a property of the species rather than of this collection. With few lineages it is markedly wider; with many the two nearly coincide.

For 7 traits one lineage carries more than half of the between-lineage variation (nalidixic acid 87.2 %; ciprofloxacin 85.0 %; chloramphenicol 85.3 %; trimethoprim-sulfamethoxazole 86.9 %; ceftiofur 61.7 %; cefoxitin 71.6 %; and others). The species interval describes lineage effects drawn from one law, and a collection in which one lineage carries the resistance is not that: on the validation grid a carrier law of this kind took the species interval below its level with few lineages while the interval for the lineages in hand held. For them read the first interval as the statement about this collection and the species interval with that reservation.

The last column restates the share on the latent scale of a threshold model, the scale on which a mixed model with a binomial link reports an intraclass correlation. The two scales are one-to-one at a given prevalence and the map is monotone, so the interval keeps its coverage; the latent figure is larger because a call discards the part of the liability that does not cross the threshold. A share at or below zero has no latent counterpart and is printed as zero.

**Table 3.** The two intervals for every trait: for lineage membership in this collection and for the species, with the share on the latent scale.

| Trait | Share | 95 % interval, these lineages | 95 % interval, species | Latent scale (95 % interval) |
|---|---:|---:|---:|---:|
| nalidixic acid | 0.765 | 0.017 to 0.851 | 0.019 to 0.862 | 0.957 (0.053 to 0.983) |
| ciprofloxacin | 0.666 | 0.065 to 0.792 | 0.070 to 0.806 | 0.911 (0.175 to 0.965) |
| sulfisoxazole | 0.478 | 0.116 to 0.635 | 0.126 to 0.656 | 0.701 (0.198 to 0.851) |
| streptomycin | 0.365 | 0.120 to 0.476 | 0.129 to 0.496 | 0.547 (0.191 to 0.683) |
| chloramphenicol | 0.346 | 0.000 to 0.403 | 0.001 to 0.451 | 0.671 (0.002 to 0.730) |
| tetracycline | 0.335 | 0.176 to 0.532 | 0.189 to 0.554 | 0.503 (0.273 to 0.742) |
| trimethoprim-sulfamethoxazole | 0.276 | -0.003 to 0.300 | 0.000 to 0.382 | 0.644 (0.000 to 0.672) |
| ceftazidime | 0.274 | -0.153 to 0.378 | 0.000 to 0.552 | 0.444 (0.000 to 0.584) |
| ceftriaxone | 0.234 | 0.052 to 0.331 | 0.056 to 0.351 | 0.453 (0.121 to 0.588) |
| ceftiofur | 0.224 | 0.058 to 0.347 | 0.063 to 0.368 | 0.425 (0.128 to 0.597) |
| cefoxitin | 0.209 | 0.029 to 0.356 | 0.032 to 0.377 | 0.446 (0.082 to 0.646) |
| gentamicin | 0.202 | 0.091 to 0.282 | 0.099 to 0.301 | 0.393 (0.197 to 0.513) |
| sulfamethoxazole | 0.198 | -0.118 to 0.375 | 0.000 to 0.429 | 0.351 (0.000 to 0.601) |
| amoxicillin-clavulanic acid | 0.194 | 0.067 to 0.284 | 0.073 to 0.303 | 0.364 (0.140 to 0.501) |
| ampicillin | 0.193 | 0.104 to 0.252 | 0.113 to 0.269 | 0.326 (0.182 to 0.415) |
| kanamycin | 0.117 | 0.068 to 0.180 | 0.074 to 0.194 | 0.253 (0.154 to 0.363) |
| amikacin | 0.024 | 0.022 to 0.024 | 0.015 to 0.042 | 0.427 (0.420 to 0.431) |
| azithromycin | 0.021 | -0.024 to 0.031 | 0.000 to 0.048 | 0.236 (0.000 to 0.294) |

**The realised share of the same call.** The variance-component ratio of the lineages in hand, with an interval that is exact under a Gaussian within-lineage law. Its gate reads the excess kurtosis of the within-lineage residuals and withholds the interval above 0.99, which on a binary call closes wherever the prevalence is far from one half; the estimate is printed either way, the interval only where the gate opened, for 5 of 22 traits here.

**Table 3b.** The realised share per trait and the verdict of its kurtosis gate.

| Trait | Realised share | Exact 95 % interval | Residual excess kurtosis | Gate |
|---|---:|---:|---:|---|
| nalidixic acid | 0.780 | withheld | 41.20 | closed |
| ciprofloxacin | 0.684 | withheld | 25.25 | closed |
| sulfisoxazole | 0.500 | withheld | 1.93 | closed |
| streptomycin | 0.381 | 0.364 to 0.398 | -0.09 | open |
| chloramphenicol | 0.368 | withheld | 10.80 | closed |
| tetracycline | 0.352 | 0.335 to 0.369 | -0.72 | open |
| trimethoprim-sulfamethoxazole | 0.286 | withheld | 17.70 | closed |
| ceftriaxone | 0.242 | withheld | 2.83 | closed |
| ceftiofur | 0.236 | withheld | 2.10 | closed |
| sulfamethoxazole | 0.223 | 0.125 to 0.326 | 0.42 | open |
| cefoxitin | 0.220 | withheld | 5.78 | closed |
| gentamicin | 0.212 | withheld | 1.96 | closed |
| ampicillin | 0.201 | 0.184 to 0.218 | -0.39 | open |
| amoxicillin-clavulanic acid | 0.197 | withheld | 1.24 | closed |
| ceftazidime | 0.183 | 0.037 to 0.342 | -0.52 | open |
| kanamycin | 0.105 | withheld | 2.17 | closed |
| azithromycin | 0.013 | withheld | 183.50 | closed |
| amikacin | 0.000 | withheld | 2376.70 | closed |
| colistin | not computed | withheld | not computed | closed |
| imipenem | not computed | withheld | not computed | closed |
| meropenem | not computed | withheld | not computed | closed |
| spectinomycin | not computed | withheld | not computed | closed |

## 4. Evidence that survives re-reading

A p-value is a statement about one look at the data. A surveillance panel is looked at again every year, and a p-value recomputed each time loses its guarantee. The e-value is evidence on a scale made for that: this run's e-value is a statement about this cohort, and a programme that adds an intake each year multiplies the e-value of each new intake, scored against the lineage rates learned from the earlier ones, into a running product whose guarantee holds at whatever intake it is read (sequential_e_process). The e-BH procedure controls the false-discovery rate across the panel whatever the dependence between traits. Larger is stronger; 1 is no evidence.

**Table 4.** e-value per trait and the e-BH selection at level 0.05. The selection threshold on this run is 27.5; a trait at or above it is selected.

| Trait | e-value | natural log | Selected |
|---|---:|---:|---|
| sulfisoxazole | 3.6 × 10¹⁶² | 374.30 | yes |
| nalidixic acid | 1.4 × 10¹⁴⁴ | 331.92 | yes |
| tetracycline | 5.3 × 10¹³⁷ | 317.12 | yes |
| streptomycin | 1.2 × 10¹³⁰ | 299.53 | yes |
| ciprofloxacin | 2.8 × 10¹²⁷ | 293.47 | yes |
| chloramphenicol | 3.4 × 10⁷⁵ | 173.92 | yes |
| ampicillin | 1.3 × 10⁶⁹ | 159.13 | yes |
| ceftriaxone | 3.5 × 10⁶⁸ | 157.83 | yes |
| gentamicin | 9.5 × 10⁶⁵ | 151.91 | yes |
| amoxicillin-clavulanic acid | 5.2 × 10⁶¹ | 142.10 | yes |
| cefoxitin | 4.5 × 10⁵⁵ | 128.14 | yes |
| trimethoprim-sulfamethoxazole | 9.8 × 10⁴⁶ | 108.20 | yes |
| ceftiofur | 4.7 × 10³⁹ | 91.35 | yes |
| kanamycin | 1.1 × 10¹⁷ | 39.20 | yes |
| azithromycin | 351.0 | 5.86 | yes |
| sulfamethoxazole | 221.3 | 5.40 | yes |
| ceftazidime | 6.0 | 1.80 | no |
| colistin | 1.0 | 0.00 | no |
| imipenem | 1.0 | 0.00 | no |
| meropenem | 1.0 | 0.00 | no |
| amikacin | 0.6 | -0.47 | no |
| spectinomycin | not computed | not computed | no |

16 of 22 traits are selected. Note: e-value in the betting sense of Vovk and Wang, not the BLAST expectation value and not the E-value of VanderWeele and Ding.

*Figure 3 (drawn on the page).* Evidence per trait on the natural-log scale, 18 largest of 21. The dashed rule is 1/α, the evidence one trait alone needs at level 0.05; the solid rule is the e-BH selection threshold on this run, 27.5, which rises with the number of traits read together, and a trait at or beyond it is selected. Traits the e-BH procedure selected are drawn in blue; the scale is logarithmic, so equal steps are equal factors of evidence.

## 5. Reading at the recorded resolution

No recorded dilutions were supplied on this run; the shares above are read from binary calls only.

## 6. How resistance is carried across lineages

Of the **19** traits with a carriage reading, **17** depart from carriage in proportion to lineage size far enough to matter for a prevalence reading: 17 in fewer lineages than chance gives. The remaining 2 are consistent with proportional carriage.

The widest gap between the per-isolate and the per-lineage prevalence is -0.194, on `tetracycline`. Where the two differ, a small number of large lineages carry most of the resistance, and a change in prevalence may be a change in which lineages were sampled rather than a change in rate.

Two prevalences are reported for every trait, and they are two quantities rather than a biased and an unbiased one: per isolate is the clinical burden of the collection, per lineage is the share of its diversity that carries the trait. They separate whenever sampling across lineages is uneven. Beside them the effective number of carrying lineages, the reciprocal of a Herfindahl index over the carriers, is compared with the number that carriage in proportion to lineage size would give, by permutation. A trait departs from proportional carriage where that permutation p-value falls below 0.05; the direction says whether the carriers sit in fewer lineages than chance gives, in more, or in as many but in uneven shares across them.

**Table 6.** Prevalence on two scales and the concentration of carriage, per trait. The per-lineage interval is a two-stage cluster bootstrap over lineages and then isolates; the null interval for the effective number is the permutation reference for carriage in proportion to lineage size.

| Trait | Per isolate | Per lineage (95 % interval) | Carrying lineages | Effective number (null 95 %) | Carriage |
|---|---:|---:|---:|---:|---|
| colistin | 100.0 % | 100.0 % (100.0 % to 100.0 %) | 29 of 29 | 3.5 (3.5 to 3.5) | proportional |
| tetracycline | 52.6 % | 33.3 % (25.7 % to 41.4 %) | 52 of 85 | 8.5 (11.4 to 12.2) | concentrated in fewer lineages |
| streptomycin | 40.8 % | 24.7 % (17.8 % to 32.1 %) | 50 of 82 | 5.6 (12.2 to 13.2) | concentrated in fewer lineages |
| sulfisoxazole | 28.9 % | 19.8 % (13.5 % to 26.7 %) | 44 of 84 | 6.9 (10.9 to 12.2) | concentrated in fewer lineages |
| ceftazidime | 28.0 % | 11.9 % (2.4 % to 24.3 %) | 7 of 20 | 1.7 (3.1 to 6.9) | concentrated in fewer lineages |
| ampicillin | 26.4 % | 16.2 % (11.1 % to 21.8 %) | 44 of 85 | 9.9 (11.1 to 12.5) | concentrated in fewer lineages |
| sulfamethoxazole | 21.3 % | 16.2 % (5.5 % to 28.7 %) | 10 of 23 | 3.8 (4.9 to 8.8) | concentrated in fewer lineages |
| amoxicillin-clavulanic acid | 16.8 % | 8.6 % (4.9 % to 13.1 %) | 33 of 84 | 7.9 (10.7 to 12.5) | concentrated in fewer lineages |
| ceftiofur | 14.6 % | 9.9 % (5.4 % to 15.4 %) | 28 of 59 | 4.3 (9.9 to 12.0) | concentrated in fewer lineages |
| gentamicin | 14.3 % | 10.1 % (6.4 % to 14.3 %) | 38 of 85 | 6.2 (10.8 to 12.8) | concentrated in fewer lineages |
| kanamycin | 13.3 % | 9.9 % (4.6 % to 16.6 %) | 25 of 56 | 5.2 (9.5 to 11.7) | concentrated in fewer lineages |
| ceftriaxone | 13.0 % | 7.5 % (3.9 % to 11.8 %) | 32 of 84 | 5.3 (10.5 to 12.7) | concentrated in fewer lineages |
| cefoxitin | 9.4 % | 7.7 % (3.7 % to 12.5 %) | 32 of 84 | 4.6 (10.3 to 12.8) | concentrated in fewer lineages |
| ciprofloxacin | 8.7 % | 6.1 % (2.2 % to 11.0 %) | 27 of 85 | 1.9 (10.4 to 13.0) | concentrated in fewer lineages |
| nalidixic acid | 8.4 % | 4.1 % (0.6 % to 8.4 %) | 15 of 84 | 1.4 (10.2 to 12.9) | concentrated in fewer lineages |
| chloramphenicol | 6.0 % | 7.4 % (3.1 % to 12.8 %) | 31 of 85 | 1.9 (10.1 to 13.2) | concentrated in fewer lineages |
| trimethoprim-sulfamethoxazole | 3.0 % | 2.2 % (0.3 % to 5.1 %) | 13 of 84 | 1.3 (9.2 to 13.5) | concentrated in fewer lineages |
| azithromycin | 0.5 % | 0.3 % (0.1 % to 0.7 %) | 7 of 73 | 2.6 (4.4 to 12.0) | concentrated in fewer lineages |
| amikacin | 0.0 % | 0.0 % (0.0 % to 0.0 %) | 1 of 51 | 1.0 (1.0 to 1.0) | proportional |

*Figure 4 (drawn on the page).* Left: prevalence per isolate (filled) and per lineage (hollow, with its 95 % interval) for each trait, 18 largest of 19 by prevalence per isolate; a long connector means that a few large lineages carry most of the resistance, or that many small ones do. Right: the effective number of carrying lineages (point) against the permutation interval for carriage in proportion to lineage size (grey band); a point below the band is carriage concentrated in fewer lineages than chance gives, above it carriage spread over more, and a point inside it with a departure is carriage in as many lineages as chance gives but in uneven shares. Traits that do not depart from proportional carriage at level 0.05 are drawn in grey.

A prevalence difference between two collections can be split into a change in lineage composition and a change in rate within lineages. That decomposition needs two collections and was not run here: this record holds one.

## 7. Provenance and terms

- **Software:** amr-clonalshare 1.0.0, record schema 1.0
- **Seed:** 42
- **Configuration:** `sha256:046fc826d97d`
- **Record:** `sha256:0356b07181b3d760a3d1139e3da2f445b7adb34f0568f1dfbb9790d077faaed3`

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

amr-clonalshare 1.0.0 · record schema 1.0 · seed 42 · configuration sha256:046fc826d97d
Cite this run as: “amr-clonalshare 1.0.0, run ACCE2632, record sha256:0356b07181b3.”
This report supersedes any earlier report bearing the same run identifier. It is regenerated from the record and holds no value that the record does not. The symbols carry the same wording in every run of this software; no symbol against a value means only that none of the listed conditions fired. The interval for the lineages in hand is printed as computed, without clipping to the natural bounds of the quantity; clipping would leave its coverage unchanged and is left to the reader so that widths stay comparable across runs. The species interval is floored at zero, which is part of its construction.
