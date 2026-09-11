# amr-clonalshare run report

Every number below is read from `clonal_share_result.json`, written by the same run; none is recomputed here.

- **Run:** 074C4F05
- **Issued:** 2026-09-11T10:21Z
- **Software:** amr-clonalshare 1.0.0
- **Isolates:** 677
- **Lineage:** baps_cluster
- **Seed:** 42
- **Configuration:** sha256:32de5c33c6c1
- **Record digest:** sha256:2882fafd821a84134b43b7ff66dbc71c2376a01c46c672ac1767fd752ec37512

**13 of 13 antimicrobials estimable.** Every antimicrobial read cleared the support gate.

## 1. Measurement summary

**Quantity measured.** The share of the variation in resistance across this collection that lineage membership accounts for, read from a lineage label and a susceptibility result and scored on isolates the estimator did not see.

- **Cohort:** 677 isolates in 30 lineages, typed by `baps_cluster`
- **Antimicrobials read:** 13
- **Estimable:** 13 of 13

**Gates.** Every antimicrobial read cleared the support gate. A withheld value is a completed reading: it says the collection cannot identify the quantity at this typing resolution. A value that clears the gate is not thereby correct.

**Table 1.** Every trait, ordered by share. The reading in the last column is fixed by two conditions the record holds: whether the interval excludes zero and whether the e-BH procedure selected the trait at level 0.05.

| Trait | Share | 95 % interval | Control | e-value |  | Reading |
|---|---:|---:|---:|---:|---|---|
| penicillin | 0.517 | 0.237 to 0.698 | -0.06 | 5.6 × 10¹⁷ |  | evidence of a lineage effect |
| ceftiofur | 0.507 | 0.238 to 0.694 | -0.06 | 1.3 × 10¹⁷ |  | evidence of a lineage effect |
| tiamulin | 0.397 | 0.229 to 0.528 | -0.06 | 3.2 × 10¹³ |  | evidence of a lineage effect |
| trimethoprim | 0.301 | 0.169 to 0.379 | -0.06 | 1.8 × 10¹⁰ |  | evidence of a lineage effect |
| amoxicillin | 0.241 | -0.032 to 0.446 | -0.06 | 1.6 × 10⁵ | ‡ | lineage effect selected by e-BH; its size is not resolved |
| lincomycin | 0.235 | 0.111 to 0.367 | -0.06 | 6.3 × 10⁷ |  | evidence of a lineage effect |
| cefquinome | 0.185 | -0.041 to 0.419 | -0.06 | 1.2 × 10³ | ‡ | lineage effect selected by e-BH; its size is not resolved |
| erythromycin | 0.163 | 0.075 to 0.272 | -0.05 | 1.3 × 10⁵ |  | evidence of a lineage effect |
| tylosin | 0.159 | 0.066 to 0.269 | -0.06 | 4.2 × 10⁴ |  | evidence of a lineage effect |
| tilmicosin | 0.149 | 0.064 to 0.264 | -0.06 | 4.0 × 10⁴ |  | evidence of a lineage effect |
| spectinomycin | 0.130 | -0.058 to 0.336 | -0.06 | 1.5 × 10³ | ‡ | lineage effect selected by e-BH; its size is not resolved |
| tetracycline | 0.028 | -0.049 to 0.074 | -0.06 | 2.2 | ‡ | no detectable lineage effect |
| doxycycline | 0.022 | -0.071 to 0.075 | -0.06 | 3.1 | ‡ | no detectable lineage effect |

> **Not evaluated on this run, and why**
>
> Decomposition of a prevalence difference into lineage composition and within-lineage rate: needs two collections; this record holds one.
>

> **Interpretation (generated from Table 1 by fixed rules)**
>
> For penicillin and ceftiofur, resistance travels with lineage: the share is at or above one half, so it is carried by particular clones and will move when those clones move. Measures that act on the clone, such as movement control and clone-directed surveillance, act on these agents.
>
> For tiamulin, trimethoprim, lincomycin, erythromycin, tylosin and tilmicosin, a lineage effect is detected but the share is below one half: most of the variation in resistance sits within lineages, so the resistance moves largely independently of the clone, and selection pressure, dosing and choice of agent are the levers that reach most of it.
>
> For amoxicillin, cefquinome and spectinomycin, the e-value selects a lineage effect while the interval for its size still reaches zero: there is evidence that resistance is not spread evenly across the lineages, and this collection is too small, or too uneven across its lineages, to say how much of it the lineages carry. Read the e-value as the finding and the share as not yet resolved.
>
> For tetracycline and doxycycline, no lineage effect is distinguishable from none in this collection. Their resistance is not tied to any lineage the data can see; this run gives no basis for clone-directed action on them, and a change in their prevalence is more plausibly a change in selection than in which clones are present.
>
> These sentences follow from the table by fixed rules and carry no judgement beyond it.
>

## 2. Admissibility of the input

- **Lineage column:** `baps_cluster`
- **Susceptibility calls:** 13 antimicrobials on 677 isolates
- **Recorded dilutions:** 10832 rows read, 10832 joined; 677 of 677 isolates carry a recorded dilution
- **Resampling:** 5 folds, 20 repeats, 400 bootstrap draws, 200 permutations per antimicrobial

Antimicrobials with fewer than 20 isolates of the rarer outcome: 0 of 13. Below that count the estimate is still computed and its interval says how little it rests on.

Lineage groups in the input: 30, of which 3 hold a single isolate. Support 99.6 % against the 90.0 % the estimator needs: accepted.

*Figure 1 (drawn on the page).* Isolates per lineage, largest first. A lineage of one isolate, drawn in orange, cannot be predicted out of sample and counts against support; support is the share of isolates in the other lineages, 99.6 % here against the 90.0 % the estimator requires. The largest lineage holds 23.8 % of the isolates, which sets how much one lineage can weigh in the share.

**Table 2.** Conditions the estimator requires before any result is reported.

| Condition | Observed | Required | Verdict |
|---|---:|---:|---|
| Lineage support, lowest over the antimicrobials read | 99.6 % | ≥ 90.0 % | accepted |
| Lineage groups used | 30 | reported | accepted |

## 3. The clonal share, trait by trait

Read the intervals as frequencies. If cohorts like this one were drawn again and again, the interval printed for a trait would cover that trait's true share on about 95 draws in 100. For **5 of the 13 traits** shown the interval includes zero, so no lineage effect is distinguishable from none for that trait.

The control column is the same estimator run on shuffled lineage labels; it should sit near zero, and a share is read against it rather than against zero.

*Figure 2 (drawn on the page).* Clonal share by trait, point estimate with 95 % interval, 13 largest of 13. Traits whose interval crosses zero are drawn in grey and marked ‡. The thin line beneath each interval is the species interval of Table 3, for a fresh draw of lineages; it is floored at zero by construction. The intervals are drawn as computed. The quantity lies between 0 and 1, so the part of an interval below zero carries no information: reading each interval as its overlap with that range leaves the coverage unchanged, and a lower limit at or below zero means the same thing either way.

On the release's validation grid, the interval this estimator prints for a binary trait contained the truth in 0.970 of 28,000 runs over 70 simulated cohorts (by cohort, 0.825 to 1.000). That figure belongs to the release, not to this run; it is what the phrase "95 % interval" was measured to mean.

Two intervals answer two questions. The interval above is for the share the lineages in this collection carry. The second interval below is for the share a fresh draw of lineages from the species would show, stated on the scale of the realised share of Table 3b below; it adds the sampling of the lineages themselves, on 29 degrees of freedom, is widened to the envelope of the first taken on that scale (with a lower end no smaller than zero, since a species share is not negative), and is the one to quote when the figure is read as a property of the species rather than of this collection. With few lineages it is markedly wider; with many the two nearly coincide.

For 2 traits one lineage carries more than half of the between-lineage variation (cefquinome 69.2 %; spectinomycin 58.2 %). The species interval describes lineage effects drawn from one law, and a collection in which one lineage carries the resistance is not that: on the validation grid a carrier law of this kind took the species interval below its level with few lineages while the interval for the lineages in hand held. For them read the first interval as the statement about this collection and the species interval with that reservation.

The last column restates the share on the latent scale of a threshold model, the scale on which a mixed model with a binomial link reports an intraclass correlation. The two scales are one-to-one at a given prevalence and the map is monotone, so the interval keeps its coverage; the latent figure is larger because a call discards the part of the liability that does not cross the threshold. A share at or below zero has no latent counterpart and is printed as zero.

**Table 3.** The two intervals for every trait: for lineage membership in this collection and for the species, with the share on the latent scale.

| Trait | Share | 95 % interval, these lineages | 95 % interval, species | Latent scale (95 % interval) |
|---|---:|---:|---:|---:|
| penicillin | 0.517 | 0.237 to 0.698 | 0.256 to 0.719 | 0.755 (0.404 to 0.903) |
| ceftiofur | 0.507 | 0.238 to 0.694 | 0.257 to 0.715 | 0.742 (0.402 to 0.899) |
| tiamulin | 0.397 | 0.229 to 0.528 | 0.248 to 0.594 | 0.636 (0.407 to 0.776) |
| trimethoprim | 0.301 | 0.169 to 0.379 | 0.184 to 0.486 | 0.481 (0.285 to 0.585) |
| amoxicillin | 0.241 | -0.032 to 0.446 | 0.000 to 0.472 | 0.547 (0.000 to 0.776) |
| lincomycin | 0.235 | 0.111 to 0.367 | 0.122 to 0.399 | 0.369 (0.179 to 0.553) |
| cefquinome | 0.185 | -0.041 to 0.419 | 0.000 to 0.444 | 0.495 (0.000 to 0.773) |
| erythromycin | 0.163 | 0.075 to 0.272 | 0.083 to 0.299 | 0.254 (0.118 to 0.415) |
| tylosin | 0.159 | 0.066 to 0.269 | 0.073 to 0.295 | 0.248 (0.105 to 0.410) |
| tilmicosin | 0.149 | 0.064 to 0.264 | 0.070 to 0.285 | 0.232 (0.101 to 0.404) |
| spectinomycin | 0.130 | -0.058 to 0.336 | 0.000 to 0.360 | 0.288 (0.000 to 0.605) |
| tetracycline | 0.028 | -0.049 to 0.074 | 0.000 to 0.103 | 0.064 (0.000 to 0.159) |
| doxycycline | 0.022 | -0.071 to 0.075 | 0.000 to 0.099 | 0.050 (0.000 to 0.160) |

**The realised share of the same call.** The variance-component ratio of the lineages in hand, with an interval that is exact under a Gaussian within-lineage law. Its gate reads the excess kurtosis of the within-lineage residuals and withholds the interval above 0.99, which on a binary call closes wherever the prevalence is far from one half; the estimate is printed either way, the interval only where the gate opened, for 5 of 13 traits here.

**Table 3b.** The realised share per trait and the verdict of its kurtosis gate.

| Trait | Realised share | Exact 95 % interval | Residual excess kurtosis | Gate |
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
| penicillin | 5.6 × 10¹⁷ | 40.87 | yes |
| ceftiofur | 1.3 × 10¹⁷ | 39.42 | yes |
| tiamulin | 3.2 × 10¹³ | 31.10 | yes |
| trimethoprim | 1.8 × 10¹⁰ | 23.59 | yes |
| lincomycin | 6.3 × 10⁷ | 17.96 | yes |
| amoxicillin | 1.6 × 10⁵ | 11.96 | yes |
| erythromycin | 1.3 × 10⁵ | 11.78 | yes |
| tylosin | 4.2 × 10⁴ | 10.64 | yes |
| tilmicosin | 4.0 × 10⁴ | 10.61 | yes |
| spectinomycin | 1.5 × 10³ | 7.33 | yes |
| cefquinome | 1.2 × 10³ | 7.08 | yes |
| doxycycline | 3.1 | 1.14 | no |
| tetracycline | 2.2 | 0.79 | no |

11 of 13 traits are selected. Note: e-value in the betting sense of Vovk and Wang, not the BLAST expectation value and not the E-value of VanderWeele and Ding.

*Figure 3 (drawn on the page).* Evidence per trait on the natural-log scale, 13 largest of 13. The dashed rule is 1/α, the evidence one trait alone needs at level 0.05; the solid rule is the e-BH selection threshold on this run, 23.6, which rises with the number of traits read together, and a trait at or beyond it is selected. Traits the e-BH procedure selected are drawn in blue; the scale is logarithmic, so equal steps are equal factors of evidence.

This cohort also carries intakes, 25 of them read in the order of collection_year, from 1 to 177 isolates each, with 26 isolates carrying no intake and set aside. Each intake is scored against the lineage rates learned from the intakes before it, and the running product is the sequential e-value: it may be read after any intake without spending the guarantee. 2 of 13 traits are selected on the product, against 11 at this single look, which is the price of a guarantee that holds at whatever intake the programme is read at.

**Table 4b.** The sequential e-value per trait after the last intake, and the e-BH selection taken on the log scale, where a product over intakes does not overflow.

| Trait | natural log of the product | Selected |
|---|---:|---|
| ceftiofur | 61.56 | yes |
| penicillin | 53.99 | yes |
| cefquinome | -10.02 | no |
| amoxicillin | -15.70 | no |
| tiamulin | -22.87 | no |
| trimethoprim | -49.18 | no |
| lincomycin | -107.30 | no |
| doxycycline | -110.84 | no |
| tetracycline | -113.99 | no |
| tilmicosin | -124.03 | no |
| tylosin | -128.73 | no |
| erythromycin | -132.72 | no |
| spectinomycin | -146.88 | no |

*Figure 4 (drawn on the page).* The running product of e-values over the 25 intakes, one line per trait, on the natural-log scale. Each intake is scored against the lineage rates learned from the intakes before it, so the first intakes carry no evidence and a line may fall as well as rise. The dashed rule is 1/α; the solid rule is the e-BH threshold on the product, 130.0. Traits selected on the product are drawn in blue and named at the right; the guarantee holds at whatever intake the line is read. The frame is cut at −25: a line that leaves it at the bottom has accumulated evidence against a lineage effect and does not return in the intakes shown.

## 5. Reading at the recorded resolution

Where a dilution was recorded, the share is also read from the dilution itself rather than from the call derived from it. A call and a dilution are both intervals on one concentration scale, so the dilution carries more of the reading; a reading at an end well is censored and enters as an interval. This table reads the share of that scale carried by lineage, beside the share of the binary call above, which is a different quantity on the same isolates.

**Table 5.** Share of the dilution scale carried by lineage, per agent, from the interval-censored likelihood. The realised interval is for the lineages this cohort holds; the 95 % interval is for the species behind it.

| Agent | Share | 95 % interval | Realised interval | Censored readings | Estimable |
|---|---:|---:|---:|---:|---|
| amoxicillin | 0.713 | 0.605 to 0.821 | 0.678 to 0.743 | 7.5 % | yes |
| cefquinome | 0.468 | 0.347 to 0.622 | 0.414 to 0.518 | 1.5 % | yes |
| ceftiofur | 0.673 | 0.559 to 0.791 | 0.635 to 0.706 | 1.3 % | yes |
| doxycycline | 0.074 | 0.031 to 0.156 | 0.035 to 0.124 | 0.9 % | yes |
| enrofloxacin | 0.221 | 0.137 to 0.357 | 0.164 to 0.280 | 0.4 % | yes |
| erythromycin | 0.228 | 0.143 to 0.367 | 0.171 to 0.288 | 56.4 % | yes |
| florfenicol | 0.240 | 0.151 to 0.382 | 0.181 to 0.300 | 3.2 % | yes |
| lincomycin | 0.357 | 0.247 to 0.513 | 0.297 to 0.414 | 55.8 % | yes |
| marbofloxacin | 0.037 | 0.006 to 0.098 | 0.006 to 0.081 | 4.3 % | yes |
| penicillin | 0.679 | 0.562 to 0.797 | 0.632 to 0.718 | 72.7 % | yes |
| spectinomycin | 0.463 | 0.342 to 0.617 | 0.409 to 0.513 | 7.8 % | yes |
| tetracycline | 0.078 | 0.034 to 0.161 | 0.037 to 0.128 | 9.9 % | yes |
| tiamulin | 0.695 | 0.584 to 0.808 | 0.660 to 0.726 | 8.0 % | yes |
| tilmicosin | 0.172 | 0.101 to 0.295 | 0.119 to 0.230 | 1.9 % | yes |
| trimethoprim | 0.441 | 0.322 to 0.596 | 0.386 to 0.492 | 9.9 % | yes |
| tylosin | 0.270 | 0.175 to 0.416 | 0.211 to 0.329 | 52.4 % | yes |

*Figure 5 (drawn on the page).* Share of the dilution scale carried by lineage, per agent, 16 largest of 16: point estimate with the 95 % interval for the species (thin) and the realised interval for these lineages (thick); an agent whose reading is refused is drawn in grey. The hollow diamond is the share of the binary call from Table 1 for the same agent, a different quantity on the same isolates, placed here so that the two readings can be seen side by side. The figure in parentheses at the right is the share of readings at an end well, which enter as censored.

## 6. How resistance is carried across lineages

Of the **13** traits with a carriage reading, **11** depart from carriage in proportion to lineage size far enough to matter for a prevalence reading: 3 in fewer lineages than chance gives and 8 in as many carrying lineages as chance gives but not in the same shares across them. The remaining 2 are consistent with proportional carriage.

The widest gap between the per-isolate and the per-lineage prevalence is 0.178, on `ceftiofur`. Where the two differ, a small number of large lineages carry most of the resistance, and a change in prevalence may be a change in which lineages were sampled rather than a change in rate.

Two prevalences are reported for every trait, and they are two quantities rather than a biased and an unbiased one: per isolate is the clinical burden of the collection, per lineage is the share of its diversity that carries the trait. They separate whenever sampling across lineages is uneven. Beside them the effective number of carrying lineages, the reciprocal of a Herfindahl index over the carriers, is compared with the number that carriage in proportion to lineage size would give, by permutation. A trait departs from proportional carriage where that permutation p-value falls below 0.05; the direction says whether the carriers sit in fewer lineages than chance gives, in more, or in as many but in uneven shares across them.

**Table 6.** Prevalence on two scales and the concentration of carriage, per trait. The per-lineage interval is a two-stage cluster bootstrap over lineages and then isolates; the null interval for the effective number is the permutation reference for carriage in proportion to lineage size.

| Trait | Per isolate | Per lineage (95 % interval) | Carrying lineages | Effective number (null 95 %) | Carriage |
|---|---:|---:|---:|---:|---|
| tetracycline | 84.8 % | 87.6 % (80.6 % to 93.1 %) | 30 of 30 | 10.5 (9.7 to 10.8) | proportional |
| doxycycline | 84.3 % | 85.4 % (77.6 % to 91.9 %) | 30 of 30 | 10.4 (9.7 to 10.8) | proportional |
| lincomycin | 62.6 % | 73.0 % (61.4 % to 83.1 %) | 29 of 30 | 10.7 (9.3 to 11.2) | as many lineages as chance, uneven shares |
| erythromycin | 54.1 % | 60.3 % (47.9 % to 72.8 %) | 28 of 30 | 9.6 (9.1 to 11.3) | as many lineages as chance, uneven shares |
| tylosin | 54.1 % | 59.2 % (47.7 % to 71.1 %) | 28 of 30 | 9.4 (9.0 to 11.4) | as many lineages as chance, uneven shares |
| tilmicosin | 53.8 % | 58.7 % (46.2 % to 71.0 %) | 28 of 30 | 9.6 (9.1 to 11.4) | as many lineages as chance, uneven shares |
| trimethoprim | 28.1 % | 38.1 % (24.7 % to 51.5 %) | 23 of 30 | 9.7 (8.1 to 12.0) | as many lineages as chance, uneven shares |
| ceftiofur | 24.1 % | 41.9 % (27.1 % to 56.7 %) | 22 of 30 | 7.4 (8.0 to 12.2) | concentrated in fewer lineages |
| penicillin | 23.0 % | 40.3 % (25.6 % to 56.3 %) | 23 of 30 | 7.1 (7.8 to 12.2) | concentrated in fewer lineages |
| tiamulin | 19.4 % | 30.0 % (17.4 % to 43.7 %) | 20 of 30 | 7.1 (7.5 to 12.4) | concentrated in fewer lineages |
| spectinomycin | 11.5 % | 19.6 % (9.8 % to 31.0 %) | 22 of 30 | 10.7 (6.8 to 12.7) | as many lineages as chance, uneven shares |
| amoxicillin | 5.5 % | 11.2 % (3.5 % to 20.1 %) | 13 of 30 | 6.3 (5.5 to 12.3) | as many lineages as chance, uneven shares |
| cefquinome | 3.8 % | 5.8 % (1.2 % to 12.8 %) | 9 of 30 | 6.8 (4.5 to 11.7) | as many lineages as chance, uneven shares |

*Figure 6 (drawn on the page).* Left: prevalence per isolate (filled) and per lineage (hollow, with its 95 % interval) for each trait, 13 largest of 13 by prevalence per isolate; a long connector means that a few large lineages carry most of the resistance, or that many small ones do. Right: the effective number of carrying lineages (point) against the permutation interval for carriage in proportion to lineage size (grey band); a point below the band is carriage concentrated in fewer lineages than chance gives, above it carriage spread over more, and a point inside it with a departure is carriage in as many lineages as chance gives but in uneven shares. Traits that do not depart from proportional carriage at level 0.05 are drawn in grey.

A prevalence difference between two collections can be split into a change in lineage composition and a change in rate within lineages. That decomposition needs two collections and was not run here: this record holds one.

## 7. Provenance and terms

- **Software:** amr-clonalshare 1.0.0, record schema 1.0
- **Seed:** 42
- **Configuration:** `sha256:32de5c33c6c1`
- **Record:** `sha256:2882fafd821a84134b43b7ff66dbc71c2376a01c46c672ac1767fd752ec37512`

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

amr-clonalshare 1.0.0 · record schema 1.0 · seed 42 · configuration sha256:32de5c33c6c1
Cite this run as: “amr-clonalshare 1.0.0, run 074C4F05, record sha256:2882fafd821a.”
This report supersedes any earlier report bearing the same run identifier. It is regenerated from the record and holds no value that the record does not. The symbols carry the same wording in every run of this software; no symbol against a value means only that none of the listed conditions fired. The interval for the lineages in hand is printed as computed, without clipping to the natural bounds of the quantity; clipping would leave its coverage unchanged and is left to the reader so that widths stay comparable across runs. The species interval is floored at zero, which is part of its construction.
