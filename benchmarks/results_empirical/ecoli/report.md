# amr-clonalshare run report

Every number below is read from `clonal_share_result.json`, written by the same run; none is recomputed here.

- **Run:** 424A3AD1
- **Issued:** 2026-10-03T06:11Z
- **Software:** amr-clonalshare 1.0.0
- **Isolates:** 481
- **Lineage:** st
- **Seed:** 20261001
- **Configuration:** sha256:3ae0b5a290d7
- **Record digest:** sha256:2d32bfbb18333bf880cd99ad82b5a1276de3219db9fdc80e51443dff2ea797da

**6 antimicrobials; 18 analysis outcomes.** 17 computed, 1 full ranges, 0 unavailable. Each result retains its own target and data subset.

## 1. Measurement summary

**Quantity measured.** The share of the ordering of the recorded MIC readings within each stratum that lineage membership accounts for, and the lower bound on the lineage share of the latent MIC ordering established by the readings.

- **Collection:** 481 isolates in 75 lineages, typed by `st`
- **Antimicrobials read:** 6

**Reporting conditions.** 17 computed, 1 full ranges, 0 unavailable. Each result retains its own target and data subset. An unavailable value failed a declared reporting condition of its estimator; the recorded reason distinguishes input limitations from numerical failure.

**Table 0.** Available analyses and their data subsets. A full range is a completed, uninformative result; an unavailable value is not zero.

| Agent | Analysis | N | Status | Data scope |
|---|---|---:|---|---|
| cefotaxime | mic_order | 472 | computed | readable and typed MIC subset |
| cefotaxime | mic_latent_bounds | 472 | computed | readable and typed MIC subset |
| cefotaxime | mic_latent_lower_limit | 472 | computed | readable and typed MIC subset |
| ceftazidime | mic_order | 472 | computed | readable and typed MIC subset |
| ceftazidime | mic_latent_bounds | 472 | computed | readable and typed MIC subset |
| ceftazidime | mic_latent_lower_limit | 472 | computed | readable and typed MIC subset |
| ciprofloxacin | mic_order | 472 | computed | readable and typed MIC subset |
| ciprofloxacin | mic_latent_bounds | 472 | computed | readable and typed MIC subset |
| ciprofloxacin | mic_latent_lower_limit | 472 | computed | readable and typed MIC subset |
| gentamicin | mic_order | 472 | computed | readable and typed MIC subset |
| gentamicin | mic_latent_bounds | 472 | computed | readable and typed MIC subset |
| gentamicin | mic_latent_lower_limit | 472 | computed | readable and typed MIC subset |
| meropenem | mic_order | 472 | full_range | readable and typed MIC subset |
| meropenem | mic_latent_bounds | 472 | computed | readable and typed MIC subset |
| meropenem | mic_latent_lower_limit | 472 | computed | readable and typed MIC subset |
| tetracycline | mic_order | 472 | computed | readable and typed MIC subset |
| tetracycline | mic_latent_bounds | 472 | computed | readable and typed MIC subset |
| tetracycline | mic_latent_lower_limit | 472 | computed | readable and typed MIC subset |

**Table 0b.** Reasons and next checks for results requiring interpretation. These checks do not change the recorded status or justify changing reporting conditions.

| Agent | Analysis / status | Reason | Next check |
|---|---|---|---|
| cefotaxime | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.582 | Read the interval kind, assumptions and retained collection before comparing results. |
| ceftazidime | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.618 | Read the interval kind, assumptions and retained collection before comparing results. |
| ciprofloxacin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.830 | Read the interval kind, assumptions and retained collection before comparing results. |
| gentamicin | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.743 | Read the interval kind, assumptions and retained collection before comparing results. |
| meropenem | mic_order / full_range | The completed interval or bounds span the entire admissible range [0, 1] | Report the full range as uninformative; inspect repeated-lineage support and outcome variation. |
| meropenem | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.993 | Read the interval kind, assumptions and retained collection before comparing results. |
| tetracycline | mic_latent_bounds / computed | The largest share attained by an arrangement of the lineages inside the readings is 0.974 | Read the interval kind, assumptions and retained collection before comparing results. |

> **Not evaluated on this run, and why**
>
> Decomposition of a prevalence difference into lineage composition and within-lineage rate: needs two collections; this record holds one.
>
> e-values: not computed on this run.
>

## 2. Admissibility of the input

- **Lineage column:** `st`
- **Susceptibility calls:** 0 antimicrobials on 481 isolates
- **Recorded dilutions:** 2886 rows read, 2886 joined; 481 of 481 isolates carry a recorded dilution
- **Resampling:** 5 folds, 20 repeats, 999 bootstrap draws, 999 permutations per antimicrobial

Lineage groups in the input: 75, of which 27 hold a single isolate. Support, the share of isolates in lineages of at least two, is 94.3 %; the share is scored on those isolates and the singletons are set aside. At least two lineages hold two or more isolates.

*Figure 1 (drawn on the page).* Isolates per lineage, largest first, the 40 largest of 75. A lineage of one isolate, drawn in orange, cannot be predicted out of sample and is set aside; support is the share of isolates in the other lineages, 94.3 % here. The largest lineage holds 15.5 % of the isolates, which sets how much one lineage can weigh in the share.

**Table 2.** Conditions the estimator requires before any result is reported.

| Condition | Observed | Required | Verdict |
|---|---:|---:|---|
| Lineages with at least two isolates, fewest over the antimicrobials read | 48 | ≥ 2 | accepted |
| Support, lowest over the antimicrobials read | 94.3 % | reported | singletons set aside |

## 3. The lineage share of the call, trait by trait

No per-trait share is reported on this run.

## 4. Evidence that survives re-reading

No e-values were computed on this run.

## 5. Reading at the recorded resolution

Where a dilution was recorded, the lineage question is also read from the dilution itself rather than from the call derived from it. Every reading is placed by its position among the readings of its own stratum, a level of `panel`: the share of those readings below it plus half the share at the same reading, with a censored reading placed by the nonparametric maximum-likelihood distribution of the readings. The share of these positions that lineage accounts for is estimated as the share of a call is, and its control shuffles the lineage labels only within strata. No distribution is assumed for the readings, a change of the MIC scale such as mg/L to log2 leaves the share unchanged, and a constant offset or a different panel between strata drops out. The share of the call and this share are read on independently retained readable and typed subsets; their denominators may differ.

The readings also bound the lineage share of the MIC ordering they were read from, the ordering of the latent MICs within a level of `panel`: the rank intraclass correlation of the MIC itself. A reading says only that the MIC lies in its range, so the readings do not identify that share; they bound it. The lower bound is the smallest share that any arrangement of the latent MICs within their readings allows, the lower bound the readings establish, computed with a certified error; the upper end is the largest share, attained by an arrangement, or, marked ≤, a certified bound on it. Every share between them is allowed. The figure after the bounds is the one-sided 95 % lower confidence limit of the lower bound, and so of the share: the isolates of every lineage are split at random, one half chooses the lineage scores on which the lower bound rests and the other half tests them, and then the other way round; the tests of 25 such splits are averaged, and their average spread sets the limit. A finer panel never lowers the lower bound, and a call, one cut of the panel, never raises it. Resolution is the share of the ordering the panel resolves, one less the expected squared width of a reading on the scale of positions; it is the factor by which the share of the midpoint arrangement falls short of the plug-in share of the scores.

**Table 5.** Lineage share of the MIC ordering per agent, with the 95 % interval for the lineages of this collection and the permutation p-value and the Benjamini-Yekutieli q-value across agents; selected is the Benjamini-Yekutieli selection at level 0.05. Control is the share on shuffled lineage labels. The bounds are those the readings place on the share of the latent MIC ordering, with the 95 % lower confidence limit of the lower bound. Readings on an end well are tied at the panel's edge and carry less of the ordering; their share is given beside the estimate. A p-value written with ≤ is the smallest the permutations can give. The last column gives the conclusion: established, lineage structure established, when the interval lies above zero and the agent is selected, or the lower confidence limit is above zero; not established; whole range, when the interval spans it; or not estimable. The conclusion says whether this collection shows lineage structure in the measurement, with the values behind it. It is not a statement about transmission or a resistance mechanism.

| Agent | Readings | Share | 95 % interval, represented lineages | p | q | Control | Selected | Latent ordering: bounds; 95 % lower limit | Resolution | End wells | Strata | Conclusion |
|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---|
| cefotaxime | 472 | 0.117 | 0.056 to 0.185 | ≤ 0.001 | 0.004 | 0.00 | yes | 0.006 to ≤ 0.868; ≥ 0.000 | 0.821 | 89.6 % | 1 | established |
| ceftazidime | 472 | 0.119 | 0.049 to 0.197 | ≤ 0.001 | 0.004 | 0.00 | yes | 0.056 to ≤ 0.954; ≥ 0.000 | 0.738 | 69.9 % | 1 | established |
| ciprofloxacin | 472 | 0.197 | 0.115 to 0.281 | ≤ 0.001 | 0.004 | 0.00 | yes | 0.025 to ≤ 1.000; ≥ 0.000 | 0.615 | 90.5 % | 1 | established |
| gentamicin | 472 | 0.094 | 0.023 to 0.176 | ≤ 0.001 | 0.004 | 0.00 | yes | 0.009 to ≤ 1.000; ≥ 0.000 | 0.621 | 98.7 % | 1 | established |
| meropenem | 472 | 0.497 | 0.000 to 1.000 | 0.003 | 0.009 | 0.00 | yes | 0.013 to ≤ 1.000; ≥ 0.000 | 0.027 | 100.0 % | 1 | whole range |
| tetracycline | 472 | 0.111 | 0.000 to 0.318 | 0.063 | 0.154 | 0.00 | no | 0.012 to ≤ 1.000; ≥ 0.000 | 0.123 | 99.6 % | 1 | not established |

**Table 5b.** Lineage structure behind each share of the MIC ordering: the lineages with at least two readings on which the share is scored, the readings of singleton lineages set aside, the effective number of scored lineages (inverse of the sum of squared lineage shares) and the share of readings in lineages with at least two members (support). Few effective lineages mean the share rests on a few lineages.

| Agent | Lineages scored | Singletons set aside | Effective | Support |
|---|---:|---:|---:|---:|
| cefotaxime | 48 | 27 | 15.9 | 94.3 % |
| ceftazidime | 48 | 27 | 15.9 | 94.3 % |
| ciprofloxacin | 48 | 27 | 15.9 | 94.3 % |
| gentamicin | 48 | 27 | 15.9 | 94.3 % |
| meropenem | 48 | 27 | 15.9 | 94.3 % |
| tetracycline | 48 | 27 | 15.9 | 94.3 % |

*Figure 2 (drawn on the page).* Lineage share of the MIC ordering per agent, 6 largest of 6: the point estimate with its 95 % interval; agents selected across the panel are drawn in blue. The hollow diamond is the share of the binary call from Table 1 for the same agent, a different quantity whose retained isolates may differ, placed here so that the two readings can be seen side by side. The figure in parentheses at the right is the share of readings on an end well.

Each laboratory is read on the wells that laboratory tested: an end-well reading is censored at the panel edge of its own laboratory, not at the widest edge in the collection. The table below gives, per agent and laboratory, the panel the readings were taken to come from.

**Table 5a.** Panel geometry per agent and testing laboratory: the wells taken as tested, their range, whether they form a doubling series, and the shares of readings on the lowest and the highest well.

| Agent | Laboratory | Wells | Range | Doubling | On lowest well | On highest well |
|---|---|---:|---:|---|---:|---:|
| cefotaxime | BD Phoenix | 6 | 1.0000 to 32.0000 | yes | 47.2 % | 44.7 % |
| ceftazidime | BD Phoenix | 5 | 1.0000 to 16.0000 | yes | 64.2 % | 13.6 % |
| ciprofloxacin | BD Phoenix | 3 | 0.5000 to 2.0000 | yes | 18.9 % | 75.0 % |
| gentamicin | BD Phoenix | 3 | 2.0000 to 8.0000 | yes | 28.4 % | 71.0 % |
| meropenem | BD Phoenix | 4 | 1.0000 to 8.0000 | yes | 99.2 % | 0.8 % |
| tetracycline | BD Phoenix | 3 | 2.0000 to 8.0000 | yes | 3.8 % | 95.8 % |

*Figure 3 (drawn on the page).* Readings per lineage and dilution interval for 6 of 6 agents, the 14 largest of 75 lineages: each cell counts the isolates of one lineage whose reading fell in one interval of the panel, darker for more. Open intervals at the panel edges are the censored readings. A lineage whose readings sit in one or two adjacent cells contributes to the lineage share; readings spread along a row do not.

## 6. A change of lineages or a change within them

A prevalence difference between two collections can be split into a change in lineage composition and a change in rate within lineages. That decomposition needs two collections and was not run here: this record holds one.

## 7. Provenance and terms

- **Software:** amr-clonalshare 1.0.0, record schema 1.0
- **Seed:** 20261001
- **Configuration:** `sha256:3ae0b5a290d7`
- **Record:** `sha256:2d32bfbb18333bf880cd99ad82b5a1276de3219db9fdc80e51443dff2ea797da`

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

amr-clonalshare 1.0.0 · record schema 1.0 · seed 20261001 · configuration sha256:3ae0b5a290d7
Cite this run as: “amr-clonalshare 1.0.0, run 424A3AD1, record sha256:2d32bfbb1833.”
This report supersedes any earlier report bearing the same run identifier. It is regenerated from the record and holds no value that the record does not. The symbols carry the same wording in every run of this software; no symbol against a value means only that none of the listed conditions fired.
