# amr-clonalshare run report

Every number below is read from `clonal_share_result.json`, written by the same run; none is recomputed here.

- **Run:** 83923CE1
- **Issued:** 2026-10-03T05:58Z
- **Software:** amr-clonalshare 1.0.0
- **Isolates:** 96
- **Lineage:** lineage
- **Seed:** 42
- **Configuration:** sha256:9052f9b38cdf
- **Record digest:** sha256:0432cfb49c05bfb0d45db84cc98a97249c7db483254cbc452066148526d8cfdb

**1 antimicrobials; 3 analysis outcomes.** 3 computed, 0 full ranges, 0 unavailable. Each result retains its own target and data subset.

## 1. Measurement summary

**Quantity measured.** The share of the ordering of the recorded MIC readings within each stratum that lineage membership accounts for, and the lower bound on the lineage share of the latent MIC ordering established by the readings.

- **Collection:** 96 isolates in 12 lineages, typed by `lineage`
- **Antimicrobials read:** 1

**Reporting conditions.** 3 computed, 0 full ranges, 0 unavailable. Each result retains its own target and data subset. An unavailable value failed a declared reporting condition of its estimator; the recorded reason distinguishes input limitations from numerical failure.

**Table 0.** Available analyses and their data subsets. A full range is a completed, uninformative result; an unavailable value is not zero.

| Agent | Analysis | N | Status | Data scope |
|---|---|---:|---|---|
| demo_agent | mic_order | 96 | computed | readable and typed MIC subset |
| demo_agent | mic_latent_bounds | 96 | computed | readable and typed MIC subset |
| demo_agent | mic_latent_lower_limit | 96 | computed | readable and typed MIC subset |

> **Not evaluated on this run, and why**
>
> Decomposition of a prevalence difference into lineage composition and within-lineage rate: needs two collections; this record holds one.
>
> e-values: not computed on this run.
>

## 2. Admissibility of the input

- **Lineage column:** `lineage`
- **Susceptibility calls:** 0 antimicrobials on 96 isolates
- **Recorded dilutions:** 96 rows read, 96 joined; 96 of 96 isolates carry a recorded dilution
- **Resampling:** 2 folds, 2 repeats, 20 bootstrap draws, 19 permutations per antimicrobial

Lineage groups in the input: 12, of which 0 hold a single isolate. Support, the share of isolates in lineages of at least two, is 100.0 %; the share is scored on those isolates and the singletons are set aside. At least two lineages hold two or more isolates.

*Figure 1 (drawn on the page).* Isolates per lineage, largest first. A lineage of one isolate, drawn in orange, cannot be predicted out of sample and is set aside; support is the share of isolates in the other lineages, 100.0 % here. The largest lineage holds 8.3 % of the isolates, which sets how much one lineage can weigh in the share.

**Table 2.** Conditions the estimator requires before any result is reported.

| Condition | Observed | Required | Verdict |
|---|---:|---:|---|
| Lineages with at least two isolates, fewest over the antimicrobials read | 12 | ≥ 2 | accepted |
| Support, lowest over the antimicrobials read | 100.0 % | reported | singletons set aside |

## 3. The lineage share of the call, trait by trait

No per-trait share is reported on this run.

## 4. Evidence that survives re-reading

No e-values were computed on this run.

## 5. Reading at the recorded resolution

Where a dilution was recorded, the lineage question is also read from the dilution itself rather than from the call derived from it. Every reading is placed by its position among the readings of the collection: the share of those readings below it plus half the share at the same reading, with a censored reading placed by the nonparametric maximum-likelihood distribution of the readings. The share of these positions that lineage accounts for is estimated as the share of a call is, and its control shuffles the lineage labels across the collection. No distribution is assumed for the readings, a change of the MIC scale such as mg/L to log2 leaves the share unchanged, and the readings of one panel are compared with one another. The share of the call and this share are read on independently retained readable and typed subsets; their denominators may differ.

The readings also bound the lineage share of the MIC ordering they were read from, the ordering of the latent MICs within the collection: the rank intraclass correlation of the MIC itself. A reading says only that the MIC lies in its range, so the readings do not identify that share; they bound it. The lower bound is the smallest share that any arrangement of the latent MICs within their readings allows, the lower bound the readings establish, computed with a certified error; the upper end is the largest share, attained by an arrangement, or, marked ≤, a certified bound on it. Every share between them is allowed. The figure after the bounds is the one-sided 95 % lower confidence limit of the lower bound, and so of the share: the isolates of every lineage are split at random, one half chooses the lineage scores on which the lower bound rests and the other half tests them, and then the other way round; the tests of 25 such splits are averaged, and their average spread sets the limit. A finer panel never lowers the lower bound, and a call, one cut of the panel, never raises it. Resolution is the share of the ordering the panel resolves, one less the expected squared width of a reading on the scale of positions; it is the factor by which the share of the midpoint arrangement falls short of the plug-in share of the scores.

**Table 5.** Lineage share of the MIC ordering per agent, with the 95 % interval for the lineages of this collection and the permutation p-value and the Benjamini-Yekutieli q-value across agents; selected is the Benjamini-Yekutieli selection at level 0.05. Control is the share on shuffled lineage labels. The bounds are those the readings place on the share of the latent MIC ordering, with the 95 % lower confidence limit of the lower bound. Readings on an end well are tied at the panel's edge and carry less of the ordering; their share is given beside the estimate. A p-value written with ≤ is the smallest the permutations can give. The last column gives the conclusion: established, lineage structure established, when the interval lies above zero and the agent is selected, or the lower confidence limit is above zero; not established; whole range, when the interval spans it; or not estimable. The conclusion says whether this collection shows lineage structure in the measurement, with the values behind it. It is not a statement about transmission or a resistance mechanism.

| Agent | Readings | Share | 95 % interval, represented lineages | p | q | Control | Selected | Latent ordering: bounds; 95 % lower limit | Resolution | End wells | Strata | Conclusion |
|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---|
| demo_agent | 96 | -0.177 | 0.000 to 0.000 | 0.950 | 0.950 | -0.03 | no | 0.001 to 0.119; ≥ 0.000 | 0.972 | 33.3 % | 1 | not established |

**Table 5b.** Lineage structure behind each share of the MIC ordering: the lineages with at least two readings on which the share is scored, the readings of singleton lineages set aside, the effective number of scored lineages (inverse of the sum of squared lineage shares) and the share of readings in lineages with at least two members (support). Few effective lineages mean the share rests on a few lineages.

| Agent | Lineages scored | Singletons set aside | Effective | Support |
|---|---:|---:|---:|---:|
| demo_agent | 12 | 0 | 12.0 | 100.0 % |

*Figure 2 (drawn on the page).* Lineage share of the MIC ordering per agent, 1 largest of 1: the point estimate with its 95 % interval; agents selected across the panel are drawn in blue. The hollow diamond is the share of the binary call from Table 1 for the same agent, a different quantity whose retained isolates may differ, placed here so that the two readings can be seen side by side. The figure in parentheses at the right is the share of readings on an end well.

Each laboratory is read on the wells that laboratory tested: an end-well reading is censored at the panel edge of its own laboratory, not at the widest edge in the collection. The table below gives, per agent and laboratory, the panel the readings were taken to come from.

**Table 5a.** Panel geometry per agent and testing laboratory: the wells taken as tested, their range, whether they form a doubling series, and the shares of readings on the lowest and the highest well.

| Agent | Laboratory | Wells | Range | Doubling | On lowest well | On highest well |
|---|---|---:|---:|---|---:|---:|
| demo_agent | all readings | 6 | 0.2500 to 8.0000 | yes | 16.7 % | 16.7 % |

*Figure 3 (drawn on the page).* Readings per lineage and dilution interval for 1 of 1 agents: each cell counts the isolates of one lineage whose reading fell in one interval of the panel, darker for more. Open intervals at the panel edges are the censored readings. A lineage whose readings sit in one or two adjacent cells contributes to the lineage share; readings spread along a row do not.

## 6. A change of lineages or a change within them

A prevalence difference between two collections can be split into a change in lineage composition and a change in rate within lineages. That decomposition needs two collections and was not run here: this record holds one.

## 7. Provenance and terms

- **Software:** amr-clonalshare 1.0.0, record schema 1.0
- **Seed:** 42
- **Configuration:** `sha256:9052f9b38cdf`
- **Record:** `sha256:0432cfb49c05bfb0d45db84cc98a97249c7db483254cbc452066148526d8cfdb`

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

amr-clonalshare 1.0.0 · record schema 1.0 · seed 42 · configuration sha256:9052f9b38cdf
Cite this run as: “amr-clonalshare 1.0.0, run 83923CE1, record sha256:0432cfb49c05.”
This report supersedes any earlier report bearing the same run identifier. It is regenerated from the record and holds no value that the record does not. The symbols carry the same wording in every run of this software; no symbol against a value means only that none of the listed conditions fired.
