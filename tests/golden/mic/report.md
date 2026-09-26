# amr-clonalshare run report

Every number below is read from `clonal_share_result.json`, written by the same run, or from the release's validation grid where the text says so; none is recomputed here.

- **Run:** E682176B
- **Issued:** 2026-09-26T01:26Z
- **Software:** amr-clonalshare 1.0.0
- **Isolates:** 96
- **Lineage:** lineage
- **Seed:** 42
- **Configuration:** sha256:0393bc31c696
- **Record digest:** sha256:356e72bd7dc75d9e534472e1e20e811cc4ec80c9fa8c033aede394158f5d5cf7

**1 antimicrobials; 1 analysis outcomes.** 1 computed, 0 full ranges, 0 unavailable or incomplete. Each result retains its own target and data subset.

## 1. Measurement summary

**Quantity measured.** Interval-censored MIC variation associated with recorded lineage. The recorded panel and coarsening assumptions determine the interval interpretation.

- **Cohort:** 96 isolates in an unrecorded number of lineages, typed by `lineage`
- **Antimicrobials read:** 1

**Gates.** 1 computed, 0 full ranges, 0 unavailable or incomplete. Each result retains its own target and data subset. A withheld value failed a declared reporting condition of its estimator; the recorded reason distinguishes input limitations from numerical failure. Passing a gate does not verify the model.

**Table 0.** Available analyses and their data subsets. A full range is a completed, uninformative result; an unavailable value is not zero.

| Agent | Analysis | N | Status | Data scope |
|---|---|---:|---|---|
| demo_agent | censored_mic | 96 | computed | population model fitted to readable and typed MIC subset |

**Table 0b.** Reasons and next checks for results requiring interpretation. These checks do not change the recorded status or justify changing reporting gates.

| Agent | Analysis / status | Reason | Next check |
|---|---|---|---|
| demo_agent | censored_mic / computed | ok | Read the interval kind, assumptions and retained cohort before comparing results. |

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
| Lineages with at least two isolates, fewest over the antimicrobials read | not recorded | ≥ 2 | not evaluated |
| Support, lowest over the antimicrobials read | not recorded | reported | not evaluated |

## 3. The clonal share, trait by trait

No per-trait share is reported on this run.

## 4. Evidence that survives re-reading

No e-values were computed on this run.

## 5. Reading at the recorded resolution

Where a dilution was recorded, the share is also read from the dilution itself rather than from the call derived from it. A call and a dilution are both intervals on one concentration scale, so the dilution carries more of the reading; a reading at an end well is censored and enters as an interval. This table reads the share of that scale carried by lineage, beside the share of the binary call above. These are different quantities on independently retained readable and typed subsets; their denominators may differ.

**Table 5.** Share of the dilution scale carried by lineage, per agent, from the interval-censored likelihood. The share and its approximate F interval come from the moment iteration; that interval undercovered in simulation. The calibrated column gives the maximum-likelihood share with the interval validated by simulation for the Gaussian population model, subject to its sampling and coarsening assumptions. The realised interval is for the lineages this cohort holds; here it is an approximation built on the fitted mean-square ratio and not the exact interval of uncensored data.

| Agent | Share | Approximate F interval | Calibrated share and 95 % interval | Realised interval | Censored readings | Gaussian model check | Estimable |
|---|---:|---:|---:|---:|---:|---|---|
| demo_agent | 0.013 | 0.000 to 0.223 | not computed | 0.000 to 0.169 | 33.3 % | not run | yes |

*Figure 2 (drawn on the page).* Share of the dilution scale carried by lineage, per agent, 1 largest of 1: the point estimate with its 95 % interval (thin) and the realised interval for these lineages (thick); an agent whose reading is refused is drawn in grey. No calibrated interval was computed on this run, so every interval drawn is the approximate F interval beside the moment estimate; that interval undercovered in simulation. The hollow diamond is the share of the binary call from Table 1 for the same agent, a different quantity whose retained isolates may differ, placed here so that the two readings can be seen side by side. The figure in parentheses at the right is the share of readings at an end well, which enter as censored.

The calibrated interval reads each laboratory on the wells that laboratory tested. The table below gives, per agent and laboratory, the panel the readings were taken to come from and the cut points the model uses; an end-well reading is censored at the panel edge of its own laboratory, not at the widest edge in the collection.

**Table 5a.** Panel geometry per agent and testing laboratory. Cut points are the boundaries between adjacent wells on which the likelihood is written; the admissible modes say which readings of the end wells the panel supports.

| Agent | Laboratory | Wells | Range | Cut points | Doubling | Admissible modes |
|---|---|---:|---:|---:|---|---|
| demo_agent | all readings | 6 | 0.2500 to 8.0000 | none | yes | binary and interval |

**Table 5c.** Status of the calibrated interval per agent. The status names the outcome of the interval inversion; the reason is given whenever no interval is reported.

| Agent | Status | Method | Draws | Failed refits | Reason |
|---|---|---|---:|---:|---|
| demo_agent | not computed |  |  |  | calibrated_interval is off in the configuration |

*Figure 3 (drawn on the page).* Readings per lineage and dilution interval for 1 of 1 agents: each cell counts the isolates of one lineage whose reading fell in one interval of the panel, darker for more. Open intervals at the panel edges are the censored readings. A lineage whose readings sit in one or two adjacent cells contributes to the lineage share; readings spread along a row do not.

## 6. A change of lineages or a change within them

A prevalence difference between two collections can be split into a change in lineage composition and a change in rate within lineages. That decomposition needs two collections and was not run here: this record holds one.

## 7. Provenance and terms

- **Software:** amr-clonalshare 1.0.0, record schema 1.0
- **Seed:** 42
- **Configuration:** `sha256:0393bc31c696`
- **Record:** `sha256:356e72bd7dc75d9e534472e1e20e811cc4ec80c9fa8c033aede394158f5d5cf7`

> **Terms used in this report**
>
> **Share.** How much of the variation of the call between isolates is associated with their recorded lineage. Near 1: a strong association on the measured scale. Near 0: little association on that scale. This does not identify transmission or its mechanism.
>
> **Collection-bootstrap interval.** The lineage-membership interval targets the sampled collection. Separate Gaussian-model component and species intervals retain their own targets and assumptions; coverage measured on a finite simulation grid is not a universal guarantee.
>
> **Control.** The same calculation on shuffled lineage labels, which should return roughly zero; a share is read against it.
>
> **Support.** The share of isolates in lineages with at least two members. A singleton lineage cannot be predicted out of sample, so the share is scored on the other isolates and the singletons are set aside; the share then speaks for singleton lineages only if they are drawn from the same population of lineages as the repeated ones.
>
> **Gate.** A condition fixed before the run. If it does not hold, the software withholds the estimate and names the condition that failed.
>
> **e-value.** A fixed-look e-value and a sequential e-process have different uses. Repeated-look validity requires the sequential construction and its stated null-model assumptions; a fixed-look value alone does not justify repeated inspection. Larger values are stronger evidence against the specified null.
>

## Notes on reading this report

- ‡ The 95 % interval includes zero. No lineage effect is distinguishable from none for this trait, and the point estimate must not be read on its own.
- § Derived arithmetically from a quantity recorded elsewhere in the record rather than estimated in this run.

amr-clonalshare 1.0.0 · record schema 1.0 · seed 42 · configuration sha256:0393bc31c696
Cite this run as: “amr-clonalshare 1.0.0, run E682176B, record sha256:356e72bd7dc7.”
This report supersedes any earlier report bearing the same run identifier. It is regenerated from the record and holds no value that the record does not. The symbols carry the same wording in every run of this software; no symbol against a value means only that none of the listed conditions fired.
Classical collection-bootstrap endpoints are printed without clipping; the separately labeled species interval is floored at zero as part of its construction.
