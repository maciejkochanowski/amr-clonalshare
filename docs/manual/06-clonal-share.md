# 6. The clonal share, agent by agent

## The question

For one antimicrobial, take the non-wild-type call of every isolate. Some of
the variation in that call sits *between* lineages (one lineage is mostly
resistant, another mostly susceptible) and some sits *within* them (isolates
of the same lineage differ). The clonal share is the between part as a
fraction of the whole. Near 1: the resistance travels with the clone, and the
intervention is against the clone. Near 0: the resistance moves across
lineages, and the intervention is against the selecting agent.

## Three shares, three questions

The tool reports three quantities and does not let one interval stand for
another. The first two are both statements about the collection in hand and
differ by the design factor `1 - sum w_g^2`, which is `(G - 1) / G` when the
lineages are of equal size: at five lineages a realised share of 0.50 is a
lineage-membership share of 0.44, at fifty it is 0.49.

| | question | interval | precision grows with |
|---|---|---|---|
| **lineage-membership share** | in *this* collection, how much of the variance does knowing the lineage carry | lineage bootstrap (`clonal_share`) | isolates per lineage |
| **realised share** | in *this* collection, holding *these* lineages, how much sits between them as a variance component | exact, from the noncentral F (`realised_share`) | isolates per lineage |
| **superpopulation share** | if new lineages were drawn from the species, how much would sit between them | variance-component F (`realised_share`); within-lineage bootstrap under a chi-square layer on the number of lineages (`clonal_share`) | number of lineages |

A veterinary collection usually holds few lineages and many isolates, which is
the case where the two diverge most. The realised share is what a laboratory
holding one collection usually means.

Two cases the report names. With fewer than ten lineages the lineage
bootstrap of `clonal_share` is not the interval to read: on the validation
grid it contained the truth in only 0.83 to 0.96 of runs at five lineages,
while the species interval held its level, so the record flags the collection
(`few_lineages`) and the report says to read the species interval. And where
one lineage carries more than half of the between-lineage variation
(`dominant_lineage_share`), the species interval rests on a law the collection
does not display, a normal spread of lineage effects rather than one carrier
among many; on the grid a carrier law of that kind took the species interval,
the variance component and the mixed model alike below their level with few
lineages, while the interval for the lineages in hand held. The report names
the traits concerned and says to read the first interval as the statement
about the collection.

## Why the estimate is out of sample

A lineage label with many levels can "explain" a trait by memorising it. The
clonal share is therefore estimated on isolates the lineage model has not
seen, in repeated cross-validation, with the fold means debiased. A label
that carries no information scores zero, not the number of its levels. Two
intervals follow: one drawn by resampling whole lineages, for the share these
lineages carry, and one that adds the draw of the lineages themselves, for
the species; with ten lineages the second is markedly wider, with a hundred
the two nearly coincide. On the validation grid the first holds its level
from ten lineages up, and the second covers 0.81 to 1.00 under the normal law,
lowest on a binary call at 8 % prevalence.

## The gates, and what they mean for you

**Support below 0.90.** Too many isolates sit in single-member lineages for
the within-lineage variation to be measured. The cell is reported as not
estimable. Use a coarser lineage definition, or accept that this collection
cannot answer the question at this resolution.

**Residual kurtosis above the limit.** The exact interval for the realised
share rests on the within-lineage departures having tails no heavier than the
Gaussian; the gate is one-sided and was read off the measured coverage curve.
A binary call passes it only where its prevalence lies between about 17 and
83 per cent, because the excess kurtosis of a Bernoulli variable rises without
bound as its prevalence leaves one half; outside that band the gate closes,
which is the usual case for a rare resistance, and the record prints the
kurtosis beside the refusal. The interval's coverage on binary calls inside
the band is reported separately in the validation ledger.

## Reading the per-trait table

Table 1 of `report.md`, on its first page, lists every trait with its clonal
share, the 95 % interval, the permuted-label control, the e-value and a
reading, sorted by share, and Table 3 sets the two intervals side by side. The
interval is the statement: a share whose interval sits below 0.5 says the
trait moves mostly across lineages in this collection; one whose interval sits
above 0.5 says it is mostly a lineage trait here; one that spans 0.5 says the
collection cannot yet tell, and the input check says which of the two binding
quantities would sharpen it.

![Section 3 of the report on the shipped cohort: the intervals drawn](../img/06_clonal_share.png)

The record also carries a permutation p-value per trait, for readers who want
one. It answers a narrower question than the interval: whether the lineage
carries any information about the trait at all. The statistic on both sides
of the comparison is the same function of the labels, a mean over the same
number of fold draws, so the p-value is exact under exchangeability rather
than conservative; `benchmarks/null_uniformity.py` checks that it is uniform
on the permutation grid when there is nothing to find, and records the power
against a graded lineage effect. It cannot fall below one over the number of
permutations plus one, and the record carries that floor (`p_floor`); a
p-value at the floor means no permutation reached the observed share, and the
report says so instead of printing the number as if it were an estimate.

For a binary trait the record also gives the share and its two intervals on
the latent scale of a threshold model (`latent_share` and the four endpoints
beside it), which is the scale on which a mixed model with a binomial link,
as `rptR` or `lme4` fit it, reports an intraclass correlation. The table in
section 3 of the report prints it in its last column. The two scales are
one-to-one at a given prevalence and the transformation is monotone, so the
interval keeps the coverage that was measured for it; the latent figure is
the larger of the two, because a susceptibility call discards the part of the
underlying liability that does not cross the threshold. Read the observed
scale for a surveillance table and the latent one to compare with a mixed
model.

## The e-value

Beside each share is an e-value: the evidence that the lineage carries *any*
information about the trait, from this run read as one look. An e-value of 20
corresponds to a decision at the 0.05 level, 100 to the 0.01 level.
False-discovery control across the panel uses the e-BH procedure, which holds
under any dependence between antimicrobials.

A programme that re-reads the panel at every intake should not recompute this
value on the growing cohort and read the sequence as one test: each value is
valid at its own look, but the sequence is not a test martingale. The rule for
a programme is `evalues.sequential_e_process`, which the run applies when
`dataset.batch_column` names the intake: keep the intakes as separate
batches in arrival order, and at each new batch it fits the lineage
probabilities on every earlier batch and scores the new one, multiplying the
factors. The running product may be inspected after any intake and stopped at
any time; an intake's isolates must not appear in a later batch.

```python
from amr_clonalshare.evalues import sequential_e_process, e_bh
seq = sequential_e_process([calls_2019, calls_2020, calls_2021],
                           [types_2019, types_2020, types_2021])
seq.e_value[-1]        # the running evidence after the last intake
seq.reject_05          # whether it has ever reached 20
```

## What the share is not

It is not heritability, which splits variance against a relatedness
matrix estimated from the genome; the random effect here is a discrete label
whose resolution is an analytical choice. It is not the accuracy of a
lineage-based classifier, which is bounded below by the majority class: a cell
at 90 % prevalence scores 0.90 with no lineage information at all, while the
clonal share returns zero.
