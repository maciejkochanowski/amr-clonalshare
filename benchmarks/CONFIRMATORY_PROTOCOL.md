# Confirmatory campaign: protocol

Written and committed before any run of the campaign. It fixes what every
study measures, on which datasets and with how many of them, and the rule that
decides it. A result that fails its rule is reported as a failure, and the
claim it would have supported is withdrawn or restricted to the designs that
pass. `confirmatory_verdicts.py` applies the rules to the outputs of the runs
and names every design that fails one. The runs use the package at the commit
that carries this file; the commit is recorded in the receipt of every run
(`campaign/run_logged.py`).

## 1. What is tested

The package estimates the lineage share of a call or of the ordering of MIC
readings among the lineages in hand, tests it by permutation, and bounds the
lineage share of the unobserved MIC ordering. Five procedures carry claims
that can fail on data, and each is tested here against a truth written from
its definition, in code that shares nothing with the package:

1. **The interval for the lineages in hand** (`observed_low`,
   `observed_high`): the studentized bootstrap of
   `attribution._smoothed_bootstrap`, lineages and their sizes fixed, the
   isolates of every lineage drawn again, studentized by the jackknife
   variance plus the variance of the share without lineage differences,
   its quantiles taken no closer to the estimate than those of the
   studentized shares of the permuted labellings; its upper end raised to
   the largest share one lineage can carry within the empirical-likelihood
   region of its readings (`attribution._single_lineage_gain`). Target: the
   share of the lineages of the dataset, computed from their laws.
2. **The permutation test** (`p_value`): lineage labels exchanged within
   strata, and within sampling units where they are declared. Target: its
   level where no lineage effect exists within strata and units.
3. **The bounds on the share of the unobserved ordering** (Theorem 1): the
   share of the ordering of the latent values of the isolates in hand lies
   within `[latent_order_lower, latent_order_upper_bound]` in every dataset.
   This is a theorem; one exception refutes the implementation.
4. **The lower confidence limit** (`latent_order_lower_limit`): below the
   share the readings guarantee in the population of the lineages in hand
   (`rho_min`), and so below the latent share `rho`.
5. **The paired intervals of the comparison of two lineage definitions**
   (`paired_intervals`): records and labels fixed, outcomes drawn again.

The prevalence-difference decomposition, the e-values and the uniformity of
the permutation p-value are unchanged from release 1.0.0 in every path the
defaults take and are not rerun; their results folders are kept with their
receipts.

## 2. Studies, designs and replicates

Every replicate of every design is drawn from its own stream, a
`numpy.random.SeedSequence` of the study's root seed, the design and the
replicate; the root seeds below have not been used before. No development
pilot of this revision used them; the pilots ran on seeds of their own and no
pilot output is shipped or supports a claim.

**Study 1, calls and traits (`estimator_benchmark.py`, root seed
20261004101).** The grid of 175 cells: binary calls and normal traits; 5 to
1,000 lineages of 2 to 20 isolates; balanced, Poisson, half-pairs and
all-pairs sizes; normal and two-point (carrier) lineage effects; shares 0 to
0.7 of the latent variance; prevalence 0.25 and 0.08. Lineage effects are
drawn in every replicate and the targets are those of the lineages drawn.
2,000 replicates per cell. Measured: procedures 1 and, for calls, 3 and 4;
bias and root mean squared error of the package's estimate and of three
alternatives (in-sample R², the ANOVA variance component, REML) against
their own targets.

**Study 2, the MIC ordering share (`order_calibration.py`, root seed
20261002207).** The 23 designs of that file: Gaussian lineage effects at
shares 0.1, 0.5 and 0.9 of the latent variance; heavy censoring; eight
lineages; singleton-rich; a two-point clone; t4, contaminated and two-mode
residuals; the *S. suis* lineage sizes with two modes (wild type and non-wild
type); sizes ordered as the effects; two laboratories on different panels;
lineages of two isolates; and six designs without a lineage effect, with
laboratories that differ in offset, panel or spread and lineages aligned with
them. 5,000 replicates per design with an effect, 20,000 per design without.
Measured: procedures 1, 2 (within laboratories; the test that ignores the
laboratory is reported beside it), 3 and 4 against `rho`.

**Study 3, fixed laws (`conditional/`, root seed 20261003311).** Lineages,
their sizes and their laws fixed, only the isolates drawn again, so every
target is exact:

- *observed* (32 designs, 10,000 replicates each): near-deterministic
  lineages on two to ten wells (the audit's design among them), two
  laboratories with the laboratory drawn with the reading, no effect, weak
  and moderate effects on eight wells, uneven sizes, pairs, pairs with a
  weak effect on twelve wells and on a practically continuous panel of
  1,000 wells (5, 10 and 30 pairs, and 100 pairs with a moderate effect),
  and the sparse regime: a call rare in all lineages but one or three,
  which carry it at a rate of one half, on lineages of 5, 10 and 20
  isolates. Procedure 1.
- *latent* (21 designs, 10,000 replicates each): normal, t4 and
  contaminated latent laws; the audit's null designs (two lineages of two on
  a median call; 30 and 100 pairs on 9 and on 1,000 wells); near-null and
  effect designs; a call; two laboratories with confined lineages; the rare
  calls; and two designs whose laboratory reads on two overlapping panels.
  Procedures 3 and 4 (4 against `rho_min` where one panel per laboratory
  makes it computable, against `rho` everywhere).
- *units* (5 designs, 20,000 replicates each): farms that shift the
  readings of their isolates, lineages mostly or wholly on one farm, no
  lineage effect within farms; panels and calls, nested and crossed with
  laboratories, and pairs. Procedure 2 within units (the test that ignores
  the units is reported beside it).
- *comparison* (5 designs, 5,000 replicates each; a sixth added by the amendment of Section 5): a second definition that
  merges lineages in pairs or fives, or splits every lineage in two, with a
  tenth of the records missing the second label; shares 0, 0.4 and 0.6;
  prevalence 0.25 and 0.08. Procedure 5, every term whose target is not zero
  by construction.

**Study 4, published procedures (`competitors/compare_order.py`, the
generator of Study 2 from replicate 2,000,000 on).** Seven published
procedures for a lineage effect on MIC readings, run in R on the same
datasets as the test of the ordering share: 8 designs without an effect
(5,000 replicates each) and 8 with a weak one (2,000 each).

## 3. Decision rules

`D` below is the number of designs (cells) of one study and one procedure.
All Clopper–Pearson (CP) bounds are one-sided at level `0.05 / D`
(Bonferroni over the designs of the study).

**Coverage of a 95 % interval (procedures 1 and 5).** A design *meets* the
rule if the CP lower bound of its coverage is at least 0.925, the liberal
criterion of Bradley (1978) for a nominal 0.95; it *falls short* if the CP
upper bound is below 0.925; otherwise it is *undecided*: a cell of Study 1 is then
run again with 10,000 replicates on the same stream, whose first 2,000 are
those already scored, and decided on them; a design of Studies 2 and 3,
whose replicate counts already resolve a difference of 0.01, is reported as
undecided. Coverage at or
above 0.945 by the same CP bound is reported as meeting the stringent
criterion. Coverage is counted over the datasets on which the interval was
reported; the share of datasets without an interval, and why, is reported.
Over-coverage is not a failure; the median width is reported.

**Level of the test (procedure 2).** A design without an effect *meets* the
rule if the CP upper bound of the rejection rate at the 5 % level is at most
0.075 (Bradley's liberal criterion); at most 0.055 is reported as stringent.

**Theorem 1 (procedure 3).** Every dataset of every design must have the
share of the latent values in hand within the bounds, up to 1e-9. Any
exception is a failure of the implementation, whatever the design.

**The lower limit (procedure 4).** A design *meets* the rule if the CP upper
bound of the rate at which the limit exceeds its target is at most 0.075; at
most 0.055 is reported as stringent. The rate at which the limit is positive
and the median ratio of the limit to its target are reported as its power.

**Study 4.** H1: the test of the ordering share meets the level rule above
in every design without an effect. H2: in every design with a weak effect,
its rejection rate is not below that of the stratified Wilcoxon test by more
than one percentage point (lower limit of the two-sided 95 % interval of the
paired difference at least -0.01). Differences against every procedure are
reported whichever way they fall.

**Study 1, the estimate.** Reported, not decided: bias and root mean squared
error of every estimator against its own target and the other. What would
mean the grid measured the wrong thing, checked and reported: that the
in-sample R² is not the most biased estimator at high lineage counts (the
generator would not produce the regime of the comparison), or that the
estimators never differ by more than 0.05 (the grid would have no leverage).

## 4. Reporting

Every design is reported with its coverage or rate, its CP bound, the number
of replicates and the Monte Carlo standard error; the designs that fail are
named in the manuscript and in the supplement, and the claims are restricted
to the designs that pass. A change to the package after this protocol is
committed voids the campaign for the procedures it touches; a deviation in
how the runs are carried out is reported beside the results it affects.

## 5. Amendment after the runs

Written after the runs of Studies 1 to 4 had started and before any of their
outputs or verdicts was read. On the *Salmonella* example
(`empirical/salmonella_analysis.py`), the paired interval of the relabelling
term did not contain its own estimate. The cause: the centre of the drawn world
(`comparison._arm_truth`) and the truth of the comparison designs
(`conditional/run.py`, `_arm`) gave every lineage the mean of the Bernoulli
variances of its records, leaving out the spread of the rates of the cells
within the lineage, where the share (procedure 1) gives a lineage the variance
of its own law, that of a record drawn from it. Where the cells of a lineage
share one rate the two agree, so procedures 1 to 4 are untouched; where a
lineage holds cells of different rates, as a serovar of several SNP clusters
does, the drawn world was centred on a share that its estimate does not
estimate, and in the designs above the truth was that same share. Both are
corrected to the definition; arms that score the same records under the same
labels are now one analysis, so that a term between them is zero in every
draw. Procedure 5 is therefore run again from the commit that carries this
amendment, on the same streams of its five designs, with a sixth design added
from the example (`small_rare_merge_fives`: 200 lineages of three isolates,
most without a positive, merged in fives, prevalence 0.08; 5,000 replicates),
by `campaign/comparison_rerun.txt`. The rule is unchanged, with D the number of
scored terms of the six designs. The verdicts of procedure 5 are those of the
rerun; the outputs of the first run are kept beside them.

Bradley, J. V. (1978). Robustness? *British Journal of Mathematical and
Statistical Psychology*, 31, 144–152.

## 6. Note on wording (2026-10-02)

Added at the release, after the verdicts were read; it changes no rule, design
or seed. "Release 1.0.0" in Section 1 denotes the package as it stood before
this campaign, whose results folders are kept with their receipts; the version
number of the published package is 1.0.0. "This revision" in Section 2 denotes
the campaign described here. The commits named in the receipts are in the git
history attached to the release (`history.bundle`, see REPRODUCIBILITY.md).
