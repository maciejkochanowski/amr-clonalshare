# Methodology

Every statistic the package computes, its definition, the choice it embodies and
the citation behind it. Where a default is consequential the reasoning is given
here instead of in a docstring.

---

## 1. Two questions a clonal share can answer

`attribution.clonal_share` scores the means of the lineages in hand on
isolates they were not fitted on, so its point estimate is the
lineage-membership share of the collection, `B_w / (B_w + sigma^2)` with
`B_w` the isolate-weighted variance of those means, and its lineage bootstrap
(`ci_low`, `ci_high`) is an interval for that quantity rather than for the
design-corrected component ratio `realised_share` estimates: on the estimator
grid of `benchmarks/estimator_benchmark.py`
it covers its own target at 0.92 to 1.00 from ten lineages to three
hundred, on binary and continuous traits alike; at five lineages it covers
0.83 to 0.96, because a percentile bootstrap over five lineages has too few
distinct resamples, so the record flags a collection below `FEW_LINEAGES`
(ten) and the report names the species interval as the one to read there;
and at a thousand lineages of three isolates and a share of 0.7 it covers
0.885, a bias of -0.017 meeting an interval 0.055 wide, where the exact
realised interval is the one to read on a continuous trait. The record says
which quantity the interval targets in `interval_target`.

The share a *new* draw of lineages from the species would show is a second
question, and the same function reports a second interval for it,
`superpopulation_low` and `superpopulation_high`. Information about the
between-lineage variance carries `G - 1` degrees of freedom, where `G` is the
number of lineages, and no number of extra isolates adds to it; a lineage
bootstrap over ten lineages does not reproduce that law, which is why on the
grid it covers the superpopulation target at 0.71 to 0.96 with ten lineages,
0.82 on average. The second interval is built from what the cohort does
identify. A within-lineage bootstrap, resampling isolates inside each lineage
with the set of lineages fixed, measures the noise of the estimate around the
share these lineages carry; the draw of the lineages is then added as the
chi-square layer the one-way model implies, `S_a^2 = tau^2 * Q` with `Q` a
chi-square on `G - 1` degrees of freedom over `G - 1`, so that for a candidate
superpopulation share `rho` the realised share is `rho Q / (rho Q + 1 - rho)`.
The interval is the set of `rho` for which the observed estimate is not in
either 2.5 % tail of that predictive law, found by bisection
(`attribution._superpopulation_interval`). The within-lineage noise it adds is
the bootstrap measured at the observed share, and the same centred draws are
added at every candidate `rho`: the law of that noise is held independent of
the share, which it is not exactly, and the level of that approximation is what
the grid measures rather than something the construction guarantees. Neither source of uncertainty is
counted twice, and the construction assumes only that the lineage effects are
a draw whose dispersion follows the chi-square law, not that the trait is
Gaussian within lineages. A collection in which one lineage carries most of
the resistance is not such a draw, and there the lineage bootstrap, which
sees the influence of that lineage, is the wider statement; the reported
species interval is therefore the envelope of the two, which can only cover
more often than either. Its coverage is measured on the same grid and
reported in Table 2 of the article and section S14 of its supplement. The grid runs from five
lineages to a thousand, with three to twenty isolates a lineage, and holds
cells whose lineage effects follow a two-point carrier law instead of the
normal law the chi-square layer is derived under, so the assumption is
measured against its violation instead of being stated. A cohort of ten
lineages therefore has a wide second interval, and that is a property of the
question, not a defect of the estimator.

The permutation p-value beside the share cannot fall below
`1 / (n_perm + 1)`, and the record carries that floor as `p_floor`. A p-value
at the floor says that no permutation reached the observed share; it is not
an estimate of how far below the floor the p-value lies, and the report
prints it as such. The control to read is the permuted-label share and its
interval, which the number of permutations does not bound.

For a single binary trait the record also restates the share and both of its
intervals on the latent scale (`latent_share`, `latent_low`, `latent_high`,
`latent_species_low`, `latent_species_high`), the scale a mixed model with a
binomial link reports its intraclass correlation on. Under a threshold model
the observed-scale share `rho_obs` at prevalence `p` is
`(Phi_2(-t, -t; rho_L) - p^2) / (p (1 - p))` with `t = Phi^{-1}(1 - p)` and
`Phi_2` the bivariate normal distribution function; the map from the latent
share `rho_L` is increasing, so `latent.latent_share` inverts it by
bisection and the endpoints of an interval carry across with their coverage
unchanged (`latent.py`). This is the liability-scale transformation of
Dempster and Lerner (1950), computed exactly and not to first order. The
observed scale remains the one every interval is measured on and the one a
prevalence table uses; the latent figure is larger, because a call discards
the part of the liability that does not cross the threshold, and is offered
so that the two readings of one cohort can be set side by side.

`realised.realised_share` reports the other question: how much of the trait
variance sits between the lineages this collection actually holds. The lineage
effects are then fixed unknowns rather than a fresh draw, and the parameter is
`S_a^2 / (S_a^2 + sigma^2)` with `S_a^2 = sum_g n_g (a_g - abar_w)^2 / (n0 (G - 1))`,
`abar_w` the isolate-weighted mean of the effects and `n0` the lineage size a
balanced design of the same total would need; on a balanced design this is
`sum_g (a_g - abar)^2 / (G - 1)`.

Conditional on those effects, `SSB / sigma^2` is a noncentral chi-square on
`G - 1` degrees of freedom with noncentrality
`lambda = sum_g n_g (a_g - abar_w)^2 / sigma^2`, independent of the central
`SSW / sigma^2` on `n - G`. The ratio of mean squares is therefore a noncentral
F, the set of noncentralities it does not reject is a confidence interval for
`lambda`, and the share is a monotone function of it,
`lambda / (lambda + n0 (G - 1))`, so the endpoints carry across. The argument
does not use balance, which is why the interval keeps its level under the
unequal lineage sizes real cohorts have while the variance-component one does
not. Obtaining a confidence interval for a noncentrality parameter by
inverting the noncentral distribution in it is Venables',
[10.1111/j.2517-6161.1975.tb01554.x](https://doi.org/10.1111/j.2517-6161.1975.tb01554.x),
and the same inversion carries an exact interval to a monotone function of
that parameter, which is how exact intervals for fixed-effect analysis of
variance effect sizes are built,
[10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164).
Treating the group effects present in a collection as a finite population
with its own variance, rather than as a draw from an infinite one, is
Cornfield and Tukey's,
[10.1214/aoms/1177728067](https://doi.org/10.1214/aoms/1177728067). What is
assembled here is the pairing: the two estimands reported side by side from
one fit of a bacterial cohort, with a gate on the assumption the exact one
needs.

Both intervals come from one fit and are reported together. Which to read is a
question about the target of inference, not about precision: a laboratory
asking what its own collection shows wants the realised one, and a reader
generalising to the species wants the wider one.

They also answer the same yes-or-no question identically. Both invert the same
ratio of mean squares, so a lower endpoint stands above zero exactly when the
same central F test rejects, and the narrower interval is therefore not a more
sensitive test and is not offered as one. Sensitivity and specificity of that
decision are measured in section S12.8 of the article's supplement; the identity
is held by a contract test, not by this paragraph. What the realised
interval buys is width, and validity under the unequal lineage sizes real
cohorts have.

### 1.1 What was measured, and what was refused

The operating characteristics are in section S12 of the article supplement and
are produced by `benchmarks/realised_calibration.py`. Two conditions were fixed
before that campaign: an anchor arm that re-derives the reference law from the
model statement instead of the implementation had to reproduce the calibration
arm, and a cross-scoring arm had to break coverage when each interval was
scored against the other estimand's truth. Both held.

Coverage of the realised interval degrades as the within-lineage residuals grow
heavy, because it is the construction that extracts more from their shape. The
`KURTOSIS_LIMIT` in `realised.py` is the largest median excess kurtosis at
which coverage stayed at or above 0.93 on the measured curve, and above it the
result is marked not estimable with the measured kurtosis in the reason. The
kurtosis is read from residuals of lineages with at least two members, each
divided by `sqrt(1 - 1/n_g)`: a singleton contributes a residual of exactly
zero, and a gate reading those would be measuring the lineage size
distribution instead of the trait.

Two repairs were tried and are recorded in the campaign instead of being
deleted. A
Box-type deflation of the degrees of freedom by the estimated residual kurtosis
is no better than the exact interval on Gaussian residuals and destroys it on
heavy ones. A within-lineage percentile bootstrap, intended as an
assumption-free fallback, under-covers under the very process the method is
derived against. Neither is used.

---

## 2. Lineage-aware surveillance reading

Section 1 describes the share of a level. This one describes a change, and
what a prevalence table cannot separate.

### 2.1 Decomposing a difference

For two collections A and B, with lineage shares `w` and within-lineage rates
`p`, the Kitagawa (1955, JASA 50:1168-1194) identity splits the observed
difference exactly:

    P_A - P_B = sum_l (w_A - w_B)(p_A + p_B)/2      composition
              + sum_l (w_A + w_B)/2 (p_A - p_B)     within-lineage

This is the demographic ancestor of the Blinder-Oaxaca decomposition (1973) and
of Fairlie's binary-outcome extension (2005). The non-parametric form is used
instead of a fitted one because a sequence-type variable carries of the order
of a hundred levels on a few hundred isolates with a long singleton tail: a
logit on ST dummies is separated before it is asked anything.
`identity_residual` is emitted so exactness can be checked, and
`benchmarks/decomposition_vs_regression.py` measures the claim instead of
resting on it. Over 14,400 simulated cohorts the decomposition and a
penalised fixed-effect logistic agree when lineages are shared (root mean
squared error 5.26 against 5.31 percentage points) and separate under 40 %
lineage turnover (4.71 against 6.53). With 100 lineages and 100 isolates per
collection, a mean of 34 lineages have all outcomes equal, so the unpenalised
fit does not exist and the penalty that produces one is a choice that moves the
estimate.

A lineage seen in only one collection has no rate in the other, and the
identity holds for any value substituted, so the substitution is a definition.
Setting the absent rate equal to the observed one makes that lineage's
within-lineage term exactly zero, and therefore makes the within-lineage
component a function of shared lineages alone, independent of the convention.
Lineage turnover is charged to composition, which is the demographic reading.

### 2.2 Two estimability gates

The convention removes the dependence on an arbitrary choice but not its cost,
so the cost is gated.

**Shared support.** `shared_support_isolate_share` is the fraction of isolates
in lineages common to both collections. Below `min_shared_support` the
within-lineage component is returned with `within_lineage_estimable` false and
a signed margin. The default of 0.8 comes from
`benchmarks/decomposition_calibration.py`: over 720 scenario cells and 720,000
simulated decompositions, coverage of the nominal 95 % interval for the
within-lineage component runs 0.82 below a support of 0.4, 0.91 between 0.5 and
0.6, and 0.94 above 0.8, while the composition component keeps nominal coverage
throughout and is not gated. Above the gate the mean absolute bias is 0.26
percentage points and the worst cell 2.3; below it the bias reaches 7.8 points,
which is what the gate refuses. Under an exact null the type I error is 0.049.

**Label availability.** A lineage-resolved statistic is computed on the
isolates that carry a label and describes the collections only if labelling was
unrelated to the trait. Two conditions must hold together. If coverage differs
between the collections but labelling is unrelated to the trait, the same
selection applies to both arms and the difference is largely unaffected. If the
trait differs between labelled and unlabelled isolates but coverage is equal,
missingness costs precision rather than validity. Together they mean the
difference is computed between two differently selected subsets. Both are
tested by Fisher exact tests and reported in
`lineage_label_availability`; where both fail, `composition_estimable` and
`within_lineage_estimable` are both false, and neither component can be a
discovery, since the labelled isolates are then not the collections. A Fisher
exact test needs a binary trait, so on a trait that is not binary the second
condition cannot be tested at all, and a cohort holding any unlabelled isolate
is refused there rather than cleared on an untested gate.

On the shipped *S. suis* cohort the sequence type is missing for 219 of 677
isolates, 40 % of the later United Kingdom period carries one against 97 % of
the earlier, and untyped isolates carry about two more non-wild-type results
out of thirteen. Ceftiofur then appears to fall by 11.0 points on the labelled
subset while the collection rises by 9.5. The gate fires and names the reason.

### 2.3 The panel is a family

An antimicrobial panel is a family of tests whose members are not independent.
`decompose_panel` controls the false discovery rate within each component
family by the Benjamini-Yekutieli step-up, which is valid whatever the
dependence between agents, separately for composition and for the
within-lineage rate because they answer different questions and are not
exchangeable, and reports the effective number of independent agents as the
participation ratio of the panel correlation eigenspectrum. On the shipped
thirteen-agent panel that number is 4.8, which is what the tetracycline,
macrolide-lincosamide, beta-lactam and remaining blocks amount to.

### 2.4 Two prevalence estimands

Prevalence per isolate answers clinical burden; prevalence per lineage answers
diversity. Neither is the true one, so a report should say which it means.
The interval on the per-lineage figure is a two-stage
cluster bootstrap: lineages are resampled first, because the estimand averages
over lineages and its uncertainty is dominated by how many distinct ones were
seen, and isolates are then resampled inside each drawn lineage so that a rate
measured on two isolates is not treated as exact. Because a lineage of size `m`
with rate `r` resampled to size `m` yields `Binomial(m, r) / m`, the second
stage is drawn exactly rather than by index.

### 2.5 Concentration of carriage

The reciprocal of a Herfindahl-Hirschman index gives the effective number of
lineages carrying a trait. Concentration alone is not evidence, since an uneven
cohort produces an uneven carrier pool with no biology involved, so the tested
quantity is the departure of the carrier distribution from the
lineage-abundance distribution in bits, against a permutation null that holds
the number of carriers fixed. **The departure is a magnitude**: it fires when
carriage piles into one clone and when carriage avoids the dominant one.
`direction` is reported with it, and on the shipped *S. suis* panel it is the
second case that fires. The permutation floor, the exceedance count and an
exact Clopper-Pearson interval for the tail probability are reported with the
p-value, so that a panel run at a small budget cannot present `1 / (n_perm + 1)`
thirteen times as thirteen strong results.

## 3. Reading a susceptibility panel at its recorded resolution

### 3.1 One likelihood, three interval widths

A non-wild-type call, a recorded minimum inhibitory concentration and a
censored reading are one likelihood at three interval widths on the log2
concentration scale. A call is `(-inf, c]` or `(c, inf)`; a dilution recorded
at well *w* is `(previous tested well, w]`; a reading on the lowest or highest
tested well is unbounded on that side. `censored.intervals_from_binary` and
`censored.intervals_from_mic` build the intervals and
`censored.censored_clonal_share` consumes nothing else, so the input mode
changes how much an observation carries and not what is estimated.

Which wells were tested is a property of the panel, not of the readings. A
laboratory that records its panel passes it as `wells`, and every reading is
placed on the nearest tested well. Without it, readings that all sit within a
quarter of a doubling of an integer power of two are taken to come from a
doubling panel whose tested wells are every doubling from the lowest to the
highest recorded value: a dilution no isolate landed on was still tested, so
it bounds its neighbours, and 0.06 and 0.064 are one well recorded two ways
rather than two wells a twentieth of a doubling apart. Readings off any such
lattice are taken as their own wells, and `panel_geometry` reports which
reading was made (`lattice`), how many wells the panel is taken to have and
how many distinct values were recorded. On the shipped cohort this reading
folds rounding variants on seven agents and fills unoccupied wells on six;
it moves the penicillin share read from the dilution by 0.024 and no other
agent's by more than 0.01.

The conditional moments of a normal restricted to an interval are computed by
reflecting an interval that lies above the mean, because `log_ndtr` is
accurate for large negative arguments and saturates at zero for large
positive ones. Before the reflection an interval more than about seven
residual standard deviations above the lineage mean returned the mean itself
as its conditional expectation, so an isolate far beyond the panel top was
read as sitting at its lineage mean; on the shipped cohort that moved the
amoxicillin share read from the dilution from 0.758 to 0.713 and no other
agent by more than 0.005. `tests/test_censored.py` checks both tails against
an independent implementation.

Two tests in `tests/test_censored.py` hold that claim to account. With
zero-width intervals the estimate must equal the classical one-way
variance-component ratio computed in closed form; and a concentration
dichotomised at a cut-off must produce the intervals the call constructor
produces, and therefore the same share.

### 3.2 The estimand, and how it differs from the clonal share

The quantity returned is `tau2 / (tau2 + sigma2)`, the share of the latent
log2 concentration variance carried by lineage. This is what a mixed model
reports and what the heritability literature compares against. It is **not**
the same quantity as `attribution.clonal_share`, which reports achievable
out-of-sample predictive skill and falls when lineages are small because a
group mean estimated from three isolates predicts badly however much variance
the grouping truly carries. A cohort should carry both; their ratio is a
property of the study design, not of the biology.

### 3.3 The maximisation step

Fitting is expectation maximisation on the interval likelihood with a
restricted-likelihood maximisation step. Three choices in that step were each
adopted because it moved a measured bias, and all three are checked in
`benchmarks/censored_calibration.py`.

* **The posterior variance of a lineage mean stays in the between-lineage
  variance.** A mean estimated from a handful of censored readings is
  uncertain, and a moment estimator that keeps only squared deviations
  discards that uncertainty.
* **The noise in a lineage average of conditional expectations is
  `sigma2` less the mean conditional variance**, by the law of total variance,
  and not `sigma2`. An interval has already resolved part of an observation, so
  a censored reading contributes less noise to the average than an exact one;
  subtracting the whole of `sigma2` over-subtracts, and the error grows with
  the censored fraction.
* **The divisors are restricted**: populated lineages less one for the
  between-lineage variance, `n` less the number of populated lineages for the
  residual scale. On a balanced exact design this reproduces the classical
  estimator; on an unbalanced one it beats it.

### 3.4 Identification, and two intervals

A single cut point recovers only the standardised distance from the cut, so
the residual scale is not identified. `censored.scale_is_identified` detects
that and the fit holds the residual standard deviation at one, which puts the
result on the liability scale of Dempster and Lerner (1950) and makes it
comparable with published heritabilities instead of with a raw log2 variance.

Three widths are reported and only one of them is an interval.

**The interval** inverts the variance ratio. For a one-way random model the
ratio of the between-lineage to the within-lineage mean square is a scaled F,
and inverting it gives an interval for the intraclass correlation. The
expectation-maximisation fit returns the same two components, so the same
inversion applies to interval data once the within-lineage degrees of freedom
are discounted by the share of each observation the interval left unresolved.
For a recorded value or a dilution interval its coverage is 0.93 to 0.98 by
true share and 0.85 to 1.00 cell by cell over 40 cells of 200 simulated
cohorts, the lowest at a hundred lineages of six isolates, where the
within-lineage degrees of freedom are fewest. Read from a single cut point,
the same interval holds its level only where the cut sits near the middle of
the distribution and loses it towards either tail (the 120-cell grid of the
supplement reads the same cohorts six ways); a binary call is therefore read
by `clonal_share`, and this function is for the recorded dilution.

**The cluster bootstrap** draws whole lineages, with bias-corrected and
accelerated limits. It is reported beside the interval and is not the interval:
on the first run of the same grid, with 200 draws a cohort, it covered a
median 0.61 against a nominal 0.95. Resampling
thirty lineages does not reproduce the sampling distribution of a variance
component built from those thirty lineages, and the bias correction does not
repair that. It is kept because it assumes no distribution and therefore fails
differently from the F inversion, which is what a check is for.

**The likelihood-ratio width**, from the Gauss-Hermite marginal likelihood in
`censored.marginal_loglik`, is not a confidence statement at all. It measures
how sharply this cohort pins the share down, which is the quantity that
improves when a dilution is read instead of a call, and its coverage is
reported so that no reader mistakes it: 0.68 on an exact reading and 0.74 on a
dilution at thirty lineages, 0.91 to 0.93 on a single cut point. The
chi-squared calibration is optimistic exactly where the residual scale is well
determined and the whole of the uncertainty sits in the group means.

### 3.5 The coarsening assumption, and the refusal

Treating a reading on an end well as censored is an assumption about why the
reading is there, not an observation: it is the coarsened-at-random condition
of Heitjan and Rubin (1991), under which the coarsening is ignorable for
likelihood inference. Three things follow. `censored.panel_geometry` reports
the mass on each end well and refuses the point mode above five per cent. A
recorded operator column overrides the assumption wherever the laboratory
supplied one, because an operator is an observation and the heuristic is not.
And `censored.sensitivity_endpoints` computes the share under both readings,
so the width of that bracket can be set beside the bootstrap interval.

Two refusals follow from the same grid, not from judgement.

A lineage every one of whose readings is one-sided has no identified latent
mean. Counting such lineages is the wrong gate, since on a susceptible agent a
few are normal; what moves the answer is how much of the cohort sits there.
`censored.CENSORED_GROUP_LIMIT` is set at one half, read off the censoring
sweep.

A single cut point in a tail leaves almost no contrast to divide. At a
non-wild-type prevalence near 0.08 the share reads 0.23 above a true zero and
its interval covers 0.08 of the time; `censored.SINGLE_CUT_PREVALENCE`
therefore refuses a call outside 0.10 to 0.90. Inside the window the estimate
is defined and its bias stays under 0.05, but the interval is calibrated only
from a prevalence of about 0.24: one doubling from the median it covers a
true zero 0.58 of the time (supplement S10.5), and the record's `notes` say
so. A single-cut reading is therefore a reading of the call for comparison
with the dilution, which is how the article uses it, and the pipeline reads a
call by `clonal_share`. A dilution at the same prevalence is unaffected,
which is the practical argument for reading the panel: on the shipped cohort
two agents are refused as calls and estimable as dilutions.

## 4. Anytime-valid evidence for a lineage effect

Surveillance re-reads one panel whenever a year of isolates arrives. A
false-discovery procedure recomputed at an unplanned number of looks controls
nothing, because the number of looks is not fixed in advance and the looks are
not independent.

`evalues.e_process` returns an e-value in the betting sense of Vovk and Wang
(2021). This is neither the BLAST expectation value nor the E-value of
VanderWeele and Ding for sensitivity to unmeasured confounding; both are common
in this literature and neither is meant. The construction is the split
likelihood ratio of Wasserman, Ramdas and Balakrishnan (2020): lineage
probabilities are fitted on a training fold, shrunk towards the training grand
mean by a one-way empirical-Bayes factor estimated on that fold, and scored on
the held-out fold against the null probability maximised on that same fold.
Because the numerator parameters never saw the held-out rows, the expectation
under the null is at most one with no regularity conditions. By Ville's
inequality an e-value of at least `1/alpha` is the counterpart of a p-value
below `alpha`, and a product of such e-values over independent intakes may be
inspected as often as wanted (section 4.1); the e-value of one look, on its
own, is a statement about that look.

Values from independent batches multiply, forming a test martingale, which is
what makes stopping at any point legitimate; folds of one cohort share data and
are averaged instead, which is valid because any convex combination of e-values
is an e-value. Across a panel, `evalues.e_bh` applies the procedure of Wang and
Ramdas (2022), which controls the false discovery rate under arbitrary
dependence between the hypotheses. That is the property the panel needs:
cross-resistance makes a macrolide block behave as one trait, so the agents are
neither independent nor reliably positively dependent, and the usual
justification for Benjamini-Hochberg does not apply.

### 4.1 One look, and a programme of looks

`e_process` is the e-value of one look. Recomputing it on the accumulated
cohort after every intake gives a valid e-value at each look, and e-BH on
those values controls the false discovery rate at each look; but the sequence
of recomputed values is not a test martingale, because the training and
held-out rows of one look are the training rows of the next, so Ville's
inequality does not bound the chance that the running value ever crosses
`1/alpha` over the programme. The benchmark of `benchmarks/repeated_looks.py`
measures that rule (`ebh_accumulated`) and records what it did on 144,000
simulated panels; what it did is a measurement, not a guarantee.

The construction that carries a guarantee over the programme is
`evalues.sequential_e_process`. Batches are the intakes in arrival order.
For batch `t` the lineage probabilities are fitted, with the same shrinkage,
on batches `1` to `t - 1` pooled, and the batch is scored against the null
probability maximised on the batch alone:

    E_t = prod_i q_{t-1}(y_ti) / sup_p prod_i f_p(y_ti).

Given the past, the numerator is a fixed probability law on the batch and the
denominator dominates every null law, so `E[E_t | past] <= 1` for every null
`p` and the product over batches is a test supermartingale with initial value
one. It may be inspected after any batch and stopped at any time, and e-BH on
the running products of a panel controls the false discovery rate at any
stopping time a programme chooses. It is the fixed-look split likelihood ratio
with the training fold replaced by the past, which is why its power grows
with the programme rather than with the batch: on the same simulated panels it
declares the true agents at least as often as the recomputed rule and more
often than the product of per-intake e-values, which trains each factor on
one intake alone. The first batch contributes nothing, because nothing
precedes it to train on. The batches must be disjoint isolates in the order
they were collected; the function cannot tell a re-used isolate or a reordered
batch from a new one, and the guarantee rests on the caller keeping them apart.
Isolates with no lineage label or no finite reading are set aside within their
batch and counted in `n_dropped_untyped` and `n_dropped_non_finite`, so the
anytime validity is conditional on that missingness being ignorable: where a
well is left untested, or a typing left unattempted, because of what the result
was expected to be, the batches scored are not the batches collected and no
supermartingale property is claimed for them.

A run reaches this construction without leaving the command line by declaring
`dataset.batch_column`, the metadata column naming the intake an isolate
arrived in. The batches are its levels in ascending order, which is why the
column has to sort into arrival order; an isolate with no intake recorded is
set aside and counted, not assigned to one. The record then carries
`lineage_evidence.sequential` beside the single-look block: the batches and
their sizes, the running product per trait, and a second e-BH decision. That
decision is taken by `evalues.e_bh_log` on the logarithms, because a product
over many intakes overflows a double long before its logarithm does, and an
e-value of positive infinity is the strongest evidence there is, not a
missing value; `e_bh` ranks an infinite e-value first and delegates
to the log-scale routine.

## 5. The estimator on other species

`benchmarks/atlas_cross_species.py` runs `attribution.clonal_share` on NCBI
Pathogen Detection, which supplies both variables from one public release: the
lineage is the SNP cluster (`PDS_acc`) and the trait is non-susceptibility,
R or I against S, read from the `AST_phenotypes` metadata field. Nothing in it
depends on a cut-off derived in this package, so what it exercises is the
estimator and not the choices made for the shipped cohort.

A species enters if 60 isolates join between the metadata and the cluster file,
and an agent enters if 50 of them carry a call and at least 2 per cent fall in
the minority class. Every cell is run twice, once as recorded and once with the
lineage labels permuted within the analysed subset. Permuting the labels of the
whole cohort and subsetting afterwards would draw from a larger pool and change
the lineage size distribution, so the control would not be comparable.

The condition set before the run was that the permuted arm must be
indistinguishable from zero; the measured values are in section S11 of the
article supplement, with the release accession, isolate count and row count
recorded per species in the evidence receipt beside them.

Two limits are part of the design. The phenotypes are contributed by many
laboratories, under different standards and over different years, so a share
measures how the estimator behaves on heterogeneous data rather than the
epidemiology of the species. And NCBI recomputes the SNP clusters at every
release, so a rerun against a later release reads a different lineage variable;
that is why the release accession is recorded and not only the date.

The same source is read a second way for the veterinary question
(`benchmarks/vet_atlas.py`): the isolates of one organism are first cut by
host and by where in the chain the sample was taken
(`benchmarks/vet_source_taxonomy.py`), because a laboratory holds a species
in a host rather than a species, and every cell is run at two definitions of
a lineage, the SNP cluster and the serovar (`benchmarks/resolution_atlas.py`),
with the permuted-label control at both. The pair is reported rather than one
number, because the share a cohort carries is a function of the resolution at
which a lineage is defined. `benchmarks/period_split.py` splits one cell at a
year fixed before the run and reports the share on each side of the cut,
which is how a rise in prevalence is placed in one lineage of the collection
or across the lineages it already held; a frame with fewer than twenty
non-susceptible isolates is reported as void and not estimated.

## 6. Reproducibility

Every stochastic stage is seeded from `--seed`: the share, the e-values and
the surveillance readings draw from one generator spawned off it via
`numpy.random.SeedSequence`, in agent order, and the dilution reading takes
the seed itself, so two runs of one configuration on one stack write one
record. Every agent's censored bootstrap is therefore given the same integer
seed, so the per-agent interval widths are common random numbers rather than
independent draws, and two agents on one cohort move together more than two
independent runs would.
`clonal_share_result.json` is serialised with `allow_nan=False` and an explicit
numpy encoder, so it is valid RFC 8259 JSON that a non-Python parser can read.

## 7. What the package assumes, and where its guarantees stop

Every estimator above rests on conditions the collection may not meet, and
the gates catch some of them but not all. The ones a reader should hold in
mind are these.

Isolates are treated as exchangeable within a lineage. Farm, herd, outbreak
or repeat sampling nested inside a lineage inflates the between-lineage
component and narrows every interval built by resampling isolates: the
lineage bootstrap, the within-lineage bootstrap of the species interval, the
two-stage prevalence bootstrap and the Kitagawa bootstrap. The package has no
column for a sampling unit below the lineage and does not correct for one; a
collection with strong within-lineage clustering should be thinned to one
isolate per unit, or its shares read as upper bounds.

Results are taken as missing at random within an agent. An isolate without a
result for an agent leaves that agent's estimate, and the record counts it
(`n_dropped_non_finite`); testing that depends on lineage or on the year of
intake biases the share and the sequential e-process in a direction the
package cannot tell from the data.

The species interval adds the sampling of lineages to a within-lineage
resample taken at the observed share. It is a plug-in approximation whose
level was measured on the grid, not derived, and under a two-point
carrier law with few lineages it falls below its level (supplement S14);
the record flags a dominant lineage so that the reader knows when that law is
the likelier one.

The permutation debiasing of the clonal share is exact under a linear
relation between the raw out-of-sample score and the share; the residual bias
at a thousand lineages of three isolates, −0.017, is the measured departure.

The dilution reading assumes the coarsening is ignorable: a value on an end
well is treated as censored at that well unless the laboratory recorded an
operator. Read from a single cut point instead of a dilution, its interval
holds its level only where the cut sits near the middle of the distribution,
which is why a binary call is read by `clonal_share` and not by this
function. The latent-scale restatement of a binary share assumes a Gaussian
liability with a fixed prevalence.

The exact interval of the realised share is available only where the
within-lineage residuals pass the kurtosis gate, which on a binary call
closes wherever the prevalence sits far from one half; on the shipped
cohorts it opens for 5 of 13 and 4 to 6 of 22 agents. Its refusal is a
statement about the interval, not about the estimate, which is printed
either way.

The sequential e-process holds its guarantee at any intake only if the
intakes are disjoint sets of isolates in the order they were collected;
within an intake the order is ignored, and a year is a coarse intake. The
e-BH decision on the running products is valid at a stopping time common to
the panel, not at a different intake for each agent.

A lineage is one flat label. Hierarchical structure, such as sequence types
nested in clusters, and phylogenetic distance between lineages are not
modelled; reading a cohort at two resolutions, as the article does, is the
package's way of showing what the choice of label is worth. The concentration
of carriage reports a direction and a permutation p-value that come from two
different statistics and can disagree; the report reads the direction only
where the p-value resolves.

Ten of the thirteen cut-offs of the shipped *S. suis* example were derived
from the cohort itself, where no published cut-off exists; they are
collection-specific, and the supplement of the article measures what a
one-dilution shift of each is worth.
