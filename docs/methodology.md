# Methodology

Every statistic the package computes, its definition, the choice it embodies and
the citation behind it. Where a default is consequential the reasoning is given
here rather than in a docstring.

---

## 1. Inputs and gating

**Layers.** One or more wide binary matrices, strains × features, declared in
the config with a `kind`:

* `wide_binary`: independent 0/1 presence calls. Jaccard-family distance.
* `one_hot`: a categorical variable expanded to indicators (capsule K locus, O
  antigen, MLST). Validated at load time to have at most one positive per row,
  given a matching coefficient rather than Jaccard, excluded from count
  aggregation, and counted as **one** hypothesis in the FDR set rather than one
  per level. Jaccard over a one-hot block takes only a handful of distinct
  values and is not a usable similarity.

**Feature groups.** `files.<layer>.groups` names co-inherited blocks (operons,
mobile elements). With `collapse_feature_groups: true` each group becomes a
single presence call *before* distances are computed. This matters twice: an
uncollapsed encoding weights a biological unit by its gene count, and it also
breaks the post-clustering inference (§4.2).

**Prevalence gate.** Keep features with `lo < prevalence < hi` and at least
`min_count` carriers. The default `lo = 0.0, min_count = 2` deliberately keeps
rare determinants: a 2 % floor removes exactly the emergent colistin,
tigecycline and carbapenemase columns that AMR surveillance exists to detect.
`files.<layer>.protected` exempts named features from the gate. A feature
present in **zero** isolates is always removed, because it cannot inform a
distance, cannot be tested, and would occupy a slot in the FDR denominator, and
is
reported with that reason. Everything removed appears in
`layer_reports[].gated_out`.

---

## 2. Distances

| metric | definition |
|---|---|
| `jaccard` | `1 - a/(a+b+c)` |
| `dice` | `1 - 2a/(2a+b+c)` |
| `simple_matching` / `hamming` | `(b+c)/p`, counts shared absence |

with `a` the shared presences and `b`, `c` the exclusive counts.

**The empty-union convention is a declared choice, not a default.**
Jaccard is undefined when two isolates carry none of a layer's features. The
textbook convention sets the similarity to 1, on the reading that two genomes
sharing no detected
determinant are declared maximally similar. On accessory-genome data this is
consequential: absence is the default state, it is confounded with assembly
quality, and trait-absent isolates are usually the majority. In the *Klebsiella*
cohort 1004 of 1500 isolates had an all-zero virulence row; under the classical
convention they formed one clique of 503 506 pairs at distance zero, and the
all-zero indicator alone reproduced the resulting partition at ARI 0.83.

`distance.undefined_pair` takes `identical` (D = 0, classical), `distinct`
(D = 1) or `nan`. Whichever is chosen is echoed into the result, and
`artifact_diagnostics.per_layer[].empty_stratum_ari` reports the ARI between the
partition and each layer's all-zero indicator. Above 0.5 the run is flagged.

---

## 3. Fusion

Similarity network fusion, Wang et al. (2014), *Nat Methods* 11:333-337,
[10.1038/nmeth.2810](https://doi.org/10.1038/nmeth.2810).

**Kernel.** `W(i,j) = exp(-d²/(mu·eps_ij))` with

```
eps_ij = ( mean_{k in N_i} d(i,k) + mean_{k in N_j} d(j,k) + d(i,j) ) / 3
```

where `N_i` is `i`'s **K nearest neighbours**, self excluded. Using the global
mean distance instead, as this package did before v2.0, removes the local
scaling the kernel exists to provide. Both reference implementations (SNFtool's
`affinityMatrix.R`, snfpy's `make_affinity`) use the K-nearest-neighbour mean.
Note that those implementations use a Gaussian *density* rather than the
exponential form above; the two differ by a scale factor absorbed into `mu`.

**Cross-diffusion.** `P_v <- normalise( S_v · mean(P_{u != v}) · S_vᵀ )` for `T`
iterations, with each view's transition matrix re-symmetrised **inside** every
iteration (SNFtool does this; symmetrising only the final average lets asymmetry
accumulate over `T` steps).

**Sparsification and ties.** `S_v` keeps each row's `K` largest affinities.
`argsort(-W[i])[:K]` is ill-defined when more than `K` candidates tie at the
cut, and binary trait data produces exactly that: on the *Klebsiella* virulence
layer the median row had ~1000 candidates tied at `K = 30`. Which of them was
kept depended on the input file's row order, so the partition was not
equivariant under relabelling of samples (ARI 0.59-0.98 against the unpermuted
run). `tie_policy: inclusive` keeps **every** neighbour tied with the `K`-th
largest, which is order-invariant by construction. The test suite asserts
ARI = 1.0000 across row permutations, and the run reports
`layer_reports[].knn_ties.tie_inflation` = mean degree / `K`. A value far above
1 means the affinity has too few distinct values for a KNN graph to be
meaningful.

`snf.K`, `snf.mu`, `snf.T`, `snf.alpha` and `snf.tie_policy` are all config
keys. A single-layer fusion is legal and returns that layer's normalised
transition matrix.

---

## 4. Post-clustering inference

### 4.1 The problem, and why data thinning does not solve it here

Testing which features distinguish clusters estimated from those same features
is double dipping: Gao, Bien & Witten (2024), *JASA* 119:332-342,
[10.1080/01621459.2022.2116331](https://doi.org/10.1080/01621459.2022.2116331);
Chen & Witten (2023), *JMLR* 24(152).

Data thinning, Neufeld, Dharamshi, Gao and Witten (2024), *JMLR* 25(57):1-35,
splits an observation into independent parts from the same family. Its Remark 8
states that infinite divisibility "prevents us from thinning the Bernoulli or
categorical distributions". Coarse-graining to counts and thinning those fails
twice on real panels: the aggregates are often rescaled indicators rather than
counts (a yersiniabactin tally taking only {0, 8}; one-hot capsule "counts"
in {0,1} with an exact linear dependence among them), and the negative-binomial
recipe needs the size parameter **known**. Proposition 11 of the same paper
gives the damage when it is not:

```
cov(X1, X2) = (1 - eps) · r · ((1-p)/p)² · (1 - (r+1)/(r~+1))
```

negative when `r~ < r`, positive when `r~ > r`. With a pooled estimate
`r̂ = 0.43` the realised within-column correlation on the *Klebsiella* aggregates
was −0.174, against +0.003 for a genuine NB matrix with `r` known. Hivert et al.
(2024, [arXiv:2405.13591](https://arxiv.org/abs/2405.13591)) reach the same
conclusion empirically.

Where a count model *is* defensible the machinery is retained and **checked**:
per-feature dispersion, `screen_thinnable` refusing under-dispersed,
near-binary and exactly collinear columns with a stated reason, and the realised
`corr(X1, X2)` reported and gated. `binomial_thin` implements the hypergeometric
split for bounded row-sums over `m` known members, where nothing is estimated.

### 4.2 Feature splitting, and the block requirement

For binary panels the split is over **features**: cluster on a random half,
test the held-out half, repeat. The group contrast within a split is the largest
*discovery-block* cluster against the rest, so no grouping or direction is ever
chosen from the values being tested.

**The split must respect the dependence structure of the features.** The null
actually tested is

> H0: the held-out features are independent of the labels learned from the
> discovery features

which is implied by "features are mutually independent", **not** by "there are
no clusters". Measured on a no-cluster null with features in correlated blocks
of five, a per-feature split has family-wise error **1.000** at nominal 0.05,
because a block straddling the two halves makes the discovery labels predict the
held-out features for reasons unrelated to clustering. Passing `groups`, which
`core.run` does automatically from the config's declared feature groups, with
each `one_hot` layer as a single unit, assigns whole blocks to one side and
restores validity: **0.042** at within-block correlation 0.5, **0.017** at 0.9.

Each split re-runs the *same* procedure the reported partition came from
(per-layer distances with each layer's declared kind, cross-diffusion fusion,
spectral clustering) on the discovery features, and
`ari_discovery_vs_reported` records whether the per-split labels resemble the
reported partition. If they do not, the p-value is about detectable structure in
general and not about this partition, and `concerns_reported_partition` says so.

### 4.3 Merging the splits

One split gives one p-value, and that p-value is random even with the data
fixed: forty draws of a single-split implementation returned p between
3.8 × 10⁻⁶³ and 0.97. Repeated splits of one fixed dataset are **exchangeable**
but arbitrarily dependent, so Fisher and Stouffer combination are invalid and
the minimum is catastrophically so. A multiplicative constant is unavoidable.

| `tva.merge` | rule | reference |
|---|---|---|
| `exchangeable_ruger` (default) | `(K/k)·min_{l≤K} p^(l)_(⌈lk/K⌉)` | Gasparin, Wang & Ramdas (2025), *PNAS* 122:e2410849122, Thm 4.1 |
| `ruger` | `(K/k)·p_(k)` | Rüger; see Vovk & Wang (2020) §2 |
| `twice_mean` | `2·mean(p)` | Vovk & Wang (2020), *Biometrika* 107:791-808 |
| `mmb` | `min{1, (1−log γ_min)·inf_γ q_γ(p/γ)}` | Meinshausen, Meier & Bühlmann (2009), *JASA* 104:1671-1681, Eq. 2.3 |

All four hold their level on an equicorrelated Gaussian-copula exchangeable null
at ρ ∈ {0, 0.5, 0.9}; the exchangeable Rüger rule is the least conservative,
which is why it is the default.

### 4.4 Structure is not the same as discreteness

The split test rejects on a single continuous gradient as readily as on discrete
groups. Comparing a `k`-component Bernoulli mixture with one component by BIC
does not settle it either: two components fit a gradient better than one.

`inference.continuum_null_test` fits a *q*-factor logistic latent-trait model.
It evaluates `q = 1..3` first and, only when `q = 3` is the provisional boundary
optimum, extends to the supported `q = 4`; the default maximum is therefore 4.
It then simulates from the selected model and compares the
observed `BIC(q-factor) − BIC(k-component mixture)` against what the continuum
itself produces. The null therefore preserves the feature dependence that
breaks the naive tests. On planted discrete data it rejects; on data simulated
from a sufficiently dimensioned continuum it does not. If BIC is still falling
at `q_max`, the function stops before bootstrap and returns
`withheld_under_dimensioned`, `p_value = null` and a null discreteness verdict:
rejection against a known-under-dimensioned continuum is not evidence of
archetypes. With an estimable interior optimum the verdict is explicitly
`discrete` or `continuum_compatible`. The p-value is Phipson-Smyth corrected, so
its floor is `1/(n_boot+1)`, reported alongside; at least 20 bootstrap draws are
required so a test at alpha 0.05 is attainable.

### 4.5 Small cohorts: exact phylogenetic convergence

Splitting throws away information and has no power at *n* of order 100.
`amr_clonalshare.archephy` asks a different question with an exact answer:
did a candidate archetype's `k` traits arise *repeatedly and independently*
across a phylogeny?

Fitch parsimony (Fitch 1971, *Syst Zool* 20:406-416) gives per trait the set
`S_t` of change edges, `c_t = |S_t|`, on a tree with `E` non-root edges. The
`k`-way joint homoplasy count is `m = |S_1 ∩ … ∩ S_k|`. Under independent
evolution with exchangeable edges the binomial moments are
`E[C(M,j)] = C(E,j)·∏_t C(c_t,j)/C(E,j)`, so

```
P(M >= m) = sum_{j>=m} (-1)^(j-m) · C(j-1, m-1) · E[C(M,j)]
```

evaluated in exact rational arithmetic. At `k = 2` this equals the
hypergeometric survival function (verified in the test suite); at `k = 3` it
agrees with brute-force Monte Carlo. No resampling and no asymptotics means
calibration at small `n` by construction.

`m >= 2` is required for a positive call: one joint change edge is a single
clonal co-origin, the confounder the test exists to reject. The legacy Poisson
approximation (`method="poisson"`) is *conservative*, not anticonservative: at
`E = 20` with three traits changing on ten edges each it rejects at 0.5 % where
the exact test rejects at 4.0 %.

---

## 5. Choosing k, including k = 1

Every criterion sweeps `[1] + k_range`. A criterion whose sweep starts at 2
cannot report "there are no clusters", and on noise it will not: measured, such
a sweep returns k = 2 in 100 % of replicates where a sweep including 1 returns
k = 1 in 100 %.

| `k_select_primary` | rule |
|---|---|
| `mdl` | largest description-length gain; **k = 1 unless some k has positive gain** |
| `prediction_strength` | largest k with mean PS ≥ threshold, else 1 (Tibshirani & Walther 2005) |
| `gap` | 1-SE rule over a sweep including k = 1 (Tibshirani, Walther & Hastie 2001) |
| `bic_mixture` | BIC over Bernoulli-mixture fits including k = 1 |

**MDL.** `gain = L_null − (L_assign + L_profiles + L_residuals)` with
`L_assign = n log₂ k`, residuals coded under the per-cluster Krichevsky-Trofimov
estimate, and parametric complexity `(d/2) log₂ n` charged to **both** models
(Rissanen 1983; Grünwald 2007). Charging the clustered model twice that rate and the
null model nothing, which is the obvious shortcut, is not a code-length
comparison. `mdl_gain_fraction` puts the gain on
an interpretable scale: 0.63 on the planted control, 0.003 on the *Klebsiella*
cohort, because a large bit count can still be a negligible fraction.

**Gap statistic.** Observed and reference samples are clustered by the **same**
procedure, so estimator bias cancels; `null_labeler` is a required argument.
Clustering the two with different algorithms accounts for 92-105 % of the
reported Gap on null data. The reference is
independent Bernoulli matched to each feature's prevalence, a stated deviation
from the paper's uniform reference, which is not defined for binary data.

**Permutation calibration.** `mdl_null_calibration` re-runs the whole
distance → fusion → spectral chain on column-permuted data. It costs
O(n_perm·T·m·n³) and is therefore **opt-in** (`run(..., n_perm=199)`). The
p-value uses the Phipson & Smyth (2010) `(b+1)/(B+1)` correction and reports its
minimum attainable value; a permutation p-value can never legitimately be 0.

---

## 6. Diagnostics

**Layer influence** (`influence.layer_influence`). `delta_loo[l] = clip(1 −
ARI(fuse(all), fuse(all \ l)), 0, 1)`, which is interventional rather than
merely descriptive.
Normalised to weights whose Hill number of order 1 (Hill 1973, *Ecology*
54:427-432) is the **effective number of contributing layers**. `n_eff` alone
cannot separate "every layer agrees, so removing one changes nothing" from "no
layer matters": both give `delta_loo ≈ 0`. The `regime` field uses the solo
agreements to resolve it: `complementary`, `collapsed`, `redundant` and
`uninformative`. Only the last two set the gate. An optional permutation
test (row-permute one layer, re-fuse) gives each layer a p-value.

**Population structure** (`lineage.lineage_concordance`). The fraction of
same-lineage pairs that are also same-cluster, against a permutation null of the
cluster labels. AMI is the wrong instrument: with hundreds of singleton lineages
the chance correction dominates, and on the *Klebsiella* cohort AMI(cluster, ST)
= 0.057 while the pair concordance was 0.914 against a null of 0.766
(z = 11.5). `cluster_composition` reports the maximum single-lineage share per
cluster; `dereplicate_index` supports the one-isolate-per-lineage sensitivity
run.

**Baselines** (`baselines`). Every candidate partition, meaning the fusion, the
naive concatenation, each single layer and a Bernoulli-mixture latent class
model, is
scored on the same external criteria and on MDL gain. Agreement with the fused
partition cannot rank candidates, because it measures agreement with the tool's
own answer. If a single layer or a domain rule scores better on the external
criteria, that is the result.

**Collinearity.** `stats.effective_dimension` is the participation ratio of a
layer's correlation eigenspectrum: how many mutually uncorrelated columns the
block is worth. A set of `p` perfectly collinear columns has effective dimension
1 however large `p` is.

---

## 7. Descriptive profiles, labelled as such

`archetype_profiles.tsv` gives per (cluster, feature): Δp with a bootstrap
interval, Cohen's *h*, a Haldane-corrected log-odds-ratio, a Fisher exact
p-value with Benjamini-Hochberg FDR, and `is_defining_descriptive`
(|h| ≥ 0.5, BH reject, support ≥ 3; all three configurable).

These are **effect sizes, not inference**. The partition was estimated from the
same matrix, so the p-values do not have their nominal meaning under the null.
They rank features; §4 tests them. The output keys carry the `descriptive_`
prefix for that reason.

`label_recoverability` reports how well a partial feature vector recovers the
label the pipeline itself assigned. It has no external target and no null; a
value near 1.0 is the expected, uninformative result. It is reported with that
caveat attached rather than as predictive validation.

---

## 8. Lineage-aware surveillance reading

Sections 4 to 7 all describe a partition. This one describes an *agent*, and it
is the only part of the package a surveillance laboratory can run without ever
clustering anything.

### 8.1 Decomposing a difference

For two collections A and B, with lineage shares `w` and within-lineage rates
`p`, the Kitagawa (1955, JASA 50:1168-1194) identity splits the observed
difference exactly:

    P_A - P_B = sum_l (w_A - w_B)(p_A + p_B)/2      composition
              + sum_l (w_A + w_B)/2 (p_A - p_B)     within-lineage

This is the demographic ancestor of the Blinder-Oaxaca decomposition (1973) and
of Fairlie's binary-outcome extension (2005). The non-parametric form is used
rather than a fitted one because a sequence-type variable carries of the order
of a hundred levels on a few hundred isolates with a long singleton tail: a
logit on ST dummies is separated before it is asked anything.
`identity_residual` is emitted so exactness can be checked rather than trusted,
and `benchmarks/decomposition_vs_regression.py` measures the claim rather than
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

### 8.2 Two estimability gates

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
throughout and is not gated. Bias never exceeds 0.3 percentage points anywhere
in the grid, and under an exact null the type I error is 0.049.

**Label availability.** A lineage-resolved statistic is computed on the
isolates that carry a label and describes the collections only if labelling was
unrelated to the trait. Two conditions must hold together. If coverage differs
between the collections but labelling is unrelated to the trait, the same
selection applies to both arms and the difference is largely unaffected. If the
trait differs between labelled and unlabelled isolates but coverage is equal,
missingness costs precision rather than validity. Together they mean the
difference is computed between two differently selected subsets. Both are
tested by Fisher exact tests and reported in
`lineage_label_availability`.

On the shipped *S. suis* cohort the sequence type is missing for 219 of 677
isolates, 40 % of the later United Kingdom period carries one against 97 % of
the earlier, and untyped isolates carry about two more non-wild-type results
out of thirteen. Ceftiofur then appears to fall by 11.0 points on the labelled
subset while the collection rises by 9.5. The gate fires and names the reason.

### 8.3 The panel is a family

An antimicrobial panel is a family of tests whose members are not independent.
`decompose_panel` controls the false discovery rate within each component
family by Benjamini-Hochberg, separately for composition and for the
within-lineage rate because they answer different questions and are not
exchangeable, and reports the effective number of independent agents as the
participation ratio of the panel correlation eigenspectrum. On the shipped
thirteen-agent panel that number is 4.8, which is what the tetracycline,
macrolide-lincosamide, beta-lactam and remaining blocks amount to.

### 8.4 Two prevalence estimands

Prevalence per isolate answers clinical burden; prevalence per lineage answers
diversity. Neither is the true one, and a report that does not say which it
means is the problem. The interval on the per-lineage figure is a two-stage
cluster bootstrap: lineages are resampled first, because the estimand averages
over lineages and its uncertainty is dominated by how many distinct ones were
seen, and isolates are then resampled inside each drawn lineage so that a rate
measured on two isolates is not treated as exact. Because a lineage of size `m`
with rate `r` resampled to size `m` yields `Binomial(m, r) / m`, the second
stage is drawn exactly rather than by index.

### 8.5 Concentration of carriage

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

## 9. Reading a susceptibility panel at its recorded resolution

### 9.1 One likelihood, three interval widths

A non-wild-type call, a recorded minimum inhibitory concentration and a
censored reading are one likelihood at three interval widths on the log2
concentration scale. A call is `(-inf, c]` or `(c, inf)`; a dilution recorded
at well *w* is `(previous tested well, w]`; a reading on the lowest or highest
tested well is unbounded on that side. `censored.intervals_from_binary` and
`censored.intervals_from_mic` build the intervals and
`censored.censored_clonal_share` consumes nothing else, so the input mode
changes how much an observation carries and not what is estimated.

Two tests in `tests/test_censored.py` hold that claim to account. With
zero-width intervals the estimate must equal the classical one-way
variance-component ratio computed in closed form; and a concentration
dichotomised at a cut-off must produce the intervals the call constructor
produces, and therefore the same share.

### 9.2 The estimand, and how it differs from the clonal share

The quantity returned is `tau2 / (tau2 + sigma2)`, the share of the latent
log2 concentration variance carried by lineage. This is what a mixed model
reports and what the heritability literature compares against. It is **not**
the same quantity as `attribution.clonal_share`, which reports achievable
out-of-sample predictive skill and falls when lineages are small because a
group mean estimated from three isolates predicts badly however much variance
the grouping truly carries. A cohort should carry both; their ratio is a
property of the study design rather than of the biology.

### 9.3 The maximisation step

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

### 9.4 Identification, and two intervals

A single cut point recovers only the standardised distance from the cut, so
the residual scale is not identified. `censored.scale_is_identified` detects
that and the fit holds the residual standard deviation at one, which puts the
result on the liability scale of Dempster and Lerner (1950) and makes it
comparable with published heritabilities rather than with a raw log2 variance.

Three widths are reported and only one of them is an interval.

**The interval** inverts the variance ratio. For a one-way random model the
ratio of the between-lineage to the within-lineage mean square is a scaled F,
and inverting it gives an interval for the intraclass correlation. The
expectation-maximisation fit returns the same two components, so the same
inversion applies to interval data once the within-lineage degrees of freedom
are discounted by the share of each observation the interval left unresolved.
Over a grid of 24,000 simulated cohorts its coverage is 0.92 to 0.99.

**The cluster bootstrap** draws whole lineages, with bias-corrected and
accelerated limits. It is reported beside the interval and is not the interval:
on the same grid it covers a median 0.61 against a nominal 0.95. Resampling
thirty lineages does not reproduce the sampling distribution of a variance
component built from those thirty lineages, and the bias correction does not
repair that. It is kept because it assumes no distribution and therefore fails
differently from the F inversion, which is what a check is for.

**The likelihood-ratio width**, from the Gauss-Hermite marginal likelihood in
`censored.marginal_loglik`, is not a confidence statement at all. It measures
how sharply this cohort pins the share down, which is the quantity that
improves when a dilution is read rather than a call, and its coverage is
reported so that no reader mistakes it: 0.68 on an exact reading and 0.74 on a
dilution at thirty lineages, 0.91 to 0.93 on a single cut point. The
chi-squared calibration is optimistic exactly where the residual scale is well
determined and the whole of the uncertainty sits in the group means.

### 9.5 The coarsening assumption, and the refusal

Treating a reading on an end well as censored is an assumption about why the
reading is there, not an observation: it is the coarsened-at-random condition
of Heitjan and Rubin (1991), under which the coarsening is ignorable for
likelihood inference. Three things follow. `censored.panel_geometry` reports
the mass on each end well and refuses the point mode above five per cent. A
recorded operator column overrides the assumption wherever the laboratory
supplied one, because an operator is an observation and the heuristic is not.
And `censored.sensitivity_endpoints` computes the share under both readings,
so the width of that bracket can be set beside the bootstrap interval.

Two refusals follow from the same grid rather than from judgement.

A lineage every one of whose readings is one-sided has no identified latent
mean. Counting such lineages is the wrong gate, since on a susceptible agent a
few are normal; what moves the answer is how much of the cohort sits there.
`censored.CENSORED_GROUP_LIMIT` is set at one half, read off the censoring
sweep.

A single cut point in a tail leaves almost no contrast to divide. At a
non-wild-type prevalence near 0.08 the share reads 0.23 above a true zero and
its interval covers 0.08 of the time, while from 0.24 upward both behave.
`censored.SINGLE_CUT_PREVALENCE` therefore refuses a call outside 0.10 to 0.90.
A dilution at the same prevalence is unaffected, which is the practical
argument for reading the panel: on the shipped cohort two agents are refused as
calls and estimable as dilutions.

## 10. Anytime-valid evidence for a lineage effect

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
mean by the same empirical-Bayes factor the clonal share uses, and scored on
the held-out fold against the null probability maximised on that same fold.
Because the numerator parameters never saw the held-out rows, the expectation
under the null is at most one with no regularity conditions, and the running
value may be inspected as often as wanted. By Ville's inequality an e-value of
at least `1/alpha` is the anytime-valid counterpart of a p-value below `alpha`.

Values from independent batches multiply, forming a test martingale, which is
what makes stopping at any point legitimate; folds of one cohort share data and
are averaged instead, which is valid because any convex combination of e-values
is an e-value. Across a panel, `evalues.e_bh` applies the procedure of Wang and
Ramdas (2022), which controls the false discovery rate under arbitrary
dependence between the hypotheses. That is the property the panel needs:
cross-resistance makes a macrolide block behave as one trait, so the agents are
neither independent nor reliably positively dependent, and the usual
justification for Benjamini-Hochberg does not apply.

## 11. Reproducibility

Every stochastic stage draws from its own generator spawned off `--seed` via
`numpy.random.SeedSequence`, so adding a diagnostic cannot perturb the
clustering. `spectral_from_similarity` takes an explicit `random_state`.

Byte-identity of a floating-point spectral pipeline is not achievable across
BLAS and scikit-learn builds and is not claimed. What is claimed, and tested, is
invariance to input row order, stability of the selected k, and reproducibility
of the whole result within one environment at a fixed seed. `cluster_result.json`
is serialised with `allow_nan=False` and an explicit numpy encoder, so it is
valid RFC 8259 JSON that a non-Python parser can read.

## 12. The estimator on other species

`benchmarks/atlas_cross_species.py` runs `attribution.clonal_share` on NCBI
Pathogen Detection, which supplies both variables from one public release: the
lineage is the SNP cluster (`PDS_acc`) and the trait is non-susceptibility,
R or I against S, read from the `AST_phenotypes` metadata field. Nothing in it
depends on a cut-off derived in this package, and no clustering of ours takes
part, so what it exercises is the estimator and not the pipeline around it.

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
that is why the release accession is recorded rather than only the date.

## 13. Two questions a clonal share can answer

`attribution.clonal_share` and `censored.censored_clonal_share` report the
share a *new* draw of lineages from the species would show. Their information
about the between-lineage variance carries `G - 1` degrees of freedom, where
`G` is the number of lineages, and no number of extra isolates adds to it. A
cohort of thirty lineages therefore has a wide interval, and that is a property
of the question rather than a defect of the estimator.

`realised.realised_share` reports the other question: how much of the trait
variance sits between the lineages this collection actually holds. The lineage
effects are then fixed unknowns rather than a fresh draw, and the parameter is
`S_a^2 / (S_a^2 + sigma^2)` with `S_a^2 = sum_g (a_g - abar)^2 / (G - 1)`.

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
decision are measured in section S12.8 of the article supplement; the identity
is held by a contract test rather than by this paragraph. What the realised
interval buys is width, and validity under the unequal lineage sizes real
cohorts have.

### 13.1 What was measured, and what was refused

The operating characteristics are in section S12 of the article supplement and
are produced by `benchmarks/realised_calibration.py`. Two conditions were fixed
before that campaign: an anchor arm that re-derives the reference law from the
model statement rather than the implementation had to reproduce the calibration
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
distribution rather than the trait.

Two repairs were tried and are recorded in the campaign rather than deleted. A
Box-type deflation of the degrees of freedom by the estimated residual kurtosis
is no better than the exact interval on Gaussian residuals and destroys it on
heavy ones. A within-lineage percentile bootstrap, intended as an
assumption-free fallback, under-covers under the very process the method is
derived against. Neither is used.
