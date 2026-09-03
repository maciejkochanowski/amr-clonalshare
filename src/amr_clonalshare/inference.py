"""inference.py — valid post-clustering inference for binary trait panels.

The problem
-----------
Discovering clusters and then testing, on the same data, which features
"define" them is double dipping: under the null of no clusters the resulting
p-values are not p-values (Gao, Bien & Witten 2024, JASA 119:332-342,
doi:10.1080/01621459.2022.2116331; Chen & Witten 2023, JMLR 24(152):1-41).

What is and is not available for binary data
--------------------------------------------
The exact selective-inference results for clustering are Gaussian: Gao-Bien-
Witten for hierarchical clustering and Chen-Witten for k-means both condition
on a Gaussian model with known or estimated variance. They do not transfer to
0/1 presence-absence panels.

Data thinning (Neufeld, Dharamshi, Gao & Witten 2024, JMLR 25(57):1-35) splits
an observation into independent parts drawn from the same family, and would be
the natural tool - except that its Remark 8 states plainly that the requirement
of infinite divisibility "prevents us from thinning the Bernoulli or
categorical distributions". A single 0/1 call cannot be split. Coarse-graining
a layer to counts and thinning those is only valid if the counts really follow
a thinnable family with a **known** parameter, and Proposition 11 of the same
paper quantifies what happens when the parameter is wrong:

    cov(X1, X2) = (1 - eps) * r * ((1-p)/p)^2 * (1 - (r+1)/(r~+1)),

non-zero whenever ``r~ != r``: **negative** when ``r~ < r`` and positive
when ``r~ > r``. Estimating ``r`` from the data being split is
therefore not a technicality. Hivert, Agniel, Thiebaut & Hejblum (2024,
arXiv:2405.13591) reach the same conclusion empirically for post-clustering
differential analysis: thinning and fission are useful only at extreme
signal-to-noise, because the nuisance parameters must be estimated from the
structure being tested.

This module therefore offers three routes and is explicit about the assumption
each one carries.

``feature_split_test`` (default for binary panels)
    Partition the **features** at random into a discovery block and a test
    block; cluster on the discovery block only; test the held-out features
    against those labels. Nothing is assumed about marginal distributions and
    0/1 data needs no coarse-graining. The cost is power: each split uses half
    the features to find structure.

    **What this tests, exactly.** The null is

        H0: the held-out features are independent of the labels learned from
            the discovery features,

    which is implied by "the features are mutually independent" - **not** by
    "there are no clusters". The distinction is not academic on genotype data.
    Two measured failures, both shipped as controls in
    ``benchmarks/calibration_study.py``:

    * *Blocked features.* With features in correlated blocks of five and no
      clusters at all, a per-feature split rejects with family-wise error
      **1.000** at nominal 0.05, because a block straddling the two halves makes
      the discovery labels predict the held-out features for reasons unrelated
      to clustering. Passing ``groups`` (one block label per feature) assigns
      whole blocks to one side and restores validity exactly: **0.042** at
      within-block correlation 0.5 and **0.017** at 0.9.

    * *Latent continuum.* With a continuous latent gradient driving all features
      - a unimodal population with no discrete groups - the test rejects at
      essentially any sample size, and grouping does not help, because the
      dependence is global rather than block-local.

    So a small p-value here means "there is structure in this panel that
    clustering the discovery features detects", which is weaker than "there are
    discrete archetypes". The result key is named ``p_value_structure`` and
    ``structure_detected`` for that reason, and discreteness is a *separate*
    question answered by :func:`continuum_null_test` below, not by this test -
    and specifically **not** by ``discreteness_evidence``, which compares a
    mixture against independence and therefore calls a continuum discrete every
    time (see its docstring).

    Blocks should be the units of co-inheritance - operons, mobile elements,
    one-hot categorical variables - which are the same groups the distance
    computation needs collapsed. Declared groups are not enough on their own:
    they cover the dependence the analyst thought of, and on the *Klebsiella*
    config that left all fourteen acquired-resistance classes as singleton units
    despite pairwise correlations up to 0.81. :func:`correlation_blocks` unions
    the declared groups with the connected components of the thresholded
    correlation graph, and :func:`split_unit_diagnostics` reports the strongest
    dependence the design can still leak across, so ``block_aware`` means "the
    units are adequate for this matrix" rather than "the analyst declared
    something".

    Between-unit dependence is *not* checkable from the data being tested:
    between-block dependence and clustering are the same observable. What is
    reported is the measured cross-unit correlation, not a guarantee.

``nb_thin`` / ``binomial_thin`` (counts, when the family is checkable)
    Count splitting for over-dispersed counts (negative binomial, ``r``
    per feature) and for bounded row-sums over ``m`` known members
    (binomial, split by the hypergeometric of JMLR 25(57) Table 2). The
    binomial route is the correct one for "how many of these ``m`` loci does
    this isolate carry", and unlike the NB route its size parameter ``m`` is
    known by construction rather than estimated. :func:`thinning_dependence`
    reports the realised ``corr(X1, X2)`` so the guarantee is checked rather
    than assumed.

Not offered: sample splitting with a classifier trained on the same features.
It is the intuitive thing to reach for and it is catastrophically invalid here
(family-wise error rate 1.00 in the calibration study shipped in
``benchmarks/``), because the classifier transfers the discovery-half labels'
dependence on the features straight into the test half.

Merging the splits
------------------
One random split gives one p-value, and that p-value is a random variable even
with the data fixed. Reporting it is not reproducible: on the *Klebsiella*
example the previous single-split implementation returned p between 3.8e-63 and
0.97 across seeds. Repeated splits of one fixed dataset are **exchangeable**
but arbitrarily dependent, so Fisher and Stouffer combination are invalid and
taking the minimum is catastrophically so. A multiplicative constant is
unavoidable. Three valid rules are provided:

``"exchangeable_ruger"`` (default)
    Gasparin, Wang & Ramdas (2025), PNAS 122(11):e2410849122, Theorem 4.1:
    ``F(p) = (K/k) * min_{l<=K} p^{(l)}_{(ceil(l k / K))}``, the exchangeable
    improvement of "twice the median". Strictly dominates classical Ruger.
``"ruger"``
    ``(K/k) p_(k)``; valid under arbitrary dependence.
``"twice_mean"``
    ``2 * mean(p)``; sharp under arbitrary dependence (Vovk & Wang 2020,
    Biometrika 107:791-808, doi:10.1093/biomet/asaa027).
``"mmb"``
    Meinshausen, Meier & Buhlmann (2009), JASA 104:1671-1681,
    doi:10.1198/jasa.2009.tm08647, Eq. 2.3: the quantile-adaptive rule
    ``min{1, (1 - log gamma_min) inf_gamma q_gamma(p/gamma)}``. Valid under
    arbitrary dependence; the most conservative of the four.

Measured on a correct null (see ``benchmarks/calibration_study.py``), the
default rule holds its level and is *more* powerful than a single split at the
detection boundary, because merging removes the randomisation noise.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import beta, fisher_exact, hypergeom, kruskal, mannwhitneyu

from .stats import benjamini_hochberg

__all__ = [
    "continuum_null_test",
    "one_factor_loglik",
    "latent_trait_loglik",
    "fit_one_factor",
    "fit_latent_trait",
    "select_latent_dimension",
    "merge_pvalues",
    "merge_exchangeable",
    "merge_order_sensitivity",
    "ruger_merge",
    "twice_mean_merge",
    "mmb_merge",
    "nb_thin",
    "binomial_thin",
    "poisson_thin",
    "thinning_dependence",
    "feature_split_test",
    "feature_split_report",
    "correlation_blocks",
    "split_unit_diagnostics",
]


# --------------------------------------------------------------------------- #
# 1. Merging exchangeable p-values
# --------------------------------------------------------------------------- #
def ruger_merge(pvals: Sequence[float], k: Optional[int] = None) -> float:
    """Classical Ruger bound ``(K/k) * p_(k)``; valid under arbitrary dependence."""
    p = np.sort(np.asarray(list(pvals), dtype=float))
    K = p.size
    if K == 0:
        return 1.0
    k = int(np.clip(k if k is not None else math.ceil(K / 2), 1, K))
    return float(min(1.0, (K / k) * p[k - 1]))


def twice_mean_merge(pvals: Sequence[float]) -> float:
    """``2 * mean(p)``; sharp under arbitrary dependence (Vovk & Wang 2020)."""
    p = np.asarray(list(pvals), dtype=float)
    return float(min(1.0, 2.0 * p.mean())) if p.size else 1.0


def merge_exchangeable(pvals: Sequence[float], k: Optional[int] = None) -> float:
    """Exchangeable Ruger improvement (Gasparin, Wang & Ramdas 2025, Thm 4.1).

    **Order-dependent by construction.** ``p^(l)_(...)`` is a quantile of the
    *first* ``l`` p-values in the order supplied, so the rule is not a function
    of the multiset: re-seeding the split loop permutes the order and moves the
    value. Over 500 random vectors of 9 exchangeable p-values evaluated under 30
    orderings each, the median max/min ratio is 5.75 and the 95th percentile
    42.5. Validity is unaffected -- the theorem holds for any *fixed* order --
    but reproducibility is, so the reported value is the one for the seeded
    split order and :func:`merge_order_sensitivity` reports the spread.
    Minimising over orderings would break validity and must not be done.
    """
    p = np.asarray(list(pvals), dtype=float)
    K = p.size
    if K == 0:
        return 1.0
    bad = p[~np.isfinite(p) | (p < 0.0) | (p > 1.0)]
    if bad.size:
        # This is the default merge rule feeding Benjamini-Hochberg. A value
        # outside [0, 1] came back through the rule unchanged -- [-1.0, 0.1,
        # 0.1] merged to -1.5 -- and a NaN sorted to the end of every prefix
        # and was silently dropped from the minimum.
        raise ValueError(
            f"p-value {float(bad[0])!r} is not a probability; "
            f"merge_exchangeable needs every input finite and in [0, 1]"
        )
    if K == 1:
        return float(np.clip(p[0], 0.0, 1.0))
    k = int(np.clip(k if k is not None else math.ceil(K / 2), 1, K))
    best = np.inf
    for l in range(1, K + 1):
        j = int(np.clip(math.ceil(l * k / K), 1, l))
        best = min(best, float(np.sort(p[:l])[j - 1]))
    return float(np.clip((K / k) * best, 0.0, 1.0))


def merge_order_sensitivity(pvals: Sequence[float], *, merge: str = "exchangeable_ruger",
                            n_orders: int = 50, seed: int = 0) -> Dict[str, float]:
    """Spread of a merged p-value over random orderings of the same split p-values.

    ``exchangeable_ruger`` is valid for any fixed order but is not symmetric in
    its arguments, so the reported value depends on the order the splits happened
    to be drawn in. This reports what that costs: the merged value under the
    order given, and the range over ``n_orders`` random permutations. A large
    ratio does not invalidate the reported p-value; it means the reader should
    read it as one draw from that range rather than as a fixed number.
    """
    p = np.asarray(list(pvals), dtype=float)
    rng = np.random.default_rng(seed)
    vals = [merge_pvalues(p, merge)]
    for _ in range(int(n_orders)):
        vals.append(merge_pvalues(rng.permutation(p), merge))
    v = np.asarray(vals, dtype=float)
    lo, hi = float(v.min()), float(v.max())
    return {"as_reported": float(vals[0]), "min_over_orders": lo,
            "max_over_orders": hi, "n_orders": int(n_orders),
            "max_over_min_ratio": float(hi / lo) if lo > 0 else float("inf")}


def mmb_merge(pvals: Sequence[float], gamma_min: float = 0.05) -> float:
    """Meinshausen-Meier-Buhlmann (2009) quantile-adaptive merge, Eq. 2.3.

    ``min{1, (1 - log gamma_min) * inf_{gamma in [gamma_min, 1]}
    q_gamma({p_b / gamma})}`` where ``q_gamma`` is the empirical gamma-quantile.
    Valid for arbitrarily dependent, marginally super-uniform p-values; the most
    conservative rule offered here.
    """
    p = np.asarray(list(pvals), dtype=float)
    K = p.size
    if K == 0:
        return 1.0
    grid = np.unique(np.clip(np.linspace(gamma_min, 1.0, 200), gamma_min, 1.0))
    best = np.inf
    for g in grid:
        q = float(np.quantile(p / g, g, method="higher"))
        best = min(best, q)
    return float(min(1.0, (1.0 - math.log(gamma_min)) * best))


_MERGERS = {
    "exchangeable_ruger": merge_exchangeable,
    "ruger": ruger_merge,
    "twice_mean": twice_mean_merge,
    "mmb": mmb_merge,
}


def merge_pvalues(pvals: Sequence[float], how: str = "exchangeable_ruger") -> float:
    """Merge exchangeable p-values by the named rule."""
    try:
        return float(_MERGERS[how](pvals))
    except KeyError:
        raise ValueError(
            f"unknown merge rule {how!r}; choose from {sorted(_MERGERS)}") from None


# --------------------------------------------------------------------------- #
# 2. Count splitting (only where the family is checkable)
# --------------------------------------------------------------------------- #
def poisson_thin(X, *, eps: float = 0.5,
                 rng: Optional[np.random.Generator] = None
                 ) -> Tuple[np.ndarray, np.ndarray]:
    """Poisson splitting ``X1 ~ Binomial(X, eps)``. Valid only for Poisson counts."""
    rng = rng or np.random.default_rng()
    Xi = np.asarray(X, dtype=np.int64)
    X1 = rng.binomial(Xi, eps)
    return X1, Xi - X1


def nb_thin(X, *, r, eps: float = 0.5,
            rng: Optional[np.random.Generator] = None
            ) -> Tuple[np.ndarray, np.ndarray]:
    """Negative-binomial splitting (JMLR 25(57) Table 2).

    ``X1 | X = x ~ BetaBinomial(x, eps*r, (1-eps)*r)``. ``r`` may be scalar or
    per-column; **the theorem assumes it is known**, so an estimated ``r``
    yields only approximate independence - use :func:`thinning_dependence` to
    measure the residual dependence rather than assuming it away.
    """
    rng = rng or np.random.default_rng()
    if not (0.0 < eps < 1.0):
        raise ValueError(f"eps must lie in (0, 1); got {eps}")
    Xi = np.asarray(X, dtype=np.int64)
    r_arr = np.atleast_1d(np.asarray(r, dtype=float))
    if r_arr.size == 1:
        r_arr = np.full(Xi.shape[1], float(r_arr[0]))
    if r_arr.size != Xi.shape[1]:
        raise ValueError(f"r must be scalar or length {Xi.shape[1]}")
    if np.any(r_arr <= 0):
        raise ValueError("all dispersions r must be positive")
    p_split = rng.beta(np.broadcast_to(eps * r_arr, Xi.shape),
                       np.broadcast_to((1 - eps) * r_arr, Xi.shape))
    X1 = rng.binomial(Xi, p_split)
    return X1, Xi - X1


def binomial_thin(X, *, m, eps: float = 0.5,
                  rng: Optional[np.random.Generator] = None
                  ) -> Tuple[np.ndarray, np.ndarray]:
    """Binomial splitting by the hypergeometric (JMLR 25(57) Table 2).

    For ``X ~ Binomial(m, p)`` with ``m`` **known**, drawing
    ``X1 | X = x ~ Hypergeometric(m, x, eps*m)`` gives
    ``X1 ~ Binomial(eps*m, p)`` independent of ``X2 ~ Binomial((1-eps)*m, p)``.

    This is the right split for "how many of these ``m`` loci does the isolate
    carry": ``m`` is the number of member columns, known by construction, so
    unlike the negative-binomial route nothing is estimated from the data being
    split. It assumes the ``m`` members are exchangeable Bernoullis within an
    isolate, which is false for a tightly co-inherited operon - collapse such a
    group to a single presence call first.

    ``m`` may be scalar or per-column.
    """
    rng = rng or np.random.default_rng()
    if not (0.0 < eps < 1.0):
        raise ValueError(f"eps must lie in (0, 1); got {eps}")
    Xi = np.asarray(X, dtype=np.int64)
    m_arr = np.atleast_1d(np.asarray(m, dtype=np.int64))
    if m_arr.size == 1:
        m_arr = np.full(Xi.shape[1], int(m_arr[0]), dtype=np.int64)
    if m_arr.size != Xi.shape[1]:
        raise ValueError(f"m must be scalar or length {Xi.shape[1]}")
    if np.any(Xi > m_arr[None, :]):
        raise ValueError("counts exceed the stated binomial size m")
    m1 = np.maximum(1, np.minimum(m_arr - 1, np.round(eps * m_arr).astype(np.int64)))
    X1 = np.empty_like(Xi)
    for j in range(Xi.shape[1]):
        X1[:, j] = hypergeom.rvs(int(m_arr[j]), Xi[:, j], int(m1[j]),
                                 random_state=rng)
    return X1, Xi - X1


def thinning_dependence(X, X1, X2) -> np.ndarray:
    """Per-column ``corr(X1, X2)``: the direct check of the thinning guarantee.

    Should be ~0. A systematically non-zero value means the assumed family or
    parameter is wrong (JMLR 25(57) Prop. 11) and the downstream p-values are
    not valid.
    """
    X1 = np.asarray(X1, float)
    X2 = np.asarray(X2, float)
    out = np.full(X1.shape[1], np.nan)
    for j in range(X1.shape[1]):
        if X1[:, j].std() > 0 and X2[:, j].std() > 0:
            out[j] = float(np.corrcoef(X1[:, j], X2[:, j])[0, 1])
    return out


# --------------------------------------------------------------------------- #
# 3. Feature splitting (the default route for binary panels)
# --------------------------------------------------------------------------- #
def _cluster_block(block: np.ndarray, k: int, cluster_fn, seed: int,
                   cols: Optional[np.ndarray] = None) -> Optional[np.ndarray]:
    """Cluster the discovery block. ``cluster_fn`` may accept a ``cols=`` kwarg
    naming the original column indices, so it can apply per-layer distances."""
    if block.shape[1] < 1 or block.shape[0] < 2 * k:
        return None
    try:
        try:
            labels = cluster_fn(block, k, seed, cols=cols)
        except TypeError:
            labels = cluster_fn(block, k, seed)
    except Exception:
        return None
    labels = np.asarray(labels)
    uniq, cnt = np.unique(labels, return_counts=True)
    if uniq.size < 2 or cnt.min() < 2:
        return None
    return labels


def _abs_corr(X: np.ndarray) -> np.ndarray:
    """Absolute feature-feature correlation, with constant columns set to zero.

    For 0/1 data the Pearson correlation is the phi coefficient.
    """
    A = np.asarray(X, dtype=float)
    sd = A.std(axis=0)
    ok = sd > 0
    R = np.zeros((A.shape[1], A.shape[1]))
    if ok.sum() > 1:
        sub = np.corrcoef(A[:, ok], rowvar=False)
        sub = np.nan_to_num(np.abs(sub), nan=0.0)
        R[np.ix_(ok, ok)] = sub
    np.fill_diagonal(R, 1.0)
    return R


def correlation_blocks(X: np.ndarray, *, threshold: float = 0.5,
                       groups: Optional[Sequence] = None) -> np.ndarray:
    """Split units within which features may be dependent and between which they are not.

    Multi-split feature splitting is only valid if the held-out features are
    independent of the discovery features *under the null*. Splitting
    per-feature on a panel whose features are correlated leaks the labels into
    the held-out half; on the calibration in ``benchmarks/calibration_study.py``
    that costs family-wise error 1.000 against a nominal 0.05. Declared
    co-inheritance groups fix the part of the dependence the analyst knew
    about; this function fixes the rest.

    Two features are joined when ``|corr| >= threshold``, and each split unit is
    a **connected component** of the resulting graph, seeded with ``groups`` so
    that declared blocks are never broken up. Single linkage is the right
    linkage here and not an arbitrary choice: it is the only one that makes the
    property the test actually needs -- every *between*-unit ``|corr|`` is below
    ``threshold`` -- true by construction. Its cost is chaining, which shows up
    as a small number of large units and is reported rather than hidden
    (:func:`split_unit_diagnostics`).

    Choosing the blocks uses the feature-by-feature correlation matrix only. It
    never touches the cluster labels, so it does not spend the selective
    inference the split exists to protect.

    Returns one integer block label per column of ``X``.
    """
    A = np.asarray(X)
    p = A.shape[1]
    parent = list(range(p))

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    if groups is not None:
        if len(list(groups)) != p:
            raise ValueError(f"groups must have one entry per feature ({p})")
        first: Dict[object, int] = {}
        for j, g in enumerate(groups):
            if g is None:
                continue
            if g in first:
                union(first[g], j)
            else:
                first[g] = j

    if threshold is not None and 0.0 < float(threshold) <= 1.0 and p > 1:
        R = _abs_corr(A)
        for a, b in zip(*np.triu_indices(p, k=1)):
            if R[a, b] >= float(threshold):
                union(int(a), int(b))

    return pd.factorize(pd.Series([find(j) for j in range(p)]))[0]


def split_unit_diagnostics(X: np.ndarray, key: Sequence[int]) -> Dict[str, float]:
    """Is this set of split units actually adequate for the dependence in ``X``?

    ``max_abs_corr_between_units`` is the quantity that decides validity: it is
    the strongest dependence the split can leak across. A design whose units are
    all singletons reports the panel's largest off-diagonal correlation here,
    which is what makes an all-singleton design visibly inadequate rather than
    silently so.
    """
    k = np.asarray(list(key))
    p = k.size
    sizes = np.bincount(k) if p else np.zeros(0, dtype=int)
    R = _abs_corr(X) if p > 1 else np.ones((p, p))
    cross = R[np.triu_indices(p, k=1)][
        (k[:, None] != k[None, :])[np.triu_indices(p, k=1)]] if p > 1 else np.zeros(0)
    within = R[np.triu_indices(p, k=1)][
        (k[:, None] == k[None, :])[np.triu_indices(p, k=1)]] if p > 1 else np.zeros(0)
    return {
        "n_split_units": int(sizes.size),
        "n_multi_feature_units": int((sizes > 1).sum()),
        "largest_unit_size": int(sizes.max()) if sizes.size else 0,
        "max_abs_corr_between_units": float(cross.max()) if cross.size else 0.0,
        "mean_abs_corr_within_units": float(within.mean()) if within.size else 0.0,
    }


#: Cross-unit |corr| above this makes a split design inadequate, and the default
#: ``tva.block_threshold``. Swept in ``benchmarks/calibration_study.py``
#: experiment G. Under a no-cluster null with correlated blocks the type-I error
#: of the merged test is 0.9998 with no blocking, 0.0163 at 0.5, 0.0195 at 0.7
#: and 0.7897 at 0.85 (blocks stop forming), over 4000 replicates. Counts,
#: Clopper-Pearson intervals and the campaign receipt are in
#: ``benchmarks/results_n4000_2026-08-14/``.
#:
#: That sweep does **not** separate 0.5 from 0.7. Both hold the nominal 0.05,
#: and 0.5 is marginally the tighter of the two. A power justification for
#: this default was attributed to experiment G and is withdrawn: experiment G
#: returns type-I error, unit counts and between-unit correlation
#: and does not measure power at all, so those figures are dropped rather
#: than restated. What protects an analysis is ``split_design_adequate`` below,
#: which is evaluated on the user's own design on every run and enforced in
#: ``cli.py``.
#:
#: The choice of 0.7 over 0.5 is made on the shipped cohorts rather than by
#: simulation. ``benchmarks/results_adequacy_2026-08-14/`` sweeps this value
#: from 0.30 to 0.99 and records where ``split_design_adequate`` holds: 0.30
#: to 0.75 on the Klebsiella case study, 0.49 to 0.72 on the clonal control,
#: everywhere on the null control. Both level-verified values sit inside the
#: intersection 0.49-0.72; 0.7 gives the finer units (15 against 10, and 20
#: against 11) and more margin, since at 0.5 the largest unit on the clonal
#: control holds exactly the 60 % the chaining bound allows. Power is the one
#: axis on which the two are still not compared.
BLOCK_ADEQUACY_MAX_CROSS_CORR = 0.7

#: Nodes per dimension for the matched-grid sensitivity check in
#: :func:`select_latent_dimension`. Small enough that q = 3 stays affordable.
_MATCHED_NODES = 9

#: A split needs at least two units on each side to be worth drawing.
MIN_SPLIT_UNITS = 4
#: Single linkage chains. When one unit swallows most of the panel the design
#: stays valid but can no longer answer much, so it is refused rather than
#: reported. This is a geometric condition on the units the split produced. It
#: is not a measured power loss, and nothing in this package measures power as
#: a function of ``tva.block_threshold``.
MAX_FRAC_IN_ONE_UNIT = 0.6


def _block_report(X: np.ndarray, key: Sequence[int],
                  units: np.ndarray) -> Dict[str, object]:
    """The split-unit fields every split-based result must carry.

    ``block_aware`` is not ``groups is not None``, which is true even when
    every unit is a singleton, that is, true in exactly the case it is meant to
    rule out. It reports whether the strongest correlation between two
    *different* units is below the calibrated bound.

    That is necessary but, on its own, close to vacuous: when
    ``tva.block_threshold > 0``, :func:`correlation_blocks` builds the units by
    single linkage at that threshold, so every cross-unit correlation is below
    it **by construction**. ``block_aware`` can therefore only fail when blocks
    are declared by hand, or when the threshold is set above the calibrated
    bound. A flag that cannot fail in the default configuration is a restatement
    of the configuration, not a diagnostic.

    The falsifiable companion is ``split_design_adequate``, which asks whether
    blocking left a design that can still answer anything. Single linkage
    chains, and when one unit swallows the panel what is left is valid but
    close to uninformative. The failure is real and reachable at the default
    settings: the synthetic planted control in ``examples/synthetic/`` builds
    only three split units, is refused there, and has its inference withheld
    with ``inference_status = "withheld_inadequate_split_design"``. ``cli.py``
    enforces that refusal rather than reporting it as advice.

    A 60-feature, 4-unit chaining example at threshold 0.5 was attributed to
    experiment G and is withdrawn. No generator, fixture or result file in this
    package produces it: experiment G runs at p = 30 and reports six units at
    that threshold, and the shipped planted control has 48 columns. It is
    dropped rather than re-attributed.
    """
    diag = split_unit_diagnostics(X, key)
    corr_ok = diag["max_abs_corr_between_units"] < BLOCK_ADEQUACY_MAX_CROSS_CORR
    n_units = int(diag["n_split_units"])
    p_feat = int(np.asarray(list(key)).size)
    frac_largest = (diag["largest_unit_size"] / p_feat) if p_feat else 1.0
    enough_units = n_units >= MIN_SPLIT_UNITS
    no_chaining = frac_largest <= MAX_FRAC_IN_ONE_UNIT
    adequate = bool(corr_ok and enough_units and no_chaining)

    reasons = []
    if not corr_ok:
        reasons.append(
            f"two features in different split units correlate at "
            f"{diag['max_abs_corr_between_units']:.2f}, at or above the "
            f"calibrated bound {BLOCK_ADEQUACY_MAX_CROSS_CORR}: the held-out "
            "half is not independent of the discovery half under the null, and "
            "per-feature splitting on correlated features has family-wise error "
            "1.000 (calibration study, experiment B)")
    if not enough_units:
        reasons.append(
            f"only {n_units} split units (need >= {MIN_SPLIT_UNITS}); a split "
            "cannot put two units on each side, so the discovery half is too "
            "small to find structure and the held-out half too small to test it")
    if not no_chaining:
        reasons.append(
            f"{frac_largest:.0%} of the features are in a single unit "
            f"(bound {MAX_FRAC_IN_ONE_UNIT:.0%}); single-linkage chaining has "
            "collapsed the design, and whichever side that unit lands on "
            "carries essentially all the signal")

    return {
        **diag,
        "frac_features_in_largest_unit": float(frac_largest),
        "block_aware": bool(corr_ok),
        "split_design_adequate": adequate,
        "split_design_problems": reasons,
        "block_adequacy_threshold": float(BLOCK_ADEQUACY_MAX_CROSS_CORR),
        "validity_note": (
            "the split is drawn over whole feature blocks; the strongest "
            f"correlation between two blocks is "
            f"{diag['max_abs_corr_between_units']:.2f} and the largest block "
            f"holds {frac_largest:.0%} of the features. Block-aware splitting "
            "has family-wise error 0.000-0.067 against a nominal 0.05 under a "
            "no-cluster null with correlated features"
            if adequate else
            "INVALID OR UNUSABLE SPLIT DESIGN: " + "; ".join(reasons)),
    }


def _split_units(p: int, groups: Optional[Sequence]) -> Tuple[np.ndarray, np.ndarray]:
    """Map each feature to a split unit. Without ``groups`` each feature is its own."""
    if groups is None:
        key = np.arange(p)
    else:
        key = pd.factorize(pd.Series(list(groups), dtype="object")
                           .fillna("__ungrouped__").astype(str))[0]
        if key.size != p:
            raise ValueError(f"groups must have one entry per feature ({p})")
    return key, np.unique(key)


def _draw_split(key: np.ndarray, units: np.ndarray, frac: float,
                rng: np.random.Generator) -> Tuple[np.ndarray, np.ndarray]:
    perm = rng.permutation(units)
    n_disc = int(np.clip(round(frac * units.size), 1, units.size - 1))
    disc_units = set(perm[:n_disc].tolist())
    disc = np.flatnonzero(np.isin(key, list(disc_units)))
    held = np.flatnonzero(~np.isin(key, list(disc_units)))
    return disc, held


def feature_split_test(X: np.ndarray, k: int, *, cluster_fn,
                       rng: Optional[np.random.Generator] = None,
                       n_splits: int = 25, frac_discovery: float = 0.5,
                       merge: str = "exchangeable_ruger",
                       binary: bool = True,
                       groups: Optional[Sequence] = None) -> Dict[str, object]:
    """Multi-split feature-split test of "there is real k-group structure".

    Each split: draw a random discovery block of features, cluster on it, and
    test the **held-out** features against those labels. For binary features the
    per-feature statistic is a Fisher exact test of the largest cluster against
    the rest (grouping chosen from the discovery block, never from the held-out
    block); the split-level statistic is the minimum held-out p-value corrected
    by Bonferroni over the held-out features, which is valid under arbitrary
    dependence among them. Split p-values are then merged by ``merge``.

    Parameters
    ----------
    X : (n, p) data matrix. Binary unless ``binary=False``.
    k : number of groups to look for.
    cluster_fn : ``(block, k, seed) -> labels`` applied to the discovery block.
    n_splits : number of random feature splits.
    frac_discovery : fraction of features used for discovery.
    merge : merging rule, see :func:`merge_pvalues`.
    """
    rng = rng or np.random.default_rng()
    A = np.asarray(X)
    n, p = A.shape
    key, units = _split_units(p, groups)
    if p < 4 or units.size < 2:
        return {"status": "skipped",
                "reason": f"feature splitting needs >= 4 features in >= 2 split "
                          f"units; got p={p}, units={units.size}",
                "p_value": None}
    per_split = []
    n_used = 0
    for s in range(int(n_splits)):
        disc, held = _draw_split(key, units, frac_discovery, rng)
        if disc.size < 1 or held.size < 1:
            per_split.append(1.0)
            continue
        labels = _cluster_block(A[:, disc], k, cluster_fn, seed=s, cols=disc)
        if labels is None:
            per_split.append(1.0)
            continue
        n_used += 1
        pv = _held_out_pvalues(A[:, held], labels, binary=binary)
        finite = pv[np.isfinite(pv)]
        per_split.append(float(min(1.0, finite.min() * finite.size))
                         if finite.size else 1.0)
    if n_used == 0:
        # Every split failed to produce usable labels - most often because
        # `cluster_fn` has the wrong signature, which `_cluster_block` catches
        # and turns into a `None`. Returning p = 1.0 here would be a clean
        # "no structure" verdict for a test that never ran.
        return {"status": "failed",
                "reason": "no split produced usable cluster labels; check that "
                          "cluster_fn accepts (block, k, seed) and that each "
                          "discovery block has >= 2k rows",
                "n_splits": int(n_splits), "n_splits_usable": 0,
                "p_value": None}
    merged = merge_pvalues(per_split, merge)
    arr = np.asarray(per_split, float)
    return {
        "status": "ok",
        "p_value": merged,
        "merge_rule": merge,
        "n_splits": int(n_splits),
        "n_splits_usable": int(n_used),
        **_block_report(X, key, units),
        "frac_discovery": float(frac_discovery),
        "p_per_split": [float(v) for v in per_split],
        "p_split_median": float(np.median(arr)),
        "p_split_min": float(arr.min()),
        "p_split_max": float(arr.max()),
        "k": int(k),
    }


def _held_out_pvalues(block: np.ndarray, labels: np.ndarray, *, binary: bool
                      ) -> np.ndarray:
    """Per-feature p-values on the held-out block given discovery-block labels.

    The grouping is a function of the discovery block alone, so no selection
    effect leaks into the held-out test.
    """
    uniq, cnt = np.unique(labels, return_counts=True)
    top = int(uniq[np.argmax(cnt)])          # chosen from the discovery labels
    mask = labels == top
    out = np.ones(block.shape[1])
    if mask.sum() < 2 or (~mask).sum() < 2:
        return out
    for j in range(block.shape[1]):
        col = block[:, j]
        if binary:
            n11 = int(col[mask].sum())
            n10 = int(mask.sum() - n11)
            n01 = int(col[~mask].sum())
            n00 = int((~mask).sum() - n01)
            if (n11 + n01) == 0 or (n10 + n00) == 0:
                continue
            out[j] = float(fisher_exact([[n11, n10], [n01, n00]],
                                        alternative="two-sided")[1])
        else:
            groups = [col[labels == c] for c in uniq if (labels == c).sum() >= 2]
            if len(groups) < 2 or np.ptp(np.concatenate(groups)) == 0:
                continue
            try:
                out[j] = float(kruskal(*groups).pvalue if len(groups) > 2
                               else mannwhitneyu(groups[0], groups[1],
                                                 alternative="two-sided").pvalue)
            except ValueError:
                out[j] = 1.0
    return out


def feature_split_report(X: pd.DataFrame, k: int, *, cluster_fn,
                         rng: Optional[np.random.Generator] = None,
                         n_splits: int = 25, frac_discovery: float = 0.5,
                         merge: str = "exchangeable_ruger",
                         q_fdr: float = 0.05, binary: bool = True,
                         groups: Optional[Sequence] = None,
                         fdr_features: Optional[Sequence[str]] = None,
                         reference_labels: Optional[Sequence[int]] = None
                         ) -> Dict[str, object]:
    """Global test plus per-feature FDR-controlled defining features.

    A feature is tested only on the splits in which it was held out; its merged
    p-value uses those splits, and Benjamini-Hochberg is applied across
    features. Features held out too rarely to be testable are reported with
    ``n_tested`` so the reader can see the coverage.
    """
    if not isinstance(X, pd.DataFrame):
        raise TypeError("feature_split_report expects a DataFrame")
    rng = rng or np.random.default_rng()
    A = X.values
    n, p = A.shape
    key, units = _split_units(p, groups)
    if p < 4 or n < 2 * k or units.size < 2:
        return {"status": "skipped",
                "reason": f"need >= 4 features in >= 2 split units and >= "
                          f"{2 * k} isolates; got p={p}, units={units.size}, n={n}",
                "structure_detected": None}

    split_p: List[float] = []
    disc_ari: List[float] = []
    dir_acc = np.zeros(p)
    dir_cnt = np.zeros(p)
    per_feature: Dict[int, List[float]] = {j: [] for j in range(p)}
    n_used = 0
    for s in range(int(n_splits)):
        disc, held = _draw_split(key, units, frac_discovery, rng)
        labels = (_cluster_block(A[:, disc], k, cluster_fn, seed=s, cols=disc)
                  if disc.size and held.size else None)
        if labels is None:
            split_p.append(1.0)
            continue
        n_used += 1
        if reference_labels is not None:
            from sklearn.metrics import adjusted_rand_score as _ari
            disc_ari.append(float(_ari(np.asarray(reference_labels), labels)))
        pv = _held_out_pvalues(A[:, held], labels, binary=binary)
        for idx, j in enumerate(held):
            per_feature[int(j)].append(float(pv[idx]))
        # Direction, from the DISCOVERY labels, in the same pass: no extra cost
        # and no selection effect, since it never touches the tested values.
        uq, cn = np.unique(labels, return_counts=True)
        top = labels == int(uq[np.argmax(cn)])
        if top.sum() >= 2 and (~top).sum() >= 2:
            dir_acc += np.sign(A[top].mean(axis=0) - A[~top].mean(axis=0))
            dir_cnt += 1
        finite = pv[np.isfinite(pv)]
        split_p.append(float(min(1.0, finite.min() * finite.size))
                       if finite.size else 1.0)

    if n_used == 0:
        # Same guard as feature_split_test: every split failed to produce usable
        # labels, so merging the 1.0 placeholders would report "no structure"
        # for a test that never ran.
        return {"status": "failed",
                "reason": "no split produced usable cluster labels; check that "
                          "cluster_fn accepts (block, k, seed) and that each "
                          "discovery block has >= 2k rows",
                "n_splits": int(n_splits), "n_splits_usable": 0,
                "structure_detected": None}

    merged_global = merge_pvalues(split_p, merge)
    feat_p = np.ones(p)
    n_tested = np.zeros(p, dtype=int)
    for j in range(p):
        vals = per_feature[j]
        n_tested[j] = len(vals)
        feat_p[j] = merge_pvalues(vals, merge) if vals else 1.0
    # Benjamini-Hochberg over the hypotheses that actually exist. Mutually
    # exclusive indicator columns (a one-hot block) are one hypothesis, not one
    # per level, so the caller can restrict the FDR set with `fdr_features`.
    cols_all = list(X.columns)
    if fdr_features is None:
        in_fdr = np.ones(p, dtype=bool)
    else:
        keep = set(map(str, fdr_features))
        in_fdr = np.array([str(c) in keep for c in cols_all], dtype=bool)
    adj = np.ones(p)
    reject = np.zeros(p, dtype=bool)
    if in_fdr.any():
        a_sub, r_sub = benjamini_hochberg(feat_p[in_fdr], q=q_fdr)
        adj[in_fdr] = a_sub
        reject[in_fdr] = r_sub
    arr = np.asarray(split_p, float)
    # Direction accumulated in the loop above from the DISCOVERY half only, so a
    # validated feature says which side of the split it belongs to without
    # carrying any selection effect into the test.
    with np.errstate(invalid="ignore", divide="ignore"):
        direction = np.sign(np.where(dir_cnt > 0,
                                     dir_acc / np.maximum(dir_cnt, 1),
                                     0.0)).astype(int)
    table = pd.DataFrame({
        "feature": cols_all,
        "p_value": feat_p, "p_value_bh": adj, "is_defining": reject,
        "n_splits_tested": n_tested,
        "in_fdr_set": in_fdr,
        "enriched_in_larger_group": direction,
    })
    return {
        "status": "ok",
        "method": "multi-split feature splitting; split p-values merged by "
                  "Gasparin, Wang & Ramdas (2025) PNAS 122:e2410849122 Thm 4.1",
        "merge_rule": merge,
        "n_splits": int(n_splits),
        "n_splits_usable": int(n_used),
        **_block_report(A, key, units),
        "frac_discovery": float(frac_discovery),
        "k_tested": int(k),
        "p_value_structure": merged_global,
        "p_value_structure_order_sensitivity": merge_order_sensitivity(
            split_p, merge=merge),
        "structure_detected": bool(merged_global < 0.05),
        "what_this_tests": "H0: held-out features are independent of the labels "
                           "learned from the discovery features. Rejection means "
                           "detectable structure, NOT necessarily discrete "
                           "archetypes - see 'discreteness'",
        "p_per_split": [float(v) for v in split_p],
        "p_split_median": float(np.median(arr)),
        "p_split_min": float(arr.min()),
        "p_split_max": float(arr.max()),
        "per_feature": table.to_dict("records"),
        "n_defining": int(reject.sum()),
        # Does the test concern the partition the pipeline reports? Each split
        # re-derives labels from its own discovery block; if those labels do not
        # resemble the reported ones, the p-value is about detectable structure
        # in general, not about this partition.
        "ari_discovery_vs_reported": (
            {"median": float(np.median(disc_ari)),
             "min": float(np.min(disc_ari)),
             "max": float(np.max(disc_ari)),
             "n": len(disc_ari)} if disc_ari else None),
        "concerns_reported_partition": (
            None if not disc_ari else bool(np.median(disc_ari) >= 0.5)),
    }


# --------------------------------------------------------------------------- #
# 4. Is the structure discrete, or a continuum?
# --------------------------------------------------------------------------- #
def discreteness_evidence(X, k: int, *, rng: Optional[np.random.Generator] = None,
                          n_init: int = 8) -> Dict[str, object]:
    """**This does not test discreteness.** It compares a mixture against
    independence, and calls a smooth continuum "strongly discrete" every time.

    Use :func:`continuum_null_test` for the discreteness question. This function
    is kept only because the comparison it makes -- ``BIC(1 component)`` against
    ``BIC(k components)`` -- is a useful descriptive statistic and because §4 of
    the manuscript reports its failure as a result. On 20 one-factor continuum
    datasets (n = 200, p = 30, no groups at all) it returns verdict "strong"
    20/20 with median ``delta_bic`` 586, and on two-factor continua 20/20 with
    median 727; it only says "against" on independent Bernoulli noise. A product
    of independent Bernoullis is a bad model for *any* correlated data, so a
    mixture beats it on a gradient as easily as on archetypes.

    The mechanics, for the descriptive use: comparing a ``k``-component
    Bernoulli mixture against a single component by BIC, on the same binary
    matrix.

    ``delta_bic = BIC(1) - BIC(k)`` is positive when the discrete model is
    preferred. Following the usual reading of the Bayes-factor scale
    (Kass & Raftery 1995, JASA 90:773-795, doi:10.1080/01621459.1995.10476572),
    values above 10 are strong evidence and below 2 are not worth mentioning.

    This is a model comparison, not a test: a mixture of Bernoullis can also
    approximate a continuum with enough components, so a positive ``delta_bic``
    supports discreteness only relative to the specific one-component
    alternative. Reported next to the split p-value so the two questions are not
    conflated.
    """
    from .baselines import bernoulli_mixture

    rng = rng or np.random.default_rng()
    A = np.asarray(X, dtype=float)
    if k <= 1 or A.shape[0] < 4 * k:
        return {"status": "skipped",
                "reason": f"needs k >= 2 and n >= {4 * max(k, 2)}"}
    one = bernoulli_mixture(A, 1, rng=rng, n_init=1)
    many = bernoulli_mixture(A, int(k), rng=rng, n_init=n_init)
    delta = float(one["bic"] - many["bic"])
    if delta > 10:
        verdict = "strong"
    elif delta > 2:
        verdict = "positive"
    elif delta > -2:
        verdict = "inconclusive"
    else:
        verdict = "against"
    return {
        "status": "ok",
        "k": int(k),
        "bic_one_component": one["bic"],
        "bic_k_components": many["bic"],
        "delta_bic": delta,
        "delta_bic_per_isolate": delta / A.shape[0],
        "verdict": verdict,
        "note": "delta_bic = BIC(1) - BIC(k); > 0 favours discrete groups over a "
                "single population. A significant split p-value with delta_bic "
                "<= 0 means structure without evidence of discreteness, i.e. a "
                "gradient rather than archetypes",
    }


# --------------------------------------------------------------------------- #
# 5. Discrete groups, or one gradient? A dependence-preserving null
# --------------------------------------------------------------------------- #
def fit_one_factor(X, *, n_iter: int = 30, ridge: float = 1.0,
                   max_abs: float = 6.0):
    """One-factor logistic latent trait: ``fit_latent_trait(X, q=1)``.

    ``P(x_ij = 1) = sigmoid(a_j + b_j z_i)`` with ``z_i`` a single continuous
    latent score. This is the canonical *unimodal* alternative to a mixture: one
    population, one gradient, features correlated through it. Kept as a named
    entry point because it is the ``q = 1`` special case the continuum test used
    exclusively before multi-dimensional continua were calibrated.

    ``ridge`` and ``max_abs`` are not cosmetic. Binary features are frequently
    near-separable along the latent score, and an unpenalised fit then sends the
    slope to infinity: bootstrap samples drawn from such a fit have
    astronomically strong dependence and the null distribution of any downstream
    statistic becomes useless (we measured a null standard deviation of 2e6
    before adding them). Coefficients are ridge-penalised and clipped to
    ``[-max_abs, max_abs]``, which corresponds to probabilities within about
    ``[0.002, 0.998]``.

    Returns ``(a, b, z)`` with ``b`` and ``z`` one-dimensional.
    """
    a, B, Z = fit_latent_trait(X, 1, n_iter=n_iter, ridge=ridge,
                               max_abs=max_abs)
    return a, B[:, 0], Z[:, 0]


def _quadrature(q: int, n_nodes: Optional[int] = None):
    """Product Gauss-Hermite grid over ``R^q`` against a standard normal.

    The node count per dimension shrinks with ``q`` so the grid stays around a
    thousand points: ``q = 1`` keeps 21 nodes, which reproduces the one-factor
    quadrature exactly, and ``q = 3`` uses 9 (729 points). Deterministic, so the
    observed statistic and every bootstrap replicate are integrated the same
    way and their difference carries no Monte-Carlo noise.
    """
    if n_nodes is None:
        n_nodes = {1: 21, 2: 15, 3: 9}.get(q, 7)
    nodes, weights = np.polynomial.hermite_e.hermegauss(int(n_nodes))
    weights = weights / weights.sum()
    grids = np.meshgrid(*([nodes] * q), indexing="ij")
    pts = np.stack([g.ravel() for g in grids], axis=1)       # (G, q)
    wg = np.ones(pts.shape[0])
    for d in range(q):
        wg *= weights[np.indices([len(nodes)] * q)[d].ravel()]
    return pts, wg / wg.sum()


def fit_latent_trait(X, q: int = 1, *, n_iter: int = 30, ridge: float = 1.0,
                     max_abs: float = 6.0):
    """Fit a ``q``-factor logistic latent-trait model to a binary matrix.

    ``P(x_ij = 1) = sigmoid(a_j + b_j . z_i)`` with ``z_i`` in R^q. ``q = 1`` is
    :func:`fit_one_factor`.

    Why more than one factor is not optional. The one-factor version of this
    model is the null of :func:`continuum_null_test`, and the *only* continuum
    it can represent is a single ordering of isolates. An accessory genome is
    not that: resistance load and virulence load are separate gradients, and the
    *Klebsiella* cohort's own layer reports give effective dimensions 5.2 and
    2.5. Against cluster-free data generated from a two- or three-factor
    continuum, a one-factor null rejects 95-100 % of the time (calibration
    study, experiment F) - i.e. it calls a multi-dimensional gradient
    "discrete". The fix is to let the null have as many dimensions as the data
    ask for.

    Fitted by **marginal** maximum likelihood (Bock & Aitkin 1981,
    *Psychometrika* 46:443-459, doi:10.1007/BF02293801): an EM whose E-step is
    the posterior over a fixed Gauss-Hermite grid and whose M-step is a
    quadrature-weighted logistic regression per feature. This is not a stylistic
    choice. Fitting ``(a, B)`` and ``Z`` by alternating Newton steps, which is
    the obvious alternative, maximises the *joint* likelihood
    with ``n q`` incidental parameters, which inflates the loadings badly: on
    data simulated with mean ``|b| = 0.98`` the joint fit returned 1.99, and the
    parametric bootstrap then drew replicates with mean absolute feature
    correlation 0.25 against the observed 0.16. A null model more strongly
    dependent than the data it is calibrating is not a null model, and the
    resulting continuum test rejected on cluster-free two-factor data every
    time. Fitting the marginal likelihood instead removes that bias: on the same
    simulation the loadings come back at 1.08 and the bootstrap replicates at
    correlation 0.18 against the observed 0.16, and the false-positive rate on
    two-factor continuum data falls from 1.00 to 0.08.

    One caveat on the word "EM": the M-step here is two ridge-penalised, clipped
    Newton steps rather than an exact maximisation of the Q-function, so the
    marginal likelihood is not guaranteed to increase monotonically and
    ``n_iter`` is a fixed budget rather than a convergence criterion. What is
    claimed for it is the measured calibration above, not a monotonicity
    theorem.

    ``ridge`` and ``max_abs`` control the near-separation blow-up documented in
    :func:`fit_one_factor`; they are much less load-bearing under marginal ML
    than under the joint fit.

    Returns ``(a, B, Z)`` with ``a`` of shape ``(p,)``, ``B`` ``(p, q)`` and
    ``Z`` ``(n, q)`` the posterior mean scores (reported for convenience; they
    are not parameters of the model).
    """
    A = np.asarray(X, dtype=float)
    n, p = A.shape
    q = int(max(1, q))
    pts, wg = _quadrature(q)                          # (G, q), (G,)
    Z1 = np.hstack([np.ones((pts.shape[0], 1)), pts])  # (G, q+1)
    logw = np.log(np.clip(wg, 1e-300, None))

    # start from the leading q principal components, rescaled to unit variance
    C = A - A.mean(axis=0, keepdims=True)
    try:
        _, _, Vt = np.linalg.svd(C, full_matrices=False)
        Z0 = C @ Vt[:q].T
    except np.linalg.LinAlgError:
        Z0 = np.tile(C.sum(axis=1)[:, None], (1, q))
    if Z0.shape[1] < q:
        Z0 = np.hstack([Z0, np.zeros((n, q - Z0.shape[1]))])
    sd = np.where(Z0.std(axis=0) > 0, Z0.std(axis=0), 1.0)
    Z0 = (Z0 - Z0.mean(axis=0)) / sd
    a = np.clip(np.log(np.clip(A.mean(axis=0), 1e-3, 1 - 1e-3)
                       / np.clip(1 - A.mean(axis=0), 1e-3, 1)), -max_abs, max_abs)
    B = np.clip(np.linalg.lstsq(Z0, A - A.mean(axis=0), rcond=None)[0].T * 4.0,
                -max_abs, max_abs)
    eye = np.eye(q + 1) * ridge

    # The M-step is clipped, ridge-penalised Newton rather than an exact
    # maximisation of the Q-function, so the marginal likelihood is NOT
    # guaranteed to increase -- measured on the case-study matrix it decreases
    # on 25 of 29 iterations at q = 2, peaking around iteration 5. Keeping the
    # best iterate rather than the last makes the returned fit insensitive to
    # `n_iter`, which is what `select_latent_dimension` needs in order to
    # compare across q at all.
    R = None
    best = (-np.inf, a.copy(), B.copy())
    for _ in range(int(n_iter)):
        # E step: posterior over the grid
        eta = a[None, :] + pts @ B.T                  # (G, p)
        lp = -np.logaddexp(0.0, -eta)
        lq = -np.logaddexp(0.0, eta)
        acc = A @ lp.T + (1.0 - A) @ lq.T + logw[None, :]      # (n, G)
        acc -= acc.max(axis=1, keepdims=True)
        R = np.exp(acc)
        R /= R.sum(axis=1, keepdims=True)
        # M step: quadrature-weighted logistic regression, two Newton steps
        Ng = R.sum(axis=0)                            # (G,)
        Sg = R.T @ A                                  # (G, p) expected successes
        for _ in range(2):
            mu = 1.0 / (1.0 + np.exp(-np.clip(a[None, :] + pts @ B.T, -30, 30)))
            w = np.clip(mu * (1 - mu), 1e-9, None) * Ng[:, None]     # (G, p)
            H = np.einsum("gp,ga,gb->pab", w, Z1, Z1) + eye
            g = np.einsum("gp,ga->pa", Sg - Ng[:, None] * mu, Z1)
            step = np.linalg.solve(H, g[..., None])[..., 0]
            a = np.clip(a + step[:, 0], -max_abs, max_abs)
            B = np.clip(B + step[:, 1:], -max_abs, max_abs)
        ll = latent_trait_loglik(A, a, B)
        if ll > best[0]:
            best = (ll, a.copy(), B.copy())
    a, B = best[1], best[2]
    eta = a[None, :] + pts @ B.T
    acc = (A @ (-np.logaddexp(0.0, -eta)).T
           + (1.0 - A) @ (-np.logaddexp(0.0, eta)).T + logw[None, :])
    acc -= acc.max(axis=1, keepdims=True)
    R = np.exp(acc)
    R /= R.sum(axis=1, keepdims=True)
    return a, B, R @ pts


def latent_trait_loglik(X, a, B, *, n_nodes: Optional[int] = None) -> float:
    """Marginal log-likelihood of the ``q``-factor model, integrating out ``z``.

    ``log P(x_i) = log int prod_j Bern(x_ij | sigmoid(a_j + b_j . z)) N(z; 0, I) dz``
    by a product Gauss-Hermite grid. Integrating rather than plugging in the
    fitted scores is what makes the comparison with a mixture fair: the latent
    scores are not free parameters, so the model costs ``p (1 + q)`` parameters
    rather than ``p (1 + q) + n q``.
    """
    A = np.asarray(X, dtype=float)
    B = np.atleast_2d(np.asarray(B, dtype=float))
    if B.shape[0] != A.shape[1]:
        B = B.T
    q = B.shape[1]
    pts, wg = _quadrature(q, n_nodes)
    logw = np.log(np.clip(wg, 1e-300, None))
    acc = np.empty((A.shape[0], pts.shape[0]))
    for g in range(pts.shape[0]):
        eta = a[None, :] + (B @ pts[g])[None, :]
        lp = -np.logaddexp(0.0, -eta)
        lq = -np.logaddexp(0.0, eta)
        acc[:, g] = (A * lp + (1 - A) * lq).sum(axis=1) + logw[g]
    m = acc.max(axis=1, keepdims=True)
    return float((m[:, 0] + np.log(np.exp(acc - m).sum(axis=1))).sum())


def select_latent_dimension(X, *, q_max: int = 4, n_iter: int = 30
                            ) -> Dict[str, object]:
    """Choose the number of continuum dimensions by BIC.

    The routine evaluates ``q = 1..3`` first. It evaluates the supported fourth
    dimension only when ``q = 3`` is the provisional boundary optimum. This is
    adaptive in the literal gate-relevant sense: if BIC has already turned up,
    the optimum is interior and a costly fourth-dimensional product quadrature
    cannot change the boundary verdict; if it has not turned up, ``q = 4`` is
    required. ``q_max`` is the maximum permitted dimension, not a requirement
    to fit every dimension unconditionally.

    Returns the selected ``q``, every BIC actually evaluated, whether the
    adaptive extension was needed, and whether the curve was still falling at
    the final permitted dimension (in which case the continuum is at least that
    rich and the selected value is a floor, not an optimum).
    """
    A = np.asarray(X, dtype=float)
    n, p = A.shape
    logn = math.log(n)
    bics, bics_matched, fits = {}, {}, {}

    q_limit = int(q_max)
    if q_limit < 1:
        raise ValueError("q_max must be >= 1")
    initial_max = min(q_limit, 3)

    def fit_dimension(q: int) -> None:
        a, B, _ = fit_latent_trait(A, q, n_iter=n_iter)
        fits[q] = (a, B)
        bics[q] = float(-2 * latent_trait_loglik(A, a, B) + p * (1 + q) * logn)

    for q in range(1, initial_max + 1):
        fit_dimension(q)
    provisional_best = min(bics, key=bics.get)
    adaptive_extension_used = bool(
        provisional_best == initial_max and initial_max < q_limit
    )
    if adaptive_extension_used:
        for q in range(initial_max + 1, q_limit + 1):
            fit_dimension(q)
    # The default grid uses fewer nodes per dimension as q grows, so the
    # log-likelihood is least accurate exactly where the model is richest --
    # which biases the curve AGAINST large q, not for it. The matched-grid
    # column uses the same node count at every q so a reader can check that the
    # selected q is not an artefact of the integration.
    matched = min(_MATCHED_NODES, 21)
    for q, (a, B) in fits.items():
        bics_matched[q] = float(-2 * latent_trait_loglik(A, a, B, n_nodes=matched)
                                + p * (1 + q) * logn)
    best = min(bics, key=bics.get)
    best_matched = min(bics_matched, key=bics_matched.get)
    return {"q_selected": int(best),
            "bic_by_q": {str(k): v for k, v in bics.items()},
            "bic_by_q_matched_grid": {str(k): v for k, v in bics_matched.items()},
            "matched_grid_nodes_per_dim": int(matched),
            "q_selected_matched_grid": int(best_matched),
            "q_max": q_limit,
            "q_evaluated_max": int(max(fits)),
            "adaptive_extension_used": adaptive_extension_used,
            "at_boundary": bool(best == q_limit and q_limit > 1)}


def one_factor_loglik(X, a, b, *, n_nodes: int = 21) -> float:
    """Marginal log-likelihood of the one-factor model: ``latent_trait_loglik``
    with a single column of loadings."""
    return latent_trait_loglik(X, a, np.asarray(b, dtype=float)[:, None],
                               n_nodes=n_nodes)


def continuum_null_test(X, k: int, *, rng: Optional[np.random.Generator] = None,
                        n_boot: int = 99, n_init: int = 5,
                        n_nodes: Optional[int] = None,
                        q: Optional[int] = None, q_max: int = 4
                        ) -> Dict[str, object]:
    """Test discrete groups against a **unimodal continuum** null.

    A significant feature-split p-value says there is structure; it says so just
    as loudly for one continuous gradient as for discrete archetypes, because
    both make the held-out features depend on the discovery labels. This
    function answers the other question.

    The statistic is ``BIC(q-factor) - BIC(k-component mixture)``, positive
    when a mixture describes the data better than a continuum does. Both models
    are fitted to the same binary matrix; the latent-trait likelihood is
    *marginal*, integrating the latent score out by Gauss-Hermite quadrature, so
    it is charged ``2p`` parameters rather than ``2p + n`` and the comparison is
    fair.

    Comparing instead against a **single independent-Bernoulli component**,
    which is the obvious thing to do, does not work: a product of independent Bernoullis is a bad model for any
    correlated data, so a mixture beats it on a smooth gradient as easily as on
    discrete groups, and the comparison calls everything discrete.

    Calibration comes from a parametric bootstrap: simulate ``n_boot`` datasets
    from the *fitted latent-trait model* - one population, no groups, features
    correlated through ``q`` continuous scores - and compare the observed
    statistic with what the continuum itself produces. The null therefore
    preserves exactly the feature-feature dependence that breaks the naive
    tests. The p-value is Phipson-Smyth corrected, so its floor is
    ``1/(n_boot + 1)``, which is reported. The observed number of bootstrap
    exceedances and its exact 95 % Clopper-Pearson interval are reported too.
    This separates a resolvable decision at alpha=0.05 from the misleading
    impression that a floor value is a precisely estimated tail probability.

    **How many dimensions the continuum gets is not a free choice.** A one-factor
    null can only represent a single ordering of the isolates, and against
    cluster-free data drawn from a two- or three-factor continuum it rejects
    95-100 % of the time (calibration study, experiment F): it calls a
    multi-dimensional gradient "discrete", which is the one error this test
    exists to prevent. ``q`` is therefore selected by BIC over ``1..q_max`` on
    the observed matrix (:func:`select_latent_dimension`) unless it is given
    explicitly. Dimensions 1--3 are fitted first; the supported fourth
    dimension is fitted only when q=3 is the provisional boundary optimum.
    The selected value, every fitted BIC and whether the optimum remains at
    ``q_max`` are reported. A boundary optimum yields the explicit third state
    ``withheld_under_dimensioned`` before bootstrapping; it is never serialized
    as a negative result or read as evidence of discreteness.
    """
    from .baselines import bernoulli_mixture
    from .stats import permutation_pvalue

    rng = rng or np.random.default_rng()
    A = np.asarray(X, dtype=float)
    n, p = A.shape
    if k <= 1 or n < 4 * k or p < 3:
        return {"status": "skipped",
                "reason": f"needs k >= 2, n >= {4 * max(k, 2)} and p >= 3"}

    logn = math.log(n)
    dim = (select_latent_dimension(A, q_max=int(q_max)) if q is None else
           {"q_selected": int(q), "bic_by_q": {}, "q_max": int(q_max),
            "at_boundary": False})
    q_use = int(dim["q_selected"])

    # A boundary optimum means that the null is known to be too small. Running
    # the bootstrap anyway produces an attractive but scientifically unusable
    # p-value and used to serialize the withheld verdict as ``false``. Stop
    # before simulation and expose an actual third state instead.
    if dim["at_boundary"]:
        return {
            "status": "withheld_under_dimensioned",
            "k": int(k),
            "statistic": None,
            "observed": None,
            "null_mean": None,
            "null_sd": None,
            "p_value": None,
            "p_value_floor": float(1.0 / (int(n_boot) + 1)),
            "n_boot": 0,
            "n_boot_requested": int(n_boot),
            "latent_dimension": dim,
            "discrete_beyond_a_gradient": None,
            "verdict": "withheld_under_dimensioned",
            "null_model": None,
            "note": "BIC selected the largest permitted latent dimension; "
                    "increase tva.continuum_q_max within the supported range "
                    "or report the discreteness verdict as not estimable.",
        }

    def stat(M):
        aa, BB, _ = fit_latent_trait(M, q_use)
        ll1 = latent_trait_loglik(M, aa, BB, n_nodes=n_nodes)
        bic1 = -2 * ll1 + p * (1 + q_use) * logn
        mix = bernoulli_mixture(M, int(k), rng=rng, n_init=n_init)
        return float(bic1 - mix["bic"])

    obs = stat(A)
    a, B, _ = fit_latent_trait(A, q_use)
    null = np.empty(int(n_boot))
    for t in range(int(n_boot)):
        Z = rng.standard_normal((n, q_use))
        eta = a[None, :] + Z @ B.T
        P = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
        null[t] = stat((rng.random((n, p)) < P).astype(float))
    exceedances = int(np.sum(null >= obs))
    pval = permutation_pvalue(null, obs, tail="greater")
    # Exact binomial uncertainty for the latent-null tail probability. This is
    # intentionally attached to the raw exceedance count, not to the
    # Phipson-Smyth corrected Monte Carlo p-value.
    alpha_ci = 0.05
    ci_low = (0.0 if exceedances == 0 else
              float(beta.ppf(alpha_ci / 2, exceedances,
                             int(n_boot) - exceedances + 1)))
    ci_high = (1.0 if exceedances == int(n_boot) else
               float(beta.ppf(1 - alpha_ci / 2, exceedances + 1,
                              int(n_boot) - exceedances)))
    return {
        "status": "ok",
        "k": int(k),
        "statistic": f"BIC({q_use}-factor latent trait) - BIC(k-component mixture)",
        "observed": obs,
        "null_mean": float(null.mean()),
        "null_sd": float(null.std(ddof=1)) if null.size > 1 else 0.0,
        "p_value": pval,
        "p_value_floor": float(1.0 / (int(n_boot) + 1)),
        "n_boot": int(n_boot),
        "bootstrap_exceedances": exceedances,
        "bootstrap_tail_probability_ci95": {
            "method": "Clopper-Pearson exact binomial",
            "confidence_level": 0.95,
            "low": ci_low,
            "high": ci_high,
        },
        "resolved_at_alpha_0_05": bool(ci_high < 0.05 or ci_low > 0.05),
        "latent_dimension": dim,
        "discrete_beyond_a_gradient": bool(pval < 0.05),
        "verdict": "discrete" if pval < 0.05 else "continuum_compatible",
        "null_model": f"{q_use}-factor logistic latent trait (one population, "
                      f"{q_use} continuous gradients, features correlated "
                      "through them), marginal likelihood by Gauss-Hermite "
                      "quadrature; q chosen by BIC over 1.." + str(q_max),
        "note": "a significant feature-split p-value with a non-significant "
                "value here means structure without evidence of discreteness. "
                "A boundary optimum is returned before bootstrap with status "
                "withheld_under_dimensioned and a null verdict.",
    }
