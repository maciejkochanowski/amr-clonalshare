"""tva.py — Thinning-Validated Archetypes (TVA), multi-split edition.

Post-clustering inference for multi-layer trait archetypes. Clustering the data
and then testing, on the same data, which features "define" the clusters is
double dipping: under the null of no clusters the resulting p-values are not
p-values at all (Gao, Bien & Witten 2024, JASA 119:332-342,
doi:10.1080/01621459.2022.2116331; Chen & Witten 2023, JMLR 24(152)).

Data thinning (Neufeld, Dharamshi, Gao & Witten 2024, JMLR 25(57):1-35) splits
each observation into two independent parts drawn from the same family, so that
clustering may use one part and testing the other. For a negative binomial
count with **known** size parameter ``r`` (their Table 2):

    X1 | X = x  ~  BetaBinomial(x, eps*r, (1-eps)*r),    X2 = X - X1
    =>  X1 ~ NB(eps*r, p)  independent of  X2 ~ NB((1-eps)*r, p).

Three things go wrong when this recipe is applied naively to bacterial trait
panels, and this module fixes all three.

1. **The verdict is a coin flip.** One thinning draw gives one p-value, and
   that p-value is a random variable even with the data held fixed. On the
   *Klebsiella* example shipped with this package, 40 draws of the previous
   single-split implementation produced p between 3.8e-63 and 0.97 (median
   0.43). A tool that answers "are the archetypes real?" differently depending
   on an internal seed is not usable.

   *Fix:* draw ``n_splits`` independent thinnings and merge the resulting
   p-values with a rule that is valid for **exchangeable** p-values. Repeated
   splits of one fixed dataset are exchangeable but arbitrarily dependent, so
   ordinary Fisher/Stouffer combination is invalid. We use the exchangeable
   Ruger improvement of Gasparin, Wang & Ramdas (2025, PNAS 122(11):e2410849122,
   Theorem 4.1),

       F_ER(p) = (K/k) * min_{l=1..K} p^{(l)}_{(ceil(l*k/K))}

   with ``k = ceil(K/2)`` (their recommended "twice the median" variant for
   strongly dependent p-values). ``F_ER`` strictly dominates the classical
   Ruger bound ``(K/k) p_(k)`` and requires no independence. The classical
   Ruger rule and the sharp Vovk-Wang ``2 * mean`` rule (Vovk & Wang 2020,
   Biometrika 107:791-808) are available via ``merge=``.

2. **The per-feature test double-dipped inside X2.** The previous
   implementation chose, for each feature, the cluster with the largest mean of
   ``X2[:, j]`` and then ran a one-sided Mann-Whitney test in that direction on
   the same ``X2[:, j]``. Under a global null with no clusters and ``r`` known,
   the marginal rejection rate of that test measured here is 0.145 at nominal
   0.05, and the probability of at least one BH-"defining" feature is 0.095 at
   nominal 0.05.

   *Fix:* a direction-free k-sample Kruskal-Wallis test of ``X2[:, j]`` across
   the labels learned from ``X1``. Nothing about the test depends on ``X2``
   beyond the values being tested. Measured under the same null: marginal
   0.044, FWER 0.010.

3. **The negative binomial model was wrong for most columns.** ``r`` was
   estimated as one pooled scalar from the very data being split. Neufeld et
   al. (2024) Proposition 11 gives the exact damage: with an assumed ``r~`` and
   a true ``r``,

       cov(X1, X2) = (1-eps) * r * ((1-p)/p)^2 * (1 - (r+1)/(r~+1)),

   which is non-zero whenever ``r~ != r``: negative when ``r~ < r``, positive
   when ``r~ > r``. On the
   *Klebsiella* example the shipped code produced a mean within-column
   corr(X1, X2) of -0.174, against +0.003 for a genuine NB matrix with known
   ``r`` - i.e. the independence the method rests on did not hold.

   *Fix:* ``screen_thinnable`` refuses to thin columns that are not plausibly
   negative binomial - under-dispersed columns (``var <= mean``), columns with
   fewer than ``min_distinct`` distinct values (a rescaled indicator is not a
   count), constant columns, and columns that are exact linear combinations of
   columns already retained. Dispersion is estimated **per feature**, not
   pooled. Everything that is dropped is reported, with the reason, in
   ``screening``.

Small cohorts. Thinning removes information: at ``n`` of order 100 the split
halves are too small for the discovery step to find anything, and the test has
essentially no power. ``amr_clonalshare.archephy`` provides a
phylogeny-aware alternative whose null is closed-form and therefore calibrated
at small ``n`` by construction.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.special import gammaln
from scipy.stats import kruskal, mannwhitneyu
from sklearn.cluster import KMeans

from .stats import benjamini_hochberg

__all__ = [
    "estimate_dispersion_mom",
    "estimate_dispersion_mle",
    "estimate_dispersion_per_feature",
    "nb_thin",
    "poisson_thin",
    "screen_thinnable",
    "merge_exchangeable",
    "ruger_merge",
    "twice_mean_merge",
    "tva_test_separation",
    "tva_defining_features",
    "tva_report",
    "binary_layer_to_counts",
    "auto_prefix_groups",
]


# --------------------------------------------------------------------------- #
# 1. p-value merging for exchangeable, arbitrarily dependent p-values
# --------------------------------------------------------------------------- #
def ruger_merge(pvals: Sequence[float], k: Optional[int] = None) -> float:
    """Classical Ruger bound ``(K/k) * p_(k)``.

    Valid for arbitrarily dependent p-values (Ruger 1978; see Vovk & Wang 2020,
    Biometrika 107:791-808, doi:10.1093/biomet/asaa027, Section 2).
    """
    p = np.sort(np.asarray(list(pvals), dtype=float))
    K = p.size
    if K == 0:
        return 1.0
    if k is None:
        k = int(math.ceil(K / 2))
    k = int(np.clip(k, 1, K))
    return float(min(1.0, (K / k) * p[k - 1]))


def twice_mean_merge(pvals: Sequence[float]) -> float:
    """``2 * mean(p)``. Sharp under arbitrary dependence (Vovk & Wang 2020)."""
    p = np.asarray(list(pvals), dtype=float)
    if p.size == 0:
        return 1.0
    return float(min(1.0, 2.0 * p.mean()))


def merge_exchangeable(pvals: Sequence[float], k: Optional[int] = None) -> float:
    """Exchangeable Ruger improvement (Gasparin, Wang & Ramdas 2025, Thm 4.1).

    ``F_ER(p) = (K/k) * min_{l=1..K} p^{(l)}_{(ceil(l*k/K))}`` where
    ``p^{(l)}_{(j)}`` is the j-th smallest of the first ``l`` p-values. Requires
    the p-values to be exchangeable (they are: repeated random splits of one
    fixed dataset) but allows arbitrary dependence. Strictly dominates
    :func:`ruger_merge` at the same ``k``.
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
    if k is None:
        k = int(math.ceil(K / 2))
    k = int(np.clip(k, 1, K))
    best = np.inf
    for l in range(1, K + 1):
        j = int(np.clip(math.ceil(l * k / K), 1, l))
        best = min(best, float(np.sort(p[:l])[j - 1]))
    return float(np.clip((K / k) * best, 0.0, 1.0))


_MERGERS = {
    "exchangeable_ruger": merge_exchangeable,
    "ruger": ruger_merge,
    "twice_mean": twice_mean_merge,
}


def _merge(pvals: "Sequence[float] | np.ndarray", how: str) -> float:
    try:
        fn = _MERGERS[how]
    except KeyError:
        raise ValueError(
            f"unknown merge rule {how!r}; choose from {sorted(_MERGERS)}"
        ) from None
    return fn(list(pvals))


# --------------------------------------------------------------------------- #
# 2. Dispersion estimation
# --------------------------------------------------------------------------- #
def estimate_dispersion_mom(X: np.ndarray, *, floor: float = 1e-3,
                            ceil: float = 1e4) -> float:
    """Pooled method-of-moments dispersion (median over over-dispersed columns).

    Kept for backwards compatibility and for the ``dispersion="pooled"`` option.
    ``estimate_dispersion_per_feature`` is the default and is preferred: a
    single pooled ``r`` is a misspecification for every column whose true
    dispersion differs, and Neufeld et al. (2024) Prop. 11 quantifies the
    resulting dependence between the two halves.
    """
    X = np.asarray(X, dtype=float)
    if X.size == 0:
        return ceil
    mean = X.mean(axis=0)
    var = X.var(axis=0, ddof=1) if X.shape[0] > 1 else np.zeros(X.shape[1])
    over = var > mean + 1e-9
    if not over.any():
        return ceil
    with np.errstate(divide="ignore", invalid="ignore"):
        r_per = mean ** 2 / (var - mean)
    return float(np.clip(np.median(r_per[over]), floor, ceil))


def _nb_neg_log_lik(r: float, X: np.ndarray) -> float:
    if r <= 0:
        return np.inf
    means = X.mean(axis=0)
    p = r / (r + np.maximum(means, 1e-9))
    term1 = gammaln(X + r) - gammaln(r) - gammaln(X + 1.0)
    term2 = r * np.log(p)[None, :] + X * np.log1p(-p)[None, :]
    return float(-(term1 + term2).sum())


def estimate_dispersion_mle(X: np.ndarray, *, floor: float = 1e-3,
                            ceil: float = 1e4) -> float:
    """Pooled profile-likelihood dispersion. Falls back to MoM on failure."""
    X = np.asarray(X, dtype=float)
    if X.size == 0 or X.shape[0] < 2:
        return ceil
    try:
        res = minimize_scalar(
            lambda lr: _nb_neg_log_lik(math.exp(lr), X),
            bounds=(math.log(floor), math.log(ceil)),
            method="bounded", options={"xatol": 1e-4},
        )
        if res.success:
            return float(np.clip(math.exp(res.x), floor, ceil))
    except Exception:
        pass
    return estimate_dispersion_mom(X, floor=floor, ceil=ceil)


def estimate_dispersion_per_feature(X: np.ndarray, *, method: str = "mom",
                                    floor: float = 1e-3,
                                    ceil: float = 1e4) -> np.ndarray:
    """Per-feature dispersion vector ``r_j`` (the default for thinning)."""
    X = np.asarray(X, dtype=float)
    n, p = X.shape
    if method == "mle":
        return np.array([estimate_dispersion_mle(X[:, j:j + 1], floor=floor, ceil=ceil)
                         for j in range(p)], dtype=float)
    mean = X.mean(axis=0)
    var = X.var(axis=0, ddof=1) if n > 1 else np.zeros(p)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(var > mean, mean ** 2 / (var - mean), np.inf)
    return np.clip(np.nan_to_num(r, nan=ceil, posinf=ceil), floor, ceil)


# --------------------------------------------------------------------------- #
# 3. Thinnability screening (condition of validity)
# --------------------------------------------------------------------------- #
def screen_thinnable(X: pd.DataFrame, *, min_distinct: int = 3,
                     rank_tol: float = 1e-9) -> Dict[str, Any]:
    """Decide which count columns may be NB-thinned, and say why not otherwise.

    A column is retained only if all of the following hold:

    * ``var > mean``  - NB requires over-dispersion; ``var <= mean`` means the
      column is Poisson-like or binomial-like and the NB split is misspecified
      (Neufeld et al. 2024 Prop. 11 then makes X1 and X2 dependent).
    * ``>= min_distinct`` distinct values - a column taking two values is an
      indicator, possibly rescaled, and an indicator is its own sufficient
      statistic: it is not thinnable in any non-trivial sense.
    * ``mean > 0`` - all-zero columns carry nothing.
    * it is not an exact linear combination of the columns already retained -
      otherwise the same information is tested repeatedly and the multiplicity
      correction counts hypotheses that do not exist.

    Returns ``{"keep": [names], "dropped": [{feature, reason, ...}],
    "n_kept": int, "n_dropped": int}``.
    """
    if not isinstance(X, pd.DataFrame):
        raise TypeError("screen_thinnable expects a DataFrame of counts")
    A = X.values.astype(float)
    n, p = A.shape
    mean = A.mean(axis=0)
    var = A.var(axis=0, ddof=1) if n > 1 else np.zeros(p)
    distinct = np.array([np.unique(A[:, j]).size for j in range(p)])

    keep_idx: List[int] = []
    dropped: List[Dict[str, Any]] = []
    basis = np.zeros((n, 0))

    for j, name in enumerate(X.columns):
        rec = {"feature": str(name), "mean": float(mean[j]),
               "var": float(var[j]), "n_distinct": int(distinct[j])}
        if mean[j] <= 0:
            dropped.append({**rec, "reason": "all-zero column"})
            continue
        if distinct[j] < min_distinct:
            dropped.append({**rec, "reason":
                            f"only {int(distinct[j])} distinct values "
                            f"(< {min_distinct}); indicator-like, not a count"})
            continue
        if var[j] <= mean[j]:
            dropped.append({**rec, "reason":
                            f"not over-dispersed (var/mean = "
                            f"{var[j] / mean[j]:.3f} <= 1); NB model invalid"})
            continue
        cand = np.column_stack([basis, A[:, j]]) if basis.shape[1] else A[:, j:j + 1]
        if np.linalg.matrix_rank(cand, tol=rank_tol) <= basis.shape[1]:
            dropped.append({**rec, "reason":
                            "exact linear combination of retained columns"})
            continue
        keep_idx.append(j)
        basis = cand

    return {
        "keep": [str(X.columns[j]) for j in keep_idx],
        "keep_idx": keep_idx,
        "dropped": dropped,
        "n_kept": len(keep_idx),
        "n_dropped": len(dropped),
    }


# --------------------------------------------------------------------------- #
# 4. Count splitting
# --------------------------------------------------------------------------- #
def poisson_thin(X: np.ndarray, *, eps: float = 0.5,
                 rng: Optional[np.random.Generator] = None
                 ) -> Tuple[np.ndarray, np.ndarray]:
    """Poisson count splitting: ``X1 ~ Binomial(X, eps)``, ``X2 = X - X1``.

    Valid only for genuinely Poisson counts. On over-dispersed data X1 and X2
    remain dependent and type-I control is lost.
    """
    if rng is None:
        rng = np.random.default_rng()
    Xi = np.asarray(X, dtype=np.int64)
    X1 = rng.binomial(Xi, eps)
    return X1, Xi - X1


def nb_thin(X: np.ndarray, *, r, eps: float = 0.5,
            rng: Optional[np.random.Generator] = None
            ) -> Tuple[np.ndarray, np.ndarray]:
    """Negative-binomial count splitting (Neufeld et al. 2024, Table 2).

    ``X1 | X = x ~ BetaBinomial(x, eps*r, (1-eps)*r)``, ``X2 = X - X1``,
    implemented as the beta-binomial compound (draw ``p ~ Beta``, then
    ``Binomial``).

    ``r`` may be a scalar (one dispersion for every column) or a length-``p``
    array (per-column dispersion, recommended).
    """
    if rng is None:
        rng = np.random.default_rng()
    if not (0.0 < eps < 1.0):
        raise ValueError(f"eps must lie in (0, 1); got {eps}")
    Xi = np.asarray(X, dtype=np.int64)
    r_arr = np.atleast_1d(np.asarray(r, dtype=float))
    if r_arr.size == 1:
        r_arr = np.full(Xi.shape[1], float(r_arr[0]))
    if r_arr.size != Xi.shape[1]:
        raise ValueError(
            f"r must be scalar or length {Xi.shape[1]}; got size {r_arr.size}")
    if np.any(r_arr <= 0):
        raise ValueError("all dispersions r must be positive")
    a = np.broadcast_to(eps * r_arr, Xi.shape)
    b = np.broadcast_to((1.0 - eps) * r_arr, Xi.shape)
    p_split = rng.beta(a, b)
    X1 = rng.binomial(Xi, p_split)
    return X1, Xi - X1


def split_informativeness(r, eps: float = 0.5) -> float:
    """Fraction of entries whose Beta split proportion falls in [0.05, 0.95].

    The beta-binomial split proportion is ``Beta(eps*r, (1-eps)*r)``, which is
    U-shaped for ``r < 2``: most entries then go almost entirely to one half and
    the split carries little information. Values near 0 mean the thinning is
    effectively all-or-nothing. Reported so a low-power result can be told apart
    from a genuine null.
    """
    r_arr = np.atleast_1d(np.asarray(r, dtype=float))
    from scipy.stats import beta as _beta
    vals = [float(_beta.cdf(0.95, eps * rj, (1 - eps) * rj)
                  - _beta.cdf(0.05, eps * rj, (1 - eps) * rj)) for rj in r_arr]
    return float(np.mean(vals))


# --------------------------------------------------------------------------- #
# 5. Tests on a single thinned pair
# --------------------------------------------------------------------------- #
def _labels_from_X1(X1: np.ndarray, k: int, seed: int, n_init: int
                    ) -> Optional[tuple]:
    if X1.shape[0] < 2 * k or X1.shape[1] < 1:
        return None
    try:
        km = KMeans(n_clusters=k, n_init=n_init, random_state=seed).fit(
            X1.astype(float))
    except Exception:
        return None
    return km, km.labels_


def _separation_pvalue(X1: np.ndarray, X2: np.ndarray, k: int, seed: int,
                       n_init: int) -> float:
    """Cluster X1; test X2 for separation along an X1-derived direction.

    The projection direction comes from the leading principal axis of the X1
    centroids, so it is a function of X1 only. The test is two-sided
    (Mann-Whitney for k=2, Kruskal-Wallis for k>2): no direction or grouping is
    ever selected using X2.
    """
    got = _labels_from_X1(X1, k, seed, n_init)
    if got is None:
        return 1.0
    km, labels = got
    uniq, cnt = np.unique(labels, return_counts=True)
    if uniq.size < 2 or cnt.min() < 2:
        return 1.0
    C = km.cluster_centers_
    Cc = C - C.mean(axis=0, keepdims=True)
    try:
        _, _, Vt = np.linalg.svd(Cc, full_matrices=False)
    except np.linalg.LinAlgError:
        return 1.0
    d = Vt[0]
    nrm = float(np.linalg.norm(d))
    if not np.isfinite(nrm) or nrm < 1e-12:
        return 1.0
    with np.errstate(over="ignore", invalid="ignore"):
        proj = X2.astype(float) @ (d / nrm)
    if not np.isfinite(proj).all():
        return 1.0
    groups = [proj[labels == c] for c in uniq if (labels == c).sum() >= 2]
    if len(groups) < 2:
        return 1.0
    if all(np.ptp(g) == 0 for g in groups) and np.ptp(np.concatenate(groups)) == 0:
        return 1.0
    try:
        if len(groups) == 2:
            return float(mannwhitneyu(groups[0], groups[1],
                                      alternative="two-sided").pvalue)
        return float(kruskal(*groups).pvalue)
    except ValueError:
        return 1.0


def _per_feature_pvalues(X1: np.ndarray, X2: np.ndarray, k: int, seed: int,
                         n_init: int) -> Tuple[np.ndarray, np.ndarray]:
    """Kruskal-Wallis of each X2 column across X1 labels; plus X1 group means.

    Returns ``(pvalues, top_cluster_from_X1)``. The reported "top" cluster is
    read off X1, never off X2, so it carries no selection effect into the test.
    """
    p = X2.shape[1]
    pv = np.ones(p, dtype=float)
    top = np.full(p, -1, dtype=int)
    got = _labels_from_X1(X1, k, seed, n_init)
    if got is None:
        return pv, top
    _, labels = got
    uniq = [int(c) for c in np.unique(labels) if (labels == c).sum() >= 2]
    if len(uniq) < 2:
        return pv, top
    for j in range(p):
        col2 = X2[:, j].astype(float)
        groups = [col2[labels == c] for c in uniq]
        if np.ptp(np.concatenate(groups)) == 0:
            continue
        try:
            pv[j] = float(kruskal(*groups).pvalue)
        except ValueError:
            pv[j] = 1.0
        means1 = {c: float(X1[labels == c, j].mean()) for c in uniq}
        top[j] = max(means1, key=lambda c: means1[c])
    return pv, top


# --------------------------------------------------------------------------- #
# 6. Public multi-split API
# --------------------------------------------------------------------------- #
#: Largest spread of the *estimated* per-feature dispersions for which a single pooled
#: value may still be used to thin. Set from the null calibration sweep in
#: ``audit_checks/calib_pooled_heterogeneity_sweep.py``: with a pooled estimate the
#: rejection rate on structure-free data is 0.033, 0.050, 0.250, 0.567, 0.800 and 0.950
#: at true dispersion ratios of 1, 2, 4, 8, 16 and 100 against a nominal 0.05. The last
#: safe level (true ratio 2) has a median estimated ratio of 2.25 and the first unsafe
#: level (true ratio 4) a median of 3.58, so the guard sits between them. Per-feature
#: estimation stayed at nominal across the whole sweep and is unaffected.
POOLED_DISPERSION_MAX_ESTIMATED_RATIO = 3.0

_DISPERSION_CHOICES = ("mom", "mle", "pooled_mom", "pooled_mle")


def pooled_dispersion_spread(X: np.ndarray, *, method: str = "mom") -> float:
    """Ratio of the largest to the smallest estimated per-feature dispersion."""
    per_feature = estimate_dispersion_per_feature(X, method=method)
    smallest = float(np.min(per_feature))
    if not np.isfinite(smallest) or smallest <= 0:
        return float("inf")
    return float(np.max(per_feature) / smallest)


def _resolve_r(X: np.ndarray, r, dispersion: str) -> np.ndarray:
    if r is not None:
        r_arr = np.atleast_1d(np.asarray(r, dtype=float))
        if r_arr.size == 1:
            r_arr = np.full(X.shape[1], float(r_arr[0]))
        return r_arr
    if dispersion not in _DISPERSION_CHOICES:
        # Silently falling back to the per-feature default would let a typo look like a
        # deliberate choice and would hide which estimator actually produced a result.
        raise ValueError(
            f"unknown dispersion {dispersion!r}; expected one of {_DISPERSION_CHOICES}"
        )
    if dispersion in ("pooled_mle", "pooled_mom"):
        method = "mle" if dispersion == "pooled_mle" else "mom"
        spread = pooled_dispersion_spread(X, method=method)
        if spread > POOLED_DISPERSION_MAX_ESTIMATED_RATIO:
            # Thinning is only valid when the r used for the split is the r that
            # generated the counts. One pooled value cannot be that for every feature
            # once they genuinely differ, and the leftover dependence is read as
            # cluster structure: on structure-free data the null rejection rate reaches
            # 0.95. Refusing is the only honest option; per-feature estimation works.
            raise ValueError(
                f"dispersion={dispersion!r} pools one value across features whose "
                f"estimated dispersions span a factor of {spread:.2f}, above the "
                f"calibrated limit of {POOLED_DISPERSION_MAX_ESTIMATED_RATIO}. A pooled "
                f"split on such data reports structure that is not present. Use "
                f"dispersion='mom' or 'mle', which estimate per feature."
            )
        estimator = estimate_dispersion_mle if dispersion == "pooled_mle" else estimate_dispersion_mom
        return np.full(X.shape[1], estimator(X))
    return estimate_dispersion_per_feature(
        X, method="mle" if dispersion == "mle" else "mom")


def tva_test_separation(X: np.ndarray, *, r=None, k: int = 2, eps: float = 0.5,
                        rng: Optional[np.random.Generator] = None,
                        n_splits: int = 9, n_init: int = 10,
                        dispersion: str = "mom",
                        merge: str = "exchangeable_ruger") -> Dict[str, Any]:
    """Multi-split thinning test of "there is real k-group structure".

    Returns a dict with the merged p-value plus the per-split p-values, so the
    randomisation spread is visible rather than hidden.
    """
    if rng is None:
        rng = np.random.default_rng()
    Xi = np.asarray(X, dtype=np.int64)
    # Thinning is valid only on columns that are plausibly negative binomial,
    # and ``tva_report`` screens before it tests. This entry point did not, so
    # a matrix of identical rows -- no variance to split at all -- came back at
    # p = 3e-14 from the very screen the module docstring credits with taking
    # the marginal rejection rate from 0.145 to 0.044.
    screening = screen_thinnable(
        pd.DataFrame(Xi, columns=[str(j) for j in range(Xi.shape[1])]))
    if screening["n_kept"] == 0:
        nan = float("nan")
        return {
            "p_value": nan,
            "status": "not_thinnable",
            "reason": f"none of {Xi.shape[1]} count columns are NB-thinnable "
                      f"(need >= 1); see 'screening.dropped' for the reason on "
                      f"each",
            "screening": screening,
            "merge_rule": merge,
            "n_splits": int(n_splits),
            "p_per_split": [],
            "p_split_median": nan,
            "p_split_min": nan,
            "p_split_max": nan,
            "k": int(k),
            "eps": float(eps),
            "r": [],
            "split_informativeness": nan,
        }
    keep_idx = np.asarray(screening["keep_idx"], dtype=int)
    r_for_kept = r
    if r is not None:
        supplied = np.atleast_1d(np.asarray(r, dtype=float))
        if supplied.size == Xi.shape[1]:
            r_for_kept = supplied[keep_idx]
    Xk = Xi[:, keep_idx]
    r_arr = _resolve_r(Xk, r_for_kept, dispersion)
    per_split = []
    for m in range(int(n_splits)):
        X1, X2 = nb_thin(Xk, r=r_arr, eps=eps, rng=rng)
        per_split.append(_separation_pvalue(X1, X2, k=k, seed=m, n_init=n_init))
    merged = _merge(per_split, merge)
    arr = np.asarray(per_split, dtype=float)
    return {
        "p_value": merged,
        "status": "ok",
        "screening": screening,
        "merge_rule": merge,
        "n_splits": int(n_splits),
        "p_per_split": [float(v) for v in per_split],
        "p_split_median": float(np.median(arr)),
        "p_split_min": float(arr.min()),
        "p_split_max": float(arr.max()),
        "k": int(k),
        "eps": float(eps),
        "r": [float(v) for v in r_arr],
        "split_informativeness": split_informativeness(r_arr, eps),
    }


def tva_defining_features(X: pd.DataFrame, *, r=None, k: int = 2,
                          eps: float = 0.5,
                          rng: Optional[np.random.Generator] = None,
                          n_splits: int = 9, n_init: int = 10,
                          q_fdr: float = 0.05, dispersion: str = "mom",
                          merge: str = "exchangeable_ruger") -> pd.DataFrame:
    """Per-feature multi-split thinning test, BH-corrected across features.

    Per split: cluster ``X1``, then Kruskal-Wallis each ``X2`` column across
    those labels. Per feature: merge the ``n_splits`` p-values with the
    exchangeable rule, then apply Benjamini-Hochberg across features.
    """
    if not isinstance(X, pd.DataFrame):
        raise TypeError("tva_defining_features expects a DataFrame of counts")
    if rng is None:
        rng = np.random.default_rng()
    Xi = X.values.astype(np.int64)
    p = Xi.shape[1]
    r_arr = _resolve_r(Xi, r, dispersion)

    if Xi.shape[0] < 2 * k or p == 0:
        return pd.DataFrame({
            "feature": list(X.columns),
            "p_value": np.ones(p), "p_value_bh": np.ones(p),
            "is_defining": np.zeros(p, dtype=bool),
            "top_cluster": np.full(p, -1, dtype=int),
            "n_splits": np.full(p, int(n_splits), dtype=int),
            "k": np.full(p, int(k), dtype=int),
        })

    P = np.ones((int(n_splits), p), dtype=float)
    T = np.full((int(n_splits), p), -1, dtype=int)
    for m in range(int(n_splits)):
        X1, X2 = nb_thin(Xi, r=r_arr, eps=eps, rng=rng)
        P[m], T[m] = _per_feature_pvalues(X1, X2, k=k, seed=m, n_init=n_init)

    merged = np.asarray([_merge(P[:, j], merge) for j in range(p)], dtype=float)
    adj, reject = benjamini_hochberg(merged, q=q_fdr)
    # modal top cluster across splits (from X1 only)
    top = np.array([np.bincount(T[:, j][T[:, j] >= 0]).argmax()
                    if (T[:, j] >= 0).any() else -1 for j in range(p)], dtype=int)
    return pd.DataFrame({
        "feature": list(X.columns),
        "p_value": merged,
        "p_value_bh": adj,
        "is_defining": reject,
        "p_split_median": np.median(P, axis=0),
        "p_split_min": P.min(axis=0),
        "p_split_max": P.max(axis=0),
        "top_cluster": top,
        "r": r_arr,
        "n_splits": np.full(p, int(n_splits), dtype=int),
        "k": np.full(p, int(k), dtype=int),
    })


def tva_report(counts: pd.DataFrame, *, k: int, rng: np.random.Generator,
               eps: float = 0.5, n_splits: int = 9, n_init: int = 10,
               q_fdr: float = 0.05, dispersion: str = "mom",
               merge: str = "exchangeable_ruger",
               min_distinct: int = 3, min_columns: int = 3,
               r=None, r_source: str | None = None) -> Dict[str, Any]:
    """Screen, thin, test and assemble the full TVA block.

    This is the entry point ``core.run`` uses. It never silently proceeds on
    non-thinnable data: if screening leaves fewer than ``min_columns`` usable
    count columns the block reports ``status="skipped"`` with the reasons.

    Negative-binomial thinning is nominally valid only when its size parameter
    is known independently of the matrix being split. Therefore the default,
    which estimates ``r`` from ``counts``, is explicitly exploratory. Supplying
    ``r`` is not enough by itself: callers must also declare its provenance as
    ``r_source="simulation_oracle"`` or ``r_source="external_validated"``.
    Exploratory p-values and feature calls are retained only under keys prefixed
    with ``diagnostic_`` and are never converted into an archetype verdict.
    """
    valid_r_sources = {"simulation_oracle", "external_validated"}
    if r_source is not None and r_source not in valid_r_sources:
        raise ValueError(
            "r_source must be one of 'simulation_oracle', "
            "'external_validated', or None"
        )
    if r is None and r_source is not None:
        raise ValueError("r_source requires an explicitly supplied r")

    nominal_inference = r is not None and r_source in valid_r_sources
    if r is None:
        parameter_source = "estimated_from_same_matrix"
    elif r_source is None:
        parameter_source = "supplied_without_validated_source"
    else:
        parameter_source = r_source

    screening = screen_thinnable(counts, min_distinct=min_distinct)
    base: Dict[str, Any] = {
        "method": "multi-split NB data thinning (Neufeld et al. 2024 JMLR 25(57); "
                  "merged by Gasparin, Wang & Ramdas 2025 PNAS 122:e2410849122)",
        "screening": screening,
        "counts_columns_offered": [str(c) for c in counts.columns],
        "dispersion_parameter_source": parameter_source,
        "nominal_inference_valid": nominal_inference,
    }
    if screening["n_kept"] < min_columns:
        base.update({
            "status": "skipped",
            "reason": f"only {screening['n_kept']} of {counts.shape[1]} count "
                      f"columns are NB-thinnable (need >= {min_columns}); see "
                      f"'screening.dropped' for the reason on each",
            "archetypes_real": None,
        })
        return base
    if len(counts) < 4 * max(k, 2):
        base.update({
            "status": "skipped",
            "reason": f"n = {len(counts)} is too small for a {k}-group thinned "
                      f"split (need >= {4 * max(k, 2)}); use "
                      f"amr_clonalshare.archephy for small cohorts",
            "archetypes_real": None,
        })
        return base

    kept = counts.loc[:, screening["keep"]]
    Xi = kept.values.astype(np.int64)
    r_for_kept = r
    if r is not None:
        supplied = np.atleast_1d(np.asarray(r, dtype=float))
        if supplied.size == counts.shape[1]:
            r_for_kept = supplied[np.asarray(screening["keep_idx"], dtype=int)]
    r_arr = _resolve_r(Xi, r_for_kept, dispersion)

    sep = tva_test_separation(Xi, r=r_arr, k=k, eps=eps, rng=rng,
                              n_splits=n_splits, n_init=n_init, merge=merge)
    per_feat = tva_defining_features(kept, r=r_arr, k=k, eps=eps, rng=rng,
                                     n_splits=n_splits, n_init=n_init,
                                     q_fdr=q_fdr, merge=merge)
    common = {
        "counts_columns_used": list(screening["keep"]),
        "dispersion_estimator": dispersion,
        "r_per_feature": {c: float(v) for c, v in zip(screening["keep"], r_arr)},
        "eps": float(eps),
        "k_tested": int(k),
        "merge_rule": merge,
        "n_splits": int(n_splits),
        "split_informativeness": sep["split_informativeness"],
    }
    base.update(common)
    if nominal_inference:
        base.update({
            "status": "ok",
            "separation": sep,
            "p_value_separation": sep["p_value"],
            "archetypes_real": bool(sep["p_value"] < 0.05),
            "per_feature": per_feat.to_dict("records"),
            "n_defining_thinned": int(per_feat["is_defining"].sum()),
        })
    else:
        base.update({
            "status": "exploratory_only",
            "reason": (
                "negative-binomial size parameters were estimated from the "
                "matrix being split or were supplied without a validated source; "
                "the resulting p-values are diagnostic, not nominal"
            ),
            "separation": None,
            "p_value_separation": None,
            "archetypes_real": None,
            "per_feature": None,
            "n_defining_thinned": None,
            "diagnostic_separation": sep,
            "diagnostic_p_value_separation": sep["p_value"],
            "diagnostic_per_feature": per_feat.to_dict("records"),
            "diagnostic_n_defining_thinned": int(per_feat["is_defining"].sum()),
        })
    return base


# --------------------------------------------------------------------------- #
# 7. Binary -> count helpers
# --------------------------------------------------------------------------- #
def auto_prefix_groups(columns: Sequence[str], *, min_group: int = 2,
                       prefix_len: int = 3) -> Dict[str, List[str]]:
    """Group columns by their leading alphabetic prefix.

    Groups columns by the first ``prefix_len`` alphabetic characters and keeps
    groups with at least ``min_group`` members, e.g. ``abc_1, abc_2, xyz_1,
    xyz_2`` gives ``{"abc": [...], "xyz": [...]}``.

    This is a convenience for panels whose feature names encode locus families.
    It is a heuristic, not a biological grouping: verify the groups it returns
    before relying on them, and pass explicit ``groups=`` to
    :func:`binary_layer_to_counts` when the panel's structure is known. In
    particular, a one-hot layer (exactly one member positive per row, e.g.
    capsule typing) aggregates to a constant and is not a count.
    """
    import re

    out: Dict[str, List[str]] = {}
    pat = re.compile(rf"^([A-Za-z]{{1,{prefix_len}}})")
    for col in columns:
        m = pat.match(str(col))
        if not m:
            continue
        out.setdefault(m.group(1).lower(), []).append(col)
    return {k: v for k, v in out.items() if len(v) >= min_group}


def binary_layer_to_counts(
    df: pd.DataFrame, *, layer_name: str,
    groups: Optional[Mapping[str, Sequence[str]]] = None,
) -> pd.DataFrame:
    """Aggregate a binary presence/absence layer into count columns.

    Binary features are not thinnable: a Bernoulli variable is its own
    sufficient statistic, so no non-trivial independent split exists. Counts
    are, which is why layers must be coarse-grained first. ``groups`` maps a
    group name to member columns (each becomes a row-sum); with ``groups=None``
    a single ``n_<layer>`` row-sum column is produced.

    Coarse-graining does not guarantee the result is a *count* in the modelling
    sense - see :func:`screen_thinnable`, which rejects the aggregates that are
    still indicators in disguise.
    """
    if df.empty:
        return pd.DataFrame(index=df.index)
    X = df.astype(int)
    if groups is None:
        return pd.DataFrame({f"n_{layer_name}": X.sum(axis=1).astype(int)},
                            index=df.index)
    cols = {}
    for grp, members in groups.items():
        members = [m for m in members if m in X.columns]
        if members:
            cols[grp] = X.loc[:, members].sum(axis=1).astype(int)
    return pd.DataFrame(cols, index=df.index)
