"""The lineage share of a call or of the recorded ordering of MIC readings.

The share is the out-of-sample skill of lineage means against the overall
mean in predicting the value of held-out isolates, set against the same score
of permuted lineage labels. Isolates are assigned to folds within the lineages, so the share
refers to further isolates of the represented lineages, not to unseen lineages.

Its target is the lineage share of the represented lineages: the lineages and their sizes
held fixed and the isolates of every lineage drawn again from its distribution, the
between-lineage share of the variance. Its interval draws the isolates of
every lineage again from the smoothed empirical distribution of that lineage,
recomputes the share, and reads the studentized spread of the recomputed
shares around the share of that smoothed scenario; its upper end also covers a
share accounted for by one lineage whose readings sit close to the rest. Within
declared strata the trait is scored stratum by stratum, a call centred on
the prevalence of its stratum, so the share is that of the variance within
strata. The p-value compares the between-lineage sum of squares with its
values under labels permuted within strata (and sampling units), a test
exact under exchangeability within them.

Lineage association does not identify transmission, horizontal transfer,
heritability, clinical resistance or a causal genetic mechanism.
"""
from __future__ import annotations

import hashlib
import math
import typing
from numbers import Integral
from dataclasses import dataclass

import numpy as np

#: Lineages with at least two isolates the share needs. A singleton lineage
#: cannot be predicted out of sample -- its one isolate is in the training
#: fold or the held-out fold, never both -- so the share is scored on the
#: isolates of repeated lineages and singletons are set aside and counted.
MIN_REPEATED_LINEAGES = 2

#: Code given to an isolate with no lineage assignment while labels are
#: coded; the estimators set those isolates aside before fitting and count
#: them in ``n_dropped_untyped``, because an untyped isolate is not a member
#: of any lineage.
_MISSING = "__missing__"

#: Fewer bootstrap draws with a share than this fraction of those asked for
#: leave the interval undefined rather than read from a few draws.
MIN_USABLE_DRAWS = 0.9

__all__ = [
    "MIN_REPEATED_LINEAGES",
    "ShareResult",
    "clonal_share",
    "layer_clonal_share",
]


# --------------------------------------------------------------------- types
@dataclass(frozen=True)
class ShareResult:
    """The lineage share of one trait.

    Attributes
    ----------
    kappa : the out-of-sample skill of the lineage-mean predictor against the
        overall-mean predictor, with the cost of estimating each lineage mean
        from its training members removed lineage by lineage (see
        :func:`_skill`): its expectation is the lineage share of the represented lineages
        whatever the variance within each lineage.
    kappa_adj : the estimate to report: kappa set against the same score of
        the permuted labellings, ``(kappa - null_mean) / (1 - null_mean)``,
        which keeps the estimate at zero where the labels account for nothing.
    observed_low, observed_high : 95% interval for the lineage share of the
        represented lineages, within 0 and 1: the studentized bootstrap of
        :func:`_smoothed_bootstrap`, its quantiles taken no closer to the
        estimate than those of the studentized shares of the permuted
        labellings, and its upper end raised where one lineage alone, within
        the empirical-likelihood region of its readings, can account for more
        of the share (:func:`_single_lineage_gain`).
    observed_se : standard error of the estimate that the interval is
        studentized by: the jackknife variance of :func:`_moment_se`, the
        spread left by the fold assignments and the variance of the share where
        the lineages do not differ, under the permuted labellings.
    null_mean : mean kappa when the lineage labels are permuted, the control
        against which the share is read; near zero.
    p_value : permutation p-value, Phipson-Smyth corrected, of the
        between-lineage sum of squares of the trait (see
        :func:`_between_lineage_ss`) against its values under the permuted
        labels. For a call and no strata this is the permutation form of the
        chi-square test of the lineage-by-outcome table, and for scores the
        permutation form of the Kruskal-Wallis test on them. It cannot fall
        below ``p_floor = 1 / (n_perm + 1)``.
    n, n_groups, prevalence : collection description carried with the estimate.
    support : share of the typed isolates sitting in a lineage with at least
        two members; the share is scored on those isolates only. ``n``
        counts the retained isolates, ``n_scored`` those scored and
        ``n_singletons_set_aside`` the ones set aside.
    estimable : ``kappa_adj`` finite and at least
        :data:`MIN_REPEATED_LINEAGES` repeated lineages.
    cv_sd : spread of kappa across repeated fold assignments; a diagnostic of
        the estimator, not of the collection, and never a confidence interval.
    n_dropped_non_finite, n_dropped_untyped : isolates set aside because the
        trait was missing or the lineage label was.
    n_boot_used : bootstrap draws that gave a share.
    """

    kappa: float
    kappa_adj: float
    null_mean: float
    p_value: float
    n: int
    n_groups: int
    prevalence: float
    support: float
    missing_share: float
    estimable: bool
    cv_sd: float
    observed_low: float = float("nan")
    observed_high: float = float("nan")
    observed_se: float = float("nan")
    n_dropped_non_finite: int = 0
    n_dropped_untyped: int = 0
    #: isolates of singleton lineages, typed and read but not scored
    n_singletons_set_aside: int = 0
    #: isolates scored: those of lineages with at least two members
    n_scored: int = 0
    #: the quantity the interval is for
    interval_target: str = "lineage share of the represented lineages"
    #: The smallest p-value ``n_perm`` permutations can return.
    p_floor: float = float("nan")
    #: Lineages with at least two analysed isolates.
    n_groups_repeated: int = 0
    #: Effective number of lineages among the analysed isolates, the inverse
    #: of the sum of squared lineage shares (Hill number of order two).
    effective_groups: float = float("nan")
    #: Share of the analysed isolates in the largest lineage.
    largest_group_share: float = float("nan")
    #: Positive outcomes among isolates in repeated lineages, for a call.
    n_positive_repeated: float = float("nan")
    n_boot_used: int = 0
    #: The denominator of the skill averaged over the fold assignments,
    #: over the total sum of squares (:func:`_fold_average`).
    fold_ratio: float = float("nan")
    #: Why a number is missing.
    reason: str = ""
    #: Strata the permutations of the test respected (strata crossed with
    #: sampling units): 0 when there were none.
    n_strata: int = 0

    def as_dict(self) -> dict:
        return {f.name: getattr(self, f.name) for f in self.__dataclass_fields__.values()}


# ----------------------------------------------------------------- internals
def _rng(seed) -> np.random.Generator:
    if isinstance(seed, np.random.Generator):
        return seed
    return np.random.default_rng(seed)


def _is_untyped(v) -> bool:
    """Recognise missing scalar labels without evaluating pandas.NA as bool."""
    if v is None:
        return True
    return str(v).strip().lower() in {
        "", "nan", "na", "none", "<na>", "nat", "null", "n/a", "__missing__"}


def _codes(labels) -> np.ndarray:
    """Integer codes for any hashable label vector, missing values included.

    A missing lineage is coded as one ``__missing__`` level so that a label
    vector of any type can be coded; the estimators set those rows aside
    before fitting (see ``n_dropped_untyped``), because an untyped isolate is
    not a member of any lineage.
    """
    arr = np.asarray(labels, dtype=object)
    clean = np.empty(arr.shape, dtype=object)
    for i, v in enumerate(arr.ravel()):
        clean.ravel()[i] = _MISSING if _is_untyped(v) else str(v)
    _, codes = np.unique(clean.astype(str), return_inverse=True)
    return codes.astype(np.int64)


def _missing_share(labels) -> float:
    """Share of isolates whose lineage label is absent.

    Reported, not used as a reporting condition. Untyped isolates are set aside before the
    estimate, and if typing failed informatively -- as it does on the shipped
    *S. suis* collection at sequence-type resolution, where the untyped isolates
    carry about two more non-wild-type results out of thirteen -- the estimate
    describes a subset selected by the typing process. The number is emitted
    so a reader can see how much of the collection the estimate rests on, and so
    that the same collection typed two ways can be compared.
    """
    arr = np.asarray(labels, dtype=object).ravel()
    return float(np.mean([_is_untyped(v) for v in arr])) if arr.size else float("nan")


def _folds(n: int, k: int, rng: np.random.Generator) -> np.ndarray:
    order = rng.permutation(n)
    fold = np.empty(n, dtype=np.int64)
    fold[order] = np.arange(n) % k
    return fold


def _dealt_folds(code: np.ndarray, k: int, rng: np.random.Generator) -> np.ndarray:
    """Folds assigned lineage by lineage, so that no lineage sits in one fold.

    The isolates are taken lineage by lineage, the lineages in random order and
    the isolates of each in random order, and assigned to the folds in turn from a
    random first fold. The folds stay balanced in size, and a lineage of m
    isolates spreads over min(m, k) folds, so every held-out isolate of a
    repeated lineage keeps at least one member of its lineage in training.
    Folds drawn at random instead leave all of a small lineage in one fold now
    and then, and its held-out isolates are then predicted by the grand mean:
    a loss that grows with the share itself, which the permuted-label control,
    where there is no share, cannot measure. On lineages of three isolates it
    took the estimate about 0.016 below a share of 0.7 on the estimator grid,
    and assigning the folds this way removes it. The folds depend on the labels
    only, never on the trait, and the permuted runs assign them by the permuted
    labels.
    """
    code = np.asarray(code)
    n = code.size
    # A random priority for every isolate; a lineage takes the smallest
    # priority of its isolates, so the order of the lineages is random and
    # depends on which isolates belong together, not on what they are called.
    priority = rng.permutation(n)
    first = np.full(int(code.max()) + 1 if n else 0, n, dtype=np.int64)
    np.minimum.at(first, code, priority)
    order = np.lexsort((priority, first[code]))
    fold = np.empty(n, dtype=np.int64)
    fold[order] = (np.arange(n) + int(rng.integers(k))) % k
    return fold


def _lineage_variances(X: np.ndarray, code: np.ndarray, G: int) -> np.ndarray:
    """Unbiased within-lineage variance of every column, zero for a lineage of
    one isolate: (G, p)."""
    m = np.bincount(code, minlength=G).astype(float)
    out = np.zeros((G, X.shape[1]))
    for j in range(X.shape[1]):
        # two passes, the squares taken about the lineage means, so that an
        # origin far from the values costs no digits
        means = np.bincount(code, weights=X[:, j], minlength=G) / np.maximum(m, 1)
        square = np.bincount(code, weights=(X[:, j] - means[code]) ** 2, minlength=G)
        out[:, j] = np.where(m > 1, square / np.maximum(m - 1, 1), 0.)
    return np.maximum(out, 0.)


def _skill_parts(X: np.ndarray, code: np.ndarray, fold: np.ndarray):
    """The two sums of the score of one fold draw (see :func:`_skill`): the
    held-out squared error of the lineage means less its penalty, and the
    held-out squared error of the training grand mean less its penalty."""
    X = np.atleast_2d(np.asarray(X, dtype=float))
    if X.shape[0] != fold.shape[0]:
        X = X.T
    p = X.shape[1]
    G = int(code.max()) + 1
    s2 = _lineage_variances(X, code, G)
    sse = sst = penalty = penalty0 = 0.
    for f in np.unique(fold):
        tr, te = fold != f, fold == f
        if not tr.any() or not te.any():
            continue
        ctr = code[tr]
        counts = np.bincount(ctr, minlength=G).astype(float)
        sums = np.empty((G, p))
        for j in range(p):
            sums[:, j] = np.bincount(ctr, weights=X[tr, j], minlength=G)
        grand = X[tr].mean(axis=0)
        means = np.where(counts[:, None] > 0, sums / np.maximum(counts[:, None], 1),
                         grand[None, :])
        g = code[te]
        sse += float(((X[te] - means[g]) ** 2).sum())
        sst += float(((X[te] - grand[None, :]) ** 2).sum())
        with np.errstate(invalid="ignore", divide="ignore"):
            penalty += float(np.where(counts[g, None] > 0, s2[g] / np.maximum(counts[g, None], 1), 0.).sum())
        # a sum over lineages taken exactly, so that renaming the lineages
        # cannot move it in its last digit
        penalty0 += float(te.sum() * math.fsum((counts[:, None] * s2).ravel()) / tr.sum() ** 2)
    return sse - penalty, sst - penalty0




def _skill(X: np.ndarray, code: np.ndarray, fold: np.ndarray) -> float:
    """Out-of-sample variance share of ``code`` for columns of ``X``.

    Held-out rows are predicted by the mean of their group among the training
    rows; a group absent from training falls back to the training grand mean,
    which is what an honest predictor would do and is the mechanism by which
    extra levels stop paying for themselves. The score is

        1 - (SSE(group means) - penalty) / (SSE(marginal mean) - penalty0)

    on the pooled sums of squares over all columns. A held-out isolate of a
    lineage with t training members is predicted by a mean of t readings, so
    its expected squared error is the lineage's variance times 1 + 1/t: the
    lineage's own variance, which the share leaves to the lineage, and the
    variance of the training mean, which is the cost of estimating it. The
    penalty is that cost, the unbiased within-lineage variance over t, summed
    over the held-out isolates; penalty0 is the variance of the training grand
    mean the baseline pays. Both are exact in expectation for any lineage
    whose isolates are exchangeable, whatever the variance of each lineage,
    so a lineage whose calls or readings vary less than the others does not
    move the estimate. Uncorrected, 1 - SSE/SSE(marginal) is the Brier skill
    score for one binary column and the multivariate out-of-sample R-squared
    for several.
    """
    numerator, denominator = _skill_parts(X, code, fold)
    return float("nan") if denominator <= 0 else 1.0 - numerator / denominator


def _repeated_skill(X, code, *, folds, repeats, rng, keep_folds=None):
    draws = [_dealt_folds(code, folds, rng) for _ in range(repeats)]
    if keep_folds is not None:
        keep_folds.extend(draws)
    vals = np.array([_skill(X, code, fold) for fold in draws])
    return float(np.nanmean(vals)), float(np.nanstd(vals)), vals


def _jumped(rng: np.random.Generator):
    """A bit generator far ahead of ``rng`` on the same stream, taken without
    advancing ``rng``. Bit generators without a jump are given a fresh one
    seeded from their state, which is equally deterministic."""
    jump = getattr(rng.bit_generator, "jumped", None)
    if jump is not None:
        return jump()
    else:
        digest = hashlib.sha256(repr(rng.bit_generator.state).encode()).digest()
        return np.random.PCG64(int.from_bytes(digest[:8], "little"))


def _between_lineage_ss(X: np.ndarray, code: np.ndarray) -> float:
    """Between-lineage sum of squares of the columns of ``X``, less a constant.

    ``sum_g n_g (xbar_g - xbar)^2`` summed over columns equals
    ``sum_g S_g^2 / n_g - S^2 / n`` with ``S_g`` the lineage sums; the last
    term does not change when the labels are permuted, so the statistic of
    the permutation test is the first. ``code`` numbers the lineages from
    zero without gaps, as :func:`_codes` and its permutations do. Each
    lineage's values are summed in sorted order, so that labellings that put
    the same values in the same lineages give the same number to the last
    bit.
    """
    counts = np.bincount(code)
    total = 0.0
    for j in range(X.shape[1]):
        order = np.lexsort((X[:, j], code))
        sums = np.bincount(code[order], weights=X[order, j])
        total += float(np.sum(sums ** 2 / counts))
    return total


def _phipson_smyth(exceed: int, n_perm: int) -> float:
    """Permutation p-value that cannot be zero (Phipson & Smyth 2010)."""
    return (exceed + 1.0) / (n_perm + 1.0)


def _debias(kappa: float, null_mean: float) -> float:
    """Set the score against the score of the permuted labellings.

    ``(kappa - null_mean) / (1 - null_mean)``: zero where the labels score as
    the permuted ones do, and one at a perfect prediction. A finite sampled
    permutation control can reach one on a small design. At that boundary
    the map is undefined; values above one also violate the positive
    rescaling denominator and must not produce a reportable estimate.
    """
    if not np.isfinite(kappa) or not np.isfinite(null_mean):
        return float("nan")
    if null_mean >= 1.0:
        return float("nan")
    return float((kappa - null_mean) / (1.0 - null_mean))


def _permute_within(code: np.ndarray, strata: np.ndarray,
                    rng: np.random.Generator) -> np.ndarray:
    """The lineage labels exchanged among the isolates of each stratum only."""
    out = code.copy()
    for s in range(int(strata.max()) + 1):
        idx = np.flatnonzero(strata == s)
        out[idx] = code[idx[rng.permutation(idx.size)]]
    return out


def _share_of(values, sizes, law, cell_of, cell_lineage):
    """The between-lineage share of the variance when every isolate of cell c
    takes value v_k with probability law[c, k], the cells' sizes fixed."""
    mean = law @ values
    second = law @ values ** 2
    counts = np.bincount(cell_of, minlength=law.shape[0]).astype(float)
    n = counts.sum()
    mu = np.bincount(cell_lineage, weights=counts * mean, minlength=sizes.size) / sizes
    m2 = np.bincount(cell_lineage, weights=counts * second, minlength=sizes.size) / sizes
    w = sizes / n
    centre = math.fsum(w * mu)
    between = math.fsum(w * (mu - centre) ** 2)
    within = math.fsum(w * (m2 - mu ** 2))
    total = between + within
    return between / total if total > 0 else float("nan")


def _fold_average(y: np.ndarray, code: np.ndarray, ratio: float) -> float:
    """The skill of :func:`_skill` averaged over the fold assignments.

    Averaged over fold assignments, the corrected numerator of a lineage of m
    isolates is exactly m times its unbiased variance (the assignment is
    exchangeable within the lineage and the correction fixes the multiple),
    and the denominator is ``ratio`` times the total sum of squares, the
    ratio read from the fold assignments of the data. So the average is
    1 - sum_g m_g s_g^2 / (ratio * SST), with no fold noise."""
    sizes = np.bincount(code).astype(float)
    means = np.bincount(code, weights=y) / sizes
    # the squares about the lineage means (two passes), and the sums over
    # lineages taken exactly, so that neither the origin of the values nor
    # the names of the lineages can move the result in its last digits
    within = np.bincount(code, weights=(y - means[code]) ** 2)
    numerator = math.fsum(within * sizes / np.maximum(sizes - 1., 1.))
    total = float(np.sum((y - y.mean()) ** 2))
    return 1. - numerator / (ratio * total) if total > 0 and ratio > 0 else float("nan")


def _smoothed_laws(counts, allowed):
    """The distribution of every cell smoothed as the bootstrap draws from it: the
    Dirichlet posterior mean under Perks' prior of total mass one over the
    categories the cell may take, (count + 1/K) / (m + 1)."""
    size = allowed.sum(1, keepdims=True)
    return (counts + allowed / size) / (counts.sum(1, keepdims=True) + 1.)


def _share_derivatives(values, counts, ratio, c_null):
    """The derivative of the fold-averaged share (:func:`_fold_average`, set
    against the control ``c_null``) with respect to the weight of one isolate
    of lineage g in category k, the infinitesimal jackknife; it depends on g
    and k alone:

        D_gk = -(dN_gk T - N dT_k) / (ratio T^2) / (1 - c_null),
        dN_gk = f'(m_g) SS_g + f(m_g) (v_k - ybar_g)^2,  dT_k = (v_k - ybar)^2,

    with N = sum_g f(m_g) SS_g, f(m) = m / (m - 1), SS_g the within-lineage
    and T the total sum of squares, the values of the categories held fixed.
    ``values`` holds the value of every category and ``counts`` the counts of
    every lineage over the categories (lineage by category). None if the
    share is not defined."""
    m = counts.sum(1)
    n = m.sum()
    sums = counts @ values
    squares = counts @ (values * values)
    ybar_g = sums / m
    within = squares - sums * sums / m
    f = m / np.maximum(m - 1., 1.)
    fprime = -1. / np.maximum(m - 1., 1.) ** 2
    N = math.fsum(f * within)
    ybar = sums.sum() / n
    T = float(squares.sum() - n * ybar * ybar)
    if not (T > 0 and ratio > 0 and c_null < 1):
        return None
    dN = fprime[:, None] * within[:, None] + f[:, None] * (values[None, :] - ybar_g[:, None]) ** 2
    dT = (values - ybar) ** 2
    return -(dN * T - N * dT[None, :]) / (ratio * T * T) / (1. - c_null)


def _moment_se(own, prior, ratio, c_null):
    """Standard error of the fold-averaged share when the isolates of every
    lineage are drawn again from its distribution and the lineages are held
    fixed: sqrt(sum_g m_g^2 / (m_g - 1) Var_g(D_g.)), D the derivatives of
    :func:`_share_derivatives` and their variance over the categories taken
    under the lineage's smoothed distribution (:func:`_smoothed_laws`) rather
    than its empirical one, so that a lineage whose isolates share one value
    is not given a variance of zero. Under the empirical distribution it is
    the infinitesimal jackknife with the small-sample factor m_g / (m_g - 1).

    ``own`` holds, for every lineage, the sums of v^0 to v^4 over its
    isolates and ``prior`` the means of v^0 to v^4 over the categories its
    smoothed distribution spreads its prior mass over. D_g(x) of
    :func:`_share_derivatives` is a quadratic in the value x,
    K0 (c0_g + c1_g x + c2_g x^2) with K0 = 1 / (ratio T^2 (1 - c_null)),
    c1_g = 2 f_g T ybar_g - 2 N ybar and c2_g = N - f_g T, so its variance
    under the smoothed distribution, whose moments are (own + prior) / (m_g + 1), is
    K0^2 (c1^2 Var x + 2 c1 c2 Cov(x, x^2) + c2^2 Var x^2): the cost of a draw
    is linear in the isolates and the categories."""
    m = own[:, 0]
    n = math.fsum(m)
    s1, s2 = own[:, 1], own[:, 2]
    ybar_g = s1 / m
    within = s2 - s1 * s1 / m
    f = m / np.maximum(m - 1., 1.)
    N = math.fsum(f * within)
    ybar = math.fsum(s1) / n
    T = float(math.fsum(s2) - n * ybar * ybar)
    if not (T > 0 and ratio > 0 and c_null < 1):
        return float("nan")
    K0 = 1. / (ratio * T * T * (1. - c_null))
    c1 = 2. * f * T * ybar_g - 2. * N * ybar
    c2 = N - f * T
    mu = (own[:, 1:] + prior[:, 1:]) / (m + 1.)[:, None]
    var_x = mu[:, 1] - mu[:, 0] ** 2
    cov = mu[:, 2] - mu[:, 0] * mu[:, 1]
    var_x2 = mu[:, 3] - mu[:, 1] ** 2
    variance = np.maximum(K0 * K0 * (c1 * c1 * var_x + 2. * c1 * c2 * cov + c2 * c2 * var_x2), 0.)
    return float(np.sqrt(math.fsum(m * m / np.maximum(m - 1., 1.) * variance)))


def _standardised(values, weights):
    """The values centred and scaled by the mean and standard deviation of
    the weights over them; the share and its derivatives do not see it, and
    the moments of :func:`_moment_se` keep their digits."""
    total = weights.sum()
    centre = float(weights @ values) / total
    scale = float(np.sqrt(weights @ (values - centre) ** 2 / total))
    return (values - centre) / scale if scale > 0 else values - centre


class _Sampler:
    """Draws the isolates of every lineage again from its smoothed distribution.

    The smoothed distribution of lineage g, (counts + 1/K_g) / (m_g + 1) over its K_g
    allowed categories (:func:`_smoothed_laws`), is the mixture of its own
    isolates, each with weight 1 / (m_g + 1), and of the uniform distribution on the
    allowed categories, with the remaining weight 1 / (m_g + 1). An isolate
    is drawn from it in two steps from two uniforms: the first picks the part
    of the mixture, the second the member (by its position among the
    isolates of the lineage) or the category (in an order fixed by the count
    of the category in the collection and then by its code, so that
    renaming the lineages or recoding the values does not change a draw).
    A draw costs time linear in the isolates."""

    def __init__(self, categories, cat_stratum, code, seen, world_count):
        self.categories = categories
        self.code = code
        m = np.bincount(code).astype(float)
        self.m = m
        order = np.argsort(code, kind="stable")
        self.members = order
        self.start = np.r_[0, np.cumsum(m)[:-1]].astype(np.int64)
        patterns, self.set_of = np.unique(seen, axis=0, return_inverse=True)
        self.set_of = self.set_of.ravel()
        rank = np.lexsort((np.arange(world_count.size), -world_count))
        self.allowed = [rank[pattern[cat_stratum[rank]]] for pattern in patterns]

    def draw(self, rng):
        n = self.code.size
        u, v = rng.random(n), rng.random(n)
        m = self.m[self.code]
        pick = self.start[self.code] + np.minimum((v * m).astype(np.int64), m.astype(np.int64) - 1)
        out = self.categories[self.members[pick]]
        prior = u * (m + 1.) < 1.
        sets = self.set_of[self.code]
        for s, allowed in enumerate(self.allowed):
            idx = np.flatnonzero(prior & (sets == s))
            if idx.size:
                out[idx] = allowed[np.minimum((v[idx] * allowed.size).astype(np.int64), allowed.size - 1)]
        return out

    def se(self, values, draw, ratio, c_null, code=None):
        """:func:`_moment_se` of one sample of categories; ``code``
        replaces the lineages by a labelling permuted within strata, which
        keeps every lineage's strata and so its smoothed distribution's support."""
        weights = np.bincount(draw, minlength=values.size).astype(float)
        v = _standardised(values, weights)
        y = v[draw]
        code = self.code if code is None else code
        own = np.column_stack([np.bincount(code, weights=y ** j, minlength=self.m.size) for j in range(5)])
        prior_sets = np.array([[np.mean(v[a] ** j) for j in range(5)] for a in self.allowed])
        return _moment_se(own, prior_sets[self.set_of], ratio, c_null)


#: Directions of the scan of :func:`_single_lineage_gain` in the plane of the
#: first two moments of a lineage's distribution.
SINGLE_LINEAGE_ANGLES = 64


def _single_lineage_gain(values, code, categories, ratio, c_null):
    """How far the share can rise when one lineage alone accounts for more of it.

    The studentized bootstrap reads the spread of the share at the readings,
    and the derivative of the share in a lineage grows with the lineage's
    distance from the rest: where the share rests on a few lineages of few
    isolates and one of them happens to read close to the rest, the standard
    error at the readings is too small and the upper end falls short of the
    share. This returns the largest rise that one lineage can give the
    share while its distribution stays compatible with its own readings.

    The distribution P of lineage g ranges over its empirical-likelihood region (Owen
    2001): distributions on the categories it was read in with 2 m_g KL(Phat_g || P)
    at most the 95 % point of F(1, m_g - 1), the small-sample calibration of
    the chi-square(1) limit, so that each end of the region answers at
    97.5 %. Every other lineage keeps its readings, and the values of the
    categories are those of the data. The share of such a scenario is written
    on the scale of the fold-averaged share (:func:`_fold_average`): the
    lineages kept contribute f(m_h) SS_h, their unbiased sum of squares, and
    lineage g m_g Var_P, the variance of its distribution, since P stands for the distribution
    itself; the total sum of squares is that of the mixture. The share
    depends on P through its first two moments (mu, M2) only, and its
    sublevel sets are convex there (the epigraph M2 >= Q(mu) of a convex
    quadratic, or the hypograph of a concave one) at every level above
    1 - n / (m_g ratio), which is below -1 whenever the lineage holds at
    most half of the isolates; the region's image in that plane is convex,
    so the largest share is reached at an extreme point of the image, which
    its support function exposes. The support function is taken in
    :data:`SINGLE_LINEAGE_ANGLES` directions, the best refined by golden
    section, each a one-dimensional convex problem solved by bisection: the
    distribution maximising <c, P> is P_k proportional to Phat_k / (t - c_k), t above
    the largest c_k read, and the end of the bisection kept lies inside the
    region, so every share found is reached by a distribution of the region.
    Returns (share of the best such scenario - share of the readings) / (1 -
    c_null), on the scale of the debiased estimate, the largest over the
    lineages and directions; nan where the share is not defined."""
    nan = float("nan")
    v = np.asarray(values, dtype=float)
    code = np.asarray(code, dtype=np.int64)
    categories = np.asarray(categories, dtype=np.int64)
    n = code.size
    y = v[categories]
    centre = math.fsum(y) / n
    total = float(np.sum((y - centre) ** 2))
    if not (total > 0 and ratio > 0 and c_null < 1):
        return nan
    z = (v - centre) / np.sqrt(total / n)
    # the cells read: (lineage, category) with their counts, by lineage
    K = v.size
    cells, counts = np.unique(code * K + categories, return_counts=True)
    lin, cat = cells // K, cells % K
    m = np.bincount(code).astype(float)
    G = m.size
    start = np.r_[0, np.flatnonzero(np.diff(lin)) + 1]
    ph = counts / m[lin]
    zc = z[cat]
    s1 = np.bincount(lin, weights=counts * zc, minlength=G)
    s2 = np.bincount(lin, weights=counts * zc * zc, minlength=G)
    ss = np.maximum(s2 - s1 * s1 / m, 0.)
    f = m / np.maximum(m - 1., 1.)
    N = math.fsum(f * ss)
    S1, S2 = math.fsum(s1), math.fsum(s2)
    base = 1. - N / (ratio * (S2 - S1 * S1 / n))
    from scipy import stats
    radius = stats.f.ppf(.95, 1, np.maximum(m - 1., 1.)) / (2. * m)
    logph = np.log(ph)

    def boundary_share(theta):
        """The share when lineage g takes the distribution of its region that
        maximises <cos theta z + sin theta z^2, P>, theta (rows, G)."""
        c = np.cos(theta)[:, lin] * zc + np.sin(theta)[:, lin] * (zc * zc)
        top = np.maximum.reduceat(c, start, axis=1)
        lo = np.full(top.shape, -30.)
        hi = np.full(top.shape, 30.)

        def law(s):
            q = ph / ((top + np.exp(s))[:, lin] - c)
            return q, np.add.reduceat(q, start, axis=1)

        for _ in range(40):
            mid = .5 * (lo + hi)
            q, norm = law(mid)
            # KL(Phat || P) at P proportional to Phat / (t - c), decreasing
            # in t = top + exp(s); the end kept is always inside the region
            kl = np.add.reduceat(ph * (logph - np.log(q)), start, axis=1) + np.log(norm)
            far = kl > radius
            lo = np.where(far, mid, lo)
            hi = np.where(far, hi, mid)
        q, norm = law(hi)
        P = q / norm[:, lin]
        mu = np.add.reduceat(P * zc, start, axis=1)
        m2 = np.add.reduceat(P * zc * zc, start, axis=1)
        Nc = N - f * ss + m * (m2 - mu * mu)
        S1c = S1 - s1 + m * mu
        S2c = S2 - s2 + m * m2
        share = 1. - Nc / (ratio * (S2c - S1c * S1c / n))
        return np.where(np.isfinite(share), share, -np.inf)

    step = 2. * np.pi / SINGLE_LINEAGE_ANGLES
    scan = boundary_share(np.repeat(np.arange(SINGLE_LINEAGE_ANGLES)[:, None] * step, G, axis=1))
    best = scan.max(0)
    # the best direction of every lineage refined by golden section between
    # its neighbours on the scan
    a = scan.argmax(0) * step - step
    b = a + 2. * step
    ratio_g = (np.sqrt(5.) - 1.) / 2.
    for _ in range(20):
        x1, x2 = b - ratio_g * (b - a), a + ratio_g * (b - a)
        both = boundary_share(np.vstack([x1, x2]))
        best = np.maximum(best, both.max(0))
        left = both[0] >= both[1]
        b = np.where(left, x2, b)
        a = np.where(left, a, x1)
    top_share = float(best.max())
    if not np.isfinite(top_share):
        return nan
    return (top_share - base) / (1. - c_null)


def _smoothed_bootstrap(categories, cat_stratum, scorer, code, stratum, *,
                        n_boot, c_null, ratio, fold_sd, rng, null_var=0., perms=(),
                        permuted=()):
    """Studentized bootstrap for the lineage share of the represented lineages.

    The isolates of every lineage are drawn again from the smoothed distribution of
    its readings over the K categories of the strata its isolates were read
    in (:func:`_smoothed_laws`, drawn by :class:`_Sampler`; a category is a
    reading within a stratum, so the stratum of every isolate is drawn with
    its reading, in the proportions of the lineage): a lineage whose readings
    all fall in one category is drawn with a small chance of every other
    category of its strata. The draws are scored as the data were (for MIC
    readings, within stratum at the distribution of the draw), the share is
    computed as its average over the fold assignments
    (:func:`_fold_average`), set against the permuted control of the data,
    and given a normal draw of standard deviation ``fold_sd``, the spread
    that averaging a finite number of fold assignments leaves in the estimate. Every
    draw is studentized by its own standard error, t = (share - eta*) / se,
    eta* the share of the smoothed distribution itself, computed exactly, and se the
    root of the jackknife variance of :func:`_moment_se`, the fold
    variance and ``null_var``, the variance the share has where the lineages
    do not differ (read from the permuted labellings): the jackknife is the
    first-order term, which vanishes where the lineages read alike, and the
    interval would collapse to a point there without the second. The
    interval is (share - se q_975(t), share - se q_025(t)) with se the
    standard error of the data, within 0 and 1: the bootstrap-t interval,
    whose studentization carries the change of the spread of the share with
    the share itself, which the plain bootstrap misreads where a few lineages
    of rare readings account for the share. The distribution of the readings in every cell
    is the only thing resampled, so the interval answers for the lineages and
    their sizes held fixed, and the scores are recomputed in every draw, so
    the uncertainty of the scores, which depend on the distribution of the
    readings, is part of it.

    The labellings ``perms``, permuted within strata, are studentized in the
    same way, t = share / se with ``permuted`` the share of each averaged
    over the fold assignments (:func:`_fold_average`), given the same fold spread,
    and its standard error taken as the data's: under exchangeable labels this is the distribution
    of the data's own studentized share, and a studentized permutation test
    keeps its level asymptotically where the lineages share only their mean
    (Janssen 1997; Chung and Romano 2013). Returns the list of t over the
    draws that gave a share, the standard error of the data and the list of
    t over the permuted labellings."""
    nan = float("nan")
    C = int(code.max()) + 1
    K = int(categories.max()) + 1
    counts = np.zeros((C, K))
    np.add.at(counts, (code, categories), 1.)
    seen = np.zeros((C, int(stratum.max()) + 1), dtype=bool)
    seen[code, stratum] = True
    allowed = seen[:, cat_stratum].astype(float)
    law = _smoothed_laws(counts, allowed)
    sizes = np.bincount(code).astype(float)
    world = law.T @ sizes
    eta_star = _share_of(scorer(world), sizes, law, code, np.arange(C))
    world_count = counts.sum(0)
    sampler = _Sampler(categories, cat_stratum, code, seen, world_count)
    added = fold_sd * fold_sd + null_var
    se_data = float(np.sqrt(sampler.se(scorer(world_count), categories, ratio, c_null) ** 2 + added))
    if not (np.isfinite(eta_star) and np.isfinite(se_data) and se_data > 0):
        return [], nan, []
    values = scorer(world_count)
    t_perm = []
    for perm, share in zip(perms, permuted):
        se = float(np.sqrt(sampler.se(values, categories, ratio, c_null, code=perm) ** 2 + added))
        share = _debias(share, c_null) + float(rng.normal(0., fold_sd))
        if np.isfinite(share) and np.isfinite(se) and se > 0:
            t_perm.append(share / se)
    t = []
    for _ in range(n_boot):
        draw = sampler.draw(rng)
        weights = np.bincount(draw, minlength=K).astype(float)
        drawn_values = scorer(weights)
        y = drawn_values[draw]
        noise = float(rng.normal(0., fold_sd))
        if np.ptp(y) <= 0:
            continue
        share = _debias(_fold_average(y, code, ratio), c_null) + noise
        se = float(np.sqrt(sampler.se(drawn_values, draw, ratio, c_null) ** 2 + added))
        if np.isfinite(share) and np.isfinite(se) and se > 0:
            t.append((share - eta_star) / se)
    return t, se_data, t_perm


def layer_clonal_share(X, lineage, *, folds: int = 5, repeats: int = 20,
                       n_boot: int = 999, n_perm: int = 999, seed=0, strata=None,
                       units=None, categories=None, scorer=None) -> ShareResult:
    """Out-of-sample share of a trait's variance accounted for by lineage.

    Parameters
    ----------
    X : (n,) values of one trait: calls (1 = the positive outcome) or scores.
    lineage : (n,) label vector; MLST sequence type, a BAPS cluster, a clonal
        complex, anything categorical. It is used as a label and never as a
        tree.
    folds, repeats : cross-validation design. ``repeats`` averages over fold
        draws so the estimate does not depend on one fold split.
    n_boot : draws of the smoothed bootstrap that gives the interval for the
        represented lineages; 0 gives the estimate and the test without it.
    n_perm : permutations of the lineage label. Each permuted labelling is
        scored once by cross-validation, and the mean of those scores is the
        control the estimate is set against; the same labellings give the
        p-value.
    strata : optional (n,) labels, such as the testing laboratory. The trait
        is then scored within each stratum: without a scorer, every value is
        centred on the mean of its stratum (a call on the prevalence of its
        stratum, the reading with one cut point of :mod:`.mic_order`), and
        the scores are recomputed in every bootstrap draw. Lineage labels are
        exchanged only among isolates of one stratum, which keeps the test
        exact, and the bootstrap draws the stratum of every isolate with its
        reading. The share is then the share of the variance within strata
        that lies between the lineages: a stratum that follows lineage
        exactly removes the lineage difference with it.
    units : optional (n,) labels of the sampling unit, such as the farm or
        host. The permutations of the test then exchange labels only within
        units (and strata): a test of lineage differences within units, exact
        when isolates of one unit are exchangeable, whatever the units share.
        The estimate, its control and its interval describe the lineages and
        the represented units and do not change.
    categories : optional (n,) labels of the reading behind each value, such
        as the well of an MIC reading, read within its stratum; without them
        the distinct values of the trait are the categories.
    scorer : optional function of the weights of the categories (their counts
        in a sample, coded as the categories are) returning the value of
        every category; with it the values are recomputed from the readings
        of every bootstrap draw. Without it every category keeps its value.

    Returns
    -------
    ShareResult
    """
    for name, value, minimum in (("folds", folds, 2), ("repeats", repeats, 1),
                                  ("n_boot", n_boot, 0), ("n_perm", n_perm, 0)):
        if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
            raise ValueError(f"{name} must be an integer >= {minimum}; at least two folds are required" if name == "folds" else f"{name} must be an integer >= {minimum}")
    rng = _rng(seed)
    X = np.asarray(X, dtype=float)
    if X.ndim == 2 and 1 in X.shape:
        X = X.ravel()
    if X.ndim != 1:
        raise ValueError("the share is computed for one trait at a time; X must be one column")
    lineage = np.asarray(lineage, dtype=object).ravel()
    if X.size != lineage.size:
        raise ValueError(f"X has {X.size} rows and lineage has "
                         f"{lineage.size} entries; they must match")
    extras: typing.Dict[str, typing.Any] = {}
    for name, given in (("strata", strata), ("units", units), ("categories", categories)):
        vector = None if given is None else np.asarray(given, dtype=object).ravel()
        if vector is not None and vector.size != lineage.size:
            raise ValueError(f"{name} has {vector.size} entries and lineage has "
                             f"{lineage.size}; they must match")
        extras[name] = vector
    missing_share = _missing_share(lineage) if lineage.size else float("nan")
    finite = np.isfinite(X)
    n_dropped_non_finite = int((~finite).sum())
    typed = np.array([not _is_untyped(v) for v in lineage], dtype=bool)
    n_dropped_untyped = int((~typed).sum())
    keep = finite & typed
    X, lineage = X[keep], lineage[keep]
    extras = {k: (None if v is None else v[keep]) for k, v in extras.items()}
    # Singleton lineages are set aside before anything is scored; their count
    # and the support describe the typed collection.
    typed_sizes = np.bincount(_codes(lineage)) if lineage.size else np.zeros(0, dtype=int)
    typed_code = _codes(lineage) if lineage.size else np.zeros(0, dtype=int)
    repeated = typed_sizes[typed_code] >= 2 if lineage.size else np.zeros(0, dtype=bool)
    support = float(repeated.mean()) if lineage.size else float("nan")
    n_singletons_set_aside = int((~repeated).sum())
    n_retained = int(lineage.size)
    n_groups_all = int(typed_sizes.size)
    X, lineage = X[repeated], lineage[repeated]
    extras = {k: (None if v is None else v[repeated]) for k, v in extras.items()}
    code = _codes(lineage)
    n = X.size
    stratum = (_codes(np.asarray([str(v) for v in extras["strata"]], dtype=object))
               if extras["strata"] is not None and n else np.zeros(n, dtype=np.int64))
    unit = (_codes(np.asarray([str(v) for v in extras["units"]], dtype=object))
            if extras["units"] is not None and n else np.zeros(n, dtype=np.int64))
    trait = X
    if scorer is None:
        # a category is a reading within its stratum, and keeps its value
        raw = extras["categories"] if extras["categories"] is not None else X
        cat_labels = np.asarray([f"{s}|{c!r}" for s, c in zip(stratum, raw)], dtype=object)
        categories = _codes(cat_labels) if n else np.zeros(0, dtype=np.int64)
        K = int(categories.max()) + 1 if n else 0
        values = np.zeros(K)
        values[categories] = X
        if extras["strata"] is None:
            def scorer(_weights, values=values):
                return values
        else:
            # scored within strata: every value centred on the mean of its
            # stratum under the weights, a call on the prevalence of its
            # stratum; recomputed in every bootstrap draw, as the
            # mid-distribution scores of the readings are
            of_category = np.zeros(K, dtype=np.int64)
            of_category[categories] = stratum
            n_strata_scored = int(stratum.max()) + 1 if n else 0

            def scorer(weights, values=values, of_category=of_category, n_strata=n_strata_scored):
                w = np.asarray(weights, dtype=float)
                scores = values.copy()
                for s in range(n_strata):
                    idx = of_category == s
                    total = w[idx].sum()
                    if total > 0:
                        scores[idx] = values[idx] - (w[idx] @ values[idx]) / total
                return scores
            X = scorer(np.bincount(categories, minlength=K).astype(float))[categories]
    else:
        # the categories are the scorer's, numbered from zero across strata
        if extras["categories"] is None:
            raise ValueError("a scorer needs the category of every value")
        categories = np.asarray(extras["categories"], dtype=np.int64)
        K = int(categories.max()) + 1 if n else 0
        if n:
            X = scorer(np.bincount(categories, minlength=K).astype(float))[categories]
    cat_stratum = np.zeros(K, dtype=np.int64)
    cat_stratum[categories] = stratum
    permute_by = None if extras["strata"] is None else stratum
    test_by = (permute_by if extras["units"] is None else
               _codes(np.asarray([f"{s}|{u}" for s, u in zip(stratum, unit)], dtype=object)))

    def unavailable(**more):
        nan = float("nan")
        return ShareResult(kappa=nan, kappa_adj=nan, null_mean=nan, p_value=nan, n=n_retained,
                           n_groups=n_groups_all, prevalence=float(trait.mean()) if trait.size else nan,
                           support=support, missing_share=missing_share, estimable=False,
                           cv_sd=nan, n_dropped_non_finite=n_dropped_non_finite,
                           n_dropped_untyped=n_dropped_untyped,
                           n_singletons_set_aside=n_singletons_set_aside, n_scored=int(n), **more)
    if n < 2:
        # Fewer than two isolates leave nothing to hold out and nothing to
        # compare, so there is no share to estimate.
        return unavailable()
    sizes = np.bincount(code)
    shares = sizes[sizes > 0] / n
    structure: typing.Dict[str, typing.Any] = dict(
        n_groups_repeated=int((sizes >= 2).sum()),
        effective_groups=float(1.0 / np.sum(np.sort(shares) ** 2)),
        largest_group_share=float(shares.max()),
        n_positive_repeated=(float(trait.sum()) if np.isin(trait, (0.0, 1.0)).all() else float("nan")))
    if sizes.size < 2:
        # One lineage leaves no contrast between lineages, so the share is not
        # defined on this collection; returning 0 would read as "lineage explains
        # nothing" when nothing was compared.
        return unavailable(**structure)
    if float(np.ptp(X)) <= 0:
        # A trait alike in every isolate, or alike within every stratum, has
        # no variance to attribute.
        return unavailable(**structure)
    Xc = X.reshape(-1, 1)
    fold_draws: typing.List[np.ndarray] = []
    kappa, cv_sd, per_draw = _repeated_skill(Xc, code, folds=folds, repeats=repeats, rng=rng,
                                             keep_folds=fold_draws)
    if permute_by is None:
        perms = [rng.permutation(code) for _ in range(n_perm)]
    else:
        perms = [_permute_within(code, permute_by, rng) for _ in range(n_perm)]
    null = np.array([_skill(Xc, perm, _dealt_folds(perm, folds, rng)) for perm in perms],
                    dtype=float)
    # The control is the mean score of the permuted labellings, one fold draw
    # each.
    c_null = float(np.nanmean(null)) if null.size else float("nan")
    kappa_adj = _debias(kappa, c_null)
    reason = ""
    if np.isfinite(c_null) and c_null >= 1.0:
        reason = ("permuted null mean is at least one; the null correction denominator "
                  "is nonpositive, so the corrected share is not estimable")
    # The interval, on a stream of its own jumped far ahead, so that the
    # estimate and the test do not depend on whether it is computed.
    obs_lo = obs_hi = obs_se = float("nan")
    used = 0
    total = float(np.sum((X - X.mean()) ** 2))
    ratio = float(np.mean([_skill_parts(Xc, code, f)[1] for f in fold_draws])) / total
    if n_boot and np.isfinite(kappa_adj):
        boot_rng = np.random.Generator(_jumped(rng))
        finite_draws = per_draw[np.isfinite(per_draw)]
        fold_sd = (float(np.std(finite_draws, ddof=1)) / np.sqrt(finite_draws.size) / (1. - c_null)
                   if finite_draws.size > 1 else 0.)
        # the variance of the share where the lineages do not differ: the
        # fold-averaged share of the permuted labellings of the control
        permuted = np.array([_fold_average(X, perm, ratio) for perm in perms], dtype=float)
        finite = permuted[np.isfinite(permuted)]
        null_var = float(np.var(finite, ddof=1)) / (1. - c_null) ** 2 if finite.size > 1 else 0.
        t, se, t_perm = _smoothed_bootstrap(categories, cat_stratum, scorer, code, stratum,
                                            n_boot=n_boot, c_null=c_null, ratio=ratio,
                                            fold_sd=fold_sd, rng=boot_rng, null_var=null_var,
                                            perms=perms, permuted=permuted)
        t = np.asarray(t, dtype=float)
        used = int(t.size)
        if used >= MIN_USABLE_DRAWS * n_boot and used >= 2:
            q025, q975 = np.quantile(t, [.025, .975])
            if len(t_perm) >= 2:
                # neither end is read closer to the estimate than the
                # studentized permutation distribution allows: an interval
                # that leaves out zero is then one the studentized
                # permutation test would reject, and where the lineages are
                # small the bootstrap scenario understates how far the share
                # can be misread
                p025, p975 = np.quantile(np.asarray(t_perm, dtype=float), [.025, .975])
                q025, q975 = min(q025, p025), max(q975, p975)
            # the upper end also answers for a share accounted for by one lineage
            # read close to the rest (_single_lineage_gain); both ends within
            # 0 and 1, so that an interval wholly below zero is reported as
            # [0, 0] and the ends never cross
            gain = _single_lineage_gain(scorer(np.bincount(categories, minlength=K).astype(float)),
                                        code, categories, ratio, c_null)
            high = kappa_adj - se * q025
            if np.isfinite(gain):
                high = max(high, kappa_adj + gain)
            obs_lo = float(np.clip(kappa_adj - se * q975, 0., 1.))
            obs_hi = float(np.clip(high, 0., 1.))
            obs_se = float(se)
        else:
            reason = reason or (f"{used} of {n_boot} bootstrap draws gave a share; the "
                                "interval is not read from fewer than "
                                f"{MIN_USABLE_DRAWS:.0%} of them")
    # The test: the between-lineage sum of squares against its values under
    # the same permuted labellings; exact under exchangeability within
    # strata, and a function of the labels alone, so it carries no fold
    # noise. A permuted value that equals the observed one up to the rounding
    # of the sums counts as reaching it.
    observed = _between_lineage_ss(Xc, code)
    if test_by is not None and test_by is not permute_by:
        # within units: a set of permutations of its own, on the stream of
        # the interval's, so the estimate and its control are unchanged
        unit_rng = np.random.Generator(_jumped(np.random.Generator(_jumped(rng))))
        perms = [_permute_within(code, test_by, unit_rng) for _ in range(n_perm)]
    null_ss = np.array([_between_lineage_ss(Xc, perm) for perm in perms])
    exceed = int((null_ss >= observed - 1e-11 * abs(observed)).sum())
    return ShareResult(
        kappa=kappa, kappa_adj=kappa_adj, null_mean=c_null,
        observed_low=obs_lo, observed_high=obs_hi, observed_se=obs_se,
        p_value=_phipson_smyth(exceed, n_perm),
        n=n_retained, n_groups=n_groups_all, prevalence=float(trait.mean()),
        support=support, missing_share=missing_share,
        estimable=bool(np.isfinite(kappa_adj) and sizes.size >= MIN_REPEATED_LINEAGES),
        cv_sd=cv_sd, n_dropped_non_finite=n_dropped_non_finite,
        n_dropped_untyped=n_dropped_untyped,
        n_singletons_set_aside=n_singletons_set_aside, n_scored=int(n),
        p_floor=(1.0 / (n_perm + 1.0)) if n_perm else float("nan"),
        n_boot_used=used, fold_ratio=ratio,
        n_strata=int(test_by.max()) + 1 if test_by is not None else 0,
        reason=reason, **structure)


def clonal_share(y, lineage, **kwargs) -> ShareResult:
    """The lineage share of one call: :func:`layer_clonal_share` for a
    binary trait, the quantity a surveillance table reports."""
    return layer_clonal_share(np.asarray(y, dtype=float).ravel(), lineage, **kwargs)
