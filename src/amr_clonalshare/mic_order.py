"""Lineage share of the recorded MIC ordering.

Every reading is replaced by its mid-distribution score within its stratum,
normally the testing laboratory: the share of that stratum's readings that lie
below it plus half the share that fall in the same reading. For interval-
censored readings the distribution is the nonparametric maximum-likelihood
estimate of the stratum (Turnbull 1976), and the score is then the Wilcoxon
score of Peto and Peto (1972) and Prentice (1978) for censored data; for exact
readings, and for a stratum read on one panel, it is the mid-rank less one
half, divided by the number of readings. The share of these scores accounted
for by lineage is estimated by
:func:`amr_clonalshare.attribution.layer_clonal_share`, the estimator of the
binary share, with the lineage labels permuted only within strata; the
p-value compares the between-lineage sum of squares of the scores with its
values under those permutations, the permutation form of the Kruskal-Wallis
test on the scores within strata. The score
is the mid-distribution transform on which the rank intraclass correlation of
Tu, Li, Zeng and Shepherd (2023) is defined, so the share is the analogue of
that correlation for the observed lineages, extended to censored readings.

What this buys, and what it does not:

* No distributional form is assumed. A reading with one mode or two, a heavy
  tail or outlying wells gives a score all the same, and the permutation test
  is exact under exchangeability within strata whatever that form is.
* The score is invariant to any increasing transformation of the MIC scale,
  so mg/L and log2 give the same share.
* Scored within the laboratory, a constant offset between laboratories drops
  out, and a difference in their panels does not appear as a difference
  between lineages, though a coarser panel resolves less of its laboratory's
  ordering. A difference between lineages that were read only in different
  laboratories drops out as well: the share is the within-laboratory
  ordering.
* With a single cut point the score is 1/2 + (call - prevalence)/2 within
  each stratum, so with one stratum the share equals the share of the binary
  call exactly, and with several it is the share of the call centred on its
  prevalence within each stratum.
* The share describes the ordering the panel resolves. Readings piled on an
  end well are tied, so a panel with most readings at its ends shows less of
  the lineage ordering than there is, as a call shows less than the dilution
  it was cut from. ``share_end_wells`` is reported beside the share for that
  reason.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np

from .attribution import layer_clonal_share

__all__ = ["order_scores", "mic_order_share"]

#: Iterations of the constrained Newton method, and the tolerance on its
#: optimality condition: no cell's directional derivative of the mean log
#: likelihood exceeds one by more than TOLERANCE, which bounds the gap to the
#: maximum.
MAX_ITERATIONS = 500
TOLERANCE = 1e-9


def _cell_ranges(lo: np.ndarray, hi: np.ndarray):
    """The contiguous run of elementary cells each reading covers.

    With the finite endpoints E_0 < ... < E_{m-1}, the cells in order are the
    gap below E_0, the point E_0, the gap (E_0, E_1), ..., the point E_{m-1}
    and the gap above it. A reading (a, b] covers the gaps and points between
    a (excluded) and b (included); an exact reading covers its point.
    """
    ends = np.unique(np.r_[lo[np.isfinite(lo)], hi[np.isfinite(hi)]])
    m = ends.size
    ia = np.where(np.isfinite(lo), np.searchsorted(ends, lo), -1)
    ib = np.where(np.isfinite(hi), np.searchsorted(ends, hi), m)
    exact = lo == hi
    first = np.where(exact, 2 * ia + 1, 2 * ia + 2)
    last = np.where(exact, 2 * ia + 1, np.where(np.isfinite(hi), 2 * ib + 1, 2 * m))
    return first.astype(int), last.astype(int), 2 * m + 1


def _npmle(cover: np.ndarray, weight: np.ndarray) -> np.ndarray:
    """Masses of the cells at the nonparametric maximum-likelihood distribution.

    ``cover`` holds one row per distinct reading and one column per cell (1
    where the reading covers the cell), ``weight`` the share of the readings
    each row stands for. The mean log likelihood sum_i w_i log(cover_i . p) is
    maximised over the simplex by the constrained Newton method of Wang
    (2007): the cell of largest directional derivative joins the support, the
    masses on the support take the Newton step of the quadratic approximation,
    and the step is halved until the slope of the likelihood at its end is not
    negative. The step is that of
    sum_i w_i log(cover_i . q) - sum(q) - (sum(q) - 1)^2 / 2 over nonnegative
    q, which has the same maximum and the same optimality condition, so the
    step leaves the maximum where it is; its quadratic approximation is the
    nonnegative least-squares problem solved here, with one row of ones and
    target zero for the last two terms. The likelihood is concave along the
    step, so it cannot fall; and the slope, unlike a difference of two
    likelihoods, is resolved in floating point down to the tolerance.
    It stops when no cell's directional derivative exceeds one by more than
    the tolerance, the condition that characterises the maximum.
    """
    from scipy.optimize import nnls
    cells = cover.shape[1]
    # start from mass spread evenly over the cells some reading covers, so
    # that every reading has a positive likelihood
    p = cover.max(0)
    p /= p.sum()

    def slope(r, direction):
        # derivative of the mean log likelihood at r along a direction whose
        # entries sum to zero, written as sum_j (d_j - 1) direction_j
        mass = cover @ r
        if (mass <= 0).any():
            return -np.inf
        return float(((weight / mass) @ cover - 1.) @ direction)

    for _ in range(MAX_ITERATIONS):
        mass = cover @ p
        gradient = (weight / mass) @ cover
        if gradient.max() - 1. <= TOLERANCE:
            return p
        support = (p > 0) | (np.arange(cells) == np.argmax(gradient))
        scaled = cover[:, support] * (np.sqrt(weight) / mass)[:, None]
        q = np.zeros(cells)
        q[support], _ = nnls(np.vstack([scaled, np.ones((1, int(support.sum())))]),
                             np.r_[2. * np.sqrt(weight), 0.])
        if q.sum() <= 0:
            raise ArithmeticError("the Newton step of the nonparametric maximum-likelihood "
                                  "distribution returned no mass")
        q /= q.sum()
        step = 1.
        while step > 1e-8 and slope((1. - step) * p + step * q, q - p) < 0.:
            step /= 2.
        p = (1. - step) * p + step * q
        p[p < 1e-15] = 0.
        p /= p.sum()
    raise ArithmeticError("the nonparametric maximum-likelihood distribution of the "
                          "readings did not reach its optimality condition")


def _stratum_scores(lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    first, last, cells = _cell_ranges(lo, hi)
    pairs, inverse, counts = np.unique(np.column_stack([first, last]), axis=0,
                                       return_inverse=True, return_counts=True)
    inverse = inverse.ravel()
    if np.all(pairs[1:, 0] > pairs[:-1, 1]):
        # one panel: the readings do not overlap, and the distribution puts
        # on every reading its share of the readings
        share = counts / counts.sum()
        return (np.cumsum(share) - share / 2.)[inverse]
    cover = ((np.arange(cells)[None, :] >= pairs[:, :1]) & (np.arange(cells)[None, :] <= pairs[:, 1:])).astype(float)
    p = _npmle(cover, counts / counts.sum())
    cum = np.r_[0., np.cumsum(p)]
    below = cum[pairs[:, 0]]
    inside = cum[pairs[:, 1] + 1] - cum[pairs[:, 0]]
    return (below + inside / 2.)[inverse]


def order_scores(lo, hi, strata=None) -> np.ndarray:
    """Mid-distribution score of every reading within its stratum.

    ``lo`` and ``hi`` are the log2 bounds of each reading, (lo, hi], with
    -inf and inf at the ends of the panel and lo == hi for an exact reading.
    Readings without a value (NaN) or unbounded on both sides are not
    accepted; the caller removes them first.
    """
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    if lo.shape != hi.shape or lo.ndim != 1:
        raise ValueError("lo and hi must be one-dimensional and of one length")
    if np.isnan(lo).any() or np.isnan(hi).any():
        raise ValueError("readings without a value must be removed before scoring")
    if (np.isneginf(lo) & np.isposinf(hi)).any():
        raise ValueError("a reading unbounded on both sides carries no ordering")
    if (lo > hi).any():
        raise ValueError("every reading needs lo <= hi")
    out = np.empty(lo.size)
    if strata is None:
        if lo.size:
            out[:] = _stratum_scores(lo, hi)
        return out
    strata = np.asarray([str(s) for s in strata], dtype=object)
    if strata.size != lo.size:
        raise ValueError("strata must give one label per reading")
    for name in np.unique(strata):
        idx = strata == name
        out[idx] = _stratum_scores(lo[idx], hi[idx])
    return out


class _Scorer:
    """Mid-distribution scores of the readings of every stratum from the
    weights of the readings (counts in a sample, or a distribution): the categories
    are the distinct readings of every stratum, numbered from zero across the
    strata in the order of the strata and, within each, of their cells."""

    def __init__(self, lo, hi, stratum_labels):
        self.names = sorted(set(stratum_labels.tolist()), key=str)
        self.category = np.empty(lo.size, dtype=np.int64)
        self.parts = []
        start = 0
        for name in self.names:
            idx = np.flatnonzero(stratum_labels == name)
            first, last, cells = _cell_ranges(lo[idx], hi[idx])
            pairs, inverse = np.unique(np.column_stack([first, last]), axis=0, return_inverse=True)
            self.category[idx] = start + inverse.ravel()
            disjoint = bool(np.all(pairs[1:, 0] > pairs[:-1, 1]))
            cover = ((np.arange(cells)[None, :] >= pairs[:, :1])
                     & (np.arange(cells)[None, :] <= pairs[:, 1:])).astype(float)
            self.parts.append((slice(start, start + len(pairs)), pairs, cover, disjoint))
            start += len(pairs)
        self.K = start

    def __call__(self, weights):
        scores = np.zeros(self.K)
        for where, pairs, cover, disjoint in self.parts:
            w = np.asarray(weights[where], dtype=float)
            total = w.sum()
            if total <= 0:
                continue
            share = w / total
            if disjoint:
                scores[where] = np.cumsum(share) - share / 2.
                continue
            present = share > 0
            mass = _npmle(cover[present], share[present])
            cum = np.r_[0., np.cumsum(mass)]
            below = cum[pairs[:, 0]]
            scores[where] = below + (cum[pairs[:, 1] + 1] - below) / 2.
        return scores


def mic_order_share(lo, hi, lineage, *, strata=None, units=None,
                    stratified_by: Optional[str] = None,
                    folds: int = 5, repeats: int = 20, n_boot: int = 999,
                    n_perm: int = 999, seed=0) -> Dict[str, Any]:
    """The lineage share of the within-stratum MIC ordering, as a record.

    The readings of repeated, typed lineages are scored within their stratum
    at the distribution of those readings, and the share, its interval for
    the represented lineages and its permutation p-value are those of
    :func:`layer_clonal_share`, with the permutations kept within strata
    (and sampling units) and the scores recomputed in every bootstrap draw.
    The record adds what the share rests on: the strata, the readings per
    stratum, the distinct readings and the share of readings on an end well,
    and the bounds on the latent ordering (:mod:`.latent_order`).
    """
    from .attribution import _is_untyped
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    lineage = np.asarray(lineage, dtype=object)
    keep = ~(np.isnan(lo) | np.isnan(hi)) & ~(np.isneginf(lo) & np.isposinf(hi))
    lo, hi, lineage = lo[keep], hi[keep], lineage[keep]
    labels = (np.zeros(lo.size, dtype=object) if strata is None else
              np.asarray([str(s) for s in strata], dtype=object)[keep])
    units = None if units is None else np.asarray([str(u) for u in units], dtype=object)[keep]
    typed = np.array([not _is_untyped(v) for v in lineage], dtype=bool)
    names_all = np.asarray([str(v) for v in lineage], dtype=object)
    _, code, sizes = np.unique(names_all[typed], return_inverse=True, return_counts=True)
    repeated = np.zeros(lo.size, dtype=bool)
    repeated[np.flatnonzero(typed)] = sizes[code.ravel()] >= 2
    counts = dict(n=int(typed.sum()), n_groups=int(sizes.size),
                  support=float(repeated.sum() / typed.sum()) if typed.any() else float("nan"),
                  missing_share=float(np.mean(~typed)) if lo.size else float("nan"),
                  n_dropped_untyped=int((~typed).sum()),
                  n_singletons_set_aside=int((typed & ~repeated).sum()))
    lo_r, hi_r, lin_r, lab_r = lo[repeated], hi[repeated], lineage[repeated], labels[repeated]
    unit_r = None if units is None else units[repeated]
    if lo_r.size:
        scorer = _Scorer(lo_r, hi_r, lab_r)
        scores = scorer(np.bincount(scorer.category, minlength=scorer.K).astype(float))[scorer.category]
        categories = scorer.category
    else:
        scorer, scores, categories = None, np.zeros(0), None
    result = layer_clonal_share(scores, lin_r, folds=folds, repeats=repeats,
                                n_boot=n_boot, n_perm=n_perm, seed=seed,
                                strata=None if strata is None else lab_r, units=unit_r,
                                categories=categories, scorer=scorer)
    record = result.as_dict()
    record.pop("prevalence", None)
    record.pop("n_positive_repeated", None)
    record.update(counts)
    names = ["all"] if strata is None else sorted(set(labels))
    record.update(
        method="lineage share of the within-stratum mid-distribution score",
        stratified_by=stratified_by if strata is not None else None,
        strata={name: int(lo.size if strata is None else (labels == name).sum()) for name in names},
        distinct_readings={name: int(len(set(zip(
            (lo if strata is None else lo[labels == name]).tolist(),
            (hi if strata is None else hi[labels == name]).tolist())))) for name in names},
        share_end_wells=float(np.mean(np.isinf(lo) | np.isinf(hi))) if lo.size else float("nan"),
        n_readings=int(lo.size),
        sampling_units=None if units is None else int(len(set(units.tolist()))),
    )
    # The same readings bound the lineage share of the latent MIC ordering
    # (latent_order.py), on a stream of its own when the seed is a number, so
    # the share is unchanged by it.
    from .latent_order import order_bounds
    bounds_rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng([int(seed), 1])
    record.update(order_bounds(lo, hi, lineage, strata=None if strata is None else labels,
                               seed=bounds_rng))
    return record
