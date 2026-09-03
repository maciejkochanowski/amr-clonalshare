"""stats.py — statistical and distance primitives for amr-clonalshare.

Contents
--------
``binary_distance_matrix``
    Jaccard / Dice / simple-matching / Hamming distances for binary matrices,
    with the **empty-union case made explicit** rather than decided by a
    numerical convention. See below.
``jaccard_distance_matrix``
    Backwards-compatible wrapper.
``uninformative_rows``
    Which isolates carry none of a layer's retained features.
``effective_dimension``
    Participation ratio of a correlation eigenspectrum — how many independent
    features a block of collinear columns is really worth.
``fisher_exact_p``, ``benjamini_hochberg``
    Two-sided Fisher exact test and the Benjamini-Hochberg step-up.
``permutation_pvalue``
    Phipson-Smyth corrected permutation p-value.

The empty-union problem
-----------------------
Jaccard similarity ``|A & B| / |A | B|`` is undefined when both isolates carry
none of the features in a layer. The textbook convention sets it to 1, i.e.
declares two isolates that share *no detected determinant* to be maximally
similar. On accessory-genome presence/absence data that convention is not
harmless: absence is the default state, it is confounded with assembly quality,
and trait-absent isolates are usually the majority. In the *Klebsiella* example
shipped with this package, 1004 of 1500 isolates have an all-zero virulence row
after gating; under the ``"identical"`` convention they form a single
503,506-pair clique at distance 0, and the indicator "is this row all zero"
alone reproduces the pipeline's headline partition at ARI 0.83. The largest
discovered "archetype" is then the trait-absent stratum wearing a cluster's
clothes.

This module therefore requires the caller to choose, and records the choice:

``undefined_pair="identical"``
    ``D = 0``. The classical convention. Reproduces prior behaviour.
``undefined_pair="distinct"``
    ``D = 1``. Trait-absent isolates are not merged, at the cost of making them
    mutually maximally distant.
``undefined_pair="nan"``
    ``D = NaN``. For callers that handle missingness explicitly.

Whichever is chosen, :func:`uninformative_rows` should be reported alongside
the partition so a reader can see how much of a cluster is the empty stratum.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy import stats

__all__ = [
    "binary_distance_matrix",
    "jaccard_distance_matrix",
    "uninformative_rows",
    "effective_dimension",
    "fisher_exact_p",
    "benjamini_hochberg",
    "permutation_pvalue",
]

_UNDEFINED = ("identical", "distinct", "nan")
_METRICS = ("jaccard", "dice", "simple_matching", "hamming")


def uninformative_rows(X: np.ndarray) -> np.ndarray:
    """Boolean mask of rows carrying none of the layer's features."""
    return np.asarray(X).astype(bool).sum(axis=1) == 0


def binary_distance_matrix(X: np.ndarray, *, metric: str = "jaccard",
                           undefined_pair: str = "identical") -> np.ndarray:
    """Pairwise distance matrix for a binary (n, p) matrix.

    Parameters
    ----------
    X : (n, p) binary matrix.
    metric :
        ``"jaccard"``          ``1 - a / (a + b + c)``   (presence-only)
        ``"dice"``             ``1 - 2a / (2a + b + c)`` (presence-only)
        ``"simple_matching"``  ``1 - (a + d) / p``       (counts shared absence)
        ``"hamming"``          ``(b + c) / p``           (= simple matching)
        with ``a = |A & B|``, ``b, c`` the exclusive counts and ``d`` the shared
        absences.
    undefined_pair : how to score pairs whose union is empty; only the
        presence-only metrics have such pairs. One of ``"identical"`` (D=0,
        classical), ``"distinct"`` (D=1), ``"nan"``.

    Returns
    -------
    (n, n) symmetric distance matrix with a zero diagonal.
    """
    if metric not in _METRICS:
        raise ValueError(f"metric must be one of {_METRICS}; got {metric!r}")
    if undefined_pair not in _UNDEFINED:
        raise ValueError(
            f"undefined_pair must be one of {_UNDEFINED}; got {undefined_pair!r}")

    Xb = np.asarray(X).astype(np.float64)
    n, p = Xb.shape
    a = Xb @ Xb.T                                   # shared presences
    pop = Xb.sum(axis=1)
    union = pop[:, None] + pop[None, :] - a         # b + c + a

    if metric in ("simple_matching", "hamming"):
        # b + c = pop_i + pop_j - 2a ; distance = (b + c) / p
        D = (pop[:, None] + pop[None, :] - 2.0 * a) / max(p, 1)
        np.fill_diagonal(D, 0.0)
        return (D + D.T) / 2.0

    with np.errstate(divide="ignore", invalid="ignore"):
        if metric == "jaccard":
            sim = np.where(union > 0, a / np.where(union > 0, union, 1.0), np.nan)
        else:  # dice
            denom = pop[:, None] + pop[None, :]
            sim = np.where(denom > 0, 2.0 * a / np.where(denom > 0, denom, 1.0), np.nan)
        D = 1.0 - sim

    if undefined_pair == "identical":
        D = np.where(np.isnan(D), 0.0, D)
    elif undefined_pair == "distinct":
        D = np.where(np.isnan(D), 1.0, D)
    # "nan" -> leave as NaN

    np.fill_diagonal(D, 0.0)
    return (D + D.T) / 2.0


def jaccard_distance_matrix(X: np.ndarray, *,
                            undefined_pair: str = "identical") -> np.ndarray:
    """Jaccard distance matrix. Thin wrapper over :func:`binary_distance_matrix`.

    The default ``undefined_pair="identical"`` reproduces the classical
    convention; see the module docstring for why that choice is consequential
    and must be reported.
    """
    return binary_distance_matrix(X, metric="jaccard",
                                  undefined_pair=undefined_pair)


def effective_dimension(X: np.ndarray) -> float:
    """Participation ratio of the correlation eigenspectrum of ``X``.

    ``(sum lambda_i)^2 / sum lambda_i^2`` for the eigenvalues of the column
    correlation matrix: the number of mutually uncorrelated columns the block
    is worth. A set of ``p`` perfectly collinear columns has effective
    dimension 1 however large ``p`` is. Reported so that a layer made of
    co-inherited operon members is not treated as if it carried ``p``
    independent hypotheses.

    Standard participation-ratio / effective-rank construction; see e.g. Roy &
    Vetterli (2007), "The effective rank: a measure of effective
    dimensionality", EUSIPCO.
    """
    A = np.asarray(X, dtype=float)
    if A.ndim != 2 or A.shape[1] == 0:
        return 0.0
    keep = A.std(axis=0) > 0
    if keep.sum() == 0:
        return 0.0
    C = np.corrcoef(A[:, keep], rowvar=False)
    C = np.atleast_2d(C)
    lam = np.linalg.eigvalsh(C)
    lam = np.clip(lam, 0.0, None)
    denom = float(np.sum(lam ** 2))
    if denom <= 0:
        return 0.0
    return float(np.sum(lam) ** 2 / denom)


def fisher_exact_p(n11, n10, n01, n00) -> float:
    """Two-sided Fisher exact p-value for a 2x2 table."""
    return float(stats.fisher_exact([[n11, n10], [n01, n00]],
                                    alternative="two-sided")[1])


def benjamini_hochberg(pvals, q: float = 0.05, *,
                       nan_policy: str = "raise",
                       dependence: str = "independent") -> Tuple[np.ndarray, np.ndarray]:
    """Benjamini-Hochberg step-up. Returns ``(adjusted_p, reject_mask)``.

    Benjamini & Hochberg (1995), JRSS-B 57:289-300,
    doi:10.1111/j.2517-6161.1995.tb02031.x. The rejection mask is
    ``adjusted_p <= q``; callers must use this mask rather than re-thresholding,
    so that a p-value exactly equal to ``q`` is treated consistently.

    A non-finite p-value is refused. ``np.argsort`` puts NaN at the largest
    ranks and the step-up runs a running minimum down from there, so a single
    NaN turns the whole family into NaN, every comparison against ``q`` into
    ``False``, and a family holding real discoveries into "nothing significant"
    with nothing printed. That is the same output a true null gives, which is
    why it stops the run instead: ``nan_policy="raise"`` by default, or
    ``"omit"`` to drop the non-finite entries from the family, which shrinks
    ``m`` and is a claim the caller has to make deliberately.

    ``dependence`` selects the assumption the family is willing to make about
    its own p-values. ``"independent"`` is the 1995 step-up and is valid under
    independence or positive regression dependence. ``"arbitrary"`` is the
    Benjamini-Yekutieli step-up, which runs the same procedure at ``q`` divided
    by the harmonic sum of the family size and is valid whatever the dependence
    (Benjamini & Yekutieli 2001, Annals of Statistics 29(4),
    doi:10.1214/aos/1013699998). A panel of antimicrobials is the second case:
    cross-resistance makes agents sharing a target site rise together, while a
    resistance trade-off can make a pair move apart, so the sign of the
    dependence is not known in advance and cannot be assumed positive.
    """
    pvals = np.asarray(pvals, dtype=float)
    m = len(pvals)
    if m == 0:
        return np.array([]), np.array([], dtype=bool)
    bad = ~np.isfinite(pvals)
    n_bad = int(bad.sum())
    if n_bad:
        if nan_policy == "raise":
            raise ValueError(
                f"{n_bad} of {m} p-values are not finite (indices "
                f"{np.flatnonzero(bad)[:10].tolist()}); Benjamini-Hochberg "
                f"would report no discoveries for the whole family. Fix the "
                f"tests that produced them, or pass nan_policy='omit' and "
                f"report how many were dropped")
        if nan_policy != "omit":
            raise ValueError(f"nan_policy must be 'raise' or 'omit', "
                             f"not {nan_policy!r}")
    if dependence not in ("independent", "arbitrary"):
        raise ValueError(f"dependence must be 'independent' or 'arbitrary', "
                         f"not {dependence!r}")
    adj = np.full(m, np.nan)
    keep = np.flatnonzero(~bad)
    k = len(keep)
    if k:
        good = pvals[keep]
        order = np.argsort(good)
        ranked = good[order]
        scale = (float(np.sum(1.0 / np.arange(1, k + 1)))
                 if dependence == "arbitrary" else 1.0)
        adj_sorted = np.minimum.accumulate(
            (ranked * scale * k / np.arange(1, k + 1))[::-1])[::-1]
        adj_sorted = np.minimum(adj_sorted, 1.0)
        adj_good = np.empty(k)
        adj_good[order] = adj_sorted
        adj[keep] = adj_good
    with np.errstate(invalid="ignore"):
        reject = adj <= q
    return adj, np.where(np.isnan(adj), False, reject)


def permutation_pvalue(null_stats, observed: float, *,
                       tail: str = "greater") -> float:
    """Phipson-Smyth corrected permutation p-value: ``(b + 1) / (B + 1)``.

    A p-value estimated from ``B`` randomly drawn permutations can never
    legitimately be 0; the unbiased and valid estimator adds one to numerator
    and denominator. Phipson & Smyth (2010), "Permutation p-values should never
    be zero", Stat Appl Genet Mol Biol 9(1):39, doi:10.2202/1544-6115.1585.
    """
    s = np.asarray(list(null_stats), dtype=float)
    B = s.size
    if B == 0:
        return 1.0
    if tail == "greater":
        b = int((s >= observed).sum())
    elif tail == "less":
        b = int((s <= observed).sum())
    else:
        raise ValueError(f"tail must be 'greater' or 'less'; got {tail!r}")
    return float((b + 1) / (B + 1))
