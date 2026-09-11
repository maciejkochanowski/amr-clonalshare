"""stats.py — statistical primitives for amr-clonalshare.

Contents
--------
``effective_dimension``
    Participation ratio of a correlation eigenspectrum: how many independent
    agents a block of correlated resistance columns is really worth.
``fisher_exact_p``, ``benjamini_hochberg``
    Two-sided Fisher exact test and the step-up, Benjamini-Hochberg under
    independence and Benjamini-Yekutieli under arbitrary dependence.
``permutation_pvalue``
    Phipson-Smyth corrected permutation p-value.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy import stats

__all__ = [
    "effective_dimension",
    "fisher_exact_p",
    "benjamini_hochberg",
    "permutation_pvalue",
]


def effective_dimension(X: np.ndarray) -> float:
    """Participation ratio of the correlation eigenspectrum of ``X``.

    ``(sum lambda_i)^2 / sum lambda_i^2`` for the eigenvalues of the column
    correlation matrix: the number of mutually uncorrelated columns the block
    is worth. A set of ``p`` perfectly collinear columns has effective
    dimension 1 however large ``p`` is. Reported so that a panel of
    co-selected agents is not treated as if it carried ``p`` independent
    hypotheses.

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
