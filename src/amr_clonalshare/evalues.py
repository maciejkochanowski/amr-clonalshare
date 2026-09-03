"""Anytime-valid evidence for a lineage effect, and e-BH across a panel.

WHY THIS EXISTS. Antimicrobial resistance surveillance is sequential by
construction: the same panel is re-analysed every time a year of isolates
arrives. Recomputing a Benjamini-Hochberg procedure on the accumulated data at
each look has no error guarantee, because the number of looks is not fixed in
advance and the looks are not independent. An e-process does have one. Its
expectation under the null is at most one however much data has accrued and
however many times it has been inspected, so a surveillance programme may
inspect the running value whenever it likes and stop when it likes.

NAMING. "e-value" here is the betting or martingale sense of Vovk and Wang
(2021) and Wang and Ramdas (2022). It is not the BLAST expectation value, and
it is not the E-value of VanderWeele and Ding for sensitivity to unmeasured
confounding. Both of those are common in this literature and neither is meant.
The names ``e_process``, ``test martingale`` and ``anytime_valid`` are used in
preference to the bare word wherever the code allows it. One call of
:func:`e_process` on one cohort returns one e-value, the first term of the
process; the process itself is the product that :func:`combine_independent`
forms as batches arrive.

CONSTRUCTION. The e-value is the split likelihood ratio of Wasserman, Ramdas
and Balakrishnan (2020), often called universal inference: fit the alternative
on a training fold, evaluate its likelihood on the held-out fold, and divide by
the null likelihood maximised on that same held-out fold. Because the numerator
parameters never saw the held-out rows, the ratio has expectation at most one
under the null with no regularity conditions on the model. That is exactly the
fold structure the clonal share already uses, so the evidence costs one extra
pass over data that has already been split.

Group probabilities are shrunk toward the training grand mean by the same
empirical-Bayes factor the clonal share uses. Shrinkage is not a convenience
here: an unshrunk group of one training isolate gives a probability of zero or
one, the log likelihood diverges, and the e-value becomes infinite for a reason
that has nothing to do with lineage.

COMBINING. Within one cohort the folds share data, so their e-values are
averaged, which is valid because any convex combination of e-values is an
e-value. Across independent batches, for instance successive years of
collection, they are multiplied, which forms a test martingale and is what
makes the procedure anytime-valid. Across a panel of antimicrobials the e-BH
procedure of Wang and Ramdas (2022) controls the false discovery rate under
arbitrary dependence between agents, which matters because cross-resistance
makes a macrolide block behave as one trait.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Dict, Iterable, List, Sequence

import numpy as np

from .attribution import _codes, _folds, _rng

__all__ = ["EResult", "e_process", "e_bh", "combine_independent",
           "combine_within_cohort"]

# Ville's inequality gives P(sup_t E_t >= 1/alpha) <= alpha, so an e-value of
# at least 1/alpha is the anytime-valid analogue of a p-value below alpha.
REJECT_AT = {0.05: 20.0, 0.01: 100.0}


@dataclass(frozen=True)
class EResult:
    """Evidence against the hypothesis that a trait is independent of lineage.

    e_value : the averaged split likelihood ratio. One is no evidence; larger
        is more. Under the null its expectation is at most one.
    log_e : the same on the log scale, which is the number that stays finite
        when evidence accumulates over many years of surveillance.
    reject_05, reject_01 : whether the value clears 1/alpha, the counterpart
        of a significance threshold by Ville's inequality. For one cohort this
        is a single e-value and a single decision. The anytime property
        belongs to the running product across independent batches formed by
        :func:`combine_independent`: that product may be inspected after any
        batch and as often as wanted without inflating the error rate. Adding
        isolates to the same cohort and recomputing is not a new batch.
    n, n_groups, prevalence : the design the evidence came from.
    n_splits : how many training and held-out splits were averaged.
    """

    e_value: float
    log_e: float
    reject_05: bool
    reject_01: bool
    n: int
    n_groups: int
    prevalence: float
    n_splits: int

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def _shrunk_group_p(y_tr: np.ndarray, code_tr: np.ndarray, G: int):
    """Empirical-Bayes group probabilities from the training rows.

    The shrinkage factor is the same one the clonal share uses,
    ``tau2 / (tau2 + sigma2 / n_g)``, with both components from a one-way
    moment estimator on the training rows only. A group absent from training
    falls back to the training grand mean, which is what an honest predictor
    would do.
    """
    cnt = np.bincount(code_tr, minlength=G).astype(float)
    present = cnt > 0
    tot = np.bincount(code_tr, weights=y_tr, minlength=G)
    gm = np.divide(tot, cnt, out=np.zeros(G), where=present)
    grand = float(y_tr.mean())
    n, gp = len(y_tr), int(present.sum())
    if gp < 2 or n - gp < 1:
        return np.full(G, grand), grand
    msw = float(((y_tr - gm[code_tr]) ** 2).sum() / (n - gp))
    msb = float((cnt[present] * (gm[present] - grand) ** 2).sum() / (gp - 1))
    n0 = (n - (cnt[present] ** 2).sum() / n) / (gp - 1)
    tau2 = max(0.0, (msb - msw) / max(n0, 1e-12))
    denom = tau2 + msw / np.maximum(cnt, 1.0)
    b = np.where(present & (denom > 0.0),
                 tau2 / np.where(denom > 0.0, denom, 1.0), 0.0)
    p = grand + b * (gm - grand)
    p = np.where(present, p, grand)
    return np.clip(p, 1e-6, 1 - 1e-6), grand


def _log_bernoulli(y: np.ndarray, p: np.ndarray) -> float:
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return float((y * np.log(p) + (1.0 - y) * np.log1p(-p)).sum())


def e_process(y, lineage, *, folds: int = 5, repeats: int = 20,
              seed=None) -> EResult:
    """Anytime-valid evidence that ``y`` depends on ``lineage``.

    The split likelihood ratio is computed once per held-out fold: group
    probabilities are fitted on the training rows and scored on the held-out
    rows, against the null probability maximised on those same held-out rows.
    Fold values are averaged, and repeats are averaged, because folds drawn
    from one cohort share data and a convex combination of e-values is an
    e-value while a product of dependent ones is not.

    Returns one for a trait carrying no lineage information, and grows without
    bound as evidence accumulates. It never needs a null distribution, a
    permutation, or an asymptotic approximation.
    """
    y = np.asarray(y, dtype=float).ravel()
    code = _codes(lineage)
    n = y.size
    if code.shape[0] != n:
        raise ValueError(f"lineage has {code.shape[0]} entries for {n} traits")
    rng = _rng(seed)
    G = int(code.max()) + 1
    prevalence = float(np.nanmean(y))

    per_split: List[float] = []
    for _ in range(int(repeats)):
        fold = _folds(n, folds, rng)
        for f in np.unique(fold):
            tr, te = fold != f, fold == f
            if not tr.any() or not te.any():
                continue
            p_alt, _ = _shrunk_group_p(y[tr], code[tr], G)
            # The null is given its best fit on the held-out rows themselves,
            # so the ratio cannot be won by the alternative merely knowing the
            # overall prevalence better.
            p_null = float(np.clip(y[te].mean(), 1e-6, 1 - 1e-6))
            ll_alt = _log_bernoulli(y[te], p_alt[code[te]])
            ll_null = _log_bernoulli(y[te], np.full(int(te.sum()), p_null))
            per_split.append(ll_alt - ll_null)

    if not per_split:
        return EResult(float("nan"), float("nan"), False, False, n, G,
                       prevalence, 0)
    logs = np.asarray(per_split, dtype=float)
    # Average on the natural scale, in a way that does not overflow.
    m = float(logs.max())
    e = float(np.exp(m) * np.exp(logs - m).mean())
    return EResult(e_value=e, log_e=float(np.log(e)) if e > 0 else float("-inf"),
                   reject_05=bool(e >= REJECT_AT[0.05]),
                   reject_01=bool(e >= REJECT_AT[0.01]),
                   n=n, n_groups=G, prevalence=prevalence,
                   n_splits=len(per_split))


def combine_independent(values: Iterable[float]) -> float:
    """Multiply e-values from independent batches, for example successive
    years of collection. The running product is a test martingale, which is
    what makes it legitimate to inspect the total after every batch and stop
    at any point without spending a significance budget."""
    out = 1.0
    for v in values:
        if not np.isfinite(v) or v < 0:
            raise ValueError(f"not an e-value: {v}")
        out *= float(v)
    return out


def combine_within_cohort(values: Iterable[float]) -> float:
    """Average e-values computed on overlapping splits of one cohort. A convex
    combination of e-values is an e-value; a product of dependent ones is
    not."""
    v = np.asarray(list(values), dtype=float)
    if v.size == 0 or np.any(v < 0):
        raise ValueError("e-values must be non-negative and non-empty")
    return float(v.mean())


def e_bh(evalues: Sequence[float], alpha: float = 0.05) -> Dict[str, object]:
    """The e-BH procedure of Wang and Ramdas (2022).

    Sort the e-values in decreasing order and reject the largest ``k`` for
    which the ``k``-th largest is at least ``m / (alpha * k)``. This controls
    the false discovery rate at ``alpha`` under *arbitrary* dependence between
    the hypotheses, with no correction. That property is what makes it the
    right instrument for an antimicrobial panel, where cross-resistance means
    the agents are neither independent nor reliably positively dependent, so
    the usual justification for Benjamini-Hochberg does not apply.

    Returns the indices rejected and the effective threshold.
    """
    e = np.asarray(list(evalues), dtype=float)
    m = e.size
    if m == 0:
        return {"rejected": [], "n_rejected": 0, "threshold": float("inf"),
                "alpha": alpha, "m": 0}
    finite = np.where(np.isfinite(e), e, 0.0)
    order = np.argsort(-finite)
    sorted_e = finite[order]
    k_star = 0
    for k in range(1, m + 1):
        if sorted_e[k - 1] >= m / (alpha * k):
            k_star = k
    rejected = sorted(int(i) for i in order[:k_star])
    thr = m / (alpha * k_star) if k_star else float("inf")
    return {"rejected": rejected, "n_rejected": k_star, "threshold": float(thr),
            "alpha": float(alpha), "m": int(m)}
