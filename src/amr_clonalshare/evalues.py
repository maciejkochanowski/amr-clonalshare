"""Anytime-valid evidence for a lineage effect, and e-BH across a panel.

Why this module exists
----------------------
Antimicrobial resistance surveillance is sequential by construction: the same
panel is re-analysed every time a year of isolates arrives. Recomputing a
Benjamini-Hochberg procedure on the accumulated data at each look has no error
guarantee, because the number of looks is not fixed in advance and the looks
are not independent. An e-process does have one. Its expectation under the null
is at most one however much data has accrued and however many times it has been
inspected, so a surveillance programme may inspect the running value whenever
it likes and stop when it likes.

Naming
------
"e-value" here is the betting or martingale sense of Vovk and Wang
(2021) and Wang and Ramdas (2022). It is not the BLAST expectation value, and
it is not the E-value of VanderWeele and Ding for sensitivity to unmeasured
confounding. Both of those are common in this literature and neither is meant.
The names ``e_process``, ``test martingale`` and ``anytime_valid`` are used in
preference to the bare word wherever the code allows it. One call of
:func:`e_process` on one cohort returns one e-value, the first term of the
process; the process itself is the product that :func:`combine_independent`
forms as batches arrive.

Construction
------------
The e-value is the split likelihood ratio of Wasserman, Ramdas
and Balakrishnan (2020), often called universal inference: fit the alternative
on a training fold, evaluate its likelihood on the held-out fold, and divide by
the null likelihood maximised on that same held-out fold. Because the numerator
parameters never saw the held-out rows, the ratio has expectation at most one
under the null with no regularity conditions on the model. That is exactly the
fold structure the clonal share already uses, so the evidence costs one extra
pass over data that has already been split.

Group probabilities are shrunk toward the training grand mean by a one-way
empirical-Bayes factor estimated on the training rows. Shrinkage is not a
convenience here: an unshrunk group of one training isolate gives a probability of zero or
one, the log likelihood diverges, and the e-value becomes infinite for a reason
that has nothing to do with lineage.

Combining
---------
Within one cohort the folds share data, so their e-values are
averaged, which is valid because any convex combination of e-values is an
e-value. Across independent batches, for instance successive years of
collection, they are multiplied, which forms a test martingale and is what
makes the procedure anytime-valid. Across a panel of antimicrobials the e-BH
procedure of Wang and Ramdas (2022) controls the false discovery rate under
arbitrary dependence between agents, which matters because cross-resistance
makes a macrolide block behave as one trait.

One look and many
-----------------
:func:`e_process` is the e-value of one look: recomputing it on the accumulated
cohort after every batch gives a valid e-value at each look, and e-BH on those
values controls the false discovery rate at each look, but the sequence of
recomputed values is not a test martingale, so no guarantee attaches to the
programme of looks as a whole. The anytime-valid construction for a programme
is :func:`sequential_e_process`: at each new batch the alternative is fitted on
every batch before it, the batch is scored against the null maximised on the
batch alone, and the factors are multiplied. Each factor has conditional
expectation at most one given the past, so the running product is a test
supermartingale and may be inspected after any batch and stopped at any time.
It is the fixed-look split likelihood ratio with the training fold replaced by
the past, which is why its power grows with the programme rather than with the
batch.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any, Dict, Iterable, List, Sequence

import numpy as np

from .attribution import _codes, _folds, _is_untyped, _rng

__all__ = ["EResult", "SequentialEResult", "e_process", "sequential_e_process",
           "e_bh", "e_bh_log", "combine_independent", "combine_within_cohort"]

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
    n_dropped_untyped : isolates set aside because their lineage cell was
        empty; an untyped isolate is not a lineage, and counting it as one
        would turn the pattern of typing into evidence about the trait.
    n_dropped_non_finite : isolates set aside because the trait value was not
        finite. A missing well is normal in a susceptibility panel; the
        isolate leaves this trait's evidence rather than becoming a NaN.
        Counted on the cohort as supplied, as ``n_dropped_untyped`` is, so an
        isolate that is both untyped and unread is counted in both.
    """

    e_value: float
    log_e: float
    reject_05: bool
    reject_01: bool
    n: int
    n_groups: int
    prevalence: float
    n_splits: int
    n_dropped_untyped: int = 0
    n_dropped_non_finite: int = 0

    def as_dict(self) -> Dict[str, Any]:
        return asdict(self)


def _shrunk_group_p(y_tr: np.ndarray, code_tr: np.ndarray, G: int):
    """Empirical-Bayes group probabilities from the training rows.

    The shrinkage factor is ``tau2 / (tau2 + sigma2 / n_g)``, with both
    components from a one-way moment estimator on the training rows only. A group absent from training
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
    labels = np.asarray(lineage, dtype=object).ravel()
    if labels.shape[0] != y.size:
        raise ValueError(f"lineage has {labels.shape[0]} entries for {y.size} traits")
    # A missing well is normal in a susceptibility panel; the isolate leaves
    # this trait's evidence rather than turning it into a NaN that the
    # step-up would read as no evidence at all. An untyped isolate leaves it
    # too, as it does in every other estimator of the package.
    typed = np.array([not _is_untyped(v) for v in labels], dtype=bool)
    finite = np.isfinite(y)
    keep = finite & typed
    n_dropped_untyped = int((~typed).sum())
    n_dropped_non_finite = int((~finite).sum())
    y, code = y[keep], _codes(labels[keep])
    n = y.size
    rng = _rng(seed)
    G = int(code.max()) + 1 if n else 0
    prevalence = float(y.mean()) if n else float("nan")

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
                       prevalence, 0, n_dropped_untyped, n_dropped_non_finite)
    logs = np.asarray(per_split, dtype=float)
    # Average on the natural scale, in a way that does not overflow.
    m = float(logs.max())
    e = float(np.exp(m) * np.exp(logs - m).mean())
    return EResult(e_value=e, log_e=float(np.log(e)) if e > 0 else float("-inf"),
                   reject_05=bool(e >= REJECT_AT[0.05]),
                   reject_01=bool(e >= REJECT_AT[0.01]),
                   n=n, n_groups=G, prevalence=prevalence,
                   n_splits=len(per_split),
                   n_dropped_untyped=n_dropped_untyped,
                   n_dropped_non_finite=n_dropped_non_finite)


@dataclass(frozen=True)
class SequentialEResult:
    """The running evidence of a programme of looks.

    log_e : the running log e-value after each batch, the first batch
        contributing nothing because nothing preceded it to train on.
    e_value : the same on the natural scale, infinite where the log overflows.
    reject_05, reject_01 : whether the running value has ever cleared
        ``1/alpha``; by Ville's inequality that happens under the null with
        probability at most ``alpha`` however many batches are inspected.
    n_batches, n_per_batch : the programme the evidence came from, after
        the untyped isolates of each batch were set aside (their count is
        ``n_dropped_untyped``).
    n_dropped_non_finite : isolates whose trait value was not finite, summed
        over the batches as supplied. They enter neither the training nor the
        held-out rows of any look. Counted on the batches as supplied, as
        ``n_dropped_untyped`` is, so an isolate that is both untyped and
        unread is counted in both.
    """

    log_e: List[float]
    e_value: List[float]
    reject_05: bool
    reject_01: bool
    n_batches: int
    n_per_batch: List[int]
    n_dropped_untyped: int = 0
    n_dropped_non_finite: int = 0

    def as_dict(self) -> Dict[str, Any]:
        return asdict(self)


def sequential_e_process(batches, lineage_batches, *,
                         min_training: int = 2) -> SequentialEResult:
    """Anytime-valid evidence over a programme of batches.

    ``batches`` is a sequence of trait vectors, one per look in arrival order,
    and ``lineage_batches`` the matching lineage labels. For batch ``t`` the
    lineage probabilities are fitted, with the same empirical-Bayes shrinkage
    as :func:`e_process`, on batches ``1`` to ``t - 1`` pooled, and the batch
    is scored against the null probability maximised on the batch itself:

        E_t = prod_i q_{t-1}(y_ti) / sup_p prod_i f_p(y_ti).

    Given the past, the numerator is a fixed probability law on the batch and
    the denominator dominates every null law, so ``E[E_t | past] <= 1`` for
    every null ``p``, and the product over batches is a test supermartingale
    with initial value one. A batch with fewer than ``min_training`` isolates
    before it, or whose training rows hold a single lineage, contributes a
    factor of one: no evidence is claimed where none could be earned.

    The batches must be disjoint sets of isolates in the order they were
    collected. Re-using an isolate in two batches, or reordering batches after
    seeing them, breaks the conditioning the guarantee rests on; the function
    cannot detect either and does not try to.
    """
    ys = [np.asarray(b, dtype=float).ravel() for b in batches]
    ls = [np.asarray(b, dtype=object).ravel() for b in lineage_batches]
    if len(ys) != len(ls):
        raise ValueError(f"{len(ys)} trait batches for {len(ls)} lineage batches")
    for t, (y, lab) in enumerate(zip(ys, ls)):
        if y.size != lab.size:
            raise ValueError(f"batch {t}: {y.size} traits for {lab.size} labels")
    # An untyped isolate is set aside within its batch, as in every other
    # estimator of the package; it is not a lineage of its own.
    typed = [np.array([not _is_untyped(v) for v in lab], dtype=bool) for lab in ls]
    n_dropped_untyped = int(sum(int((~t).sum()) for t in typed))
    n_dropped_non_finite = int(sum(int((~np.isfinite(y)).sum()) for y in ys))
    ys = [y[t] for y, t in zip(ys, typed)]
    ls = [lab[t] for lab, t in zip(ls, typed)]
    code_all = _codes(np.concatenate(ls)) if ys else np.zeros(0, dtype=np.int64)
    G = int(code_all.max()) + 1 if code_all.size else 0
    sizes = [int(y.size) for y in ys]
    bounds = np.cumsum([0] + sizes)
    y_all = np.concatenate(ys) if ys else np.zeros(0)
    running = 0.0
    logs: List[float] = []
    for t in range(len(ys)):
        lo, hi = int(bounds[t]), int(bounds[t + 1])
        y_tr, c_tr = y_all[:lo], code_all[:lo]
        y_te, c_te = y_all[lo:hi], code_all[lo:hi]
        keep = np.isfinite(y_te)
        y_te, c_te = y_te[keep], c_te[keep]
        ok_tr = np.isfinite(y_tr)
        y_tr, c_tr = y_tr[ok_tr], c_tr[ok_tr]
        if (y_tr.size >= int(min_training) and y_te.size
                and np.unique(c_tr).size >= 2):
            p_alt, _ = _shrunk_group_p(y_tr, c_tr, G)
            p_null = float(np.clip(y_te.mean(), 1e-6, 1 - 1e-6))
            running += (_log_bernoulli(y_te, p_alt[c_te])
                        - _log_bernoulli(y_te, np.full(y_te.size, p_null)))
        logs.append(float(running))
    e_vals = [float(np.exp(v)) if v < 700 else float("inf") for v in logs]
    peak = max(logs) if logs else float("-inf")
    return SequentialEResult(
        log_e=logs, e_value=e_vals,
        reject_05=bool(peak >= np.log(REJECT_AT[0.05])),
        reject_01=bool(peak >= np.log(REJECT_AT[0.01])),
        n_batches=len(ys), n_per_batch=sizes,
        n_dropped_untyped=n_dropped_untyped,
        n_dropped_non_finite=n_dropped_non_finite)


def combine_independent(values: Iterable[float]) -> float:
    """Multiply e-values from independent batches, for example successive
    years of collection. The running product is a test martingale, which is
    what makes it legitimate to inspect the total after every batch and stop
    at any point without spending a significance budget.

    An e-value of ``+inf`` is evidence beyond the floating-point range and is
    admissible here as it is in :func:`e_bh`: the product saturates at
    ``+inf`` and stays there, which is the right reading, since no later batch
    can take evidence away. A NaN is a missing e-value and a negative value is
    not an e-value; both raise."""
    out = 1.0
    for v in values:
        if np.isnan(v) or v < 0:
            raise ValueError(f"not an e-value: {v}")
        out *= float(v)
    if np.isnan(out):
        raise ValueError("an infinite e-value multiplied by a zero one; their "
                         "product is not defined and must not be read as one")
    return out


def combine_within_cohort(values: Iterable[float]) -> float:
    """Average e-values computed on overlapping splits of one cohort. A convex
    combination of e-values is an e-value; a product of dependent ones is
    not."""
    v = np.asarray(list(values), dtype=float)
    if v.size == 0 or np.any(v < 0):
        raise ValueError("e-values must be non-negative and non-empty")
    return float(v.mean())


def e_bh(evalues: Sequence[float], alpha: float = 0.05) -> Dict[str, Any]:
    """The e-BH procedure of Wang and Ramdas (2022).

    Sort the e-values in decreasing order and reject the largest ``k`` for
    which the ``k``-th largest is at least ``m / (alpha * k)``. This controls
    the false discovery rate at ``alpha`` under *arbitrary* dependence between
    the hypotheses, with no correction. That property is what makes it the
    right instrument for an antimicrobial panel, where cross-resistance means
    the agents are neither independent nor reliably positively dependent, so
    the usual justification for Benjamini-Hochberg does not apply.

    An e-value of ``+inf`` is evidence beyond the floating-point range, as a
    running product over many intakes can be, and ranks above every finite
    value; a NaN is a missing e-value and counts as no evidence while keeping
    its place in ``m``, so that dropping a hypothesis cannot loosen the
    threshold for the others. A negative value is not an e-value and raises.

    Returns the indices rejected and the effective threshold.
    """
    e = np.asarray(list(evalues), dtype=float)
    if np.any(e < 0.0):
        raise ValueError("an e-value is non-negative; a negative value was passed")
    with np.errstate(divide="ignore"):
        return e_bh_log(np.log(np.where(np.isnan(e), 0.0, e)), alpha)


def e_bh_log(log_evalues: Iterable[float] | np.ndarray,
             alpha: float = 0.05) -> Dict[str, Any]:
    """:func:`e_bh` on the log scale, where a product over many intakes
    stays finite. ``-inf`` is an e-value of zero, ``+inf`` evidence beyond
    the floating-point range, NaN a missing e-value counted as zero."""
    le = np.asarray(list(log_evalues), dtype=float)
    m = le.size
    if m == 0:
        return {"rejected": [], "n_rejected": 0, "threshold": float("inf"),
                "alpha": float(alpha), "m": 0}
    le = np.where(np.isnan(le), -np.inf, le)
    order = np.argsort(-le, kind="stable")
    sorted_le = le[order]
    k_star = 0
    for k in range(1, m + 1):
        if sorted_le[k - 1] >= np.log(m / (alpha * k)):
            k_star = k
    rejected = sorted(int(i) for i in order[:k_star])
    thr = m / (alpha * k_star) if k_star else float("inf")
    return {"rejected": rejected, "n_rejected": k_star, "threshold": float(thr),
            "alpha": float(alpha), "m": int(m)}
