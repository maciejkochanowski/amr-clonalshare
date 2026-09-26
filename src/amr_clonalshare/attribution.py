"""Predictive association of recorded phenotypes with supplied lineage labels.

The score compares held-out squared prediction error from lineage-specific
rates with a marginal predictor. Random isolate-level folds target another
isolate from the observed lineage distribution, not an unseen lineage.
Permutation calibration adjusts finite-sample prediction penalties under the
stated exchangeability null; it does not remove sampling bias or confounding.

Intervals and support gates are empirically assessed under documented designs.
A cluster bootstrap requires a sampling model and is not assumption-free.
Lineage association does not identify transmission, horizontal transfer,
heritability, clinical resistance or a causal genetic mechanism.
"""
from __future__ import annotations

import hashlib
import typing
from numbers import Integral
from dataclasses import dataclass

import numpy as np

#: Lineages with at least two isolates the share needs. A singleton lineage
#: cannot be predicted out of sample -- its one isolate is in the training
#: fold or the held-out fold, never both -- so the share is scored on the
#: isolates of repeated lineages and singletons are set aside and counted.
#: Cell D of ``benchmarks/attribution_calibration.py`` checks its level as the
#: support falls, under the condition the record states: that singleton
#: lineages are drawn from the same population of lineages as the repeated
#: ones.
MIN_REPEATED_LINEAGES = 2

#: Below this many lineages the lineage bootstrap is not the interval to
#: read: a percentile bootstrap over five lineages has too few distinct
#: resamples, and on the estimator grid (``validation_grid.json``, entry
#: ``few_lineages``) it lost its level with five lineages and a non-zero share
#: while from ten lineages up it held. The estimate is returned either way;
#: the record flags the case and the report names it.
FEW_LINEAGES = 10

#: Code given to an isolate with no lineage assignment while labels are
#: coded; the estimators set those isolates aside before fitting and count
#: them in ``n_dropped_untyped``, because an untyped isolate is not a member
#: of any lineage.
_MISSING = "__missing__"

__all__ = [
    "MIN_REPEATED_LINEAGES",
    "FEW_LINEAGES",
    "ShareResult",
    "clonal_share",
    "layer_clonal_share",
]


# --------------------------------------------------------------------- types
@dataclass(frozen=True)
class ShareResult:
    """Out-of-sample variance share of a categorical predictor.

    Attributes
    ----------
    kappa : the raw out-of-sample skill of the group-mean predictor against the
        marginal-prevalence predictor (see :func:`_skill`). It is biased low by
        the cost of estimating group means from a training fold.
    kappa_adj : the bias-corrected estimate, and the one to report.
        Cross-validation charges a penalty ``c`` per group that is paid whether
        or not the grouping carries signal, so ``kappa_hat = kappa - (1-kappa)c``
        while the permuted-label run measures ``kappa_null = -c`` on the same
        design and the same folds. Therefore

            (kappa_hat - kappa_null) / (1 - kappa_null) = kappa

        exactly under that model, and the correction needs no extra assumption
        about the number of levels or the cohort size, because both are already
        in ``kappa_null``.
    ci_low, ci_high : percentile bootstrap interval that draws lineages
        whole, on the corrected scale. It is an interval for the
        lineage-membership share of this collection, B_w / (B_w + sigma^2)
        with B_w the isolate-weighted variance of the lineage means: the
        predictor scores those means, and the resample never draws a lineage
        the cohort does not hold. That quantity is not the design-corrected
        component ratio ``realised_share`` estimates, whose between-lineage
        component is divided by 1 - sum w_g^2. The two part company where lineages are
        few. ``validation_grid.json`` records how often this interval covers
        its own target and how often it covers the superpopulation target,
        the share a fresh draw of lineages would show, which it misses more
        often where lineages are few. ``interval_target`` says so in the
        record, and ``superpopulation_low``/``superpopulation_high`` carry
        the interval for that second question.
    null_mean : mean kappa when the predictor's labels are permuted. Negative
        rather than zero, and its magnitude is the penalty ``c`` above.
    p_value : permutation p-value, Phipson-Smyth corrected. It cannot fall
        below ``p_floor = 1 / (n_perm + 1)``; a value at the floor says no
        permutation reached the observed share, and the permuted-label share
        with its interval, not the p-value, is the control to read.
    n, n_groups, prevalence : cohort description carried with the estimate.
    support : share of the typed isolates sitting in a lineage with at least
        two members. The share is scored on those isolates only; a singleton
        lineage can never be predicted out of sample, and scored with the
        others it would pull the share towards zero. ``n`` counts the retained
        isolates, ``n_scored`` those scored and ``n_singletons_set_aside`` the
        ones set aside. A low support
        does not bias the share, provided singleton lineages are drawn from the
        same population of lineages as the repeated ones; it does make the
        share describe fewer isolates, which the record says.
    estimable : ``kappa_adj`` finite and at least
        :data:`MIN_REPEATED_LINEAGES` repeated lineages; a share that did not
        come back is not an estimate.
    cv_sd : spread of kappa across repeated fold assignments; a diagnostic of
        the estimator, not of the cohort, and never a confidence interval.
    n_dropped_non_finite : isolates removed before the estimate because a trait
        column was not finite. A missing well is normal in a susceptibility
        panel, and every comparison against a NaN is false, so leaving those
        rows in made the permuted runs never exceed the observed skill and
        drove the p-value to its floor. The count is carried because the
        estimate then describes a subset of the cohort.
    n_dropped_untyped : isolates removed before the estimate because they
        carry no lineage label. An untyped isolate is not a lineage, and a
        level made of every isolate whose typing failed measures the typing
        process rather than lineage membership; ``missing_share`` reports how
        much of the cohort was set aside, and the estimate describes the typed
        isolates.
    """

    kappa: float
    kappa_adj: float
    ci_low: float
    ci_high: float
    null_mean: float
    p_value: float
    n: int
    n_groups: int
    prevalence: float
    support: float
    missing_share: float
    estimable: bool
    cv_sd: float
    #: Defaulted so that constructions written before the drop existed still
    #: build a valid record.
    n_dropped_non_finite: int = 0
    n_dropped_untyped: int = 0
    #: isolates of singleton lineages, typed and read but not scored
    n_singletons_set_aside: int = 0
    #: isolates scored: those of lineages with at least two members
    n_scored: int = 0
    #: The quantity ``ci_low``/``ci_high`` is calibrated for. Always the
    #: lineage-membership share of this collection, ``B_w / (B_w + sigma^2)``
    #: with ``B_w`` the isolate-weighted variance of the lineage means, which
    #: is what the out-of-sample estimate converges to and is not the
    #: design-corrected component ratio ``realised_share`` reports. It is
    #: carried in the record so that a reader of ``clonal_share_result.json``
    #: does not have to know the estimator to know which question the
    #: interval answers.
    interval_target: str = "lineage_membership"
    #: The interval for the superpopulation share, the share a fresh draw of
    #: lineages from the species would show. Built by :func:`_superpopulation_
    #: interval` from a within-lineage bootstrap and a chi-square layer on
    #: ``n_groups - 1`` degrees of freedom, widened to the envelope of that
    #: and the lineage bootstrap; NaN where the estimate is. It is stated on
    #: the design-corrected component scale ``S_a^2 / (S_a^2 + sigma^2)``,
    #: the scale of a population intraclass correlation, not on the
    #: lineage-membership scale of ``kappa_adj``; ``superpopulation_share`` is
    #: the estimate on that scale.
    superpopulation_low: float = float("nan")
    superpopulation_high: float = float("nan")
    superpopulation_share: float = float("nan")
    #: The smallest p-value ``n_perm`` permutations can return,
    #: ``1 / (n_perm + 1)``. A ``p_value`` equal to it says every permutation
    #: fell below the observed share, not that the p-value is that number.
    p_floor: float = float("nan")
    #: For a single binary trait, the same share and its two intervals on the
    #: model-equivalent scale of a Gaussian probit threshold model, using
    #: plug-in prevalence (:mod:`amr_clonalshare.latent`). Collection and
    #: population endpoints retain their different targets; coverage is not
    #: guaranteed by the transformation. NaN for a continuous trait or a
    #: block of several traits. A share at or below zero maps to zero.
    latent_share: float = float("nan")
    latent_low: float = float("nan")
    latent_high: float = float("nan")
    latent_species_low: float = float("nan")
    latent_species_high: float = float("nan")
    #: For a single trait, the largest share of the between-lineage sum of
    #: squares that one lineage carries. The chi-square layer of the species
    #: interval describes lineage effects drawn from one law; a collection in
    #: which one lineage carries most of the between-lineage variation is not
    #: that, and on the estimator grid a two-point carrier law can take the
    #: species interval below its level at few lineages while the interval for
    #: the lineages in hand holds. NaN for a block of several traits.
    dominant_lineage_share: float = float("nan")
    #: True below :data:`FEW_LINEAGES`, where the lineage bootstrap loses its
    #: level; with so few lineages the species interval is the one to read,
    #: and the report says so.
    few_lineages: bool = False
    #: Lineages with at least two analysed isolates, the only ones a held-out
    #: prediction can learn from.
    n_groups_repeated: int = 0
    #: Effective number of lineages among the analysed isolates, the inverse
    #: of the sum of squared lineage shares (Hill number of order two).
    effective_groups: float = float("nan")
    #: Share of the analysed isolates in the largest lineage.
    largest_group_share: float = float("nan")
    #: Positive outcomes among isolates in repeated lineages, for a single
    #: binary trait; NaN for a block of traits.
    n_positive_repeated: float = float("nan")
    #: Explicit diagnostic when the null correction cannot be defined.
    reason: str = ""

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

    Reported rather than gated. Untyped isolates are set aside before the
    estimate, and if typing failed informatively -- as it does on the shipped
    *S. suis* cohort at sequence-type resolution, where the untyped isolates
    carry about two more non-wild-type results out of thirteen -- the estimate
    describes a subset selected by the typing process. The number is emitted
    so a reader can see how much of the cohort the estimate rests on, and so
    that the same cohort typed two ways can be compared.
    """
    arr = np.asarray(labels, dtype=object).ravel()
    return float(np.mean([_is_untyped(v) for v in arr])) if arr.size else float("nan")


def _folds(n: int, k: int, rng: np.random.Generator) -> np.ndarray:
    order = rng.permutation(n)
    fold = np.empty(n, dtype=np.int64)
    fold[order] = np.arange(n) % k
    return fold


def _skill(X: np.ndarray, code: np.ndarray, fold: np.ndarray) -> float:
    """Out-of-sample variance share of ``code`` for columns of ``X``.

    Held-out rows are predicted by the mean of their group among the training
    rows; a group absent from training falls back to the training grand mean,
    which is what an honest predictor would do and is the mechanism by which
    extra levels stop paying for themselves. The score is

        1 - SSE(group means) / SSE(marginal mean)

    on the pooled sum of squares over all columns. For a single binary column
    this is exactly the Brier skill score against the prevalence baseline, and
    for several columns it is the multivariate out-of-sample R-squared. Both
    denominators are evaluated out of sample as well, so the score is not
    flattered by an in-sample baseline.
    """
    X = np.atleast_2d(np.asarray(X, dtype=float))
    if X.shape[0] != fold.shape[0]:
        X = X.T
    p = X.shape[1]
    pred = np.empty_like(X)
    base = np.empty_like(X)
    for f in np.unique(fold):
        tr, te = fold != f, fold == f
        if not tr.any() or not te.any():
            continue
        c = code
        G = int(c.max()) + 1
        ctr = c[tr]
        counts = np.bincount(ctr, minlength=G).astype(float)
        sums = np.empty((G, p))
        for j in range(p):
            sums[:, j] = np.bincount(ctr, weights=X[tr, j], minlength=G)
        grand = X[tr].mean(axis=0)
        means = np.where(counts[:, None] > 0, sums / np.maximum(counts[:, None], 1),
                         grand[None, :])
        pred[te] = means[c[te]]
        base[te] = grand
    sse = float(((X - pred) ** 2).sum())
    sst = float(((X - base) ** 2).sum())
    return float("nan") if sst <= 0 else 1.0 - sse / sst


def _repeated_skill(X, code, *, folds, repeats, rng):
    n = np.asarray(code).shape[0]
    vals = np.array([_skill(X, code, _folds(n, folds, rng)) for _ in range(repeats)])
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


def _phipson_smyth(exceed: int, n_perm: int) -> float:
    """Permutation p-value that cannot be zero (Phipson & Smyth 2010)."""
    return (exceed + 1.0) / (n_perm + 1.0)


def _group_index(code: np.ndarray):
    """Sort isolates by lineage once, so a bootstrap replicate is O(n).

    Rebuilding the membership lists inside the resampling loop costs O(G n) per
    replicate, which on 768 sequence types, 1500 isolates, 400 replicates and
    seventeen agents is billions of operations and turned a one-minute run into
    an unfinished one. The sorted order and the group boundaries do not change
    between replicates, so they are computed once.
    """
    order = np.argsort(code, kind="stable")
    counts = np.bincount(code)
    bounds = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)
    return order, bounds


def _cluster_resample_fast(order: np.ndarray, bounds: np.ndarray,
                           rng: np.random.Generator):
    """Draw lineages with replacement from a prepared index.

    Lineages are taken whole because the isolates of one lineage are not
    exchangeable with the isolates of another: a resample that split them
    would put rows of the same lineage on both sides of a fold and score a
    group mean on the isolates that formed it."""
    G = len(bounds) - 1
    picks = rng.integers(0, G, G)
    rows, new = [], []
    for j, gi in enumerate(picks):
        idx = order[bounds[gi]:bounds[gi + 1]]
        if idx.size:
            rows.append(idx)
            new.append(np.full(idx.size, j, dtype=np.int64))
    if not rows:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    return np.concatenate(rows), np.concatenate(new)


def _debias(kappa: float, null_mean: float) -> float:
    """Remove the cross-validation penalty measured on the permuted design.

    ``kappa_hat = kappa - (1 - kappa) c`` and ``null_mean = -c``, so the map
    below returns ``kappa`` exactly under that model. The denominator is
    ``1 + c >= 1`` under that penalty model. A finite sampled permutation
    control can instead reach one on a small design. At that boundary the
    correction is undefined; values above one also violate the positive
    rescaling denominator and must not produce a reportable estimate.
    """
    if not np.isfinite(kappa) or not np.isfinite(null_mean):
        return float("nan")
    if null_mean >= 1.0:
        return float("nan")
    return float((kappa - null_mean) / (1.0 - null_mean))


# ------------------------------------------------------------- public: agent
def _stratified_resample(code: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Isolate indices drawn with replacement inside each lineage, sizes kept.

    The within-lineage resample holds the set of lineages fixed and moves only
    the isolates, so its spread is the uncertainty of the share those lineages
    carry and nothing else; the between-lineage sampling is added separately
    by :func:`_superpopulation_interval`, so neither source is counted twice.
    """
    order = np.argsort(code, kind="stable")
    sizes = np.bincount(code)
    out = np.empty(code.size, dtype=np.int64)
    start = 0
    for g, size in enumerate(sizes):
        block = order[start:start + size]
        out[start:start + size] = block[rng.integers(0, size, size)]
        start += int(size)
    return out


def _to_component_scale(x, design: float):
    """Move a lineage-membership share onto the component-ratio scale.

    The out-of-sample share estimates ``B_w / (B_w + sigma^2)`` with ``B_w``
    the isolate-weighted variance of the lineage means, while the chi-square
    law of the species interval is written for the design-corrected ratio
    ``S_a^2 / (S_a^2 + sigma^2)`` with ``S_a^2 = B_w / (1 - sum w_g^2)``. The
    two differ by that factor on the odds scale and by nothing else, so the
    odds are divided by ``design = 1 - sum w_g^2`` before the law is inverted.
    Values at or above one, and values whose transformed odds fall at or below
    minus one, are returned as they are.
    """
    x = np.asarray(x, dtype=float)
    out = x.copy()
    ok = np.isfinite(x) & (x < 1.0)
    odds = np.where(ok, x / np.where(ok, 1.0 - x, 1.0), 0.0) / design
    good = ok & (odds > -1.0)
    out[good] = odds[good] / (1.0 + odds[good])
    return out if out.ndim else float(out)


def _superpopulation_interval(kappa_hat: float, within_draws: np.ndarray,
                              n_groups: int, rng: np.random.Generator,
                              n_sim: int = 4000, alpha: float = 0.05):
    """Interval for the share a fresh draw of lineages would show.

    The lineage effects in hand are one draw of ``n_groups`` from the species,
    so the realised between-lineage dispersion is the superpopulation value
    times a chi-square on ``n_groups - 1`` degrees of freedom over its degrees
    of freedom, and the realised share for a candidate superpopulation share
    ``rho`` is ``rho Q / (rho Q + 1 - rho)``. That law is written on the scale
    of the design-corrected component ratio, so the caller hands in the
    estimate and its within-lineage draws on that scale (see
    :func:`_to_component_scale`). The estimate adds the
    within-lineage noise the stratified bootstrap measured. The interval is
    the set of ``rho`` for which the observed estimate is not in either
    ``alpha / 2`` tail of that predictive law, found by bisection. With ten
    lineages the chi-square layer is what a lineage bootstrap misses; with a
    hundred it is negligible and the two intervals nearly coincide. The
    caller widens the result to the envelope of this interval and the
    lineage bootstrap, so that a collection dominated by one lineage, whose
    effects do not follow the chi-square law, keeps the wider statement.

    The within-lineage noise is resampled from the draws measured at the
    observed share, and the same centred draws are added at every candidate
    ``rho``: the law of that noise is held independent of the share, which it
    is not exactly, since a share near either end has less room for it. That
    is an approximation, and its level is what the grid measures rather than
    something the construction guarantees.
    """
    draws = within_draws[np.isfinite(within_draws)]
    if not np.isfinite(kappa_hat) or draws.size < 20 or n_groups < 2:
        return float("nan"), float("nan")
    q = rng.chisquare(n_groups - 1, n_sim) / (n_groups - 1)
    e = rng.choice(draws - draws.mean(), n_sim)

    def realised(rho):
        return rho * q / (rho * q + 1.0 - rho)

    def below(rho):  # P(estimate <= observed | rho)
        return float(((realised(rho) + e) <= kappa_hat).mean())

    def above(rho):  # P(estimate >= observed | rho)
        return float(((realised(rho) + e) >= kappa_hat).mean())

    lo, hi = 0.0, 1.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if below(mid) >= alpha / 2:
            lo = mid
        else:
            hi = mid
    upper = lo
    if above(0.0) >= alpha / 2:
        lower = 0.0
    else:
        lo, hi = 0.0, 1.0
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if above(mid) >= alpha / 2:
                hi = mid
            else:
                lo = mid
        lower = hi
    # below() falls and above() rises with rho, and the two sum to at least
    # one, so the set is an interval and never empty; min() only absorbs the
    # last bisection step.
    return float(min(lower, upper)), float(upper)


def _reason(null_mean: float, low: float, high: float, n_boot: int,
            repeated_groups: int) -> str:
    """Why a number is missing, in the words the report prints.

    An estimate without an interval is reported as computed, so the record has
    to say what happened to the interval; the point-only route of ``n_boot=0``
    is not a failure and says nothing.
    """
    if null_mean >= 1.0:
        return ("permuted null mean is at least one; the null correction denominator "
                "is nonpositive, so the corrected share is not estimable")
    if n_boot and not (np.isfinite(low) and np.isfinite(high)):
        return (f"the lineage bootstrap produced no usable draw: {repeated_groups} "
                "lineage(s) hold repeated isolates, so resampling whole lineages "
                "cannot form a percentile; the estimate stands without an interval")
    return ""


def layer_clonal_share(X, lineage, *, folds: int = 5, repeats: int = 20,
                       n_boot: int = 400, n_perm: int = 200,
                       seed=0, null_repeats: int = 5) -> ShareResult:
    """Out-of-sample share of a trait block's variance explained by lineage.

    Parameters
    ----------
    X : (n,) or (n, p) array of binary trait calls, 1 = non-wild-type.
    lineage : (n,) label vector; MLST sequence type, a BAPS cluster, a clonal
        complex, anything categorical. It is used as a label and never as a
        tree, which is the point: a laboratory that types by MLST has one.
    folds, repeats : cross-validation design. ``repeats`` averages over fold
        draws so the estimate does not depend on one fold split of the cohort.
    n_boot : percentile bootstrap replicates drawing whole lineages; also
        used for the separate resampling of isolates within lineages.
    n_perm : permutations of the lineage label for the null and the p-value.
    null_repeats : fold draws averaged on each side of the permutation test.
        The reported estimate is a mean over ``repeats`` fold draws; a
        permuted statistic computed from one draw carries fold noise that
        the estimate has averaged away, and a p-value that compares the two
        is conservative. The test therefore compares the mean of the first
        ``null_repeats`` observed draws with the mean of ``null_repeats``
        draws under each permutation, the same function of the labels on
        both sides, which is exact under exchangeability at any value. More
        draws buy power at a cost of ``n_perm`` skill evaluations each;
        ``null_repeats`` larger than ``repeats`` is read as ``repeats``.

    Returns
    -------
    ShareResult
    """
    for name, value, minimum in (("folds", folds, 2), ("repeats", repeats, 1),
                                  ("n_boot", n_boot, 0), ("n_perm", n_perm, 0),
                                  ("null_repeats", null_repeats, 1)):
        if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
            # One fold would hold nothing out, so every score would be taken
            # on an empty held-out set; the minimum of two folds is enforced
            # here with the other budgets.
            raise ValueError(f"{name} must be an integer >= {minimum}; at least two folds are required" if name == "folds" else f"{name} must be an integer >= {minimum}")
    rng = _rng(seed)
    X = np.atleast_2d(np.asarray(X, dtype=float))
    lineage = np.asarray(lineage, dtype=object).ravel()
    if X.shape[0] != lineage.size:
        if X.shape[1] == lineage.size:
            X = X.T
        else:
            # Transposing on any mismatch turned 60 traits against 59 labels
            # into a cohort of n = 1 and reported a number for it.
            raise ValueError(f"X has {X.shape[0]} rows and lineage has "
                             f"{lineage.size} entries; they must match")
    # The share of the cohort with no label is read before any row is set
    # aside, so it describes the cohort as supplied.
    missing_share = _missing_share(lineage) if lineage.size else float("nan")
    finite = np.isfinite(X).all(axis=1)
    n_dropped_non_finite = int((~finite).sum())
    # Both counts describe the cohort as supplied, as the other modules count
    # them; an isolate that is unread and untyped is counted in both.
    typed = np.array([not _is_untyped(v) for v in lineage], dtype=bool)
    n_dropped_untyped = int((~typed).sum())
    keep = finite & typed
    if not keep.all():
        # Every comparison against a NaN is false, so `null >= kappa` counted
        # no exceedances and the permutation p-value collapsed to its floor
        # 1 / (n_perm + 1) -- the most significant value obtainable -- beside a
        # NaN kappa_adj and estimable true. A missing well is normal in a
        # susceptibility panel, so the row goes and is counted.
        X, lineage = X[keep], lineage[keep]
    # Singleton lineages are set aside before anything is scored (see
    # ``support`` above); their count and the support describe the typed cohort.
    typed_sizes = np.bincount(_codes(lineage)) if lineage.size else np.zeros(0, dtype=int)
    typed_code = _codes(lineage) if lineage.size else np.zeros(0, dtype=int)
    repeated = typed_sizes[typed_code] >= 2 if lineage.size else np.zeros(0, dtype=bool)
    support = float(repeated.mean()) if lineage.size else float("nan")
    n_singletons_set_aside = int((~repeated).sum())
    # ``n`` and ``n_groups`` in the record describe the retained cohort; the
    # scored isolates are ``n_scored``.
    n_retained = int(lineage.size)
    n_groups_all = int(typed_sizes.size)
    if not repeated.all():
        X, lineage = X[repeated], lineage[repeated]
    code = _codes(lineage)
    n = X.shape[0]
    if n < 2:
        # Fewer than two isolates leave nothing to hold out and nothing to
        # compare, so there is no share to estimate.
        nan = float("nan")
        return ShareResult(kappa=nan, kappa_adj=nan, ci_low=nan, ci_high=nan,
                           null_mean=nan, p_value=nan, n=n_retained,
                           n_groups=n_groups_all,
                           prevalence=float(X.mean()) if X.size else nan,
                           support=support,
                           missing_share=missing_share,
                           estimable=False, cv_sd=nan,
                           n_dropped_non_finite=n_dropped_non_finite,
                           n_dropped_untyped=n_dropped_untyped,
                           n_singletons_set_aside=n_singletons_set_aside, n_scored=int(n))
    sizes = np.bincount(code)
    shares = sizes[sizes > 0]/n
    structure: typing.Dict[str, typing.Any] = dict(
        n_groups_repeated=int((sizes >= 2).sum()),
        effective_groups=float(1.0/np.sum(np.sort(shares)**2)),
        largest_group_share=float(shares.max()),
        n_positive_repeated=(float(X[:, 0].sum())
                             if X.shape[1] == 1 and np.isin(X, (0.0, 1.0)).all() else float("nan")),
        n_singletons_set_aside=n_singletons_set_aside, n_scored=int(n))

    if sizes.size < 2:
        # One lineage group leaves no contrast between groups, so the share is
        # not defined on this cohort. Without this the run returns kappa 0, a
        # zero-width interval and p = 1, which reads as "lineage explains
        # nothing" when in fact nothing was compared -- the same false reading
        # the trait-variance guard below exists to prevent. A public release
        # that assigns a whole species to one cluster produces exactly this.
        nan = float("nan")
        return ShareResult(kappa=nan, kappa_adj=nan, ci_low=nan, ci_high=nan,
                           null_mean=nan, p_value=nan, n=n_retained,
                           n_groups=n_groups_all,
                           prevalence=float(X.mean()), support=support,
                           missing_share=missing_share,
                           estimable=False, cv_sd=nan,
                           n_dropped_non_finite=n_dropped_non_finite,
                           n_dropped_untyped=n_dropped_untyped, **structure)

    if float(np.ptp(X, axis=0).max()) <= 0:
        # A trait present in every isolate, or in none, has no variance to
        # attribute. Returning NaN says so; returning 0 would read as "not
        # clonal", which is a different and false statement.
        nan = float("nan")
        return ShareResult(kappa=nan, kappa_adj=nan, ci_low=nan, ci_high=nan,
                           null_mean=nan, p_value=nan, n=n_retained,
                           n_groups=n_groups_all,
                           prevalence=float(X.mean()), support=support,
                           missing_share=missing_share,
                           estimable=False, cv_sd=nan,
                           n_dropped_non_finite=n_dropped_non_finite,
                           n_dropped_untyped=n_dropped_untyped, **structure)

    kappa, cv_sd, draws = _repeated_skill(X, code, folds=folds, repeats=repeats,
                                          rng=rng)

    perms = [rng.permutation(code) for _ in range(n_perm)]
    null = np.array([_skill(X, perm, _folds(n, folds, rng)) for perm in perms],
                    dtype=float)
    # The penalty is read from the first fold draw of each permutation, so the
    # debiased estimate and the interval below do not depend on how many
    # draws the p-value averages over.
    c_null = float(np.nanmean(null)) if null.size else float("nan")
    kappa_adj = _debias(kappa, c_null)
    # A second stream for the extra permuted draws, jumped far ahead of the
    # main one and taken here, before the bootstrap, so that it depends on
    # the isolate positions alone and not on which lineages the bootstrap
    # happens to pick.
    more = np.random.Generator(_jumped(rng))

    boot = []
    order, bounds = _group_index(code)
    for _ in range(n_boot):
        idx, newcode = _cluster_resample_fast(order, bounds, rng)
        if idx.size == 0 or np.ptp(X[idx], axis=0).max() <= 0:
            continue
        boot.append(_debias(
            _skill(X[idx], newcode, _folds(len(idx), folds, rng)), c_null))
    draws_b = np.asarray([b for b in boot if np.isfinite(b)], dtype=float)
    lo, hi = (float(np.percentile(draws_b, 2.5)), float(np.percentile(draws_b, 97.5))) \
        if draws_b.size >= 20 else (float("nan"), float("nan"))

    # The second question, on its own stream so that every number above is
    # unchanged by its presence: a within-lineage resample for the noise of
    # the estimate, then the chi-square layer for the draw of the lineages.
    within_rng = np.random.Generator(_jumped(more))
    within = []
    for _ in range(n_boot):
        idx = _stratified_resample(code, within_rng)
        if np.ptp(X[idx], axis=0).max() <= 0:
            continue
        within.append(_debias(
            _skill(X[idx], code[idx], _folds(n, folds, within_rng)), c_null))
    # The estimate and its draws are lineage-membership shares; the chi-square
    # law is written for the component ratio, one design factor away on the
    # odds scale, so both are moved onto that scale before the inversion.
    design = float(1.0 - np.sum((sizes / n) ** 2))
    sp_share = float(_to_component_scale(kappa_adj, design))
    sp_lo, sp_hi = _superpopulation_interval(
        sp_share,
        _to_component_scale(np.asarray(within, dtype=float), design),
        int(code.max()) + 1, within_rng)
    # The chi-square layer describes lineage effects drawn from one law; a
    # collection in which one lineage carries most of the resistance is not
    # that, and there the lineage bootstrap, which sees the influence of that
    # lineage, is the wider statement. The species interval is therefore the
    # envelope of the two, the lineage bootstrap taken on the same scale: it
    # can only cover more often than either.
    if np.isfinite(sp_lo) and np.isfinite(lo):
        lo_c = float(_to_component_scale(max(lo, 0.0), design))
        hi_c = float(_to_component_scale(min(hi, 1.0), design))
        sp_lo = float(min(sp_lo, lo_c))
        sp_hi = float(max(sp_hi, hi_c))

    # Extra fold draws for the permuted statistics come last and from their
    # own stream, so that every quantity above is the same number whatever
    # ``null_repeats`` is. The observed side of the comparison is the mean
    # of as many draws as the permuted side.
    extra = max(1, min(int(null_repeats), repeats))
    if extra > 1 and n_perm:
        null = np.array([
            np.nanmean([null[i]] + [_skill(X, perm, _folds(n, folds, more))
                                    for _ in range(extra - 1)])
            for i, perm in enumerate(perms)], dtype=float)
    observed = float(np.nanmean(draws[:extra]))
    # A permuted statistic that could not be scored (no variance held out)
    # is not a non-exceedance; it leaves the comparison, and the floor of
    # the p-value follows the permutations that remain.
    null = null[np.isfinite(null)]
    n_perm_scored = int(null.size)
    exceed = int((null >= observed).sum())

    prevalence = float(np.asarray(X, dtype=float).mean())
    dominant = float("nan")
    if X.shape[1] == 1:
        # The largest single-lineage share of the between-lineage sum of
        # squares, read from the whole cohort: the diagnostic for a carrier
        # structure the chi-square layer does not describe.
        y = X[:, 0]
        sums = np.bincount(code, weights=y)
        means = sums / sizes
        contrib = sizes * (means - y.mean()) ** 2
        # Summed in sorted order so that the figure does not move by a unit
        # in the last place when the lineages are renamed.
        total = float(np.sort(contrib).sum())
        dominant = float(contrib.max() / total) if total > 0 else float("nan")
    latent: typing.Dict[str, typing.Any] = dict(
        latent_share=float("nan"), latent_low=float("nan"),
        latent_high=float("nan"), latent_species_low=float("nan"),
        latent_species_high=float("nan"))
    if X.shape[1] == 1 and np.isin(X, (0.0, 1.0)).all():
        # One binary trait: the same numbers on the liability scale, a
        # monotone transformation at the observed prevalence, treated as fixed.
        from .latent import latent_interval, latent_share as _latent
        latent["latent_share"] = _latent(kappa_adj, prevalence)
        latent["latent_low"], latent["latent_high"] = latent_interval(
            lo, hi, prevalence)
        latent["latent_species_low"], latent["latent_species_high"] = (
            latent_interval(sp_lo, sp_hi, prevalence))

    return ShareResult(
        kappa=kappa, kappa_adj=kappa_adj, ci_low=lo, ci_high=hi,
        null_mean=c_null,
        p_value=_phipson_smyth(exceed, n_perm_scored),
        n=n_retained, n_groups=n_groups_all,
        prevalence=prevalence,
        support=support, missing_share=missing_share,
        estimable=bool(np.isfinite(kappa_adj) and sizes.size >= MIN_REPEATED_LINEAGES),
        cv_sd=cv_sd,
        n_dropped_non_finite=n_dropped_non_finite,
        n_dropped_untyped=n_dropped_untyped,
        superpopulation_low=sp_lo, superpopulation_high=sp_hi,
        superpopulation_share=sp_share if np.isfinite(sp_lo) else float("nan"),
        p_floor=(1.0 / (n_perm_scored + 1.0)) if n_perm_scored else float("nan"),
        dominant_lineage_share=dominant,
        few_lineages=bool(sizes.size < FEW_LINEAGES),
        **structure,
        reason=_reason(c_null, lo, hi, n_boot, structure["n_groups_repeated"]),
        **latent,
    )


def clonal_share(y, lineage, **kwargs) -> ShareResult:
    """Per-agent clonal share. Thin alias of :func:`layer_clonal_share` for a
    single binary trait, kept separate because it is the quantity a
    surveillance table reports and the one the documentation names."""
    y = np.asarray(y, dtype=float).reshape(-1, 1)
    return layer_clonal_share(y, lineage, **kwargs)
