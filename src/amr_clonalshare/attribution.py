"""attribution.py — how much of a resistance pattern is the clone?

Contents
--------
``clonal_share``
    Out-of-sample share of one binary trait's variance that a lineage label
    explains. The per-agent number a surveillance report can carry.
``layer_clonal_share``
    The same quantity for a block of traits taken together.
``attribute_partition``
    Commonality and Shapley decomposition of a partition's explanatory power
    between the lineage factor and the partition factor.
``concordance_z``
    The same-lineage pair concordance z the interpretation gate uses, kept here
    so the two can be reported side by side.

Why this module exists
----------------------
A partition of resistance profiles is checked for lineage confounding by
permuting labels and asking whether same-lineage pairs co-cluster more often
than chance. That statistic answers "is there any lineage signal", and on a
real cohort the answer is always yes, because resistance is partly inherited.
It is computed over the O(n^2) pairs of the cohort, so its null standard
deviation shrinks with n while the effect it measures does not: on the shipped
*S. suis* example the same structure gives z = 4.7 at n = 67 and z = 23.2 at
n = 677. A gate on that statistic refuses every cohort large enough to be worth
analysing, and refusing is not a result.

The quantity that does not scale with n is a variance share. This module
estimates it out of sample, so a lineage variable with many levels is not paid
for its levels, and reports it with an interval and against a permutation null.

Reading the number
------------------
``kappa`` near 1: non-wild-type status for this agent travels with the clone.
The lever is biosecurity, movement control, cleaning and mixing.
``kappa`` near 0: it travels independently of the clone, which in practice
means a mobile element or repeated independent selection. The lever is
selection pressure, dosing and choice of agent.

That is the same fork the prevalence decomposition draws between a composition
change and a within-lineage rate change, applied to the level rather than to a
difference between two collections.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

#: Share of isolates that must sit in a lineage with at least two members
#: before a clonal share is flagged estimable. Read off the coverage curve in
#: ``benchmarks/attribution_calibration.py`` (cell D), not chosen: over 180
#: cells the 95% interval covers the truth 0.00 to 0.40 of the time below a
#: support of 0.75, 0.65 at 0.80, and 1.00 at 0.90 and above, while the bias
#: falls from -0.18 to -0.02 across the same range. The curve is built on the
#: worst case, every lineage outside the large ones being a singleton, so a
#: real cohort with a graded size distribution clears it more easily than this
#: threshold implies. See ``SUPPORT_THRESHOLD_EVIDENCE``.
SUPPORT_THRESHOLD = 0.90
SUPPORT_THRESHOLD_EVIDENCE = "benchmarks/results_attribution_calibration/support_curve.json"

#: Label used for isolates with no lineage assignment. They form one level
#: rather than being dropped, because dropping them changes the estimand.
_MISSING = "__missing__"

__all__ = [
    "SUPPORT_THRESHOLD",
    "ShareResult",
    "AttributionResult",
    "clonal_share",
    "layer_clonal_share",
    "attribute_partition",
    "concordance_z",
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
    ci_low, ci_high : percentile bootstrap interval over isolates, on the
        corrected scale.
    null_mean : mean kappa when the predictor's labels are permuted. Negative
        rather than zero, and its magnitude is the penalty ``c`` above.
    p_value : permutation p-value, Phipson-Smyth corrected.
    n, n_groups, prevalence : cohort description carried with the estimate.
    support : share of isolates sitting in a lineage with at least two members.
        A singleton lineage can never be predicted out of sample -- its one
        isolate is either in the training fold or in the held-out fold, never
        both -- so it contributes only to the denominator. A cohort typed at a
        resolution that makes most lineages singletons therefore carries a
        kappa that is mostly a statement about the typing scheme.
    estimable : ``support >= SUPPORT_THRESHOLD``. Below it the estimate is
        returned, because suppressing a number is not the same as reporting an
        honest one, but the flag is false and callers are expected to gate on
        it. The threshold is set from the coverage curve in
        ``benchmarks/attribution_calibration.py``, not chosen.
    cv_sd : spread of kappa across repeated fold assignments; a diagnostic of
        the estimator, not of the cohort, and never a confidence interval.
    n_dropped_non_finite : isolates removed before the estimate because a trait
        column was not finite. A missing well is normal in a susceptibility
        panel, and every comparison against a NaN is false, so leaving those
        rows in made the permuted runs never exceed the observed skill and
        drove the p-value to its floor. The count is carried because the
        estimate then describes a subset of the cohort.
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

    def as_dict(self) -> dict:
        return {f.name: getattr(self, f.name) for f in self.__dataclass_fields__.values()}


@dataclass(frozen=True)
class AttributionResult:
    """Split of a partition's explanatory power against a lineage label.

    All three R-squared terms are out-of-sample and are scored on the *same*
    fold assignment, so the commonality identity holds to within Monte Carlo
    error rather than being disturbed by three different fold draws.

    Attributes
    ----------
    r2_lineage, r2_partition, r2_joint : out-of-sample variance shares of the
        lineage design, the partition design, and the additive model of both.
        All three are debiased on their own permutation null, because the three
        designs have different numbers of levels and an undebiased split would
        report how finely the cohort was typed rather than what the partition
        carries. ``r2_joint`` is the best member of the additive family on the
        same folds, i.e. at least the better marginal.
    shared : ``r2_lineage + r2_partition - r2_joint``. The part of the
        partition's power that the lineage label could equally have supplied.
        Negative when the lineage label predicts held-out traits worse than the
        cohort prevalence does, which is a real state for a typing scheme whose
        levels are mostly singletons.
    lam : ``shared / r2_partition``, clipped to [0, 1]; the lineage-attributable
        share of the partition. ``lam_clipped`` records whether clipping bit.
    shapley_lineage, shapley_partition : Shapley shares of ``r2_joint`` between
        the two factors (Shorrocks 1982; Israeli 2007). They sum to
        ``r2_joint`` by construction and, unlike ``lam``, need no clipping.
        This is *not* SHAP: it is the two-player Shapley value of the
        variance-explained game, computed exactly, not approximated.
    partition_refit : whether the partition was rebuilt inside every fold
        rather than taken as one labelling of the whole cohort. With a fixed
        labelling the fold split holds isolates out of the group means but not
        out of the clustering that drew the groups, so a held-out isolate had a
        hand in deciding which group it would then be scored against. On the
        shipped *S. suis* cohort that is worth 44 per cent of ``r2_partition``
        and 0.10 of ``lam``, which is most of the distance to the gate; on
        *Klebsiella* it is worth 4 per cent and moves ``lam`` the other way,
        because a two-cluster partition of 1500 isolates barely changes when a
        fifth of them are withheld. The size of it tracks how far the partition
        moves when a fold is dropped, so it cannot be read off the cohort size
        and has to be measured on each cohort.
    lam_cv_sd : spread of ``lam`` across the repeated fold assignments. A
        diagnostic of the estimator, not of the cohort, and never a confidence
        interval, on the same reading as ``ShareResult.cv_sd``. Under a refit
        it is the only dispersion reported, because the lineage bootstrap
        cannot be run there for the reason :func:`attribute_partition` gives.
    r2_partition_assignment_control, lam_assignment_control : the same two
        quantities with the clustering left fitted on every isolate but the
        held-out isolates placed by the refit's own assignment rule. Choosing a
        held-out isolate's cluster from its traits and then asking how well
        that cluster predicts those traits is itself worth something: it raised
        ``r2_partition`` by 28 per cent on *Klebsiella* with the partition held
        fixed. The control separates that channel, so a reader can tell how
        much of the change under a refit came from the clustering moving and
        how much from the placement. Both are NaN when no refit was done.
    """

    r2_lineage: float
    r2_partition: float
    r2_joint: float
    shared: float
    lam: float
    lam_clipped: bool
    shapley_lineage: float
    shapley_partition: float
    ci_low: float
    ci_high: float
    n: int
    missing_share: float = 0.0
    partition_refit: bool = False
    lam_cv_sd: float = float("nan")
    r2_partition_assignment_control: float = float("nan")
    lam_assignment_control: float = float("nan")

    def as_dict(self) -> dict:
        return {f.name: getattr(self, f.name) for f in self.__dataclass_fields__.values()}


# ----------------------------------------------------------------- internals
def _rng(seed) -> np.random.Generator:
    if isinstance(seed, np.random.Generator):
        return seed
    return np.random.default_rng(seed)


def _codes(labels) -> np.ndarray:
    """Integer codes for any hashable label vector, missing values included.

    Missing lineage is coded as its own level rather than dropped. A cohort
    where typing failed non-randomly is a cohort where dropping the untyped
    isolates changes the estimand, and the caller is told to decide that
    explicitly by passing a subset if that is what they mean.
    """
    arr = np.asarray(labels, dtype=object)
    clean = np.empty(arr.shape, dtype=object)
    for i, v in enumerate(arr.ravel()):
        clean.ravel()[i] = _MISSING if (v is None or v != v
                                        or str(v).strip() == ""
                                        or str(v).lower() in ("nan", "na", "none")
                                        ) else str(v)
    _, codes = np.unique(clean.astype(str), return_inverse=True)
    return codes.astype(np.int64)


def _missing_share(labels) -> float:
    """Share of isolates whose lineage label is absent.

    Reported rather than gated. An untyped isolate joins a single
    ``__missing__`` level, and if typing failed informatively -- as it does on
    the shipped *S. suis* cohort, where the untyped isolates carry about two
    more non-wild-type results out of thirteen -- that level carries real
    signal and the attribution partly measures the typing process. The number
    is emitted so a reader can see how much of the estimate rests on it, and so
    that the same cohort typed two ways can be compared.
    """
    arr = np.asarray(labels, dtype=object).ravel()
    return float(np.mean([(v is None or v != v or str(v).strip() == ""
                           or str(v).lower() in ("nan", "na", "none"))
                          for v in arr]))


def _cross(designs: Sequence[np.ndarray]) -> np.ndarray:
    """Codes for the crossing of several code vectors.

    Retained for callers that want the saturated model. It is **not** what the
    commonality decomposition uses: crossing 299 sequence types with 2 clusters
    asks for 598 cells on 1031 isolates, and out of sample that predicts worse
    than either factor alone, which then appears in the arithmetic as shared
    variance and drove one Shapley share negative. See :func:`_skill_additive`.
    """
    key = np.zeros(len(designs[0]), dtype=np.int64)
    for g in designs:
        key = key * (int(g.max()) + 1) + g
    _, out = np.unique(key, return_inverse=True)
    return out.astype(np.int64)


def _group_means(X, code, G):
    counts = np.bincount(code, minlength=G).astype(float)
    sums = np.empty((G, X.shape[1]))
    for j in range(X.shape[1]):
        sums[:, j] = np.bincount(code, weights=X[:, j], minlength=G)
    return np.where(counts[:, None] > 0, sums / np.maximum(counts[:, None], 1), 0.0), counts


def _folds(n: int, k: int, rng: np.random.Generator) -> np.ndarray:
    order = rng.permutation(n)
    fold = np.empty(n, dtype=np.int64)
    fold[order] = np.arange(n) % k
    return fold


def _code_for(code, f: int) -> np.ndarray:
    """The partition design that applies to fold ``f``.

    A single array is one design fixed for the whole cohort, which is what a
    caller who already holds a partition passes. A sequence of arrays is one
    design per fold, which is what refitting the clustering inside each fold
    produces, and there the held-out fold has to be scored against the design
    that was built without it. Resolving the two shapes here rather than at
    every call site keeps the scorers identical in every other respect, so a
    refitted run and a fixed-partition run differ only in where the labels came
    from and the commonality identity still holds on one fold assignment.
    """
    return code if isinstance(code, np.ndarray) else code[int(f)]


def _permute_code(code, rng: np.random.Generator):
    """Permute a partition design, whether it is one design or one per fold.

    A single array is permuted directly, which is what the fixed-partition path
    has always done and is why that path still draws the same numbers from the
    generator in the same order. A per-fold design is permuted through one
    shared index instead, because the folds label the same isolates and
    permuting each of them on its own would scramble that correspondence and
    measure a penalty for a design nothing is scored against.
    """
    if isinstance(code, np.ndarray):
        return rng.permutation(code)
    pi = rng.permutation(len(code[0]))
    return [c[pi] for c in code]


def _skill(X: np.ndarray, code, fold: np.ndarray) -> float:
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

    ``code`` is either one design for the whole cohort or one design per fold;
    see :func:`_code_for` for why the second shape exists. The number of groups
    is therefore read inside the fold loop, which for a single design gives the
    same value on every fold and leaves the score exactly as it was.
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
        c = _code_for(code, f)
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


def _skill_additive(X, code_a, code_b, fold, n_iter: int = 4) -> float:
    """Out-of-sample skill of the **additive** two-factor model.

    Predicts ``mu + a[A_i] + b[B_i]`` rather than a mean per crossed cell,
    fitted by backfitting on the training rows: take group means for A, take
    group means of the residual for B, repeat. Four passes are enough here
    because two factors converge geometrically and the folds add more noise
    than the fifth pass removes.

    The additive model is the right joint term for a commonality split. The
    saturated crossing is a different question -- whether the partition means
    something *different inside each lineage* -- and on cohorts with many
    lineages it cannot be estimated at all.

    ``code_b`` follows :func:`_code_for` and may be one design per fold, which
    is what the partition becomes when it is refitted inside the fold. The
    lineage factor ``code_a`` never varies that way, because a lineage label is
    an observation about the isolate rather than something the run fits.
    """
    X = np.atleast_2d(np.asarray(X, dtype=float))
    if X.shape[0] != code_a.shape[0]:
        X = X.T
    Ga = int(code_a.max()) + 1
    pred = np.empty_like(X)
    base = np.empty_like(X)
    for f in np.unique(fold):
        tr, te = fold != f, fold == f
        if not tr.any() or not te.any():
            continue
        cb = _code_for(code_b, f)
        Gb = int(cb.max()) + 1
        Xtr, atr, btr = X[tr], code_a[tr], cb[tr]
        mu = Xtr.mean(axis=0)
        a = np.zeros((Ga, X.shape[1]))
        b = np.zeros((Gb, X.shape[1]))
        for _ in range(n_iter):
            a, _c = _group_means(Xtr - mu - b[btr], atr, Ga)
            b, _c = _group_means(Xtr - mu - a[atr], btr, Gb)
        pred[te] = mu + a[code_a[te]] + b[cb[te]]
        base[te] = mu
    sse = float(((X - pred) ** 2).sum())
    sst = float(((X - base) ** 2).sum())
    return float("nan") if sst <= 0 else 1.0 - sse / sst


def _repeated_skill(X, code, *, folds, repeats, rng):
    n = np.asarray(code).shape[0]
    vals = np.array([_skill(X, code, _folds(n, folds, rng)) for _ in range(repeats)])
    return float(np.nanmean(vals)), float(np.nanstd(vals))


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
    """Draw lineages with replacement from a prepared index. See
    :func:`_cluster_resample` for why lineages are taken whole."""
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


def _cluster_resample(code: np.ndarray, rng: np.random.Generator):
    """Two-stage bootstrap: draw lineages with replacement, then take their
    isolates.

    The isolate bootstrap is wrong here and undercovers badly. A variance share
    across lineages is uncertain mainly because a cohort contains few lineages,
    not because it contains few isolates: resampling isolates leaves every
    lineage present at close to its original size and so holds fixed the very
    thing the estimate depends on. A lineage drawn twice becomes two distinct
    groups, which is the standard convention and correctly reduces the
    effective information rather than duplicating it.

    Lineages are taken **whole**. Resampling isolates within a drawn lineage as
    well puts duplicate rows into the cohort, and a duplicated row that lands in
    a training fold and in the held-out fold leaks its own value into its own
    prediction. On a cohort with many small lineages that leak is not small: it
    lifted every bootstrap replicate above the point estimate, so the interval
    sat entirely above the quantity it was supposed to cover.

    Returns ``(row_index, new_codes)``.
    """
    groups = np.unique(code)
    members = [np.flatnonzero(code == g) for g in groups]
    picks = rng.integers(0, len(groups), len(groups))
    rows, new = [], []
    for j, gi in enumerate(picks):
        idx = members[gi]
        rows.append(idx)
        new.append(np.full(idx.size, j, dtype=np.int64))
    return np.concatenate(rows), np.concatenate(new)


def _debias(kappa: float, null_mean: float) -> float:
    """Remove the cross-validation penalty measured on the permuted design.

    ``kappa_hat = kappa - (1 - kappa) c`` and ``null_mean = -c``, so the map
    below returns ``kappa`` exactly under that model. The denominator is
    ``1 + c >= 1``, so the correction can never blow up; it is a rescaling, not
    an extrapolation.
    """
    if not np.isfinite(kappa) or not np.isfinite(null_mean):
        return float("nan")
    return float((kappa - null_mean) / (1.0 - null_mean))


# ------------------------------------------------------------- public: agent
def layer_clonal_share(X, lineage, *, folds: int = 5, repeats: int = 20,
                       n_boot: int = 400, n_perm: int = 200,
                       seed=None) -> ShareResult:
    """Out-of-sample share of a trait block's variance explained by lineage.

    Parameters
    ----------
    X : (n,) or (n, p) array of binary trait calls, 1 = non-wild-type.
    lineage : (n,) label vector; MLST sequence type, a BAPS cluster, a clonal
        complex, anything categorical. It is used as a label and never as a
        tree, which is the point: a laboratory that types by MLST has one.
    folds, repeats : cross-validation design. ``repeats`` averages over fold
        draws so the estimate does not depend on one partition of the cohort.
    n_boot : percentile bootstrap replicates over isolates.
    n_perm : permutations of the lineage label for the null and the p-value.

    Returns
    -------
    ShareResult
    """
    if int(folds) < 2:
        # One fold holds nothing out, so every score is taken on an empty
        # held-out set and the run returns kappa 0, a zero-width interval,
        # p = 1 and estimable true -- a clean "lineage explains nothing"
        # verdict that nothing was ever tested for. Zero folds reaches a
        # modulo by zero first.
        raise ValueError(f"folds={folds}; cross-validation needs at least two "
                         f"folds, because one fold holds nothing out")
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
    finite = np.isfinite(X).all(axis=1)
    n_dropped_non_finite = int((~finite).sum())
    if n_dropped_non_finite:
        # Every comparison against a NaN is false, so `null >= kappa` counted
        # no exceedances and the permutation p-value collapsed to its floor
        # 1 / (n_perm + 1) -- the most significant value obtainable -- beside a
        # NaN kappa_adj and estimable true. A missing well is normal in a
        # susceptibility panel, so the row goes and is counted.
        X, lineage = X[finite], lineage[finite]
    code = _codes(lineage)
    n = X.shape[0]
    if n < 2:
        # Fewer than two isolates leave nothing to hold out and nothing to
        # compare, so there is no share to estimate.
        nan = float("nan")
        return ShareResult(kappa=nan, kappa_adj=nan, ci_low=nan, ci_high=nan,
                           null_mean=nan, p_value=nan, n=n,
                           n_groups=int(code.max()) + 1 if code.size else 0,
                           prevalence=float(X.mean()) if X.size else nan,
                           support=nan,
                           missing_share=(_missing_share(lineage) if n
                                          else nan),
                           estimable=False, cv_sd=nan,
                           n_dropped_non_finite=n_dropped_non_finite)
    sizes = np.bincount(code)
    support = float((sizes[code] >= 2).mean())

    if float(np.ptp(X, axis=0).max()) <= 0:
        # A trait present in every isolate, or in none, has no variance to
        # attribute. Returning NaN says so; returning 0 would read as "not
        # clonal", which is a different and false statement.
        nan = float("nan")
        return ShareResult(kappa=nan, kappa_adj=nan, ci_low=nan, ci_high=nan,
                           null_mean=nan, p_value=nan, n=n,
                           n_groups=int(code.max()) + 1,
                           prevalence=float(X.mean()), support=support,
                           missing_share=_missing_share(lineage),
                           estimable=False, cv_sd=nan,
                           n_dropped_non_finite=n_dropped_non_finite)

    kappa, cv_sd = _repeated_skill(X, code, folds=folds, repeats=repeats, rng=rng)

    null = np.array([_skill(X, rng.permutation(code), _folds(n, folds, rng))
                     for _ in range(n_perm)], dtype=float)
    exceed = int((null >= kappa).sum())
    c_null = float(np.nanmean(null))
    kappa_adj = _debias(kappa, c_null)

    boot = []
    order, bounds = _group_index(code)
    for _ in range(n_boot):
        idx, newcode = _cluster_resample_fast(order, bounds, rng)
        if idx.size == 0 or np.ptp(X[idx], axis=0).max() <= 0:
            continue
        boot.append(_debias(
            _skill(X[idx], newcode, _folds(len(idx), folds, rng)), c_null))
    boot = np.asarray([b for b in boot if np.isfinite(b)], dtype=float)
    lo, hi = (float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))) \
        if boot.size >= 20 else (float("nan"), float("nan"))

    return ShareResult(
        kappa=kappa, kappa_adj=kappa_adj, ci_low=lo, ci_high=hi,
        null_mean=c_null,
        p_value=_phipson_smyth(exceed, n_perm),
        n=n, n_groups=int(code.max()) + 1,
        prevalence=float(np.asarray(X, dtype=float).mean()),
        support=support, missing_share=_missing_share(lineage),
        estimable=bool(support >= SUPPORT_THRESHOLD),
        cv_sd=cv_sd,
        n_dropped_non_finite=n_dropped_non_finite,
    )


def clonal_share(y, lineage, **kwargs) -> ShareResult:
    """Per-agent clonal share. Thin alias of :func:`layer_clonal_share` for a
    single binary trait, kept separate because it is the quantity a
    surveillance table reports and the one the documentation names."""
    y = np.asarray(y, dtype=float).reshape(-1, 1)
    return layer_clonal_share(y, lineage, **kwargs)


# --------------------------------------------------------- public: partition
def _attribution_penalties(X, lab, lin, *, folds, n_perm, rng):
    """Cross-validation penalty of each of the three designs, from permutation.

    Each design pays for its own levels, and the three designs here do not have
    the same number: a 299-level sequence type against a 2-cluster partition.
    Forming a commonality split out of the raw scores would therefore make the
    answer a function of how finely the cohort was typed -- the same defect the
    module exists to remove from the gate statistic. Every term is debiased on
    its own null before any arithmetic is done with it.

    ``lab`` may be one design per fold, in which case it is permuted through
    one shared index so that the folds go on describing the same isolates. The
    penalty is a property of the level count and the fold geometry rather than
    of the labels, so a refitted partition and a fixed one measure nearly the
    same penalty and the difference between the two is not an artefact of this
    correction.
    """
    n = len(lin)
    p_l, p_c, p_j = [], [], []
    for _ in range(n_perm):
        f = _folds(n, folds, rng)
        p_l.append(_skill(X, rng.permutation(lin), f))
        p_c.append(_skill(X, _permute_code(lab, rng), f))
        p_j.append(_skill_additive(X, rng.permutation(lin),
                                   _permute_code(lab, rng), f))
    return (float(np.nanmean(p_l)), float(np.nanmean(p_c)), float(np.nanmean(p_j)))


def _one_attribution(X, lab, lin, fold, pen=(0.0, 0.0, 0.0)):
    r_l = _debias(_skill(X, lin, fold), pen[0])
    r_c = _debias(_skill(X, lab, fold), pen[1])
    r_fit = _debias(_skill_additive(X, lin, lab, fold), pen[2])
    # The additive family contains both single-factor models, so the best it can
    # do is at least the better marginal. A fitted backfit can still land below
    # that when one factor is noise, and letting it would turn estimation noise
    # into apparent shared variance; selecting within the family on the same
    # folds is the fix, and it also bounds lambda above by one.
    r_lc = float(np.nanmax([r_fit, r_l, r_c]))
    shared = r_l + r_c - r_lc
    raw = shared / r_c if r_c > 0 else float("nan")
    lam = float(np.clip(raw, 0.0, 1.0)) if np.isfinite(raw) else float("nan")
    clipped = bool(np.isfinite(raw) and (raw < 0.0 or raw > 1.0))
    return dict(r2_lineage=r_l, r2_partition=r_c, r2_joint=r_lc,
                r2_joint_fitted=r_fit, shared=shared,
                lam=lam, lam_clipped=clipped,
                shapley_lineage=0.5 * (r_l + (r_lc - r_c)),
                shapley_partition=0.5 * (r_c + (r_lc - r_l)))


def _fold_codes(refit, fold, folds, labels=None):
    """One partition of the whole cohort per fold, built by the caller's refit.

    Each fold gets its own labelling because the clustering is rebuilt without
    the isolates that fold holds out, and those isolates are then placed
    against the clusters the training rows produced. Passing ``labels`` asks
    the callable to skip the rebuild and place the held-out isolates against a
    partition that already exists, which is how the assignment control is
    built: it changes the placement and nothing else.
    """
    out = []
    for f in range(folds):
        train = np.flatnonzero(fold != f)
        held_out = np.flatnonzero(fold == f)
        out.append(_codes(refit(train, held_out) if labels is None
                          else refit(train, held_out, labels=labels)))
    return out


def attribute_partition(X, labels, lineage, *, folds: int = 5, repeats: int = 20,
                        n_boot: int = 400, n_perm: int = 60, refit=None,
                        seed=None) -> AttributionResult:
    """Split a partition's explanatory power between lineage and itself.

    The partition is taken as given. This function does not ask whether it is a
    good partition; it asks how much of what it explains about the traits a
    lineage label already explained. A ``lam`` near 1 says the partition is a
    lineage relabelling and should not be read as an archetype. A ``lam`` near 0
    says it is carrying something the clone does not.

    ``lam_clipped`` true means crossing lineage with the partition predicted
    held-out traits *worse* than lineage alone. That is not a numerical
    nuisance: it says the partition adds nothing to lineage on these traits,
    and the honest reading is ``lam = 1``.

    ``refit`` makes the partition fold-internal, and it is what makes the
    out-of-sample claim in ``r2_partition`` true. Without it the partition is
    one labelling of the whole cohort, so the fold split holds isolates out of
    the group means but not out of the clustering that defined the groups, and
    every held-out isolate had a hand in deciding which group it would be
    scored against. The callable is handed the training and held-out row
    indices and must return a label for every isolate in the cohort, so it
    carries the assignment rule as well as the clustering; when it is given
    ``labels`` as well, it must place the held-out isolates against that
    partition without rebuilding one, which is what the assignment control
    asks of it.

    The placement is not free, and the result says so rather than folding it
    away. Choosing a held-out isolate's cluster from its own traits and then
    scoring how well that cluster predicts those same traits raised
    ``r2_partition`` by 28 per cent on the shipped *Klebsiella* cohort with the
    clustering held fixed, so ``r2_partition_assignment_control`` is computed
    beside the refitted figure and the two are meant to be read together.

    ``n_boot`` must be zero under a refit. The lineage bootstrap draws lineages
    whole, so a lineage drawn twice puts the same isolates into the cohort
    twice, and refitting inside such a replicate would train a cluster on rows
    that also sit in the held-out fold. That is the leak
    :func:`_cluster_resample` already refuses to create within a lineage.
    ``lam_cv_sd`` is reported in its place; it is a Monte Carlo spread and not
    an interval, but the gate has always been decided by the point estimate, so
    nothing that carries a threshold is lost.
    """
    if int(folds) < 2:
        # As in layer_clonal_share: one fold holds nothing out and zero folds
        # divides by zero, and neither refusal was ever raised.
        raise ValueError(f"folds={folds}; cross-validation needs at least two "
                         f"folds, because one fold holds nothing out")
    rng = _rng(seed)
    X = np.atleast_2d(np.asarray(X, dtype=float))
    if X.shape[0] != len(labels):
        if X.shape[1] == len(labels):
            X = X.T
        else:
            # Transposing on any mismatch analyses a cohort of one or two rows
            # rather than saying that the two inputs do not describe the same
            # isolates.
            raise ValueError(f"X has {X.shape[0]} rows and labels has "
                             f"{len(labels)} entries; they must match")
    lab, lin = _codes(labels), _codes(lineage)
    n = X.shape[0]

    if refit is None:
        # Left exactly as it was, so a caller who already holds a partition
        # draws the same numbers from the generator in the same order.
        pen = _attribution_penalties(X, lab, lin, folds=folds, n_perm=n_perm,
                                     rng=rng)
        reps = [_one_attribution(X, lab, lin, _folds(n, folds, rng), pen)
                for _ in range(repeats)]
        control = {}
    else:
        if n_boot:
            raise ValueError(
                f"n_boot={n_boot} together with refit; the lineage bootstrap "
                f"draws lineages whole, so a lineage drawn twice puts the same "
                f"isolates in the cohort twice and refitting inside the "
                f"replicate would train a cluster on rows that are also held "
                f"out. Pass n_boot=0 and read lam_cv_sd instead")
        draws = [_folds(n, folds, rng) for _ in range(repeats)]
        codes = [_fold_codes(refit, fd, folds) for fd in draws]
        placed = [_fold_codes(refit, fd, folds, labels=lab) for fd in draws]
        pen = _attribution_penalties(X, codes[0], lin, folds=folds,
                                     n_perm=n_perm, rng=rng)
        reps = [_one_attribution(X, c, lin, fd, pen)
                for c, fd in zip(codes, draws)]
        ctrl = [_one_attribution(X, c, lin, fd, pen)
                for c, fd in zip(placed, draws)]
        control = {
            "r2_partition_assignment_control":
                float(np.nanmean([r["r2_partition"] for r in ctrl])),
            "lam_assignment_control":
                float(np.nanmean([r["lam"] for r in ctrl])),
        }
    agg = {k: float(np.nanmean([r[k] for r in reps]))
           for k in reps[0] if k != "lam_clipped"}
    agg["lam_clipped"] = bool(np.mean([r["lam_clipped"] for r in reps]) > 0.5)
    agg.pop("r2_joint_fitted", None)
    agg["lam_cv_sd"] = float(np.nanstd([r["lam"] for r in reps]))
    agg["partition_refit"] = refit is not None
    agg.update(control)

    boot = []
    order, bounds = _group_index(lin)
    for _ in range(n_boot):
        idx, newlin = _cluster_resample_fast(order, bounds, rng)
        if idx.size == 0:
            continue
        r = _one_attribution(X[idx], _codes(lab[idx]), newlin,
                             _folds(len(idx), folds, rng), pen)
        if np.isfinite(r["lam"]):
            boot.append(r["lam"])
    boot = np.asarray(boot, dtype=float)
    lo, hi = (float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))) \
        if boot.size >= 20 else (float("nan"), float("nan"))

    return AttributionResult(ci_low=lo, ci_high=hi, n=n,
                             missing_share=_missing_share(lineage), **agg)


def concordance_z(labels, lineage, *, n_perm: int = 500, seed=None):
    """Same-lineage pair concordance and its permutation z.

    Reported alongside the attribution so that a reader can see both the
    statistic the gate used and the magnitude it failed to carry. Returns
    ``(observed, null_mean, null_sd, z, p_value)``.
    """
    rng = _rng(seed)
    lab, lin = _codes(labels), _codes(lineage)
    iu = np.triu_indices(len(lab), 1)
    same_lin = lin[iu[0]] == lin[iu[1]]
    if same_lin.sum() == 0:
        return (float("nan"),) * 5
    obs = float((lab[iu[0]] == lab[iu[1]])[same_lin].mean())
    null = np.empty(n_perm)
    for b in range(n_perm):
        pl = rng.permutation(lab)
        null[b] = (pl[iu[0]] == pl[iu[1]])[same_lin].mean()
    sd = float(null.std())
    z = (obs - float(null.mean())) / sd if sd > 0 else float("inf")
    return obs, float(null.mean()), sd, z, _phipson_smyth(int((null >= obs).sum()), n_perm)
