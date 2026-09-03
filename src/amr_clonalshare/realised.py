"""The share of trait variance carried by the lineages a cohort actually holds.

ESTIMAND, BEFORE ESTIMATOR. Two different questions share one point estimate and
must not share one interval.

*Superpopulation.* If a new set of lineages were drawn from this species, what
share of the trait variance would sit between them? The parameter is
``tau^2 / (tau^2 + sigma^2)``, and the information about ``tau^2`` carries
``G - 1`` degrees of freedom whatever the number of isolates. A cohort with
thirty lineages therefore has a wide interval and no estimator repairs that:
`attribution` and `censored` report this question and say so.

*Realised.* In this collection, holding these lineages, what share of the trait
variance sits between them? The lineage effects are then fixed unknowns rather
than a fresh draw, and the parameter is ``S_a^2 / (S_a^2 + sigma^2)`` with
``S_a^2`` their realised dispersion. This module reports that question.

The second is what a laboratory holding one collection usually means, and it is
identified far more sharply, because nothing is being extrapolated to lineages
that were never seen. Under normality the between-lineage sum of squares is a
noncentral chi-square whose noncentrality is exactly the quantity wanted, so the
interval is obtained by inverting the noncentral F and is exact, including under
unequal lineage sizes, where the variance-component interval is only
approximate.

WHAT WOULD MEAN THIS MEASURED THE WRONG THING. Each interval must be scored
against its own truth. Measured at thirty lineages of twenty isolates and a true
share of 0.30 over 3000 replicates: the variance-component interval covers its
own target 0.951 of the time and the realised target 1.000; the realised
interval covers its own target 0.950 and the superpopulation target 0.658. Had
the second row not fallen away from nominal, the two estimands would not be
distinguishable and this module would have no reason to exist.

WHAT WAS TRIED AND DISCARDED, recorded here rather than deleted. A Box-type
deflation of the degrees of freedom by the estimated residual kurtosis was meant
to repair heavy tails; it made coverage worse, not better, at every tail weight
tested (0.872 to 0.040 at t(4), 0.941 to 0.487 at t(10)), and is not used. A
within-lineage percentile bootstrap was meant to provide an assumption-free
fallback; it covered 0.786 under the very process the method is derived against,
so it is not a fallback either. What remains is the exact interval together with
a gate on the residual kurtosis, whose limit is read off the measured coverage
curve rather than chosen.

WHAT THE GATE IS AND IS NOT, measured conditional on its own decision (90
cells, 10000 replicates each; validation ledger V38). Given an open gate the
interval covers 0.947 to 0.957 under the Gaussian law and 0.929 to 0.961 on
binary traits at prevalence 0.25, the shortfall sitting at high shares. On a
heavy-tailed continuous process (t on four degrees of freedom) the gate
refuses 88 to 100 percent of cohorts with twenty or more isolates per lineage,
and the few it admits cover as low as 0.36: a sample-kurtosis gate selects the
light-tailed draws of a heavy-tailed process. The gate is therefore a refusal
rule, not a certificate, and an open gate on a continuous trait says the
sample looked Gaussian, not that the interval holds. A Bernoulli trait passes
the gate only for prevalence between 0.173 and 0.827; outside that band the
refusal is by construction and the withheld interval would have covered 0.86
to 0.96.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, Sequence, Tuple

import numpy as np
from scipy.stats import f as _f_dist
from scipy.stats import ncf as _ncf
from scipy.stats import kurtosis as _kurtosis

__all__ = ["KURTOSIS_LIMIT", "MIN_GROUPS", "RealisedShare", "realised_share",
           "realised_interval", "superpopulation_interval"]

#: Above this estimated excess kurtosis of the within-lineage residuals the
#: exact interval loses its level, and the run is marked not estimable. The
#: limit is the largest median excess kurtosis at which coverage of the
#: realised share stayed at or above 0.93 on the curve measured in
#: ``benchmarks/realised_calibration.py``, and is read off it rather
#: than chosen. Below it the measured coverage runs 0.93 to 0.95;
#: at t(4) residuals, whose fourth moment does not exist, it is 0.87.
KURTOSIS_LIMIT = 0.99

#: Two populated lineages give one between-lineage degree of freedom, which is
#: the least that identifies anything at all.
MIN_GROUPS = 2

_MAX_NONCENTRALITY = 1.0e7

#: Label used for isolates with no lineage assignment, spelled as
#: ``attribution._codes`` spells it so that the two modules read one cohort the
#: same way.
_MISSING = "__missing__"


def _label_strings(labels: np.ndarray) -> np.ndarray:
    """Lineage labels in the string form every other module compares them by.

    ``np.unique`` sorts, and sorting a vector that mixes an integer ``1`` with
    a text ``"1"`` -- which is what a cohort assembled from two pipelines looks
    like -- raises a bare ``TypeError``, as does a NaN sitting among strings.
    Missing values (None, NaN, empty, and the spellings ``nan``/``na``/``none``)
    collapse to one level here, exactly as :func:`attribution._codes`
    documents it, rather than becoming several levels or an error.
    """
    out = np.empty(labels.shape, dtype=object)
    for i, v in enumerate(labels.ravel()):
        out.ravel()[i] = _MISSING if (v is None or v != v
                                      or str(v).strip() == ""
                                      or str(v).lower() in ("nan", "na", "none")
                                      ) else str(v)
    return out.astype(str)


@dataclass
class RealisedShare:
    """The realised share, its exact interval, and the wider question beside it."""

    kappa: float
    ci_low: float
    ci_high: float
    #: The variance-component interval for the superpopulation question, on the
    #: same data. Reported beside the realised one so that a reader can see what
    #: the extrapolation to unseen lineages costs, and never in place of it.
    superpopulation_low: float
    superpopulation_high: float
    n: int
    n_groups: int
    #: The unbalance-corrected effective lineage size, equal to the common size
    #: on a balanced design.
    effective_group_size: float
    f_ratio: float
    #: Excess kurtosis of the within-lineage residuals, which is what the gate
    #: reads. Zero for a Gaussian.
    residual_excess_kurtosis: float
    estimable: bool
    reason: str

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def _gate_residuals(residual: np.ndarray, counts: np.ndarray,
                    code: np.ndarray) -> np.ndarray:
    """Residuals the kurtosis gate may read.

    A lineage of one contributes a residual of exactly zero, carrying no
    information about the shape of the within-lineage law while pulling a
    kurtosis towards a spike at the origin: on a cohort with many singletons the
    gate would then be reading the lineage size distribution rather than the
    trait. Singletons are therefore dropped, and what remains is divided by
    ``sqrt(1 - 1/n_g)``, the factor by which fitting a group mean shrinks the
    residual of a group of that size, so that groups of different sizes are on
    one scale before their shapes are compared.
    """
    keep = counts[code] >= 2
    if not np.any(keep):
        return np.empty(0)
    sizes = counts[code][keep]
    return residual[keep] / np.sqrt(1.0 - 1.0 / sizes)


def _one_way(y: np.ndarray, code: np.ndarray, n_groups: int):
    """Mean squares of the one-way layout, with unequal group sizes allowed."""
    counts = np.bincount(code, minlength=n_groups).astype(float)
    sums = np.bincount(code, weights=y, minlength=n_groups)
    present = counts > 0
    means = np.zeros(n_groups)
    means[present] = sums[present] / counts[present]
    n = int(y.size)
    populated = int(present.sum())
    grand = float(y.mean())
    ssb = float(np.sum(counts[present] * (means[present] - grand) ** 2))
    residual = y - means[code]
    ssw = float(np.sum(residual ** 2))
    return counts, present, populated, n, ssb, ssw, residual


def _effective_size(counts: np.ndarray, present: np.ndarray, n: int,
                    populated: int) -> float:
    """Satterthwaite's ``n0``: the common size a balanced design would need."""
    if populated < 2:
        return float("nan")
    return (n - float(np.sum(counts[present] ** 2)) / n) / (populated - 1)


def realised_interval(f_ratio: float, df_between: float, df_within: float,
                      effective_size: float, n_groups: int,
                      alpha: float = 0.05) -> Tuple[float, float]:
    """Invert the noncentral F for the realised share.

    Holding the lineage effects fixed, ``SSB / sigma^2`` is a noncentral
    chi-square on ``G - 1`` degrees of freedom whose noncentrality is
    ``lambda = sum_g n_g (a_g - abar)^2 / sigma^2``, independent of the central
    ``SSW / sigma^2``. Their ratio of mean squares is therefore a noncentral F,
    and the set of noncentralities the observed ratio does not reject is a
    confidence interval for ``lambda``. The share is a monotone function of it,
    so the endpoints carry across.
    """
    if not np.isfinite(f_ratio) or f_ratio <= 0.0 or n_groups < MIN_GROUPS:
        return float("nan"), float("nan")
    if df_between < 1.0 or df_within < 1.0 or not np.isfinite(effective_size):
        return float("nan"), float("nan")

    def solve(target: float) -> float:
        """Smallest noncentrality whose distribution puts ``target`` below the
        observed ratio. Zero when even a central F is already below it."""
        if _ncf.cdf(f_ratio, df_between, df_within, 0.0) < target:
            return 0.0
        high = 1.0
        while (_ncf.cdf(f_ratio, df_between, df_within, high) > target
               and high < _MAX_NONCENTRALITY):
            high *= 2.0
        if _ncf.cdf(f_ratio, df_between, df_within, high) > target:
            return float(high)
        low = 0.0
        for _ in range(80):
            mid = 0.5 * (low + high)
            if _ncf.cdf(f_ratio, df_between, df_within, mid) > target:
                low = mid
            else:
                high = mid
        return 0.5 * (low + high)

    lam_low = solve(1.0 - alpha / 2.0)
    lam_high = solve(alpha / 2.0)
    scale = effective_size * (n_groups - 1)
    to_share = lambda lam: lam / (lam + scale) if scale > 0 else float("nan")
    low, high = sorted((to_share(lam_low), to_share(lam_high)))
    return float(np.clip(low, 0.0, 1.0)), float(np.clip(high, 0.0, 1.0))


def superpopulation_interval(f_ratio: float, df_between: float,
                             df_within: float, effective_size: float,
                             alpha: float = 0.05) -> Tuple[float, float]:
    """The variance-component interval, by inverting the central F.

    Kept here so that the two questions can be reported side by side from one
    fit. It is the interval `censored` reports, restated for uncensored input.
    """
    if not np.isfinite(f_ratio) or f_ratio <= 0.0:
        return float("nan"), float("nan")
    if df_between < 1.0 or df_within < 1.0 or not np.isfinite(effective_size):
        return float("nan"), float("nan")
    out = []
    for q in (1.0 - alpha / 2.0, alpha / 2.0):
        adjusted = f_ratio / _f_dist.ppf(q, df_between, df_within)
        out.append((adjusted - 1.0) / (adjusted - 1.0 + effective_size))
    low, high = sorted(out)
    return float(np.clip(low, 0.0, 1.0)), float(np.clip(high, 0.0, 1.0))


def realised_share(y: Sequence[float], lineage: Sequence,
                   *, alpha: float = 0.05,
                   kurtosis_limit: float = KURTOSIS_LIMIT) -> RealisedShare:
    """The share of trait variance carried by the lineages this cohort holds.

    ``y`` may be binary or continuous; a binary trait is read on its own scale,
    so the share is the intraclass correlation of the indicator rather than a
    liability-scale quantity. For a susceptibility panel read at its recorded
    resolution use :mod:`amr_clonalshare.censored`, which reports the same
    two intervals from the interval likelihood.

    Lineage labels are compared by their string form, so that a label arriving
    as an integer from one pipeline and as text from another is one lineage.
    """
    if not 0.0 < float(alpha) < 1.0:
        # alpha = 1 returned a zero-width interval, alpha = 2 a reversed one
        # and alpha = nan the pair (0, 0), each of them flagged estimable.
        raise ValueError(f"alpha={alpha}; the interval needs 0 < alpha < 1")
    y = np.asarray(y, dtype=float).ravel()
    labels = np.asarray(lineage, dtype=object).ravel()
    if y.size != labels.size:
        raise ValueError(f"y has {y.size} entries and lineage has "
                         f"{labels.size}; they must match")
    keep = np.isfinite(y)
    y, labels = y[keep], labels[keep]
    if y.size == 0:
        return RealisedShare(float("nan"), float("nan"), float("nan"),
                             float("nan"), float("nan"), 0, 0, float("nan"),
                             float("nan"), float("nan"), False,
                             "no finite observations")
    _, code = np.unique(_label_strings(labels), return_inverse=True)
    n_groups = int(code.max()) + 1
    counts, present, populated, n, ssb, ssw, residual = _one_way(y, code,
                                                                 n_groups)
    n0 = _effective_size(counts, present, n, populated)

    def refuse(reason: str) -> RealisedShare:
        return RealisedShare(float("nan"), float("nan"), float("nan"),
                             float("nan"), float("nan"), n, populated,
                             n0, float("nan"), float("nan"), False, reason)

    if populated < MIN_GROUPS:
        return refuse(f"{populated} populated lineage(s); at least "
                      f"{MIN_GROUPS} are needed for one between-lineage "
                      f"degree of freedom")
    if n <= populated:
        return refuse(f"{n} isolates over {populated} lineages leaves no "
                      f"within-lineage degrees of freedom")
    df_between = float(populated - 1)
    df_within = float(n - populated)
    msb = ssb / df_between
    msw = ssw / df_within
    if msw <= 0.0:
        return refuse("every lineage is constant, so the within-lineage scale "
                      "is zero and no share is identified")
    f_ratio = msb / msw
    between = max((msb - msw) / n0, 0.0)
    kappa = between / (between + msw)
    scaled = _gate_residuals(residual, counts, code)
    excess = (float(_kurtosis(scaled, fisher=True, bias=False))
              if scaled.size >= 4 else float("nan"))

    low, high = realised_interval(f_ratio, df_between, df_within, n0,
                                  populated, alpha)
    slow, shigh = superpopulation_interval(f_ratio, df_between, df_within, n0,
                                           alpha)
    if not (np.isfinite(low) and np.isfinite(high)):
        # f_ratio is zero when every lineage mean coincides. The inversion
        # then has no interval to return, and a result with no interval must
        # not be read as estimable whatever the gate says.
        estimable = False
        reason = ("the between-lineage sum of squares is zero, so the "
                  "noncentral-F inversion has no interval to return")
    elif not np.isfinite(excess):
        estimable = False
        reason = ("fewer than four residuals come from lineages of two or "
                  "more, so the shape of the within-lineage law cannot be "
                  "read and the gate cannot be applied")
    elif excess > kurtosis_limit:
        estimable = False
        reason = (f"within-lineage residuals have excess kurtosis "
                  f"{excess:.2f} against a limit of {kurtosis_limit:.2f}, "
                  f"above which the exact interval was measured to lose its "
                  f"level")
    else:
        estimable, reason = True, ""
    return RealisedShare(float(kappa), low, high, slow, shigh, n, populated,
                         float(n0), float(f_ratio), excess, estimable, reason)
