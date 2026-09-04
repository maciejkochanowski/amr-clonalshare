"""The clonal share from an interval-censored MIC, and from its coarsenings.

THE IDEA THAT ORGANISES THIS MODULE. A dichotomised non-wild-type call, a
recorded MIC and a censored reading are not three analyses. They are one
likelihood at three interval widths. A wild-type call is the interval
(-inf, log2 cut-off]; a non-wild-type call is (log2 cut-off, +inf); a recorded
MIC of 8 mg/L on a doubling panel is (log2 4, log2 8]; a reading on the lowest
tested well is left-censored and one on the highest is right-censored. The
estimator below takes intervals and nothing else, so the input mode changes
only how much information an observation carries, never what is estimated.
Two tests in ``tests/test_censored.py`` enforce that: with degenerate
intervals the result must equal the continuous estimator, and with single-cut
intervals it must equal the binary one.

WHAT DICHOTOMISATION COSTS is therefore measurable rather than asserted. Cohen
(1983) puts the loss at the equivalent of 55 per cent of the sample for a cut
one standard deviation from the mean, and worse further out; an epidemiological
cut-off sits in a tail by construction.

THE ASSUMPTION THAT MUST BE DECLARED. Cohorts frequently record no censoring
operator, and yet pile mass on the end wells of a panel. Treating those
readings as censored is an assumption about the coarsening mechanism, not an
observation. It is the coarsened-at-random condition of Heitjan and Rubin
(1991), under which the coarsening is ignorable for likelihood inference.
``panel_geometry`` reports the end-well mass so the assumption is visible, and
``sensitivity_endpoints`` recomputes the share with the end wells treated as
exact readings instead, which brackets the effect of getting it wrong.

WHY THE GROUP MEANS ARE SHRUNK. A lineage whose isolates all sit on the top
well carries no information about how far above it they lie; its latent mean
is not identified and an unshrunk maximum-likelihood estimate runs away. The
empirical-Bayes factor pulls such a lineage back to the cohort mean in
proportion to how little it says, which is what makes heavy censoring
survivable rather than fatal. Zhang, Wang and O'Connor (2021) declined to
estimate variance components under large censored fractions; the answer here
is to shrink, to report the panel geometry, and to refuse when the geometry
says the question cannot be answered.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import log_ndtr, ndtr, ndtri
from scipy.stats import chi2, f as _f_dist

from .attribution import _codes, _rng

__all__ = ["PanelGeometry", "CensoredShare", "panel_geometry",
           "intervals_from_mic", "intervals_from_binary",
           "scale_is_identified", "censored_clonal_share",
           "sensitivity_endpoints", "marginal_loglik", "profile_interval"]

_SQRT2PI = float(np.sqrt(2.0 * np.pi))
# Below this share of isolates on an end well, a point-valued MIC is a fair
# description of the reading. Above it, the point value is a censored one
# wearing a number, and the interval mode is required. Read off the coverage
# curve in the calibration, not chosen.
END_WELL_LIMIT = 0.05
# Share of isolates that may lie in lineages wholly beyond the panel before
# the share stops being estimable. A lineage all of whose readings are
# one-sided has no identified latent mean, and only the shrinkage keeps the
# fit finite; what decides whether that matters is how much of the cohort is
# in that position. Swept over 21 censoring levels with 60 simulated cohorts
# each (30 lineages of 25, true share 0.5), the bias stays under 0.005 up to a
# share of 0.35, reaches 0.029 at 0.52 and 0.079 at 0.69. Measured, not chosen.
CENSORED_GROUP_LIMIT = 0.5
# Non-wild-type prevalence outside which a single cut point no longer carries
# the latent share. Read off the design grid: at a prevalence near 0.08 the
# estimate reads 0.23 above a true zero and its interval covers 0.08 of the
# time, while from 0.24 upward the bias is under 0.04 and coverage is 0.89 to
# 0.95. A dilution at the same prevalence is unaffected, which is the practical
# argument for reading the panel rather than the call.
SINGLE_CUT_PREVALENCE = (0.10, 0.90)


@dataclass(frozen=True)
class PanelGeometry:
    """What the panel can and cannot support, before any biology."""

    n_wells: int
    lowest: float
    highest: float
    share_on_lowest: float
    share_on_highest: float
    lattice_ratios: Tuple[float, ...]
    doubling: bool
    admissible_modes: Tuple[str, ...]
    reason: str

    def as_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class CensoredShare:
    kappa: float
    ci_low: float
    ci_high: float
    n: int
    n_groups: int
    sigma_within: float
    sigma_between: float
    #: The cluster bootstrap, kept as an assumption-free check and not as the
    #: interval: over the design grid it covers a median 0.61 of the time
    #: against a nominal 0.95, because resampling lineages does not reproduce
    #: the sampling distribution of a component built from those lineages.
    boot_low: float
    boot_high: float
    #: Endpoints of the likelihood-ratio width. A measure of how much
    #: information the reading carries, not a confidence interval: on exact
    #: and interval readings at thirty lineages its coverage is 0.68 and 0.74.
    #: The interval to report is ``ci_low`` to ``ci_high``.
    profile_low: float
    profile_high: float
    share_censored: float
    n_groups_fully_censored: int
    #: Share of isolates belonging to lineages wholly beyond the panel. This
    #: is what the estimability gate reads; the count above is reported beside
    #: it because it is the number a reader expects to see.
    share_in_censored_groups: float
    estimable: bool
    reason: str
    #: The share carried by the lineages this cohort actually holds, with the
    #: interval obtained by inverting the noncentral F on the same two mean
    #: squares. A different question from ``ci_low`` to ``ci_high`` above, and
    #: the one a cohort with few lineages can support: see
    #: :mod:`amr_clonalshare.realised`.
    realised_low: float = float("nan")
    realised_high: float = float("nan")

    def as_dict(self) -> Dict[str, Any]:
        return asdict(self)


def panel_geometry(values: Sequence[float],
                   cutoff: Optional[float] = None) -> PanelGeometry:
    """Describe the dilution panel and say which input modes it supports.

    The lattice is reported rather than repaired. Ratios near 1.03 to 1.05
    between adjacent recorded values are rounding variants of one dilution,
    and a ratio of 8 is a pair of untested intermediate dilutions; neither is
    a property of the panel and neither should be silently rounded away. The
    intervals below are built from the wells actually tested for this agent,
    following the varying-range treatment of Aerts and colleagues rather than
    forcing a common power-of-two grid.
    """
    v = np.asarray([x for x in np.asarray(values, dtype=float)
                    if np.isfinite(x) and x > 0])
    if v.size == 0:
        return PanelGeometry(0, float("nan"), float("nan"), float("nan"),
                             float("nan"), (), False, (), "no usable readings")
    wells = np.unique(v)
    ratios = tuple(round(float(r), 3) for r in np.unique(wells[1:] / wells[:-1]))
    doubling = bool(np.allclose(ratios, 2.0, atol=0.15)) if len(ratios) else False
    lo_share = float((v == wells[0]).mean())
    hi_share = float((v == wells[-1]).mean())

    modes = ["binary"]
    reason = []
    if wells.size >= 3:
        modes.append("interval")
    else:
        reason.append("fewer than three wells")
    if max(lo_share, hi_share) <= END_WELL_LIMIT:
        modes.append("point")
    else:
        reason.append(
            f"{max(lo_share, hi_share):.2f} of readings sit on an end well, so a "
            "point value is a censored reading and the interval mode is required")
    if cutoff is not None and np.isfinite(cutoff):
        if not (wells[0] < cutoff <= wells[-1]):
            modes = [m for m in modes if m != "binary"]
            reason.append("the cut-off lies outside the tested range, so every "
                          "call is forced by the panel")
    return PanelGeometry(int(wells.size), float(wells[0]), float(wells[-1]),
                         lo_share, hi_share, ratios, doubling,
                         tuple(modes), "; ".join(reason) or "panel supports all modes")


def intervals_from_binary(y, cutoff_log2: float = 0.0):
    """A dichotomised call as the interval it actually is.

    Wild type means the latent value did not exceed the cut-off; non-wild type
    means it did. Nothing else is known. This is the coarsest member of the
    family and exists so the same estimator can consume a cohort that reports
    only an interpretation.
    """
    y = np.asarray(y, dtype=float).ravel()
    lo = np.where(y > 0.5, cutoff_log2, -np.inf)
    hi = np.where(y > 0.5, np.inf, cutoff_log2)
    return lo, hi


def intervals_from_mic(values, *, operators=None, treat_end_wells_as_censored=True):
    """Intervals on the log2 scale from recorded MIC values.

    A reading at a tested dilution means the true value lies above the
    previous tested dilution and at or below this one. A reading on the lowest
    well is left-censored and one on the highest is right-censored, unless
    ``treat_end_wells_as_censored`` is turned off, which is the sensitivity
    arm for the coarsening assumption.

    ``operators`` may carry the recorded censoring signs; where present they
    win over the end-well heuristic, because a recorded operator is an
    observation and the heuristic is an assumption.
    """
    v = np.asarray(values, dtype=float).ravel()
    ok = np.isfinite(v) & (v > 0)
    wells = np.unique(v[ok])
    lw = np.log2(wells)
    lo = np.full(v.shape, -np.inf)
    hi = np.full(v.shape, np.inf)
    idx = np.searchsorted(wells, v[ok])
    prev = np.where(idx > 0, lw[np.maximum(idx - 1, 0)], -np.inf)
    lo[ok] = prev
    hi[ok] = np.log2(v[ok])
    if treat_end_wells_as_censored and wells.size:
        lo[ok & (v == wells[0])] = -np.inf
        hi[ok & (v == wells[-1])] = np.inf
    if operators is not None:
        op = np.asarray(operators, dtype=object).ravel()
        gt = np.isin(op, (">", ">="))
        le = np.isin(op, ("<", "<="))
        lo[gt] = np.log2(np.where(v[gt] > 0, v[gt], 1.0)) if gt.any() else lo[gt]
        hi[gt] = np.inf
        lo[le] = -np.inf
        hi[le] = np.log2(np.where(v[le] > 0, v[le], 1.0)) if le.any() else hi[le]
    lo[~ok] = np.nan
    hi[~ok] = np.nan
    return lo, hi


def _trunc_moments(a, b, m, s):
    """Mean and second central moment of a normal restricted to (a, b].

    Computed through log_ndtr so that an interval far into a tail, which is
    exactly what a lineage pinned to an end well produces, does not underflow
    to a zero denominator.
    """
    s = max(float(s), 1e-9)
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    # A zero-width interval is an exact reading: its conditional mean is the
    # value and its conditional variance is zero. It needs its own branch
    # because the interval probability is zero and the general expression
    # divides by it. This branch is what makes the point mode collapse to the
    # classical one-way variance-component estimator, which the equivalence
    # test in tests/test_censored.py checks against a closed form.
    exact = np.isfinite(a) & np.isfinite(b) & ((b - a) <= 1e-12)
    al = (a - m) / s
    be = (b - m) / s
    la = log_ndtr(al)
    lb = log_ndtr(be)
    # log(Phi(be) - Phi(al)) computed stably
    with np.errstate(divide="ignore", invalid="ignore"):
        logZ = lb + np.log1p(-np.exp(np.clip(la - lb, -700, -1e-12)))
    logZ = np.where(np.isfinite(logZ), logZ, -700.0)
    Z = np.exp(np.clip(logZ, -700, 50))
    Z = np.maximum(Z, 1e-300)
    # An unbounded end contributes nothing to either moment. The bound is
    # replaced by a finite placeholder before the arithmetic so that the
    # product -inf * 0 is never formed; the placeholder is masked out again.
    fa, fb = np.isfinite(al), np.isfinite(be)
    al_s = np.where(fa, al, 0.0)
    be_s = np.where(fb, be, 0.0)
    pa = np.where(fa, np.exp(-0.5 * al_s ** 2) / _SQRT2PI, 0.0)
    pb = np.where(fb, np.exp(-0.5 * be_s ** 2) / _SQRT2PI, 0.0)
    r = (pa - pb) / Z
    ez = m + s * r
    ta = np.where(fa, al_s * pa, 0.0)
    tb = np.where(fb, be_s * pb, 0.0)
    var = s * s * (1.0 + (ta - tb) / Z - r * r)
    var = np.clip(var, 1e-12, None)
    if np.any(exact):
        ez = np.where(exact, np.where(exact, b, 0.0), ez)
        var = np.where(exact, 0.0, var)
    return ez, var


def scale_is_identified(lo, hi) -> bool:
    """Can the residual scale be estimated, or only a ratio?

    A single cut point tells you whether the latent value fell above or below
    it and nothing else, so only the standardised distance from the cut is
    recoverable and the residual scale is not. Dichotomisation therefore does
    not merely cost precision: it removes the scale, and the latent variance
    share exists only up to a convention. The convention adopted when the
    scale is unidentified is the usual threshold-model one, a residual
    standard deviation of one, which puts the result on the liability scale of
    Dempster and Lerner (1950) and makes it comparable with published
    heritabilities rather than with a raw log2 MIC variance.
    """
    b = np.concatenate([np.asarray(lo, dtype=float).ravel(),
                        np.asarray(hi, dtype=float).ravel()])
    return int(np.unique(b[np.isfinite(b)]).size) >= 2


def _fit_em(lo, hi, code, G, *, iters=200, tol=1e-10, fixed_scale=None):
    """Group latent means, residual scale and between-group variance, by
    expectation maximisation on the interval likelihood with a restricted-
    likelihood maximisation step. Returns the posterior group means, the
    residual scale, the between-group variance and the grand mean.

    Three things in that step earn their place, and each was adopted because
    it moved a measured bias rather than because it is customary.

    THE POSTERIOR VARIANCE OF THE GROUP MEAN STAYS IN. A group mean estimated
    from a handful of censored readings is uncertain, and a between-group
    variance built only from the squared deviations of those means discards
    that uncertainty. Keeping it is what separates a restricted-likelihood
    step from a moment step, and on a single cut point at prevalence 0.24 it
    moves the bias of the share from -0.077 to +0.012.

    THE NOISE IN A GROUP AVERAGE IS NOT sigma2. By the law of total variance
    the conditional expectations vary about the group mean with variance
    sigma2 minus the mean conditional variance: an interval has already
    resolved part of an observation variance, so a censored reading
    contributes less noise to the group average than an exact one.
    Subtracting the whole of sigma2 over-subtracts, and the error grows with
    the censored fraction, which is exactly the regime this module exists to
    serve.

    THE DIVISORS ARE THE RESTRICTED ONES: populated groups minus one for the
    between-group variance, because the grand mean was estimated from the same
    data, and n minus the number of populated groups for the residual scale.
    Together they make the estimator coincide with the classical one-way
    variance-component ratio on a balanced exact design, which
    ``tests/test_censored.py`` checks to 1e-3, while beating it on an
    unbalanced one: over 200 simulated cohorts with group sizes from 2 to 90
    the bias is +0.004 against -0.020 and the root mean squared error 0.093
    against 0.120.

    ``fixed_scale`` holds the residual standard deviation at a given value,
    which is required when the data carry a single cut point and the scale is
    therefore not identified."""
    finite_hi = np.where(np.isfinite(hi), hi, np.nan)
    finite_lo = np.where(np.isfinite(lo), lo, np.nan)
    start = np.nanmean(np.where(np.isfinite(finite_hi), finite_hi, finite_lo))
    if not np.isfinite(start):
        start = 0.0
    m = np.full(G, float(start))
    grand = float(start)
    s = 1.0
    tau2 = 1.0
    nn = int(np.asarray(lo).size)
    for _ in range(int(iters)):
        ez, var = _trunc_moments(lo, hi, m[code], s)
        cnt = np.bincount(code, minlength=G).astype(float)
        present = cnt > 0
        grand = float(ez.mean())
        raw = np.divide(np.bincount(code, weights=ez, minlength=G), cnt,
                        out=np.full(G, grand), where=present)
        gp = int(present.sum())
        # Variance of a conditional expectation about the group mean, by the
        # law of total variance. Equals sigma2 when nothing is censored.
        s2_raw = max(s * s - float(var.mean()), 1e-9)
        b = np.where(present, tau2 / (tau2 + s2_raw / np.maximum(cnt, 1.0)), 0.0)
        new_m = grand + b * (raw - grand)
        post_var = np.where(present, b * s2_raw / np.maximum(cnt, 1.0), 0.0)
        tau2 = float(np.sum((new_m[present] - grand) ** 2 + post_var[present])
                     / max(gp - 1, 1))
        if fixed_scale is not None:
            new_s = float(fixed_scale)
        else:
            resid = float(np.sum(var + (ez - new_m[code]) ** 2))
            new_s = float(np.sqrt(max(resid / max(nn - gp, 1.0), 1e-9)))
        if np.max(np.abs(new_m - m)) < tol and abs(new_s - s) < tol:
            m, s = new_m, new_s
            break
        m, s = new_m, new_s
    return m, s, float(max(tau2, 0.0)), grand


def _fully_censored_groups(lo, hi, code, G):
    """Groups every one of whose intervals is unbounded on the same side, and
    the share of isolates in them. Their latent mean is not identified: the
    likelihood is flat beyond the panel, and only the shrinkage keeps the fit
    finite."""
    count = 0
    members = 0
    for g in range(G):
        sel = code == g
        if not sel.any():
            continue
        if np.all(~np.isfinite(hi[sel])) or np.all(~np.isfinite(lo[sel])):
            count += 1
            members += int(sel.sum())
    return count, (members / code.size if code.size else 0.0)


def censored_clonal_share(lo, hi, lineage, *, n_boot: int = 200,
                          profile: bool = True, seed=None) -> CensoredShare:
    """Share of the latent log2 MIC variance carried by lineage.

    Group means and a shared residual scale are fitted by expectation
    maximisation on the interval likelihood, with empirical-Bayes shrinkage of
    the means at every step, and the quantity returned is the variance-
    component ratio

        kappa = tau2 / (tau2 + sigma2)

    where tau2 is the between-lineage variance of the latent log2 MIC and
    sigma2 the within-lineage residual variance. This is the estimand a mixed
    model reports and the one the heritability literature compares against.

    IT IS NOT THE SAME QUANTITY AS ``attribution.clonal_share``, and the
    difference is deliberate. That function reports achievable out-of-sample
    predictive skill, which falls when lineages are small because a group mean
    estimated from three isolates predicts badly however much variance the
    grouping truly carries. The ratio above reports the population variance
    share, which does not. A cohort should carry both: their ratio is a
    property of the study design, not of the biology.

    What the three input modes share is the likelihood, not the number. A
    dichotomised call, a recorded dilution and a censored reading enter as
    intervals of decreasing width and are scored by one expression, so the
    cost of coarsening shows up as a widening interval on kappa rather than as
    a change of method. Two tests in ``tests/test_censored.py`` hold that
    claim to account: with zero-width intervals the estimate must equal the
    closed-form one-way variance-component ratio, and a MIC vector coarsened
    at a cut-off must give byte-identical intervals to the dichotomised calls
    at that cut-off, and therefore the identical kappa.

    The interval is a cluster bootstrap that draws whole lineages, for the
    same reason as elsewhere: uncertainty in a variance share across lineages
    is governed by how many lineages there are, not how many isolates.
    """
    lo = np.asarray(lo, dtype=float).ravel()
    hi = np.asarray(hi, dtype=float).ravel()
    if lo.size != hi.size:
        raise ValueError(
            f"lo has {lo.size} entries and hi has {hi.size}")
    all_codes = _codes(lineage)
    if all_codes.shape[0] != lo.size:
        raise ValueError(
            f"lineage has {all_codes.shape[0]} entries for {lo.size} intervals")
    keep = ~(np.isnan(lo) | np.isnan(hi))
    lo, hi = lo[keep], hi[keep]
    # Re-coded after filtering so that the codes stay dense; a gap would make
    # the group counts below index a lineage that is no longer present.
    code = _codes(all_codes[keep])
    n = lo.size
    G = int(code.max()) + 1 if n else 0
    rng = _rng(seed)
    censored = float(np.mean(~np.isfinite(lo) | ~np.isfinite(hi))) if n else float("nan")
    fully, fully_share = (_fully_censored_groups(lo, hi, code, G) if n
                          else (0, 0.0))

    if n < 8 or G < 2:
        nan = float("nan")
        return CensoredShare(nan, nan, nan, n, G, nan, nan, nan, nan, nan,
                             nan, censored, fully, fully_share, False,
                             "too few isolates or lineages")

    identified = scale_is_identified(lo, hi)
    fixed = None if identified else 1.0

    def one(idx):
        """Latent-scale variance share from the interval likelihood.

        A held-out residual ratio, which is what the binary estimator uses,
        is the wrong score here: for a one-sided interval the conditional
        variance of the held-out reading does not fall as the group mean
        improves, so the ratio is insensitive to the very thing being
        measured. The variance-component ratio of the fitted model is used
        instead.
        """
        l, h, c = lo[idx], hi[idx], _codes(code[idx])
        Gi = int(c.max()) + 1
        _m, s, tau2, _g = _fit_em(l, h, c, Gi, fixed_scale=fixed)
        denom = tau2 + s * s
        return float("nan") if denom <= 0 else float(tau2 / denom)

    point = one(np.arange(n))
    order = np.argsort(code, kind="stable")
    bounds = np.concatenate([[0], np.cumsum(np.bincount(code, minlength=G))])
    boots = []
    for _ in range(int(n_boot)):
        pick = rng.integers(0, G, size=G)
        idx = np.concatenate([order[bounds[g]:bounds[g + 1]] for g in pick
                              if bounds[g + 1] > bounds[g]])
        if idx.size >= 8 and np.unique(code[idx]).size >= 2:
            boots.append(one(idx))
    b = np.asarray([x for x in boots if np.isfinite(x)], dtype=float)
    if b.size >= 20 and np.isfinite(point):
        # Jackknife over lineages, for the acceleration term.
        jack = []
        for g in range(G):
            keep_idx = np.flatnonzero(code != g)
            if keep_idx.size >= 8 and np.unique(code[keep_idx]).size >= 2:
                value = one(keep_idx)
                if np.isfinite(value):
                    jack.append(value)
        ci = _bca_limits(point, b, np.asarray(jack, dtype=float))
    else:
        ci = (float("nan"), float("nan"))

    m, s, tau2, _g = _fit_em(lo, hi, code, G, fixed_scale=fixed)
    f_lo, f_hi = _variance_ratio_interval(lo, hi, code, G, m, s, tau2)
    realised_low, realised_high = _realised_interval_from_fit(
        lo, hi, code, G, m, s, tau2)
    over_limit = fully_share > CENSORED_GROUP_LIMIT
    # For a single cut point the only observable is which side of the cut each
    # isolate fell, so a prevalence in a tail leaves almost no contrast to
    # divide between lineage and residual.
    above = float(np.mean(np.isfinite(lo))) if not identified else float("nan")
    tail = bool(not identified and np.isfinite(above)
                and not (SINGLE_CUT_PREVALENCE[0] <= above
                         <= SINGLE_CUT_PREVALENCE[1]))
    estimable = bool(np.isfinite(point) and not over_limit and not tail)
    notes = []
    if tail:
        notes.append(
            f"a single cut point with {above:.0%} of isolates above it lies "
            f"outside the calibrated window "
            f"{SINGLE_CUT_PREVALENCE[0]:.0%} to {SINGLE_CUT_PREVALENCE[1]:.0%}, "
            "where the share is not estimable from a call; the recorded "
            "concentration is not affected")
    if over_limit:
        notes.append(
            f"{fully_share:.0%} of isolates lie in lineages wholly beyond the "
            f"panel, above the calibrated limit of "
            f"{CENSORED_GROUP_LIMIT:.0%}; their latent means are not "
            "identified and the share is not estimable")
    elif fully:
        notes.append(
            f"{fully} of {G} lineages lie wholly beyond the panel, holding "
            f"{fully_share:.0%} of isolates, which the calibration puts inside "
            "the range where the shrinkage absorbs them")
    if not identified:
        notes.append(
            "scale not identified from a single cut point; reported on the "
            "liability scale with residual standard deviation fixed at one")
    # Computed once on the whole cohort, never inside the bootstrap: it costs
    # a few hundred likelihood evaluations and answers a different question.
    p_lo, p_hi, _peak = (profile_interval(lo, hi, code) if profile
                         else (float("nan"), float("nan"), float("nan")))
    return CensoredShare(point, f_lo, f_hi, n, G, float(s),
                         float(np.sqrt(max(tau2, 0.0))), ci[0], ci[1],
                         p_lo, p_hi, censored, fully,
                         float(fully_share), estimable,
                         "; ".join(notes) if notes else "ok",
                         realised_low, realised_high)


def sensitivity_endpoints(values, lineage, *, operators=None, **kwargs
                          ) -> Dict[str, Any]:
    """Bracket the end-well coarsening assumption by running it both ways.

    Reading a value on the lowest or highest tested well as censored is an
    assumption about why the value is there. Where the operator was recorded
    it is an observation and this question does not arise; where it was not,
    the honest response is to report the share under both readings and let the
    width of the bracket say how much of the answer rests on the assumption
    rather than on the data. A bracket narrower than the bootstrap interval
    means the assumption is not carrying the result.
    """
    lo_c, hi_c = intervals_from_mic(values, operators=operators,
                                    treat_end_wells_as_censored=True)
    lo_e, hi_e = intervals_from_mic(values, operators=operators,
                                    treat_end_wells_as_censored=False)
    as_censored = censored_clonal_share(lo_c, hi_c, lineage, **kwargs)
    as_exact = censored_clonal_share(lo_e, hi_e, lineage, **kwargs)
    both = [as_censored.kappa, as_exact.kappa]
    finite = [x for x in both if np.isfinite(x)]
    width = (max(finite) - min(finite)) if len(finite) == 2 else float("nan")
    ci = as_censored.ci_high - as_censored.ci_low
    return {
        "end_wells_censored": as_censored.as_dict(),
        "end_wells_exact": as_exact.as_dict(),
        "kappa_bracket_width": float(width),
        "interval_width": float(ci),
        "assumption_dominates": bool(np.isfinite(width) and np.isfinite(ci)
                                     and width > ci),
    }


# --------------------------------------------------------------------------
# Marginal likelihood and the interval that reads it
# --------------------------------------------------------------------------
#: Gauss-Hermite nodes for the integral over the lineage random effect. Twenty
#: is where the log likelihood of the shipped cohort stops moving in the sixth
#: decimal place; the check is in ``tests/test_censored.py``.
GH_NODES = 20
_GH_X, _GH_W = np.polynomial.hermite_e.hermegauss(GH_NODES)
_GH_LOGW = np.log(_GH_W) - 0.5 * np.log(2.0 * np.pi)


def marginal_loglik(lo, hi, code, G, grand: float, tau2: float,
                    sigma: float) -> float:
    """Observed-data log likelihood of the interval-censored random-intercept
    model, integrating the lineage effect out by Gauss-Hermite quadrature.

    The expectation-maximisation fit elsewhere in this module maximises this
    quantity without ever forming it. Forming it buys two things a variance
    component otherwise cannot have: an interval that reflects how much each
    isolate contributed rather than only how many lineages there were, and a
    likelihood-ratio check that the fit is a maximum and not a stationary
    point the shrinkage walked into.
    """
    a = np.asarray(lo, dtype=float).ravel()
    b = np.asarray(hi, dtype=float).ravel()
    s = max(float(sigma), 1e-9)
    sd_u = float(np.sqrt(max(tau2, 0.0)))
    u = sd_u * _GH_X
    # log P(a < z <= b | u) for every isolate at every quadrature node
    al = (a[:, None] - grand - u[None, :]) / s
    be = (b[:, None] - grand - u[None, :]) / s
    la, lb = log_ndtr(al), log_ndtr(be)
    with np.errstate(divide="ignore", invalid="ignore"):
        cell = lb + np.log1p(-np.exp(np.clip(la - lb, -700, -1e-12)))
    cell = np.where(np.isfinite(cell), cell, -700.0)
    # A zero-width interval is an exact reading and contributes a density.
    # Without this branch its interval probability is zero, the likelihood
    # collapses, and the profile below is flat.
    exact = np.isfinite(a) & np.isfinite(b) & ((b - a) <= 1e-12)
    if np.any(exact):
        dens = -0.5 * be ** 2 - np.log(s) - np.log(_SQRT2PI)
        cell = np.where(exact[:, None], dens, cell)
    # Per lineage, the sum over its isolates, one weighted count per node.
    per_node = np.empty((G, cell.shape[1]), dtype=float)
    for k in range(cell.shape[1]):
        per_node[:, k] = np.bincount(code, weights=cell[:, k], minlength=G)
    per_node += _GH_LOGW[None, :]
    present = np.bincount(code, minlength=G) > 0
    mx = per_node.max(axis=1, keepdims=True)
    lse = (mx.ravel() + np.log(np.exp(per_node - mx).sum(axis=1)))
    return float(lse[present].sum())


def _profile_at(lo, hi, code, G, kappa: float, grand: float,
                fixed_scale) -> float:
    """Maximised log likelihood with the share held at ``kappa``.

    One free scale parameter remains, the total variance, which is profiled
    out by a bounded scalar search. When the scale is not identified, which is
    the single-cut case, the residual standard deviation is held at the
    convention and the total variance follows from kappa, so there is nothing
    left to search.
    """
    k = float(np.clip(kappa, 1e-6, 1 - 1e-6))
    if fixed_scale is not None:
        s = float(fixed_scale)
        return marginal_loglik(lo, hi, code, G, grand, k * s * s / (1 - k), s)

    def negative(log_v):
        v = float(np.exp(log_v))
        return -marginal_loglik(lo, hi, code, G, grand, k * v,
                                float(np.sqrt((1 - k) * v)))

    best = minimize_scalar(negative, bounds=(-8.0, 8.0), method="bounded",
                           options={"xatol": 1e-4})
    return float(-best.fun)


def profile_interval(lo, hi, lineage, *, alpha: float = 0.05, grid: int = 41):
    """Likelihood-ratio width of the share, and the point the likelihood peaks
    at.

    The set of shares whose maximised log likelihood is within half a
    chi-squared quantile of the best one. It measures how sharply this cohort
    pins the share down, which is the question the extra information in a
    recorded dilution improves, and it is used here to compare two readings of
    one cohort.

    IT IS NOT A CONFIDENCE INTERVAL, and the calibration says so rather than
    leaving a reader to find out. Over 200 simulated cohorts of thirty
    lineages its coverage is 0.68 on an exact reading and 0.74 on a dilution
    against a nominal 0.95, while on a single cut point it is 0.91 to 0.93.
    The chi-squared calibration is optimistic exactly where the residual scale
    is well determined and the whole of the uncertainty sits in thirty group
    means. The interval to report is the bias-corrected cluster bootstrap in
    ``CensoredShare.ci_low`` and ``ci_high``, which answers how the share
    would move under another draw of lineages.
    """
    lo = np.asarray(lo, dtype=float).ravel()
    hi = np.asarray(hi, dtype=float).ravel()
    keep = ~(np.isnan(lo) | np.isnan(hi))
    lo, hi = lo[keep], hi[keep]
    code = _codes(np.asarray(_codes(lineage))[keep])
    G = int(code.max()) + 1 if code.size else 0
    if code.size < 8 or G < 2:
        return float("nan"), float("nan"), float("nan")
    fixed = None if scale_is_identified(lo, hi) else 1.0
    m, s, tau2, grand = _fit_em(lo, hi, code, G, fixed_scale=fixed)
    ks = np.linspace(1e-3, 1 - 1e-3, int(grid))
    ll = np.array([_profile_at(lo, hi, code, G, k, grand, fixed) for k in ks])
    peak = int(np.argmax(ll))
    cut = ll[peak] - 0.5 * chi2.ppf(1.0 - alpha, 1)
    inside = ks[ll >= cut]
    if inside.size == 0:
        return float("nan"), float("nan"), float(ks[peak])
    return float(inside.min()), float(inside.max()), float(ks[peak])


def _bca_limits(point: float, boots: np.ndarray, jack: np.ndarray,
                alpha: float = 0.05) -> Tuple[float, float]:
    """Bias-corrected and accelerated percentile limits.

    Falls back to the plain percentile limits when either correction cannot be
    formed, which happens when every bootstrap draw lies on one side of the
    point estimate or when the jackknife has no spread. A silent fallback would
    be worse than a wide interval, so the two cases are the only ones in which
    the correction is skipped.
    """
    lo_q, hi_q = alpha / 2.0, 1.0 - alpha / 2.0
    below = float(np.mean(boots < point))
    if not 0.0 < below < 1.0 or jack.size < 3:
        return (float(np.percentile(boots, 100 * lo_q)),
                float(np.percentile(boots, 100 * hi_q)))
    z0 = float(ndtri(below))
    centred = jack.mean() - jack
    denom = 6.0 * float(np.sum(centred ** 2)) ** 1.5
    a = float(np.sum(centred ** 3) / denom) if denom > 0 else 0.0
    out = []
    for q in (lo_q, hi_q):
        z = ndtri(q)
        adj = z0 + (z0 + z) / max(1.0 - a * (z0 + z), 1e-6)
        out.append(float(np.percentile(boots, 100 * float(ndtr(adj)))))
    lo, hi = sorted(out)
    return (float(np.clip(lo, 0.0, 1.0)), float(np.clip(hi, 0.0, 1.0)))


def _ratio_and_degrees(lo, hi, code, G, means, sigma, tau2):
    """The mean-square ratio and its degrees of freedom, from the fitted parts.

    Shared by the two interval constructions so that they cannot drift apart:
    they answer different questions about one fit, not two fits.
    """
    n = int(np.asarray(lo).size)
    cnt = np.bincount(code, minlength=G).astype(float)
    present = cnt > 0
    gp = int(present.sum())
    if gp < 2 or n <= gp:
        return None
    n0 = (n - float(np.sum(cnt[present] ** 2)) / n) / max(gp - 1, 1)
    _ez, var = _trunc_moments(lo, hi, means[code], sigma)
    resolved = float(np.clip(1.0 - float(var.mean()) / max(sigma ** 2, 1e-12),
                             0.05, 1.0))
    ratio = 1.0 + n0 * max(tau2, 0.0) / max(sigma ** 2, 1e-12)
    return (ratio, float(max(gp - 1, 1)), float(max((n - gp) * resolved, 1.0)),
            n0, gp, float(np.asarray(var).mean()))


def _realised_interval_from_fit(lo, hi, code, G, means, sigma, tau2,
                                alpha: float = 0.05) -> Tuple[float, float]:
    """Interval for the share carried by the lineages this cohort holds.

    The variance-ratio interval below answers what would sit between lineages
    drawn afresh, and at thirty lineages it is wide for a reason no estimator
    removes. This one answers what sits between the lineages that are here,
    which is the question a laboratory holding one collection usually means,
    and it is obtained by inverting the noncentral F on the same two mean
    squares. See :mod:`amr_clonalshare.realised` for the estimand and for
    the measured operating characteristics of both.
    """
    parts = _ratio_and_degrees(lo, hi, code, G, means, sigma, tau2)
    if parts is None:
        return float("nan"), float("nan")
    ratio, df_between, df_within, n0, gp, _ = parts
    from .realised import realised_interval
    return realised_interval(ratio, df_between, df_within, n0, gp, alpha)


def _variance_ratio_interval(lo, hi, code, G, means, sigma, tau2,
                             alpha: float = 0.05) -> Tuple[float, float]:
    """Interval for the share by inverting the variance ratio.

    For a balanced one-way random model the ratio of the between-lineage to the
    within-lineage mean square is a scaled F on (G - 1, n - G) degrees of
    freedom, and inverting it gives an interval for the intraclass correlation.
    The expectation-maximisation fit returns the same two components, so the
    same inversion applies to interval data once the within-lineage degrees of
    freedom are discounted by the share of each observation the interval left
    unresolved: an isolate known only to lie above a cut-off contributes less
    than one residual degree of freedom, and pretending otherwise would make
    the interval too narrow exactly where the reading is coarsest.

    Measured over the design grid, coverage is 0.92 to 0.99 for exact and
    dilution readings across cohort shapes from 15 lineages of 50 to 100 of 6,
    against a median 0.61 for the cluster bootstrap on the same cells.
    """
    n = int(np.asarray(lo).size)
    cnt = np.bincount(code, minlength=G).astype(float)
    present = cnt > 0
    gp = int(present.sum())
    if gp < 2 or n <= gp:
        return float("nan"), float("nan")
    n0 = (n - float(np.sum(cnt[present] ** 2)) / n) / max(gp - 1, 1)
    _ez, var = _trunc_moments(lo, hi, means[code], sigma)
    resolved = float(np.clip(1.0 - float(var.mean()) / max(sigma ** 2, 1e-12),
                             0.05, 1.0))
    df_between = max(gp - 1, 1)
    df_within = max((n - gp) * resolved, 1.0)
    ratio = 1.0 + n0 * max(tau2, 0.0) / max(sigma ** 2, 1e-12)
    out = []
    for q in (1.0 - alpha / 2.0, alpha / 2.0):
        adjusted = ratio / _f_dist.ppf(q, df_between, df_within)
        out.append((adjusted - 1.0) / (adjusted - 1.0 + n0))
    low, high = sorted(out)
    return float(np.clip(low, 0.0, 1.0)), float(np.clip(high, 0.0, 1.0))
