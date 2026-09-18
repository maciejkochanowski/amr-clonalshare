"""Lineage-associated variation under an interval-censored Gaussian model.

Exact log2 values, MIC intervals and thresholded calls can be represented in
the same latent model. Their information and identifiability differ. The
latent variance fraction is not the observed binary-scale attribution score.

Endpoint censoring is observed only when a recorded operator or known panel
supports it. Inferred endpoints are an explicit working assumption.

The point estimate uses an empirical-Bayes moment iteration, not an exact EM,
ML or REML fit. The F-based interval uses an approximate information-fraction
adjustment. Neither shrinkage nor an eligibility threshold proves calibration.
Historical calibration applies only to its recorded source and scenarios.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np
from scipy.stats import f as _f_dist

from .attribution import _codes, _is_untyped

__all__ = ["PanelGeometry", "CensoredShare", "panel_geometry",
           "intervals_from_mic", "intervals_by_panel", "intervals_from_binary",
           "scale_is_identified", "censored_clonal_share", "marginal_loglik"]

_SQRT2PI = float(np.sqrt(2.0 * np.pi))
# Below this share of isolates on an end well, a point-valued MIC is a fair
# description of the reading. Above it, the point value is a censored one
# wearing a number, and the interval mode is required. Read off the coverage
# curve in the calibration.
END_WELL_LIMIT = 0.05
# Share of isolates that may lie in lineages wholly beyond the panel before
# the share stops being estimable. A lineage all of whose readings are
# one-sided has no identified latent mean, and only the shrinkage keeps the
# fit finite; what decides whether that matters is how much of the cohort is
# in that position. Swept over 21 censoring levels with 60 simulated cohorts
# each (30 lineages of 25, true share 0.5), the bias stays under 0.005 up to a
# share of 0.35, reaches 0.029 at 0.52 and 0.079 at 0.69.
CENSORED_GROUP_LIMIT = 0.5
# Non-wild-type prevalence outside which a single cut point no longer carries
# the latent share. Read off the design grid: at a prevalence near 0.08 the
# estimate reads 0.23 above a true zero and its interval covers 0.08 of the
# time. Inside the window the estimate is defined and its bias stays under
# 0.05 at every true share, but the interval is calibrated only from a
# prevalence of about 0.24 (coverage 0.89 to 0.95); one doubling from the
# median, near 0.16 or 0.84, it covers a true zero 0.58 of the time. The
# record says so in ``notes``. A dilution at the same prevalence is
# unaffected, which is the practical argument for reading the panel rather
# than the call.
SINGLE_CUT_PREVALENCE = (0.10, 0.90)


@dataclass(frozen=True)
class PanelGeometry:
    """What the panel can and cannot support, before any biology.

    ``n_wells`` counts the wells the panel is taken to have tested: every
    doubling from the lowest to the highest recorded value when the recorded
    values sit on a doubling lattice (``lattice == "doubling"``), the distinct
    recorded values otherwise (``lattice == "recorded"``). ``n_wells_recorded``
    counts the distinct recorded values before rounding variants of one
    dilution, such as 0.06 and 0.064, are folded onto the same well.
    """

    n_wells: int
    lowest: float
    highest: float
    share_on_lowest: float
    share_on_highest: float
    lattice_ratios: Tuple[float, ...]
    doubling: bool
    admissible_modes: Tuple[str, ...]
    reason: str
    n_wells_recorded: int = 0
    lattice: str = "recorded"

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
    #: Isolates set aside because they carry no lineage label, under the same
    #: rule as :func:`amr_clonalshare.attribution.clonal_share`: an untyped
    #: isolate is not a lineage, and a level made of every isolate whose
    #: typing failed would measure the typing process.
    n_dropped_untyped: int = 0
    estimator_method: str = "empirical_Bayes_moment_iteration"
    interval_status: str = "approximate_F_interval; scenario-specific calibration only"
    n_dropped_uninformative: int = 0
    fit_converged: bool = False
    fit_iterations: int = 0
    fit_max_scaled_change: float = float("nan")
    complete: bool = True

    def as_dict(self) -> Dict[str, Any]:
        return asdict(self)


LATTICE_TOLERANCE = 0.25


def _lattice(v: np.ndarray, wells=None):
    """The tested wells on the log2 scale and the index of each reading.

    With ``wells`` given, every reading must match a tested concentration.
    Without them, readings that all sit within a quarter of a doubling of an
    integer power of two are taken to come from a doubling panel, whose
    tested wells are every doubling from the lowest to the highest recorded
    value: a dilution no isolate landed on was still tested, and 0.06 and
    0.064 are one well recorded two ways. Readings off any such lattice are
    taken as their own wells, which is the reading of a panel whose steps are
    not known.
    """
    lv = np.log2(v)
    if wells is not None:
        w = np.asarray(wells, dtype=float)
        if (w.ndim != 1 or w.size == 0 or not np.isfinite(w).all()
                or np.any(w <= 0) or np.any(np.diff(w) <= 0)):
            raise ValueError("wells must be strictly increasing positive finite concentrations")
        matches = np.isclose(v[:, None], w[None, :], rtol=1e-10, atol=0.0)
        if not matches.any(axis=1).all():
            raise ValueError("MIC values are outside configured wells; no panel snapping is allowed")
        lattice = np.log2(w)
        idx = matches.argmax(axis=1)
        return lattice, idx, "given"
    k = np.round(lv)
    if v.size and np.all(np.abs(lv - k) <= LATTICE_TOLERANCE):
        lo_k, hi_k = int(k.min()), int(k.max())
        lattice = np.arange(lo_k, hi_k + 1, dtype=float)
        return lattice, (k - lo_k).astype(int), "doubling"
    lattice = np.unique(lv)
    return lattice, np.searchsorted(lattice, lv), "recorded"


def panel_geometry(values: Sequence[float],
                   cutoff: Optional[float] = None, wells=None) -> PanelGeometry:
    """Describe the dilution panel and say which input modes it supports.

    The lattice ratios between adjacent recorded values are reported as they
    are: ratios near 1.03 to 1.05 are rounding variants of one dilution and a
    ratio of 8 is a pair of dilutions no isolate landed on. Neither is a
    property of the panel, and :func:`_lattice` decides what the tested wells
    are; this function reports both readings so that a reader can see what
    was recorded and what was inferred from it.
    """
    v = np.asarray([x for x in np.asarray(values, dtype=float)
                    if np.isfinite(x) and x > 0])
    if v.size == 0:
        return PanelGeometry(0, float("nan"), float("nan"), float("nan"),
                             float("nan"), (), False, (), "no usable readings",
                             0, "recorded")
    recorded = np.unique(v)
    ratios = tuple(round(float(r), 3)
                   for r in np.unique(recorded[1:] / recorded[:-1]))
    lattice, idx, kind = _lattice(v, wells)
    doubling = kind == "doubling" or (
        kind == "given" and lattice.size > 1
        and bool(np.allclose(np.diff(lattice), 1.0)))
    lo_share = float((idx == 0).mean())
    hi_share = float((idx == lattice.size - 1).mean())

    modes = ["binary"]
    reason = []
    if lattice.size >= 3:
        modes.append("interval")
    else:
        reason.append("fewer than three wells")
    if max(lo_share, hi_share) <= END_WELL_LIMIT:
        modes.append("point")
    else:
        reason.append(
            f"{max(lo_share, hi_share):.2f} of readings sit on an end well, so a "
            "point value is a censored reading and the interval mode is required")
    lowest, highest = float(2.0 ** lattice[0]), float(2.0 ** lattice[-1])
    if cutoff is not None and np.isfinite(cutoff):
        if not (lowest < cutoff <= highest):
            modes = [m for m in modes if m != "binary"]
            reason.append("the cut-off lies outside the tested range, so every "
                          "call is forced by the panel")
    return PanelGeometry(int(lattice.size), lowest, highest,
                         lo_share, hi_share, ratios, doubling,
                         tuple(modes), "; ".join(reason) or "panel supports all modes",
                         int(recorded.size), kind)


def intervals_from_binary(y, cutoff_log2: float = 0.0):
    """A dichotomised call as the interval it actually is.

    Wild type means the latent value did not exceed the cut-off; non-wild type
    means it did. Nothing else is known. This is the coarsest member of the
    family and exists so the same estimator can consume a cohort that reports
    only an interpretation.
    """
    y = np.asarray(y, dtype=float).ravel()
    if not np.isfinite(cutoff_log2):
        raise ValueError("cutoff_log2 must be finite")
    missing = np.isnan(y)
    if not np.isin(y[~missing], (0.0, 1.0)).all():
        raise ValueError("values must be binary numeric 0/1 or NaN")
    lo = np.where(y == 1.0, cutoff_log2, -np.inf)
    hi = np.where(y == 1.0, np.inf, cutoff_log2)
    lo[missing] = np.nan
    hi[missing] = np.nan
    return lo, hi


def _normalise_operators(operators):
    """Recorded signs or an empty cell; an unknown sign is not an assumption."""
    op = np.asarray(operators, dtype=object).ravel()
    cleaned = np.array(["" if v is None or str(v).strip().lower() in
                        ("", "nan", "none", "na", "<na>") else str(v).strip()
                        for v in op], dtype=str)
    invalid = sorted(set(cleaned) - {"", "=", "<", "<=", ">", ">="})
    if invalid:
        raise ValueError(f"unknown MIC operator(s): {invalid}; use =, <, <=, >, >= or an empty cell")
    return cleaned


def intervals_from_mic(values, *, operators=None, treat_end_wells_as_censored=True,
                       wells=None):
    """Intervals on the log2 scale from recorded MIC values.

    A reading at a tested dilution means the true value lies above the
    previous tested dilution and at or below this one. A reading on the lowest
    well is left-censored and one on the highest is right-censored, unless
    ``treat_end_wells_as_censored`` is turned off, which is the sensitivity
    arm for the coarsening assumption.

    ``wells`` are the concentrations the panel tested, when the laboratory
    recorded them; every reading is then placed on the nearest tested well.
    Without them the tested wells are inferred by :func:`_lattice`: every
    doubling from the lowest to the highest recorded value when the readings
    sit on a doubling lattice, so that a dilution no isolate landed on still
    bounds its neighbours and a rounding variant of one dilution is not read
    as a well of its own, and the recorded values themselves otherwise.

    ``operators`` may carry the recorded censoring signs; where present they
    win over the end-well heuristic, because a recorded operator is an
    observation and the heuristic is an assumption.
    """
    v = np.asarray(values, dtype=float).ravel()
    ok = np.isfinite(v) & (v > 0)
    lo = np.full(v.shape, -np.inf)
    hi = np.full(v.shape, np.inf)
    if ok.any():
        lattice, idx, _kind = _lattice(v[ok], wells)
        # The lowest well: left-censored under the assumption, or one
        # doubling wide, the width of every other well, in the sensitivity
        # arm that reads the end wells as exact; the highest well likewise.
        bottom_lo = -np.inf if treat_end_wells_as_censored else lattice[0] - 1.0
        prev = np.where(idx > 0, lattice[np.maximum(idx - 1, 0)], bottom_lo)
        lo[ok] = prev
        hi[ok] = lattice[idx]
        if treat_end_wells_as_censored:
            top = np.zeros(v.shape, dtype=bool)
            top[ok] = idx == lattice.size - 1
            hi[top] = np.inf
    if operators is not None:
        op = _normalise_operators(operators)
        if op.size != v.size:
            raise ValueError(
                f"operators has {op.size} entries and values has {v.size}; "
                f"they must match")
        gt = np.isin(op, (">", ">="))
        le = np.isin(op, ("<", "<="))
        exact = (op == "=") & ok
        if exact.any():
            # An exact MIC call is a finite dilution interval, including at
            # an end well. Only missing signs invoke the end-well heuristic.
            exact_lo = np.where(idx > 0, lattice[np.maximum(idx - 1, 0)], lattice[0] - 1.0)
            full_lo = np.full(v.shape, np.nan)
            full_hi = np.full(v.shape, np.nan)
            full_lo[ok], full_hi[ok] = exact_lo, lattice[idx]
            lo[exact], hi[exact] = full_lo[exact], full_hi[exact]
        lo[gt] = np.log2(np.where(v[gt] > 0, v[gt], 1.0)) if gt.any() else lo[gt]
        hi[gt] = np.inf
        lo[le] = -np.inf
        hi[le] = np.log2(np.where(v[le] > 0, v[le], 1.0)) if le.any() else hi[le]
    lo[~ok] = np.nan
    hi[~ok] = np.nan
    return lo, hi


def _normal_interval_logmass(a, b):
    """Stable log standard-normal interval mass; see _normal_numerics."""
    from ._normal_numerics import normal_interval_logmass
    return normal_interval_logmass(a, b)


def _trunc_moments(a, b, m, s):
    """Stable truncated-normal moments, retaining every positive interval width."""
    from ._normal_numerics import normal_truncated_moments
    return normal_truncated_moments(a, b, m, s)


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


def _fit_em(lo, hi, code, G, *, iters=1000, tol=1e-10, fixed_scale=None, diagnostics=None):
    """Approximate empirical-Bayes moment iteration for the interval model.

    The legacy function name is retained for compatibility. These updates
    are not an exact EM or REML maximisation of the marginal likelihood.
    Returns group means, residual scale and observed-group mean variance.
    """
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
    converged = False
    max_scaled_change = float("inf")
    iterations = 0
    for _ in range(int(iters)):
        iterations += 1
        previous_tau2, previous_grand = tau2, grand
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
        max_scaled_change = float(max(np.max(np.abs(new_m-m)/(1+np.abs(m))),
            abs(new_s-s)/(1+abs(s)), abs(tau2-previous_tau2)/(1+abs(previous_tau2)),
            abs(grand-previous_grand)/(1+abs(previous_grand))))
        if max_scaled_change < tol:
            m, s = new_m, new_s
            converged = True
            break
        m, s = new_m, new_s
    if diagnostics is not None:
        diagnostics.update(converged=converged, iterations=iterations, max_scaled_change=max_scaled_change)
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


def censored_clonal_share(lo, hi, lineage) -> CensoredShare:
    """Estimate a latent concentration variance ratio with a moment iteration.

    The empirical-Bayes updates are not exact EM, ML or REML. The F limits
    use an approximate information adjustment; nominal coverage is not
    guaranteed. Observed binary predictive association is a different
    target."""
    lo = np.asarray(lo, dtype=float).ravel()
    hi = np.asarray(hi, dtype=float).ravel()
    if lo.size != hi.size:
        raise ValueError(
            f"lo has {lo.size} entries and hi has {hi.size}")
    labels = np.asarray(lineage, dtype=object).ravel()
    if labels.shape[0] != lo.size:
        raise ValueError(
            f"lineage has {labels.shape[0]} entries for {lo.size} intervals")
    if np.any(lo > hi):
        bad = int(np.sum(lo > hi))
        raise ValueError(
            f"{bad} interval(s) have lo > hi; an interval is [lo, hi] on the "
            "log2 concentration scale, with lo == hi for an exact reading")
    if np.isposinf(lo).any() or np.isneginf(hi).any():
        raise ValueError("interval endpoints cannot be infinite exact observations")
    typed = np.array([not _is_untyped(v) for v in labels], dtype=bool)
    n_dropped_untyped = int((~typed).sum())
    uninformative = np.isneginf(lo) & np.isposinf(hi)
    n_dropped_uninformative = int(uninformative.sum())
    keep = typed & ~(np.isnan(lo) | np.isnan(hi)) & ~uninformative
    lo, hi = lo[keep], hi[keep]
    # Coded after filtering so that the codes stay dense; a gap would make
    # the group counts below index a lineage that is no longer present.
    code = _codes(labels[keep])
    n = lo.size
    G = int(code.max()) + 1 if n else 0
    censored = float(np.mean(~np.isfinite(lo) | ~np.isfinite(hi))) if n else float("nan")
    fully, fully_share = (_fully_censored_groups(lo, hi, code, G) if n
                          else (0, 0.0))

    if n < 8 or G < 2 or n <= G:
        nan = float("nan")
        return CensoredShare(nan, nan, nan, n, G, nan, nan,
                             censored, fully, fully_share, False,
                             f"{n} readable isolate(s) over {G} lineage(s), "
                             f"against the 8 isolates and 2 lineages the fit "
                             f"needs before a between-lineage component is "
                             f"identified; within-lineage replication is required",
                             n_dropped_untyped=n_dropped_untyped,
                             n_dropped_uninformative=n_dropped_uninformative)

    if np.all(lo == hi) and float(np.ptp(lo)) == 0.0:
        nan = float("nan")
        return CensoredShare(nan, nan, nan, n, G, nan, nan,
                             censored, fully, fully_share, False,
                             "every reading is the same value; there is no "
                             "variance to split",
                             n_dropped_untyped=n_dropped_untyped,
                             n_dropped_uninformative=n_dropped_uninformative)
    identified = scale_is_identified(lo, hi)
    fixed = None if identified else 1.0

    # The variance-component ratio of the fitted model. A held-out residual
    # ratio, which is what the binary estimator uses, is the wrong score here:
    # for a one-sided interval the conditional variance of the held-out reading
    # does not fall as the group mean improves, so the ratio is insensitive to
    # the very thing being measured.
    fit_diagnostics: Dict[str, Any] = {}
    m, s, tau2, _g = _fit_em(lo, hi, code, G, fixed_scale=fixed, diagnostics=fit_diagnostics)
    denom = tau2 + s * s
    point = float("nan") if denom <= 0 else float(tau2 / denom)
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
    converged = bool(fit_diagnostics["converged"])
    estimable = bool(np.isfinite(point) and not over_limit and not tail and converged)
    notes = []
    if not converged:
        notes.append("moment iteration limit reached; point retained for diagnostics only")
    if tail:
        notes.append(
            f"a single cut point with {above:.0%} of isolates above it lies "
            f"outside the historical reporting window "
            f"{SINGLE_CUT_PREVALENCE[0]:.0%} to {SINGLE_CUT_PREVALENCE[1]:.0%}, "
            "where the share is not estimable from a call; the recorded "
            "concentration is not affected")
    if over_limit:
        notes.append(
            f"{fully_share:.0%} of isolates lie in lineages wholly beyond the "
            f"panel, above the reporting threshold of "
            f"{CENSORED_GROUP_LIMIT:.0%}; their latent means are not "
            "identified and the share is not estimable")
    elif fully:
        notes.append(
            f"{fully} of {G} lineages lie wholly beyond the panel, holding "
            f"{fully_share:.0%} of isolates, which the historical calibration placed inside "
            "the range where the shrinkage absorbs them")
    if not identified:
        notes.append(
            "scale not identified from a single cut point; reported on the "
            "liability scale with residual standard deviation fixed at one, "
            "and the interval is a reading of the call rather than a "
            "calibrated interval: on the design grid it covers a true share "
            "of zero 0.58 of the time at a prevalence one doubling from the "
            "median")
    return CensoredShare(point, f_lo, f_hi, n, G, float(s),
                         float(np.sqrt(max(tau2, 0.0))), censored, fully,
                         float(fully_share), estimable,
                         "; ".join(notes) if notes else "ok",
                         realised_low, realised_high, n_dropped_untyped,
                         n_dropped_uninformative=n_dropped_uninformative,
                         fit_converged=converged, fit_iterations=int(fit_diagnostics["iterations"]),
                         fit_max_scaled_change=float(fit_diagnostics["max_scaled_change"]),
                         complete=converged)


def intervals_by_panel(values, panel, *, operators=None,
                       treat_end_wells_as_censored=True, wells=None):
    """:func:`intervals_from_mic` applied within each panel.

    ``panel`` names the panel every reading was made on. The tested wells,
    and with them the end wells, are inferred from the readings of one panel
    at a time, so that the lowest well of one laboratory is not taken for an
    interior dilution because another laboratory tested lower. Recorded
    ``wells`` apply to every panel alike.
    """
    v = np.asarray(values, dtype=float).ravel()
    labels = np.asarray([str(x) for x in np.asarray(panel, dtype=object).ravel()])
    if labels.size != v.size:
        raise ValueError("panel must name the panel of every reading")
    ops = None if operators is None else np.asarray(operators, dtype=object).ravel()
    lo, hi = np.empty(v.size), np.empty(v.size)
    for name in np.unique(labels):
        m = labels == name
        lo[m], hi[m] = intervals_from_mic(
            v[m], operators=None if ops is None else ops[m],
            treat_end_wells_as_censored=treat_end_wells_as_censored, wells=wells)
    return lo, hi


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
    """Gaussian random-intercept log likelihood for exact or interval data.

    Exact Gaussian observations use the analytic rank-one covariance formula.
    Interval observations use adaptive Gauss-Hermite integration with node
    refinement and a checked adaptive-integration fallback. Numerical agreement
    validates likelihood evaluation, not point-estimator or interval coverage.
    """
    a = np.asarray(lo, dtype=float).ravel()
    b = np.asarray(hi, dtype=float).ravel()
    code = np.asarray(code)
    if (not np.isfinite(grand) or not np.isfinite(tau2) or tau2 < 0
            or not np.isfinite(sigma) or sigma <= 0):
        raise ValueError("grand and tau2 must be finite, tau2 >= 0 and sigma > 0")
    if (a.shape != b.shape or code.ndim != 1 or code.size != a.size
            or not np.issubdtype(code.dtype, np.integer)
            or isinstance(G, bool) or not isinstance(G, (int, np.integer)) or G < 1
            or np.any(code < 0) or np.any(code >= G)):
        raise ValueError("intervals and integer group codes must have matching lengths and 0 <= code < G")
    if (np.isnan(a).any() or np.isnan(b).any() or np.any(a > b)
            or np.isposinf(a).any() or np.isneginf(b).any()):
        raise ValueError("interval endpoints must be ordered and not NaN or infinite exact readings")
    code = code.astype(np.int64, copy=False)
    s = float(sigma)
    if np.all(np.isfinite(a) & (a == b)):
        # A Gaussian random intercept adds a rank-one covariance. Evaluate
        # its exact likelihood rather than integrating a narrow peak on fixed nodes.
        counts = np.bincount(code, minlength=G).astype(float)
        present = counts > 0
        counts = counts[present]
        means = np.bincount(code, weights=b, minlength=G)[present] / counts
        group_means = np.zeros(G)
        group_means[present] = means
        within = np.bincount(code, weights=(b-group_means[code])**2, minlength=G)[present]
        variance = s*s
        between_variance = variance + counts*max(tau2, 0.0)
        logdet = (counts-1)*np.log(variance) + np.log(between_variance)
        quadratic = within/variance + counts*(means-grand)**2/between_variance
        return float(-0.5*np.sum(counts*np.log(2*np.pi)+logdet+quadratic))
    if tau2 <= 0.0:
        cell = _normal_interval_logmass((a-grand)/s, (b-grand)/s)
        exact = np.isfinite(a) & (a == b)
        density = -0.5*((b-grand)/s)**2 - np.log(s) - np.log(_SQRT2PI)
        return float(np.where(exact, density, cell).sum())
    from ._mic_quadrature import adaptive_loglik
    return adaptive_loglik(a, b, np.asarray(code, dtype=int), G, grand,
                           tau2, s, _normal_interval_logmass, _trunc_moments)


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
    # The within-lineage degrees of freedom are scaled by the share of the
    # residual variance the readings resolve, one for exact values and less
    # for wide intervals. The floor of 0.05 keeps a cohort read almost
    # entirely from end wells from reaching zero degrees of freedom; it is a
    # guard, not a calibrated constant.
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
    """Approximate F-based limits from fitted variance components.
    
    A balanced Gaussian mean-square identity motivates this calculation.
    Censored moment estimates and an information-fraction adjustment do
    not yield an exact pivot. Scenario-specific calibration is required."""
    parts = _ratio_and_degrees(lo, hi, code, G, means, sigma, tau2)
    if parts is None:
        return float("nan"), float("nan")
    ratio, df_between, df_within, n0, _gp, _ = parts
    out = []
    for q in (1.0 - alpha / 2.0, alpha / 2.0):
        adjusted = ratio / _f_dist.ppf(q, df_between, df_within)
        out.append((adjusted - 1.0) / (adjusted - 1.0 + n0))
    low, high = sorted(out)
    return float(np.clip(low, 0.0, 1.0)), float(np.clip(high, 0.0, 1.0))
