"""Descriptive lineage composition and within-lineage rate contrasts.

For two analysed collections, the Kitagawa midpoint identity splits the
prevalence difference exactly into composition and within-lineage terms.
For a lineage observed in only one collection, the absent rate is set equal
to its observed rate. This declared convention allocates its entire change
to composition; another convention can change the full decomposition.

The shared-support statistic is the unweighted mean of collection-specific
overlap proportions, which are also returned separately. Eligibility is an
empirical reporting rule, not proof of identifiability or absence of selection.
Isolate resampling assumes the corresponding within-collection sampling model;
additional farm, patient, laboratory or temporal dependence is not accounted for.

Multiplicity correction uses BY within each component family. Its theoretical
FDR guarantee requires valid input p-values. Bootstrap-tail values are
approximate, so arbitrary-dependence correction alone does not prove coverage
or error control. Neither component identifies a biological mechanism.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import beta as _beta

from ._seeding import estimator_rng
from .attribution import _is_untyped

from .stats import benjamini_hochberg, effective_dimension, fisher_exact_p

__all__ = ["decompose_prevalence_difference", "decompose_panel"]

#: Default lower bound on the isolate share held by lineages common to both
#: collections. Below it both components are reported but marked not
#: estimable.
#:
#: A lineage that one collection holds and the other sampled no isolate of is
#: charged wholly to composition, so where many lineages are seen in one
#: collection only, both components are biased, in opposite directions. The
#: value is read from the grid of ``benchmarks/decomposition_calibration.py``
#: by ``benchmarks/summarize_model_checks.py gate``: the lowest shared
#: support, in steps of 0.05, at and above which no cell covers either
#: component below 0.89. The same summary records four other bootstrap
#: intervals (basic, bias-corrected, bias-shifted and studentised) computed
#: from the same draws.
DEFAULT_MIN_SHARED_SUPPORT = 0.9

#: Cap on the number of resampled isolate labels held in memory at once. The
#: bootstrap is drawn in chunks of this many cells so that a large replicate
#: count does not scale memory with it.
_CHUNK_CELLS = 4_000_000


def _clean_with_mask(y: Sequence,
                     lineage: Sequence) -> Tuple[np.ndarray, np.ndarray,
                                                 np.ndarray]:
    """The same cleaning, returning the full arrays and the keep mask.

    The dropped isolates are the evidence for whether dropping them was safe,
    so they have to survive the step that removes them.
    """
    arr = np.asarray(list(y), dtype=float)
    lin = pd.Series(list(lineage), dtype="object")
    if arr.size != lin.size:
        raise ValueError("y and lineage must have the same length")
    ok = np.asarray(lin.notna().to_numpy(), dtype=bool, copy=True)
    lin_s = np.asarray(lin.fillna("__missing__").astype(str).to_numpy(),
                       dtype=object)
    ok &= np.array([not _is_untyped(x) for x in lin_s], dtype=bool)
    ok &= np.isfinite(arr)
    return arr, lin_s, ok


def _lineage_present_mask(lineage: Sequence) -> np.ndarray:
    """Which isolates carry a usable lineage label, whatever the trait is.

    ``_clean_with_mask`` returns "label present **and** trait finite", which is
    the right mask to run the analysis on and the wrong one for the numbers
    reported beside it: a NaN in the trait was published as a missing lineage
    label under ``n_dropped_missing_lineage`` and as a fall in
    ``label_coverage_a``. The coverage
    and the dropped-label counts read this mask instead; the trait-missing
    count is reported separately.
    """
    lin = pd.Series(list(lineage), dtype="object")
    ok = np.asarray(lin.notna().to_numpy(), dtype=bool, copy=True)
    lin_s = np.asarray(lin.fillna("__missing__").astype(str).to_numpy(),
                       dtype=object)
    ok &= np.array([not _is_untyped(x) for x in lin_s], dtype=bool)
    return ok


def _shared_index(lin_a: np.ndarray,
                  lin_b: np.ndarray) -> Tuple[np.ndarray, np.ndarray, list]:
    """Integer codes for both collections over one sorted lineage index.

    Coding once and counting with ``bincount`` replaces a pandas ``groupby``
    per bootstrap replicate. The arithmetic is unchanged - a lineage absent
    from a replicate has zero share and an undefined rate in both collections,
    so it contributes nothing to either component - and the resampling becomes
    fast enough for the calibration study that justifies the intervals.
    """
    lineages = sorted(set(lin_a.tolist()) | set(lin_b.tolist()))
    index = {value: position for position, value in enumerate(lineages)}
    codes_a = np.fromiter((index[v] for v in lin_a), dtype=np.intp,
                          count=lin_a.size)
    codes_b = np.fromiter((index[v] for v in lin_b), dtype=np.intp,
                          count=lin_b.size)
    return codes_a, codes_b, lineages


def _shares_and_rates(codes: np.ndarray, y: np.ndarray,
                      n_lineages: int) -> Tuple[np.ndarray, np.ndarray]:
    """Lineage shares ``w`` and within-lineage rates ``p`` over a fixed index.

    ``p`` is NaN for a lineage absent from this collection; the caller applies
    the declared convention.
    """
    counts = np.bincount(codes, minlength=n_lineages).astype(float)
    sums = np.bincount(codes, weights=y, minlength=n_lineages)
    total = counts.sum()
    w = counts / total if total > 0 else counts
    with np.errstate(invalid="ignore", divide="ignore"):
        p = np.where(counts > 0, sums / np.where(counts > 0, counts, 1.0),
                     np.nan)
    return w, p


def _shares_and_rates_batch(codes: np.ndarray, y: np.ndarray,
                            n_lineages: int) -> Tuple[np.ndarray, np.ndarray]:
    """``_shares_and_rates`` for a stack of resampled collections."""
    n_rep, _ = codes.shape
    offset = (np.arange(n_rep, dtype=np.intp)[:, None] * n_lineages
              + codes).ravel()
    size = n_rep * n_lineages
    counts = np.bincount(offset, minlength=size).reshape(n_rep,
                                                         n_lineages).astype(float)
    sums = np.bincount(offset, weights=y.ravel(),
                       minlength=size).reshape(n_rep, n_lineages)
    total = counts.sum(axis=1, keepdims=True)
    w = np.divide(counts, total, out=np.zeros_like(counts), where=total > 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = np.where(counts > 0, sums / np.where(counts > 0, counts, 1.0),
                     np.nan)
    return w, p


def _kitagawa(wa: np.ndarray, pa: np.ndarray, wb: np.ndarray,
              pb: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Composition and within-lineage components, absent rates filled.

    Works on a single pair of collections and on a stack of them; the sum runs
    over the last axis either way.
    """
    pa_f = np.nan_to_num(np.where(np.isnan(pa), pb, pa))
    pb_f = np.nan_to_num(np.where(np.isnan(pb), pa, pb))
    composition = ((wa - wb) * (pa_f + pb_f) / 2.0).sum(axis=-1)
    within = ((wa + wb) / 2.0 * (pa_f - pb_f)).sum(axis=-1)
    return composition, within


def _bootstrap_pvalue(draws: np.ndarray) -> float:
    """Two-sided bootstrap p-value for a component being zero.

    The percentile interval and this p-value are close to the same statement
    read two ways: the p-value is about the smallest ``alpha`` at which the
    ``1 - alpha`` interval still excludes zero. The ``+1`` follows Phipson and
    Smyth for the same reason it does in a permutation test - a resampled
    tail count of zero is not evidence that the tail is empty - and it makes
    the p-value a little larger than the interval alone would say, so at the
    boundary the two can disagree by one resample.
    """
    finite = draws[np.isfinite(draws)]
    n_rep = finite.size
    if n_rep == 0:
        return float("nan")
    below = float((np.sum(finite <= 0.0) + 1) / (n_rep + 1))
    above = float((np.sum(finite >= 0.0) + 1) / (n_rep + 1))
    return float(min(1.0, 2.0 * min(below, above)))


def _clopper_pearson(successes: int, trials: int,
                     confidence: float = 0.95) -> Tuple[float, float]:
    """Exact binomial interval for a Monte Carlo tail probability."""
    alpha = 1.0 - confidence
    low = (0.0 if successes == 0 else
           float(_beta.ppf(alpha / 2, successes, trials - successes + 1)))
    high = (1.0 if successes == trials else
            float(_beta.ppf(1 - alpha / 2, successes + 1, trials - successes)))
    return low, high


def _label_availability(y_a: np.ndarray, labelled_a: np.ndarray,
                        y_b: np.ndarray, labelled_b: np.ndarray,
                        alpha: float) -> Dict[str, Any]:
    """Diagnose selection into the observed, lineage-labelled subset.

    Fisher tests report whether the observed data detect an association; a
    nonsignificant result and equal label coverage do not establish absence of
    selection. ``selection_status`` therefore distinguishes complete data,
    detected association, an unresolved diagnostic, and a trait the binary
    diagnostic cannot assess. It makes no MCAR or MAR claim.

    ``collection_generalization_supported`` is true only for complete outcome
    and label data, and concerns only the supplied finite collections and
    never a broader population.
    """
    binary = bool(np.all(np.isin(y_a[np.isfinite(y_a)], (0.0, 1.0)))
                  and np.all(np.isin(y_b[np.isfinite(y_b)], (0.0, 1.0))))
    coverage_a = float(labelled_a.mean()) if labelled_a.size else None
    coverage_b = float(labelled_b.mean()) if labelled_b.size else None
    coverage_p = fisher_exact_p(
        int(labelled_a.sum()), int((~labelled_a).sum()),
        int(labelled_b.sum()), int((~labelled_b).sum()))

    gaps: Dict[str, Any] = {}
    trait_p: Dict[str, Optional[float]] = {}
    for side, y, labelled_mask in (("a", y_a, labelled_a),
                                    ("b", y_b, labelled_b)):
        unlabelled = ~labelled_mask & np.isfinite(y)
        labelled = labelled_mask & np.isfinite(y)
        if not binary or unlabelled.sum() == 0 or labelled.sum() == 0:
            gaps[side] = None
            trait_p[side] = None
            continue
        gaps[side] = float(y[unlabelled].mean() - y[labelled].mean())
        trait_p[side] = fisher_exact_p(
            int(y[unlabelled].sum()), int((1 - y[unlabelled]).sum()),
            int(y[labelled].sum()), int((1 - y[labelled]).sum()))

    informative = any(p is not None and p < alpha for p in trait_p.values())
    differential = bool(coverage_p < alpha)
    any_missing_label = bool((~labelled_a).any() or (~labelled_b).any())
    any_missing_trait = bool((~np.isfinite(y_a)).any()
                             or (~np.isfinite(y_b)).any())
    nonempty = bool(y_a.size and y_b.size)
    complete = nonempty and not any_missing_label and not any_missing_trait
    if not nonempty:
        selection_status = "unassessable"
    elif complete:
        selection_status = "complete"
    elif binary and informative:
        selection_status = "detected association"
    elif not binary and any_missing_label:
        selection_status = "unassessable"
    else:
        selection_status = "undetermined"
    note = (
        "the decomposition describes isolates with an observed trait and a "
        "lineage label. Selection tests are diagnostics: significance detects "
        "an association, while equal coverage or nonsignificance cannot "
        "establish representativeness. Complete outcome and label data alone "
        "support generalization to the supplied finite collections; no field "
        "supports extrapolation beyond those collections")
    return {
        "label_coverage_a": coverage_a,
        "label_coverage_b": coverage_b,
        "label_coverage_differs_p": float(coverage_p),
        "label_coverage_differs": differential,
        "trait_gap_unlabelled_minus_labelled_a": gaps["a"],
        "trait_gap_unlabelled_minus_labelled_b": gaps["b"],
        "trait_gap_p_a": trait_p["a"],
        "trait_gap_p_b": trait_p["b"],
        "missingness_informative": bool(informative),
        "selection_status": selection_status,
        "collection_generalization_supported": bool(complete),
        "assessable": bool(binary),
        "note": note,
    }


def decompose_prevalence_difference(
    y_a: Sequence, lineage_a: Sequence,
    y_b: Sequence, lineage_b: Sequence,
    *,
    n_boot: int = 2000,
    rng: Optional[np.random.Generator] = None,
    min_shared_support: float = DEFAULT_MIN_SHARED_SUPPORT,
    label_alpha: float = 0.05,
) -> Dict[str, Any]:
    """Split a prevalence difference into composition and within-lineage parts.

    Parameters
    ----------
    y_a, y_b : per-isolate trait indicator in each collection. Designed for
        0/1 non-wild-type calls; any finite numeric works and the components
        are then differences of means.
    lineage_a, lineage_b : lineage label per isolate. Missing labels are
        dropped and counted.
    n_boot : bootstrap replicates. Isolates are resampled within each
        collection, which is the sampling that produced the difference.
    min_shared_support : both components are marked not estimable when the
        isolate share held by lineages common to both collections falls below
        this. Set to 0 to report them unconditionally.

    Returns
    -------
    dict carrying the two components with percentile intervals, two-sided
    bootstrap p-values, the exactness residual of the identity, the lineage
    coverage of each collection, and ``turnover_share`` - the fraction of the
    moving mass held by lineages seen in only one collection.
    ``within_lineage`` is identified on the shared lineages alone;
    ``shared_support_isolate_share`` is the coverage it rests on and
    ``within_lineage_estimable`` states whether that coverage clears the gate.
    The arithmetic describes only isolates with both an observed trait and a
    lineage label. ``collection_generalization_supported`` is separately true
    only when neither supplied finite collection has a missing trait or label;
    it does not support extrapolation beyond those collections.
    """
    if rng is None:
        rng = np.random.default_rng(0)
    all_ya, all_lina, keep_a = _clean_with_mask(y_a, lineage_a)
    all_yb, all_linb, keep_b = _clean_with_mask(y_b, lineage_b)
    label_a, label_b = (_lineage_present_mask(lineage_a),
                        _lineage_present_mask(lineage_b))
    ya, lina, drop_a = all_ya[keep_a], all_lina[keep_a], int((~label_a).sum())
    yb, linb, drop_b = all_yb[keep_b], all_linb[keep_b], int((~label_b).sum())
    drop_trait = int((~np.isfinite(all_ya)).sum() + (~np.isfinite(all_yb)).sum())
    availability = _label_availability(all_ya, label_a, all_yb, label_b,
                                       label_alpha)
    finite_a, finite_b = np.isfinite(all_ya), np.isfinite(all_yb)
    analysis_scope = {
        "target": "observed trait among lineage-labelled isolates",
        "frame": "supplied finite collections",
        "n_total_a": int(all_ya.size),
        "n_total_b": int(all_yb.size),
        "n_trait_observed_a": int(finite_a.sum()),
        "n_trait_observed_b": int(finite_b.sum()),
        "n_lineage_labelled_a": int(label_a.sum()),
        "n_lineage_labelled_b": int(label_b.sum()),
        "n_analyzed_a": int(keep_a.sum()),
        "n_analyzed_b": int(keep_b.sum()),
        "n_missing_trait_a": int((~finite_a).sum()),
        "n_missing_trait_b": int((~finite_b).sum()),
        "n_missing_lineage_a": int((~label_a).sum()),
        "n_missing_lineage_b": int((~label_b).sum()),
    }
    observed_prevalence_a = (float(all_ya[finite_a].mean())
                             if finite_a.any() else None)
    observed_prevalence_b = (float(all_yb[finite_b].mean())
                             if finite_b.any() else None)
    observed_difference = (
        float(observed_prevalence_a - observed_prevalence_b)
        if observed_prevalence_a is not None and observed_prevalence_b is not None
        else None)
    if ya.size == 0 or yb.size == 0:
        return {"status": "skipped",
                "reason": "no isolate with a usable lineage label in one or "
                          "both collections",
                "analysis_scope": analysis_scope,
                "selection_status": availability["selection_status"],
                "collection_generalization_supported": bool(
                    availability["collection_generalization_supported"]),
                "lineage_label_availability": availability,
                "n_a": int(ya.size), "n_b": int(yb.size),
                "n_trait_observed_a": int(finite_a.sum()),
                "n_trait_observed_b": int(finite_b.sum()),
                "n_typed_a": int(ya.size), "n_typed_b": int(yb.size),
                "observed_prevalence_a": observed_prevalence_a,
                "observed_prevalence_b": observed_prevalence_b,
                "observed_difference": observed_difference,
                "typed_prevalence_a": float(ya.mean()) if ya.size else None,
                "typed_prevalence_b": float(yb.mean()) if yb.size else None,
                "typed_difference": None,
                "n_dropped_missing_lineage": int(drop_a + drop_b),
                "n_dropped_non_finite_trait": int(drop_trait)}

    codes_a, codes_b, lineages = _shared_index(lina, linb)
    n_lin = len(lineages)
    wa, pa = _shares_and_rates(codes_a, ya, n_lin)
    wb, pb = _shares_and_rates(codes_b, yb, n_lin)
    comp_arr, within_arr = _kitagawa(wa, pa, wb, pb)
    comp, within = float(comp_arr), float(within_arr)
    diff = float(ya.mean() - yb.mean())

    pa_f = np.nan_to_num(np.where(np.isnan(pa), pb, pa))
    pb_f = np.nan_to_num(np.where(np.isnan(pb), pa, pb))
    contrib = wa * pa_f - wb * pb_f
    only_a = np.isnan(pb) & ~np.isnan(pa)
    only_b = np.isnan(pa) & ~np.isnan(pb)
    unshared = only_a | only_b
    shared = ~unshared & (wa + wb > 0)
    mass = np.abs(contrib).sum()
    share = float(np.abs(contrib[unshared]).sum() / mass) if mass > 0 else None
    support = float((wa[shared].sum() + wb[shared].sum()) / 2.0)
    # The part of the difference carried by lineages seen in one collection
    # only. Under the convention below it sits wholly inside ``composition``,
    # which is therefore the sum of a shared-lineage term and this one; the
    # two are reported apart so that a composition finding can be read for
    # how much of it rests on lineages the other collection never held.
    nonshared = float(contrib[unshared].sum()) if unshared.any() else 0.0

    boot = _bootstrap_components(codes_a, ya, codes_b, yb, n_lin,
                                 int(n_boot), rng)
    if boot.shape[0]:
        ci = np.percentile(boot, [2.5, 97.5], axis=0)
        se = boot.std(axis=0, ddof=1) if boot.shape[0] > 1 else np.zeros(2)
        p_comp = _bootstrap_pvalue(boot[:, 0])
        p_within = _bootstrap_pvalue(boot[:, 1])
        # The floor is set by the draws that came back finite, which is what
        # :func:`_bootstrap_pvalue` divides by; a replicate that produced no
        # finite component bought no resolution and must not be counted as if
        # it had.
        n_finite = int(min(np.isfinite(boot[:, 0]).sum(),
                           np.isfinite(boot[:, 1]).sum()))
    else:
        ci = np.full((2, 2), np.nan)
        se = np.full(2, np.nan)
        p_comp = p_within = float("nan")
        n_finite = 0

    estimable = support >= float(min_shared_support)
    return {
        "status": "ok",
        "analysis_scope": analysis_scope,
        "selection_status": availability["selection_status"],
        "collection_generalization_supported": bool(
            availability["collection_generalization_supported"]),
        "lineage_label_availability": availability,
        "n_a": int(ya.size), "n_b": int(yb.size),
        "n_trait_observed_a": int(finite_a.sum()),
        "n_trait_observed_b": int(finite_b.sum()),
        "n_typed_a": int(ya.size), "n_typed_b": int(yb.size),
        "n_dropped_missing_lineage": int(drop_a + drop_b),
        "n_dropped_non_finite_trait": int(drop_trait),
        "n_lineages": int(n_lin),
        "n_lineages_shared": int(shared.sum()),
        "n_lineages_only_a": int(only_a.sum()),
        "n_lineages_only_b": int(only_b.sum()),
        "shared_support_isolate_share": support,
        "shared_support_a": float(wa[shared].sum()),
        "shared_support_b": float(wb[shared].sum()),
        "shared_support_definition": "unweighted mean of the two collection-specific proportions",
        "prevalence_a": float(ya.mean()),
        "prevalence_b": float(yb.mean()),
        "difference": diff,
        "typed_prevalence_a": float(ya.mean()),
        "typed_prevalence_b": float(yb.mean()),
        "typed_difference": diff,
        "observed_prevalence_a": observed_prevalence_a,
        "observed_prevalence_b": observed_prevalence_b,
        "observed_difference": observed_difference,
        "composition": comp,
        "composition_shared": float(comp - nonshared),
        "nonshared": nonshared,
        "within_lineage": within,
        "composition_ci95": [float(ci[0, 0]), float(ci[1, 0])],
        "within_lineage_ci95": [float(ci[0, 1]), float(ci[1, 1])],
        "composition_se": float(se[0]),
        "within_lineage_se": float(se[1]),
        "composition_p": p_comp,
        "within_lineage_p": p_within,
        "n_boot_finite": n_finite,
        "p_value_floor": float(2.0 / (n_finite + 1)) if n_finite else None,
        "identity_residual": float(comp + within - diff),
        "turnover_share": share,
        "within_lineage_estimable": bool(estimable),
        "composition_estimable": bool(estimable),
        "not_estimable_because": (
            [] if support >= float(min_shared_support)
            else ["shared lineage support below the threshold"]),
        "shared_support_threshold": float(min_shared_support),
        "shared_support_margin": float(support - float(min_shared_support)),
        "n_boot": int(n_boot),
        "convention": "a lineage absent from one collection is given that "
                      "collection's rate equal to the observed one; this makes "
                      "within_lineage a function of shared lineages only, and "
                      "under this declared convention, and allocates the whole "
                      "contribution of a lineage seen in one collection only "
                      "in composition, where it is reported apart as "
                      "nonshared, with composition_shared the remainder",
        "note": "the typed difference is described, not explained: it applies "
                "to isolates with observed traits and lineage labels. The "
                "observed collection difference is reported separately; "
                "neither component is a causal effect",
    }


def _bootstrap_components(codes_a: np.ndarray, ya: np.ndarray,
                          codes_b: np.ndarray, yb: np.ndarray,
                          n_lineages: int, n_boot: int,
                          rng: np.random.Generator) -> np.ndarray:
    """Percentile-bootstrap draws of the two components.

    Isolates are resampled with replacement inside each collection, which is
    the sampling that produced the observed difference. Replicates are drawn
    in chunks so that memory is set by the cohort rather than by ``n_boot``.
    """
    if n_boot <= 0:
        return np.empty((0, 2), dtype=float)
    na, nb = ya.size, yb.size
    per_replicate = max(na + nb, 1)
    chunk = max(1, min(n_boot, _CHUNK_CELLS // per_replicate))
    out = np.empty((n_boot, 2), dtype=float)
    done = 0
    while done < n_boot:
        size = min(chunk, n_boot - done)
        ia = rng.integers(0, na, (size, na))
        ib = rng.integers(0, nb, (size, nb))
        wa, pa = _shares_and_rates_batch(codes_a[ia], ya[ia], n_lineages)
        wb, pb = _shares_and_rates_batch(codes_b[ib], yb[ib], n_lineages)
        composition, within = _kitagawa(wa, pa, wb, pb)
        out[done:done + size, 0] = composition
        out[done:done + size, 1] = within
        done += size
    return out


# The stream the panel decomposition draws from, one generator per agent.
_DECOMPOSITION_STREAM = 0


def decompose_panel(
    y_a: "pd.DataFrame", lineage_a: Sequence,
    y_b: "pd.DataFrame", lineage_b: Sequence,
    *,
    n_boot: int = 2000,
    rng: Optional[np.random.Generator] = None,
    q: float = 0.05,
    min_shared_support: float = DEFAULT_MIN_SHARED_SUPPORT,
    label_alpha: float = 0.05,
) -> Dict[str, Any]:
    """Decompose every agent in a panel and control the family error rate.

    An antimicrobial panel is a family of tests and its members are not
    independent: cross-resistance makes a macrolide block behave as one trait,
    and a tetracycline pair can be the same column twice. Reporting thirteen
    nominal 95 % intervals therefore overstates both the number of findings and
    their number of chances.

    Two corrections are applied and both are reported. The false discovery rate
    is controlled **within each component family**, since a composition finding
    and a within-lineage finding answer different questions and are not
    exchangeable. The procedure is the Benjamini-Yekutieli step-up, which is
    valid whatever the dependence between the agents. Plain Benjamini-Hochberg
    is not used here because it needs the dependence to be positive, and the
    paragraph above gives the reason it cannot be assumed to be: agents sharing
    a target site rise together, while a resistance trade-off makes a pair move
    apart. The discovery count under the independence assumption is reported
    beside it, so the price of the weaker assumption is visible rather than
    silent. The effective number of independent
    agents is the participation ratio of the panel's correlation eigenspectrum
    (:func:`stats.effective_dimension`): a panel of
    thirteen agents carrying four distinct resistance phenotypes has an
    effective size near four, and that number belongs beside the count of
    discoveries.

    Parameters
    ----------
    y_a, y_b : isolate-by-agent indicator frames for the two collections. They
        must carry the same agent columns.
    lineage_a, lineage_b : lineage label per isolate in each frame.
    q : target false discovery rate within each component family.
    label_alpha : level of the two label-availability tests that decide
        whether a lineage-resolved reading describes the collections. Both are
        Fisher exact, so on a trait that is not binary the comparison cannot
        be made and an agent with any unlabelled isolate is refused rather
        than cleared.

    Returns
    -------
    dict with ``per_agent`` results, each carrying ``composition_q`` and
    ``within_lineage_q``, and a ``family`` block recording the panel size, its
    effective size, and the discovery counts before and after control.
    """
    if rng is None:
        rng = np.random.default_rng(0)
    agents = [str(c) for c in y_a.columns]
    if [str(c) for c in y_b.columns] != agents:
        raise ValueError("both collections must carry the same agent columns")

    # Each agent resamples from a fresh generator on the same stream, drawn
    # from the caller's generator once: the bootstrap of one antimicrobial
    # then does not move when another is added to the panel or the columns
    # arrive in another order.
    root = int(rng.integers(0, 2 ** 63 - 1))
    per_agent: Dict[str, dict] = {}
    for agent in agents:
        per_agent[agent] = decompose_prevalence_difference(
            y_a[agent].to_numpy(), lineage_a,
            y_b[agent].to_numpy(), lineage_b,
            n_boot=n_boot, rng=estimator_rng(root, _DECOMPOSITION_STREAM),
            min_shared_support=min_shared_support, label_alpha=label_alpha)

    usable = [a for a in agents if per_agent[a].get("status") == "ok"]
    family: Dict[str, Any] = {
        "n_agents": len(agents),
        "n_agents_decomposed": len(usable),
        "q": float(q),
        "discovery_scope": "observed trait among lineage-labelled isolates",
        "collection_generalization_rule": (
            "requires complete trait and lineage data for the supplied finite "
            "collections"),
        "method": "Benjamini-Yekutieli within each component family, "
                  "valid under arbitrary dependence",
    }
    for key, label in (("composition", "composition"),
                       ("within_lineage", "within_lineage")):
        pvals = np.array([per_agent[a][f"{key}_p"] for a in usable], dtype=float)
        if usable and np.all(np.isfinite(pvals)):
            adjusted, reject = benjamini_hochberg(pvals, q=q,
                                                  dependence="arbitrary")
            _, reject_pos = benjamini_hochberg(pvals, q=q,
                                               dependence="independent")
        else:
            adjusted = np.full(len(usable), np.nan)
            reject = np.zeros(len(usable), dtype=bool)
            reject_pos = np.zeros(len(usable), dtype=bool)
            if usable:
                family[f"{key}_note"] = (
                    "no false-discovery control: the bootstrap p-values are "
                    "not finite (no finite draw, or surveillance.n_boot is 0), "
                    "so no component was selected")
        # A within-lineage component the shared-support gate refused cannot be
        # a subset-scope discovery whatever its p-value. Label missingness does
        # not erase the analysed-subset result; its separate selection and
        # collection-generalization fields limit the interpretation.
        admitted = np.array([
            (key != "within_lineage"
             or bool(per_agent[a]["within_lineage_estimable"]))
            and (key != "composition"
                 or bool(per_agent[a]["composition_estimable"]))
            for a in usable], dtype=bool)
        reject = reject & admitted
        reject_pos = reject_pos & admitted
        for position, agent in enumerate(usable):
            per_agent[agent][f"{key}_q"] = float(adjusted[position])
            per_agent[agent][f"{key}_discovery"] = bool(reject[position])
        family[f"n_{label}_nominal"] = int(sum(
            1 for a, ok in zip(usable, admitted)
            if ok and per_agent[a][f"{key}_ci95"][0] * per_agent[a][f"{key}_ci95"][1] > 0))
        family[f"n_{label}_discoveries"] = int(reject.sum())
        family[f"n_{label}_discoveries_under_independence"] = int(reject_pos.sum())
        family[f"n_{label}_refused"] = int((~admitted).sum())

    # The BY step-up can reject several hypotheses jointly. Its attainable
    # minimum is not the rank-one threshold m*H_m*p_min. Apply the actual
    # adjustment to the vector of attainable per-agent floors.
    m = len(usable)
    realised = [int(per_agent[a].get("n_boot_finite") or 0) for a in usable]
    n_finite = min(realised) if realised else 0
    if m and n_finite:
        floors = np.array([min(1.0, 2.0/(n+1)) for n in realised])
        adjusted_floors, _ = benjamini_hochberg(floors, q=q, dependence="arbitrary")
        floor_q = float(np.min(adjusted_floors))
        family["smallest_attainable_q"] = floor_q
        family["rank_one_q_floor"] = float(min(1.0, m*sum(1.0/k for k in range(1,m+1))*floors.min()))
        if floor_q > float(q):
            family["warning"] = (
                f"with at least {n_finite} finite draw(s), the best attainable "
                f"BY-adjusted value over {m} agents is {floor_q:.3f}, above "
                f"the target {float(q):g}; increase surveillance.n_boot before "
                "interpreting an absence of discoveries")

    stacked = pd.concat([y_a[agents], y_b[agents]], axis=0).to_numpy(dtype=float)
    family["effective_independent_agents"] = float(effective_dimension(stacked))
    family["note"] = (
        "discoveries have the stated labelled-subset scope. The effective "
        "count is the participation ratio of the panel correlation "
        "eigenspectrum; a panel whose agents share a resistance mechanism is "
        "worth fewer independent tests than it has columns")
    return {"status": "ok", "per_agent": per_agent, "family": family}
