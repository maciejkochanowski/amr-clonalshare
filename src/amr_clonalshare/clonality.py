"""clonality.py - separating a change in lineage *mix* from a change in *rate*.

Surveillance asks one question of a resistance trend that the reported number
cannot answer. When non-wild-type prevalence for an agent differs between two
collections - two years, two countries, two hosts - the difference has two
mechanisms, and they call for opposite interventions:

* the **composition** of the population changed, because a lineage that already
  carried the trait became a larger share of what was sampled. The response is
  transmission control: biosecurity, movement restriction, cleaning.
* the **within-lineage rate** changed, because the same lineages became more
  often non-wild-type. The response is selection-pressure control:
  stewardship, dosing, agent choice.

The distinction is recognised in the genomic-surveillance literature - reports
of "widespread acquisition" alongside "clonal expansion of resistant lineages"
are describing exactly these two components - but it is read qualitatively off
a phylogeny. There is no estimator that splits an observed prevalence
difference into the two parts with an interval attached.

Demography and labour economics have had one since Kitagawa (1955, JASA
50(272):1168-1194, "Components of a Difference Between Two Rates"), which
became the Blinder-Oaxaca decomposition (Blinder 1973, J Human Resources
8(4):436-455; Oaxaca 1973, Int Econ Rev 14(3):693-709) and, for binary
outcomes, Fairlie's non-linear extension (Fairlie 2005, J Econ Soc Meas
30(4):305-316).

This module uses the **non-parametric Kitagawa form** rather than a fitted
Blinder-Oaxaca or Fairlie decomposition, and the reason is the shape of
bacterial data rather than convenience. A regression decomposition needs the
grouping variable as design columns; a sequence-type variable here has of the
order of a hundred levels on a few hundred isolates, dominated by one clone and
with a long singleton tail, so a logit on ST dummies is separated and
overfitted before it is asked anything. ``benchmarks/decomposition_vs_glm.py``
measures that head to head. Working directly with lineage shares ``w`` and
within-lineage rates ``p`` avoids the model entirely:

    P_A - P_B = sum_l (w_A - w_B) * (p_A + p_B) / 2      <- composition
              + sum_l (w_A + w_B) / 2 * (p_A - p_B)      <- within-lineage

The identity is exact by construction and ``identity_residual`` is reported so
a reader can confirm it rather than trust it.

**A lineage seen in only one collection needs a convention, and the choice of
convention is what makes the within-lineage component interpretable.** Its rate
in the other collection does not exist, and the identity holds for any value
substituted for it, so the substitution is a definition rather than an
estimate. Setting the absent rate equal to the observed one makes that
lineage's within-lineage term exactly zero, which has a consequence worth
stating plainly:

    the within-lineage component is then a function only of lineages observed
    in **both** collections, and does not depend on the convention at all.

Everything a lineage contributes by appearing or disappearing is charged to
composition, which is the demographic reading - a lineage that is not there is
a fact about the mix, not about a rate. That is the argument for this
convention over any other, and it is why the component a laboratory acts on is
the identified one.

That identification has a price, and the price is reported rather than
described. ``shared_support_isolate_share`` is the fraction of isolates sitting
in lineages common to both collections: the support the within-lineage
component is actually estimated on. Below ``min_shared_support`` the component
is returned with ``within_lineage_estimable`` false and a signed margin, in the
same form as the other gates in this package, because two collections that
barely overlap are describing lineage turnover and a within-lineage rate read
off them is an artefact of the convention rather than a measurement.
``benchmarks/decomposition_calibration.py`` measures interval coverage against
that share and is the evidence for the default.

``turnover_share`` reports the complementary quantity: the fraction of the
moving mass carried by lineages seen in only one collection. High turnover
means the population was replaced; low turnover means the same lineages
shifted.

The second function answers a smaller question that surveillance reporting
conflates. ``prevalence_per_isolate`` and ``prevalence_per_lineage`` are two
different estimands, not a biased and an unbiased version of one:

* per isolate - what fraction of the isolates you encounter carry the trait.
  Correct for clinical burden, and the number surveillance reports today.
* per lineage - what fraction of the distinct lineages carry it. Correct for
  diversity and emergence, and not reported anywhere.

They separate whenever sampling is uneven across lineages, which for a
submission-driven collection is always. Neither is the true one; a report that
does not say which it means is the problem.

Both functions need a lineage label and nothing else - no phylogeny, no
alignment, no core genome. A sequence type is enough, which is what a
diagnostic laboratory has. For a tree-aware test of independent trait
acquisition see :mod:`amr_clonalshare.archephy`; for whether a *partition*
is a clonal relabelling see :mod:`amr_clonalshare.lineage`.

An agent panel is a family of tests, and the agents in it are not independent:
cross-resistance makes a macrolide block move as one trait. ``decompose_panel``
runs the decomposition across a panel, controls the false discovery rate within
each component family, and reports the effective number of independent agents
so the size of the family is visible next to its count.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import beta as _beta

from .stats import (benjamini_hochberg, effective_dimension, fisher_exact_p,
                    permutation_pvalue)

__all__ = ["decompose_prevalence_difference", "decompose_panel",
           "lineage_resolved_prevalence", "trait_concentration"]

_MISSING = ("nan", "None", "", "__missing__", "NA", "<NA>")

#: Default lower bound on the isolate share held by lineages common to both
#: collections. Below it the within-lineage component is reported but marked
#: not estimable.
#:
#: The value is measured rather than chosen. Across 720 scenario cells and
#: 720,000 simulated decompositions in
#: ``benchmarks/decomposition_calibration.py``, coverage of the nominal 95 %
#: interval for the within-lineage component runs 0.82 below a shared support
#: of 0.4, 0.89 from 0.4 to 0.5, 0.91 from 0.5 to 0.6, 0.93 from 0.6 to 0.8 and
#: 0.94 above 0.8. Above 0.8 the fifth percentile across cells is 0.92; below
#: 0.5 the worst cell is 0.47. The composition component keeps nominal coverage
#: throughout and is not gated.
DEFAULT_MIN_SHARED_SUPPORT = 0.8

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
    ok &= ~np.isin(lin_s, _MISSING)
    ok &= np.isfinite(arr)
    return arr, lin_s, ok


def _lineage_present_mask(lineage: Sequence) -> np.ndarray:
    """Which isolates carry a usable lineage label, whatever the trait is.

    ``_clean_with_mask`` returns "label present **and** trait finite", which is
    the right mask to run the analysis on and the wrong one for the numbers
    reported beside it: a NaN in the trait was published as a missing lineage
    label under ``n_dropped_missing_lineage`` and as a fall in
    ``label_coverage_a``, which are claims SWX-057 and SWX-071. The coverage
    and the dropped-label counts read this mask instead; the trait-missing
    count is reported separately.
    """
    lin = pd.Series(list(lineage), dtype="object")
    ok = np.asarray(lin.notna().to_numpy(), dtype=bool, copy=True)
    lin_s = np.asarray(lin.fillna("__missing__").astype(str).to_numpy(),
                       dtype=object)
    ok &= ~np.isin(lin_s, _MISSING)
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

    The percentile interval and this p-value are the same statement read two
    ways: the p-value is the smallest ``alpha`` at which the ``1 - alpha``
    interval still excludes zero. The ``+1`` follows Phipson and Smyth for the
    same reason it does in a permutation test - a resampled tail count of zero
    is not evidence that the tail is empty.
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


def _label_availability(y_a: np.ndarray, keep_a: np.ndarray,
                        y_b: np.ndarray, keep_b: np.ndarray,
                        alpha: float) -> Dict[str, Any]:
    """Is the labelled subset representative of the collection it came from?

    A lineage-resolved statistic is computed on the isolates that carry a
    lineage label, and it describes the two collections only if labelling was
    unrelated to what is being measured. Two things have to go wrong together
    before that fails, and only together:

    * **coverage differs between the collections.** If both arms are labelled
      to the same degree, the same selection applies to both and a difference
      between them is largely unaffected.
    * **the trait differs between labelled and unlabelled isolates.** If
      labelling is unrelated to the trait, missingness costs precision and not
      validity, however uneven it is.

    When both hold, the difference computed on labelled isolates is a
    difference between two differently selected subsets and can carry the
    opposite sign to the collections it is drawn from. That is not a
    hypothetical: on the shipped *S. suis* cohort, ceftiofur non-wild-type
    prevalence rises by 9 points between two United Kingdom periods and falls
    by 11 points on the sequence-typed isolates within them, because only 40 %
    of the later period carries a sequence type and the untyped isolates are
    the more resistant ones.

    Both tests are Fisher exact on the trait, which requires a binary trait; a
    non-binary trait is reported as not assessable rather than tested with a
    rule that does not apply to it.
    """
    binary = bool(np.all(np.isin(y_a[np.isfinite(y_a)], (0.0, 1.0)))
                  and np.all(np.isin(y_b[np.isfinite(y_b)], (0.0, 1.0))))
    coverage_a = float(keep_a.mean()) if keep_a.size else float("nan")
    coverage_b = float(keep_b.mean()) if keep_b.size else float("nan")
    coverage_p = fisher_exact_p(int(keep_a.sum()), int((~keep_a).sum()),
                                int(keep_b.sum()), int((~keep_b).sum()))

    gaps: Dict[str, Any] = {}
    trait_p: Dict[str, Optional[float]] = {}
    for side, y, keep in (("a", y_a, keep_a), ("b", y_b, keep_b)):
        unlabelled = ~keep & np.isfinite(y)
        labelled = keep & np.isfinite(y)
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
    representative = not (differential and informative)
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
        "labels_representative": bool(representative),
        "assessable": bool(binary),
        "note": ("the decomposition describes the labelled isolates. It "
                 "describes the collections themselves only while "
                 "labels_representative holds: coverage differing between the "
                 "collections and the trait differing between labelled and "
                 "unlabelled isolates are harmless apart and not together"),
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
    min_shared_support : the within-lineage component is marked not estimable
        when the isolate share held by lineages common to both collections
        falls below this. Set to 0 to report the component unconditionally.

    Returns
    -------
    dict carrying the two components with percentile intervals, two-sided
    bootstrap p-values, the exactness residual of the identity, the lineage
    coverage of each collection, and ``turnover_share`` - the fraction of the
    moving mass held by lineages seen in only one collection.
    ``within_lineage`` is identified on the shared lineages alone;
    ``shared_support_isolate_share`` is the coverage it rests on and
    ``within_lineage_estimable`` states whether that coverage clears the gate.
    """
    if rng is None:
        rng = np.random.default_rng()
    all_ya, all_lina, keep_a = _clean_with_mask(y_a, lineage_a)
    all_yb, all_linb, keep_b = _clean_with_mask(y_b, lineage_b)
    label_a, label_b = (_lineage_present_mask(lineage_a),
                        _lineage_present_mask(lineage_b))
    ya, lina, drop_a = all_ya[keep_a], all_lina[keep_a], int((~label_a).sum())
    yb, linb, drop_b = all_yb[keep_b], all_linb[keep_b], int((~label_b).sum())
    drop_trait = int((~np.isfinite(all_ya)).sum() + (~np.isfinite(all_yb)).sum())
    availability = _label_availability(all_ya, label_a, all_yb, label_b,
                                       label_alpha)
    if ya.size == 0 or yb.size == 0:
        return {"status": "skipped",
                "reason": "no isolate with a usable lineage label in one or "
                          "both collections",
                "n_a": int(ya.size), "n_b": int(yb.size),
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

    boot = _bootstrap_components(codes_a, ya, codes_b, yb, n_lin,
                                 int(n_boot), rng)
    if boot.shape[0]:
        ci = np.percentile(boot, [2.5, 97.5], axis=0)
        se = boot.std(axis=0, ddof=1) if boot.shape[0] > 1 else np.zeros(2)
        p_comp = _bootstrap_pvalue(boot[:, 0])
        p_within = _bootstrap_pvalue(boot[:, 1])
    else:
        ci = np.full((2, 2), np.nan)
        se = np.full(2, np.nan)
        p_comp = p_within = float("nan")

    estimable = (support >= float(min_shared_support)
                 and bool(availability["labels_representative"]))
    return {
        "status": "ok",
        "lineage_label_availability": availability,
        "n_a": int(ya.size), "n_b": int(yb.size),
        "n_dropped_missing_lineage": int(drop_a + drop_b),
        "n_dropped_non_finite_trait": int(drop_trait),
        "n_lineages": int(n_lin),
        "n_lineages_shared": int(shared.sum()),
        "n_lineages_only_a": int(only_a.sum()),
        "n_lineages_only_b": int(only_b.sum()),
        "shared_support_isolate_share": support,
        "prevalence_a": float(ya.mean()),
        "prevalence_b": float(yb.mean()),
        "difference": diff,
        "composition": comp,
        "within_lineage": within,
        "composition_ci95": [float(ci[0, 0]), float(ci[1, 0])],
        "within_lineage_ci95": [float(ci[0, 1]), float(ci[1, 1])],
        "composition_se": float(se[0]),
        "within_lineage_se": float(se[1]),
        "composition_p": p_comp,
        "within_lineage_p": p_within,
        "p_value_floor": float(2.0 / (int(n_boot) + 1)) if n_boot else None,
        "identity_residual": float(comp + within - diff),
        "turnover_share": share,
        "within_lineage_estimable": bool(estimable),
        "not_estimable_because": sorted(
            ([] if support >= float(min_shared_support)
             else ["shared lineage support below the threshold"])
            + ([] if availability["labels_representative"]
               else ["lineage labels are differentially missing and the "
                     "missingness is associated with the trait"])),
        "shared_support_threshold": float(min_shared_support),
        "shared_support_margin": float(support - float(min_shared_support)),
        "n_boot": int(n_boot),
        "convention": "a lineage absent from one collection is given that "
                      "collection's rate equal to the observed one; this makes "
                      "within_lineage a function of shared lineages only, and "
                      "so independent of the convention",
        "note": "the difference is described, not explained: two collections "
                "separated by country or period also differ in sampling frame "
                "and laboratory method, and neither component is a causal "
                "effect",
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


def decompose_panel(
    y_a: "pd.DataFrame", lineage_a: Sequence,
    y_b: "pd.DataFrame", lineage_b: Sequence,
    *,
    n_boot: int = 2000,
    rng: Optional[np.random.Generator] = None,
    q: float = 0.05,
    min_shared_support: float = DEFAULT_MIN_SHARED_SUPPORT,
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
    agents is the participation ratio of the panel's correlation eigenspectrum,
    the same construction the package uses for contributing layers: a panel of
    thirteen agents carrying four distinct resistance phenotypes has an
    effective size near four, and that number belongs beside the count of
    discoveries.

    Parameters
    ----------
    y_a, y_b : isolate-by-agent indicator frames for the two collections. They
        must carry the same agent columns.
    lineage_a, lineage_b : lineage label per isolate in each frame.
    q : target false discovery rate within each component family.

    Returns
    -------
    dict with ``per_agent`` results, each carrying ``composition_q`` and
    ``within_lineage_q``, and a ``family`` block recording the panel size, its
    effective size, and the discovery counts before and after control.
    """
    if rng is None:
        rng = np.random.default_rng()
    agents = [str(c) for c in y_a.columns]
    if [str(c) for c in y_b.columns] != agents:
        raise ValueError("both collections must carry the same agent columns")

    per_agent: Dict[str, dict] = {}
    for agent in agents:
        per_agent[agent] = decompose_prevalence_difference(
            y_a[agent].to_numpy(), lineage_a,
            y_b[agent].to_numpy(), lineage_b,
            n_boot=n_boot, rng=rng, min_shared_support=min_shared_support)

    usable = [a for a in agents if per_agent[a].get("status") == "ok"]
    family: Dict[str, Any] = {
        "n_agents": len(agents),
        "n_agents_decomposed": len(usable),
        "q": float(q),
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
        for position, agent in enumerate(usable):
            per_agent[agent][f"{key}_q"] = float(adjusted[position])
            per_agent[agent][f"{key}_discovery"] = bool(reject[position])
        family[f"n_{label}_nominal"] = int(sum(
            1 for a in usable
            if per_agent[a][f"{key}_ci95"][0] * per_agent[a][f"{key}_ci95"][1] > 0))
        family[f"n_{label}_discoveries"] = int(reject.sum())
        family[f"n_{label}_discoveries_under_independence"] = int(reject_pos.sum())

    stacked = pd.concat([y_a[agents], y_b[agents]], axis=0).to_numpy(dtype=float)
    family["effective_independent_agents"] = float(effective_dimension(stacked))
    family["note"] = (
        "the effective count is the participation ratio of the panel "
        "correlation eigenspectrum; a panel whose agents share a resistance "
        "mechanism is worth fewer independent tests than it has columns")
    return {"status": "ok", "per_agent": per_agent, "family": family}


def lineage_resolved_prevalence(
    y: Sequence, lineage: Sequence,
    *,
    n_boot: int = 2000,
    rng: Optional[np.random.Generator] = None,
) -> Dict[str, Any]:
    """Report trait prevalence per isolate and per lineage, with an interval.

    Parameters
    ----------
    y : per-isolate trait indicator.
    lineage : lineage label per isolate; missing labels are dropped.
    n_boot : replicates of a **two-stage cluster bootstrap**. Lineages are
        resampled first, because the per-lineage estimand averages over
        lineages and its uncertainty is dominated by how many distinct ones
        were seen; isolates are then resampled inside each drawn lineage, so
        that a rate measured on two isolates is not treated as if it were
        measured exactly. Resampling isolates alone gives an interval that is
        far too narrow; resampling lineages alone omits the second stage
        entirely. Because a lineage of size ``m`` with rate ``r`` resampled to
        size ``m`` yields ``Binomial(m, r) / m``, the second stage is drawn
        exactly rather than by index.

    Returns
    -------
    dict with both estimands, their difference, the share held by the largest
    lineage, and the fraction of lineages carrying the trait at all.
    """
    if rng is None:
        rng = np.random.default_rng()
    all_y, all_lin, keep = _clean_with_mask(y, lineage)
    arr, lin = all_y[keep], all_lin[keep]
    dropped = int((~_lineage_present_mask(lineage)).sum())
    dropped_trait = int((~np.isfinite(all_y)).sum())
    if arr.size == 0:
        return {"status": "skipped",
                "reason": "no isolate with a usable lineage label",
                "n_dropped_missing_lineage": int(dropped),
                "n_dropped_non_finite_trait": int(dropped_trait)}

    lineages = sorted(set(lin.tolist()))
    index = {value: position for position, value in enumerate(lineages)}
    codes = np.fromiter((index[v] for v in lin), dtype=np.intp, count=lin.size)
    n_lin = len(lineages)
    counts = np.bincount(codes, minlength=n_lin).astype(float)
    sums = np.bincount(codes, weights=arr, minlength=n_lin)
    rates = sums / counts

    per_isolate = float(arr.mean())
    per_lineage = float(rates.mean())
    binary = bool(np.all((arr == 0) | (arr == 1)))
    if n_boot > 0 and n_lin > 0:
        # Chunked for the same reason the other resamplers are: a cohort with
        # thousands of lineages would otherwise size the working array by the
        # replicate count as well as by the cohort.
        chunk = max(1, min(int(n_boot), _CHUNK_CELLS // max(n_lin, 1)))
        draws = np.empty(int(n_boot), dtype=float)
        done = 0
        while done < int(n_boot):
            size = min(chunk, int(n_boot) - done)
            pick = rng.integers(0, n_lin, (size, n_lin))
            drawn_sizes = counts[pick]
            drawn_rates = rates[pick]
            if binary:
                resampled = rng.binomial(
                    drawn_sizes.astype(np.int64),
                    np.clip(drawn_rates, 0.0, 1.0)) / drawn_sizes
            else:
                # A non-binary trait has no exact second stage; the
                # between-lineage stage is reported alone and the note says so.
                resampled = drawn_rates
            draws[done:done + size] = resampled.mean(axis=1)
            done += size
        lo, hi = (float(np.percentile(draws, 2.5)),
                  float(np.percentile(draws, 97.5)))
    else:
        lo = hi = float("nan")

    return {
        "status": "ok",
        "n_isolates": int(arr.size),
        "n_lineages": int(n_lin),
        "n_dropped_missing_lineage": int(dropped),
        "n_dropped_non_finite_trait": int(dropped_trait),
        "largest_lineage_share": float(counts.max() / counts.sum()),
        "prevalence_per_isolate": per_isolate,
        "prevalence_per_lineage": per_lineage,
        "prevalence_per_lineage_ci95": [lo, hi],
        "fraction_of_lineages_with_any": float((rates > 0).mean()),
        "difference_per_lineage_minus_per_isolate": per_lineage - per_isolate,
        "n_boot": int(n_boot),
        "bootstrap": ("two-stage cluster bootstrap: lineages, then isolates "
                      "within each drawn lineage" if binary else
                      "lineage-stage bootstrap only; the trait is not binary "
                      "so the within-lineage stage has no exact form"),
        "note": "two estimands, not a biased and an unbiased one: per isolate "
                "is clinical burden, per lineage is diversity. They separate "
                "whenever sampling across lineages is uneven",
    }


def trait_concentration(
    y: Sequence, lineage: Sequence,
    *,
    n_perm: int = 500,
    rng: Optional[np.random.Generator] = None,
) -> Dict[str, Any]:
    """How concentrated in lineages is a trait, and by how much beyond chance.

    Two instruments, both borrowed rather than invented, and neither in use on
    resistance data.

    The **Herfindahl-Hirschman index** is the competition-authority measure of
    market concentration: the sum of squared shares, here the shares of the
    carrier pool held by each lineage. Its reciprocal is the *effective number
    of lineages carrying the trait* - one when a single clone holds it, ``L``
    when ``L`` lineages hold it equally. The package already reports an
    effective number of contributing layers, so the reading carries over: a
    determinant with an effective 2 carriers is a clonal object whatever its
    prevalence, and one with an effective 40 is dispersed.

    Concentration alone is not evidence, because an uneven cohort produces an
    uneven carrier pool with no biology involved. The tested quantity is
    therefore the **departure from proportional carriage**: the divergence of
    the carrier distribution from the lineage-abundance distribution,
    ``sum_l s_l log2(s_l / w_l)``, in bits - the Theil form of the
    inequality-economics decomposition. It is exactly zero when carriage is
    proportional to abundance, however uneven the lineages are, and the
    permutation null reshuffles carriage across isolates with the number of
    carriers held fixed, so only its lineage placement varies.

    **The departure is a magnitude and says nothing about which way.** It is
    large both when a trait piles into one clone and when a trait avoids the
    dominant one, and on the shipped *S. suis* panel it is the second case that
    fires: beta-lactam and tiamulin non-wild-type carriage sits in more
    lineages than chance allows, because the dominant sequence type is nearly
    free of it. ``direction`` therefore accompanies every result and must be
    read with it - a significant departure labelled ``dispersed`` is the
    opposite finding to one labelled ``concentrated``, and the p-value alone
    cannot tell them apart.

    The permutation p-value has a floor at ``1 / (n_perm + 1)`` and a panel run
    at a small budget reports that floor for every agent, which reads as
    thirteen strong results and is one statement about the budget.
    ``tail_probability_ci95`` is the exact binomial interval for the tail the
    permutations estimate and ``resolved_at_alpha_0_05`` states whether that
    interval settles the decision, in the same form the continuum test uses.

    Parameters
    ----------
    y : per-isolate trait indicator; non-zero counts as carrying.
    lineage : lineage label per isolate; missing labels are dropped.
    n_perm : permutations of carriage across isolates.

    Returns
    -------
    dict with the effective number of carrying lineages, the raw index, the
    excess in bits, and the permutation null against which the excess is read.
    """
    if rng is None:
        rng = np.random.default_rng()
    all_y, all_lin, keep = _clean_with_mask(y, lineage)
    arr, lin = all_y[keep], all_lin[keep]
    dropped = int((~_lineage_present_mask(lineage)).sum())
    dropped_trait = int((~np.isfinite(all_y)).sum())
    if arr.size == 0:
        return {"status": "skipped",
                "reason": "no isolate with a usable lineage label",
                "n_dropped_missing_lineage": int(dropped),
                "n_dropped_non_finite_trait": int(dropped_trait)}
    carrier = arr > 0
    n_carriers = int(carrier.sum())
    if n_carriers == 0:
        return {"status": "skipped",
                "reason": "no isolate carries the trait",
                "n_isolates": int(arr.size),
                "n_dropped_missing_lineage": int(dropped),
                "n_dropped_non_finite_trait": int(dropped_trait)}

    codes, lineages = pd.factorize(pd.Series(lin, dtype="object"))
    codes = np.asarray(codes, dtype=np.intp)
    n_lin = len(lineages)
    n_iso = arr.size
    w = np.bincount(codes, minlength=n_lin).astype(float) / n_iso

    def statistics(mask: np.ndarray) -> Tuple[float, float]:
        s = np.bincount(codes[mask], minlength=n_lin).astype(float)
        s /= s.sum()
        hhi = float((s ** 2).sum())
        nz = s > 0
        bits = float((s[nz] * np.log2(s[nz] / w[nz])).sum())
        return hhi, bits

    hhi, bits = statistics(carrier)
    null = _concentration_null(codes, w, n_iso, n_lin, n_carriers,
                               int(n_perm), rng)
    sd = float(null[:, 1].std(ddof=1)) if n_perm > 1 else 0.0
    eff, null_eff = 1.0 / hhi, float(null[:, 0].mean())
    lo, hi = np.percentile(null[:, 0], [2.5, 97.5])
    direction = ("concentrated" if eff < lo else
                 "dispersed" if eff > hi else "proportional")
    exceedances = int(np.sum(null[:, 1] >= bits))
    tail_low, tail_high = _clopper_pearson(exceedances, int(n_perm))

    return {
        "status": "ok",
        "n_isolates": int(n_iso),
        "n_lineages": int(n_lin),
        "n_dropped_missing_lineage": int(dropped),
        "n_dropped_non_finite_trait": int(dropped_trait),
        "n_carriers": n_carriers,
        "prevalence": float(carrier.mean()),
        "n_lineages_carrying": int((np.bincount(codes[carrier],
                                                minlength=n_lin) > 0).sum()),
        "herfindahl_index": hhi,
        "effective_number_of_lineages": eff,
        "null_effective_number_ci95": [float(lo), float(hi)],
        "departure_from_proportional_bits": bits,
        "null_mean_effective_number": null_eff,
        "null_mean_departure_bits": float(null[:, 1].mean()),
        "null_sd_departure_bits": sd,
        "z": float((bits - null[:, 1].mean()) / sd) if sd > 0 else float("nan"),
        "p_value": permutation_pvalue(null[:, 1], bits, tail="greater"),
        "p_value_floor": float(1.0 / (int(n_perm) + 1)),
        "permutation_exceedances": exceedances,
        "tail_probability_ci95": {
            "method": "Clopper-Pearson exact binomial",
            "confidence_level": 0.95,
            "low": tail_low,
            "high": tail_high,
        },
        "resolved_at_alpha_0_05": bool(tail_high < 0.05 or tail_low > 0.05),
        "n_perm": int(n_perm),
        "direction": direction,
        "note": "the departure in bits is a magnitude: it fires on carriage "
                "piling into few lineages and on carriage avoiding the "
                "dominant one alike. Read `direction` with the p-value or the "
                "sign of the finding is lost",
    }


def _concentration_null(codes: np.ndarray, w: np.ndarray, n_iso: int,
                        n_lineages: int, n_carriers: int, n_perm: int,
                        rng: np.random.Generator) -> np.ndarray:
    """Permutation null for carriage placement, drawn in chunks.

    Carriage is reassigned to ``n_carriers`` isolates chosen without
    replacement, so prevalence is held and only the lineage placement varies.
    ``argpartition`` selects that subset in linear time per replicate.
    """
    out = np.empty((max(n_perm, 1), 2), dtype=float)
    if n_perm <= 0:
        out[0] = (1.0 / max((w ** 2).sum(), 1e-300), 0.0)
        return out
    chunk = max(1, min(n_perm, _CHUNK_CELLS // max(n_iso, 1)))
    log_w = np.where(w > 0, w, 1.0)
    done = 0
    while done < n_perm:
        size = min(chunk, n_perm - done)
        keys = rng.random((size, n_iso))
        picked = np.argpartition(keys, n_carriers - 1, axis=1)[:, :n_carriers]
        picked_codes = codes[picked]
        offset = (np.arange(size, dtype=np.intp)[:, None] * n_lineages
                  + picked_codes).ravel()
        counts = np.bincount(offset, minlength=size * n_lineages
                             ).reshape(size, n_lineages).astype(float)
        s = counts / counts.sum(axis=1, keepdims=True)
        out[done:done + size, 0] = 1.0 / (s ** 2).sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            terms = np.where(s > 0, s * np.log2(s / log_w[None, :]), 0.0)
        out[done:done + size, 1] = terms.sum(axis=1)
        done += size
    return out
