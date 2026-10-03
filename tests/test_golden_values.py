"""Every public function of the modules that compute a reported number, on
fixtures that reach its branches, compared in full against a committed record.

The other tests check what a result must be: an oracle, an identity, a bound,
an invariance. These pin what it is. Each case calls one public function on a
small fixture chosen to pass through a branch the others do not -- one
stratum and several, one panel and overlapping panels, few lineages and more
than the exact search reaches, singletons, sampling units, missing values,
censored and uncensored readings -- and every field of the result, numbers,
labels and messages alike, is compared with the record in
`tests/golden/values/`, numbers to 1e-8 relative. A change to any reported
value therefore fails here, and the diff of the record describes it.

Regenerate after an intended change with `python scripts/update_golden.py`,
read the diff, and commit it with the change that caused it.
"""
from __future__ import annotations

import dataclasses
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from amr_clonalshare import (attribution, censored, clonality, comparison, evalues, latent_order,
                             mic_order, missingness, phenotype, qc, stats)

VALUES = Path(__file__).resolve().parent / "golden" / "values"
SMALL = dict(folds=5, repeats=5, n_boot=199, n_perm=99)


def plain(x):
    """A JSON-ready copy of a result: dataclasses, arrays, frames and numpy
    scalars as plain values; non-finite numbers as their names."""
    if dataclasses.is_dataclass(x) and not isinstance(x, type):
        return plain(dataclasses.asdict(x))
    if isinstance(x, dict):
        return {str(k): plain(v) for k, v in x.items()}
    if isinstance(x, pd.DataFrame):
        return {"columns": [str(c) for c in x.columns], "index": plain(list(x.index)),
                "rows": plain(x.astype(object).where(x.notna(), None).values.tolist()),
                "attrs": plain(dict(x.attrs))}
    if isinstance(x, pd.Series):
        return plain(x.astype(object).where(x.notna(), None).tolist())
    if isinstance(x, (list, tuple, set, frozenset)):
        items = sorted(x, key=str) if isinstance(x, (set, frozenset)) else x
        return [plain(v) for v in items]
    if isinstance(x, np.ndarray):
        return plain(x.tolist())
    if isinstance(x, (np.bool_, bool)):
        return bool(x)
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, (np.floating, float)):
        v = float(x)
        return v if math.isfinite(v) else ("nan" if v != v else ("inf" if v > 0 else "-inf"))
    if x is None or isinstance(x, (int, str)):
        return x
    if x is pd.NA:
        return None
    return repr(x)


# --------------------------------------------------------------------- fixtures

def _panel(y, cuts):
    k = np.searchsorted(cuts, y)
    ext = np.r_[-np.inf, cuts, np.inf]
    return ext[k], ext[k + 1]


def _latent(G, m, effect, seed, cuts):
    rng = np.random.default_rng(seed)
    sizes = np.full(G, m) if np.isscalar(m) else np.asarray(m)
    code = np.repeat(np.arange(G), sizes)
    y = effect * rng.normal(0, 1, G)[code] + rng.normal(0, 1, code.size)
    lo, hi = _panel(y, cuts)
    return lo, hi, code, y


def _two_labs(seed):
    rng = np.random.default_rng(seed)
    code = np.repeat(np.arange(10), 6)
    lab = np.where(np.arange(code.size) % 3 == 0, "x", "y")
    y = rng.normal(0, 1, 10)[code] + .8 * (lab == "y") + rng.normal(0, 1, code.size)
    lo, hi = np.empty(code.size), np.empty(code.size)
    for name, cuts in (("x", np.arange(-2., 3.)), ("y", np.arange(-1.5, 3., 1.5))):
        m = lab == name
        lo[m], hi[m] = _panel(y[m], cuts)
    return lo, hi, code, lab


def _overlapping(seed):
    """One stratum read on two panels whose wells overlap."""
    rng = np.random.default_rng(seed)
    code = np.repeat(np.arange(8), 6)
    y = rng.normal(0, 1, 8)[code] + rng.normal(0, 1, code.size)
    lo, hi = _panel(y, np.arange(-2., 3.))
    odd = np.arange(code.size) % 2 == 1
    lo[odd], hi[odd] = _panel(y[odd], np.arange(-1.5, 2.6, 1.))
    return lo, hi, code


def _calls(seed, G=12, prevalence=.3, effect=1.5, singletons=2, units=False):
    rng = np.random.default_rng(seed)
    sizes = rng.integers(2, 9, G)
    code = np.r_[np.repeat(np.arange(G), sizes), np.arange(G, G + singletons)]
    z = effect * rng.normal(0, 1, G + singletons)[code] + rng.normal(0, 1, code.size)
    y = (z > np.quantile(z, 1 - prevalence)).astype(float)
    labels = np.array([f"L{c:02d}" for c in code], dtype=object)
    farm = np.array([f"farm{int(v)}" for v in rng.integers(0, 4, code.size)], dtype=object)
    return y, labels, farm, z


def _majority():
    lineage = np.repeat(np.arange(10), 4)
    return np.repeat(np.arange(10) % 2, 4).astype(float), lineage


def _strong_labs(seed):
    """Two laboratories on different panels, a strong lineage effect, and
    lineages with a single isolate in one laboratory."""
    rng = np.random.default_rng(seed)
    code = np.repeat(np.arange(12), 9)
    lab = np.where(rng.random(code.size) < .3, "x", "y")
    y = 2. * rng.normal(0, 1, 12)[code] + .8 * (lab == "y") + rng.normal(0, 1, code.size)
    lo, hi = np.empty(code.size), np.empty(code.size)
    for name, cuts in (("x", np.arange(-3., 3.5, 1.)), ("y", np.arange(-2.5, 3.5, .5))):
        m = lab == name
        lo[m], hi[m] = _panel(y[m], cuts)
    return lo, hi, code, lab


def _strong_overlapping(seed):
    rng = np.random.default_rng(seed)
    code = np.repeat(np.arange(9), 10)
    y = 2. * rng.normal(0, 1, 9)[code] + rng.normal(0, 1, code.size)
    lo, hi = _panel(y, np.arange(-3., 3.5, 1.))
    odd = np.arange(code.size) % 2 == 1
    lo[odd], hi[odd] = _panel(y[odd], np.arange(-2.5, 3., 1.))
    return lo, hi, code


def _eight_lineages():
    counts = [[5, 0, 11, 11, 0, 5], [2, 2, 11, 11, 2, 2], [48, 1, 16, 16, 1, 48], [6, 1, 0, 0, 1, 6],
              [28, 3, 25, 25, 3, 28], [6, 41, 12, 12, 41, 6], [4, 16, 0, 0, 16, 4], [5, 29, 0, 0, 29, 5]]
    ext = np.r_[-np.inf, np.arange(5.), np.inf]
    lo, hi, lineage = [], [], []
    for g, row in enumerate(counts):
        for k, c in enumerate(row):
            lo += [ext[k]] * c
            hi += [ext[k + 1]] * c
            lineage += [g] * c
    return np.array(lo), np.array(hi), np.array(lineage)


def _collections(seed):
    rng = np.random.default_rng(seed)
    lin_a = rng.choice(list("ABCDEFG"), 160, p=[.25, .2, .15, .15, .1, .1, .05])
    lin_b = rng.choice(list("ABCDEFH"), 140, p=[.1, .1, .2, .2, .15, .15, .1])
    rate = dict(zip("ABCDEFGH", [.1, .6, .3, .5, .2, .8, .4, .7]))
    y_a = (rng.random(160) < [rate[l] for l in lin_a]).astype(float)
    y_b = (rng.random(140) < [rate[l] + .1 * (l in "AB") for l in lin_b]).astype(float)
    return y_a, lin_a, y_b, lin_b, rng


def _matched(seed):
    rng = np.random.default_rng(seed)
    a = np.repeat(np.arange(24), rng.integers(2, 6, 24))
    b = (a // 2).astype(object)
    b[rng.random(a.size) < .1] = None
    y = (rng.random(a.size) < np.where(a % 3 == 0, .7, .2)).astype(float)
    return y, a.astype(object), b


def _sir(values):
    return pd.DataFrame({"Strain_ID": [f"{i:03}" for i in range(len(values))], "antibiotic": "a",
                         "resistant_phenotype": values})


def _qc_frames():
    panel = pd.DataFrame({"complete": [0, 1] * 15,
                          "sparse": [0, np.nan, 1, np.nan, 0, np.nan] * 5,
                          "rare": [0] * 28 + [1, 1]})
    meta = pd.DataFrame({"lineage": [None, "a"] + list("aabbbccddeeffgghhiijjkkllmmn")})
    return panel, meta


# ------------------------------------------------------------------------ cases

def _cases():
    c = {}
    # latent_order: the bounds, the certificate and the exact upper end
    for name, (G, m, eff, seed, cuts) in {
            "bounds_six_lineages": (6, 8, 1., 1, np.arange(-2., 3.)),
            "bounds_sixteen_lineages": (16, 5, .8, 2, np.arange(-2., 2.5, .5)),
            "bounds_twenty_six_lineages": (26, 4, .6, 3, np.arange(-1.5, 2., .5)),
            "bounds_null_thirty_pairs": (30, 2, 0., 4, np.linspace(-1.2, 1.2, 8)),
            "bounds_uneven_sizes": (7, [2, 3, 12, 5, 2, 9, 4], 1.2, 5, np.arange(-2., 3.))}.items():
        c[name] = (lambda G=G, m=m, eff=eff, seed=seed, cuts=cuts:
                   latent_order.order_bounds(*_latent(G, m, eff, seed, cuts)[:3], seed=7))
    c["bounds_two_laboratories"] = lambda: latent_order.order_bounds(*_two_labs(6)[:3], strata=_two_labs(6)[3],
                                                                     seed=8)
    c["bounds_overlapping_panels"] = lambda: latent_order.order_bounds(*_overlapping(9), seed=3)
    c["bounds_alpha_ten_percent"] = lambda: latent_order.order_bounds(
        *_latent(9, 6, 1., 10, np.arange(-2., 3.))[:3], alpha=.1, starts=3, seed=4)
    c["bounds_audit_eight_lineages"] = lambda: latent_order.order_bounds(*_eight_lineages(), seed=0)
    c["bounds_exact_readings_with_singleton"] = lambda: latent_order.order_bounds(
        np.r_[0., 1, 2, 3, 4, 5, 6], np.r_[0., 1, 2, 3, 4, 5, 6], [0, 0, 1, 1, 2, 2, 3], seed=1)
    c["call_bounds_strata"] = lambda: latent_order.call_bounds(
        _calls(11)[0], _calls(11)[1], strata=_calls(11)[2], seed=2)
    c["call_bounds_majority"] = lambda: latent_order.call_bounds(*_majority(), seed=0)
    c["bounds_strong_ten_lineages"] = lambda: latent_order.order_bounds(
        *_latent(10, 12, 2., 61, np.arange(-3., 3.5, .75))[:3], seed=9)
    c["bounds_strong_eighteen_lineages"] = lambda: latent_order.order_bounds(
        *_latent(18, 8, 1.8, 62, np.arange(-3., 3.5, .75))[:3], seed=9)
    c["bounds_sixty_lineages"] = lambda: latent_order.order_bounds(
        *_latent(60, 4, 1.2, 63, np.arange(-2., 2.5, .5))[:3], seed=9)
    c["bounds_strong_two_laboratories"] = lambda: latent_order.order_bounds(*_strong_labs(64)[:3],
                                                                            strata=_strong_labs(64)[3], seed=10)
    c["bounds_strong_overlapping_panels"] = lambda: latent_order.order_bounds(*_strong_overlapping(65), seed=11)
    c["call_bounds_strong_strata"] = lambda: latent_order.call_bounds(
        _calls(66, G=14, prevalence=.4, effect=3.)[0], _calls(66, G=14, prevalence=.4, effect=3.)[1],
        strata=_calls(66, G=14, prevalence=.4, effect=3.)[2], alpha=.1, seed=12)
    c["bounds_two_lineages"] = lambda: latent_order.order_bounds(*_latent(2, 9, 1.5, 67, np.arange(-2., 3.))[:3],
                                                                 seed=13)
    c["bounds_one_reading"] = lambda: latent_order.order_bounds(np.zeros(12), np.ones(12), np.repeat(np.arange(4), 3),
                                                                seed=14)
    c["bounds_refused_alpha"] = lambda: latent_order.order_bounds(np.zeros(4), np.ones(4), [0, 0, 1, 1], alpha=1.)
    # mic_order: the share of the recorded ordering
    c["order_scores_one_stratum"] = lambda: mic_order.order_scores(*_latent(8, 5, 1., 12, np.arange(-2., 3.))[:2])
    c["order_scores_two_laboratories"] = lambda: mic_order.order_scores(*_two_labs(13)[:2], strata=_two_labs(13)[3])
    c["order_scores_overlapping"] = lambda: mic_order.order_scores(*_overlapping(14)[:2])
    c["mic_share_one_stratum"] = lambda: mic_order.mic_order_share(
        *_latent(10, 6, 1., 15, np.arange(-2., 3.))[:3], seed=1, **SMALL)
    c["mic_share_two_laboratories"] = lambda: mic_order.mic_order_share(
        *_two_labs(16)[:3], strata=_two_labs(16)[3], stratified_by="laboratory", seed=2, **SMALL)
    c["mic_share_units"] = lambda: mic_order.mic_order_share(
        *_latent(12, 5, .7, 17, np.arange(-2., 3.))[:3], units=np.tile(["u1", "u2", "u3"], 20), seed=3, **SMALL)
    c["mic_share_end_wells"] = lambda: mic_order.mic_order_share(
        *_latent(8, 6, 1., 18, np.array([-.5, .5]))[:3], seed=4, **SMALL)
    c["mic_share_pairs"] = lambda: mic_order.mic_order_share(
        *_latent(20, 2, .5, 19, np.arange(-2., 3.))[:3], seed=5, **SMALL)
    c["mic_share_majority_two_wells"] = lambda: mic_order.mic_order_share(
        np.where(_majority()[0] == 0, -np.inf, 0.), np.where(_majority()[0] == 0, 0., np.inf), _majority()[1],
        seed=1, **SMALL)
    c["mic_share_strong_two_laboratories_units"] = lambda: mic_order.mic_order_share(
        *_strong_labs(68)[:3], strata=_strong_labs(68)[3], units=np.arange(108) % 6, seed=6, **SMALL)
    c["mic_share_all_one_well"] = lambda: mic_order.mic_order_share(np.zeros(12), np.ones(12),
                                                                     np.repeat(np.arange(4), 3), seed=7, **SMALL)
    c["mic_share_missing_labels"] = lambda: mic_order.mic_order_share(
        *_latent(8, 6, 1., 69, np.arange(-2., 3.))[:2],
        np.where(np.arange(48) % 7 == 0, None, np.repeat(np.arange(8), 6).astype(object)), seed=8, **SMALL)
    # attribution: the share of a call or a trait
    c["share_call"] = lambda: attribution.clonal_share(_calls(20)[0], _calls(20)[1], seed=1, **SMALL)
    c["share_call_strata"] = lambda: attribution.clonal_share(_calls(21)[0], _calls(21)[1],
                                                              strata=_calls(21)[2], seed=2, **SMALL)
    c["share_call_units"] = lambda: attribution.clonal_share(_calls(22)[0], _calls(22)[1],
                                                             units=_calls(22)[2], seed=3, **SMALL)
    c["share_rare_call"] = lambda: attribution.clonal_share(_calls(23, G=15, prevalence=.06)[0],
                                                            _calls(23, G=15, prevalence=.06)[1], seed=4, **SMALL)
    c["share_trait"] = lambda: attribution.layer_clonal_share(_calls(24)[3], _calls(24)[1], seed=5, **SMALL)
    c["share_trait_categories"] = lambda: attribution.layer_clonal_share(
        np.round(_calls(25)[3]), _calls(25)[1], categories=np.round(_calls(25)[3]).astype(int), seed=6, **SMALL)
    c["share_majority"] = lambda: attribution.clonal_share(*_majority(), seed=1, **SMALL)
    c["share_null"] = lambda: attribution.clonal_share(_calls(26, effect=0.)[0], _calls(26, effect=0.)[1],
                                                       seed=7, **SMALL)
    c["share_without_interval"] = lambda: attribution.clonal_share(_calls(27)[0], _calls(27)[1], seed=8,
                                                                   folds=5, repeats=5, n_boot=0, n_perm=99)
    c["share_one_lineage"] = lambda: attribution.clonal_share([0., 1, 0, 1, 1], ["a"] * 5, seed=1, **SMALL)
    c["share_constant_trait"] = lambda: attribution.layer_clonal_share(np.full(30, 2.5), np.repeat(np.arange(10), 3),
                                                                       seed=9, **SMALL)
    c["share_ten_folds_on_pairs"] = lambda: attribution.clonal_share(
        _calls(28, G=20)[0], _calls(28, G=20)[1], seed=10, folds=10, repeats=5, n_boot=199, n_perm=99)
    c["share_missing_values_strata_units"] = lambda: attribution.layer_clonal_share(
        np.where(np.arange(_calls(29)[3].size) % 9 == 0, np.nan, _calls(29)[3]),
        np.where(np.arange(_calls(29)[3].size) % 11 == 0, None, _calls(29)[1]),
        strata=np.arange(_calls(29)[3].size) % 2, units=_calls(29)[2], seed=11, **SMALL)
    # comparison of two lineage definitions
    c["compare_merged_pairs"] = lambda: comparison.compare_lineage_definitions(
        _matched(30)[0], _matched(30)[1], _matched(30)[2], seed=1, folds=5, repeats=5, n_perm=49, n_boot=99)
    c["compare_split"] = lambda: comparison.compare_lineage_definitions(
        _matched(31)[0], _matched(31)[2], _matched(31)[1], name_a="coarse", name_b="fine", seed=2,
        folds=5, repeats=5, n_perm=49, n_boot=99)
    c["compare_missing_outcomes"] = lambda: comparison.compare_lineage_definitions(
        np.where(np.arange(_matched(32)[0].size) % 8 == 0, np.nan, _matched(32)[0]), _matched(32)[1], _matched(32)[2],
        ids=np.array([f"i{k:03d}" for k in range(_matched(32)[0].size)]), seed=3, folds=5, repeats=5, n_perm=49,
        n_boot=99)
    # clonality: the prevalence decomposition
    c["decompose_two_collections"] = lambda: clonality.decompose_prevalence_difference(
        *_collections(40)[:4], n_boot=199, rng=np.random.default_rng(1))
    c["decompose_units"] = lambda: clonality.decompose_prevalence_difference(
        *_collections(41)[:4], units_a=np.arange(160) % 20, units_b=np.arange(140) % 14, n_boot=199,
        rng=np.random.default_rng(2), min_shared_support=.5)
    c["decompose_panel"] = lambda: clonality.decompose_panel(
        pd.DataFrame({"x": _collections(42)[0], "y": 1. - _collections(42)[0]}), _collections(42)[1],
        pd.DataFrame({"x": _collections(42)[2], "y": _collections(42)[2]}), _collections(42)[3],
        n_boot=199, rng=np.random.default_rng(3), min_shared_support=.5)
    c["decompose_missing_values"] = lambda: clonality.decompose_prevalence_difference(
        np.where(np.arange(160) % 10 == 0, np.nan, _collections(43)[0]),
        np.where(np.arange(160) % 7 == 0, None, _collections(43)[1]).astype(object),
        _collections(43)[2], np.where(np.arange(140) % 5 == 0, None, _collections(43)[3]).astype(object),
        n_boot=199, rng=np.random.default_rng(4), min_shared_support=.5, label_alpha=.1)
    c["decompose_panel_q_ten_percent"] = lambda: clonality.decompose_panel(
        pd.DataFrame({"x": _collections(44)[0], "y": np.where(np.arange(160) % 6 == 0, np.nan, 1. - _collections(44)[0]),
                      "z": (np.arange(160) % 3 == 0).astype(float)}), _collections(44)[1],
        pd.DataFrame({"x": _collections(44)[2], "y": _collections(44)[2], "z": (np.arange(140) % 4 == 0).astype(float)}),
        _collections(44)[3], n_boot=199, rng=np.random.default_rng(5), q=.1, min_shared_support=.5)
    # evalues
    c["e_process"] = lambda: evalues.e_process(_calls(50)[0], _calls(50)[1], seed=1, repeats=5)
    c["e_process_null"] = lambda: evalues.e_process(_calls(51, effect=0.)[0], _calls(51, effect=0.)[1], seed=2,
                                                    repeats=5)
    c["sequential_e_process"] = lambda: evalues.sequential_e_process(
        [_calls(52)[0][:25], _calls(52)[0][25:45], _calls(52)[0][45:]],
        [_calls(52)[1][:25], _calls(52)[1][25:45], _calls(52)[1][45:]])
    c["e_bh"] = lambda: evalues.e_bh([40., 1., 3., 250., .5, 21.], alpha=.1)
    c["e_bh_log"] = lambda: evalues.e_bh_log(np.log([40., 1., 3., 250., .5, 21.]))
    c["anytime_bonferroni"] = lambda: evalues.anytime_bonferroni(
        [[0., 1., 3.5], [0., np.nan, 5.], [np.inf], [-1., -2.]], alpha=.05)
    c["combine"] = lambda: [evalues.combine_independent([2., 3., .5]), evalues.combine_within_cohort([2., 3., .5])]
    # missingness
    c["finite_collection_bounds"] = lambda: missingness.difference_bounds(
        missingness.finite_collection_bounds([1, 0, None, 1, np.nan, 0, 1], total_count=9),
        missingness.finite_collection_bounds([0, 0, 1, None]))
    # stats
    c["stats"] = lambda: {
        "effective_dimension": stats.effective_dimension(np.c_[np.arange(10.), np.arange(10.) ** 2,
                                                                np.sin(np.arange(10.))]),
        "fisher": [stats.fisher_exact_p(8, 2, 1, 5), stats.fisher_exact_p(0, 0, 3, 4)],
        "bh": stats.benjamini_hochberg([.001, .04, .03, .5, np.nan], q=.05, nan_policy="omit"),
        "by": stats.benjamini_hochberg([.001, .04, .03, .5], dependence="arbitrary"),
        "permutation": [stats.permutation_pvalue([1., 2., 3., 4.], 3.5),
                        stats.permutation_pvalue([1., 2., 3., 4.], 1.5, tail="less"),
                        stats.permutation_pvalue([-3., 1., 2., -1.], 2.)]}
    # censored: readings as intervals
    c["censored"] = lambda: {
        "geometry": censored.panel_geometry([.06, .125, .25, .25, .5, 1., 2., 2., 4.]),
        "geometry_wells": censored.panel_geometry([.5, 1., 1., 8.], wells=[.25, .5, 1., 2., 4., 8.]),
        "from_mic": censored.intervals_from_mic([.25, .5, 1., 2., 8.], operators=["<=", "", "", ">", ""]),
        "from_mic_uncensored_ends": censored.intervals_from_mic([.25, .5, 1., 8.],
                                                                 treat_end_wells_as_censored=False),
        "by_panel": censored.intervals_by_panel([.25, 1., 4., .5, 2.], ["p", "p", "p", "q", "q"]),
        "from_binary": censored.intervals_from_binary([0, 1, 1, 0], cutoff_log2=1.)}
    # phenotype and qc: the reading of calls and the input diagnostics
    c["phenotype_clinical"] = lambda: phenotype.to_non_susceptible(_sir(["S", "I", "R", pd.NA, "", "R"]),
                                                                   kind="clinical_sir")
    c["phenotype_undeclared"] = lambda: phenotype.to_non_susceptible(
        pd.concat([_sir(["resistant", "susceptible", "intermediate", "R", "x"]),
                   _sir(["S", "R"]).assign(Strain_ID=["000", "001"])]),
        kind="undeclared", intermediate="drop", duplicate_policy="positive_wins")
    c["phenotype_wt_nwt"] = lambda: phenotype.to_non_susceptible(_sir(["WT", "NWT", "NWT", "WT"]), kind="wt_nwt")
    c["qc"] = lambda: {"input": qc.input_qc(_qc_frames()[0], metadata=_qc_frames()[1], lineage_column="lineage"),
                       "groups": qc.group_adequacy(_qc_frames()[1]["lineage"]),
                       "traits": qc.trait_adequacy(_qc_frames()[0]),
                       "markdown": qc.render_markdown(qc.input_qc(_qc_frames()[0], metadata=_qc_frames()[1],
                                                                  lineage_column="lineage"))}
    return c


CASES = _cases()


def compute(name):
    """The result of a case, or the refusal it ends in: its type and message."""
    try:
        result = CASES[name]()
    except Exception as refusal:  # a refusal is a result too, and is pinned
        result = {"refused": type(refusal).__name__, "message": str(refusal)}
    return json.loads(json.dumps(plain(result)))


@pytest.mark.parametrize("name", sorted(CASES))
def test_every_reported_value_matches_the_committed_record(name):
    from test_golden_artefacts import compare_json
    path = VALUES / f"{name}.json"
    assert path.exists(), f"{path} is missing; run python scripts/update_golden.py"
    problems = compare_json(compute(name), json.loads(path.read_text(encoding="utf-8")))
    assert not problems, f"{name} differs from its record:\n  " + "\n  ".join(problems[:12])
