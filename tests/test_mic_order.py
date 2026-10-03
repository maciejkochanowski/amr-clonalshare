"""The lineage share of the MIC ordering: the scores, the stratified
permutation and the record."""
from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import rankdata

from amr_clonalshare.attribution import _permute_within, clonal_share, layer_clonal_share
from amr_clonalshare.latent_order import order_bounds
from amr_clonalshare.mic_order import TOLERANCE, _Scorer, _cell_ranges, _npmle, mic_order_share, order_scores


def observe_panel(values, edges):
    """Intervals (-inf, e0], (e0, e1], ..., (e_last, inf) of each value."""
    e = np.asarray(edges, dtype=float)
    y = np.asarray(values, dtype=float)
    k = np.searchsorted(e, y, side='left')
    return np.r_[-np.inf, e][k], np.r_[e, np.inf][k]


def _turnbull_reference(lo, hi, iterations=20000):
    """Self-consistency (EM) on the Turnbull cells written from the
    definition, by an explicit membership matrix; exact readings are not
    accepted here. EM approaches the maximum slowly where a cell's mass
    vanishes, so it agrees with the package to a few decimals only."""
    points = np.unique(np.r_[lo[np.isfinite(lo)], hi[np.isfinite(hi)]])
    edges = np.r_[-np.inf, points, np.inf]
    member = (edges[None, :-1] >= lo[:, None]) & (edges[None, 1:] <= hi[:, None])
    w = np.full(member.shape[1], 1 / member.shape[1])
    for _ in range(iterations):
        w = w * (member / (member @ w)[:, None]).mean(0)
    cdf = np.r_[0, np.cumsum(w)]
    below = cdf[np.searchsorted(edges, lo)]
    upto = cdf[np.searchsorted(edges, hi)]
    return (below + upto) / 2


def test_exact_readings_score_as_mid_ranks():
    rng = np.random.default_rng(1)
    y = rng.integers(0, 6, 40).astype(float)
    assert np.allclose(order_scores(y, y), (rankdata(y) - .5) / y.size, atol=1e-9)


def test_readings_of_one_panel_score_as_mid_ranks_of_their_wells():
    rng = np.random.default_rng(2)
    lo, hi = observe_panel(rng.normal(size=60), np.arange(-2., 3.))
    well = np.searchsorted(np.arange(-2., 3.), np.where(np.isfinite(hi), hi, np.inf))
    assert np.allclose(order_scores(lo, hi), (rankdata(well) - .5) / lo.size, atol=1e-9)


def test_overlapping_panels_follow_the_nonparametric_maximum_likelihood_distribution():
    rng = np.random.default_rng(3)
    y = rng.normal(size=80)
    lab = np.arange(80) % 2
    lo, hi = np.empty(80), np.empty(80)
    for k, edges in ((0, np.arange(-3., 4.)), (1, np.arange(-1.5, 2.))):
        lo[lab == k], hi[lab == k] = observe_panel(y[lab == k], edges)
    assert np.allclose(order_scores(lo, hi), _turnbull_reference(lo, hi), atol=1e-4)
    # the package's distribution meets the optimality condition of the maximum:
    # no cell's directional derivative exceeds one, and it is one on the support
    first, last, cells = _cell_ranges(lo, hi)
    cover = ((np.arange(cells) >= first[:, None]) & (np.arange(cells) <= last[:, None])).astype(float)
    weight = np.full(lo.size, 1 / lo.size)
    p = _npmle(cover, weight)
    derivative = (weight / (cover @ p)) @ cover
    assert derivative.max() <= 1 + 1e-8
    assert np.allclose(derivative[p > 1e-9], 1., atol=1e-8)


@pytest.mark.parametrize("readings", [
    # a wide and a three-well panel scored as one stratum
    ((-np.inf, -1., 14), (-4., -3., 1), (-3., -2., 17), (-2., -1., 54), (-1., 0., 169),
     (0., 1., 163), (1., 2., 12), (1., np.inf, 56), (2., 3., 1), (3., 4., 1)),
    # one panel and no reading on an end well
    ((-4., -3., 2), (-3., -2., 7), (-2., -1., 74), (-1., 0., 173), (0., 1., 134), (1., 2., 41),
     (2., 3., 11), (3., 4., 4)),
])
def test_the_distribution_reaches_its_optimality_condition(readings):
    # on these readings a Newton step that held the masses to one by a
    # weighted row, and a step halved on differences of the likelihood, each
    # stalled short of the optimality condition and raised
    lo = np.repeat([r[0] for r in readings], [r[2] for r in readings])
    hi = np.repeat([r[1] for r in readings], [r[2] for r in readings])
    first, last, cells = _cell_ranges(lo, hi)
    pairs, counts = np.unique(np.column_stack([first, last]), axis=0, return_counts=True)
    cover = ((np.arange(cells) >= pairs[:, :1]) & (np.arange(cells) <= pairs[:, 1:])).astype(float)
    weight = counts / counts.sum()
    p = _npmle(cover, weight)
    assert ((weight / (cover @ p)) @ cover).max() <= 1 + TOLERANCE
    assert np.allclose(order_scores(lo, hi), _turnbull_reference(lo, hi), atol=1e-4)


def test_exact_and_interval_readings_share_the_cells_they_both_cover():
    # two exact readings of 1 and a reading (0, 1] put their mass at 1; the
    # readings (1, 2] above them share no cell with the exact readings
    lo = np.array([1., 1., 0., 1., 1., 1.])
    hi = np.array([1., 1., 1., 2., 2., 2.])
    assert np.allclose(order_scores(lo, hi), [.25, .25, .25, .75, .75, .75], atol=1e-9)


def test_scores_are_read_within_each_stratum():
    # three laboratories on shifted panels; a reading is placed among the
    # readings of its own laboratory, not among all of them
    rng = np.random.default_rng(9)
    lab = np.repeat(np.array(["x", "y", "z"], dtype=object), 40)
    lo, hi = np.empty(120), np.empty(120)
    for name, shift in (("x", 0.), ("y", 1.), ("z", 3.)):
        idx = lab == name
        lo[idx], hi[idx] = observe_panel(rng.normal(shift, 1., 40), np.arange(-3., 4.) + shift)
    within = order_scores(lo, hi, lab)
    for name in ("x", "y", "z"):
        idx = lab == name
        assert np.allclose(within[idx], order_scores(lo[idx], hi[idx]), atol=1e-12)
    assert not np.allclose(within, order_scores(lo, hi), atol=1e-3)


def test_the_share_is_read_on_the_scores_within_strata():
    rng = np.random.default_rng(10)
    lineage = np.repeat(np.arange(12), 8)
    lab = (np.arange(12) % 3)[lineage]
    y = 1.5 * lab + .8 * rng.normal(size=12)[lineage] + rng.normal(size=lineage.size)
    lo, hi = observe_panel(y, np.arange(-3., 7.))
    lo[0], hi[0] = np.nan, 1.    # a reading with one bound missing is set aside
    kw = dict(folds=5, repeats=2, n_boot=30, n_perm=19, seed=3)
    record = mic_order_share(lo, hi, lineage, strata=lab, **kw)
    keep = np.arange(lineage.size) > 0
    labels = np.array([str(s) for s in lab[keep]], dtype=object)
    scorer = _Scorer(lo[keep], hi[keep], labels)
    scores = order_scores(lo[keep], hi[keep], labels)
    direct = layer_clonal_share(scores, lineage[keep], strata=labels, categories=scorer.category,
                                scorer=scorer, **kw)
    plain = layer_clonal_share(scores, lineage[keep], strata=labels, **kw)
    # the estimate and the test read the scores alone (centred within their
    # strata, which moves the share by rounding only, since it is invariant
    # to a shift within a stratum); the interval draws the readings again and
    # scores every draw
    assert plain.kappa_adj == pytest.approx(direct.kappa_adj, abs=1e-12)
    assert plain.p_value == direct.p_value
    assert record["n_readings"] == int(keep.sum())
    for key in ("kappa_adj", "observed_low", "observed_high", "p_value"):
        assert record[key] == pytest.approx(getattr(direct, key), abs=1e-12)
    # the bounds on the unobserved ordering from the same readings, within the
    # same strata, with the bootstrap budget and a stream of their own
    bounds = order_bounds(lo[keep], hi[keep], lineage[keep], strata=labels,
                          seed=np.random.default_rng([3, 1]))
    assert {k: record[k] for k in bounds} == bounds
    assert mic_order_share(lo, hi, lineage, stratified_by="lab", **kw)["stratified_by"] is None


def test_scores_do_not_depend_on_the_scale_of_the_readings():
    rng = np.random.default_rng(4)
    lo, hi = observe_panel(rng.normal(size=50), np.arange(-3., 4.))
    as_mg = order_scores(2. ** lo, 2. ** hi)
    assert np.allclose(order_scores(lo, hi), as_mg, atol=1e-12)


def test_one_cut_point_gives_the_share_of_the_binary_call():
    rng = np.random.default_rng(5)
    lineage = np.repeat(np.arange(12), 6)
    call = (rng.normal(size=12)[lineage] + rng.normal(size=72) > .3).astype(float)
    lo, hi = np.where(call == 1, 0., -np.inf), np.where(call == 1, np.inf, 0.)
    order = mic_order_share(lo, hi, lineage, n_boot=50, n_perm=49, repeats=4, seed=9)
    binary = clonal_share(call, lineage, n_boot=50, n_perm=49, repeats=4, seed=9)
    for key in ("kappa_adj", "observed_low", "observed_high", "p_value"):
        assert order[key] == pytest.approx(getattr(binary, key), abs=1e-10)


def test_permutation_within_strata_keeps_every_stratum_whole():
    rng = np.random.default_rng(6)
    code = rng.integers(0, 9, 200)
    strata = rng.integers(0, 3, 200)
    moved = _permute_within(code, strata, rng)
    for s in range(3):
        assert sorted(moved[strata == s]) == sorted(code[strata == s])
        assert not np.array_equal(moved[strata == s], code[strata == s])


def test_laboratory_offset_aligned_with_lineage_is_not_read_as_a_lineage_effect():
    # every lineage read in one of two laboratories, which differ by two wells
    # and in panel; there is no lineage effect within a laboratory
    rng = np.random.default_rng(7)
    lineage = np.repeat(np.arange(16), 8)
    lab = (np.arange(16) % 2)[lineage]
    y = 2. * lab + rng.normal(size=lineage.size)
    lo, hi = np.empty(y.size), np.empty(y.size)
    for k, edges in ((0, np.arange(-3., 4.)), (1, np.arange(-1., 6.))):
        lo[lab == k], hi[lab == k] = observe_panel(y[lab == k], edges)
    pooled = mic_order_share(lo, hi, lineage, n_perm=199, n_boot=50, repeats=4, seed=1)
    within = mic_order_share(lo, hi, lineage, strata=lab, stratified_by="laboratory",
                             n_perm=199, n_boot=50, repeats=4, seed=1)
    assert pooled["p_value"] <= .01
    assert within["p_value"] > .05 and within["n_strata"] == 2


def test_the_record_says_what_the_share_rests_on():
    rng = np.random.default_rng(8)
    lineage = np.repeat(np.array(["a", "b", "c", "d", "e"], dtype=object), 6)
    lo, hi = observe_panel(rng.normal(size=30), np.array([-1., 0., 1.]))
    lab = np.where(np.arange(30) < 20, "x", "y")
    lo[0] = hi[0] = np.nan
    r = mic_order_share(lo, hi, lineage, strata=lab, stratified_by="lab", n_boot=20,
                        n_perm=19, repeats=2)
    assert r["n_readings"] == 29 and r["strata"] == {"x": 19, "y": 10}
    assert r["stratified_by"] == "lab"
    assert r["share_end_wells"] == pytest.approx(np.mean(np.isinf(lo[1:]) | np.isinf(hi[1:])))
    # the probit restatement of a call has no place here; the bounds on the
    # share of the unobserved ordering do
    assert not any(k.startswith("latent_") and not k.startswith("latent_order_") for k in r)
    assert "prevalence" not in r and "latent_order_lower" in r
    assert mic_order_share(lo, hi, lineage, n_boot=20, n_perm=19, repeats=2)["stratified_by"] is None


@pytest.mark.parametrize("lo,hi,message", [
    ([0., 1.], [1.], "one length"),
    ([np.nan, 0.], [1., 1.], "without a value"),
    ([-np.inf, 0.], [np.inf, 1.], "unbounded"),
    ([2., 0.], [1., 1.], "lo <= hi"),
])
def test_scores_refuse_readings_they_cannot_place(lo, hi, message):
    with pytest.raises(ValueError, match=message):
        order_scores(np.array(lo), np.array(hi))


def test_strata_are_checked_against_the_readings():
    with pytest.raises(ValueError, match="one label per reading"):
        order_scores(np.zeros(3), np.ones(3), strata=["a", "b"])
    lineage = np.repeat(np.arange(4), 3)
    x = np.arange(12.)
    same = layer_clonal_share(x, lineage, n_boot=0, n_perm=9, repeats=2, strata=np.zeros(12))
    assert same.n_strata == 1


def test_the_record_bounds_the_unobserved_ordering_within_the_same_strata():
    """The bounds come from the same readings, within the same strata, with
    the bootstrap budget of the share and a stream of their own."""
    rng = np.random.default_rng(11)
    lineage = np.repeat(np.arange(10), 8)
    lab = np.tile(["x", "y"], 40)
    y = 2. * rng.normal(size=10)[lineage] + 1.5 * (lab == "y") + rng.normal(size=80)
    lo, hi = observe_panel(y, np.arange(-3., 4.))
    record = mic_order_share(lo, hi, lineage, strata=lab, n_boot=30, n_perm=19, repeats=2, seed=4)
    labels = np.array([str(s) for s in lab], dtype=object)
    bounds = order_bounds(lo, hi, lineage, strata=labels, seed=np.random.default_rng([4, 1]))
    assert record["latent_order_lower_limit"] > 0.
    assert {k: record[k] for k in bounds} == bounds
