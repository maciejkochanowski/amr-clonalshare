"""The bounds on the lineage share of the unobserved MIC ordering.

The lower end is checked against a quadratic programme over the set of
lineage means written out subset by subset, which shares no code with the
minimum-norm-point algorithm; the upper end against every order of the
lineages, and its certified bound against the same; the characterisation
against a direct model of the allocations inside every reading; the dual
certificate against the lower end it bounds; the split limit against its
definition; and the identities, the nesting and the invariances the module
states against their definitions.
"""
from __future__ import annotations

import itertools

import numpy as np
import pytest
from scipy.optimize import linprog, minimize

from amr_clonalshare import latent_order as lo_mod
from amr_clonalshare.latent_order import (Readings, _Polytope, _certificate, _lower, _minimum_norm, _split,
                                          _upper, _upper_certificate, _upper_exact, call_bounds, order_bounds)
from amr_clonalshare.mic_order import order_scores


def _instance(rng, G, B):
    """Lineage-by-reading counts of a location model read on B readings."""
    shift = rng.normal(0, 1, G)
    cuts = np.sort(rng.normal(0, 1, B - 1))
    counts = np.zeros((G, B))
    for g in range(G):
        y = rng.normal(shift[g], 1, int(rng.integers(3, 25)))
        counts[g] = np.bincount(np.searchsorted(cuts, y), minlength=B)
    counts = counts[:, counts.sum(0) > 0]
    share = counts.sum(0) / counts.sum()
    cum = np.r_[0., np.cumsum(share)]
    return _Polytope(counts, cum[:-1], cum[1:])


def _subset_programme(poly):
    """min sum x^2 / w with y_gr constrained by f_r(A) for every set A."""
    G, B = poly.p.shape
    cons = []
    for r in range(B):
        present = [g for g in range(G) if poly.pi[g, r] > 0]
        for size in range(1, len(present) + 1):
            for A in itertools.combinations(present, size):
                s = sum(poly.pi[g, r] for g in A)
                f = s * poly.H[r] - poly.D[r] * s * s / 2
                idx = [g * B + r for g in A]
                if size == len(present):
                    cons.append(dict(type="eq", fun=lambda v, idx=idx, f=f: np.sum(v[idx]) - f))
                else:
                    cons.append(dict(type="ineq", fun=lambda v, idx=idx, f=f: f - np.sum(v[idx])))
        for g in range(G):
            if poly.pi[g, r] == 0:
                cons.append(dict(type="eq", fun=lambda v, i=g * B + r: v[i]))

    def objective(v):
        x = v.reshape(G, B) @ poly.Pr
        return np.sum(x ** 2 / poly.w)

    start = (poly.pi * ((poly.L + poly.H) / 2)).ravel()
    res = minimize(objective, start, constraints=cons, method="SLSQP",
                   options=dict(ftol=1e-14, maxiter=2000))
    return 12 * (res.fun - 0.25)


@pytest.mark.parametrize("seed", range(6))
def test_the_lower_end_is_the_minimum_over_the_set_written_subset_by_subset(seed):
    rng = np.random.default_rng(seed)
    poly = _instance(rng, int(rng.integers(3, 5)), int(rng.integers(2, 4)))
    value, gap = _lower(poly)
    assert gap <= 1e-9
    assert value == pytest.approx(_subset_programme(poly), abs=1e-7)


def _every_order(poly):
    return max(poly.rho(poly.vertex(o)) for o in itertools.permutations(range(poly.w.size)))


@pytest.mark.parametrize("seed", range(8))
def test_the_upper_end_is_the_largest_value_over_every_order(seed):
    rng = np.random.default_rng(100 + seed)
    poly = _instance(rng, int(rng.integers(3, 7)), int(rng.integers(2, 6)))
    every = _every_order(poly)
    attained, bound, exact = _upper(poly, 10, np.random.default_rng(0))
    assert exact and attained == bound
    assert attained == pytest.approx(every, abs=1e-12)
    assert 12 * (_upper_exact(poly) - .25) == pytest.approx(every, abs=1e-12)
    low, _ = _lower(poly)
    assert low <= poly.rho(poly.midpoint()) + 1e-12 <= every + 2e-12


@pytest.mark.parametrize("seed", range(6))
def test_the_certified_upper_bound_holds_and_closes_when_the_search_may_finish(monkeypatch, seed):
    """Past EXACT_LINEAGES the upper end is an insertion search with a
    branch-and-bound certificate; forced onto small instances, the attained
    share never exceeds the largest over every order and the bound never falls
    below it, and with enough nodes the two meet."""
    rng = np.random.default_rng(200 + seed)
    poly = _instance(rng, int(rng.integers(4, 7)), int(rng.integers(2, 6)))
    every = _every_order(poly)
    monkeypatch.setattr(lo_mod, "EXACT_LINEAGES", 0)
    monkeypatch.setattr(lo_mod, "DP_FIRST", 0)
    for nodes in (1, 3, 10 ** 6):
        monkeypatch.setattr(lo_mod, "UPPER_NODES", nodes)
        attained, bound, exact = _upper(poly, 2, np.random.default_rng(1))
        assert attained <= every + 1e-12 <= bound + 2e-12
        assert exact == (bound - attained <= 12e-12)
    assert exact and bound == pytest.approx(every, abs=1e-11)


def test_a_gap_left_by_the_search_is_closed_by_the_dynamic_programme(monkeypatch):
    """Up to EXACT_LINEAGES lineages the upper end is exact whatever the
    branch-and-bound manages: with it cut to one bound, the dynamic
    programme is run."""
    poly = _instance(np.random.default_rng(31), 14, 5)
    monkeypatch.setattr(lo_mod, "UPPER_NODES", 1)
    attained, bound, exact = _upper(poly, 1, np.random.default_rng(0))
    assert exact and attained == bound == pytest.approx(12 * (_upper_exact(poly) - .25), abs=1e-12)


def test_interchangeable_lineages_let_the_certificate_close_on_many_calls():
    """Calls over forty lineages of five isolates fall into six kinds of
    lineage; placing only one of each kind next, the branch-and-bound closes
    on the share the search attains."""
    rng = np.random.default_rng(0)
    code = np.repeat(np.arange(40), 5)
    y = (rng.random(code.size) < np.r_[[.02] * 39, [.5]][code]).astype(float)
    reads = Readings(np.where(y == 1, 0., -np.inf), np.where(y == 1, np.inf, 0.), code,
                     np.zeros(code.size, dtype=object))
    attained, bound, exact = _upper(reads.polytope(), 10, np.random.default_rng(0))
    assert exact and attained == bound


@pytest.mark.parametrize("G, B, seed", [(2, 3, 1), (5, 4, 2), (9, 5, 3), (14, 6, 4)])
def test_the_insertion_search_returns_a_locally_best_order_and_its_share(G, B, seed):
    """The search returns an order and the sum at its vertex; no lineage moved
    to another place in that order raises the sum, the order is no worse than
    the first start (lineages by midpoint mean) and no better than the
    largest over every order."""
    from amr_clonalshare.latent_order import _upper_search
    poly = _instance(np.random.default_rng(300 + seed), G, B)
    G = poly.w.size
    value, order = _upper_search(poly, 3, np.random.default_rng(seed))
    total = lambda o: float(np.sum(poly.vertex(o) ** 2 / poly.w))
    assert sorted(order.tolist()) == list(range(G))
    assert value == pytest.approx(total(order), rel=1e-12)
    for i in range(G):
        rest = np.delete(order, i)
        for j in range(G):
            assert total(np.insert(rest, j, order[i])) <= value + 1e-12
    assert value >= total(np.argsort(-poly.midpoint() / poly.w, kind="stable")) - 1e-12
    assert value <= _upper_exact(poly) + 1e-12


@pytest.mark.parametrize("seed", range(4))
def test_the_certificate_from_no_incumbent_reaches_the_maximum_and_certifies_it(seed):
    """With no share attained beforehand, the branch and bound builds orders
    down to their last lineage: the share it attains and the bound it
    certifies are both the largest over every order."""
    rng = np.random.default_rng(400 + seed)
    poly = _instance(rng, int(rng.integers(3, 7)), int(rng.integers(2, 6)))
    top = _upper_exact(poly)
    bound, attained = _upper_certificate(poly, 0.)
    assert attained == pytest.approx(top, abs=1e-12)
    assert bound == pytest.approx(top, abs=1e-12)


def test_the_certificate_of_an_incumbent_above_the_maximum_returns_it():
    """A node is opened only if its bound beats the incumbent, so an incumbent
    already at the maximum closes every node."""
    poly = _instance(np.random.default_rng(3), 6, 4)
    top = _upper_exact(poly)
    bound, attained = _upper_certificate(poly, top)
    # The sum of an order is accumulated in a different sequence from the
    # dynamic programme's, so the two agree to round-off, not bit for bit.
    assert bound == pytest.approx(top, abs=1e-12)
    assert attained == pytest.approx(top, abs=1e-12)


def test_the_set_of_lineage_means_is_the_set_of_allocations_inside_the_readings():
    """A direct model: every reading's range cut into M slices, each lineage
    given a mass in every slice, the slices of a reading filled uniformly.
    Its lineage means form the polytope up to the slice width, so the linear
    minimum along any direction agrees with the greedy vertex."""
    rng = np.random.default_rng(7)
    poly = _instance(rng, 4, 3)
    G, B = poly.p.shape
    M = 400
    for _ in range(5):
        c = rng.normal(size=G)
        greedy = c @ poly.vertex(np.argsort(c))          # min of c.x over the polytope
        # variables q[g, r, j] >= 0: mass of lineage g in slice j of reading r
        n = G * B * M
        cost = np.zeros(n)
        A_eq, b_eq = [], []
        for r in range(B):
            u = poly.L[r] + poly.D[r] * (np.arange(M) + .5) / M
            for g in range(G):
                cost[(g * B + r) * M:(g * B + r + 1) * M] = c[g] * poly.Pr[r] * u
                row = np.zeros(n)
                row[(g * B + r) * M:(g * B + r + 1) * M] = 1.
                A_eq.append(row)
                b_eq.append(poly.pi[g, r])
            for j in range(M):
                row = np.zeros(n)
                row[[(g * B + r) * M + j for g in range(G)]] = 1.
                A_eq.append(row)
                b_eq.append(1. / M)
        res = linprog(cost, A_eq=np.asarray(A_eq), b_eq=b_eq, bounds=(0, None), method="highs")
        assert res.status == 0
        assert res.fun == pytest.approx(greedy, abs=np.abs(c).sum() * np.max(poly.D) / M)


def test_the_midpoint_is_the_resolution_times_the_share_of_the_scores():
    rng = np.random.default_rng(3)
    code = np.repeat(np.arange(9), rng.integers(2, 12, 9))
    y = rng.normal(0, 1, code.size) + rng.normal(0, .8, 9)[code]
    lab = rng.integers(0, 2, code.size)
    cuts = [np.arange(-2., 3.), np.arange(-1., 2.)]
    lo = np.empty(code.size)
    hi = np.empty(code.size)
    for s in (0, 1):
        m = lab == s
        k = np.searchsorted(cuts[s], y[m])
        ext = np.r_[-np.inf, cuts[s], np.inf]
        lo[m], hi[m] = ext[k], ext[k + 1]
    r = order_bounds(lo, hi, code, strata=lab)
    w = order_scores(lo, hi, lab)
    sizes = np.bincount(code)
    between = np.sum(sizes / code.size * (np.bincount(code, weights=w) / sizes - w.mean()) ** 2)
    # the scores come from the maximum-likelihood distribution of the readings,
    # which is computed to mic_order.TOLERANCE, the bounds from exact shares
    assert r["latent_order_midpoint"] == pytest.approx(12 * between, rel=1e-8)
    assert r["latent_order_midpoint"] == pytest.approx(r["panel_resolution"] * between / w.var(), rel=1e-8)
    assert r["latent_order_lower"] <= r["latent_order_midpoint"] <= r["latent_order_upper"]


def _feasible_shares(lo, hi, lineage):
    """Every order of the isolates the readings allow, and the lineage share of
    its mid-ranks: an isolate whose reading lies wholly below another's must
    come first; any two isolates whose readings overlap may stand either way."""
    n = lo.size
    must = hi[:, None] <= lo[None, :]
    names = sorted(set(lineage))
    shares = []
    for order in itertools.permutations(range(n)):
        rank = np.empty(n, dtype=int)
        rank[list(order)] = np.arange(n)
        if np.any(must & (rank[:, None] > rank[None, :])):
            continue
        u = (rank + .5) / n
        shares.append(12. * sum(np.mean(lineage == g) * (np.mean(u[lineage == g]) - .5) ** 2 for g in names))
    return np.asarray(shares)


def test_overlapping_readings_in_one_stratum_give_bounds_that_hold_for_every_allowed_order():
    """The example of the review: readings on one range whose intervals
    overlap. Every order the readings allow, enumerated, lies within the
    reported bounds, and the certified upper bound is above the largest."""
    from amr_clonalshare.censored import intervals_from_mic
    values = np.array([4., 16., 8., 4., 16., 8.])
    operators = np.array([">", "<=", "<=", "<=", "=", ">"])
    lineage = np.array(["A", "A", "A", "B", "B", "B"])
    lo, hi = intervals_from_mic(values, operators=operators, wells=[1, 2, 4, 8, 16])
    shares = _feasible_shares(lo, hi, lineage)
    assert shares.size == 84
    out = order_bounds(lo, hi, lineage, seed=0)
    assert not out["latent_order_sharp"]
    assert out["latent_order_lower"] <= shares.min() + 1e-12
    assert out["latent_order_upper_bound"] >= shares.max() - 1e-12
    assert out["latent_order_lower_limit"] <= out["latent_order_lower"]


@pytest.mark.parametrize("seed", range(40))
def test_the_bounds_contain_every_allowed_order_whether_or_not_the_readings_overlap(seed):
    """Small random readings, some on one range and some on two ranges mixed
    in one stratum: the bounds hold for every order the readings allow,
    found by enumeration, and are reported sharp only when the distinct
    readings are disjoint."""
    rng = np.random.default_rng(seed)
    n = int(rng.integers(4, 8))
    lineage = np.array(["A", "B", "C"])[rng.integers(0, 2 + (seed % 2), n)]
    while min(np.sum(lineage == g) for g in set(lineage)) < 2 or len(set(lineage)) < 2:
        lineage = np.array(["A", "B", "C"])[rng.integers(0, 2 + (seed % 2), n)]
    y = rng.normal(size=n)
    cuts_a = np.array([-1., 0., 1.])
    cuts_b = np.array([-0.5, 0.5]) if seed % 3 else cuts_a
    use_a = rng.random(n) < .6
    lo, hi = np.empty(n), np.empty(n)
    for m, cuts in ((use_a, cuts_a), (~use_a, cuts_b)):
        if m.any():
            lo[m], hi[m] = _panel(y[m], cuts)
    shares = _feasible_shares(lo, hi, lineage)
    out = order_bounds(lo, hi, lineage, seed=seed)
    assert out["latent_order_lower"] <= shares.min() + 1e-9
    assert out["latent_order_upper_bound"] >= shares.max() - 1e-9
    assert out["latent_order_lower_limit"] <= out["latent_order_lower"] + 1e-12
    pairs = sorted(set(zip(lo, hi)))
    disjoint = all(pairs[k][1] <= pairs[k + 1][0] for k in range(len(pairs) - 1))
    assert out["latent_order_sharp"] == disjoint


def test_overlapping_readings_are_merged_into_ordered_cells_that_take_their_share_of_the_stratum():
    """With two panels in one stratum the readings overlap; each run of
    overlapping readings is one cell, the cells are ordered, and their ranges
    partition (0, 1) in proportion to the isolates they hold."""
    rng = np.random.default_rng(4)
    y = rng.normal(size=300)
    code = rng.integers(0, 12, 300)
    a, b = np.arange(-2., 3.), np.array([-1., 1.])
    use_a = rng.random(300) < .6
    lo = np.empty(300)
    hi = np.empty(300)
    for m, cuts in ((use_a, a), (~use_a, b)):
        k = np.searchsorted(cuts, y[m])
        ext = np.r_[-np.inf, cuts, np.inf]
        lo[m], hi[m] = ext[k], ext[k + 1]
    reads = Readings(lo, hi, code, np.zeros(300, dtype=object))
    assert not reads.groups[0][2] and not reads.sharp     # the readings overlap
    assert reads.n_cells < reads.n_bins
    counts = np.bincount(reads.cell[reads.bins], minlength=reads.n_cells).astype(float)
    low, high = reads.ranges(counts)
    assert low[0] == 0. and high[-1] == pytest.approx(1.)
    assert np.all(high[:-1] == pytest.approx(low[1:]))
    assert np.all(high - low == pytest.approx(counts / counts.sum()))


def _panel(y, cuts):
    k = np.searchsorted(cuts, y)
    ext = np.r_[-np.inf, cuts, np.inf]
    return ext[k], ext[k + 1]


def test_a_coarser_panel_and_a_call_allow_a_wider_interval():
    rng = np.random.default_rng(5)
    code = np.repeat(np.arange(7), rng.integers(3, 15, 7))
    y = rng.normal(0, 1, code.size) + rng.normal(0, 1, 7)[code]
    fine = order_bounds(*_panel(y, np.arange(-3., 3.5, .5)), code)
    coarse = order_bounds(*_panel(y, np.arange(-3., 4., 1.)), code)
    call = call_bounds((y > 0.).astype(float), code)
    for wide, narrow in ((coarse, fine), (call, coarse)):
        assert wide["latent_order_lower"] <= narrow["latent_order_lower"] + 1e-12
        assert wide["latent_order_upper"] >= narrow["latent_order_upper"] - 1e-9


def test_the_bounds_do_not_depend_on_how_the_end_wells_or_the_scale_are_read():
    rng = np.random.default_rng(6)
    code = np.repeat(np.arange(6), rng.integers(3, 10, 6))
    y = rng.normal(0, 1, code.size) + rng.normal(0, 1, 6)[code]
    lo, hi = _panel(y, np.arange(-2., 3.))
    base = order_bounds(lo, hi, code)
    exact_ends = order_bounds(np.where(np.isneginf(lo), hi - 1., lo), np.where(np.isposinf(hi), lo + 1., hi),
                              code)
    rescaled = order_bounds(3. * lo + 7., 3. * hi + 7., code)
    perm = rng.permutation(code.size)
    relabelled = order_bounds(lo[perm], hi[perm], np.asarray([f"L{9 - c}" for c in code[perm]]))
    for other in (exact_ends, rescaled, relabelled):
        for key in ("latent_order_lower", "latent_order_upper", "latent_order_upper_bound",
                    "latent_order_midpoint", "panel_resolution"):
            assert other[key] == pytest.approx(base[key], abs=1e-12)


def test_lineages_alone_in_their_readings_fix_the_share_and_alike_lineages_allow_zero():
    lo = np.repeat([-np.inf, 0., 1.], 4)
    hi = np.repeat([0., 1., np.inf], 4)
    alone = order_bounds(lo, hi, np.repeat([0, 1, 2], 4))
    # the lower end is certified to the tolerance of the minimum-norm point
    tolerance = 2 * lo_mod.GAP_TOLERANCE
    assert alone["latent_order_lower"] == pytest.approx(alone["latent_order_upper"], abs=tolerance)
    assert alone["latent_order_lower"] == pytest.approx(alone["latent_order_midpoint"], abs=tolerance)
    assert alone["latent_order_lower"] <= alone["latent_order_upper"]
    alike = order_bounds(lo, hi, np.tile([0, 1, 2, 3], 3))
    assert alike["latent_order_lower"] == pytest.approx(0., abs=1e-12)


def test_the_lower_limit_lies_below_the_lower_end_and_nothing_is_reported_without_two_lineages():
    rng = np.random.default_rng(8)
    code = np.repeat(np.arange(8), 10)
    y = rng.normal(0, 1, 80) + rng.normal(0, 1, 8)[code]
    r = order_bounds(*_panel(y, np.arange(-2., 3.)), code, seed=1)
    assert 0. < r["latent_order_lower_limit"] <= r["latent_order_lower"]
    one = order_bounds(*_panel(y[:10], np.arange(-2., 3.)), code[:10])
    assert np.isnan(one["latent_order_lower"]) and one["latent_order_readings"] == 0
    with pytest.raises(ValueError):
        call_bounds(np.array([0., 2.]), np.array([0, 0]))


@pytest.mark.parametrize("cap", [1, 2, 5])
def test_a_run_cut_short_still_returns_a_lower_end_and_a_gap_that_brackets_the_minimum(monkeypatch, cap):
    rng = np.random.default_rng(9)
    poly = _instance(rng, 12, 6)
    full, full_gap = _lower(poly)
    assert full_gap <= lo_mod.GAP_TOLERANCE
    monkeypatch.setattr(lo_mod, "MAX_MAJOR", cap)
    value, gap = _lower(poly)
    # the minimum lies in [full, full + tolerance] and in [value, value + gap]
    assert value <= full + lo_mod.GAP_TOLERANCE
    assert value + max(gap, lo_mod.GAP_TOLERANCE) >= full - 1e-12


def test_lineages_of_one_profile_are_merged_without_changing_the_minimum():
    """Two pairs of lineages spread over the readings in the same proportions
    (one pair at twice the size); the lower end is the minimum over the set
    written subset by subset for the lineages as they are."""
    rng = np.random.default_rng(11)
    counts = rng.integers(1, 6, (3, 3)).astype(float)
    counts = np.vstack([counts, 2 * counts[0], counts[1]])
    share = counts.sum(0) / counts.sum()
    cum = np.r_[0., np.cumsum(share)]
    poly = _Polytope(counts, cum[:-1], cum[1:])
    merged, group = poly.merged()
    assert merged.w.size == 3 and group[0] == group[3] and group[1] == group[4]
    value, gap = _lower(poly)
    assert gap <= lo_mod.GAP_TOLERANCE
    assert value == pytest.approx(_subset_programme(poly), abs=1e-7)


def test_a_call_over_many_small_lineages_does_not_depend_on_splitting_alike_lineages():
    """Many small lineages read on one cut point share a few rates. Every
    lineage split into two of its own rate, one of them twice the size of the
    other, leaves the reading shares and the lower end as they were."""
    rng = np.random.default_rng(12)
    code = np.repeat(np.arange(400), 4)
    y = (rng.random(code.size) < rng.beta(.6, .9, 400)[code]).astype(float)
    whole = call_bounds(y, code)
    split = call_bounds(np.r_[y, y, y], np.r_[2 * code, 2 * code + 1, 2 * code + 1])
    assert whole["latent_order_gap"] <= lo_mod.GAP_TOLERANCE
    assert split["latent_order_lower"] == pytest.approx(whole["latent_order_lower"], abs=2 * lo_mod.GAP_TOLERANCE)


def test_a_split_halves_every_lineage_in_every_stratum():
    """Every cell of a lineage and a stratum is dealt at random between the
    two halves, as evenly as its size allows, and the halves partition the
    isolates."""
    rng = np.random.default_rng(13)
    code = np.repeat(np.arange(5), [2, 3, 4, 6, 9])
    lab = rng.integers(0, 2, code.size)
    y = rng.normal(0, 1, code.size) + rng.normal(0, 1, 5)[code]
    reads = Readings(*_panel(y, np.arange(-2., 3.)), code, lab.astype(object))
    seen = np.zeros(code.size, dtype=int)
    for r in range(50):
        first, second = _split(reads, np.random.default_rng(r))
        assert np.array_equal(np.sort(np.r_[first, second]), np.arange(code.size))
        for g in range(5):
            for s in (0, 1):
                size = int(np.sum((code == g) & (lab == s)))
                k = int(np.sum((code[first] == g) & (lab[first] == s)))
                assert k in (size // 2, (size + 1) // 2)
        seen[first] += 1
    assert (seen > 0).all() and (seen < 50).all()


def _dual(reads, z):
    """G(z) of the readings in hand, taken as the population, with Q(z)."""
    n_cell = np.zeros((reads.G, len(reads.names)))
    np.add.at(n_cell, (reads.code, reads.stratum), 1.)
    w = n_cell.sum(1) / n_cell.sum()
    z = z - np.sum(w * z)
    full = reads.cell_counts(np.arange(reads.code.size))
    value, _, _ = _certificate(reads, z, full, full, n_cell)
    return value, float(np.sum(w * z * z))


@pytest.mark.parametrize("seed", range(5))
def test_the_dual_certificate_bounds_the_lower_end_and_meets_it_at_the_minimum_norm_point(seed):
    """Theorem 2 on the readings in hand: 12 G(z)_+^2 / Q(z) never exceeds the
    guaranteed share, and at z = E[U | g] of the minimum-norm point it equals
    it, one panel in every stratum."""
    rng = np.random.default_rng(300 + seed)
    G = int(rng.integers(3, 8))
    code = np.repeat(np.arange(G), rng.integers(2, 12, G))
    lab = rng.integers(0, 2, code.size)
    y = rng.normal(0, 1, code.size) + rng.normal(0, 1, G)[code] + .5 * lab
    lo, hi = np.empty(code.size), np.empty(code.size)
    for s, cuts in ((0, np.arange(-2., 3.)), (1, np.array([-.5, .5, 1.5]))):
        m = lab == s
        lo[m], hi[m] = _panel(y[m], cuts)
    reads = Readings(lo, hi, code, lab.astype(object))
    lower, gap, means = _minimum_norm(reads.polytope())
    assert gap <= lo_mod.GAP_TOLERANCE
    value, Q = _dual(reads, means - .5)
    assert 12 * max(value, 0.) ** 2 / Q == pytest.approx(lower, abs=1e-8)
    for _ in range(20):
        value, Q = _dual(reads, rng.normal(size=G))
        assert 12 * max(value, 0.) ** 2 / Q <= lower + 1e-9


def test_the_estimate_of_the_certificate_is_the_average_over_pairs_of_isolates():
    """For a fixed z the estimate from one half is the U-statistic: over the
    pairs of isolates of two lineages in one stratum, the weight times whether
    the reading of the lower-scored lineage lies wholly below the other, less
    one half."""
    rng = np.random.default_rng(21)
    code = np.repeat(np.arange(4), [3, 4, 5, 6])
    lab = rng.integers(0, 2, code.size)
    y = rng.normal(0, 1, code.size) + rng.normal(0, 1, 4)[code]
    lo, hi = _panel(y, np.arange(-1., 2.))
    reads = Readings(lo, hi, code, lab.astype(object))
    n_cell = np.zeros((4, 2))
    np.add.at(n_cell, (code, lab), 1.)
    z = np.array([.3, -.2, .1, -.4])
    rows = np.flatnonzero(rng.random(code.size) < .6)
    value, _, _ = _certificate(reads, z, reads.cell_counts(rows), reads.cell_counts(np.arange(code.size)), n_cell)
    oracle = 0.
    for s in (0, 1):
        ns = n_cell[:, s].sum()
        for g in range(4):
            for h in range(4):
                if z[g] <= z[h] or n_cell[g, s] == 0 or n_cell[h, s] == 0:
                    continue
                weight = n_cell[g, s] * n_cell[h, s] / (code.size * ns) * (z[g] - z[h])
                a = [i for i in rows if code[i] == g and lab[i] == s]
                b = [j for j in rows if code[j] == h and lab[j] == s]
                if not a or not b:
                    oracle -= .5 * weight            # the worst case, P = 0
                    continue
                below = np.mean([[hi[j] <= lo[i] for i in a] for j in b])
                oracle += weight * (below - .5)
    assert value == pytest.approx(oracle, abs=1e-12)


def test_the_variance_of_the_certificate_is_the_variance_of_its_estimate():
    """Over datasets drawn from one law, with the counts of every lineage
    fixed, the variance the certificate reports for its estimate at a fixed
    z is, on average, the variance of that estimate across the datasets: the
    projection of every reading and the pair terms together, neither left
    out nor counted twice."""
    rng = np.random.default_rng(31)
    sizes = np.array([20, 25, 30, 25])
    code = np.repeat(np.arange(4), sizes)
    laws = np.array([[.4, .3, .2, .1, 0.], [.2, .3, .3, .1, .1], [.1, .2, .3, .2, .2], [0., .1, .2, .3, .4]])
    z = np.array([-.3, -.1, .1, .3])
    n_cell = sizes.astype(float)[:, None]
    estimates, variances = [], []
    for _ in range(600):
        reading = np.concatenate([rng.choice(5, m, p=laws[g]) for g, m in enumerate(sizes)])
        lo, hi = reading - 1., reading.astype(float)
        lo[reading == 0] = -np.inf
        hi[reading == 4] = np.inf
        reads = Readings(lo, hi, code, np.zeros(code.size, dtype=object))
        everyone = reads.cell_counts(np.arange(code.size))
        value, variance, _ = _certificate(reads, z, everyone, everyone, n_cell)
        estimates.append(value)
        variances.append(variance)
    ratio = np.mean(variances) / np.var(estimates, ddof=1)
    assert .8 < ratio < 1.35, ratio


def test_a_lineage_absent_from_the_evaluation_half_adds_nothing_to_the_variance():
    """A pair with a lineage that has no isolate in the evaluation half enters
    the estimate at its worst case, a constant, so the readings of that lineage
    in the whole sample move neither the estimate nor its variance."""
    rng = np.random.default_rng(7)
    sizes = np.array([12, 15, 10, 8])
    code = np.repeat(np.arange(4), sizes)
    z = np.array([-.2, .3, .1, -.4])
    n_cell = sizes.astype(float)[:, None]
    base = np.concatenate([np.r_[np.arange(5), rng.integers(0, 5, m - 5)] for m in sizes[:3]])
    rows = np.flatnonzero(code < 3)
    results = []
    for last in (np.zeros(sizes[3], dtype=int), np.full(sizes[3], 4)):
        reading = np.r_[base, last]
        lo, hi = reading - 1., reading.astype(float)
        lo[reading == 0] = -np.inf
        hi[reading == 4] = np.inf
        reads = Readings(lo, hi, code, np.zeros(code.size, dtype=object))
        value, variance, _ = _certificate(reads, z, reads.cell_counts(rows),
                                          reads.cell_counts(np.arange(code.size)), n_cell)
        results.append((value, variance))
    assert results[0][0] == pytest.approx(results[1][0], abs=1e-12)
    assert results[0][1] == pytest.approx(results[1][1], rel=1e-12)
    assert results[0][1] > 0


def _certificate_by_hand(lo, hi, code, stratum, part_rows, z):
    """The estimate of G(z), its variance and degrees of freedom written out
    isolate by isolate from the module docstring: pairs of isolates for the
    estimate; the variance of the first-order projection of a reading of each
    cell under the cell's law smoothed by Perks' prior over the readings of
    its stratum, over the isolates of the cell in the evaluation half; the pair
    terms bounded by theta (1 - theta) with theta the smoothed share of pairs of
    all isolates; Satterthwaite's degrees of freedom, a cell's term carrying
    the number of its isolates and the pair terms infinitely many."""
    n = code.size
    part = set(part_rows.tolist())
    estimate, terms, dfs = 0., [], []
    for s in np.unique(stratum):
        rows = np.flatnonzero(stratum == s)
        ns = rows.size
        lineages = np.unique(code[rows])
        n_gs = {g: np.sum(code[rows] == g) for g in lineages}
        full = {g: [i for i in rows if code[i] == g] for g in lineages}
        half = {g: [i for i in full[g] if i in part] for g in lineages}
        readings = sorted({(lo[i], hi[i]) for i in rows})
        K = len(readings)

        def A(g, h):
            return n_gs[g] * n_gs[h] / (n * ns) * (z[g] - z[h]) if z[g] > z[h] else 0.

        def below(j, i):                                    # reading j wholly below reading i
            return float(hi[j] <= lo[i])

        both = [(g, h) for g in lineages for h in lineages if A(g, h) > 0 and half[g] and half[h]]
        for g in lineages:
            for h in lineages:
                if A(g, h) == 0.:
                    continue
                if (g, h) in both:
                    p = np.mean([below(j, i) for i in half[g] for j in half[h]])
                    estimate += A(g, h) * (p - .5)
                else:
                    estimate -= .5 * A(g, h)
        for g in lineages:
            if not half[g]:
                continue
            size = len(full[g])
            f, law = [], []
            for (a, b) in readings:
                value = 0.
                for (x, y) in both:
                    if x == g:          # g above y: a reading of y below this one
                        value += A(g, y) * np.mean([float(hi[j] <= a) for j in full[y]])
                    if y == g:          # x above g: this reading below one of x
                        value += A(x, g) * np.mean([float(b <= lo[i]) for i in full[x]])
                f.append(value)
                law.append((sum(1 for i in full[g] if (lo[i], hi[i]) == (a, b)) + 1. / K) / (size + 1.))
            f, law = np.array(f), np.array(law)
            spread = np.sum(law * (f - np.sum(law * f)) ** 2) * size / max(size - 1., 1.)
            terms.append(spread / len(half[g]))
            dfs.append(size)
        second = 0.
        for (g, h) in both:
            pairs = len(full[g]) * len(full[h])
            theta = (sum(below(j, i) for i in full[g] for j in full[h]) + .5) / (pairs + 1.)
            second += A(g, h) ** 2 * theta * (1. - theta) / (len(half[g]) * len(half[h]))
        if second > 0:
            terms.append(second)
            dfs.append(np.inf)
    variance = float(np.sum(terms))
    spread_df = float(np.sum(np.square(terms) / np.array(dfs, dtype=float))) if terms else 0.
    return estimate, variance, (variance ** 2 / spread_df if spread_df > 0 else np.inf)


def test_the_certificate_is_its_written_form():
    """Two strata on different panels; a lineage with one isolate in a
    stratum, a cell whose isolates share one reading, lineages with one
    isolate in the evaluation half and lineages with none there."""
    #            lineage: 0  1  2  3  4
    sizes = {"a": [6, 5, 1, 4, 3], "b": [3, 0, 2, 1, 4]}
    in_half = {"a": [3, 2, 1, 1, 0], "b": [2, 0, 1, 0, 2]}
    cuts = {"a": np.array([-1., 0., 1.]), "b": np.array([-.5, .5])}
    rng = np.random.default_rng(12)
    code, stratum, lo, hi, part = [], [], [], [], []
    for s in ("a", "b"):
        for g, m in enumerate(sizes[s]):
            if m == 0:
                continue
            y = np.full(m, .3) if (s, g) == ("a", 1) else rng.normal(.4 * (g - 2), 1., m)
            l, h = _panel(y, cuts[s])
            for k in range(m):
                if k < in_half[s][g]:
                    part.append(len(code))
                code.append(g)
                stratum.append(s)
                lo.append(l[k])
                hi.append(h[k])
    code, lo, hi, part = np.array(code), np.array(lo), np.array(hi), np.array(part)
    stratum = np.array(stratum, dtype=object)
    for s in ("a", "b"):         # every reading of the panel is read in its stratum
        assert len({(a, b) for a, b, t in zip(lo, hi, stratum) if t == s}) == cuts[s].size + 1
    reads = Readings(lo, hi, code, stratum)
    n_cell = np.zeros((5, 2))
    np.add.at(n_cell, (code, reads.stratum), 1.)
    full = reads.cell_counts(np.arange(code.size))
    z = np.array([.4, -.3, .1, -.1, .25])
    got = _certificate(reads, z, reads.cell_counts(part), full, n_cell)
    want = _certificate_by_hand(lo, hi, code, stratum, part, z)
    assert got[0] == pytest.approx(want[0], abs=1e-14)
    assert got[1] == pytest.approx(want[1], rel=1e-12)
    assert got[2] == pytest.approx(want[2], rel=1e-12) and np.isfinite(got[2])
    # an evaluation half without an isolate: every pair at its worst case,
    # nothing random, so no variance and no finite degrees of freedom
    empty = _certificate(reads, z, reads.cell_counts(np.array([], dtype=int)), full, n_cell)
    assert empty[0] == pytest.approx(_certificate_by_hand(lo, hi, code, stratum, np.array([], dtype=int), z)[0],
                                     abs=1e-14)
    assert empty[1] == 0. and empty[2] == np.inf


def test_the_certificate_is_its_written_form_on_readings_that_overlap_in_a_stratum():
    """One stratum read on two panels whose wells interleave, so that some
    readings overlap and a pair of isolates is ordered only where one
    reading lies wholly below the other."""
    rng = np.random.default_rng(31)
    code = np.repeat(np.arange(5), [7, 5, 6, 4, 6])
    y = rng.normal(0, 1, code.size) + .6 * (code - 2)
    first = rng.random(code.size) < .5
    lo, hi = np.empty(code.size), np.empty(code.size)
    lo[first], hi[first] = _panel(y[first], np.array([-1., 0., 1.]))
    lo[~first], hi[~first] = _panel(y[~first], np.array([-.5, .5]))
    stratum = np.zeros(code.size, dtype=object)
    reads = Readings(lo, hi, code, stratum)
    assert not reads.sharp
    part = np.flatnonzero(rng.random(code.size) < .5)
    n_cell = np.zeros((5, 1))
    np.add.at(n_cell, (code, reads.stratum), 1.)
    z = np.array([.5, -.2, .1, -.4, .3])
    got = _certificate(reads, z, reads.cell_counts(part), reads.cell_counts(np.arange(code.size)), n_cell)
    want = _certificate_by_hand(lo, hi, code, stratum, part, z)
    assert got[0] == pytest.approx(want[0], abs=1e-14)
    assert got[1] == pytest.approx(want[1], rel=1e-12)
    assert got[2] == pytest.approx(want[2], rel=1e-12)


def test_the_limit_averages_the_scaled_estimates_of_every_split_both_ways(monkeypatch):
    """12 max(0, mean T - t sd)^2 over the two directions of SPLITS splits,
    t at the smallest of their degrees of freedom."""
    from scipy.stats import t as student
    rng = np.random.default_rng(14)
    code = np.repeat(np.arange(8), 6)
    y = 3. * rng.normal(0, 1, 8)[code] + rng.normal(0, 1, code.size)
    lo, hi = _panel(y, np.arange(-6., 7.))
    seen = []
    real = lo_mod._certified_terms

    def spy(*args):
        out = real(*args)
        assert len(out) == 2
        seen.extend(out)
        return out
    monkeypatch.setattr(lo_mod, "_certified_terms", spy)
    reads = Readings(lo, hi, code, np.zeros(code.size, dtype=object))
    limit = lo_mod._lower_limit(reads, .1, np.random.default_rng(3))
    assert len(seen) == 2 * lo_mod.SPLITS
    T, S, D = map(np.array, zip(*seen))
    by_hand = 12 * max(T.mean() - student.ppf(.9, D.min()) * S.mean(), 0.) ** 2
    assert limit == pytest.approx(by_hand, rel=1e-12) and limit > 0.


def test_one_direction_is_chosen_on_one_half_and_estimated_on_the_other():
    rng = np.random.default_rng(22)
    code = np.repeat(np.arange(6), 8)
    y = 2. * rng.normal(0, 1, 6)[code] + rng.normal(0, 1, code.size)
    lo, hi = _panel(y, np.arange(-4., 5.))
    reads = Readings(lo, hi, code, np.zeros(code.size, dtype=object))
    n_cell = np.bincount(code).astype(float)[:, None]
    full = reads.cell_counts(np.arange(code.size))
    terms = lo_mod._certified_terms(reads, np.random.default_rng(4), full, n_cell)
    first, second = _split(reads, np.random.default_rng(4))
    w = n_cell[:, 0] / n_cell.sum()
    for (T, sd, df), (choose, evaluate) in zip(terms, ((first, second), (second, first))):
        _, _, means = _minimum_norm(reads.polytope(choose, lineages=np.arange(6)))
        z = means - .5 - np.sum(w * (means - .5))
        estimate, variance, dof = _certificate(reads, z, reads.cell_counts(evaluate), full, n_cell)
        Q = np.sum(w * z * z)
        assert T == pytest.approx(estimate / np.sqrt(Q), rel=1e-12)
        assert sd == pytest.approx(np.sqrt(variance / Q), rel=1e-12) and df == dof


def test_a_cell_whose_isolates_share_one_reading_keeps_a_variance():
    """Two lineages, each read wholly in one reading: the pair estimate is
    the largest it can be, and its variance is not zero."""
    lo = np.r_[np.full(6, -np.inf), np.zeros(6)]
    hi = np.r_[np.zeros(6), np.full(6, np.inf)]
    code = np.repeat([0, 1], 6)
    reads = Readings(lo, hi, code, np.zeros(12, dtype=object))
    n_cell = np.bincount(code).astype(float)[:, None]
    full = reads.cell_counts(np.arange(12))
    estimate, variance, _ = _certificate(reads, np.array([-.5, .5]), full, full, n_cell)
    assert estimate == pytest.approx(.25 * .5) and variance > 0.


@pytest.mark.parametrize("sizes, bins", [([2, 2], 2), ([2] * 30, 9), ([2] * 30, 1000)])
def test_the_audit_nulls_give_no_positive_limit(sizes, bins):
    """The designs of the audit under which a lower limit of the old bootstrap
    was positive in most draws: lineages of one law, read on a fixed panel."""
    rng = np.random.default_rng(5)
    code = np.repeat(np.arange(len(sizes)), sizes)
    y = rng.random(code.size)
    lo, hi = _panel(y, np.arange(1, bins) / bins)
    r = order_bounds(lo, hi, code, seed=2)
    assert r["latent_order_lower_limit"] == 0.


def test_calls_and_readings_pass_every_setting_through():
    rng = np.random.default_rng(15)
    code = np.repeat(np.arange(10), 8)
    lab = np.tile(["x", "y"], 40)
    y = (2. * rng.normal(0, 1, 10)[code] + 1.5 * (lab == "y") + rng.normal(0, 1, 80) > .5).astype(float)
    settings = dict(strata=lab, alpha=.2, starts=3, seed=5)
    call = call_bounds(y, code, **settings)
    assert call["latent_order_lower_limit"] > 0.
    lo, hi = np.where(y == 1., 0., -np.inf), np.where(y == 1., np.inf, 0.)
    assert call == order_bounds(lo, hi, code, **settings)
    reads = Readings(lo, hi, code, np.asarray(lab, dtype=object))
    assert call == lo_mod.bounds_from_bins(reads, alpha=.2, starts=3, seed=5)
    assert call["latent_order_alpha"] == .2 and call["latent_order_upper_starts"] == 4
    assert call["latent_order_splits"] == lo_mod.SPLITS
    assert call != call_bounds(y, code, **{**settings, "strata": None})
    # the defaults, as the signatures state them
    defaults = dict(strata=None, alpha=.05, starts=10, seed=0)
    assert call_bounds(y, code)["latent_order_lower_limit"] > 0.
    assert call_bounds(y, code) == call_bounds(y, code, **defaults)
    assert order_bounds(lo, hi, code) == order_bounds(lo, hi, code, **defaults)
    pooled = Readings(lo, hi, code, np.zeros(80, dtype=object))
    assert lo_mod.bounds_from_bins(pooled) == lo_mod.bounds_from_bins(pooled, alpha=.05, starts=10, seed=0)


def test_two_lineages_are_bounded_and_fewer_leave_every_field_present_and_not_computed():
    lo, hi = np.repeat([-np.inf, 0.], 4), np.repeat([0., np.inf], 4)
    two = order_bounds(lo, hi, np.r_[0, 0, 0, 1, 0, 1, 1, 1])
    assert 0. < two["latent_order_lower"] < two["latent_order_upper"] and two["latent_order_readings"] == 2
    assert two["latent_order_upper_exact"] and two["latent_order_sharp"]
    one = order_bounds(lo, hi, np.zeros(8, dtype=int), alpha=.1, starts=4)
    assert one.keys() == two.keys()
    fixed = ("latent_order_readings", "latent_order_alpha", "latent_order_upper_starts", "latent_order_splits")
    flags = ("latent_order_upper_exact", "latent_order_sharp")
    assert all(np.isnan(one[k]) for k in one if k not in fixed + flags)
    assert tuple(one[k] for k in fixed) == (0, .1, 5, lo_mod.SPLITS)
    assert not any(one[k] for k in flags)


def test_strata_are_bounded_on_their_own_beside_the_pooled_share():
    """Opposite lineage effects in two strata cancel in the pooled share; the
    guaranteed share of every stratum on its own is reported beside it and is
    the lower end of that stratum's readings alone."""
    code = np.tile(np.repeat(np.arange(4), 5), 2)
    lab = np.repeat(["a", "b"], 20).astype(object)
    effect = np.where(lab == "a", 1., -1.) * np.array([-1.5, -.5, .5, 1.5])[code]
    y = effect + np.tile(np.linspace(-.3, .3, 5), 8)
    lo, hi = _panel(y, np.arange(-2., 3.))
    r = order_bounds(lo, hi, code, strata=lab)
    within = r["latent_order_lower_within_strata"]
    assert set(within) == {"a", "b"}
    for name in ("a", "b"):
        m = lab == name
        alone = order_bounds(lo[m], hi[m], code[m])["latent_order_lower"]
        assert within[name] == pytest.approx(alone, abs=2 * lo_mod.GAP_TOLERANCE)
        assert within[name] > .5
    assert r["latent_order_lower"] < .05
    assert "latent_order_lower_within_strata" not in order_bounds(lo, hi, code)


def test_overlapping_readings_of_one_stratum_are_marked_as_not_sharp():
    rng = np.random.default_rng(19)
    code = np.repeat(np.arange(6), 8)
    y = rng.normal(0, 1, code.size) + rng.normal(0, 1, 6)[code]
    one = order_bounds(*_panel(y, np.arange(-2., 3.)), code)
    assert one["latent_order_sharp"]
    lo, hi = _panel(y, np.arange(-2., 3.))
    lo2, hi2 = _panel(y, np.arange(-1.5, 3.))
    half = np.arange(code.size) % 2 == 0
    two = order_bounds(np.where(half, lo, lo2), np.where(half, hi, hi2), code)
    assert not two["latent_order_sharp"]
    assert two["latent_order_lower"] <= two["latent_order_upper"]


def test_a_reading_with_one_end_missing_is_left_out():
    rng = np.random.default_rng(16)
    code = np.repeat(np.arange(5), 4)
    lo, hi = _panel(rng.normal(0, 1, 20) + rng.normal(0, 1, 5)[code], np.arange(-2., 3.))
    base = order_bounds(lo[1:], hi[1:], code[1:])
    keys = ("latent_order_lower", "latent_order_upper", "latent_order_midpoint", "panel_resolution")
    for missing in (lo, hi):
        kept = missing[0]
        missing[0] = np.nan
        other = order_bounds(lo, hi, code)
        missing[0] = kept
        assert [other[k] for k in keys] == [base[k] for k in keys]


def test_a_resampled_polytope_counts_the_rows_drawn_and_reads_its_ranges_from_them():
    rng = np.random.default_rng(17)
    code = np.repeat(np.arange(4), 5)
    lo, hi = _panel(rng.normal(0, 1, 20), np.array([-1., 0., 1.]))
    reads = Readings(lo, hi, code, np.zeros(20, dtype=object))
    # a resample that holds no isolate of the lowest reading, some twice
    rows = np.flatnonzero(reads.bins > 0)
    rows = np.r_[rows, rows[:3]]
    poly = reads.polytope(rows)
    counts = np.zeros((4, reads.n_bins))
    np.add.at(counts, (code[rows], reads.bins[rows]), 1.)
    present = counts.sum(0) > 0
    assert not present.all()
    assert np.array_equal(poly.counts, counts[:, present])
    cum = np.r_[0., np.cumsum(counts.sum(0) / counts.sum())]
    assert np.allclose(poly.L, cum[:-1][present], atol=1e-15) and np.allclose(poly.H, cum[1:][present], atol=1e-15)


def test_readings_get_their_ranges_whatever_their_counts_and_whichever_strata_are_absent():
    """Two strata, each on one panel. In the first every reading holds one
    isolate; a resample may hold no isolate of the first stratum, and the
    second still gets its ranges."""
    lo = np.array([-np.inf, 0., 1., -np.inf, 0.])
    hi = np.array([0., 1., np.inf, 0., np.inf])
    reads = Readings(lo, hi, np.array([0, 0, 1, 1, 0]), np.array(["a", "a", "a", "b", "b"], dtype=object))
    low, high = reads.ranges(np.ones(5))
    assert np.allclose(low, [0., 1 / 3, 2 / 3, 0., .5]) and np.allclose(high, [1 / 3, 2 / 3, 1., .5, 1.])
    low, high = reads.ranges(np.array([0., 0., 0., 1., 3.]))
    assert np.allclose(low[3:], [0., .25]) and np.allclose(high[3:], [.25, 1.])


def test_a_run_of_one_cycle_returns_the_certificate_of_its_first_point_and_longer_runs_no_less(monkeypatch):
    """The first point is the greedy vertex of the order of the lineages'
    midpoint means; a run cut after one cycle returns rho there less its gap,
    with that gap. Every point visited certifies a lower end and the best is
    kept, so a longer run never returns less."""
    poly = _instance(np.random.default_rng(1), 20, 8)
    merged, _ = poly.merged()
    w = merged.w
    x = merged.vertex(np.argsort(merged.midpoint() / w, kind="stable"))
    q = merged.vertex(np.argsort(x / w, kind="stable"))
    gap = 24 * (np.sum(x * x / w) - np.sum(x * q / w))
    monkeypatch.setattr(lo_mod, "MAX_MAJOR", 1)
    value, returned = _lower(poly)
    assert returned == pytest.approx(gap, rel=1e-12)
    assert value == pytest.approx(max(0., merged.rho(x) - max(gap, lo_mod.GAP_TOLERANCE)), rel=1e-12)
    values = []
    for cap in range(1, 30):
        monkeypatch.setattr(lo_mod, "MAX_MAJOR", cap)
        values.append(_lower(poly)[0])
    assert all(later >= earlier for earlier, later in zip(values, values[1:]))


def test_a_run_that_needs_more_cycles_than_the_patience_still_reaches_the_tolerance(monkeypatch):
    """Twenty lineages on eight readings: every cycle lowers sum x^2 / w, and
    the gap reaches the tolerance only after more than PATIENCE cycles."""
    poly = _instance(np.random.default_rng(1), 20, 8)
    assert _lower(poly)[1] <= lo_mod.GAP_TOLERANCE
    monkeypatch.setattr(lo_mod, "MAX_MAJOR", lo_mod.PATIENCE)
    assert _lower(poly)[1] > lo_mod.GAP_TOLERANCE


def test_the_audit_collection_on_which_the_search_stopped_short_gets_its_exact_upper_end():
    """Eight lineages on six wells, 544 isolates: the collection of the audit
    of release 1.0.0 on which the insertion search alone returned 0.0256073,
    below the largest share an order attains, 0.0256500 (all 8! orders). The
    upper end is now exact and certified."""
    counts = [[5, 0, 11, 11, 0, 5], [2, 2, 11, 11, 2, 2], [48, 1, 16, 16, 1, 48], [6, 1, 0, 0, 1, 6],
              [28, 3, 25, 25, 3, 28], [6, 41, 12, 12, 41, 6], [4, 16, 0, 0, 16, 4], [5, 29, 0, 0, 29, 5]]
    ext = np.r_[-np.inf, np.arange(5.), np.inf]
    lo, hi, lineage = [], [], []
    for g, row in enumerate(counts):
        for k, c in enumerate(row):
            lo += [ext[k]] * c
            hi += [ext[k + 1]] * c
            lineage += [g] * c
    r = order_bounds(np.array(lo), np.array(hi), np.array(lineage), seed=0)
    assert r["latent_order_upper_exact"]
    assert r["latent_order_upper"] == pytest.approx(0.025650009991467826, abs=1e-12)
    assert r["latent_order_upper_bound"] >= r["latent_order_upper"] - 1e-15
