"""Clonal share: properties, planted truth, gates."""
import math

import numpy as np
import pytest

from amr_clonalshare.attribution import (MIN_REPEATED_LINEAGES, clonal_share,
                                            layer_clonal_share)


def _two_lineage_cohort(n=400, sep=1.0, seed=0):
    """Half the isolates in lineage A, half in B, prevalence separated by
    ``sep`` (0 = lineage carries nothing, 1 = lineage determines the trait)."""
    rng = np.random.default_rng(seed)
    lin = np.repeat(["A", "B"], n // 2)
    q = 0.5 + sep * np.where(lin == "A", -0.45, 0.45)
    y = (rng.random(n) < q).astype(float)
    return y, lin


def test_clonal_share_is_near_zero_when_lineage_carries_nothing():
    y, lin = _two_lineage_cohort(sep=0.0, seed=1)
    r = clonal_share(y, lin, n_boot=100, n_perm=100, seed=1)
    assert abs(r.kappa_adj) < 0.08
    assert r.observed_low <= 0.0 <= r.observed_high
    assert r.p_value > 0.05


def test_clonal_share_is_near_one_when_lineage_determines_the_trait():
    y, lin = _two_lineage_cohort(sep=1.0, seed=2)
    r = clonal_share(y, lin, n_boot=100, n_perm=100, seed=2)
    assert r.kappa_adj > 0.75
    assert r.p_value <= 0.05


def test_the_cost_of_estimating_lineage_means_is_removed_lineage_by_lineage():
    """Scored out of sample, a lineage mean estimated from a few training
    isolates pushes the plain score down; the lineage's own variance over its
    training members measures that cost, and removed it leaves the score of
    labels that carry nothing, and of their permutations, near zero. With many
    small lineages the gap must be visible."""
    from amr_clonalshare.attribution import _dealt_folds
    rng = np.random.default_rng(3)
    lin = rng.integers(0, 60, 300)
    y = (rng.random(300) < 0.4).astype(float)
    r = clonal_share(y, lin, n_boot=60, n_perm=120, seed=3)
    assert abs(r.null_mean) < 0.02 and abs(r.kappa_adj) < 0.12
    code = np.unique(lin, return_inverse=True)[1]
    fold = _dealt_folds(code, 5, np.random.default_rng(0))
    sse = sst = 0.
    for k in range(5):
        tr, te = fold != k, fold == k
        means = np.bincount(code[tr], weights=y[tr], minlength=60) / np.maximum(np.bincount(code[tr], minlength=60), 1)
        sse += np.sum((y[te] - means[code[te]]) ** 2)
        sst += np.sum((y[te] - y[tr].mean()) ** 2)
    assert 1 - sse / sst < -0.05


def test_singleton_lineages_are_set_aside_and_counted():
    """A singleton lineage cannot be predicted out of sample. The share is
    scored on the isolates of repeated lineages: the same number as the share
    of those isolates alone, with the singletons counted and the support
    reported, and no share at all when fewer than two lineages repeat."""
    rng = np.random.default_rng(4)
    y = (rng.random(300) < 0.4).astype(float)
    r = clonal_share(y, np.arange(300), n_boot=40, n_perm=40, seed=4)
    assert r.support == 0.0 and r.n_singletons_set_aside == 300
    assert r.estimable is False and math.isnan(r.kappa_adj)

    lin = np.r_[np.repeat(np.arange(20), 10), 100 + np.arange(100)]
    mixed = clonal_share(y, lin, n_boot=40, n_perm=40, seed=4)
    alone = clonal_share(y[:200], lin[:200], n_boot=40, n_perm=40, seed=4)
    assert mixed.n == 300 and mixed.n_scored == 200 and mixed.n_singletons_set_aside == 100
    assert mixed.support == pytest.approx(200 / 300)
    assert mixed.kappa_adj == alone.kappa_adj
    assert (mixed.observed_low, mixed.observed_high) == (alone.observed_low, alone.observed_high)
    assert mixed.estimable is True and MIN_REPEATED_LINEAGES == 2


def test_missing_lineage_labels_are_set_aside_and_counted():
    """An untyped isolate is not a lineage. The estimate describes the typed
    isolates, and the record says how many were set aside and what share of
    the cohort they were."""
    y, lin = _two_lineage_cohort(sep=0.6, seed=5)
    lin = np.array(lin, dtype=object)
    lin[:50] = None
    r = clonal_share(y, lin, n_boot=40, n_perm=40, seed=5)
    assert r.n == len(y) - 50
    assert r.n_dropped_untyped == 50
    assert r.missing_share == pytest.approx(50 / len(y))
    assert r.n_groups == 2                  # A and B; no __missing__ level
    typed = clonal_share(y[50:], lin[50:], n_boot=40, n_perm=40, seed=5)
    assert r.kappa == pytest.approx(typed.kappa)


def test_constant_trait_returns_nan_rather_than_a_number():
    lin = np.repeat(np.arange(10), 20)
    r = clonal_share(np.ones(200), lin, n_boot=20, n_perm=20, seed=12)
    assert np.isnan(r.kappa) or np.isnan(r.kappa_adj)


def test_a_non_finite_trait_value_is_dropped_rather_than_pinning_the_p_value_to_its_floor():
    """A missing well made every ``null >= kappa`` comparison false, so the
    permutation test counted no exceedances and reported 1 / (n_perm + 1) --
    the most significant value obtainable -- beside a NaN kappa_adj and
    estimable true, on a cohort where the lineage carries almost nothing."""
    rng = np.random.default_rng(21)
    lin = np.array([f"ST{i % 6}" for i in range(60)])
    y = (rng.random(60) < 0.4).astype(float)
    y[7] = np.nan
    r = clonal_share(y, lin, n_boot=50, n_perm=200, repeats=5, seed=1)
    assert r.n_dropped_non_finite == 1
    assert r.n == 59
    # Dropping the isolate is the whole of the treatment: the run must be the
    # run on the cohort with that isolate absent, number for number.
    keep = np.isfinite(y)
    clean = clonal_share(y[keep], lin[keep], n_boot=50, n_perm=200, repeats=5,
                         seed=1)
    assert clean.n_dropped_non_finite == 0
    for name in ("kappa", "kappa_adj", "observed_low", "observed_high", "null_mean",
                 "p_value"):
        assert getattr(r, name) == getattr(clean, name), name
    assert np.isfinite(r.kappa_adj)
    assert r.as_dict()["n_dropped_non_finite"] == 1


def test_a_length_mismatch_is_refused_rather_than_transposed_into_a_cohort_of_one():
    """Transposing on any mismatch turned 60 traits against 59 labels into
    n = 1 and reported a number for it. A transpose is the answer only when
    the transpose is what resolves the mismatch."""
    rng = np.random.default_rng(22)
    lin = np.array([f"ST{i % 6}" for i in range(60)])
    y = (rng.random(60) < 0.4).astype(float)
    with pytest.raises(ValueError, match="60 rows and lineage has 59"):
        clonal_share(y, lin[:59], n_boot=20, n_perm=20, repeats=2, seed=1)
    # a row vector is one trait, not a cohort of one
    assert layer_clonal_share(y.reshape(1, -1), lin, n_boot=20,
                              n_perm=20, repeats=2, seed=1).n == 60
    with pytest.raises(ValueError, match="one trait at a time"):
        layer_clonal_share(np.column_stack([y, y]), lin, n_boot=20, n_perm=20, repeats=2, seed=1)


def test_fewer_than_two_folds_is_refused_rather_than_read_as_no_lineage_signal():
    """One fold holds nothing out, so every score was taken on an empty
    held-out set and the run returned kappa 0, a zero-width interval, p = 1 and
    estimable true: a clean acquittal of the lineage that nothing was tested
    for. Zero folds reached a modulo by zero first."""
    rng = np.random.default_rng(23)
    lin = np.array([f"ST{i % 6}" for i in range(60)])
    y = (rng.random(60) < 0.4).astype(float)
    for folds in (1, 0):
        with pytest.raises(ValueError, match="at least two folds"):
            clonal_share(y, lin, folds=folds, n_boot=20, n_perm=20, repeats=2,
                         seed=1)


def test_a_cohort_with_one_lineage_is_refused_rather_than_scored_zero():
    """One lineage leaves nothing to contrast, so the share is undefined.

    Answering such a cohort with kappa 0, an interval of exactly zero width
    and p = 1, all flagged estimable, would be defensible arithmetic number by
    number and a false statement as a set:
    it reads as "lineage explains nothing" when nothing was compared. A public
    release that types a whole species as a single lineage produces exactly
    this shape, so the case is reachable on real data and not a curiosity.

    The design can fail rather than merely pass: the two-lineage arm on the
    same trait vector must still be estimable, so a guard that refused every
    cohort would be caught here.
    """
    rng = np.random.default_rng(0)
    y = (rng.random(120) < 0.25).astype(float)

    one = np.array(["PDS000000001"] * 120, dtype=object)
    refused = clonal_share(y, one, folds=5, repeats=4, n_boot=40, n_perm=40,
                           seed=1)
    assert refused.n_groups == 1
    assert not refused.estimable
    assert math.isnan(refused.kappa_adj)
    assert math.isnan(refused.observed_low) and math.isnan(refused.observed_high)

    two = np.array(["PDS000000001"] * 60 + ["PDS000000002"] * 60, dtype=object)
    scored = clonal_share(y, two, folds=5, repeats=4, n_boot=40, n_perm=40,
                          seed=1)
    assert scored.n_groups == 2
    assert math.isfinite(scored.kappa_adj)


def _many_lineage_cohort(n_groups, size, share, seed):
    rng = np.random.default_rng(seed)
    tau = math.sqrt(share / (1.0 - share))
    effects = rng.normal(0.0, tau, n_groups)
    lin = np.repeat(np.arange(n_groups), size)
    y = effects[lin] + rng.normal(0.0, 1.0, lin.size)
    return y, lin


def test_the_p_value_floor_is_carried_and_never_undercut():
    y, lin = _two_lineage_cohort(sep=0.9, seed=5)
    r = clonal_share(y, lin, n_boot=40, n_perm=30, seed=5)
    assert r.p_floor == pytest.approx(1.0 / 31.0)
    assert r.p_value >= r.p_floor
    assert math.isnan(clonal_share(y, lin, n_boot=40, n_perm=0, seed=5).p_floor)


def test_lineage_structure_diagnostics_are_reported():
    from amr_clonalshare.attribution import clonal_share
    y = np.r_[np.ones(10), np.zeros(10), [1, 0, 1]]
    lineage = np.array(['a']*10 + ['b']*10 + ['c', 'd', 'e'], dtype=object)
    r = clonal_share(y, lineage, folds=2, repeats=2, n_boot=10, n_perm=9, seed=1).as_dict()
    assert r['n_groups'] == 5 and r['n_groups_repeated'] == 2
    # the structure of the scored lineages; the three singletons are set aside
    assert r['n_singletons_set_aside'] == 3 and r['n_scored'] == 20
    assert r['effective_groups'] == pytest.approx(2.0)
    assert r['largest_group_share'] == pytest.approx(10/20)
    assert r['n_positive_repeated'] == 10


# --------------------------------------------------------------------------- #
# The interval for the lineages in hand
# --------------------------------------------------------------------------- #

def _majority_cohort():
    """Ten lineages of four, alternating rates of 0.05 and 0.95, every isolate
    in the majority state of its lineage: the sample of the audit."""
    lineage = np.repeat(np.arange(10), 4)
    call = np.repeat(np.arange(10) % 2, 4).astype(float)
    return call, lineage


def test_a_call_and_the_same_call_read_as_two_wells_give_one_interval():
    """The interval depends on the readings, not on how they are coded: a
    call, the call coded 3/7, and the call read as two wells of a panel give
    the same estimate and the same interval."""
    from amr_clonalshare.mic_order import mic_order_share
    call, lineage = _majority_cohort()
    a = clonal_share(call, lineage, n_perm=99, seed=1)
    b = clonal_share(3. + 4. * call, lineage, n_perm=99, seed=1)
    lo, hi = np.where(call == 0, -np.inf, 0.), np.where(call == 0, 0., np.inf)
    c = mic_order_share(lo, hi, lineage, n_perm=99, seed=1)
    for other in (b.as_dict(), c):
        for key in ("kappa_adj", "observed_low", "observed_high"):
            assert other[key] == pytest.approx(getattr(a, key), abs=1e-12)


def test_lineages_whose_calls_all_agree_still_leave_an_interval_of_width():
    """Every lineage homogeneous: resampled as they stand, the lineages would
    give one value and an interval of zero width at one. Drawn from Jeffreys'
    posterior mean, they give an interval that reaches the true share of the
    design, 0.81."""
    call, lineage = _majority_cohort()
    r = clonal_share(call, lineage, n_perm=99, seed=1)
    assert r.kappa_adj == pytest.approx(1.0)
    assert r.observed_low < 0.81 <= r.observed_high == 1.0
    assert r.observed_se > 0.05


def test_the_share_of_the_smoothed_world_is_the_share_by_its_definition():
    """eta* from the cell laws against the between-lineage share computed
    from the values and probabilities one isolate at a time."""
    from amr_clonalshare.attribution import _share_of
    rng = np.random.default_rng(8)
    values = np.array([.1, .4, .9])
    cell_lineage = np.array([0, 0, 1, 2])
    law = rng.dirichlet(np.ones(3), 4)
    cell_of = np.repeat(np.arange(4), [3, 2, 4, 5])
    sizes = np.bincount(cell_lineage[cell_of]).astype(float)
    mean_i = (law @ values)[cell_of]
    var_i = (law @ values ** 2)[cell_of] - mean_i ** 2
    code = cell_lineage[cell_of]
    mu = np.bincount(code, weights=mean_i) / sizes
    m2 = np.bincount(code, weights=var_i + mean_i ** 2) / sizes
    w = sizes / sizes.sum()
    between = w @ (mu - w @ mu) ** 2
    expected = between / (between + w @ (m2 - mu ** 2))
    assert _share_of(values, sizes, law, cell_of, cell_lineage) == pytest.approx(expected, rel=1e-12)


def test_the_average_of_the_skill_over_dealings_is_its_closed_form():
    """Averaged over the dealings of the folds, the corrected numerator of the
    skill is exactly sum_g m_g s_g^2; with the denominator ratio read from the
    dealings, the fold average is the closed form of _fold_average."""
    from amr_clonalshare.attribution import _dealt_folds, _fold_average, _skill_parts
    rng = np.random.default_rng(12)
    code = np.repeat(np.arange(6), [2, 3, 4, 5, 3, 7])
    y = rng.random(code.size)
    X = y.reshape(-1, 1)
    parts = np.array([_skill_parts(X, code, _dealt_folds(code, 5, rng)) for _ in range(4000)])
    sizes = np.bincount(code)
    s2 = np.array([np.var(y[code == g], ddof=1) for g in range(6)])
    assert parts[:, 0].mean() == pytest.approx(float(np.sum(sizes * s2)),
                                               abs=4 * parts[:, 0].std() / np.sqrt(4000))
    ratio = parts[:, 1].mean() / float(np.sum((y - y.mean()) ** 2))
    assert _fold_average(y, code, ratio) == pytest.approx(1 - parts[:, 0].mean() / parts[:, 1].mean(),
                                                          abs=0.01)


def test_the_estimate_and_the_test_do_not_depend_on_the_interval():
    y, lin = _two_lineage_cohort(sep=0.5, seed=21)
    with_interval = clonal_share(y, lin, n_boot=200, n_perm=99, seed=4)
    without = clonal_share(y, lin, n_boot=0, n_perm=99, seed=4)
    assert with_interval.kappa_adj == without.kappa_adj
    assert with_interval.p_value == without.p_value
    assert math.isnan(without.observed_low) and without.n_boot_used == 0
    assert with_interval.n_boot_used == 200


def test_the_interval_holds_the_estimate_in_a_planted_cohort():
    rng = np.random.default_rng(11)
    lin = np.repeat(np.arange(40), 8)
    y = (rng.random(320) < np.where(lin < 20, 0.25, 0.75)).astype(float)
    r = clonal_share(y, lin, n_perm=100, seed=11)
    assert r.observed_low <= r.kappa_adj <= r.observed_high
    assert 0.1 < r.observed_high - r.observed_low < 0.4


def test_sampling_units_confine_the_permutations_of_the_test():
    """With every lineage confined to its own farm there is nothing to compare
    within a farm: the permutations cannot move a label, so the p-value is
    one. The estimate and its interval describe the lineages in hand and do
    not change."""
    y, lin = _two_lineage_cohort(sep=0.8, seed=13)
    farms = np.where(lin == "A", "f1", "f2")
    r = clonal_share(y, lin, units=farms, n_boot=50, n_perm=49, seed=2)
    assert r.p_value == 1.0 and r.n_strata == 2
    free = clonal_share(y, lin, n_boot=50, n_perm=49, seed=2)
    assert free.p_value < 0.05
    for key in ("kappa_adj", "observed_low", "observed_high", "null_mean"):
        assert getattr(free, key) == getattr(r, key)


def test_the_derivatives_of_the_share_are_those_of_its_weighted_closed_form():
    """The infinitesimal jackknife of the fold-averaged share against central
    differences of the same closed form with one isolate's weight moved, the
    values held fixed."""
    from amr_clonalshare.attribution import _debias, _share_derivatives
    rng = np.random.default_rng(3)
    values = np.array([-1., .2, .5, 2.])
    code = np.repeat(np.arange(5), [2, 3, 4, 3, 6])
    category = rng.integers(0, 4, code.size)
    counts = np.zeros((5, 4))
    np.add.at(counts, (code, category), 1.)
    ratio, c_null = .93, -.04

    def share(weights):
        y = values[category]
        m = np.bincount(code, weights=weights)
        sums = np.bincount(code, weights=weights * y)
        squares = np.bincount(code, weights=weights * y * y)
        numerator = np.sum((squares - sums * sums / m) * m / (m - 1.))
        ybar = np.sum(weights * y) / weights.sum()
        return _debias(1. - numerator / (ratio * np.sum(weights * (y - ybar) ** 2)), c_null)

    D = _share_derivatives(values, counts, ratio, c_null)
    h = 1e-6
    for i in range(code.size):
        up, down = np.ones(code.size), np.ones(code.size)
        up[i] += h
        down[i] -= h
        assert D[code[i], category[i]] == pytest.approx((share(up) - share(down)) / (2 * h), rel=1e-6, abs=1e-9)


def test_the_standard_error_under_the_empirical_law_is_the_jackknife_and_never_zero_smoothed():
    from amr_clonalshare.attribution import _moment_se, _share_derivatives, _standardised
    values = np.array([0., 1.])
    counts = np.array([[4., 0.], [0., 4.], [2., 2.]])
    allowed = np.ones_like(counts)
    # the moments of every lineage's own isolates and of the categories its
    # smoothed distribution may take, as _Sampler.se forms them
    v = _standardised(values, counts.sum(0))
    powers = v[:, None] ** np.arange(5)[None, :]
    se = _moment_se(counts @ powers, (allowed @ powers) / allowed.sum(1, keepdims=True), 1., 0.)
    D = _share_derivatives(values, counts, 1., 0.)
    # the lineages whose isolates agree carry variance only through smoothing
    law = (counts + .5) / 5.
    by_hand = np.sqrt(np.sum(16. / 3. * np.sum(law * (D - np.sum(law * D, 1, keepdims=True)) ** 2, 1)))
    assert se == pytest.approx(by_hand, rel=1e-12) and se > 0.


def _share_with_one_lineage_replaced(values, code, category, g, law, ratio, c_null):
    """The share of _single_lineage_gain written out by hand: lineage g at
    ``law`` (its variance for its sum of squares, the law's moments in the
    total), every other lineage at its readings with its unbiased sum of
    squares."""
    from amr_clonalshare.attribution import _debias
    y = values[category]
    n = y.size
    num, s1, s2 = 0., 0., 0.
    for h in np.unique(code):
        yh = y[code == h]
        mh = yh.size
        if h == g:
            mu, m2 = law @ values, law @ values ** 2
            num += mh * (m2 - mu * mu)
            s1 += mh * mu
            s2 += mh * m2
        else:
            num += mh / (mh - 1.) * np.sum((yh - yh.mean()) ** 2)
            s1 += yh.sum()
            s2 += np.sum(yh ** 2)
    return _debias(1. - num / (ratio * (s2 - s1 * s1 / n)), c_null)


def _readings_share(values, code, category, ratio, c_null):
    from amr_clonalshare.attribution import _debias, _fold_average
    return _debias(_fold_average(values[category], code, ratio), c_null)


def test_the_single_lineage_end_of_a_call_is_the_best_end_of_a_binomial_region():
    """For a call the law of a lineage is its rate, and its region the rates
    p with 2 m KL(phat || p) at most the 95 % point of F(1, m - 1). The gain
    is checked against a fine grid of rates over every lineage's region."""
    from scipy import optimize, stats
    from amr_clonalshare.attribution import _single_lineage_gain
    rng = np.random.default_rng(5)
    k = np.r_[3, rng.binomial(10, .05, 11)]
    code = np.repeat(np.arange(k.size), 10)
    category = np.concatenate([np.r_[np.ones(j), np.zeros(10 - j)] for j in k]).astype(int)
    values = np.array([0., 1.])
    ratio, c_null = .96, .03
    base = _readings_share(values, code, category, ratio, c_null)
    best = -np.inf
    for g, j in enumerate(k):
        r = stats.f.ppf(.95, 1, 9) / 2.
        ph = j / 10.

        def kl(p, ph=ph):
            out = 0.
            if ph > 0:
                out += ph * np.log(ph / p)
            if ph < 1:
                out += (1 - ph) * np.log((1 - ph) / (1 - p))
            return 10. * out - r
        if 0 < ph < 1:
            lo = optimize.brentq(kl, 1e-12, ph)
            hi = optimize.brentq(kl, ph, 1 - 1e-12)
        else:
            lo = hi = ph
        for p in np.linspace(lo, hi, 2001):
            best = max(best, _share_with_one_lineage_replaced(values, code, category, g, np.array([1 - p, p]),
                                                              ratio, c_null))
    gain = _single_lineage_gain(values, code, category, ratio, c_null)
    assert gain == pytest.approx(best - base, abs=1e-7)
    assert gain > 0


def test_the_single_lineage_end_on_many_categories_is_the_maximum_over_the_region():
    """Against a direct maximisation over the laws of each lineage's region
    from many starts: the gain is never above it (every law used lies in the
    region) and reaches it."""
    from scipy import optimize, stats
    from amr_clonalshare.attribution import _single_lineage_gain
    rng = np.random.default_rng(11)
    code = np.repeat(np.arange(5), [6, 4, 8, 5, 7])
    category = rng.integers(0, 4, code.size)
    category[code == 2] = rng.integers(2, 4, 8)
    values = np.array([-1.3, -.2, .4, 1.9])
    ratio, c_null = .94, .01
    base = _readings_share(values, code, category, ratio, c_null)
    best = -np.inf
    for g in range(5):
        cats, counts = np.unique(category[code == g], return_counts=True)
        m = counts.sum()
        ph = counts / m
        r = stats.f.ppf(.95, 1, m - 1) / (2. * m)

        def law(x, cats=cats):
            p = np.zeros(4)
            p[cats] = np.exp(x) / np.exp(x).sum()
            return p

        def neg(x, g=g):
            return -_share_with_one_lineage_replaced(values, code, category, g, law(x), ratio, c_null)

        def inside(x, ph=ph, cats=cats, r=r):
            return r - np.sum(ph * np.log(ph / law(x)[cats]))
        for _ in range(25):
            x0 = np.log(ph) + rng.normal(0, .3, ph.size)
            if inside(x0) < 0:
                x0 = np.log(ph)
            fit = optimize.minimize(neg, x0, constraints=[{"type": "ineq", "fun": inside}], method="SLSQP",
                                    options={"ftol": 1e-12, "maxiter": 500})
            if inside(fit.x) >= -1e-9:
                best = max(best, -fit.fun)
    gain = _single_lineage_gain(values, code, category, ratio, c_null)
    assert gain <= best - base + 1e-7
    assert gain == pytest.approx(best - base, abs=1e-5)


def _two_category_oracle(values, code, category, ratio, c_null):
    """The largest share one lineage can give when every lineage was read in
    at most two categories: its law is (1 - p, p) on them, its region an
    interval of p found by root finding on the logit scale, and the share is
    maximised over a fine grid of that interval with its ends."""
    from scipy import optimize, stats
    best = -np.inf
    for g in np.unique(code):
        cats, counts = np.unique(category[code == g], return_counts=True)
        m = counts.sum()
        if cats.size == 1:
            law = np.zeros(values.size)
            law[cats[0]] = 1.
            best = max(best, _share_with_one_lineage_replaced(values, code, category, g, law, ratio, c_null))
            continue
        ph = counts[1] / m
        r = stats.f.ppf(.95, 1, max(m - 1, 1)) / 2.

        def kl(s, ph=ph, m=m, r=r):
            p = 1. / (1. + np.exp(-s))
            return m * (ph * np.log(ph / p) + (1 - ph) * np.log((1 - ph) / (1 - p))) - r
        s0 = np.log(ph / (1 - ph))
        lo = optimize.brentq(kl, -30., s0) if kl(-30.) > 0 else -30.
        hi = optimize.brentq(kl, s0, 30.) if kl(30.) > 0 else 30.
        for s in np.r_[np.linspace(lo, hi, 4001), lo, hi]:
            p = 1. / (1. + np.exp(-s))
            law = np.zeros(values.size)
            law[cats] = [1 - p, p]
            best = max(best, _share_with_one_lineage_replaced(values, code, category, g, law, ratio, c_null))
    return best - _readings_share(values, code, category, ratio, c_null)


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_the_single_lineage_end_of_small_lineages_is_the_best_end_of_their_regions(seed):
    """Lineages of two, three, four and ten isolates, each read in at most two
    of four categories, some of them below the rest so that the share rises
    as their mean falls: the gain against a one-dimensional oracle over each
    region, to a tolerance well inside what a coarse scan of directions or a
    changed calibration of a pair would miss."""
    from amr_clonalshare.attribution import _single_lineage_gain
    rng = np.random.default_rng(seed)
    sizes = np.r_[2, 2, 2, 3, 3, 4, 10, 10]
    code = np.repeat(np.arange(sizes.size), sizes)
    values = np.array([-1.1, .2, .7, 2.4])
    pairs = [rng.choice(4, 2, replace=False) for _ in sizes]
    category = np.concatenate([rng.choice(p, m) for p, m in zip(pairs, sizes)])
    category[code == 0] = pairs[0]
    ratio, c_null = .95, .02
    gain = _single_lineage_gain(values, code, category, ratio, c_null)
    oracle = _two_category_oracle(values, code, category, ratio, c_null)
    assert gain <= oracle + 1e-9
    assert gain == pytest.approx(oracle, abs=1e-7)


def test_a_pair_read_at_the_extremes_is_given_the_region_of_a_pair():
    """The largest rise comes from a pair that read the two outer categories
    while every other lineage read the two inner ones: its region, at the 95 %
    point of F(1, 1), reaches laws close to either outer category alone."""
    from amr_clonalshare.attribution import _single_lineage_gain
    rng = np.random.default_rng(8)
    code = np.repeat(np.arange(6), [2, 6, 6, 6, 6, 6])
    category = np.r_[0, 3, rng.integers(1, 3, 30)]
    values = np.array([-1.1, .2, .7, 2.4])
    gain = _single_lineage_gain(values, code, category, .95, .02)
    oracle = _two_category_oracle(values, code, category, .95, .02)
    assert gain == pytest.approx(oracle, abs=1e-7)


@pytest.mark.parametrize("seed", [100, 101])
def test_the_single_lineage_end_reaches_the_maximum_of_every_region_to_high_accuracy(seed):
    """Against a direct maximisation over the laws of each lineage's region
    from many starts, to 1e-8: the scan of directions and its refinement
    must reach the extreme point, not a direction near it."""
    from scipy import optimize, stats
    from amr_clonalshare.attribution import _single_lineage_gain
    rng = np.random.default_rng(seed)
    code = np.repeat(np.arange(6), rng.integers(3, 9, 6))
    category = rng.integers(0, 5, code.size)
    values = np.sort(rng.normal(size=5))
    ratio, c_null = .95, .02
    best = -np.inf
    for g in np.unique(code):
        cats, counts = np.unique(category[code == g], return_counts=True)
        m = counts.sum()
        ph = counts / m
        r = stats.f.ppf(.95, 1, m - 1) / (2. * m)
        if cats.size == 1:
            law = np.zeros(values.size)
            law[cats] = 1.
            best = max(best, _share_with_one_lineage_replaced(values, code, category, g, law, ratio, c_null))
            continue

        def law(x, cats=cats):
            p = np.zeros(values.size)
            e = np.exp(x - x.max())
            p[cats] = e / e.sum()
            return p

        def neg(x, g=g):
            return -_share_with_one_lineage_replaced(values, code, category, g, law(x), ratio, c_null)

        def inside(x, ph=ph, cats=cats, r=r):
            with np.errstate(divide="ignore"):
                return r - np.sum(ph * np.log(ph / law(x)[cats]))
        for _ in range(25):
            x0 = np.log(ph) + rng.normal(0, .5, ph.size)
            if inside(x0) < 0:
                x0 = np.log(ph)
            fit = optimize.minimize(neg, x0, constraints=[{"type": "ineq", "fun": inside}], method="SLSQP",
                                    options={"ftol": 1e-15, "maxiter": 1000})
            if inside(fit.x) >= -1e-12:
                best = max(best, -fit.fun)
    oracle = best - _readings_share(values, code, category, ratio, c_null)
    assert _single_lineage_gain(values, code, category, ratio, c_null) == pytest.approx(oracle, abs=1e-8)


def test_the_single_lineage_end_does_not_see_the_names_or_the_scale_of_the_values():
    from amr_clonalshare.attribution import _single_lineage_gain
    rng = np.random.default_rng(2)
    code = np.repeat(np.arange(9), rng.integers(2, 9, 9))
    category = rng.integers(0, 5, code.size)
    values = np.sort(rng.normal(size=5))
    gain = _single_lineage_gain(values, code, category, .97, .02)
    assert _single_lineage_gain(4. * values - 3., code, category, .97, .02) == pytest.approx(gain, rel=1e-9)
    renamed = np.array([3, 7, 0, 5, 1, 8, 2, 6, 4])[code]
    order = rng.permutation(code.size)
    assert _single_lineage_gain(values, renamed[order], category[order], .97, .02) == pytest.approx(gain, rel=1e-9)


def test_a_share_carried_by_one_lineage_read_close_to_the_rest_stays_within_the_interval():
    """Nineteen lineages of ten isolates rarely positive and one positive at
    a rate of one half, which reads three positives: the share of these
    laws is 0.143. The studentized bootstrap alone ends below it; the end
    that one lineage can carry keeps it."""
    from amr_clonalshare import attribution as A
    rng = np.random.default_rng(4)
    k = np.r_[3, rng.binomial(10, .05, 19)]
    lineage = np.repeat(np.arange(20), 10)
    y = np.concatenate([np.r_[np.ones(j), np.zeros(10 - j)] for j in k])
    p = np.r_[.5, np.full(19, .05)]
    truth = np.var(p) / (np.var(p) + np.mean(p * (1 - p)))
    assert truth == pytest.approx(.143, abs=1e-3)
    r = clonal_share(y, lineage, n_boot=199, n_perm=99, seed=3)
    assert r.observed_low <= truth <= r.observed_high
    kept = A._single_lineage_gain
    try:
        A._single_lineage_gain = lambda *a: float("nan")
        alone = clonal_share(y, lineage, n_boot=199, n_perm=99, seed=3)
    finally:
        A._single_lineage_gain = kept
    assert alone.observed_high < truth
    assert alone.observed_low == r.observed_low and alone.kappa_adj == r.kappa_adj


def test_lineages_that_read_alike_leave_an_interval_around_the_estimate():
    """Four lineages of five with two positives each: the share of the
    readings is at its minimum, where its first-order (jackknife) error is
    zero. The variance of the share without lineage differences keeps the
    interval from collapsing onto the estimate, on both sides."""
    lineage = np.repeat(np.arange(4), 5)
    y = np.tile([1., 1., 0., 0., 0.], 4)
    r = clonal_share(y, lineage, n_boot=199, n_perm=199, seed=4)
    assert r.observed_se > 0.
    assert r.observed_low < r.kappa_adj or r.observed_low == 0.
    assert r.observed_high > r.kappa_adj + .1


def test_the_interval_leaves_out_zero_only_where_the_studentized_permutation_test_rejects(monkeypatch):
    """The lower end is never above what the studentized shares of the
    permuted labellings allow: where it leaves out zero, the data's
    studentized share exceeds their 97.5 % point."""
    from amr_clonalshare import attribution as A
    seen = {}
    kept = A._smoothed_bootstrap

    def spy(*args, **kwargs):
        out = kept(*args, **kwargs)
        seen["t"], seen["se"], seen["t_perm"] = out
        return out
    monkeypatch.setattr(A, "_smoothed_bootstrap", spy)
    rng = np.random.default_rng(7)
    above = 0
    for trial in range(6):
        lineage = np.repeat(np.arange(8), 2)
        y = rng.normal(size=16) + (.8 * rng.normal(size=8))[lineage]
        r = clonal_share(y, lineage, n_boot=199, n_perm=199, seed=trial)
        assert len(seen["t_perm"]) == 199
        if r.observed_low > 0:
            above += 1
            assert r.kappa_adj / r.observed_se > np.quantile(seen["t_perm"], .975)
    assert above >= 1


@pytest.mark.slow
def test_pairs_without_a_lineage_effect_leave_zero_in_the_interval_at_the_level():
    """Six pairs of exchangeable continuous values: the lower end leaves out
    zero in at most about 2.5 % of datasets (bound at 5 % for 400 datasets)."""
    rng = np.random.default_rng(11)
    lineage = np.repeat(np.arange(6), 2)
    above = sum(clonal_share(rng.normal(size=12), lineage, folds=2, repeats=5, n_boot=199,
                             n_perm=199, seed=i).observed_low > 0 for i in range(400))
    assert above / 400 <= .05


@pytest.mark.parametrize("seed", range(6))
def test_the_ends_of_the_interval_never_cross(seed):
    """On collections without a lineage effect and small budgets the
    estimate can fall well below zero; both ends are kept within 0 and 1, so
    the interval is never reported with its ends crossed."""
    rng = np.random.default_rng(seed)
    lineage = np.repeat(np.arange(12), 4)
    y = rng.random(lineage.size)
    r = clonal_share(y, lineage, folds=2, repeats=2, n_boot=20, n_perm=19, seed=seed)
    assert 0. <= r.observed_low <= r.observed_high <= 1.


# --------------------------------------------------------------------------- #
# The fold dealing and the corrected numerator
# --------------------------------------------------------------------------- #

def _m_s2(X, code):
    """m s^2 summed over lineages and columns."""
    total = 0.
    for g in np.unique(code):
        Y = X[code == g]
        total += Y.shape[0] * Y.var(axis=0, ddof=1).sum()
    return total


@pytest.mark.parametrize("columns", [1, 3])
def test_the_corrected_numerator_of_a_lineage_is_m_times_its_variance(columns):
    """With no lineage larger than the number of folds, every held-out isolate
    is predicted by the other members of its lineage, and the corrected
    held-out error of the lineage is m s^2 on every fold draw."""
    from amr_clonalshare import attribution as A
    rng = np.random.default_rng(1)
    code = np.repeat(np.arange(12), rng.integers(2, 6, 12))
    X = rng.normal(0, 1, (12, columns))[code] + rng.normal(0, 1, (code.size, columns))
    for _ in range(5):
        numerator, _ = A._skill_parts(X, code, A._dealt_folds(code, 5, rng))
        assert numerator == pytest.approx(_m_s2(X, code), rel=1e-12)


def test_folds_are_dealt_so_that_no_lineage_sits_in_one_fold():
    from amr_clonalshare import attribution as A
    rng = np.random.default_rng(10)
    code = np.repeat(np.arange(40), rng.integers(2, 12, 40))
    for _ in range(20):
        fold = A._dealt_folds(code, 5, rng)
        sizes = np.bincount(fold, minlength=5)
        assert sizes.max() - sizes.min() <= 1
        for g in range(40):
            m = int((code == g).sum())
            assert np.unique(fold[code == g]).size == min(m, 5)


def test_dealt_folds_depend_on_which_isolates_belong_together_not_on_the_names():
    from amr_clonalshare import attribution as A
    code = np.repeat(np.arange(9), [2, 3, 4, 5, 6, 7, 8, 2, 3])
    a = A._dealt_folds(code, 5, np.random.default_rng(3))
    b = A._dealt_folds(8 - code, 5, np.random.default_rng(3))
    assert (a == b).all()


def test_one_analysis_selects_agents_by_the_permutation_p_values():
    from amr_clonalshare.core import _lineage_selection
    shares = {name: {"p_value": p, "p_floor": .001} for name, p in
              (("a", .001), ("b", .001), ("c", .001), ("d", .4), ("e", .9))}
    sel = _lineage_selection(shares, .05)
    assert sel["method"] == "benjamini_yekutieli"
    assert sel["rejected_features"] == ["a", "b", "c"]
    assert _lineage_selection({"x": {"p_value": float("nan")}}, .05)["n_tested"] == 0


def test_a_stratum_that_follows_lineage_exactly_removes_the_lineage_difference():
    """Two lineages, each in a stratum of its own, with rates 0.1 and 0.9.
    Pooled, the lineages carry most of the variance of the call; within
    strata they carry none, since each stratum holds one lineage. A call
    read within strata is centred on the prevalence of its stratum, so the
    share and its control are those of the same, within-stratum quantity and
    the estimate is near zero, as the manual says and as the score of a
    reading with one cut point makes it; pooled, the estimate is near the
    share of the pooled variance between the lineages, 0.64."""
    rng = np.random.default_rng(21)
    lin = np.repeat(["A", "B"], 200)
    y = (rng.random(400) < np.where(lin == "A", 0.1, 0.9)).astype(float)
    strata = np.where(lin == "A", "north", "south")
    within = clonal_share(y, lin, strata=strata, n_boot=100, n_perm=99, seed=3)
    pooled = clonal_share(y, lin, n_boot=100, n_perm=99, seed=3)
    assert abs(within.kappa_adj) < 0.05 and within.p_value == 1.0
    assert within.observed_low == 0.0
    assert abs(pooled.kappa_adj - 0.64) < 0.08 and pooled.p_value < 0.05
    assert within.prevalence == pooled.prevalence == y.mean()


def test_a_call_within_strata_is_the_reading_with_one_cut_point():
    """A call read within strata is scored as the MIC reading with one cut
    point of the same call is (Appendix A.2.4 of the article): the share, its control,
    its interval and its p-value of clonal_share with strata are those of
    mic_order_share on the readings (-inf, c] for 0 and (c, inf) for 1 in the
    same strata, whose mid-distribution score is an affine map of the centred
    call, stratum by stratum."""
    from amr_clonalshare.mic_order import mic_order_share
    rng = np.random.default_rng(5)
    lin = np.repeat(np.arange(12), 15)
    strata = rng.choice(["lab1", "lab2", "lab3"], size=lin.size)
    rate = np.where(lin < 6, 0.3, 0.7) + np.where(strata == "lab1", 0.15, 0.0)
    y = (rng.random(lin.size) < rate).astype(float)
    call = clonal_share(y, lin, strata=strata, n_boot=200, n_perm=199, seed=17)
    lo = np.where(y == 1.0, 0.0, -np.inf)
    hi = np.where(y == 1.0, np.inf, 0.0)
    order = mic_order_share(lo, hi, lin, strata=strata, n_boot=200, n_perm=199, seed=17)
    for key in ("kappa", "kappa_adj", "null_mean", "observed_low", "observed_high"):
        assert getattr(call, key) == pytest.approx(order[key], abs=1e-9), key
    assert call.p_value == order["p_value"]
