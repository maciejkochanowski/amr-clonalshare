"""The realised share: identities it must satisfy, and behaviour it must refuse.

The Monte-Carlo cells here are small and are a guard, not the calibration. The
operating characteristics quoted anywhere else come from
``benchmarks/realised_calibration.py``.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy import stats

from amr_clonalshare.realised import (KURTOSIS_LIMIT, MIN_GROUPS,
                                         realised_interval, realised_share,
                                         superpopulation_interval)


def balanced(n_groups, size, share, seed, noise=None):
    rng = np.random.default_rng(seed)
    tau = np.sqrt(share / (1.0 - share)) if share < 1.0 else 0.0
    effects = rng.normal(0.0, tau, n_groups)
    draw = noise(rng, (n_groups, size)) if noise else rng.normal(0.0, 1.0,
                                                                 (n_groups, size))
    y = (effects[:, None] + draw).ravel()
    lineage = np.repeat(np.arange(n_groups), size)
    realised = ((effects - effects.mean()) ** 2).sum() / (n_groups - 1)
    return y, lineage, realised / (realised + 1.0)


# ------------------------------------------------------------------ identities
def test_the_effective_size_is_the_common_size_on_a_balanced_design():
    y, lineage, _ = balanced(12, 7, 0.4, 1)
    assert realised_share(y, lineage).effective_group_size == pytest.approx(7.0)


def test_unequal_lineage_sizes_lower_the_effective_size_below_the_mean():
    rng = np.random.default_rng(2)
    sizes = np.array([2, 2, 3, 40, 40])
    lineage = np.repeat(np.arange(sizes.size), sizes)
    y = rng.normal(0.0, 1.0, sizes.sum())
    result = realised_share(y, lineage)
    assert result.effective_group_size < sizes.mean()


def test_the_point_estimate_is_the_classical_one_way_estimator():
    y, lineage, _ = balanced(15, 9, 0.35, 3)
    groups = y.reshape(15, 9)
    msb = 9 * ((groups.mean(1) - y.mean()) ** 2).sum() / 14
    msw = ((groups - groups.mean(1)[:, None]) ** 2).sum() / (135 - 15)
    expected = max((msb - msw) / 9, 0.0)
    expected = expected / (expected + msw)
    assert realised_share(y, lineage).kappa == pytest.approx(expected)


def test_the_share_is_invariant_to_shifting_and_rescaling_the_trait():
    y, lineage, _ = balanced(10, 8, 0.5, 4)
    plain = realised_share(y, lineage)
    moved = realised_share(3.5 * y - 17.0, lineage)
    assert moved.kappa == pytest.approx(plain.kappa)
    assert moved.ci_low == pytest.approx(plain.ci_low)
    assert moved.ci_high == pytest.approx(plain.ci_high)


def test_a_ratio_at_the_central_median_puts_no_lineage_effect_in_the_interval():
    df1, df2 = 20.0, 300.0
    median = stats.f.ppf(0.5, df1, df2)
    low, _ = realised_interval(median, df1, df2, 15.0, 21)
    assert low == 0.0


def test_the_endpoints_increase_with_the_observed_ratio():
    lows, highs = [], []
    for ratio in (1.5, 3.0, 6.0, 12.0):
        low, high = realised_interval(ratio, 29.0, 570.0, 20.0, 30)
        lows.append(low)
        highs.append(high)
    assert lows == sorted(lows)
    assert highs == sorted(highs)


def test_the_realised_interval_is_the_narrower_of_the_two_questions():
    narrower = 0
    for seed in range(40):
        y, lineage, _ = balanced(25, 12, 0.3, 500 + seed)
        r = realised_share(y, lineage)
        if (r.ci_high - r.ci_low) < (r.superpopulation_high
                                     - r.superpopulation_low):
            narrower += 1
    assert narrower == 40


def test_the_superpopulation_helper_matches_the_reported_field():
    y, lineage, _ = balanced(18, 6, 0.45, 9)
    r = realised_share(y, lineage)
    low, high = superpopulation_interval(r.f_ratio, r.n_groups - 1,
                                         r.n - r.n_groups,
                                         r.effective_group_size)
    assert (low, high) == pytest.approx((r.superpopulation_low,
                                         r.superpopulation_high))


def test_the_two_questions_leave_zero_together():
    """Detection is one decision, not two.

    Both intervals invert the same ratio of mean squares, so a lower endpoint
    stands above zero exactly when the same central F test rejects. The
    realised interval buys width, and validity under unequal lineage sizes,
    never power. A claim that it detects structure the other misses would be
    wrong, and this is the contract that keeps the claim out.
    """
    df_between, df_within, n0, g = 29.0, 570.0, 20.0, 30
    seen = set()
    for f_ratio in (0.4, 0.8, 1.0, 1.2, 1.3, 1.4, 1.5, 1.6, 2.0, 4.0, 10.0):
        low, _ = realised_interval(f_ratio, df_between, df_within, n0, g)
        slow, _ = superpopulation_interval(f_ratio, df_between, df_within, n0)
        assert (low > 0.0) == (slow > 0.0)
        seen.add(low > 0.0)
    # A test that only ever sees one side of the boundary is not a test.
    assert seen == {True, False}


# -------------------------------------------------------------------- refusals
def test_one_lineage_is_refused_because_it_carries_no_contrast():
    result = realised_share(np.arange(20.0), np.zeros(20))
    assert not result.estimable
    assert str(MIN_GROUPS) in result.reason


def test_a_constant_trait_is_refused_rather_than_scored():
    lineage = np.repeat(np.arange(5), 4)
    result = realised_share(np.ones(20), lineage)
    assert not result.estimable
    assert "within-lineage scale is zero" in result.reason


def test_identical_lineage_means_are_refused_rather_than_flagged_estimable():
    # Every lineage mean is 0.5: the between-lineage sum of squares is zero,
    # the inversion returns no interval, and a missing interval is not an
    # estimate however light the residual tails are.
    result = realised_share([0, 1, 0, 1, 0, 1], [0, 0, 1, 1, 2, 2])
    assert result.f_ratio == 0.0
    assert not np.isfinite(result.ci_low) and not np.isfinite(result.ci_high)
    assert not result.estimable
    assert "sum of squares is zero" in result.reason


def test_mismatched_lengths_raise_rather_than_align_silently():
    with pytest.raises(ValueError, match="must match"):
        realised_share(np.zeros(10), np.zeros(9))


def test_missing_observations_are_dropped_and_counted():
    y, lineage, _ = balanced(8, 5, 0.4, 11)
    y = y.copy()
    y[:3] = np.nan
    assert realised_share(y, lineage).n == y.size - 3


def test_heavy_tailed_residuals_are_refused_by_the_kurtosis_gate():
    result = realised_share(*balanced(30, 20, 0.3, 12,
                                      noise=lambda r, s: r.standard_t(3, s))[:2])
    assert result.residual_excess_kurtosis > KURTOSIS_LIMIT
    assert not result.estimable
    assert "excess kurtosis" in result.reason


# ------------------------------------------------------------------ behaviour
def test_destroying_the_lineage_structure_sends_the_share_to_zero():
    y, lineage, _ = balanced(30, 15, 0.6, 13)
    rng = np.random.default_rng(14)
    permuted = realised_share(y, rng.permutation(lineage))
    assert permuted.kappa < 0.05
    assert permuted.ci_low == 0.0


@pytest.mark.slow
def test_the_interval_covers_the_realised_share_it_targets():
    covered = 0
    trials = 300
    for seed in range(trials):
        y, lineage, truth = balanced(30, 20, 0.3, 900 + seed)
        r = realised_share(y, lineage)
        covered += r.ci_low <= truth <= r.ci_high
    assert 0.91 <= covered / trials <= 0.99


@pytest.mark.slow
def test_the_two_questions_agree_once_the_cohort_holds_many_lineages():
    y, lineage, _ = balanced(400, 20, 0.3, 21)
    r = realised_share(y, lineage)
    assert abs(r.ci_low - r.superpopulation_low) < 0.05
    assert abs(r.ci_high - r.superpopulation_high) < 0.05


def test_an_impossible_alpha_is_refused_and_mixed_type_labels_are_read_as_text():
    """Two ways the same function returned a confident answer to a question it
    could not answer. ``alpha=1`` gave a zero-width interval, ``alpha=2`` a
    reversed one and ``alpha=nan`` the pair (0, 0), all three flagged
    estimable. And a lineage vector mixing an integer 1 with a text "1" -- one
    cohort assembled from two pipelines -- raised a bare TypeError out of
    ``np.unique`` rather than being read as one lineage."""
    y, lineage, _ = balanced(6, 8, 0.4, 31)
    for bad in (1.0, 2.0, 0.0, float("nan")):
        with pytest.raises(ValueError, match="0 < alpha < 1"):
            realised_share(y, lineage, alpha=bad)
    assert realised_share(y, lineage, alpha=0.05).estimable

    mixed = np.array([1, "1", 2, "2", 1, 2], dtype=object)
    assert realised_share([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], mixed).n_groups == 2
    absent = np.array(["A", "B", np.nan, "A", "B", None], dtype=object)
    assert realised_share([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], absent).n_groups == 3
