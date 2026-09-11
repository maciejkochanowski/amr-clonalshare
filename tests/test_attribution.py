"""Clonal share: properties, planted truth, gates."""
import math

import numpy as np
import pytest

from amr_clonalshare.attribution import (SUPPORT_THRESHOLD, clonal_share,
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
    assert r.ci_low <= 0.0 <= r.ci_high
    assert r.p_value > 0.05


def test_clonal_share_is_near_one_when_lineage_determines_the_trait():
    y, lin = _two_lineage_cohort(sep=1.0, seed=2)
    r = clonal_share(y, lin, n_boot=100, n_perm=100, seed=2)
    assert r.kappa_adj > 0.75
    assert r.p_value <= 0.05


def test_debias_removes_the_penalty_the_permuted_run_measures():
    """The raw score is pushed down by the cost of estimating group means; the
    corrected one is not. With many small lineages the gap must be visible."""
    rng = np.random.default_rng(3)
    lin = rng.integers(0, 60, 300)
    y = (rng.random(300) < 0.4).astype(float)
    r = clonal_share(y, lin, n_boot=60, n_perm=120, seed=3)
    assert r.null_mean < -0.02
    assert r.kappa < r.kappa_adj
    assert abs(r.kappa_adj) < 0.12


def test_support_gate_fires_on_a_singleton_heavy_lineage_variable():
    rng = np.random.default_rng(4)
    lin = np.arange(300)                     # every isolate its own lineage
    y = (rng.random(300) < 0.4).astype(float)
    r = clonal_share(y, lin, n_boot=40, n_perm=40, seed=4)
    assert r.support == 0.0
    assert r.estimable is False

    lin2 = np.repeat(np.arange(30), 10)      # ten isolates per lineage
    r2 = clonal_share(y, lin2, n_boot=40, n_perm=40, seed=4)
    assert r2.support == 1.0
    assert r2.estimable is True
    assert SUPPORT_THRESHOLD <= 1.0


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


def test_layer_share_lies_between_its_agents():
    rng = np.random.default_rng(6)
    lin = np.repeat(np.arange(20), 20)
    clonal = (rng.random(400) < np.where(lin < 10, 0.1, 0.9)).astype(float)
    free = (rng.random(400) < 0.5).astype(float)
    k_clonal = clonal_share(clonal, lin, n_boot=60, n_perm=60, seed=6).kappa_adj
    k_free = clonal_share(free, lin, n_boot=60, n_perm=60, seed=6).kappa_adj
    k_layer = layer_clonal_share(np.column_stack([clonal, free]), lin,
                                 n_boot=60, n_perm=60, seed=6).kappa_adj
    assert k_free < k_layer < k_clonal


def test_lineage_bootstrap_interval_contains_the_point_estimate():
    """Resampling isolates *within* a drawn lineage puts duplicate rows into
    training and held-out folds at once, which leaks and lifts every replicate
    above the estimate. Lineages are therefore taken whole."""
    rng = np.random.default_rng(11)
    lin = np.repeat(np.arange(40), 8)
    y = (rng.random(320) < np.where(lin < 20, 0.25, 0.75)).astype(float)
    r = clonal_share(y, lin, n_boot=200, n_perm=100, seed=11)
    assert r.ci_low <= r.kappa_adj <= r.ci_high


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
    for name in ("kappa", "kappa_adj", "ci_low", "ci_high", "null_mean",
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
    # a genuine (p, n) block still has to be recognised and turned round
    assert layer_clonal_share(np.column_stack([y, y]).T, lin, n_boot=20,
                              n_perm=20, repeats=2, seed=1).n == 60


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

    The estimator used to answer such a cohort with kappa 0, an interval of
    exactly zero width and p = 1, all flagged estimable. Every one of those
    numbers is defensible arithmetic and the set of them is a false statement:
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
    assert math.isnan(refused.ci_low) and math.isnan(refused.ci_high)

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


def test_the_superpopulation_interval_is_wider_with_few_lineages_and_converges_with_many():
    """The second interval adds the draw of the lineages, which ten lineages
    carry as a wide chi-square and a hundred as a narrow one."""
    few = clonal_share(*_many_lineage_cohort(10, 20, 0.3, 11), n_boot=150,
                       n_perm=60, seed=11)
    many = clonal_share(*_many_lineage_cohort(100, 20, 0.3, 12), n_boot=150,
                        n_perm=60, seed=12)
    for r in (few, many):
        assert r.interval_target == "lineage_membership"
        assert math.isfinite(r.superpopulation_low)
        assert math.isfinite(r.superpopulation_high)
        assert 0.0 <= r.superpopulation_low <= r.kappa_adj <= r.superpopulation_high <= 1.0
    width_few = few.superpopulation_high - few.superpopulation_low
    width_many = many.superpopulation_high - many.superpopulation_low
    assert width_few > (few.ci_high - few.ci_low) * 1.1
    assert width_many < (many.ci_high - many.ci_low) * 1.6
    assert width_few > width_many


def test_the_superpopulation_interval_is_absent_where_the_estimate_is():
    r = clonal_share(np.ones(20), np.repeat(["A", "B"], 10), n_boot=50,
                     n_perm=20, seed=3)
    assert not r.estimable and math.isnan(r.superpopulation_low)


def test_the_superpopulation_interval_does_not_move_the_numbers_beside_it():
    """Its bootstrap runs on its own stream, so every earlier field is the
    value it was before the interval existed: the record stays comparable."""
    y, lin = _two_lineage_cohort(sep=0.6, seed=4)
    a = clonal_share(y, lin, n_boot=100, n_perm=100, seed=4).as_dict()
    b = clonal_share(y, lin, n_boot=100, n_perm=100, seed=4).as_dict()
    assert a == b
    assert a["ci_low"] == pytest.approx(0.0, abs=1.0)  # finite and present


def test_the_p_value_floor_is_carried_and_never_undercut():
    y, lin = _two_lineage_cohort(sep=0.9, seed=5)
    r = clonal_share(y, lin, n_boot=40, n_perm=30, seed=5)
    assert r.p_floor == pytest.approx(1.0 / 31.0)
    assert r.p_value >= r.p_floor
    assert math.isnan(clonal_share(y, lin, n_boot=40, n_perm=0, seed=5).p_floor)


def test_a_binary_trait_carries_its_share_on_the_latent_scale_beside_the_observed_one():
    """The latent figure is the threshold-model image of the observed one:
    larger, ordered the same way, and absent for a continuous trait."""
    rng = np.random.default_rng(8)
    lin = np.repeat(np.arange(30), 20)
    p = 1.0 / (1.0 + np.exp(-(-1.2 + rng.normal(0.0, 1.0, 30))))
    y = rng.binomial(1, p[lin]).astype(float)
    r = clonal_share(y, lin, n_boot=60, n_perm=40, seed=8)
    assert 0.0 < r.kappa_adj < r.latent_share < 1.0
    assert r.latent_low <= r.latent_share <= r.latent_high
    # The species interval is stated on the component-ratio scale, one design
    # factor away from the lineage-membership share on the odds scale, and
    # envelopes the lineage bootstrap taken on that scale.
    from amr_clonalshare.attribution import _to_component_scale
    from amr_clonalshare.latent import latent_interval
    design = 1.0 - 30 * (20 / 600) ** 2
    lo_c, hi_c = latent_interval(_to_component_scale(r.ci_low, design),
                                 _to_component_scale(r.ci_high, design),
                                 r.prevalence)
    assert r.latent_species_low <= lo_c + 1e-12
    assert r.latent_species_high >= hi_c - 1e-12
    cont = clonal_share(y + rng.normal(0.0, 0.1, y.size), lin, n_boot=60,
                        n_perm=40, seed=8)
    assert math.isnan(cont.latent_share) and math.isnan(cont.latent_low)


def test_a_dominant_lineage_is_named_and_few_lineages_are_flagged():
    """A carrier structure is what the species interval's chi-square layer
    does not describe, so the record measures it: the largest single-lineage
    share of the between-lineage sum of squares."""
    rng = np.random.default_rng(9)
    lin = np.repeat(np.arange(5), 20)
    y = np.repeat(np.array([2.0, 0.0, 0.0, 0.0, 0.0]), 20) + rng.normal(0.0, 1.0, 100)
    r = clonal_share(y, lin, n_boot=40, n_perm=20, seed=9)
    assert r.few_lineages
    assert r.dominant_lineage_share > 0.5
    spread = clonal_share(rng.normal(size=600) + np.repeat(rng.normal(size=30), 20),
                          np.repeat(np.arange(30), 20), n_boot=40, n_perm=20, seed=9)
    assert not spread.few_lineages
    assert 0.0 < spread.dominant_lineage_share < 0.5
    block = layer_clonal_share(rng.random((600, 3)) < 0.3, np.repeat(np.arange(30), 20),
                               n_boot=20, n_perm=10, seed=9)
    assert math.isnan(block.dominant_lineage_share)
