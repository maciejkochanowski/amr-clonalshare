"""Anytime-valid evidence: the null property, growth, combination, e-BH.

An e-value earns its name by one property: its expectation under the null is
at most one, however long the data have been accumulating and however often
the running value has been looked at. Everything else the module offers rests
on that, so it is tested first and tested as a mean over many null cohorts
rather than as a single draw.
"""
import numpy as np
import pytest

from amr_clonalshare.evalues import (REJECT_AT, combine_independent,
                                        combine_within_cohort, e_bh, e_process)


def _null_cohort(rng, n=300, G=20, prevalence=0.35):
    """Lineage labels that carry nothing about the trait."""
    lineage = rng.integers(0, G, n).astype(object)
    y = (rng.random(n) < prevalence).astype(float)
    return y, lineage


def _signal_cohort(rng, n=300, G=20, sep=0.8):
    lineage = rng.integers(0, G, n)
    base = rng.random(G)
    q = np.clip(0.5 + sep * (base[lineage] - 0.5), 0.02, 0.98)
    y = (rng.random(n) < q).astype(float)
    return y, lineage.astype(object)


def test_expectation_under_the_null_is_at_most_one():
    """The defining property. Averaged over independent null cohorts the
    e-value must not exceed one; a procedure whose null mean drifts above one
    would spend an error budget it has not been granted."""
    rng = np.random.default_rng(0)
    values = [e_process(*_null_cohort(rng), folds=5, repeats=4,
                        seed=int(rng.integers(1 << 30))).e_value
              for _ in range(120)]
    assert float(np.mean(values)) <= 1.0


def test_rejection_is_rare_under_the_null():
    """Ville's inequality caps the probability that the running value ever
    reaches 1/alpha. A single look must therefore reject far less often than
    alpha, and over this many null cohorts must reject at most a handful."""
    rng = np.random.default_rng(1)
    rejects = sum(e_process(*_null_cohort(rng), folds=5, repeats=4,
                            seed=int(rng.integers(1 << 30))).reject_05
                  for _ in range(120))
    assert rejects <= 6


def test_evidence_grows_when_the_lineage_carries_the_trait():
    rng = np.random.default_rng(2)
    res = e_process(*_signal_cohort(rng), folds=5, repeats=10, seed=2)
    assert res.e_value > REJECT_AT[0.05]
    assert res.reject_05 and res.log_e > 0


def test_more_data_gives_more_evidence():
    """Evidence accumulates. The e-value on a larger cohort drawn from the
    same generating process must be the larger of the two, which is what makes
    the running value worth re-inspecting as a surveillance year closes."""
    small = e_process(*_signal_cohort(np.random.default_rng(3), n=200),
                      folds=5, repeats=10, seed=3)
    large = e_process(*_signal_cohort(np.random.default_rng(3), n=1200),
                      folds=5, repeats=10, seed=3)
    assert large.log_e > small.log_e


def test_thresholds_are_the_reciprocal_of_alpha():
    assert REJECT_AT[0.05] == pytest.approx(20.0)
    assert REJECT_AT[0.01] == pytest.approx(100.0)


def test_independent_batches_multiply_and_dependent_ones_average():
    assert combine_independent([2.0, 3.0, 4.0]) == pytest.approx(24.0)
    assert combine_within_cohort([2.0, 3.0, 4.0]) == pytest.approx(3.0)
    with pytest.raises(ValueError):
        combine_independent([1.0, -1.0])
    with pytest.raises(ValueError):
        combine_within_cohort([])


def test_e_bh_rejects_the_largest_k_that_clears_its_threshold():
    out = e_bh([100.0, 80.0, 0.5, 0.4, 0.3], alpha=0.05)
    assert out["n_rejected"] == 2
    assert out["rejected"] == [0, 1]
    assert out["threshold"] == pytest.approx(5 / (0.05 * 2))


def test_e_bh_rejects_nothing_when_no_value_clears_the_smallest_threshold():
    out = e_bh([1.0, 1.0, 1.0], alpha=0.05)
    assert out["n_rejected"] == 0
    assert out["rejected"] == []
    assert not np.isfinite(out["threshold"])


def test_e_bh_handles_an_empty_panel():
    out = e_bh([], alpha=0.05)
    assert out["n_rejected"] == 0 and out["m"] == 0


def test_e_bh_controls_the_discovery_rate_under_dependence():
    """The panel of an antimicrobial susceptibility test is not independent:
    cross-resistance makes a macrolide block behave as one trait. e-BH holds
    under arbitrary dependence, so a panel of perfectly correlated nulls must
    still produce few false rejections."""
    rng = np.random.default_rng(5)
    false_rejections = 0
    for _ in range(60):
        y, lineage = _null_cohort(rng, n=300)
        # thirteen copies of one null trait: maximal positive dependence
        values = [e_process(y, lineage, folds=5, repeats=4,
                            seed=int(rng.integers(1 << 30))).e_value
                  for _ in range(13)]
        false_rejections += e_bh(values, alpha=0.05)["n_rejected"] > 0
    assert false_rejections <= 6


def test_lineage_length_mismatch_is_refused():
    with pytest.raises(ValueError):
        e_process(np.zeros(10), ["A"] * 9)


def test_result_reports_the_design_it_came_from():
    rng = np.random.default_rng(7)
    y, lineage = _signal_cohort(rng, n=240, G=12)
    res = e_process(y, lineage, folds=4, repeats=3, seed=7)
    assert res.n == 240
    assert res.n_groups == 12
    assert res.n_splits == 12
    assert res.prevalence == pytest.approx(float(y.mean()))
