"""Anytime-valid evidence: the null property, growth, combination, e-BH.

An e-value earns its name by one property: its expectation under the null is
at most one, however long the data have been accumulating and however often
the running value has been looked at. Everything else the module offers rests
on that, so it is tested first and tested as a mean over many null cohorts
rather than as a single draw.
"""
import numpy as np
import pytest

from amr_clonalshare.evalues import (REJECT_AT, combine_independent, e_bh_log,
                                        combine_within_cohort, e_bh, e_process,
                                        sequential_e_process)


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


def _null_programme(rng, batches=8, n=30, G=10, prevalence=0.35):
    ys, ls = [], []
    for _ in range(batches):
        y, lab = _null_cohort(rng, n=n, G=G, prevalence=prevalence)
        ys.append(y)
        ls.append(lab)
    return ys, ls


def test_the_sequential_product_is_a_supermartingale_under_the_null():
    """Averaged over null programmes the running value after the last batch
    must not exceed one, and by Ville's inequality it may ever reach 1/alpha
    in at most a fraction alpha of programmes, however many looks are taken."""
    rng = np.random.default_rng(3)
    finals, ever = [], 0
    for _ in range(400):
        res = sequential_e_process(*_null_programme(rng))
        finals.append(res.e_value[-1])
        ever += res.reject_05
    assert float(np.mean(finals)) <= 1.0
    assert ever <= 20


def test_the_first_batch_earns_nothing_and_later_batches_earn_evidence():
    rng = np.random.default_rng(4)
    ys, ls = [], []
    G = 10
    q = np.clip(0.5 + 0.8 * (rng.random(G) - 0.5), 0.02, 0.98)
    for _ in range(8):
        lineage = rng.integers(0, G, 60)
        ys.append((rng.random(60) < q[lineage]).astype(float))
        ls.append(lineage.astype(object))
    res = sequential_e_process(ys, ls)
    assert res.log_e[0] == 0.0 and res.e_value[0] == 1.0
    assert res.log_e[-1] > np.log(REJECT_AT[0.05]) and res.reject_05
    assert res.n_batches == 8 and res.n_per_batch == [60] * 8


def test_sequential_batches_must_match_and_a_lone_lineage_earns_nothing():
    with pytest.raises(ValueError):
        sequential_e_process([np.zeros(4)], [np.zeros(4), np.zeros(4)])
    with pytest.raises(ValueError):
        sequential_e_process([np.zeros(4)], [np.zeros(3)])
    ys = [np.array([0, 1, 1, 0.0]), np.array([1, 0, 1, 1.0])]
    ls = [np.array(["A"] * 4, dtype=object), np.array(["A", "B", "A", "B"],
                                                       dtype=object)]
    res = sequential_e_process(ys, ls)
    assert res.log_e == [0.0, 0.0]


def test_infinite_evidence_is_rejected_not_discarded():
    # A running product over many intakes can leave the floating-point range;
    # the procedure then works on the log scale and infinity ranks first.
    assert e_bh([float("inf")])["rejected"] == [0]
    assert e_bh([1000.0])["rejected"] == [0]
    dec = e_bh([float("inf"), 1.0, float("nan"), 50.0])
    assert dec["rejected"] == [0, 3] and dec["m"] == 4
    assert e_bh_log([1386.29, 0.0])["rejected"] == [0]
    with pytest.raises(ValueError):
        e_bh([-1.0])


def test_the_sequential_product_survives_overflow_into_the_panel_decision():
    n = 2000
    lab = np.repeat([0, 1], n // 2)
    y = lab.astype(float)
    r = sequential_e_process([y, y], [lab, lab])
    assert r.log_e[-1] > 700 and r.e_value[-1] == float("inf")
    assert e_bh([r.e_value[-1]])["n_rejected"] == 1
    assert e_bh_log([r.log_e[-1]])["n_rejected"] == 1


def test_an_untyped_isolate_is_not_a_lineage_for_the_evidence_either():
    """A trait that depends only on whether an isolate was typed must earn no
    evidence: the untyped isolates are set aside, as every other estimator
    sets them aside, and the record says how many."""
    rng = np.random.default_rng(11)
    y, lineage = _null_cohort(rng, n=400, G=10)
    lineage = lineage.astype(object)
    lineage[:170] = None
    y[:170] = (rng.random(170) < 0.9).astype(float)      # the untyped are resistant
    r = e_process(y, lineage, folds=5, repeats=4, seed=1)
    assert r.n_dropped_untyped == 170 and r.n == 230 and r.n_groups == 10
    assert r.e_value < REJECT_AT[0.05]
    same = e_process(y[170:], lineage[170:], folds=5, repeats=4, seed=1)
    assert r.e_value == pytest.approx(same.e_value)
    seq = sequential_e_process([y[:200], y[200:]], [lineage[:200], lineage[200:]])
    assert seq.n_dropped_untyped == 170 and seq.n_per_batch == [30, 200]
