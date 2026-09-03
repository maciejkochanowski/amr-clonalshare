"""Distance conventions, FDR, permutation p-values, effective dimension."""
import numpy as np
import pytest

from amr_clonalshare.stats import (benjamini_hochberg, binary_distance_matrix,
                                      effective_dimension, jaccard_distance_matrix,
                                      permutation_pvalue, uninformative_rows)


def test_jaccard_against_hand_computed_values():
    X = np.array([[1, 1, 0], [1, 0, 0], [0, 0, 1]])
    D = jaccard_distance_matrix(X)
    # rows 0,1: intersection 1, union 2 -> J = 0.5 -> D = 0.5
    assert D[0, 1] == pytest.approx(0.5)
    # rows 0,2: intersection 0, union 3 -> D = 1
    assert D[0, 2] == pytest.approx(1.0)
    assert np.allclose(np.diag(D), 0.0)
    assert np.allclose(D, D.T)


def test_empty_union_convention_is_explicit_and_changes_the_answer():
    X = np.array([[0, 0], [0, 0], [1, 1]])
    d_id = binary_distance_matrix(X, undefined_pair="identical")
    d_di = binary_distance_matrix(X, undefined_pair="distinct")
    d_nan = binary_distance_matrix(X, undefined_pair="nan")
    assert d_id[0, 1] == 0.0            # classical convention: identical
    assert d_di[0, 1] == 1.0            # refuse to merge trait-absent isolates
    assert np.isnan(d_nan[0, 1])
    with pytest.raises(ValueError):
        binary_distance_matrix(X, undefined_pair="whatever")


def test_uninformative_rows_flags_the_all_zero_stratum():
    X = np.array([[0, 0], [1, 0], [0, 0]])
    assert uninformative_rows(X).tolist() == [True, False, True]


def test_simple_matching_scores_shared_absence():
    X = np.array([[0, 0], [0, 0], [1, 1]])
    D = binary_distance_matrix(X, metric="simple_matching")
    assert D[0, 1] == pytest.approx(0.0)
    assert D[0, 2] == pytest.approx(1.0)


def test_effective_dimension_of_collinear_block_is_one():
    rng = np.random.default_rng(0)
    base = rng.integers(0, 2, size=(200, 1))
    collinear = np.repeat(base, 8, axis=1)
    assert effective_dimension(collinear) == pytest.approx(1.0, abs=1e-6)
    independent = rng.integers(0, 2, size=(2000, 8))
    assert effective_dimension(independent) > 7.0


def test_benjamini_hochberg_matches_the_step_up_definition():
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205])
    adj, rej = benjamini_hochberg(p, q=0.05)
    m = len(p)
    expected = np.minimum.accumulate((np.sort(p) * m / np.arange(1, m + 1))[::-1])[::-1]
    assert np.allclose(np.sort(adj), np.minimum(expected, 1.0))
    assert rej.sum() == int((adj <= 0.05).sum())


def test_permutation_pvalue_is_never_zero():
    null = np.zeros(5)
    p = permutation_pvalue(null, observed=10.0)
    assert p == pytest.approx(1 / 6)          # (0 + 1) / (5 + 1)
    assert p > 0


def test_benjamini_hochberg_refuses_a_failed_test_rather_than_zeroing_the_family():
    """One NaN used to make every adjusted value NaN, every rejection False, and
    a family with three clear discoveries read as "nothing significant" with no
    warning. A failed upstream test is not a null result, so it stops the run."""
    p = np.array([1e-9, 1e-8, 1e-7, 0.4, np.nan])
    with pytest.raises(ValueError, match="not finite"):
        benjamini_hochberg(p)
    adj, rej = benjamini_hochberg(p, nan_policy="omit")
    assert np.isnan(adj[4]) and not rej[4]
    assert rej[:3].all()
    assert adj[3] == pytest.approx(0.4)      # family is the four usable tests


def test_the_nan_guard_moves_no_clean_family():
    """Regression against the numbers this module already published: on families
    with no missing test the adjusted values must be bit-identical to before."""
    rng = np.random.default_rng(11)
    for _ in range(50):
        p = np.concatenate([rng.uniform(0, 1, 30), rng.uniform(0, 1e-4, 8),
                            np.full(5, 0.05)])
        adj, rej = benjamini_hochberg(p)
        m = len(p)
        order = np.argsort(p)
        want = np.minimum(np.minimum.accumulate(
            (p[order] * m / np.arange(1, m + 1))[::-1])[::-1], 1.0)
        ref = np.empty(m)
        ref[order] = want
        np.testing.assert_array_equal(adj, ref)
        np.testing.assert_array_equal(rej, ref <= 0.05)
