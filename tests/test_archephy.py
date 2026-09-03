"""Exact null of the k-way joint-homoplasy convergence test."""
import numpy as np
import pytest
from scipy import stats

dendropy = pytest.importorskip("dendropy")

from amr_clonalshare.archephy import (_exact_upper_tail, _poisson_upper_tail,
                                         archephy_cs_test, fitch_change_edges,
                                         joint_homoplasy_pvalue, load_tree)


@pytest.mark.parametrize("E,c1,c2,m", [(30, 8, 7, 2), (40, 12, 10, 3),
                                       (25, 6, 6, 1), (50, 20, 18, 7)])
def test_exact_null_equals_the_hypergeometric_at_k_equals_two(E, c1, c2, m):
    assert _exact_upper_tail(E, [c1, c2], m) == pytest.approx(
        float(stats.hypergeom.sf(m - 1, E, c1, c2)), rel=1e-9)


@pytest.mark.parametrize("E,c,m", [(30, [8, 7, 6], 2), (40, [12, 10, 9], 3),
                                   (25, [6, 6, 5], 1)])
def test_exact_null_matches_brute_force_monte_carlo_at_k_equals_three(E, c, m):
    rng = np.random.default_rng(0)
    hits = 0
    N = 40000
    for _ in range(N):
        sets = [set(rng.choice(E, size=ci, replace=False)) for ci in c]
        hits += len(set.intersection(*sets)) >= m
    assert _exact_upper_tail(E, c, m) == pytest.approx(hits / N, abs=0.01)


def test_exact_null_is_calibrated_and_the_poisson_approximation_is_conservative():
    """Where the intersection is not rare the Poisson approximation loses power.

    At E = 20 with three traits changing on 10 edges each, lambda = 2.5 and the
    Poisson tail is far too large: it rejects at 0.5% where the exact test
    rejects at 4%, i.e. it throws away most of the available power in exactly
    the dense-homoplasy regime a small tree produces.
    """
    rng = np.random.default_rng(1)
    E, c = 20, [10, 10, 10]
    exact, poisson = [], []
    for _ in range(4000):
        sets = [set(rng.choice(E, size=ci, replace=False)) for ci in c]
        m = len(set.intersection(*sets))
        exact.append(_exact_upper_tail(E, c, m) if m > 0 else 1.0)
        poisson.append(_poisson_upper_tail(E, c, m) if m > 0 else 1.0)
    exact, poisson = np.asarray(exact), np.asarray(poisson)
    assert np.mean(exact <= 0.05) <= 0.05          # valid
    assert np.mean(poisson <= 0.05) < np.mean(exact <= 0.05)   # conservative


def test_joint_homoplasy_pvalue_dispatch():
    assert joint_homoplasy_pvalue(30, [8, 7], 0) == 1.0
    p_exact = joint_homoplasy_pvalue(30, [8, 7, 6], 3, method="exact")
    p_pois = joint_homoplasy_pvalue(30, [8, 7, 6], 3, method="poisson")
    assert 0 < p_exact <= 1 and 0 < p_pois <= 1
    with pytest.raises(ValueError):
        joint_homoplasy_pvalue(30, [8, 7], 2, method="nonsense")


def _star_tree(n):
    return "(" + ",".join(f"t{i}:1.0" for i in range(n)) + ");"


def test_fitch_change_edges_on_a_star_tree():
    tree, leaves = load_tree(_star_tree(6))
    trait = np.array([1, 1, 1, 0, 0, 0])
    changed = fitch_change_edges(tree, leaves, trait)
    assert 0 < len(changed) <= 6


def test_convergence_test_requires_two_joint_change_edges():
    """One joint change edge is a single clonal co-origin, never convergence."""
    tree, leaves = load_tree(_star_tree(8))
    X = np.zeros((8, 2), dtype=int)
    X[0] = 1                                   # both traits, one leaf only
    out = archephy_cs_test(_star_tree(8), X, leaves, [0, 1])
    assert out["m_obs"] <= 1
    assert out["call"] is False


def test_convergence_test_fires_on_repeated_joint_gains():
    tree, leaves = load_tree(_star_tree(12))
    X = np.zeros((12, 2), dtype=int)
    # Two traits, not one trait written twice: identical columns are collapsed
    # before the test, so the pair has to differ somewhere off the joint edges.
    X[[0, 1, 2], 0] = 1                         # three independent co-origins
    X[[0, 1, 2, 7], 1] = 1                      # the partner, gained once more
    out = archephy_cs_test(_star_tree(12), X, leaves, [0, 1])
    assert out["m_obs"] >= 2
    assert out["p_value"] < 0.05
    assert out["call"] is True


def test_convergence_test_validates_its_inputs():
    tree_str = _star_tree(5)
    _, leaves = load_tree(tree_str)
    X = np.zeros((5, 3), dtype=int)
    with pytest.raises(ValueError):
        archephy_cs_test(tree_str, X, leaves, [0])          # k < 2
    with pytest.raises(ValueError):
        archephy_cs_test(tree_str, X, leaves, [0, 99])      # index out of range
    with pytest.raises(ValueError):
        archephy_cs_test(tree_str, np.zeros((4, 3), int), leaves, [0, 1])


def test_ctmc_null_is_the_default_and_reports_the_exact_one_alongside():
    """The tree-aware null is the default; both p-values must be visible."""
    tree_str = _star_tree(14)
    _, leaves = load_tree(tree_str)
    X = np.zeros((14, 2), dtype=int)
    X[[0, 1, 2]] = 1
    out = archephy_cs_test(tree_str, X, leaves, [0, 1], n_mc=200)
    assert out["method"] == "ctmc"
    assert "p_value_exact_uniform_null" in out
    assert 0.0 < out["p_value"] <= 1.0


def test_ctmc_null_conditions_on_the_observed_change_count():
    """It must differ from the exact null in where the changes fall, not how many.

    An unconditional simulation would let replicates carry more changes than the
    data and inflate the joint count for the wrong reason.
    """
    from amr_clonalshare.archephy import _ctmc_change_sets
    tree_str = _star_tree(20)
    tree, leaves = load_tree(tree_str)
    rng = np.random.default_rng(0)
    sets, matched = _ctmc_change_sets(tree, leaves, 5, 40, rng)
    sizes = [len(s) for s in sets]
    assert matched > 0.5, f"acceptance rate too low: {matched}"
    assert np.mean([abs(s - 5) for s in sizes]) < 1.0


def test_two_identical_archetype_columns_are_collapsed_into_one_measurement():
    """A co-inherited locus pair recorded twice is one measurement, not two.
    Counting it twice put every change edge of the trait into the intersection,
    so m_obs became the trait's own change count and the pair read as
    independent convergent evidence for itself: on this tree it called
    convergence at p = 0.018 where two distinct columns would not."""
    tree_str = _star_tree(12)
    _, leaves = load_tree(tree_str)
    X = np.zeros((12, 3), dtype=int)
    X[[0, 1, 2], 0] = 1
    X[:, 1] = X[:, 0]                           # the same locus, recorded twice
    X[[0, 1, 2, 7], 2] = 1                      # a genuinely distinct partner

    duplicated = archephy_cs_test(tree_str, X, leaves, [0, 1])
    assert duplicated["duplicate_columns_collapsed"] == [[0, 1]]
    assert duplicated["k"] == 1
    assert duplicated["p_value"] > 0.05
    assert duplicated["call"] is False

    distinct = archephy_cs_test(tree_str, X, leaves, [0, 2])
    assert distinct["duplicate_columns_collapsed"] == []
    assert distinct["k"] == 2
    assert distinct["p_value"] < 0.05
    assert distinct["call"] is True
