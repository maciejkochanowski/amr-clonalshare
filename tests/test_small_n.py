"""The gap statistic, with the same clusterer on both sides.

The published rule is a difference of two dispersions produced by the same
estimator, so the tests use one k-means labeller for the observed data and
for every reference sample. A planted three-group panel must select three and
a panel with no structure must select one; an elbow that does not fire must
be reported as not firing rather than replaced by an argmax.
"""
from __future__ import annotations

import numpy as np
from sklearn.cluster import KMeans

from amr_clonalshare.small_n import compute_Wk, gap_statistic_k
from amr_clonalshare.synthetic import synth_cluster_archetypes


def _kmeans(X, k, seed=0):
    if k == 1:
        return np.zeros(X.shape[0], dtype=int)
    return KMeans(n_clusters=k, n_init=5, random_state=seed).fit_predict(X)


def _select(X, k_range, B=15):
    return gap_statistic_k(X, k_range,
                           observed_labeler=lambda k: _kmeans(X, k),
                           null_labeler=lambda Xn, k: _kmeans(Xn, k),
                           B=B, seed=3)


def test_planted_three_groups_select_three():
    amr, vir, _ = synth_cluster_archetypes(n=150, p_amr=20, p_vir=20, k_true=3,
                                           overlap=0.05, seed=11)
    X = np.hstack([amr.to_numpy(float), vir.to_numpy(float)])
    out = _select(X, [2, 3, 4, 5])
    assert out["rule_fired"]
    assert out["k_star"] == 3
    assert [row["k"] for row in out["table"]] == [1, 2, 3, 4, 5]


def test_no_structure_selects_one():
    rng = np.random.default_rng(5)
    X = (rng.random((120, 30)) < rng.uniform(0.2, 0.6, size=30)).astype(float)
    out = _select(X, [2, 3, 4])
    assert out["k_star"] == 1
    assert out["rule_fired"]


def test_an_elbow_that_never_fires_is_reported_not_replaced():
    """With a single candidate the one-standard-error rule has no k + 1 to
    compare against, so it cannot fire; the result must say so and return
    the largest k asked for, not an argmax the paper does not define."""
    rng = np.random.default_rng(6)
    X = (rng.random((60, 12)) < 0.4).astype(float)
    out = gap_statistic_k(X, [2], observed_labeler=lambda k: _kmeans(X, k),
                          null_labeler=lambda Xn, k: _kmeans(Xn, k),
                          B=5, seed=1, include_k1=False)
    assert out["rule_fired"] is False
    assert out["k_star"] == 2


def test_within_dispersion_is_zero_for_identical_rows_and_grows_with_disagreement():
    X = np.array([[1, 0, 1, 0]] * 6, dtype=float)
    assert compute_Wk(X, np.zeros(6, dtype=int)) == 0.0
    Y = X.copy()
    Y[:3, 0] = 0
    assert compute_Wk(Y, np.zeros(6, dtype=int)) > 0.0
    # Splitting the two kinds of row into their own clusters removes the
    # disagreement inside every cluster.
    assert compute_Wk(Y, np.array([0, 0, 0, 1, 1, 1])) == 0.0
