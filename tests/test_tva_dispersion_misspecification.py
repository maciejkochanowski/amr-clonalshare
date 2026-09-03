"""Counterexamples for the pooled-dispersion failure found by the null calibration.

A pooled dispersion estimate is only valid for thinning when the features really do
share one dispersion. When they do not, the split leaves dependence behind and the
separation test reads it as cluster structure: the measured null rejection rate reaches
0.95 against a nominal 0.05 (``audit_checks/calib_pooled_heterogeneity_sweep.py``).

These tests pin the refusal, the threshold it is set at, and the fact that per-feature
estimation — the default — is not affected.
"""

from __future__ import annotations

import numpy as np
import pytest

from amr_clonalshare.tva import (
    POOLED_DISPERSION_MAX_ESTIMATED_RATIO,
    pooled_dispersion_spread,
    tva_test_separation,
)


def _homogeneous(n=200, p=10, r=5.0, mu=8.0, seed=0):
    rng = np.random.default_rng(seed)
    return rng.negative_binomial(r, r / (r + mu), size=(n, p))


def _heterogeneous(n=200, p=10, r_low=0.5, r_high=50.0, mu=8.0, seed=0):
    """Structure-free counts whose per-feature dispersion spans two orders of magnitude."""
    rng = np.random.default_rng(seed)
    r_j = np.exp(rng.uniform(np.log(r_low), np.log(r_high), size=p))
    X = np.empty((n, p), dtype=np.int64)
    for j in range(p):
        X[:, j] = rng.negative_binomial(r_j[j], r_j[j] / (r_j[j] + mu), size=n)
    return X


def test_pooling_is_refused_when_feature_dispersions_genuinely_differ():
    X = _heterogeneous(seed=1)
    assert pooled_dispersion_spread(X) > POOLED_DISPERSION_MAX_ESTIMATED_RATIO
    with pytest.raises(ValueError, match="pools one value across features"):
        tva_test_separation(X, r=None, k=2, dispersion="pooled_mom",
                            rng=np.random.default_rng(1), n_splits=3)


def test_pooling_still_works_when_the_features_share_a_dispersion():
    X = _homogeneous(seed=2)
    assert pooled_dispersion_spread(X) <= POOLED_DISPERSION_MAX_ESTIMATED_RATIO
    result = tva_test_separation(X, r=None, k=2, dispersion="pooled_mom",
                                 rng=np.random.default_rng(2), n_splits=3)
    assert 0.0 <= result["p_value"] <= 1.0
    assert len(set(result["r"])) == 1


def test_per_feature_estimation_accepts_the_data_that_pooling_refuses():
    # The refusal must be specific to pooling; the default path is valid here and the
    # sweep measured it at nominal across the whole heterogeneity range.
    X = _heterogeneous(seed=3)
    result = tva_test_separation(X, r=None, k=2, dispersion="mom",
                                 rng=np.random.default_rng(3), n_splits=3)
    assert 0.0 <= result["p_value"] <= 1.0
    assert len(set(result["r"])) > 1


def test_an_explicitly_supplied_dispersion_is_never_second_guessed():
    # Passing r is a caller assertion about the generating model, e.g. a simulation
    # oracle. The guard applies to estimation, not to a supplied truth.
    X = _heterogeneous(seed=4)
    result = tva_test_separation(X, r=np.full(X.shape[1], 5.0), k=2,
                                 dispersion="pooled_mom",
                                 rng=np.random.default_rng(4), n_splits=3)
    assert result["r"] == [5.0] * X.shape[1]


def test_an_unknown_dispersion_name_is_refused_rather_than_silently_defaulted():
    # A typo previously fell through to the per-feature default, so a result could name
    # an estimator that never ran.
    X = _homogeneous(seed=5)
    with pytest.raises(ValueError, match="unknown dispersion"):
        tva_test_separation(X, r=None, k=2, dispersion="pooled",
                            rng=np.random.default_rng(5), n_splits=3)


def test_the_guard_threshold_matches_the_calibration_it_was_taken_from():
    # Between the last safe level (true ratio 2, estimated median 2.25, rejection 0.050)
    # and the first unsafe one (true ratio 4, estimated median 3.58, rejection 0.250).
    assert 2.25 < POOLED_DISPERSION_MAX_ESTIMATED_RATIO < 3.58
