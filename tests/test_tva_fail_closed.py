"""Fail-closed contracts for public negative-binomial TVA reporting."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from amr_clonalshare.tva import tva_report, tva_test_separation


def _nb_counts(*, zero_inflation: float = 0.0) -> pd.DataFrame:
    rng = np.random.default_rng(20260811)
    n, p = 120, 4
    r = np.full(p, 5.0)
    mean = np.full((n, p), 8.0)
    counts = rng.negative_binomial(r[None, :], r[None, :] / (r[None, :] + mean))
    if zero_inflation:
        counts[rng.random(counts.shape) < zero_inflation] = 0
    return pd.DataFrame(counts, columns=[f"count_{j}" for j in range(p)])


@pytest.mark.parametrize("zero_inflation", [0.0, 0.25])
def test_estimated_dispersion_is_exploratory_even_when_screening_passes(
    zero_inflation,
):
    report = tva_report(
        _nb_counts(zero_inflation=zero_inflation),
        k=2,
        rng=np.random.default_rng(7),
        n_splits=1,
        n_init=2,
    )

    assert report["screening"]["n_kept"] == 4
    assert report["status"] == "exploratory_only"
    assert report["dispersion_parameter_source"] == "estimated_from_same_matrix"
    assert report["nominal_inference_valid"] is False
    assert report["archetypes_real"] is None
    assert report["p_value_separation"] is None
    assert report["per_feature"] is None
    assert report["n_defining_thinned"] is None
    assert 0.0 <= report["diagnostic_p_value_separation"] <= 1.0
    assert (
        report["diagnostic_separation"]["p_value"]
        == (report["diagnostic_p_value_separation"])
    )


def test_supplied_r_without_source_marker_is_not_nominal():
    report = tva_report(
        _nb_counts(),
        r=np.full(4, 5.0),
        k=2,
        rng=np.random.default_rng(8),
        n_splits=1,
        n_init=2,
    )

    assert report["status"] == "exploratory_only"
    assert report["dispersion_parameter_source"] == (
        "supplied_without_validated_source"
    )
    assert report["nominal_inference_valid"] is False
    assert report["archetypes_real"] is None
    assert report["p_value_separation"] is None


def test_oracle_marker_is_required_for_nominal_simulation_report():
    report = tva_report(
        _nb_counts(),
        r=np.full(4, 5.0),
        r_source="simulation_oracle",
        k=2,
        rng=np.random.default_rng(9),
        n_splits=1,
        n_init=2,
    )

    assert report["status"] == "ok"
    assert report["dispersion_parameter_source"] == "simulation_oracle"
    assert report["nominal_inference_valid"] is True
    assert isinstance(report["archetypes_real"], bool)
    assert 0.0 <= report["p_value_separation"] <= 1.0
    assert "diagnostic_p_value_separation" not in report


def test_source_marker_cannot_be_attached_to_implicitly_estimated_r():
    with pytest.raises(ValueError, match="requires an explicitly supplied r"):
        tva_report(
            _nb_counts(),
            r_source="external_validated",
            k=2,
            rng=np.random.default_rng(10),
            n_splits=1,
            n_init=2,
        )


def test_oracle_r_vector_tracks_columns_retained_by_screening():
    counts = _nb_counts()
    counts["count_1"] = 1
    supplied = np.array([2.0, 3.0, 4.0, 5.0])

    report = tva_report(
        counts,
        r=supplied,
        r_source="simulation_oracle",
        k=2,
        rng=np.random.default_rng(11),
        n_splits=1,
        n_init=2,
    )

    assert report["status"] == "ok"
    assert report["counts_columns_used"] == ["count_0", "count_2", "count_3"]
    assert report["r_per_feature"] == {
        "count_0": 2.0,
        "count_2": 4.0,
        "count_3": 5.0,
    }


def test_the_separation_test_screens_for_thinnability_before_it_splits_anything():
    """The module docstring credits ``screen_thinnable`` with taking the
    marginal rejection rate from 0.145 to 0.044, but ``tva_test_separation``
    never consulted it: eighty identical rows, with no variance to split at
    all, came back at p = 3e-14. It must refuse, name the reason, and still
    test the matrices that are genuinely thinnable."""
    identical = np.tile(np.array([5, 9, 2, 14, 7, 3, 11, 6]), (80, 1))
    refused = tva_test_separation(identical, k=2, n_splits=5, n_init=5,
                                  rng=np.random.default_rng(3))
    assert refused["status"] == "not_thinnable"
    assert np.isnan(refused["p_value"])
    assert refused["screening"]["n_kept"] == 0
    assert refused["p_per_split"] == []
    assert "NB-thinnable" in refused["reason"]

    draw = np.random.default_rng(11)
    counts = np.vstack([
        draw.negative_binomial(4, 4 / (4 + 12.0), size=(40, 6)),
        draw.negative_binomial(4, 4 / (4 + 40.0), size=(40, 6)),
    ]).astype(np.int64)
    real = tva_test_separation(counts, k=2, n_splits=5, n_init=5,
                               rng=np.random.default_rng(3))
    assert real["status"] == "ok"
    assert real["screening"]["n_kept"] == 6
    assert real["p_value"] < 0.01
