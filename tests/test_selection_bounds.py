"""Finite-collection prevalence bounds with no missing-data assumptions."""

import json

import numpy as np
import pytest

from amr_clonalshare.missingness import difference_bounds, finite_collection_bounds


def test_bounds_assign_every_missing_value_to_zero_or_one():
    result = finite_collection_bounds([1, 0, np.nan, 1])

    assert result == {
        "status": "ok",
        "analysis_scope": "supplied finite collection",
        "total_count": 4,
        "observed_count": 3,
        "missing_count": 1,
        "positive_count": 2,
        "observed_prevalence": pytest.approx(2 / 3),
        "lower_bound": pytest.approx(1 / 2),
        "upper_bound": pytest.approx(3 / 4),
        "prevalence_bounds": pytest.approx([1 / 2, 3 / 4]),
    }


def test_explicit_frame_counts_unlisted_members_as_missing():
    result = finite_collection_bounds([1, 0, np.nan, 1], total_count=6)

    assert result["total_count"] == 6
    assert result["observed_count"] == 3
    assert result["missing_count"] == 3
    assert result["positive_count"] == 2
    assert result["prevalence_bounds"] == pytest.approx([1 / 3, 5 / 6])


def test_all_missing_values_have_vacuous_but_defined_bounds():
    result = finite_collection_bounds([np.nan, np.nan])

    assert result["status"] == "ok"
    assert result["observed_count"] == 0
    assert result["missing_count"] == 2
    assert result["positive_count"] == 0
    assert result["observed_prevalence"] is None
    assert result["prevalence_bounds"] == [0.0, 1.0]


def test_empty_collection_uses_json_safe_undefined_values():
    result = finite_collection_bounds([])

    assert result["status"] == "empty"
    assert result["total_count"] == 0
    assert result["observed_count"] == 0
    assert result["missing_count"] == 0
    assert result["positive_count"] == 0
    assert result["observed_prevalence"] is None
    assert result["lower_bound"] is None
    assert result["upper_bound"] is None
    assert result["prevalence_bounds"] == [None, None]
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("values", [[0, 0.5, 1], [0, np.inf], [0, -np.inf]])
def test_bounds_reject_values_outside_binary_or_nan(values):
    with pytest.raises(ValueError, match="binary"):
        finite_collection_bounds(values)


@pytest.mark.parametrize("total_count", [-1, 1.5, True, 2])
def test_bounds_reject_invalid_or_contradictory_frame_size(total_count):
    values = [1, 0, np.nan]
    with pytest.raises(ValueError, match="total_count"):
        finite_collection_bounds(values, total_count=total_count)


def test_difference_bounds_take_the_extreme_cross_combinations():
    a = finite_collection_bounds([1, 1, np.nan, 0])
    b = finite_collection_bounds([1, np.nan, np.nan, 0])

    result = difference_bounds(a, b)

    assert result == {
        "status": "ok",
        "analysis_scope": "difference between supplied finite collections",
        "lower_bound": pytest.approx(-0.25),
        "upper_bound": pytest.approx(0.5),
        "difference_bounds": pytest.approx([-0.25, 0.5]),
    }


def test_difference_bounds_propagate_an_empty_collection_without_nan():
    result = difference_bounds(finite_collection_bounds([]),
                               finite_collection_bounds([0, 1]))

    assert result["status"] == "empty"
    assert result["lower_bound"] is None
    assert result["upper_bound"] is None
    assert result["difference_bounds"] == [None, None]
    json.dumps(result, allow_nan=False)


def test_none_counts_as_a_missing_outcome():
    from amr_clonalshare.missingness import finite_collection_bounds
    bounds = finite_collection_bounds([1, 0, None], total_count=5)
    assert bounds['observed_count'] == 2 and bounds['missing_count'] == 3
