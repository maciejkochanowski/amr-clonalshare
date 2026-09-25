"""Assumption-free bounds for binary outcomes in a known finite collection.

The functions in this module do not impute missing outcomes and do not attach
confidence levels. They describe only the supplied finite frame. For ``r``
observed positives among ``n`` observed members of a frame of size ``N``, the
prevalence is bounded exactly by ``r / N`` and ``(r + N - n) / N``: the two
endpoints assign every missing outcome to zero or one, respectively.
"""
from __future__ import annotations

import math
from numbers import Integral, Real
from typing import Any, Dict, Mapping, Optional, Sequence

import numpy as np

__all__ = ["finite_collection_bounds", "difference_bounds"]


def _frame_size(total_count: Optional[int], supplied_count: int) -> int:
    if total_count is None:
        return supplied_count
    if isinstance(total_count, bool) or not isinstance(total_count, Integral):
        raise ValueError("total_count must be a non-negative integer")
    frame_size = int(total_count)
    if frame_size < supplied_count:
        raise ValueError(
            "total_count cannot be smaller than the number of supplied values")
    return frame_size


def finite_collection_bounds(
    values: Sequence[Any], *, total_count: Optional[int] = None
) -> Dict[str, Any]:
    """Bound binary prevalence in the supplied finite collection.

    ``values`` contains numeric zeros, ones, and NaNs. ``total_count`` may name
    a larger known frame; members absent from ``values`` then count as missing.
    Finite values other than zero and one, infinities, invalid frame sizes, and
    frame sizes smaller than the supplied sequence are rejected.

    Returned undefined quantities are ``None`` so the record is valid strict
    JSON. These are identification bounds for this finite frame, not confidence
    intervals and not an extrapolation to a population beyond the frame.
    """
    supplied = list(values)
    frame_size = _frame_size(total_count, len(supplied))

    observed_count = 0
    positive_count = 0
    explicit_missing = 0
    for value in supplied:
        if value is None:
            explicit_missing += 1
            continue
        if not isinstance(value, (Real, np.bool_)):
            raise ValueError("values must be binary numeric 0/1 or NaN")
        numeric = float(value)
        if math.isnan(numeric):
            explicit_missing += 1
        elif not math.isfinite(numeric) or numeric not in (0.0, 1.0):
            raise ValueError("values must be binary numeric 0/1 or NaN")
        else:
            observed_count += 1
            positive_count += int(numeric)

    missing_count = explicit_missing + frame_size - len(supplied)
    if frame_size == 0:
        lower = upper = observed_prevalence = None
        status = "empty"
    else:
        lower = float(positive_count / frame_size)
        upper = float((positive_count + missing_count) / frame_size)
        observed_prevalence = (
            float(positive_count / observed_count)
            if observed_count else None
        )
        status = "ok"

    return {
        "status": status,
        "analysis_scope": "supplied finite collection",
        "total_count": int(frame_size),
        "observed_count": int(observed_count),
        "missing_count": int(missing_count),
        "positive_count": int(positive_count),
        "observed_prevalence": observed_prevalence,
        "lower_bound": lower,
        "upper_bound": upper,
        "prevalence_bounds": [lower, upper],
    }


def difference_bounds(
    bounds_a: Mapping[str, Any], bounds_b: Mapping[str, Any]
) -> Dict[str, Any]:
    """Derive exact bounds for finite-collection prevalence A minus B.

    Each argument is a record returned by :func:`finite_collection_bounds`.
    If either collection is empty, the difference is undefined and represented
    by ``None`` values rather than non-standard JSON NaNs.
    """
    if bounds_a.get("status") == "empty" or bounds_b.get("status") == "empty":
        lower = upper = None
        status = "empty"
    else:
        try:
            lower_a = float(bounds_a["lower_bound"])
            upper_a = float(bounds_a["upper_bound"])
            lower_b = float(bounds_b["lower_bound"])
            upper_b = float(bounds_b["upper_bound"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("difference_bounds requires two non-empty bounds records") from exc
        if not all(math.isfinite(value)
                   for value in (lower_a, upper_a, lower_b, upper_b)):
            raise ValueError("difference_bounds requires finite bound endpoints")
        lower = float(lower_a - upper_b)
        upper = float(upper_a - lower_b)
        status = "ok"

    return {
        "status": status,
        "analysis_scope": "difference between supplied finite collections",
        "lower_bound": lower,
        "upper_bound": upper,
        "difference_bounds": [lower, upper],
    }

