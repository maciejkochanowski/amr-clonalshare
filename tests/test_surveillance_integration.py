"""Config surface and summary reading for the lineage-aware surveillance block.

The decomposition is only run when a contrast is declared, so the config has to
refuse a half-declared one loudly rather than silently skip the analysis: a
misspelt level or a missing lineage column would otherwise look identical to
"not requested".

The summary logic is tested on hand-built records rather than on a pipeline run,
because the reading that matters - an agent whose components offset - is
precisely the one a real cohort may or may not contain.
"""
from __future__ import annotations

import pytest

from amr_clonalshare.cli import _surveillance_summary
from amr_clonalshare.config import ConfigError, DatasetConfig


def _dataset(**kw) -> DatasetConfig:
    base = dict(name="d", metadata="m.csv", lineage_column="mlst")
    base.update(kw)
    return DatasetConfig(**base)


# --- config surface ---------------------------------------------------------

def test_a_contrast_without_metadata_is_refused():
    with pytest.raises(ConfigError, match="metadata"):
        _dataset(metadata=None, lineage_column=None,
                 contrast_column="country",
                 contrast_levels=("a", "b")).validate()


def test_a_contrast_without_a_lineage_column_is_refused():
    with pytest.raises(ConfigError, match="lineage_column"):
        _dataset(lineage_column=None, contrast_column="country",
                 contrast_levels=("a", "b")).validate()


@pytest.mark.parametrize("levels", [(), ("only",), ("a", "b", "c")])
def test_a_contrast_needs_exactly_two_levels(levels):
    with pytest.raises(ConfigError, match="exactly two"):
        _dataset(contrast_column="country", contrast_levels=levels).validate()


def test_a_complete_contrast_validates():
    _dataset(contrast_column="country", contrast_levels=("a", "b")).validate()


def test_no_contrast_is_not_an_error():
    _dataset().validate()


# --- summary reading --------------------------------------------------------

def _decomp(diff, comp, within, comp_ci, within_ci):
    return {"status": "ok", "difference": diff, "composition": comp,
            "within_lineage": within, "composition_ci95": comp_ci,
            "within_lineage_ci95": within_ci}


def test_offsetting_components_are_singled_out():
    # `flat` is the case the decomposition exists for: prevalence moved by a
    # point while both components moved by twelve, in opposite directions.
    meta = {"prevalence_decomposition": {
        "contrast_column": "period", "levels": ["late", "early"],
        "per_feature": {
            "flat": _decomp(-0.006, -0.120, +0.114, [-0.21, -0.01], [0.02, 0.19]),
            "real_decline": _decomp(-0.110, -0.101, -0.008,
                                    [-0.16, -0.03], [-0.03, 0.01]),
            "nothing": _decomp(0.028, -0.039, 0.067, [-0.11, 0.04], [-0.01, 0.14]),
        }}}
    out = _surveillance_summary(meta)["decomposition"]
    assert out["n_features"] == 3
    assert out["n_composition_significant"] == 2
    assert out["n_within_lineage_significant"] == 1
    assert out["offsetting_features"] == ["flat"]


def test_a_real_decline_is_not_called_offsetting():
    # both components negative: they add up, they do not cancel
    meta = {"prevalence_decomposition": {
        "contrast_column": "period", "levels": ["late", "early"],
        "per_feature": {"x": _decomp(-0.20, -0.12, -0.08,
                                     [-0.18, -0.05], [-0.14, -0.02])}}}
    assert _surveillance_summary(meta)["decomposition"]["n_offsetting"] == 0


def test_widest_prevalence_gap_is_reported_with_its_feature():
    meta = {"lineage_resolved_prevalence": {
        "a": {"status": "ok", "difference_per_lineage_minus_per_isolate": 0.02},
        "b": {"status": "ok", "difference_per_lineage_minus_per_isolate": 0.21},
        "c": {"status": "skipped"},
    }}
    out = _surveillance_summary(meta)
    assert out["lineage_prevalence_widest_gap_feature"] == "b"
    assert out["lineage_prevalence_widest_gap"] == pytest.approx(0.21)


def test_carriage_directions_are_counted_and_significance_separated():
    meta = {"trait_concentration": {
        "a": {"status": "ok", "direction": "dispersed", "p_value": 0.0005},
        "b": {"status": "ok", "direction": "dispersed", "p_value": 0.30},
        "c": {"status": "ok", "direction": "proportional", "p_value": 0.40},
        "d": {"status": "skipped"},
    }}
    out = _surveillance_summary(meta)
    assert out["carriage_direction_counts"] == {"dispersed": 2,
                                                "proportional": 1}
    assert out["n_features_departing_from_proportional_carriage"] == 1


def test_an_empty_block_yields_an_empty_reading():
    assert _surveillance_summary({}) == {}
