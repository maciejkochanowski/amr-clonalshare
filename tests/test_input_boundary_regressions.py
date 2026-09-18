"""Regression tests for the 1.0.0 input contract."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import yaml

from amr_clonalshare import core, qc
from amr_clonalshare.censored import intervals_from_mic
from amr_clonalshare.cli import _summary, main
from amr_clonalshare.config import ConfigError, from_dict, load_config
from amr_clonalshare.io import load_dataset
from amr_clonalshare.report import render_report
from amr_clonalshare.report_html import render_html_report
from amr_clonalshare.report_model import build_report
from conftest import planted_cohort, write_config


@pytest.mark.parametrize("value", [False, 0, [], ""])
def test_non_mapping_optional_section_is_refused(tmp_path, value):
    raw = planted_cohort(tmp_path)
    raw["evidence"] = value
    with pytest.raises(ConfigError, match="mapping"):
        from_dict(raw)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), "-inf"])
def test_nonfinite_integer_budget_is_a_config_error(tmp_path, value):
    raw = planted_cohort(tmp_path)
    raw["attribution"]["n_boot"] = value
    with pytest.raises(ConfigError, match="whole number"):
        from_dict(raw)


def test_identical_contrast_levels_are_refused(tmp_path):
    raw = planted_cohort(tmp_path)
    raw["dataset"].update(contrast_column="intake", contrast_levels=["2010", "2010"])
    with pytest.raises(ConfigError, match="distinct"):
        from_dict(raw).validate(check_files_exist=False)


def test_duplicate_yaml_keys_are_refused(tmp_path):
    raw = planted_cohort(tmp_path)
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(raw) + "evidence:\n  enabled: false\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="duplicate"):
        load_config(path)


def test_skip_precheck_still_refuses_a_missing_required_table(tmp_path, capsys):
    raw = planted_cohort(tmp_path)
    raw["dataset"]["metadata"] = "absent.csv"
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    assert main(["--config", str(path), "--no-check-files", "--check-input"]) == 2
    assert "input error" in capsys.readouterr().err


def test_unreadable_calls_are_accounted_for_in_input_qc(tmp_path):
    raw = planted_cohort(tmp_path)
    path = tmp_path / "data" / "calls.csv"
    table = pd.read_csv(path)
    table.loc[table.iid == table.iid.iloc[0], "call"] = "unreadable"
    table.to_csv(path, index=False)
    _, cfg = write_config(tmp_path, raw)
    with pytest.warns(RuntimeWarning, match="unrecognized"):
        ds = load_dataset(cfg)
    counts = ds.input_qc["phenotype_reading"]
    assert counts["rows_read"] == 288
    assert counts["n_unrecognized_calls"] == 3
    assert counts["isolates_without_readable_calls"] == 1
    assert counts["unrecognized_calls"] == {"unreadable": 3}
    assert "3 unrecognized" in qc.render_markdown(ds.input_qc)


def test_single_lineage_is_not_declared_estimable_by_input_check():
    assert qc.group_adequacy(pd.Series(["A"] * 24))["estimable"] is False


def test_single_lineage_input_explains_the_missing_comparison():
    panel = pd.DataFrame({"agent": [0, 1] * 12})
    metadata = pd.DataFrame({"lineage": ["A"] * 24})
    report = qc.input_qc(panel, metadata=metadata, lineage_column="lineage")
    text = qc.render_markdown(report)
    assert "two distinct lineages" in text
    assert "nan" not in text


def test_inline_censoring_survives_an_empty_operator_cell(tmp_path):
    raw = planted_cohort(tmp_path, mic=True)
    raw["dataset"]["mic_operator_column"] = "operator"
    path = tmp_path / "data" / "mic.csv"
    table = pd.read_csv(path)
    table["measurement"] = table.measurement.astype(str)
    table.loc[0, "measurement"] = "<=8"
    table["operator"] = ""
    table.to_csv(path, index=False)
    _, cfg = write_config(tmp_path, raw)
    ds = load_dataset(cfg)
    lo, hi = intervals_from_mic(ds.mic.measurement, operators=ds.mic.operator)
    assert np.isneginf(lo[0]) and hi[0] == 3.0


@pytest.mark.parametrize("operator", ["nonsense", "=<"])
def test_invalid_censoring_operator_is_refused(operator):
    with pytest.raises(ValueError, match="operator"):
        intervals_from_mic([1, 2, 4], operators=["", operator, ""])


def test_explicit_exact_operator_overrides_end_well_heuristic():
    lo, hi = intervals_from_mic([1, 2, 4], operators=["=", "=", "="])
    np.testing.assert_array_equal(lo, [-1.0, 0.0, 1.0])
    np.testing.assert_array_equal(hi, [0.0, 1.0, 2.0])


def test_no_bootstrap_keeps_the_point_estimate_in_both_reports(tmp_path):
    raw = planted_cohort(tmp_path)
    raw["attribution"].update(n_boot=0, n_perm=5, repeats=2)
    raw["evidence"]["enabled"] = False
    raw["surveillance"] = {"enabled": False}
    _, cfg = write_config(tmp_path, raw)
    record = core.run(cfg, seed=4)
    summary = _summary(record)
    build_report(record, summary, input_qc=record["input_qc"])
    # Schema 2 reports each analysis independently; a realised-component
    # refusal must not hide a supported membership point estimate.
    assert summary["analyses"]["collection_membership"]["n_computed"] == 3
    rows = [r for r in record["results"] if r["analysis"] == "collection_membership"]
    assert len(rows) == 3
    assert all(r["status"] == "computed" and r["estimate"] is not None
               and r["lower"] is None and r["upper"] is None for r in rows)
    for render in (render_report, render_html_report):
        text = render(record, summary, input_qc=record["input_qc"])
        assert "interval not computed" in text.lower()
        for name, result in record["metadata_diagnostics"]["clonal_share"].items():
            assert name in text and f"{result['kappa_adj']:.3f}" in text


@pytest.mark.parametrize("file,column", [("calls.csv", "iid"), ("calls.csv", "antibiotic"),
                                        ("meta.csv", "iid"), ("mic.csv", "iid"),
                                        ("mic.csv", "antibiotic")])
def test_missing_table_keys_are_refused_before_joining(tmp_path, file, column):
    raw = planted_cohort(tmp_path, mic=True)
    path = tmp_path / "data" / file
    table = pd.read_csv(path)
    table.loc[0, column] = ""
    table.to_csv(path, index=False)
    _, cfg = write_config(tmp_path, raw)
    with pytest.raises(ConfigError, match="missing.*" + column):
        load_dataset(cfg)


def test_markdown_never_emits_input_labels_as_raw_html():
    payload = '<img src=x onerror="alert(1)">'
    record = {"n_isolates": 30, "n_traits": 1, "seed": 1,
              "metadata_diagnostics": {"lineage_column": payload,
                  "clonal_share": {payload: {"kappa_adj": 0.4, "ci_low": 0.2,
                      "ci_high": 0.6, "null_mean": 0.0, "support": 1.0,
                      "estimable": True, "n_groups": 6}}}}
    text = render_report(record, _summary(record))
    assert "<img" not in text
    assert "&lt;img" in text
