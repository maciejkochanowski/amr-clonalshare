"""Missing atlas diagnostics and hostile report labels remain portable text."""
import importlib
import json
from pathlib import Path
import sys

import numpy as np

from amr_clonalshare.cli import _summary
from amr_clonalshare.report import render_report

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
atlas = importlib.import_module("atlas_cross_species")


def _strict(text):
    def reject(value):
        raise ValueError(value)
    return json.loads(text, parse_constant=reject)


def test_atlas_serializes_unavailable_realised_diagnostics_as_null(tmp_path):
    record = {"organism": "synthetic", "n_isolates": 60, "n_lineages": 3,
              "support": 1., "n_drugs_available": 1, "agents": {"drug": {
                  "kappa": .4, "kappa_permuted": .0, "prevalence": .5,
                  "ci_low": .1, "ci_high": .6, "support": 1., "estimable": True,
                  "permuted_ci": [-.1, .1], "realised": {
                      "estimable": False, "kappa": np.nan,
                      "ci_low": np.nan, "ci_high": np.nan,
                      "reason": "insufficient groups"}}}}
    atlas.write_outputs(tmp_path, [record])
    table = _strict((tmp_path / "atlas_table.json").read_text(encoding="utf-8"))
    assert table[0]["realised_kappa"] is None
    assert table[0]["realised_reason"] == "insufficient groups"
    _strict((tmp_path / "atlas_summary.json").read_text(encoding="utf-8"))


def test_atlas_all_refused_can_be_aggregated(tmp_path):
    atlas.write_outputs(tmp_path, [{"organism": "synthetic", "error": "too few joined isolates"}])
    summary = _strict((tmp_path / "atlas_summary.json").read_text(encoding="utf-8"))
    assert summary["admitted"]["cells"] == 0
    assert summary["admitted"]["kappa_median"] is None


def test_one_lineage_refusal_keeps_screen_type_contract():
    screen = importlib.import_module("agent_screen").screen_agent(
        [0, 1] * 30, ["single"] * 60)
    assert not screen.analyse
    assert isinstance(screen.y, np.ndarray)
    assert screen.y.size == 0


def test_markdown_escapes_raw_html_labels_in_tables_and_identity():
    label = '<img src=x onerror="alert(1)">'
    record = {"n_isolates": 60, "seed": 1, "metadata_diagnostics": {
        "lineage_column": label, "clonal_share": {label: {
            "kappa_adj": .4, "ci_low": .2, "ci_high": .6,
            "n_groups": 3, "support": 1., "estimable": True}}}}
    report = render_report(record, _summary(record))
    assert label not in report
    assert "&lt;img" in report
    assert "| &lt;img" in report
    assert "> **" in report  # authored Markdown callouts retain their structure


def test_input_check_markdown_escapes_trait_labels():
    import pandas as pd
    from amr_clonalshare.qc import input_qc, render_markdown
    label = "<img src=x onerror=alert(1)>"
    qc = input_qc(pd.DataFrame({label: [0, 0, 0, 0]}))
    report = render_markdown(qc)
    assert label not in report
    assert "&lt;img" in report
