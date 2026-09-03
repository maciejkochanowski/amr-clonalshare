"""The run report is rendered from the record, in a fixed order, with nothing invented.

The shipped *S. suis* record is the real case; synthetic records exercise the
branches the shipped run does not reach (no structure, withheld inference,
no lineage column, gates tripped).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from amr_clonalshare.cli import _gates, _summary, main
from amr_clonalshare.report import render_report

ROOT = Path(__file__).resolve().parents[1]
SSUIS = ROOT / "out_ssuis" / "cluster_result.json"


@pytest.fixture(scope="module")
def ssuis():
    if not SSUIS.is_file():
        pytest.skip("shipped S. suis record not present")
    return json.loads(SSUIS.read_text(encoding="utf-8"))


def _sections(text):
    return [line[3:] for line in text.splitlines() if line.startswith("## ")]


def test_every_number_in_the_ssuis_report_is_read_from_the_record(ssuis):
    summary = _summary(ssuis)
    text = render_report(ssuis, summary, problems=_gates(ssuis))
    heads = _sections(text)
    assert heads[0] == "1. What was analysed" and heads[1] == "2. Input check"
    assert "Are there resistance-profile groups?" in heads[2]
    assert "How much of it is the clone?" in heads[3]
    assert heads[-1].endswith("What may be concluded")
    assert f"Isolates: {ssuis['n_isolates']}" in text
    assert f"selected {ssuis['selected_k']}" in text
    lam = summary["lineage_attributable_share"]
    if lam is not None:
        assert f"{lam:.3f}" in text
    per = (ssuis.get("metadata_diagnostics") or {}).get("clonal_share") or {}
    if per:
        first = next(iter(per))
        assert f"| {first} |" in text
        assert f"{per[first]['kappa_adj']:.3f}" in text


def test_the_no_structure_branch_and_the_absent_lineage_branch():
    result = {"layers": ["amr"], "n_isolates": 12, "seed": 1, "selected_k": 1,
              "k_selection": {"no_structure": True}}
    text = render_report(result, _summary(result))
    assert "This is a result, not a failure" in text
    assert "No lineage column was supplied" in text
    assert "No diagnostic gate tripped" in text
    assert "No input record was attached" in text


def test_the_withheld_and_confounded_branches_with_an_input_record():
    result = {"layers": ["amr", "vir"], "n_isolates": 50, "seed": 2, "selected_k": 3,
              "provenance": {"version": "1.0.0"},
              "post_clustering_inference": {"status": "withheld_inadequate_split_design"},
              "metadata_diagnostics": {
                  "lineage_attribution": {"lam": 0.8, "ci_low": 0.4, "ci_high": 0.9,
                                          "lambda_gate": 0.5},
                  "clonal_share": {f"t{i}": {"kappa_adj": i / 20, "ci_low": 0.0,
                                             "ci_high": 1.0, "support": 0.95,
                                             "estimable": True}
                                   for i in range(15)}},
              "interpretation": {"claim_level": "descriptive",
                                 "active_gate_codes": ["G1"]}}
    qc = {"missing": {"policy": "drop_rows", "layers": {"amr": {"cells": 3}}},
          "lineage": {"n_groups": 4, "n_singletons": 1, "support": 0.5,
                      "support_threshold": 0.9, "estimable": False}}
    text = render_report(result, _summary(result), input_qc=qc,
                         problems=["fusion collapsed"])
    assert "could not be tested" in text
    assert "lineage-confounded" in text
    assert "could still change this reading" in text
    assert "Version: 1.0.0" in text
    assert "Empty cells in the binary layers: 3, handled by the policy `drop_rows`" in text
    assert "not accepted at this typing resolution" in text
    assert "and 3 more" in text
    assert "Diagnostic failure: fusion collapsed" in text
    assert "Gates active on this run: G1" in text


def test_the_reproducible_and_discrete_branch_and_the_surveillance_section():
    result = {"layers": ["amr"], "n_isolates": 80, "seed": 3, "selected_k": 2,
              "post_clustering_inference": {"status": "ok", "structure_detected": True,
                                            "p_value_structure": 0.0001,
                                            "discreteness": {"verdict": "discrete",
                                                             "discrete_beyond_a_gradient": True}},
              "metadata_diagnostics": {
                  "lineage_attribution": {"lam": 0.2, "ci_low": 0.1, "ci_high": 0.3,
                                          "lambda_gate": 0.5},
                  "clonal_share_by_layer": {"all_clustering_features": {"kappa_adj": 0.31}},
                  "lineage_resolved_prevalence": {
                      "tet": {"status": "ok", "difference_per_lineage_minus_per_isolate": -0.2}},
                  "prevalence_decomposition": {
                      "contrast_column": "period", "levels": ["a", "b"],
                      "per_feature": {"tet": {"status": "ok", "composition": 0.2,
                                              "within_lineage": -0.25, "difference": -0.05,
                                              "composition_ci95": [0.1, 0.3],
                                              "within_lineage_ci95": [-0.4, -0.1]}}}}}
    text = render_report(result, _summary(result))
    assert "are reproducible (p < 0.001)" in text
    assert "separated by gaps" in text
    assert "Clonal share of the whole panel" in text and "0.310" in text
    assert "5. Mix or rate" in text and "6. What may be concluded" in text
    assert "cancel: tet" in text and "Largest gap" in text
    # the gradient branch
    result["post_clustering_inference"]["discreteness"]["discrete_beyond_a_gradient"] = False
    assert "sit on a gradient" in render_report(result, _summary(result))


def test_a_refit_partition_explains_its_missing_interval_and_gates_survive_null():
    result = {"layers": ["amr"], "n_isolates": 200, "seed": 4, "selected_k": 3,
              "metadata_diagnostics": {
                  "lineage_attribution": {"lam": 0.99, "ci_low": None, "ci_high": None,
                                          "lambda_gate": 0.5, "partition_refit": True,
                                          "lam_cv_sd": 0.0006}}}
    problems = _gates(result)  # a JSON round-trip leaves null where NaN was
    assert problems and "nan to nan" in problems[0]
    text = render_report(result, _summary(result), problems=problems)
    assert "not computed to not computed" not in text
    assert "re-fitted inside every cross-validation fold" in text and "0.001" in text


def test_formatting_never_invents_a_number():
    from amr_clonalshare.report import _ci, _f, _pct
    assert _f(None) == "not computed" and _f(float("nan")) == "not computed"
    assert _f("withheld") == "withheld" and _f(0.12345) == "0.123"
    assert _pct(None) == "not computed" and _pct(0.5) == "50.0 %"
    assert _ci(None, 0.2) == "" and _ci(0.1, 0.2) == " (95 % interval 0.100 to 0.200)"
    assert _ci(float("nan"), 0.2) == ""


def test_a_full_run_writes_the_report(planted_cfg, tmp_path):
    cfg, _ = planted_cfg
    out = tmp_path / "res"
    code = main(["--config", str(cfg.config_path), "--results-dir", str(out),
                 "--quiet", "--threads", "1"])
    assert code in (0, 3)
    text = (out / "report.md").read_text(encoding="utf-8")
    assert text.startswith("# amr-clonalshare run report")
    assert "## 2. Input check" in text and "Empty cells in the binary layers: 0" in text
