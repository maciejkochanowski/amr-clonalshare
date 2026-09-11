"""A cohort is read for the share its lineages carry, and nothing else.

The run reaches every estimator through one function and writes one record
whose shape is fixed by the schema.
"""
from __future__ import annotations

import json

import pytest

from amr_clonalshare import run
from amr_clonalshare.config import ConfigError, from_dict
from conftest import planted_cohort, write_config


def test_a_config_naming_no_phenotype_table_is_refused():
    with pytest.raises(ConfigError) as exc:
        from_dict({"dataset": {"name": "nothing"}}).validate(check_files_exist=False)
    assert "phenotype" in str(exc.value)


def test_a_config_without_a_lineage_column_is_refused(tmp_path):
    raw = planted_cohort(tmp_path)
    del raw["dataset"]["lineage_column"]
    with pytest.raises(ConfigError, match="lineage_column"):
        write_config(tmp_path, raw)


def test_a_cohort_is_read_for_the_share(share_cfg, tmp_path):
    _, cfg = share_cfg
    out = run(cfg, results_dir=tmp_path / "out", seed=7)

    assert out["n_isolates"] == 96
    assert sorted(out["traits"]) == ["agentA", "agentB", "agentC"]
    assert (tmp_path / "out" / "clonal_share_result.json").is_file()

    share = out["metadata_diagnostics"]["clonal_share"]
    assert set(share) == {"agentA", "agentB", "agentC"}
    # The two planted agents follow the lineage and the third does not, so the
    # share separates them and the permuted control stays near zero on all
    # three. This is the only claim the test makes about the numbers.
    assert share["agentA"]["kappa_adj"] > 0.25
    assert share["agentB"]["kappa_adj"] > 0.25
    assert share["agentC"]["kappa_adj"] < 0.15
    for name, d in share.items():
        # A permuted labelling scores below zero out of sample and the margin
        # grows with the number of lineages; what the control has to show is
        # that it manufactures no share, not that it lands on zero.
        assert -0.15 < d["null_mean"] < 0.05
        assert d["support"] == pytest.approx(1.0)
    for name in ("agentA", "agentB"):
        assert (share[name]["kappa_adj"] - share[name]["null_mean"]) > 0.30


def test_the_record_carries_exactly_the_documented_blocks(share_cfg):
    _, cfg = share_cfg
    out = run(cfg, results_dir=None, seed=7)
    assert set(out) == {"schema_version", "seed", "config", "versions",
                        "n_isolates", "n_traits", "traits",
                        "metadata_diagnostics", "input_qc"}
    md = out["metadata_diagnostics"]
    assert set(md) <= {"lineage_column", "clonal_share", "realised_share",
                       "lineage_evidence", "censored_share",
                       "lineage_resolved_prevalence", "trait_concentration",
                       "prevalence_decomposition"}


def test_the_run_carries_the_evidence_and_the_intakes(share_cfg):
    _, cfg = share_cfg
    out = run(cfg, results_dir=None, seed=7)
    ev = out["metadata_diagnostics"]["lineage_evidence"]
    assert set(ev["per_feature"]) == {"agentA", "agentB", "agentC"}
    assert "e_bh" in ev
    seq = ev["sequential"]
    assert seq["batch_column"] == "intake"
    assert seq["batches"] == ["2010", "2011", "2012"]
    assert "e_bh" in seq
    one = next(iter(seq["per_feature"].values()))
    # The first intake has nothing before it, so it earns a factor of one and
    # the running log starts at zero; the product may only be read forward.
    assert one["log_e"][0] == 0.0
    assert len(one["log_e"]) == len(seq["batches"])


def test_the_same_seed_gives_the_same_record(share_cfg):
    _, cfg = share_cfg
    a = run(cfg, seed=3)
    b = run(cfg, seed=3)
    assert json.dumps(a["metadata_diagnostics"], sort_keys=True) == \
        json.dumps(b["metadata_diagnostics"], sort_keys=True)


def test_the_cli_writes_the_record_and_the_two_reports(share_cfg, tmp_path):
    """The command in the article, on a configuration that names a panel and a
    lineage column, writes the record and both reports and exits zero."""
    from amr_clonalshare import cli

    cfg_path, _ = share_cfg
    out = tmp_path / "out"
    assert cli.main(["--config", str(cfg_path), "--results-dir", str(out),
                     "--seed", "7", "--quiet"]) == 0

    record = json.loads((out / "clonal_share_result.json").read_text())
    assert record["n_isolates"] == 96
    for name in ("report.md", "report.html", "input_qc.json", "input_qc.md"):
        assert (out / name).is_file(), name
    assert not (out / "cluster_result.json").exists()
    summary = cli._summary(record)
    assert summary["n_estimable"] == 3
    assert set(summary["per_agent"]) == {"agentA", "agentB", "agentC"}
    assert "sha256:" in (out / "report.md").read_text()


def test_a_dilution_table_is_read_beside_the_calls(tmp_path):
    from amr_clonalshare.cli import _summary
    from amr_clonalshare.report import render_report

    raw = planted_cohort(tmp_path, mic=True)
    _, cfg = write_config(tmp_path, raw)
    out = run(cfg, seed=7)
    cz = out["metadata_diagnostics"]["censored_share"]
    assert cz["join"]["strains_with_mic"] == 96
    assert set(cz["per_agent"]) == {"agentA", "agentB", "agentC"}
    entry = cz["per_agent"]["agentA"]
    assert entry["n"] == 96 and "panel" in entry and "share" in entry
    assert out["input_qc"]["mic_join"]["rows_joined"] == 288
    text = render_report(out, _summary(out), input_qc=out["input_qc"])
    assert "Reading at the recorded resolution" in text
    assert "| agentA |" in text


def test_a_cohort_with_dilutions_and_no_calls_is_read_on_the_dilutions(tmp_path):
    raw = planted_cohort(tmp_path, mic=True, intake=False)
    del raw["dataset"]["phenotype"]
    _, cfg = write_config(tmp_path, raw)
    out = run(cfg, seed=7)
    assert out["n_isolates"] == 96 and out["n_traits"] == 0
    md = out["metadata_diagnostics"]
    assert "clonal_share" not in md and "lineage_evidence" not in md
    assert set(md["censored_share"]["per_agent"]) == {"agentA", "agentB", "agentC"}


def test_the_record_is_strict_json():
    import numpy as np
    import pandas as pd
    import pytest as _pytest
    from pathlib import Path

    from amr_clonalshare.jsonio import dumps, to_jsonable

    obj = {"a": np.float64("nan"), "b": np.array([1, 2]), "c": np.bool_(True),
           "d": pd.Series({"x": 1.5}), "e": pd.DataFrame({"y": [1]}),
           "f": Path("p"), "g": (1, 2), "h": pd.Timestamp("2026-09-09")}
    text = dumps(obj)
    assert "NaN" not in text and '"a": null' in text
    assert to_jsonable(obj)["b"] == [1, 2] and to_jsonable(obj)["f"] == "p"
    with _pytest.raises(TypeError):
        dumps({"z": object()})
