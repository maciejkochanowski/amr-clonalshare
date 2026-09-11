"""The input check: what the loader found, and what the estimator will do with it.

Every branch is exercised on a real config and a real CSV, and the record the
run writes is checked for the numbers a reader needs before spending compute.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from amr_clonalshare import qc
from amr_clonalshare.attribution import SUPPORT_THRESHOLD
from amr_clonalshare.cli import main
from amr_clonalshare.config import ConfigError, load_config
from amr_clonalshare.io import load_dataset
from conftest import planted_cohort, write_config


def test_the_loader_reads_the_panel_and_joins_the_metadata(share_cfg):
    _, cfg = share_cfg
    ds = load_dataset(cfg)
    assert ds.n == 96 and list(ds.panel.columns) == ["agentA", "agentB", "agentC"]
    rec = ds.input_qc
    assert rec["n_isolates"] == 96 and rec["n_antimicrobials"] == 3
    assert rec["metadata_join"]["share_joined"] == 1.0
    assert rec["lineage"]["column"] == "lineage"
    assert rec["lineage"]["n_groups"] == 8 and rec["lineage"]["n_singletons"] == 0
    assert rec["lineage"]["estimable"] is True
    assert all(t["adequate"] for t in rec["traits"].values())


def test_group_adequacy_matches_the_estimator_definition_of_support():
    lin = pd.Series(["L0"] * 30 + [f"L{i}" for i in range(1, 11)])
    g = qc.group_adequacy(lin)
    assert g["n_groups"] == 11 and g["n_singletons"] == 10
    assert g["support"] == pytest.approx(30 / 40)
    assert g["estimable"] is (30 / 40 >= SUPPORT_THRESHOLD)
    n0 = (40 - (30 ** 2 + 10) / 40) / 10
    assert g["effective_group_size"] == pytest.approx(n0)
    assert g["smallest_groups"]["L1"] == 1


def test_group_adequacy_treats_untyped_isolates_exactly_as_the_estimator_does():
    # The estimator sets untyped isolates (None, NaN, "", "nan", "NA", "none")
    # aside before it fits, so the check reads support on the typed isolates
    # alone and reports the untyped count beside it.
    from amr_clonalshare.attribution import clonal_share
    labels = ["A", "A", "B", "B", "B", None, None, None, "nan", "NA", "C"]
    g = qc.group_adequacy(pd.Series(labels, dtype=object))
    y = np.array([0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1], dtype=float)
    r = clonal_share(y, np.asarray(labels, dtype=object), n_boot=5, n_perm=5,
                     repeats=2, seed=1)
    assert g["support"] == pytest.approx(r.support)
    assert g["support"] == pytest.approx(5 / 6)
    assert g["n"] == 11 and g["n_typed"] == 6 and g["n_untyped"] == 5
    assert g["n_groups"] == 3 and g["n_singletons"] == 1
    assert "__missing__" not in g["group_sizes"]
    lone = qc.group_adequacy(pd.Series(["A", "A", None, "B", "B", "B"]))
    assert lone["n_untyped"] == 1 and lone["n_groups"] == 2
    assert lone["support"] == pytest.approx(1.0)
    nothing = qc.group_adequacy(pd.Series([None, None]))
    assert nothing["n_groups"] == 0 and nothing["estimable"] is False
    assert qc.group_adequacy(pd.Series([], dtype=object))["n_groups"] == 0


def test_a_broken_yaml_config_is_a_config_error(tmp_path):
    p = tmp_path / "bad.yaml"
    p.write_text("dataset: {name: x\nfiles: [")
    with pytest.raises(ConfigError, match="not valid YAML"):
        load_config(p)


def test_trait_adequacy_uses_the_rarer_outcome():
    X = pd.DataFrame({"rare": [1] * 5 + [0] * 95, "common": [1] * 60 + [0] * 40,
                      "flat": [0] * 100})
    # an untested cell is neither outcome: the counts are over tested isolates
    sparse = pd.DataFrame({"sparse": [1.0] * 26 + [0.0] * 4 + [np.nan] * 30})
    t = qc.trait_adequacy(sparse)["sparse"]
    assert t["n"] == 30 and t["minor_count"] == 4 and not t["adequate"]
    assert t["prevalence"] == pytest.approx(26 / 30)
    t = qc.trait_adequacy(X)
    assert t["rare"]["minor_count"] == 5 and not t["rare"]["adequate"]
    assert t["common"]["minor_count"] == 40 and t["common"]["adequate"]
    assert t["flat"]["constant"] and not t["flat"]["adequate"]


def test_a_fine_lineage_definition_is_reported_as_not_estimable(tmp_path):
    raw = planted_cohort(tmp_path, n_per_lineage=1, n_lineages=40)
    _, cfg = write_config(tmp_path, raw)
    rec = load_dataset(cfg).input_qc
    assert rec["lineage"]["n_singletons"] == 40
    assert rec["lineage"]["estimable"] is False
    text = qc.render_markdown(rec)
    assert "not estimable" in text and "40 lineages hold a single isolate" in text


def test_render_markdown_covers_every_branch(share_cfg):
    _, cfg = share_cfg
    rec = load_dataset(cfg).input_qc
    rec["traits"]["agentA"]["constant"] = True
    rec["traits"]["agentB"]["adequate"] = False
    rec["metadata_join"] = {"n_joined": 1, "share_joined": 0.5,
                            "unjoined_examples": ["s9"]}
    rec["mic_join"] = {"rows_joined": 10, "strains_with_mic": 5,
                       "strains_aligned": 96, "antimicrobials": ["a", "b"]}
    text = qc.render_markdown(rec)
    assert "single value" in text and "Below the threshold: agentB" in text
    assert "s9" in text and "can be estimated" in text
    assert "Recorded dilutions: 10 rows for 5 of 96 isolates" in text


def test_a_phenotype_table_without_the_named_columns_is_refused(tmp_path):
    raw = planted_cohort(tmp_path)
    raw["dataset"]["phenotype_call_column"] = "verdict"
    _, cfg = write_config(tmp_path, raw)
    with pytest.raises(ConfigError, match="verdict"):
        load_dataset(cfg)


def test_a_metadata_table_without_the_lineage_column_is_refused(tmp_path):
    raw = planted_cohort(tmp_path)
    raw["dataset"]["lineage_column"] = "clade"
    _, cfg = write_config(tmp_path, raw)
    with pytest.raises(ConfigError, match="clade"):
        load_dataset(cfg)


# ---------------------------------------------------------------- the CLI

def test_check_input_writes_the_record_and_stops(share_cfg, tmp_path, capsys):
    cfg_path, _ = share_cfg
    out = tmp_path / "res"
    assert main(["--config", str(cfg_path), "--results-dir", str(out),
                 "--check-input"]) == 0
    assert (out / "input_qc.json").is_file() and (out / "input_qc.md").is_file()
    rec = json.loads((out / "input_qc.json").read_text())
    assert rec["n_isolates"] == 96
    assert "# Input check" in capsys.readouterr().out
    assert not (out / "clonal_share_result.json").exists()


def test_a_refused_input_leaves_with_the_config_exit_code(tmp_path, capsys):
    raw = planted_cohort(tmp_path)
    raw["dataset"]["lineage_column"] = "clade"
    cfg_path, _ = write_config(tmp_path, raw)
    assert main(["--config", str(cfg_path), "--check-input"]) == 2
    assert "input error:" in capsys.readouterr().err


def test_a_missing_file_leaves_with_the_config_exit_code(tmp_path, capsys):
    raw = planted_cohort(tmp_path)
    raw["dataset"]["phenotype"] = "absent.csv"
    cfg_path = tmp_path / "config.yaml"
    import yaml
    cfg_path.write_text(yaml.safe_dump(raw))
    assert main(["--config", str(cfg_path)]) == 2
    assert "config error:" in capsys.readouterr().err
