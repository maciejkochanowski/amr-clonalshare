"""A binary layer cannot hold "not measured", so the loader has to.

The failure this guards: an empty cell, or an isolate absent from one layer
under the union policy, was cast to integer and either crashed with a pandas
message or, once a user filled the gap with 0, was read as absence of the
trait. Every branch of the policy is exercised on a real config and a real
CSV, and the record the run writes is checked for the numbers a reader needs.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest
import yaml

from amr_clonalshare import qc
from amr_clonalshare.attribution import SUPPORT_THRESHOLD
from amr_clonalshare.cli import main
from amr_clonalshare.config import ConfigError, from_dict
from amr_clonalshare.io import load_dataset


def _cohort(tmp_path, *, hole=None, union_extra=0, value=None, policy=None,
            alignment=None, lineage=True):
    rng = np.random.default_rng(11)
    n = 40
    idx = [f"s{i}" for i in range(n)]
    amr = pd.DataFrame(rng.integers(0, 2, (n, 6)), index=idx,
                       columns=[f"a{j}" for j in range(6)]).astype(float)
    vir = pd.DataFrame(rng.integers(0, 2, (n, 5)), index=idx,
                       columns=[f"v{j}" for j in range(5)]).astype(float)
    if hole is not None:
        amr.loc[hole[0], hole[1]] = np.nan
    if value is not None:
        amr.loc["s1", "a1"] = value
    if union_extra:
        extra = pd.DataFrame(rng.integers(0, 2, (union_extra, 6)),
                             index=[f"x{i}" for i in range(union_extra)],
                             columns=amr.columns).astype(float)
        amr = pd.concat([amr, extra])
    d = tmp_path / "data"
    d.mkdir(exist_ok=True)
    amr.rename_axis("isolate").to_csv(d / "amr.csv")
    vir.rename_axis("isolate").to_csv(d / "vir.csv")
    raw = {"dataset": {"name": "t", "strain_id_column": "isolate", "data_dir": "data"},
           "files": {"amr": {"path": "amr.csv", "kind": "wide_binary"},
                     "vir": {"path": "vir.csv", "kind": "wide_binary"}},
           "trait_cluster": {"layers": ["amr", "vir"], "k_range": [2, 3],
                             "prevalence_gate": {"lo": 0.0, "hi": 1.0, "min_count": 2}},
           "snf": {"K": 10, "T": 10}, "tva": {"n_splits": 7},
           "influence": {"n_perm": 0}}
    if lineage:
        meta = pd.DataFrame({"lineage": ["L0"] * 30 + [f"L{i}" for i in range(1, 11)]},
                            index=idx)
        meta.rename_axis("isolate").to_csv(d / "metadata.csv")
        raw["dataset"].update({"metadata": "metadata.csv",
                               "lineage_column": "lineage"})
    if policy:
        raw["dataset"]["missing_policy"] = policy
    if alignment:
        raw["dataset"]["strain_alignment_policy"] = alignment
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(yaml.safe_dump(raw, sort_keys=False))
    return cfg_path, from_dict(raw, config_path=cfg_path).validate()


# ------------------------------------------------------------- the policy

def test_a_complete_cohort_loads_and_reports_no_missing(tmp_path):
    _, cfg = _cohort(tmp_path)
    ds = load_dataset(cfg)
    rec = ds.input_qc["missing"]
    assert rec["policy"] == "refuse"
    assert all(v["cells"] == 0 for v in rec["layers"].values())
    assert ds.input_qc["n_isolates_aligned"] == 40


def test_an_empty_cell_is_refused_with_its_position_named(tmp_path):
    _, cfg = _cohort(tmp_path, hole=("s3", "a2"))
    with pytest.raises(ConfigError) as caught:
        load_dataset(cfg)
    msg = str(caught.value)
    assert "1 empty cells" in msg and "'s3'" in msg and "'a2'" in msg
    assert "absence of the trait" in msg


def test_drop_rows_removes_the_isolate_from_every_layer(tmp_path):
    _, cfg = _cohort(tmp_path, hole=("s3", "a2"), policy="drop_rows")
    ds = load_dataset(cfg)
    assert "s3" not in ds.strain_ids
    assert all("s3" not in ds.binary(r).index for r in ("amr", "vir"))
    assert ds.input_qc["missing"]["dropped_rows"] == ["s3"]
    assert ds.n == 39
    assert ds.binary("amr").dtypes.eq(int).all()


def test_drop_columns_removes_the_feature_and_keeps_the_isolate(tmp_path):
    _, cfg = _cohort(tmp_path, hole=("s3", "a2"), policy="drop_columns")
    ds = load_dataset(cfg)
    assert "s3" in ds.strain_ids
    assert "a2" not in ds.binary("amr").columns
    assert ds.binary("vir").shape == (40, 5)
    assert ds.input_qc["missing"]["dropped_columns"] == {"amr": ["a2"], "vir": []}


def test_union_alignment_creates_missing_rows_that_the_policy_governs(tmp_path):
    _, cfg = _cohort(tmp_path, union_extra=3, alignment="union")
    with pytest.raises(ConfigError, match="15 empty cells in 3 rows"):
        load_dataset(cfg)
    _, cfg = _cohort(tmp_path, union_extra=3, alignment="union", policy="drop_rows")
    ds = load_dataset(cfg)
    assert ds.n == 40 and not any(str(i).startswith("x") for i in ds.strain_ids)


def test_a_value_other_than_zero_or_one_is_refused(tmp_path):
    _, cfg = _cohort(tmp_path, value=2.0)
    with pytest.raises(ConfigError, match="values other than 0 and 1"):
        load_dataset(cfg)


def test_an_unknown_missing_policy_is_a_config_error(tmp_path):
    with pytest.raises(ConfigError, match="missing_policy"):
        _cohort(tmp_path, policy="impute")


def test_missing_values_are_never_imputed():
    frames = {"a": pd.DataFrame({"x": [1.0, np.nan]}, index=["p", "q"])}
    for policy in ("drop_rows", "drop_columns"):
        out, rec = qc.apply_missing_policy(frames, policy)
        assert not out["a"].isna().any().any()
        assert rec.any_missing
    with pytest.raises(ConfigError):
        qc.apply_missing_policy(frames, "nonsense")


# ---------------------------------------------------------- the record

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
    # attribution._codes folds None, NaN, "", "nan", "NA" and "none" into one
    # __missing__ level that counts as a lineage. The check must report the
    # support the estimator will use, not a different one.
    from amr_clonalshare.attribution import _codes
    labels = ["A", "A", "B", "B", "B", None, None, None, "nan", "NA", "C"]
    g = qc.group_adequacy(pd.Series(labels, dtype=object))
    codes = _codes(np.asarray(labels, dtype=object))
    sizes = np.bincount(codes)
    assert g["support"] == pytest.approx(float((sizes[codes] >= 2).mean()))
    assert g["support"] == pytest.approx(10 / 11)
    assert g["n_groups"] == 4 and g["n_untyped"] == 5
    assert g["group_sizes"]["__missing__"] == 5 and g["n_singletons"] == 1
    lone = qc.group_adequacy(pd.Series(["A", "A", None, "B", "B", "B"]))
    assert lone["n_untyped"] == 1 and lone["n_groups"] == 3
    assert lone["support"] == pytest.approx(5 / 6)
    assert qc.group_adequacy(pd.Series([None, None]))["n_groups"] == 1
    assert qc.group_adequacy(pd.Series([], dtype=object))["n_groups"] == 0


def test_true_and_false_cells_are_refused(tmp_path):
    _, cfg = _cohort(tmp_path)
    d = tmp_path / "data" / "amr.csv"
    df = pd.read_csv(d, index_col="isolate")
    df["a0"] = df["a0"].astype(bool)
    df.to_csv(d)
    with pytest.raises(ConfigError, match="values other than 0 and 1"):
        load_dataset(cfg)


def test_an_emptied_cohort_and_missing_columns_are_refused_with_a_sentence(tmp_path):
    # every row has a hole under drop_rows: nothing remains
    _, cfg = _cohort(tmp_path, policy="drop_rows")
    d = tmp_path / "data" / "amr.csv"
    df = pd.read_csv(d, index_col="isolate")
    df["a0"] = np.nan
    df.to_csv(d)
    with pytest.raises(ConfigError, match="no isolates remain"):
        load_dataset(cfg)
    # the strain column is absent from a layer
    df = pd.read_csv(d, index_col="isolate").rename_axis("sample")
    df.to_csv(d)
    with pytest.raises(ConfigError, match="has no column 'isolate'"):
        load_dataset(cfg)
    # the lineage column is absent from the metadata
    _, cfg = _cohort(tmp_path)
    m = tmp_path / "data" / "metadata.csv"
    pd.read_csv(m).rename(columns={"lineage": "st"}).to_csv(m, index=False)
    with pytest.raises(ConfigError, match="has no column 'lineage'"):
        load_dataset(cfg)


def test_a_broken_yaml_config_is_a_config_error(tmp_path):
    from amr_clonalshare.config import load_config
    p = tmp_path / "bad.yaml"
    p.write_text("dataset: {name: x\nfiles: [")
    with pytest.raises(ConfigError, match="not valid YAML"):
        load_config(p)


def test_trait_adequacy_uses_the_rarer_outcome():
    X = pd.DataFrame({"rare": [1] * 5 + [0] * 95, "common": [1] * 60 + [0] * 40,
                      "flat": [0] * 100})
    t = qc.trait_adequacy(X)
    assert t["rare"]["minor_count"] == 5 and not t["rare"]["adequate"]
    assert t["common"]["minor_count"] == 40 and t["common"]["adequate"]
    assert t["flat"]["constant"] and not t["flat"]["adequate"]


def test_the_record_carries_the_lineage_verdict_and_the_metadata_join(tmp_path):
    _, cfg = _cohort(tmp_path)
    ds = load_dataset(cfg)
    rec = ds.input_qc
    assert rec["metadata_join"]["share_joined"] == 1.0
    assert rec["lineage"]["column"] == "lineage"
    assert rec["lineage"]["n_singletons"] == 10
    assert rec["lineage"]["estimable"] is False
    text = qc.render_markdown(rec)
    assert "not estimable" in text and "10 lineages hold a single isolate" in text


def test_render_markdown_covers_every_branch(tmp_path):
    _, cfg = _cohort(tmp_path, hole=("s3", "a2"), policy="drop_rows", lineage=False)
    ds = load_dataset(cfg)
    text = qc.render_markdown(ds.input_qc)
    assert "Isolates removed: 1." in text and "Nothing was filled in" in text
    _, cfg = _cohort(tmp_path, hole=("s3", "a2"), policy="drop_columns", lineage=False)
    text = qc.render_markdown(load_dataset(cfg).input_qc)
    assert "Features removed from layer `amr`: 1." in text
    rec = load_dataset(cfg).input_qc
    rec["layers"]["amr"]["traits"]["a0"]["constant"] = True
    assert "single value" in qc.render_markdown(rec)
    rec["metadata_join"] = {"n_joined": 1, "share_joined": 0.5,
                            "unjoined_examples": ["s9"]}
    rec["lineage"] = {"column": "l", **qc.group_adequacy(pd.Series(["A"] * 40))}
    text = qc.render_markdown(rec)
    assert "s9" in text and "can be estimated" in text


# ---------------------------------------------------------------- the CLI

def test_check_input_writes_the_record_and_stops(tmp_path, capsys):
    cfg_path, _ = _cohort(tmp_path)
    out = tmp_path / "res"
    assert main(["--config", str(cfg_path), "--results-dir", str(out),
                 "--check-input"]) == 0
    assert (out / "input_qc.json").is_file() and (out / "input_qc.md").is_file()
    rec = json.loads((out / "input_qc.json").read_text())
    assert rec["n_isolates_aligned"] == 40
    assert "# Input check" in capsys.readouterr().out
    assert not (out / "cluster_result.json").exists()


def test_a_refused_input_leaves_with_the_config_exit_code(tmp_path, capsys):
    cfg_path, _ = _cohort(tmp_path, hole=("s3", "a2"))
    assert main(["--config", str(cfg_path), "--check-input"]) == 2
    assert "input error:" in capsys.readouterr().err
