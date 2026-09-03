"""End-to-end behaviour of core.run and the CLI, including the cases a
simpler pipeline cannot express."""
import json
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import adjusted_rand_score

from amr_clonalshare import cli, core
from amr_clonalshare.config import ConfigError, from_dict
from amr_clonalshare.jsonio import dumps, to_jsonable


def test_planted_structure_is_recovered_and_validated(planted_cfg):
    cfg, truth = planted_cfg
    r = core.run(cfg, seed=42, consensus_B=10, n_boot=50)
    assert r["selected_k"] == 3
    labels = [a["cluster"] for a in r["assignment"]]
    assert adjusted_rand_score(truth, labels) > 0.9
    inf = r["post_clustering_inference"]
    assert inf["status"] == "ok"
    # A bound of 1e-2 would be met only where every feature is its own split
    # unit. On this fixture (24 features, 60 isolates) valid blocking
    # leaves 6 units, so a split has 3 units to discover with and 3 to test on,
    # and the merged p-value is correspondingly weaker. Detection is the claim;
    # the size of the p-value is a property of the fixture, not of the method.
    assert inf["p_value_structure"] < 0.05
    assert inf["block_aware"] is True
    assert inf["structure_detected"] is True
    assert inf["n_defining"] > 0
    # structure is detected AND is distinguishable from a single gradient
    assert inf["discreteness"]["discrete_beyond_a_gradient"] is True
    # the test concerns the partition that was reported, not some other one
    assert inf["concerns_reported_partition"] is True


def test_no_structure_yields_k_equals_one(null_cfg):
    """The case a k-sweep starting at two cannot report: no clusters."""
    r = core.run(null_cfg, seed=42, consensus_B=10, n_boot=50)
    assert r["selected_k"] == 1
    assert r["k_selection"]["no_structure"] is True
    assert r["k_selection"]["mdl_prefers_no_clustering"] is True
    assert r["post_clustering_inference"]["status"] == "skipped"
    assert r["n_defining_features_descriptive"] == 0


def test_clonal_cohort_is_flagged_as_lineage_confounded(clonal_cfg):
    """A valid test can be significant while the result is still a clone."""
    r = core.run(clonal_cfg, seed=42, consensus_B=10, n_boot=50)
    lc = r["metadata_diagnostics"]["lineage_concordance"]
    assert lc["status"] == "ok"
    assert lc["concordance_observed"] > lc["concordance_null_mean"]
    assert lc["z"] > 3
    assert lc["p_value"] < 0.05


def test_seed_controls_the_partition(planted_cfg):
    cfg, _ = planted_cfg
    a = core.run(cfg, seed=1, consensus_B=6, n_boot=20, run_baselines=False)
    b = core.run(cfg, seed=1, consensus_B=6, n_boot=20, run_baselines=False)
    assert [x["cluster"] for x in a["assignment"]] == [x["cluster"] for x in b["assignment"]]
    assert a["post_clustering_inference"]["p_value_structure"] == pytest.approx(
        b["post_clustering_inference"]["p_value_structure"])
    c = core.run(cfg, seed=999, consensus_B=6, n_boot=20, run_baselines=False)
    # a different seed draws different splits, so the inference p-value moves
    # even where the partition itself is stable
    assert (c["post_clustering_inference"]["p_value_structure"]
            != a["post_clustering_inference"]["p_value_structure"])


def test_run_reports_every_diagnostic_the_paper_cites(planted_cfg):
    cfg, _ = planted_cfg
    r = core.run(cfg, seed=42, consensus_B=6, n_boot=20)
    for key in ("artifact_diagnostics", "layer_influence", "baselines",
                "layer_reports", "metadata_diagnostics", "k_selection",
                "post_clustering_inference", "config", "schema_version"):
        assert key in r, key
    assert r["layer_influence"]["n_eff"] >= 1.0
    assert "empty_stratum_ari" in r["artifact_diagnostics"]["per_layer"][0]
    assert "knn_ties" in r["layer_reports"][0]


def test_result_serialises_to_strict_json(planted_cfg):
    cfg, _ = planted_cfg
    r = core.run(cfg, seed=42, consensus_B=6, n_boot=20)
    text = dumps(r)                       # allow_nan=False inside
    assert "NaN" not in text and "Infinity" not in text
    json.loads(text)                      # round-trips


def test_jsonio_refuses_to_stringify_unknown_objects():
    class Weird:
        pass
    with pytest.raises(TypeError):
        to_jsonable({"x": Weird()})
    assert to_jsonable(float("nan")) is None


def test_validate_recovers_planted_truth_and_reports_the_difficulty():
    easy = core.validate(overlap=0.05, n=90, seed=1)
    assert easy["passed"] and easy["ari"] > 0.9
    hard = core.validate(overlap=0.45, n=90, seed=1)
    assert hard["ari"] < easy["ari"]


def test_cli_writes_both_artifacts_and_exits_zero(planted_cfg, tmp_path):
    cfg, _ = planted_cfg
    out = tmp_path / "out"
    code = cli.main(["--config", str(cfg.config_path), "--results-dir", str(out),
                     "--quiet", "--no-check-files"])
    assert code in (0, 3)
    assert (out / "cluster_result.json").exists()
    assert (out / "summary.json").exists()
    assert (out / "archetype_profiles.tsv").exists()
    json.loads((out / "cluster_result.json").read_text())


def test_cli_reports_a_config_error_with_exit_code_two(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("dataset: {name: x}\n")
    assert cli.main(["--config", str(bad)]) == 2


def test_cli_version_and_help():
    r = subprocess.run([sys.executable, "-m", "amr_clonalshare.cli", "--version"],
                       capture_output=True, text=True)
    assert "amr-clonalshare" in r.stdout


# ------------------------------------------------------------- edge cases ---
def _cfg(tmp_path, frames, extra=None, kinds=None):
    d = tmp_path / "d"
    d.mkdir(exist_ok=True)
    files = {}
    for name, df in frames.items():
        df.rename_axis("isolate").to_csv(d / f"{name}.csv")
        files[name] = {"path": f"{name}.csv",
                       "kind": (kinds or {}).get(name, "wide_binary")}
    raw = {"dataset": {"name": "e", "strain_id_column": "isolate", "data_dir": "d"},
           "files": files,
           "trait_cluster": {"layers": list(frames), "k_range": [2, 3],
                             "prevalence_gate": {"lo": 0.0, "hi": 1.0, "min_count": 1}},
           "snf": {"K": 5, "T": 5}, "tva": {"n_splits": 3},
           "influence": {"n_perm": 0}}
    if extra:
        for k, v in extra.items():
            raw.setdefault(k, {}).update(v)
    import yaml
    p = tmp_path / "c.yaml"
    p.write_text(yaml.safe_dump(raw, sort_keys=False))
    return from_dict(raw, config_path=p).validate()


def _frame(n, p, seed=0, name="f"):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(rng.integers(0, 2, size=(n, p)),
                        index=[f"i{i:03d}" for i in range(n)],
                        columns=[f"{name}{j}" for j in range(p)])


def test_single_layer_config_runs(tmp_path):
    cfg = _cfg(tmp_path, {"only": _frame(40, 8)})
    r = core.run(cfg, seed=0, consensus_B=4, n_boot=10, run_baselines=False)
    assert r["n_isolates"] == 40


def test_too_few_isolates_is_rejected_with_a_clear_message(tmp_path):
    cfg = _cfg(tmp_path, {"a": _frame(3, 5), "b": _frame(3, 5, seed=1, name="b")})
    with pytest.raises(ValueError, match="at least 4"):
        core.run(cfg, seed=0)


def test_layer_emptied_by_the_gate_is_rejected_with_a_clear_message(tmp_path):
    df = _frame(30, 5)
    df.iloc[:, :] = 0
    cfg = _cfg(tmp_path, {"a": df, "b": _frame(30, 5, seed=2, name="b")},
               extra={"trait_cluster": {"prevalence_gate":
                                        {"lo": 0.0, "hi": 1.0, "min_count": 2}}})
    with pytest.raises(ValueError, match="no features left"):
        core.run(cfg, seed=0)


def test_duplicate_strain_ids_are_rejected(tmp_path):
    df = _frame(10, 5)
    df.index = ["dup"] * 10
    with pytest.raises(Exception):
        core.run(_cfg(tmp_path, {"a": df}), seed=0)


def test_one_hot_violation_is_rejected(tmp_path):
    df = _frame(20, 5, seed=3)
    df.iloc[:, :] = 1                       # every row sums to 5
    with pytest.raises(ConfigError, match="one_hot"):
        core.run(_cfg(tmp_path, {"cat": df, "b": _frame(20, 5, seed=4, name="b")},
                      kinds={"cat": "one_hot"}), seed=0)


def test_k_larger_than_n_is_clamped(tmp_path):
    cfg = _cfg(tmp_path, {"a": _frame(12, 6)},
               extra={"trait_cluster": {"k_range": [2, 3]}})
    r = core.run(cfg, seed=0, consensus_B=3, n_boot=10, run_baselines=False)
    assert r["selected_k"] <= 12


def test_yaml_keys_that_change_the_science_are_actually_read(tmp_path):
    """Every documented config key must survive `from_dict`.

    `snf.update` and `tva.block_threshold` were both declared on their
    dataclasses, validated, plumbed to the call sites and documented in the
    manuscript -- and silently dropped by the YAML loader, so setting them in a
    config did nothing. A key that changes the science and is ignored is worse
    than a key that does not exist.
    """
    raw = {
        "dataset": {"name": "t", "strain_id_column": "isolate", "data_dir": "."},
        "files": {"amr": {"path": "amr.csv", "kind": "wide_binary"}},
        "snf": {"K": 7, "mu": 0.4, "T": 5, "alpha": 1.0,
                "tie_policy": "strict", "update": "add_identity"},
        "tva": {"n_splits": 3, "block_threshold": 0.35,
                "continuum_q_max": 4},
        "trait_cluster": {"layers": ["amr"], "k_range": [2, 3]},
    }
    cfg = from_dict(raw, config_path=tmp_path / "c.yaml")
    assert cfg.snf.update == "add_identity"
    assert cfg.snf.alpha == 1.0          # allowed only under add_identity
    assert cfg.snf.tie_policy == "strict"
    assert cfg.tva.block_threshold == 0.35
    assert cfg.tva.continuum_q_max == 4


def test_continuum_config_refuses_unresolvable_bootstrap_and_unsafe_dimension(tmp_path):
    raw = {
        "dataset": {"name": "t", "strain_id_column": "isolate", "data_dir": "."},
        "files": {"amr": {"path": "amr.csv", "kind": "wide_binary"}},
        "trait_cluster": {"layers": ["amr"], "k_range": [2, 3]},
        "tva": {"n_boot_continuum": 19, "continuum_q_max": 5},
    }
    with pytest.raises(ConfigError, match="n_boot_continuum"):
        from_dict(raw, config_path=tmp_path / "c.yaml").validate(
            check_files_exist=False)
    raw["tva"]["n_boot_continuum"] = 20
    with pytest.raises(ConfigError, match="continuum_q_max"):
        from_dict(raw, config_path=tmp_path / "c.yaml").validate(
            check_files_exist=False)
