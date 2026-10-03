"""A cohort is read for the share its lineages carry, and nothing else.

The run reaches every estimator through one function and writes one record
whose shape is fixed by the schema.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

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
                        "metadata_diagnostics", "input_qc", "collection_bounds", "results"}
    assert out["schema_version"] == "1.0"
    assert len(out["collection_bounds"]) == 3
    assert {row["agent"] for row in out["results"]} == set(out["traits"])
    md = out["metadata_diagnostics"]
    assert set(md) <= {"lineage_column", "clonal_share", "realised_share",
                       "lineage_selection", "lineage_evidence", "censored_share",
                       "prevalence_decomposition", "sampling_units", "call_strata"}


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
    """The documented command, on a configuration that names a panel and a
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
    assert "sha256:" in (out / "report.md").read_text(encoding="utf-8")


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
    assert entry["n"] == 96 and "panel" in entry and "order" in entry and "share" not in entry
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


def test_sampling_units_and_concentration_units_are_kept_apart(tmp_path):
    """A farm column confines the permutations of the test of every reading,
    calls and dilutions alike, while the units of the dilutions stay a label
    of the scale; neither is read as the other."""
    raw = planted_cohort(tmp_path, mic=True)
    meta = pd.read_csv(tmp_path / "data" / "meta.csv")
    meta["farm"] = [f"F{i % 4}" for i in range(len(meta))]
    meta.to_csv(tmp_path / "data" / "meta.csv", index=False)
    raw["dataset"].update(unit_column="farm", mic_units={"agentA": "mg/L"})
    _, cfg = write_config(tmp_path, raw)
    out = run(cfg, seed=7)
    assert out["metadata_diagnostics"]["sampling_units"]["n_units"] == 4
    entry = out["metadata_diagnostics"]["censored_share"]["per_agent"]["agentA"]
    assert entry["mic_units"] == "mg/L"
    assert entry["order"]["estimable"]
    for agent in ("agentA", "agentB", "agentC"):
        order = out["metadata_diagnostics"]["censored_share"]["per_agent"][agent]["order"]
        assert order["n_strata"] == 4 and order["sampling_units"]
        assert out["metadata_diagnostics"]["clonal_share"][agent]["n_strata"] == 4


def test_units_reach_the_decomposition_and_an_isolate_without_one_is_its_own(tmp_path):
    """The farm column is read once: isolates without a recorded farm are
    units of their own and are counted, and the decomposition of a
    prevalence difference draws whole farms."""
    raw = planted_cohort(tmp_path)
    meta = pd.read_csv(tmp_path / "data" / "meta.csv")
    meta["farm"] = [f"F{i % 6}" for i in range(len(meta))]
    meta.loc[meta.index[:3], "farm"] = None
    meta["period"] = ["early" if i % 2 else "late" for i in range(len(meta))]
    meta.to_csv(tmp_path / "data" / "meta.csv", index=False)
    raw["dataset"].update(unit_column="farm", contrast_column="period",
                          contrast_levels=["early", "late"])
    raw["surveillance"] = {"n_boot": 50}
    _, cfg = write_config(tmp_path, raw)
    out = run(cfg, seed=7)
    units = out["metadata_diagnostics"]["sampling_units"]
    assert units == {"column": "farm", "n_units": 9, "n_isolates_without_unit": 3}
    decomposition = out["metadata_diagnostics"]["prevalence_decomposition"]["per_feature"]["agentA"]
    assert all(np.isfinite(decomposition["within_lineage_ci95"]))
    raw["dataset"]["unit_column"] = "no_such_column"
    with pytest.raises((ConfigError, ValueError)):
        run(write_config(tmp_path, raw)[1], seed=7)


def test_a_unit_column_sets_aside_readings_in_another_unit(tmp_path):
    """An export that keeps zone diameters beside the MICs: the rows in mm
    are set aside and counted per agent, a blank unit is read and counted,
    the unit read is recorded for every agent, and an agent tested in
    another unit only leaves no MIC analysis."""
    raw = planted_cohort(tmp_path, mic=True)
    mic = pd.read_csv(tmp_path / "data" / "mic.csv")
    mic["unit"] = "mg/L"
    mic.loc[mic.index[:5], "unit"] = "ug/mL"
    mic.loc[mic.index[5:7], "unit"] = None
    zones = mic.iloc[:4].copy()
    zones["antibiotic"] = "agentD"
    zones["measurement"] = [22, 18, 30, 25]
    zones["unit"] = "mm"
    # one isolate's agentA reading is a zone diameter, not an MIC
    swapped = mic.index[mic["antibiotic"] == "agentA"][3]
    mic.loc[swapped, ["measurement", "unit"]] = [19, "mm"]
    pd.concat([mic, zones]).to_csv(tmp_path / "data" / "mic.csv", index=False)
    raw["dataset"]["mic_unit_column"] = "unit"
    _, cfg = write_config(tmp_path, raw)
    with pytest.warns(RuntimeWarning, match="unit other than mg/L"):
        out = run(cfg, seed=7)
    join = out["input_qc"]["mic_join"]
    assert join["units_excluded"] == {"mm": 5}
    assert join["agents_in_other_units_only"] == ["agentD"]
    assert join["mic_unit_column"] == "unit"
    assert join["mic_units"] == {"agentA": "mg/L", "agentB": "mg/L", "agentC": "mg/L"}
    assert join["agents"]["agentA"]["n_unit_excluded"] == 1
    assert sum(r["n_unit_blank"] for r in join["agents"].values()) == 2
    per_agent = out["metadata_diagnostics"]["censored_share"]["per_agent"]
    assert set(per_agent) == {"agentA", "agentB", "agentC"}
    assert per_agent["agentA"]["mic_units"] == "mg/L"
    assert per_agent["agentA"]["n"] == 95
    assert any(r["analysis"] == "mic_order" and "MIC units: mg/L" in r["assumptions"] for r in out["results"])
    raw["dataset"].update(mic_unit_column="unit", mic_units={"agentB": "µg/mL"})
    out = run(write_config(tmp_path, raw)[1], seed=7)
    assert out["input_qc"]["mic_join"]["mic_units"]["agentB"] == "µg/mL"
    raw["dataset"]["mic_unit_column"] = "no_such_column"
    with pytest.raises(ConfigError):
        run(write_config(tmp_path, raw)[1], seed=7)


def test_a_phenotype_covariate_column_reads_the_calls_within_its_levels(tmp_path):
    """Two laboratories whose calls differ in level: with the column declared
    every call is centred on the prevalence of its laboratory, so the share
    is that of the calls within laboratories and differs from the pooled one
    where the levels differ; the permutations, the interval and the bounds
    are read within laboratories, the record names the strata, and an
    isolate with a call but no laboratory is refused."""
    raw = planted_cohort(tmp_path)
    meta = pd.read_csv(tmp_path / "data" / "meta.csv")
    meta["laboratory"] = ["north" if i % 2 else "south" for i in range(len(meta))]
    meta.to_csv(tmp_path / "data" / "meta.csv", index=False)
    pooled = run(write_config(tmp_path, raw)[1], seed=7)
    raw["dataset"]["phenotype_covariate_column"] = "laboratory"
    _, cfg = write_config(tmp_path, raw)
    out = run(cfg, seed=7)
    assert out["metadata_diagnostics"]["call_strata"] == {"column": "laboratory", "n_strata": 2}
    for agent in ("agentA", "agentB", "agentC"):
        share = out["metadata_diagnostics"]["clonal_share"][agent]
        assert share["n_strata"] == 2
        assert share["kappa"] != pooled["metadata_diagnostics"]["clonal_share"][agent]["kappa"]
        assert share["prevalence"] == pooled["metadata_diagnostics"]["clonal_share"][agent]["prevalence"]
        assert set(share["latent_order_lower_within_strata"]) == {"north", "south"}
    assert any(r["analysis"] == "collection_membership"
               and "within each level of laboratory" in r["assumptions"] for r in out["results"])
    from amr_clonalshare.cli import _summary
    from amr_clonalshare.report import render_report
    assert "Strata of the calls" in render_report(out, _summary(out), input_qc=out["input_qc"])
    meta.loc[meta.index[0], "laboratory"] = None
    meta.to_csv(tmp_path / "data" / "meta.csv", index=False)
    with pytest.raises(ValueError, match="phenotype covariate column"):
        run(write_config(tmp_path, raw)[1], seed=7)
    raw["dataset"]["phenotype_covariate_column"] = "no_such_column"
    with pytest.raises(ConfigError):
        run(write_config(tmp_path, raw)[1], seed=7)
