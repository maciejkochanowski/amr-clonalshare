"""The pipeline reports the validated MIC interval beside the approximate one."""
import numpy as np
import pandas as pd
import pytest

from amr_clonalshare import core
from amr_clonalshare.config import CensoredConfig, ConfigError, from_dict
from amr_clonalshare.outputs import result_rows

from conftest import planted_cohort, write_config


def small_panel():
    from amr_clonalshare._mic_panels import observe_panel
    rng = np.random.default_rng(8)
    labels = np.repeat(np.array([f"L{g}" for g in range(6)], dtype=object), 5)
    latent = rng.normal(size=6)[np.arange(6).repeat(5)] + rng.normal(size=30)
    lo, hi = observe_panel(latent, [-2., -1., 0., 1., 2.])
    return lo, hi, labels


@pytest.mark.real_calibration
def test_helper_returns_a_calibrated_interval_for_interval_readings():
    lo, hi, labels = small_panel()
    cen = CensoredConfig(calibration_n_boot=19, workers=2)
    cal = core._calibrated_mic_interval(lo, hi, labels, cen, 5)
    assert cal["method"] == "nullwise_parametric_bootstrap_LR_inversion"
    assert cal["n_boot"] == 19 and cal["seed"] == 5
    assert cal["panel_edges"] == [-2.0, -1.0, 0.0, 1.0, 2.0]
    if cal["estimable"]:
        assert 0 <= cal["low"] <= cal["estimate"] <= cal["high"] <= 1


@pytest.mark.real_calibration
def test_helper_uses_the_exact_pivot_for_exact_readings():
    rng = np.random.default_rng(3)
    labels = np.repeat(np.arange(8), 4).astype(str).astype(object)
    y = rng.normal(size=8)[np.arange(8).repeat(4)] + rng.normal(size=32)
    cal = core._calibrated_mic_interval(y, y.copy(), labels, CensoredConfig(), 5)
    assert cal["method"] == "exact_generalized_F"
    assert cal["estimable"] and 0 <= cal["low"] <= cal["high"] <= 1


def test_pipeline_carries_the_calibrated_result_to_tables_and_report(tmp_path, monkeypatch):
    fake = {"estimate": .4, "low": .2, "high": .6, "confidence": .95, "estimable": True,
            "method": "nullwise_parametric_bootstrap_LR_inversion", "status": "ok"}
    monkeypatch.setattr(core, "_calibrated_mic_interval", lambda *a, **k: dict(fake))
    raw = planted_cohort(tmp_path, mic=True)
    _, cfg = write_config(tmp_path, raw)
    out = core.run(cfg, results_dir=tmp_path / "out", seed=5)
    agents = out["metadata_diagnostics"]["censored_share"]["per_agent"]
    report = "".join(p.read_text(encoding="utf-8") for p in (tmp_path / "out").rglob("*.md"))
    assert "Calibrated share and 95 % interval" in report
    # the header alone proved nothing: the calibrated values must be printed,
    # and the figure must not present the approximate F interval as the
    # simulation-validated one.
    assert "0.400 (0.200 to 0.600)" in report
    assert "Every interval drawn is the calibrated interval" in report
    assert "interval for the Gaussian population model (thin)" not in report
    # the series itself, not only the caption: the drawn interval must be the
    # calibrated one, never the approximate F interval
    from amr_clonalshare.cli import _summary
    from amr_clonalshare.report_model import build_report
    drawn = [block for section in build_report(out, _summary(out)).sections
             for block in section.blocks
             if block.kind == "figure" and block.figure == "dilution"]
    assert drawn, "the dilution figure is missing from the report"
    assert all(row["calibrated"] and row["pt"] == .4
               and row["lo"] == .2 and row["hi"] == .6 for row in drawn[0].data)
    assert all(entry["calibrated"] == fake for entry in agents.values())
    rows = [r for r in result_rows(out) if r["analysis"] == "censored_mic_calibrated"]
    assert len(rows) == len(agents)
    assert all(r["status"] == "computed" and r["lower"] == .2 and r["upper"] == .6 for r in rows)


def test_calibrated_interval_can_be_switched_off(tmp_path):
    raw = planted_cohort(tmp_path, mic=True)
    raw["censored"] = {"calibrated_interval": False}
    _, cfg = write_config(tmp_path, raw)
    out = core.run(cfg, results_dir=None, seed=5)
    for entry in out["metadata_diagnostics"]["censored_share"]["per_agent"].values():
        assert "calibrated" not in entry


@pytest.mark.parametrize("option", [{"calibration_n_boot": 18}, {"workers": -1}])
def test_invalid_calibration_settings_are_rejected(tmp_path, option):
    raw = planted_cohort(tmp_path, mic=True)
    raw["censored"] = option
    with pytest.raises(ConfigError):
        from_dict(raw, config_path=tmp_path / "c.yaml").validate()


def test_numerical_failure_of_one_agent_is_recorded_and_the_run_completes(tmp_path, monkeypatch):
    from amr_clonalshare import censored

    def fail(*args, **kwargs):
        raise ArithmeticError("MIC quadrature failed its integration error check")
    monkeypatch.setattr(censored, "censored_clonal_share", fail)
    raw = planted_cohort(tmp_path, mic=True)
    raw["censored"] = {"calibrated_interval": False}
    _, cfg = write_config(tmp_path, raw)
    out = core.run(cfg, results_dir=tmp_path / "out", seed=5)
    agents = out["metadata_diagnostics"]["censored_share"]["per_agent"]
    assert agents
    for entry in agents.values():
        assert entry["share"]["status"] == "not_computed"
        assert "numerical failure" in entry["share"]["reason"]
    assert (tmp_path / "out" / "report.md").exists()


def test_published_results_directory_follows_the_umask(tmp_path):
    import os
    raw = planted_cohort(tmp_path)
    _, cfg = write_config(tmp_path, raw)
    mask = os.umask(0o022)
    try:
        core.run(cfg, results_dir=tmp_path / "out", seed=5)
    finally:
        os.umask(mask)
    assert (tmp_path / "out").stat().st_mode & 0o777 == 0o755
    assert (tmp_path / "out" / "report.md").stat().st_mode & 0o777 == 0o644


def test_panel_cut_points_follow_recorded_wells():
    from amr_clonalshare.censored import intervals_from_mic
    wells = [0.25, 0.5, 1, 2, 4, 8]
    lo, hi = intervals_from_mic([0.5, 0.5, 8, 0.25, 4], wells=wells)
    cuts = core._panel_cut_points(lo, hi, wells, True)
    assert cuts.tolist() == [-2.0, -1.0, 0.0, 1.0, 2.0]
    lo, hi = intervals_from_mic([0.5, 0.5, 8, 0.25, 4], wells=wells, treat_end_wells_as_censored=False)
    cuts = core._panel_cut_points(lo, hi, wells, False)
    assert cuts.tolist() == [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
    lo, hi = intervals_from_mic([0.5, 0.5, 8, 0.25, 4])
    assert core._panel_cut_points(lo, hi, None, True).tolist() == [-2.0, -1.0, 0.0, 1.0, 2.0]


def test_panel_column_reads_each_laboratory_on_its_own_wells(tmp_path, monkeypatch):
    seen = {}

    def record(lo, hi, lineage, cen, seed, wells=None, panel=None, covariate=None, n_boot=None):
        seen['panel'] = None if panel is None else sorted(set(panel))
        return {"estimable": False, "status": "not_computed", "reason": "stub"}
    monkeypatch.setattr(core, "_calibrated_mic_interval", record)
    raw = planted_cohort(tmp_path, mic=True)
    mic = pd.read_csv(tmp_path / "data" / "mic.csv")
    # laboratory B never records the lowest well, so its lowest reading is
    # an end well of B and an interior well of the pooled lattice
    mic["laboratory"] = np.where(mic["iid"].str.endswith(("0", "5")), "B", "A")
    mic.loc[mic["laboratory"].eq("B") & mic["measurement"].eq(0.25), "measurement"] = 0.5
    mic.to_csv(tmp_path / "data" / "mic.csv", index=False)
    raw["dataset"]["mic_panel_column"] = "laboratory"
    _, cfg = write_config(tmp_path, raw)
    out = core.run(cfg, results_dir=tmp_path / "out", seed=5)
    agents = out["metadata_diagnostics"]["censored_share"]["per_agent"]
    entry = agents["agentA"]
    assert seen['panel'] == ["A", "B"]
    assert entry["panel"]["panel_column"] == "laboratory"
    assert set(entry["panel"]["per_panel"]) == {"A", "B"}
    assert entry["panel"]["per_panel"]["B"]["lowest"] == 0.5
    assert entry["panel"]["per_panel"]["A"]["lowest"] == 0.25
    raw["dataset"]["mic_panel_column"] = "no_such_column"
    with pytest.raises(ConfigError):
        core.run(write_config(tmp_path, raw)[1], results_dir=tmp_path / "out2", seed=5)


def test_covariate_column_enters_the_calibrated_interval(tmp_path, monkeypatch):
    seen = {}

    def record(lo, hi, lineage, cen, seed, wells=None, panel=None, covariate=None, n_boot=None):
        seen['covariate'] = None if covariate is None else [sorted(set(c)) for c in covariate]
        return {"estimable": False, "status": "not_computed", "reason": "stub"}
    monkeypatch.setattr(core, "_calibrated_mic_interval", record)
    raw = planted_cohort(tmp_path, mic=True)
    meta = pd.read_csv(tmp_path / "data" / "meta.csv")
    meta["laboratory"] = np.where(meta["iid"].str.endswith(("0", "5")), "B", "A")
    meta.to_csv(tmp_path / "data" / "meta.csv", index=False)
    raw["dataset"]["mic_covariate_column"] = "laboratory"
    _, cfg = write_config(tmp_path, raw)
    out = core.run(cfg, results_dir=tmp_path / "out", seed=5)
    entry = out["metadata_diagnostics"]["censored_share"]["per_agent"]["agentA"]
    assert seen['covariate'] == [["A", "B"]]
    assert entry["calibrated"]["covariate_columns"] == ["laboratory"]
    raw["dataset"]["mic_covariate_column"] = "no_such_column"
    with pytest.raises(ConfigError):
        core.run(write_config(tmp_path, raw)[1], results_dir=tmp_path / "out2", seed=5)


@pytest.mark.real_calibration
def test_helper_reports_the_fixed_effects_of_a_covariate():
    lo, hi, labels = small_panel()
    laboratory = np.where(np.arange(lo.size) % 2 == 0, "lab1", "lab2")
    cen = CensoredConfig(calibration_n_boot=19, workers=2)
    cal = core._calibrated_mic_interval(lo, hi, labels, cen, 5, covariate=[laboratory])
    assert cal["method"] == "nullwise_parametric_bootstrap_LR_inversion"
    assert cal["adjustment"][0]["reference_level"] == "lab1" and set(cal["adjustment"][0]["coefficients"]) == {"lab2"}
    y = np.random.default_rng(3).normal(size=6)[np.arange(6).repeat(5)] + np.random.default_rng(4).normal(size=30)
    exact = core._calibrated_mic_interval(y, y.copy(), labels, cen, 5, covariate=[laboratory])
    assert exact["method"] == "nullwise_parametric_bootstrap_LR_inversion"


def test_two_covariate_columns_enter_the_model_together(tmp_path, monkeypatch):
    seen = {}

    def record(lo, hi, lineage, cen, seed, wells=None, panel=None, covariate=None, n_boot=None):
        seen['covariate'] = None if covariate is None else [sorted(set(c)) for c in covariate]
        return {"estimable": False, "status": "not_computed", "reason": "stub"}
    monkeypatch.setattr(core, "_calibrated_mic_interval", record)
    raw = planted_cohort(tmp_path, mic=True)
    meta = pd.read_csv(tmp_path / "data" / "meta.csv")
    meta["laboratory"] = np.where(meta["iid"].str.endswith(("0", "5")), "B", "A")
    meta.to_csv(tmp_path / "data" / "meta.csv", index=False)
    raw["dataset"]["mic_covariate_columns"] = ["laboratory", "intake"]
    _, cfg = write_config(tmp_path, raw)
    out = core.run(cfg, results_dir=tmp_path / "out", seed=5)
    entry = out["metadata_diagnostics"]["censored_share"]["per_agent"]["agentA"]
    assert seen['covariate'] == [["A", "B"], ["2010", "2011", "2012"]]
    assert entry["calibrated"]["covariate_columns"] == ["laboratory", "intake"]


def test_stratify_by_repeats_the_run_within_each_level(tmp_path, monkeypatch):
    monkeypatch.setattr(core, "_calibrated_mic_interval",
                        lambda *a, **k: {"estimable": False, "status": "not_computed", "reason": "stub"})
    raw = planted_cohort(tmp_path, mic=True)
    meta = pd.read_csv(tmp_path / "data" / "meta.csv")
    meta["laboratory"] = np.where(meta["iid"].str.endswith(("0", "5")), "B", "A")
    meta.loc[meta.index[:2], "laboratory"] = None
    meta.to_csv(tmp_path / "data" / "meta.csv", index=False)
    raw["dataset"]["stratify_by"] = "laboratory"
    _, cfg = write_config(tmp_path, raw)
    out = core.run(cfg, results_dir=tmp_path / "out", seed=5)
    strata = out["strata"]
    assert strata["column"] == "laboratory" and strata["n_without_level"] == 2
    assert set(strata["levels"]) == {"A", "B"}
    for level in strata["levels"].values():
        assert level["n_isolates"] > 0
        assert "clonal_share" in level["metadata_diagnostics"]
        assert "censored_share" in level["metadata_diagnostics"]
        assert level["results"]
    assert sum(level["n_isolates"] for level in strata["levels"].values()) + 2 == len(meta)
    csv = (tmp_path / "out" / "strata_results.csv").read_text(encoding="utf-8")
    assert csv.lstrip("\ufeff").startswith("stratum,") and "\nA," in csv and "\nB," in csv
    report = (tmp_path / "out" / "report.md").read_text(encoding="utf-8")
    assert "Table 5d" in report and "strata_results.csv" in report
    raw["dataset"]["stratify_by"] = "no_such_column"
    with pytest.raises(ConfigError):
        core.run(write_config(tmp_path, raw)[1], results_dir=tmp_path / "out2", seed=5)


def test_screening_budget_applies_to_every_agent_but_the_named_ones(tmp_path, monkeypatch):
    seen = {}

    def record(lo, hi, lineage, cen, seed, wells=None, panel=None, covariate=None, n_boot=None):
        seen[len(seen)] = n_boot
        return {"estimable": False, "status": "not_computed", "reason": "stub"}
    monkeypatch.setattr(core, "_calibrated_mic_interval", record)
    raw = planted_cohort(tmp_path, mic=True)
    raw["censored"] = {"calibration_n_boot": 199,
                       "screening_n_boot": 19, "calibration_agents": ["agentA"]}
    _, cfg = write_config(tmp_path, raw)
    core.run(cfg, results_dir=tmp_path / "out", seed=5)
    budgets = sorted(seen.values())
    assert budgets[-1] == 199 and set(budgets[:-1]) == {19}
    cen = CensoredConfig(calibration_n_boot=199, screening_n_boot=19, calibration_agents=("agentA",))
    assert cen.draws_for("agentA") == 199 and cen.draws_for("agentB") == 19
    assert CensoredConfig(calibration_n_boot=199).draws_for("agentB") == 199
    with pytest.raises(ConfigError):
        from_dict(dict(raw, censored={"screening_n_boot": 5}), tmp_path)
    with pytest.raises(ConfigError):
        from_dict(dict(raw, censored={"calibration_agents": ["agentA"]}), tmp_path)


def test_dilution_table_counts_every_reading_once():
    lo, hi, labels = small_panel()
    table = core._dilution_table(lo, hi, labels)
    counts = np.asarray(table["counts"])
    assert counts.sum() == lo.size
    assert table["lineages"] == sorted(set(labels))
    assert table["intervals"][0][0] is None
    assert all(a is None or b is None or a < b for a, b in table["intervals"])
    assert len(table["intervals"]) == counts.shape[1]


def test_report_shows_panels_coefficients_and_status(tmp_path, monkeypatch):
    from amr_clonalshare.report import render_report
    from amr_clonalshare.report_html import render_html_report
    fake = {"estimate": 0.4, "low": 0.2, "high": 0.6, "confidence": 0.95, "estimable": True,
            "method": "nullwise_parametric_bootstrap_LR_inversion", "status": "computed",
            "n": 30, "n_groups": 6, "n_boot": 19, "seed": 5, "bootstrap_failed": 0,
            "tested_values": 12, "reason": "",
            "panel_edges": {"A": [-2.0, -1.0, 0.0], "B": [-1.0, 0.0, 1.0]},
            "adjustment": [{"levels": ["A", "B"], "reference_level": "A",
                            "coefficients": {"B": 0.75}}]}
    monkeypatch.setattr(core, "_calibrated_mic_interval", lambda *a, **k: dict(fake))
    raw = planted_cohort(tmp_path, mic=True)
    mic = pd.read_csv(tmp_path / "data" / "mic.csv")
    mic["laboratory"] = np.where(mic["iid"].str.endswith(("0", "5")), "B", "A")
    mic.to_csv(tmp_path / "data" / "mic.csv", index=False)
    raw["dataset"]["mic_panel_column"] = "laboratory"
    raw["dataset"]["mic_covariate_column"] = "laboratory"
    _, cfg = write_config(tmp_path, raw)
    out = core.run(cfg, results_dir=tmp_path / "out", seed=5)
    md = render_report(out, out.get("summary") or {})
    for label in ("Table 5a", "Table 5b", "Table 5c"):
        assert label in md
    assert "| agentA | laboratory | A | B | 0.750 |" in md
    assert "| agentA | computed | nullwise_parametric_bootstrap_LR_inversion | 19 | 0 | interval reported |" in md
    html = render_html_report(out, out.get("summary") or {})
    assert "Readings per lineage and dilution interval" in html and "<rect" in html


def test_screening_agents_are_checked_against_the_mic_table(tmp_path):
    from amr_clonalshare.io import load_dataset
    raw = planted_cohort(tmp_path, mic=True)
    raw["censored"] = {"screening_n_boot": 19, "calibration_agents": "agentA"}
    with pytest.raises(ConfigError, match="list of labels"):
        write_config(tmp_path, raw)
    raw["censored"]["calibration_agents"] = ["agentX"]
    with pytest.raises(ConfigError, match="agentX"):
        load_dataset(write_config(tmp_path, raw)[1])


def test_panel_labels_are_part_of_the_reading(tmp_path):
    from amr_clonalshare.io import load_dataset
    raw = planted_cohort(tmp_path, mic=True)
    raw["dataset"]["mic_panel_column"] = "laboratory"
    mic = pd.read_csv(tmp_path / "data" / "mic.csv")
    mic["laboratory"] = "A"
    # the same reading recorded twice under two laboratories is a conflict
    twice = pd.concat([mic, mic.iloc[:1].assign(laboratory="B")], ignore_index=True)
    twice.to_csv(tmp_path / "data" / "mic.csv", index=False)
    with pytest.raises(ConfigError, match="conflict"):
        load_dataset(write_config(tmp_path, raw)[1])
    # a reading without a laboratory is refused before any estimate
    blank = mic.copy(); blank.loc[0, "laboratory"] = "NA"
    blank.to_csv(tmp_path / "data" / "mic.csv", index=False)
    with pytest.raises(ConfigError, match="without a value"):
        load_dataset(write_config(tmp_path, raw)[1])


def test_strata_list_every_requested_analysis_and_their_own_join(tmp_path, monkeypatch):
    monkeypatch.setattr(core, "_calibrated_mic_interval",
                        lambda *a, **k: {"estimable": False, "status": "not_computed", "reason": "stub"})
    raw = planted_cohort(tmp_path, mic=True)
    meta = pd.read_csv(tmp_path / "data" / "meta.csv")
    meta["laboratory"] = np.where(meta["iid"].str.endswith(("0", "5")), "B", "A")
    meta.loc[meta.index[:3], "laboratory"] = "NA"   # a missing token, not a level
    meta.to_csv(tmp_path / "data" / "meta.csv", index=False)
    raw["dataset"]["stratify_by"] = "laboratory"
    out = core.run(write_config(tmp_path, raw)[1], results_dir=tmp_path / "out", seed=5)
    strata = out["strata"]
    assert set(strata["levels"]) == {"A", "B"} and strata["n_without_level"] == 3
    pooled = {(r["agent"], r["analysis"]) for r in out["results"]}
    for level in strata["levels"].values():
        assert {(r["agent"], r["analysis"]) for r in level["results"]} == pooled
        join = level["metadata_diagnostics"]["censored_share"]["join"]
        assert join["strains_aligned"] == level["n_isolates"]
    raw["dataset"]["stratify_by"] = ["laboratory"]
    with pytest.raises(ConfigError, match="one column"):
        write_config(tmp_path, raw)
