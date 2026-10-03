"""The pipeline reads every laboratory on its own panel, scores within the
declared strata, and repeats the whole run within the levels of a column."""
import os

import numpy as np
import pandas as pd
import pytest

from amr_clonalshare import core
from amr_clonalshare.config import ConfigError

from conftest import planted_cohort, write_config


def _laboratory(frame):
    return np.where(frame["iid"].str.endswith(("0", "5")), "B", "A")


@pytest.mark.skipif(os.name == "nt", reason="POSIX permission bits")
def test_published_results_directory_follows_the_umask(tmp_path):
    raw = planted_cohort(tmp_path)
    _, cfg = write_config(tmp_path, raw)
    mask = os.umask(0o022)
    try:
        core.run(cfg, results_dir=tmp_path / "out", seed=5)
    finally:
        os.umask(mask)
    assert (tmp_path / "out").stat().st_mode & 0o777 == 0o755
    assert (tmp_path / "out" / "report.md").stat().st_mode & 0o777 == 0o644


def test_panel_column_reads_each_laboratory_on_its_own_wells(tmp_path):
    raw = planted_cohort(tmp_path, mic=True)
    mic = pd.read_csv(tmp_path / "data" / "mic.csv")
    # laboratory B never records the lowest well, so its lowest reading is
    # an end well of B and an interior well of the pooled lattice
    mic["laboratory"] = _laboratory(mic)
    mic.loc[mic["laboratory"].eq("B") & mic["measurement"].eq(0.25), "measurement"] = 0.5
    mic.to_csv(tmp_path / "data" / "mic.csv", index=False)
    raw["dataset"]["mic_panel_column"] = "laboratory"
    _, cfg = write_config(tmp_path, raw)
    out = core.run(cfg, results_dir=tmp_path / "out", seed=5)
    entry = out["metadata_diagnostics"]["censored_share"]["per_agent"]["agentA"]
    assert entry["panel"]["panel_column"] == "laboratory"
    assert set(entry["panel"]["per_panel"]) == {"A", "B"}
    assert entry["panel"]["per_panel"]["B"]["lowest"] == 0.5
    assert entry["panel"]["per_panel"]["A"]["lowest"] == 0.25
    # the readings are scored and permuted within laboratory
    assert entry["order"]["stratified_by"] == "laboratory" and entry["order"]["n_strata"] == 2
    report = (tmp_path / "out" / "report.md").read_text(encoding="utf-8")
    assert "Table 5a" in report and "| agentA | B |" in report
    raw["dataset"]["mic_panel_column"] = "no_such_column"
    with pytest.raises(ConfigError):
        core.run(write_config(tmp_path, raw)[1], results_dir=tmp_path / "out2", seed=5)


@pytest.mark.parametrize("columns, expected", [(["laboratory"], "laboratory"),
                                               (["laboratory", "intake"], "laboratory × intake")])
def test_covariate_columns_join_the_strata(tmp_path, columns, expected):
    """A declared covariate, one or more, joins the laboratory in the strata
    the readings are scored and permuted within, so that a shift between its
    levels is not read as a lineage difference."""
    raw = planted_cohort(tmp_path, mic=True)
    meta = pd.read_csv(tmp_path / "data" / "meta.csv")
    meta["laboratory"] = _laboratory(meta)
    meta.to_csv(tmp_path / "data" / "meta.csv", index=False)
    raw["dataset"]["mic_covariate_columns"] = columns
    _, cfg = write_config(tmp_path, raw)
    out = core.run(cfg, seed=5)
    order = out["metadata_diagnostics"]["censored_share"]["per_agent"]["agentA"]["order"]
    assert order["stratified_by"] == expected
    assert order["n_strata"] == (2 if len(columns) == 1 else 6)


def test_dilution_table_counts_every_reading_once():
    rng = np.random.default_rng(8)
    labels = np.repeat(np.array([f"L{g}" for g in range(6)], dtype=object), 5)
    latent = rng.normal(size=6)[np.arange(6).repeat(5)] + rng.normal(size=30)
    cuts = np.array([-2., -1., 0., 1., 2.])
    k = np.searchsorted(cuts, latent)
    ext = np.r_[-np.inf, cuts, np.inf]
    lo, hi = ext[k], ext[k + 1]
    table = core._dilution_table(lo, hi, labels)
    counts = np.asarray(table["counts"])
    assert counts.sum() == lo.size
    assert table["lineages"] == sorted(set(labels))
    assert table["intervals"][0][0] is None
    assert all(a is None or b is None or a < b for a, b in table["intervals"])
    assert len(table["intervals"]) == counts.shape[1]


def test_the_html_report_draws_the_dilution_heat_map(tmp_path):
    from amr_clonalshare.report_html import render_html_report
    raw = planted_cohort(tmp_path, mic=True)
    _, cfg = write_config(tmp_path, raw)
    out = core.run(cfg, seed=5)
    html = render_html_report(out, out.get("summary") or {})
    assert "Readings per lineage and dilution interval" in html and "<rect" in html


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
    blank = mic.copy()
    blank.loc[0, "laboratory"] = "NA"
    blank.to_csv(tmp_path / "data" / "mic.csv", index=False)
    with pytest.raises(ConfigError, match="without a value"):
        load_dataset(write_config(tmp_path, raw)[1])


def test_stratify_by_repeats_the_run_within_each_level(tmp_path):
    raw = planted_cohort(tmp_path, mic=True)
    meta = pd.read_csv(tmp_path / "data" / "meta.csv")
    meta["laboratory"] = _laboratory(meta)
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
    assert csv.lstrip("﻿").startswith("stratum,") and "\nA," in csv and "\nB," in csv
    report = (tmp_path / "out" / "report.md").read_text(encoding="utf-8")
    assert "Table 5d" in report and "strata_results.csv" in report
    raw["dataset"]["stratify_by"] = "no_such_column"
    with pytest.raises(ConfigError):
        core.run(write_config(tmp_path, raw)[1], results_dir=tmp_path / "out2", seed=5)


def test_strata_list_every_requested_analysis_and_their_own_join(tmp_path):
    raw = planted_cohort(tmp_path, mic=True)
    meta = pd.read_csv(tmp_path / "data" / "meta.csv")
    meta["laboratory"] = _laboratory(meta)
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
