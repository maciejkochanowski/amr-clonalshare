"""A shipped panel preset fills the tested concentrations as typing them does."""
import json
import math

import pytest

from amr_clonalshare.censored import panel_geometry
from amr_clonalshare.config import ConfigError, from_dict
from amr_clonalshare.core import run
from amr_clonalshare.panels import PANEL_PRESETS, preset_names, preset_wells


def test_every_preset_is_a_doubling_series_from_a_public_source():
    assert len(preset_names()) == 4
    for name, preset in PANEL_PRESETS.items():
        assert "2020/1729" in preset["source"]
        for agent, wells in preset["wells"].items():
            assert wells == sorted(wells)
            assert all(math.isclose(b / a, 2.0) for a, b in zip(wells, wells[1:])), (name, agent)
    first = preset_wells("eu-2020-1729-salmonella-ecoli")
    assert first["ampicillin"] == [1, 2, 4, 8, 16, 32]
    assert first["ciprofloxacin"][0] == pytest.approx(1 / 64) and first["ciprofloxacin"][-1] == 8


def _tables(tmp_path):
    tmp_path.mkdir(exist_ok=True)
    (tmp_path / "metadata.csv").write_text(
        "id,st\n" + "\n".join(f"A{i},ST{i % 3}" for i in range(12)) + "\n", encoding="utf-8")
    values = ["<=0.015", "0.03", "0.06", "0.12", "0.25", "0.5", "1", "2", "4", ">8", "0.12", "0.25"]
    (tmp_path / "mic.csv").write_text(
        "id,agent,mic\n" + "\n".join(f"A{i},CIP,{v}" for i, v in enumerate(values)) + "\n", encoding="utf-8")
    return {"name": "t", "data_dir": str(tmp_path), "metadata": "metadata.csv", "strain_id_column": "id",
            "lineage_column": "st", "mic": "mic.csv", "mic_id_column": "id", "mic_antibiotic_column": "agent",
            "mic_value_column": "mic"}


def test_a_preset_fills_mic_wells_exactly_as_typing_the_range_does(tmp_path):
    base = _tables(tmp_path)
    typed = from_dict({"dataset": {**base, "mic_wells": {"ciprofloxacin": [0.015625, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8]}}})
    preset = from_dict({"dataset": {**base, "mic_panel_preset": "eu-2020-1729-salmonella-ecoli"}})
    assert preset.dataset.mic_wells["ciprofloxacin"] == typed.dataset.mic_wells["ciprofloxacin"]
    assert "ampicillin" in preset.dataset.mic_wells
    # an agent written explicitly keeps its own range
    own = from_dict({"dataset": {**base, "mic_panel_preset": "eu-2020-1729-salmonella-ecoli",
                                 "mic_wells": {"ciprofloxacin": [0.25, 0.5, 1, 2]}}})
    assert own.dataset.mic_wells["ciprofloxacin"] == [0.25, 0.5, 1, 2]
    with pytest.raises(ConfigError, match="not a shipped panel"):
        from_dict({"dataset": {**base, "mic_panel_preset": "sensititre-anything"}})


def test_readings_written_as_0_12_match_the_well_at_0_125_and_the_record_names_the_preset(tmp_path):
    base = _tables(tmp_path)
    cfg = from_dict({"dataset": {**base, "mic_panel_preset": "eu-2020-1729-salmonella-ecoli"},
                     "attribution": {"folds": 2, "repeats": 2, "n_boot": 10, "n_perm": 9},
                     "evidence": {"folds": 2, "repeats": 2}})
    run(cfg, results_dir=tmp_path / "out", seed=1)
    record = json.loads((tmp_path / "out" / "clonal_share_result.json").read_text(encoding="utf-8"))
    entry = record["metadata_diagnostics"]["censored_share"]["per_agent"]["ciprofloxacin"]
    assert entry["panel_source"].startswith("preset eu-2020-1729-salmonella-ecoli")
    assert entry["panel"]["n_wells"] == 10
    assert "shipped panel preset (eu-2020-1729-salmonella-ecoli)" in (tmp_path / "out" / "report.md").read_text(encoding="utf-8")


def test_a_reading_farther_than_a_quarter_doubling_from_every_well_is_refused():
    geometry = panel_geometry([0.12, 0.25, 0.5], wells=[0.125, 0.25, 0.5])
    assert geometry.n_wells == 3
    with pytest.raises(ValueError, match="outside configured wells"):
        panel_geometry([0.18, 0.25], wells=[0.125, 0.25, 0.5])
