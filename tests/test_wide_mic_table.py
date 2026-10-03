"""A wide MIC table, one column per agent, is read exactly as its long shape."""
import csv
import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

from amr_clonalshare.config import ConfigError, from_dict, load_config
from amr_clonalshare.core import run
from amr_clonalshare.io import load_dataset

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
pytestmark = pytest.mark.skipif(not EXAMPLES.is_dir(), reason="examples folder not shipped")


def _widen(long_path: Path, wide_path: Path, *, id_column, agent_column, value_column,
           operator_column=None, keep=()):
    """One row per isolate; the cell of an agent is its sign and value."""
    rows = {}
    agents = []
    with long_path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            agent = row[agent_column]
            if agent not in agents:
                agents.append(agent)
            cell = rows.setdefault(row[id_column], {id_column: row[id_column], **{k: row[k] for k in keep}})
            sign = row.get(operator_column, "") if operator_column else ""
            cell[agent] = f"{sign}{row[value_column]}"
    with wide_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[id_column, *keep, *agents])
        writer.writeheader()
        for cell in rows.values():
            writer.writerow({k: cell.get(k, "") for k in writer.fieldnames})
    return agents


def _record(path):
    """The record without the configuration that names the shape and the
    file, and without the join fields that repeat them."""
    record = json.loads(Path(path).read_text(encoding="utf-8"))
    record.pop("config", None)
    for join in (record.get("input_qc", {}).get("mic_join"),
                 record.get("metadata_diagnostics", {}).get("censored_share", {}).get("join")):
        if join:
            join.pop("path", None)
            join.pop("wide_shape", None)
    return json.dumps(record, sort_keys=True, default=str)


def test_the_wide_and_the_long_ecoli_table_load_to_the_same_readings(tmp_path):
    source = EXAMPLES / "ecoli_swine"
    long_cfg = load_config(source / "config.yaml")
    agents = _widen(source / "data" / "mic_long.csv", tmp_path / "mic_wide.csv", id_column="isolate_id",
                    agent_column="agent", value_column="measurement", operator_column="operator",
                    keep=("unit", "panel"))
    raw = yaml.safe_load((source / "config.yaml").read_text(encoding="utf-8"))
    ds = raw["dataset"]
    ds["data_dir"] = str(source / "data")
    ds["mic"] = str(tmp_path / "mic_wide.csv")
    ds["mic_agent_columns"] = agents
    for key in ("mic_antibiotic_column", "mic_value_column", "mic_operator_column"):
        ds.pop(key, None)
    wide_cfg = from_dict(raw)
    long_ds = load_dataset(long_cfg)
    wide_ds = load_dataset(wide_cfg)
    keys = ["isolate_id", "antibiotic", "measurement", "_operator_from_value"]
    long_mic = long_ds.mic.rename(columns={"agent": "antibiotic", "operator": "_operator_from_value"})
    long_mic = long_mic.sort_values(["isolate_id", "antibiotic"]).reset_index(drop=True)[keys + ["panel"]]
    wide_mic = wide_ds.mic.sort_values(["isolate_id", "antibiotic"]).reset_index(drop=True)[keys + ["panel"]]
    pd.testing.assert_frame_equal(long_mic, wide_mic)
    assert wide_ds.mic_join["wide_shape"] is True
    assert wide_ds.mic_join["antimicrobials"] == long_ds.mic_join["antimicrobials"]
    assert wide_ds.mic_join["strains_with_mic"] == long_ds.mic_join["strains_with_mic"]


def test_the_wide_and_the_long_fictional_collection_give_one_record(tmp_path):
    source = EXAMPLES / "workflows"
    agents = _widen(source / "data" / "mic.csv", tmp_path / "mic_wide.csv", id_column="isolate_id",
                    agent_column="agent", value_column="measurement", operator_column="operator")
    raw = yaml.safe_load((source / "mic.yaml").read_text(encoding="utf-8"))
    raw["dataset"]["data_dir"] = str(source / "data")
    long_cfg = from_dict(raw)
    ds = raw["dataset"]
    ds["mic"] = str(tmp_path / "mic_wide.csv")
    ds["mic_agent_columns"] = agents
    for key in ("mic_antibiotic_column", "mic_value_column", "mic_operator_column"):
        ds.pop(key, None)
    wide_cfg = from_dict(raw)
    run(long_cfg, results_dir=tmp_path / "long", seed=7)
    run(wide_cfg, results_dir=tmp_path / "wide", seed=7)
    assert _record(tmp_path / "long" / "clonal_share_result.json") == _record(tmp_path / "wide" / "clonal_share_result.json")


def test_a_wide_table_refuses_an_operator_column_and_a_missing_agent_column(tmp_path):
    (tmp_path / "metadata.csv").write_text("id,st\nA1,ST1\nA2,ST1\nA3,ST2\nA4,ST2\n", encoding="utf-8")
    (tmp_path / "mic.csv").write_text("id,CIP,GEN\nA1,<=0.5,2\nA2,1,4\nA3,2,>8\nA4,0.5,4\n", encoding="utf-8")
    base = {"name": "t", "data_dir": str(tmp_path), "metadata": "metadata.csv", "strain_id_column": "id",
            "lineage_column": "st", "mic": "mic.csv", "mic_id_column": "id"}
    with pytest.raises(ConfigError, match="no place in a wide MIC table"):
        from_dict({"dataset": {**base, "mic_agent_columns": ["CIP", "GEN"], "mic_operator_column": "op"}})
    with pytest.raises(ConfigError, match="no column"):
        load_dataset(from_dict({"dataset": {**base, "mic_agent_columns": ["CIP", "TET"]}}))
    ds = load_dataset(from_dict({"dataset": {**base, "mic_agent_columns": ["CIP", "GEN"]}}))
    # the WHONET codes are read as the agents they stand for, and the record says so
    assert ds.mic_join["antimicrobials"] == ["ciprofloxacin", "gentamicin"]
    assert ds.mic_join["agent_names_unified"] == {"CIP": "ciprofloxacin", "GEN": "gentamicin"}


def test_one_sheet_with_lineage_and_wide_mics_gives_the_numbers_of_the_two_file_version(tmp_path):
    source = EXAMPLES / "workflows"
    agents = _widen(source / "data" / "mic.csv", tmp_path / "mic_wide.csv", id_column="isolate_id",
                    agent_column="agent", value_column="measurement", operator_column="operator")
    metadata = pd.read_csv(source / "data" / "metadata.csv", dtype=str, keep_default_na=False)
    wide = pd.read_csv(tmp_path / "mic_wide.csv", dtype=str, keep_default_na=False)
    sheet = metadata.merge(wide, on="isolate_id", how="left")
    sheet.to_csv(tmp_path / "one_sheet.csv", index=False)
    raw = yaml.safe_load((source / "mic.yaml").read_text(encoding="utf-8"))
    raw["dataset"]["data_dir"] = str(source / "data")
    long_cfg = from_dict(raw)
    ds = raw["dataset"]
    ds["metadata"] = ds["mic"] = str(tmp_path / "one_sheet.csv")
    ds["mic_agent_columns"] = agents
    for key in ("mic_antibiotic_column", "mic_value_column", "mic_operator_column"):
        ds.pop(key, None)
    one_cfg = from_dict(raw)
    run(long_cfg, results_dir=tmp_path / "two_files", seed=7)
    run(one_cfg, results_dir=tmp_path / "one_sheet", seed=7)
    assert _record(tmp_path / "two_files" / "clonal_share_result.json") == _record(tmp_path / "one_sheet" / "clonal_share_result.json")
    # the draft of that one sheet names it as the metadata and the MIC table
    from amr_clonalshare.draft import draft_config
    draft = draft_config([tmp_path / "one_sheet.csv"])
    assert draft.count('"one_sheet.csv"') == 2
    assert 'lineage_column: "lineage"' in draft
    assert 'mic_agent_columns: ["demo_agent"]' in draft
    assert '"demo_agent": [0.25, 0.5, 1, 2, 4, 8]' in draft


def test_the_draft_of_the_ecoli_tables_runs_once_the_panel_range_is_typed_in(tmp_path):
    from amr_clonalshare.cli import main
    from amr_clonalshare.draft import draft_config
    source = EXAMPLES / "ecoli_swine" / "data"
    draft = draft_config([source / "metadata.csv", source / "mic_long.csv"])
    for line in ('mic: "mic_long.csv"', 'mic_antibiotic_column: "agent"', 'mic_value_column: "measurement"',
                 'mic_operator_column: "operator"', 'mic_panel_column: "panel"', 'lineage_column: "st"',
                 '"ciprofloxacin": [0.5, 1, 2]', "OBSERVED"):
        assert line in draft
    assert "phenotype:" not in draft
    # the ranges of DATA_PROVENANCE replace the observed ones; nothing else is edited
    panel = {"cefotaxime": "[1, 2, 4, 8, 16, 32]", "ceftazidime": "[1, 2, 4, 8, 16]", "meropenem": "[1, 2, 4, 8]",
             "ciprofloxacin": "[0.5, 1, 2]", "gentamicin": "[2, 4, 8]", "tetracycline": "[2, 4, 8]"}
    lines = []
    for line in draft.splitlines():
        for agent, wells in panel.items():
            if line.startswith(f'    "{agent}":'):
                line = f'    "{agent}": {wells}'
        lines.append(line)
    (tmp_path / "draft.yaml").write_text("\n".join(lines) + "\n", encoding="utf-8")
    assert main(["--config", str(tmp_path / "draft.yaml"), "--check-input", "--results-dir", str(tmp_path / "check")]) == 0
    assert (tmp_path / "check" / "input_qc.md").is_file()
