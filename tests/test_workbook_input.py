"""A real workbook, read end to end.

The Excel path was shipped with a test that made `read_excel` raise, which
checks the refusal and nothing else. These read an `.xlsx` written by
`openpyxl` and require that a workbook and the CSV of the same sheet give the
same configuration draft and the same numbers, because that is the promise the
manual makes.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

openpyxl = pytest.importorskip("openpyxl")

ROOT = Path(__file__).resolve().parent.parent


def _write_workbook(path: Path, rows: list[list[str]]) -> None:
    book = openpyxl.Workbook()
    sheet = book.active
    for row in rows:
        sheet.append(row)
    book.save(path)


def _tables():
    metadata = [["isolate_id", "lineage"]]
    calls = [["isolate_id", "agent", "call"]]
    pattern = "RSRSS"
    for i in range(40):
        metadata.append([f"iso_{i:02d}", f"L{i // 4 + 1:02d}"])
        calls.append([f"iso_{i:02d}", "ag_alpha", pattern[i % 5]])
        calls.append([f"iso_{i:02d}", "ag_beta", pattern[(i + 2) % 5]])
    return metadata, calls


CONFIG = """\
dataset:
  name: workbook_fixture
  data_dir: data
  metadata: metadata.{ext}
  strain_id_column: isolate_id
  lineage_column: lineage
  phenotype: calls.{ext}
  phenotype_id_column: isolate_id
  phenotype_antibiotic_column: agent
  phenotype_call_column: call
  phenotype_kind: clinical_sir
attribution: {{folds: 2, repeats: 2, n_boot: 20, n_perm: 19}}
surveillance: {{n_boot: 20}}
evidence: {{folds: 2, repeats: 2}}
"""


def _case(root: Path, ext: str) -> Path:
    data = root / "data"
    data.mkdir(parents=True, exist_ok=True)
    metadata, calls = _tables()
    if ext == "csv":
        for name, rows in (("metadata", metadata), ("calls", calls)):
            (data / f"{name}.csv").write_text(
                "\n".join(",".join(row) for row in rows) + "\n", encoding="utf-8")
    else:
        _write_workbook(data / "metadata.xlsx", metadata)
        _write_workbook(data / "calls.xlsx", calls)
    config = root / "config.yaml"
    config.write_text(CONFIG.format(ext=ext), encoding="utf-8")
    return config


def _run(config: Path, out: Path) -> dict:
    from amr_clonalshare import cli
    assert cli.main(["--config", str(config), "--results-dir", str(out),
                     "--quiet"]) == 0
    return json.loads((out / "clonal_share_result.json").read_text(encoding="utf-8"))


def test_a_workbook_and_the_csv_of_the_same_sheet_give_the_same_numbers(tmp_path):
    from_csv = _run(_case(tmp_path / "csv", "csv"), tmp_path / "out_csv")
    from_book = _run(_case(tmp_path / "xlsx", "xlsx"), tmp_path / "out_xlsx")
    for block in ("clonal_share", "realised_share"):
        assert (from_book["metadata_diagnostics"][block]
                == from_csv["metadata_diagnostics"][block])
    assert (from_book["metadata_diagnostics"]["lineage_evidence"]["per_feature"]
            == from_csv["metadata_diagnostics"]["lineage_evidence"]["per_feature"])
    assert from_book["input_qc"]["phenotype_reading"]["rows_read"] == \
        from_csv["input_qc"]["phenotype_reading"]["rows_read"]


def test_the_draft_reads_a_workbook_as_it_reads_the_csv(tmp_path):
    from amr_clonalshare.draft import draft_config

    metadata, calls = _tables()
    csv_path = tmp_path / "calls.csv"
    csv_path.write_text("\n".join(",".join(row) for row in calls) + "\n",
                        encoding="utf-8")
    book_path = tmp_path / "calls.xlsx"
    _write_workbook(book_path, calls)
    from_csv = draft_config([csv_path]).replace("calls.csv", "calls.<ext>")
    from_book = draft_config([book_path]).replace("calls.xlsx", "calls.<ext>")
    assert from_book == from_csv


def test_a_workbook_with_one_column_per_agent_is_read_in_that_shape(tmp_path):
    from amr_clonalshare.draft import draft_config

    rows = [["isolate_id", "ampicillin", "tetracycline"]]
    for i in range(12):
        rows.append([f"iso_{i:02d}", "S" if i % 2 else "R", "R" if i % 3 else "S"])
    book = tmp_path / "wide.xlsx"
    _write_workbook(book, rows)
    drafted = draft_config([book])
    assert 'phenotype_agent_columns: ["ampicillin", "tetracycline"]' in drafted
    assert "phenotype_kind: clinical_sir" in drafted


def test_numeric_workbook_ids_that_lost_their_zeros_are_named(tmp_path):
    from amr_clonalshare.config import Config, DatasetConfig
    from amr_clonalshare.io import load_dataset
    _write_workbook(tmp_path / "calls.xlsx",
                    [["isolate_id", "agent", "call"], [1, "tet", "R"], [2, "tet", "S"]])
    (tmp_path / "meta.csv").write_text("isolate_id,lineage\n00001,A\n00002,B\n")
    cfg = Config(DatasetConfig(name="t", data_dir=str(tmp_path), metadata="meta.csv",
        strain_id_column="isolate_id", lineage_column="lineage", phenotype="calls.xlsx",
        phenotype_id_column="isolate_id", phenotype_antibiotic_column="agent",
        phenotype_call_column="call", phenotype_kind="clinical_sir"))
    with pytest.warns(RuntimeWarning, match="leading zeros"):
        ds = load_dataset(cfg)
    assert ds.input_qc["metadata_join"]["n_near_matches"] == 2
