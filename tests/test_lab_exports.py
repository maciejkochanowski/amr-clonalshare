"""Tables as laboratory systems and spreadsheets export them read to one table."""
import pandas as pd
import pytest

from amr_clonalshare.io import _read_table, text_separator


ROWS = [("A1", "ST10", "<=0.5", "8"), ("A2", "ST131", "1", ">16"), ("A3", "ST10", "2", "4")]
COLUMNS = ("isolate_id", "st", "ciprofloxacin", "gentamicin")


def _write(path, text, encoding="utf-8"):
    path.write_text(text, encoding=encoding)
    return path


@pytest.fixture
def exports(tmp_path):
    comma = "\n".join(",".join(c) for c in (COLUMNS, *ROWS)) + "\n"
    semicolon = "\n".join(";".join(c) for c in (COLUMNS, *ROWS)) + "\n"
    tab = "\n".join("\t".join(c) for c in (COLUMNS, *ROWS)) + "\n"
    spaced = "\n".join(",".join(c) for c in (("Isolate_ID ", " st", " ciprofloxacin", "gentamicin "), *ROWS)) + "\n"
    trailing = comma + ",,,\n,,,\n"
    return {
        "comma": _write(tmp_path / "comma.csv", comma),
        "semicolon": _write(tmp_path / "semicolon.csv", semicolon),
        "tab": _write(tmp_path / "tab.txt", tab),
        "bom": _write(tmp_path / "bom.csv", comma, encoding="utf-8-sig"),
        "spaced": _write(tmp_path / "spaced.csv", spaced),
        "trailing": _write(tmp_path / "trailing.csv", trailing),
    }


def test_six_exports_of_one_table_read_to_the_same_frame(exports):
    reference = _read_table(exports["comma"], "table", dtype=str)
    assert list(reference.columns) == list(COLUMNS)
    assert len(reference) == 3
    for name, path in exports.items():
        table = _read_table(path, "table", dtype=str)
        if name == "spaced":
            assert list(table.columns) == ["Isolate_ID", "st", "ciprofloxacin", "gentamicin"]
            table.columns = reference.columns
        pd.testing.assert_frame_equal(table, reference, check_names=False)


def test_the_separator_is_the_most_frequent_one_of_the_header(exports):
    assert text_separator(exports["comma"]) == ","
    assert text_separator(exports["semicolon"]) == ";"
    assert text_separator(exports["tab"]) == "\t"


def test_the_draft_readers_see_the_same_columns_and_rows(exports):
    from amr_clonalshare.draft import _columns, _rows
    for name, path in exports.items():
        columns = _columns(path)
        assert len(columns) == 4, name
        rows = _rows(path)
        assert len(rows) == 3, name
        assert rows[0][columns[0]] == "A1"


def test_a_decimal_comma_and_a_dash_are_read_in_a_mic_table(tmp_path):
    from amr_clonalshare.config import from_dict
    from amr_clonalshare.io import load_dataset
    (tmp_path / "metadata.csv").write_text("isolate_id;st\nA1;ST10\nA2;ST10\nA3;ST131\nA4;ST131\n", encoding="utf-8")
    (tmp_path / "mic.csv").write_text(
        "isolate_id;agent;mic\nA1;ciprofloxacin;<=0,5\nA2;ciprofloxacin;0,25\nA3;ciprofloxacin;1\n"
        "A4;ciprofloxacin;-\nA1;gentamicin;n.d.\nA2;gentamicin;2\nA3;gentamicin;4\nA4;gentamicin;8\n", encoding="utf-8")
    cfg = from_dict({"dataset": {"name": "t", "data_dir": str(tmp_path), "metadata": "metadata.csv",
                                 "strain_id_column": "isolate_id", "lineage_column": "st", "mic": "mic.csv",
                                 "mic_id_column": "isolate_id", "mic_antibiotic_column": "agent",
                                 "mic_value_column": "mic"}})
    ds = load_dataset(cfg)
    join = ds.mic_join
    assert join["n_unparseable_values"] == 0
    assert join["n_missing_values"] == 2
    values = ds.mic.set_index(["isolate_id", "agent"])["mic"]
    assert values[("A1", "ciprofloxacin")] == 0.5
    assert values[("A2", "ciprofloxacin")] == 0.25
    assert ds.mic.loc[ds.mic["isolate_id"].eq("A1") & ds.mic["agent"].eq("ciprofloxacin"), "_operator_from_value"].iloc[0] == "<="
