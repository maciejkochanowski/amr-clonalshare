"""The E. coli example: its tables follow from the receipt of build_tables.py,
its three configurations load, and the shipped records describe the tables."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from amr_clonalshare.config import load_config

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
EXAMPLE = EXAMPLES / "ecoli_swine"
# An installed wheel carries no examples folder at all, and the suite may run
# beside one; a source tree or distribution that has the folder must hold
# this example, so its absence there is a failure, not a skip.
pytestmark = pytest.mark.skipif(not EXAMPLES.is_dir(), reason="the examples folder is not beside the package")


def _rows(name: str):
    with (EXAMPLE / "data" / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def test_the_tables_match_the_build_receipt():
    receipt = json.loads((EXAMPLE / "data" / "build_receipt.json").read_text(encoding="utf-8"))
    meta = {r["isolate_id"]: r for r in _rows("metadata.csv")}
    full, two, without = _rows("mic_long.csv"), _rows("mic_two_categories.csv"), _rows("mic_long_without_st410.csv")
    assert len(meta) == receipt["isolates_in_data_1"] == 499
    assert sum(1 for r in meta.values() if not r["st"]) == 9
    isolates = {r["isolate_id"] for r in full}
    assert len(isolates) == receipt["ast_rows_in_data_13"] - len(receipt["set_aside"]) == 481
    assert not isolates & set(receipt["set_aside"])
    assert all(len([r for r in full if r["agent"] == agent]) == 481 for agent in receipt["wells"])
    # every reading sits on the declared range, and a sign marks only its ends
    for r in full:
        wells = receipt["wells"][r["agent"]]
        value = float(r["measurement"])
        assert value in wells
        assert r["operator"] in ("", "<=", ">")
        assert r["operator"] != "<=" or value == wells[0]
        assert r["operator"] != ">" or value == wells[-1]
    # the two-category table is the full table read against the cut of its agent
    by_key = {(r["isolate_id"], r["agent"]): r for r in full}
    for r in two:
        cut = receipt["two_category_cut"][r["agent"]]
        source = by_key[(r["isolate_id"], r["agent"])]
        value = float(source["measurement"])
        above = (source["operator"] == ">" and value >= cut) or (source["operator"] == "" and value > cut)
        assert float(r["measurement"]) == cut and r["operator"] == (">" if above else "<=")
    # the ST410-free table drops exactly the isolates of ST410
    st410 = {i for i, r in meta.items() if r["st"] == receipt["excluded_sequence_type"]}
    assert len(st410) == receipt["isolates_of_excluded_type"]
    assert {r["isolate_id"] for r in without} == isolates - st410
    assert len(isolates & st410) == 73


def test_the_configurations_load_and_the_records_describe_the_tables():
    for name, folder in (("config.yaml", "expected"), ("config_two_categories.yaml", "expected_two_categories"),
                         ("config_without_st410.yaml", "expected_without_st410")):
        cfg = load_config(EXAMPLE / name)
        assert cfg.dataset.lineage_column == "st" and set(cfg.dataset.mic_wells) == {
            "cefotaxime", "ceftazidime", "meropenem", "ciprofloxacin", "gentamicin", "tetracycline"}
        record = json.loads((EXAMPLE / folder / "clonal_share_result.json").read_text(encoding="utf-8"))
        assert record["seed"] == 20261001
        orders = record["metadata_diagnostics"]["censored_share"]["per_agent"]
        assert set(orders) == set(cfg.dataset.mic_wells)
        order = orders["cefotaxime"]["order"]
        expected_n = 408 if "without" in folder else 481
        assert record["n_isolates"] == expected_n
        assert order["n_singletons_set_aside"] == 27
        assert order["n"] - order["n_singletons_set_aside"] == (372 if "without" in folder else 445)


def test_the_derived_tables_follow_from_the_shipped_ones():
    receipt = json.loads((EXAMPLE / "data" / "derived_receipt.json").read_text(encoding="utf-8"))
    for name, digest in receipt["read"].items():
        assert hashlib.sha256((EXAMPLE / "data" / name).read_bytes()).hexdigest() == digest
    for name, digest in receipt["written"].items():
        assert hashlib.sha256((EXAMPLE / name).read_bytes()).hexdigest() == digest
    meta = {r["isolate_id"]: r for r in _rows("metadata_periods.csv")}
    two = {(r["isolate_id"], r["agent"]): r for r in _rows("mic_two_categories.csv")}
    calls = _rows("calls_two_categories.csv")
    assert len(calls) == len(two) == 6 * 481
    assert all(c["call"] == ("1" if two[(c["isolate_id"], c["agent"])]["operator"] == ">" else "0") for c in calls)
    for r in meta.values():
        year = r["isolation_year"]
        assert r["period"] == ("" if not year.isdigit() else "2011-2013" if int(year) <= 2013 else "2014-2017")
    assert receipt["isolates_per_period"] == {"2011-2013": 197, "2014-2017": 284}
    with (EXAMPLE / "input.csv").open(encoding="utf-8", newline="") as handle:
        compared = list(csv.DictReader(handle))
    assert len(compared) == 481
    assert all(r["ciprofloxacin_call"] == two[(r["isolate_id"], "ciprofloxacin")]["operator"].replace("<=", "0").replace(">", "1")
               and r["sequence_type"] == meta[r["isolate_id"]]["st"] and r["phylogroup"] == meta[r["isolate_id"]]["phylogroup"]
               for r in compared)
    # the second phylogroup column differs from the first only on clonal complex 23
    assert all(r["phylogroup_clermont2013"] == ("C" if r["sequence_type"] in ("ST23", "ST88", "ST410") else r["phylogroup"])
               for r in compared)


def test_the_periods_record_splits_the_change_by_phylogroup():
    cfg = load_config(EXAMPLE / "config_periods.yaml")
    assert cfg.dataset.lineage_column == "phylogroup" and cfg.dataset.contrast_column == "period"
    record = json.loads((EXAMPLE / "expected_periods" / "clonal_share_result.json").read_text(encoding="utf-8"))
    assert record["seed"] == 20261001 and record["n_isolates"] == 481
    split = record["metadata_diagnostics"]["prevalence_decomposition"]
    assert split["levels"] == ["2014-2017", "2011-2013"] and split["n_a"] == 284 and split["n_b"] == 197
    for item in split["per_feature"].values():
        assert item["shared_support_isolate_share"] == 1.0
        assert abs(item["composition"] + item["within_lineage"] - item["difference"]) < 1e-12


def test_the_comparison_record_reads_the_shipped_table():
    for folder, second in (("expected_comparison", "phylogroup"), ("expected_comparison_clermont", "phylogroup_clermont2013")):
        record = json.loads((EXAMPLE / folder / "comparison.json").read_text(encoding="utf-8"))
        assert record["input_sha256"] == hashlib.sha256((EXAMPLE / "input.csv").read_bytes()).hexdigest()
        assert record["names"] == {"a": "sequence_type", "b": second}
        assert record["counts"]["common"] == 472 and record["counts"]["common_scorable"] == 445
