"""Derive the tables of the two further E. coli analyses from the shipped tables.

    python examples/ecoli_swine/build_derived_tables.py

Reads data/metadata.csv and data/mic_two_categories.csv, which build_tables.py
writes from the source workbook, and writes:

  data/metadata_periods.csv     the metadata with a column `period`, 2011-2013
                                or 2014-2017 from the year of isolation (blank
                                for the 11 isolates without a year, none of
                                which has readings)
  data/calls_two_categories.csv one 0/1 call per isolate and agent: 1 when the
                                MIC lies above the cut of mic_two_categories.csv
  input.csv                     one row per isolate with readings: the
                                ciprofloxacin call, the sequence type, the
                                phylogroup as published and the phylogroup
                                with ST23, ST88 and ST410 in phylogroup C, as
                                the scheme of Clermont et al. (2013) places
                                them, for the comparison of lineage definitions
  data/derived_receipt.json     counts and the SHA-256 of every file read and written

No value is changed: a call is the two-category reading written as 0 or 1.
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
PERIODS = (("2011-2013", range(2011, 2014)), ("2014-2017", range(2014, 2018)))
COMPARED_AGENT = "ciprofloxacin"
#: Sequence types of clonal complex 23, which the source codes as phylogroup
#: B1 and the scheme of Clermont et al. (2013) places in phylogroup C.
PHYLOGROUP_C = ("ST23", "ST88", "ST410")


def read(path: Path) -> list:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write(path: Path, rows: list) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def period(year: str) -> str:
    """The period of a year of isolation; blank where the source gives none."""
    if not year.strip().isdigit():
        return ""
    for name, years in PERIODS:
        if int(year) in years:
            return name
    raise ValueError(f"year {year} lies in no period")


def main() -> int:
    metadata = read(DATA / "metadata.csv")
    readings = read(DATA / "mic_two_categories.csv")
    with_period = [dict(row, period=period(row["isolation_year"])) for row in metadata]
    calls = [dict(isolate_id=r["isolate_id"], agent=r["agent"],
                  call="1" if r["operator"] == ">" else "0") for r in readings]
    if any(r["operator"] not in (">", "<=") for r in readings):
        raise ValueError("a two-category reading without > or <=")
    by_id = {row["isolate_id"]: row for row in metadata}
    compared = [dict(isolate_id=c["isolate_id"], ciprofloxacin_call=c["call"],
                     sequence_type=by_id[c["isolate_id"]]["st"],
                     phylogroup=by_id[c["isolate_id"]]["phylogroup"],
                     phylogroup_clermont2013=("C" if by_id[c["isolate_id"]]["st"] in PHYLOGROUP_C
                                              else by_id[c["isolate_id"]]["phylogroup"]))
                for c in calls if c["agent"] == COMPARED_AGENT]
    outputs = {DATA / "metadata_periods.csv": with_period, DATA / "calls_two_categories.csv": calls,
               HERE / "input.csv": compared}
    for path, rows in outputs.items():
        write(path, rows)
    analysed = {c["isolate_id"] for c in calls}
    receipt = dict(
        read={p.name: sha256(p) for p in (DATA / "metadata.csv", DATA / "mic_two_categories.csv")},
        written={p.relative_to(HERE).as_posix(): sha256(p) for p in outputs},
        isolates_with_calls=len(analysed),
        isolates_per_period={name: sum(1 for r in with_period if r["period"] == name and r["isolate_id"] in analysed)
                             for name, _ in PERIODS},
        calls_above_cut={a: sum(1 for c in calls if c["agent"] == a and c["call"] == "1")
                         for a in sorted({c["agent"] for c in calls})},
        compared_agent=COMPARED_AGENT)
    (DATA / "derived_receipt.json").write_text(json.dumps(receipt, indent=1) + "\n")
    print(json.dumps(receipt, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
