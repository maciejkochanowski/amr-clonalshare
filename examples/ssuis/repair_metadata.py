#!/usr/bin/env python3
"""Repair two defects in the BV-BRC metadata dump, and derive the contrast keys.

Run once; it is idempotent and prints what it changed.

    python examples/ssuis/repair_metadata.py --check     # report, change nothing
    python examples/ssuis/repair_metadata.py             # rewrite metadata.csv

**Defect 1: serotype 1/2 arrives as the string "1-Feb".** *Streptococcus suis*
serotype 1/2 is a recognised serotype, distinct from 1 and from 2. Somewhere
upstream of the public record the string "1/2" passed through a spreadsheet,
which read it as a date and wrote it back as "1-Feb". The corruption is in the
BV-BRC dump itself rather than in anything this package did, and it affects 39
of 677 isolates. It is the textbook shape of the error Ziemann, Eren and
El-Osta documented for gene symbols (2016, Genome Biology 17:177): a value that
looks like a date to a spreadsheet stops being the value it was.

**Defect 2: free-text serotyping results.** Several rows carry human notes
rather than a serotype - "1 + 6", "8 (+2 partial)", and one row reading
"(26 partia/31 partiall)". These are left exactly as received and only
recorded, because a note about a partial reaction is information and guessing
which serotype it meant would not be.

**Derived: the contrast keys.** ``collection_period`` splits the collection at
2013, and ``country_period`` joins it to the country. Neither is a repair; both
exist so that a two-collection contrast can be named in a YAML configuration
instead of being computed inside an example script. The 2013 boundary is the
gap in the United Kingdom series, which runs 2009-2011 and then 2013-2015 with
nothing in 2012.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
METADATA = HERE / "data" / "metadata.csv"

#: Spreadsheet date coercions seen in this dump, and the serotype each was.
DATE_COERCIONS = {"1-Feb": "1/2", "2-Jan": "1/2"}

PERIOD_BOUNDARY = 2013


def repair(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    out = frame.copy()
    changes: dict = {}

    serovar = out["serovar"].astype("string")
    coerced = serovar.isin(DATE_COERCIONS)
    changes["serovar_date_coercions_repaired"] = int(coerced.sum())
    changes["serovar_repair_map"] = {k: v for k, v in DATE_COERCIONS.items()
                                     if (serovar == k).any()}
    out["serovar"] = serovar.replace(DATE_COERCIONS)

    free_text = out["serovar"].astype("string").str.contains(
        r"[+/]|partial|and ", case=False, na=False)
    free_text &= out["serovar"].astype("string") != "1/2"
    changes["serovar_free_text_rows_left_as_received"] = int(free_text.sum())

    year = pd.to_numeric(out["collection_year"], errors="coerce")
    period = pd.Series(pd.NA, index=out.index, dtype="string")
    period[year.notna() & (year < PERIOD_BOUNDARY)] = "early"
    period[year >= PERIOD_BOUNDARY] = "late"
    out["collection_period"] = period
    out["country_period"] = (
        out["isolation_country"].astype("string") + " " + period)
    changes["collection_period_counts"] = (
        period.value_counts(dropna=False).to_dict())
    changes["country_period_counts"] = (
        out["country_period"].value_counts(dropna=False).head(8).to_dict())
    return out, changes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true",
                        help="report what would change and exit non-zero if "
                             "the file is not already repaired")
    args = parser.parse_args()

    frame = pd.read_csv(METADATA, dtype=str)
    repaired, changes = repair(frame)
    for key, value in changes.items():
        print(f"{key}: {value}")

    already = (list(frame.columns) == list(repaired.columns)
               and frame.astype(str).equals(repaired.astype(str)))
    if args.check:
        print("state:", "already repaired" if already else "repair pending")
        return 0 if already else 1
    if already:
        print("nothing to write")
        return 0
    repaired.to_csv(METADATA, index=False)
    print(f"written: {METADATA}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
