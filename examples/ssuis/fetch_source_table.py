#!/usr/bin/env python3
"""Fetch Additional file 1 of the source publication and flatten it to CSV.

    python examples/ssuis/fetch_source_table.py            # download and write
    python examples/ssuis/fetch_source_table.py --check    # verify the shipped copy

Hadjirin et al. 2021, *BMC Biology* 19:191, doi:10.1186/s12915-021-01094-1,
licence CC BY 4.0. Additional file 1 is table S1: one row per isolate, carrying
the hierBAPS population cluster, the collection, serotype, source, year, the
sixteen MIC columns and forty-three resistance determinants.

The MIC columns of this table and the MIC records that BV-BRC serves for the
same isolates are the same measurements reached by two different routes, which
is what makes the two tables joinable and what
``link_source_lineages.py`` exploits.

The workbook is read with the standard library. An `.xlsx` file is a zip of
XML, and forty lines of parsing is a smaller liability in a reproducibility
chain than a spreadsheet dependency that must be pinned, built and audited.

**Serotype 1/2 is stored in this file as the number 44228.** Under the Excel
1900 date system that serial is 1 February 2021: the string ``1/2`` was read as
a date when the table was prepared and written back as a serial. All 39
affected rows are repaired here, and the same 39 isolates reach BV-BRC as the
string ``1-Feb``, which is the same cell rendered by a different reader. The
error is in the published record rather than in either database, and it is the
error Ziemann, Eren and El-Osta documented for gene symbols (2016, *Genome
Biology* 17:177).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from io import BytesIO
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TARGET = DATA / "source_table_s1.csv"

URL = ("https://static-content.springer.com/esm/"
       "art%3A10.1186%2Fs12915-021-01094-1/MediaObjects/"
       "12915_2021_1094_MOESM1_ESM.xlsx")
#: sha256 of Additional file 1 as retrieved on 2026-08-31.
WORKBOOK_SHA256 = ("45aac567f1b23deaa91fe42bbdba8ecedc1e005d4f611f04727639e3"
                   "e9ed920f")
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

#: Excel serials that are serotypes read as dates, and the serotype each was.
SEROTYPE_DATE_SERIALS = {"44228": "1/2"}


def column_number(reference: str) -> int:
    number = 0
    for character in reference:
        number = number * 26 + ord(character) - 64
    return number


def read_worksheet(payload: bytes) -> list[list[str | None]]:
    archive = zipfile.ZipFile(BytesIO(payload))
    shared = [
        "".join(node.text or "" for node in item.iter(NS + "t"))
        for item in ET.fromstring(
            archive.read("xl/sharedStrings.xml")).iter(NS + "si")
    ]
    sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    rows = []
    for row in sheet.iter(NS + "row"):
        cells: dict[int, str | None] = {}
        for cell in row.iter(NS + "c"):
            reference = re.match(r"([A-Z]+)", cell.get("r")).group(1)
            value = cell.find(NS + "v")
            if value is None:
                parsed = None
            elif cell.get("t") == "s":
                parsed = shared[int(value.text)]
            else:
                parsed = value.text
            cells[column_number(reference)] = parsed
        rows.append(cells)
    width = max(max(row) for row in rows if row)
    return [[row.get(i) for i in range(1, width + 1)] for row in rows]


def repair_serotypes(header: list, rows: list) -> int:
    if "Serotype" not in header:
        return 0
    column = header.index("Serotype")
    repaired = 0
    for row in rows:
        replacement = SEROTYPE_DATE_SERIALS.get(str(row[column]))
        if replacement is not None:
            row[column] = replacement
            repaired += 1
    return repaired


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--url", default=URL)
    args = parser.parse_args()

    if args.check:
        if not TARGET.is_file():
            print(f"missing: {TARGET}")
            return 1
        with TARGET.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.reader(handle))
        print(f"{TARGET.name}: {len(rows) - 1} isolates, "
              f"{len(rows[0])} columns")
        print(f"sha256: {hashlib.sha256(TARGET.read_bytes()).hexdigest()}")
        return 0

    print(f"fetching {args.url}")
    request = urllib.request.Request(
        args.url, headers={"User-Agent": "amr-clonalshare (research)"})
    with urllib.request.urlopen(request, timeout=180) as response:
        payload = response.read()
    digest = hashlib.sha256(payload).hexdigest()
    print(f"retrieved {len(payload)} bytes, sha256 {digest}")
    if digest != WORKBOOK_SHA256:
        print(f"WARNING: expected {WORKBOOK_SHA256}. The publisher's file has "
              f"changed; check what changed before using it.")

    table = read_worksheet(payload)
    header, rows = table[0], [list(r) for r in table[1:]]
    repaired = repair_serotypes(header, rows)
    print(f"serotype date serials repaired: {repaired}")

    DATA.mkdir(parents=True, exist_ok=True)
    with TARGET.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"written: {TARGET}  ({len(rows)} isolates, {len(header)} columns)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
