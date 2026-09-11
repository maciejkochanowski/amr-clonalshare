#!/usr/bin/env python3
"""Write the panel of this cohort in the long shape the run reads.

    python examples/ssuis/build_calls_long.py

``ribo.csv`` and ``cell.csv`` hold the 8,801 non-wild-type calls one column
per antimicrobial, as ``derive_ecoffs.py`` wrote them. The run takes the panel
in the shape a laboratory exports it in, one row per isolate and antimicrobial,
which is also the shape NCBI Pathogen Detection and BV-BRC serve. This script
turns the one into the other and adds nothing: ``tests/test_ssuis_provenance.py``
rebuilds ``data/calls_long.csv`` from the two wide files and fails on any
difference.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent / "data"
# The package's vocabulary for a binary call. The cohort's calls are
# non-wild-type against an epidemiological cut-off, not a clinical category:
# the words below are how the reader spells 1 and 0, and DATA_PROVENANCE.md
# states what was measured.
CALL = {0: "susceptible", 1: "non-susceptible"}


def build() -> pd.DataFrame:
    """The two wide call tables melted into one long panel, agents in table order."""
    ribo = pd.read_csv(DATA / "ribo.csv", dtype={"genome_id": str})
    cell = pd.read_csv(DATA / "cell.csv", dtype={"genome_id": str})
    wide = ribo.merge(cell, on="genome_id", validate="one_to_one")
    long = wide.melt(id_vars="genome_id", var_name="antibiotic",
                     value_name="call")
    long["call"] = long["call"].map(CALL)
    if long["call"].isna().any():
        raise ValueError("a call cell is neither 0 nor 1")
    return long


def main() -> int:
    long = build()
    out = DATA / "calls_long.csv"
    long.to_csv(out, index=False)
    print(f"{out}: {len(long)} rows, {long['genome_id'].nunique()} isolates, "
          f"{long['antibiotic'].nunique()} antimicrobials")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
