#!/usr/bin/env python3
"""Key the recorded MIC table to the strain identifiers the layers use.

The cohort's MIC table arrives from BV-BRC with the bare genome identifier
(``1307.7503``) while every other file in this example carries the ``g``
prefix the loader expects (``g1307.7503``). The package deliberately does not
repair identifiers during a join: a loader that stripped or added prefixes
could match the wrong isolate and would have no way to say so. The repair
therefore happens here, once, in the open, and writes a receipt.

Nothing else is changed. Values, units, method and the empty censoring-sign
column pass through as recorded.

    python examples/ssuis/build_mic_table.py
"""
from __future__ import annotations

import json
import pathlib

import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE / "data"
SOURCE = DATA / "mic_long.csv"
TARGET = DATA / "mic_panel.csv"
RECEIPT = DATA / "mic_panel_receipt.json"


def main() -> int:
    mic = pd.read_csv(SOURCE, dtype={"genome_id": str})
    layer = pd.read_csv(DATA / "ribo.csv", dtype={"genome_id": str})
    wanted = set(layer["genome_id"])

    prefixed = "g" + mic["genome_id"].astype(str)
    direct = int(mic["genome_id"].isin(wanted).sum())
    with_prefix = int(prefixed.isin(wanted).sum())
    if with_prefix <= direct:
        raise SystemExit(
            "the prefix is not the difference between the two identifier "
            f"sets (direct {direct}, prefixed {with_prefix}); inspect them "
            "before writing anything")

    out = mic.copy()
    out["genome_id"] = prefixed
    out = out[out["genome_id"].isin(wanted)]
    out.to_csv(TARGET, index=False)

    receipt = {
        "source": SOURCE.name,
        "target": TARGET.name,
        "operation": "prepend 'g' to genome_id, then restrict to the layer index",
        "rows_in": int(len(mic)),
        "rows_out": int(len(out)),
        "strains_in_layers": len(wanted),
        "strains_with_mic": int(out["genome_id"].nunique()),
        "antimicrobials": sorted(out["antibiotic"].astype(str).unique()),
        "readings_per_agent": {
            str(k): int(v) for k, v in
            out.groupby("antibiotic")["measurement"].size().items()},
        "censoring_signs_recorded": int(
            out["measurement_sign"].notna().sum()),
        "units": sorted(out["measurement_unit"].astype(str).unique()),
    }
    RECEIPT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
