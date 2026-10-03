#!/usr/bin/env python3
"""The support gate of the decomposition, read off its calibration grid.

    python benchmarks/decomposition_gate.py OUT/decomposition_calibration/cells.json OUT/gate.json

For every threshold of shared support in steps of 0.05, the lowest coverage
over the grid cells of ``decomposition_calibration.py`` of each component
among the datasets at or above it, counting cells with at least MIN_DATASETS
such datasets; the gate is the lowest threshold at and above which neither
falls below FLOOR. The coverage of four other bootstrap intervals computed
from the same draws is summarised beside it. This is the ``gate`` step that
wrote ``results_decomposition_calibration/decomposition_gate.json`` (its
receipt, GATE_RECEIPT.json, names the script it then lived in).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

MIN_DATASETS = 100
FLOOR = 0.89


def gate(cells_json: Path, output: Path) -> dict:
    cells = json.loads(cells_json.read_text())
    thresholds = [round(0.05 * i, 2) for i in range(10, 20)]
    out = dict(rule=f"lowest threshold at and above which no cell with at least {MIN_DATASETS} datasets "
                    f"at or above the threshold covers either component below {FLOOR}", thresholds=[])
    for s in thresholds:
        entry = dict(threshold=s)
        for key in ("composition", "within_lineage"):
            worst, judged = 1.0, 0
            for cell in cells:
                bins = [b for b in cell[f"{key}_coverage_by_support"] if b["support_from"] >= s - 1e-9]
                n = sum(b["n"] for b in bins)
                if n >= MIN_DATASETS:
                    judged += 1
                    worst = min(worst, sum(b["covered"] for b in bins) / n)
            entry[f"{key}_lowest_coverage"] = worst
            entry[f"{key}_cells_judged"] = judged
        out["thresholds"].append(entry)
    # The gate is the lowest threshold from which every higher one passes too.
    out["gate"] = None
    for e in reversed(out["thresholds"]):
        if min(e["composition_lowest_coverage"], e["within_lineage_lowest_coverage"]) < FLOOR:
            break
        out["gate"] = e["threshold"]
    methods = ("percentile", "basic", "bias_corrected", "bias_shifted", "studentised")
    out["intervals"] = {f"{key}_{m}": dict(
        lowest=min(c[f"{key}_{m}_coverage"] for c in cells if c[f"{key}_{m}_coverage"] is not None),
        cells_below_0_925=sum(c[f"{key}_{m}_coverage"] < 0.925 for c in cells
                              if c[f"{key}_{m}_coverage"] is not None))
        for key in ("composition", "within_lineage") for m in methods}
    output.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("cells", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    out = gate(args.cells, args.output)
    print(json.dumps({k: out[k] for k in ("gate", "intervals")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
