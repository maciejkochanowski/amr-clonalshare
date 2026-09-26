#!/usr/bin/env python3
"""Summaries of the model checks, the multi-panel and mixture designs and the
support gate of the decomposition.

    python benchmarks/summarize_model_checks.py shape OUT/shape_check OUT/shape_check.csv \
        --coverage null=OUT/mic_null --coverage extension=OUT/mic_extension \
        --coverage plasmode=OUT/mic_panels
    python benchmarks/summarize_model_checks.py mixture OUT/mic_mixture OUT/mixture.csv
    python benchmarks/summarize_model_checks.py panels OUT/mic_panels OUT/panels.csv
    python benchmarks/summarize_model_checks.py gate OUT/decomposition_calibration/cells.json OUT/gate.json

``shape`` joins every checked dataset to its row of the coverage campaign, by
cell and replicate (by agent and replicate for the S. suis design), and gives
per cell the rejection rate of the check and the coverage of the calibrated test
among the datasets the check passes. ``gate`` gives, for every threshold of
shared support in steps of 0.05, the lowest coverage over the grid cells of
each component among the datasets at or above it, counting cells with at least
MIN_DATASETS such datasets; the gate is the lowest threshold at and above which
neither falls below FLOOR.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd

Z = 1.959963984540054
MIN_DATASETS = 100
FLOOR = 0.89


def wilson(k, n):
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    centre = (p + Z * Z / (2 * n)) / (1 + Z * Z / n)
    half = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / (1 + Z * Z / n)
    return centre - half, centre + half


def _json_rows(directory: Path) -> pd.DataFrame:
    rows = []
    for path in sorted(directory.glob("*.json")):
        if not path.name.endswith("RECEIPT.json"):
            rows.extend(json.loads(path.read_text()))
    if not rows:
        raise SystemExit(f"no results in {directory}")
    return pd.DataFrame(rows)


def _rate(frame, column):
    k, n = int(frame[column].astype(bool).sum()), len(frame)
    low, high = wilson(k, n)
    return dict(n=n, rate=k / n if n else float("nan"), wilson_low=low, wilson_high=high)


def shape(directory: Path, output: Path, coverage: list[str]) -> None:
    data = _json_rows(directory)
    key = {"null": ["cell", "replicate"], "extension": ["cell", "replicate"], "plasmode": ["agent", "replicate"]}
    joined = []
    for source, frame in data.groupby("source"):
        frame = frame.copy()
        paths = dict(c.split("=", 1) for c in coverage)
        if source in paths:
            where = Path(paths[source])
            if source == "plasmode":
                cov = _json_rows(where)
                cov = cov[cov.arm == "plasmode"][["agent", "replicate", "covered"]]
            else:
                cov = pd.concat([pd.read_csv(p) for p in sorted(where.glob("cell_*.csv"))], ignore_index=True)
                cov = cov[["cell", "replicate", "interval_test_covered"]].rename(
                    columns={"interval_test_covered": "covered"})
            frame = frame.merge(cov, on=key[source], how="left", validate="one_to_one")
        joined.append(frame)
    data = pd.concat(joined, ignore_index=True)
    rows = []
    group = [c for c in ("source", "cell", "agent", "domain", "residual", "reading", "rho") if c in data]
    for labels, frame in data.groupby(group, dropna=False):
        checked = frame[frame.checked.astype(bool)]
        row = dict(zip(group, labels if isinstance(labels, tuple) else (labels,)))
        row.update({f"reject_{k}": v for k, v in _rate(checked, "rejected").items()})
        row["refused"] = int((~frame.checked.astype(bool)).sum())
        if "covered" in frame and frame.covered.notna().any():
            for name, part in (("all", checked), ("passed", checked[~checked.rejected.astype(bool)]),
                               ("rejected", checked[checked.rejected.astype(bool)])):
                part = part[part.covered.notna()]
                row.update({f"coverage_{name}_{k}": v for k, v in _rate(part, "covered").items()})
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(output, index=False)
    print(table.groupby(["source", "residual"] if "residual" in table else ["source"], dropna=False)
          .agg(cells=("reject_n", "size"), mean_rejection=("reject_rate", "mean"),
               max_rejection=("reject_rate", "max")).to_string())


def mixture(directory: Path, output: Path) -> None:
    data = _json_rows(directory)
    rows = []
    for cell, frame in data.groupby("cell"):
        first = frame.iloc[0]
        row = dict(cell=int(cell), mu_wt=first.mu_wt, pi=first.pi, rho_z=first.rho_z, rho_np=first.rho_np,
                   misclassified=float(frame.misclassified.mean()),
                   gaussian_rho_hat_mean=float(frame.g_rho_hat.mean()),
                   population_rho_hat_mean=float(frame.p_rho_hat.mean()))
        row.update({f"gaussian_coverage_{k}": v for k, v in _rate(frame, "g_covered").items()})
        checked = frame[frame.shape_checked.astype(bool)]
        row.update({f"shape_reject_{k}": v for k, v in _rate(checked, "shape_rejected").items()})
        row.update({f"population_coverage_{k}": v for k, v in _rate(frame, "p_covered").items()})
        mixed = frame[frame.p_mixing_rejected.notna()]
        row.update({f"mixing_reject_{k}": v for k, v in _rate(mixed, "p_mixing_rejected").items()})
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(output, index=False)
    print(table[["cell", "rho_np", "gaussian_coverage_rate", "shape_reject_rate", "rho_z",
                 "population_coverage_rate", "mixing_reject_rate"]].to_string())


def panels(directory: Path, output: Path) -> None:
    data = _json_rows(directory)
    group = [c for c in ("mode", "cell", "agent", "rho", "alignment", "shift", "arm") if c in data]
    rows = []
    for labels, frame in data.groupby(group, dropna=False):
        row = dict(zip(group, labels))
        row.update({f"coverage_{k}": v for k, v in _rate(frame, "covered").items()})
        row["reportable"] = float(frame.reportable.astype(bool).mean())
        if "low" in frame and frame.low.notna().any():
            width = (frame.high - frame.low).dropna()
            row.update(intervals=int(width.size), median_width=float(width.median()))
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(output, index=False)
    print(table.to_string())


def gate(cells_json: Path, output: Path) -> None:
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
    print(json.dumps({k: out[k] for k in ("gate", "intervals")}, indent=1))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("step", choices=["shape", "mixture", "panels", "gate"])
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--coverage", action="append", default=[],
                        help="SOURCE=DIRECTORY of the coverage campaign that ran the same datasets")
    args = parser.parse_args()
    if args.step == "shape":
        shape(args.source, args.output, args.coverage)
    else:
        {"mixture": mixture, "panels": panels, "gate": gate}[args.step](args.source, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
