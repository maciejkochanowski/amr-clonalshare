#!/usr/bin/env python3
"""Does the majority-class baseline have anywhere to move in these grids.

    python leverage_check.py --atlas <vet_atlas.json> --out <dir>

WHY THIS RUNS FIRST. The comparison Arm B exists to make is between a
classification accuracy and a variance share, and the whole force of it comes
from the majority-class baseline: at prevalence 0.90 an accuracy of 0.90 is
what predicting the majority for everyone already achieves, so an accuracy
number carries information about lineage only to the extent that it exceeds
that baseline. If every cell of a grid sat near prevalence 0.5 the baseline
would be a constant 0.5, accuracy and accuracy-minus-baseline would be the
same number shifted, and every correlation computed from them would agree by
construction rather than by measurement. A stable answer from such a grid
would be an artefact of the design. The prereg therefore declares a grid void
for this comparison when more than nine tenths of its cells carry a baseline
below 0.60, or when fewer than one twentieth carry a baseline of 0.70 or more,
and this script computes that verdict before any accuracy is estimated.

The real cells are read back from the published atlas rather than recomputed,
because the prevalence that matters is the one the published clonal share was
estimated at. The synthetic prevalences are generated from the benchmark's own
cohort() function, since for a binary cohort the marginal prevalence is not the
target prevalence: a lineage effect on the probit scale pulls the marginal
towards one half as the between-lineage variance grows, so the realised spread
has to be measured rather than assumed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from estimator_benchmark import SEED, cohort, design

SWEEP_PREVALENCES = (0.10, 0.25, 0.50, 0.75, 0.90)
#: The addendum grid. The parent registration set a floor of thirty cells for a
#: rank correlation and then registered a sweep of twenty five, so the sweep
#: fell below the parent's own floor. Nine prevalence levels is the smallest
#: change that clears it while keeping the axis and the geometry unchanged.
SWEEP_PREVALENCES_EXTENDED = (0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80,
                              0.90)
SWEEP_SHARES = (0.0, 0.1, 0.3, 0.5, 0.7)
SWEEP_GROUPS = 30
SWEEP_SIZE = 20
SWEEP_REPLICATES = 400


def sweep_design(prevalences=SWEEP_PREVALENCES) -> list[dict]:
    """The registered prevalence extension, in a fixed order.

    The shipped binary grid is generated at one target prevalence, so it cannot
    move the majority baseline and cannot by itself decide whether accuracy and
    the variance share order cells differently. This grid crosses prevalence
    with the true share at one cohort geometry so that the baseline sweeps from
    0.90 down to 0.50 and back up, which is the axis the objection lives on.
    """
    cells, i = [], 0
    for prevalence in prevalences:
        for share in SWEEP_SHARES:
            cells.append({"index": i, "binary": True,
                          "n_groups": SWEEP_GROUPS, "group_size": SWEEP_SIZE,
                          "share": share, "prevalence": prevalence,
                          "unbalance": "balanced",
                          "replicates": SWEEP_REPLICATES})
            i += 1
    return cells


def _summary(baselines: np.ndarray, prevalences: np.ndarray) -> dict:
    """The leverage verdict for one grid, on the baseline scale."""
    b = np.asarray(baselines, dtype=float)
    p = np.asarray(prevalences, dtype=float)
    n = int(b.size)
    return {
        "n_cells": n,
        "prevalence_min": float(p.min()), "prevalence_max": float(p.max()),
        "prevalence_quartiles": [float(x) for x in
                                 np.percentile(p, [25, 50, 75])],
        "baseline_min": float(b.min()), "baseline_max": float(b.max()),
        "baseline_range": float(b.max() - b.min()),
        "baseline_quartiles": [float(x) for x in
                               np.percentile(b, [25, 50, 75])],
        "share_baseline_at_least_0.70": float((b >= 0.70).mean()),
        "share_baseline_at_least_0.80": float((b >= 0.80).mean()),
        "share_baseline_at_least_0.85": float((b >= 0.85).mean()),
        "share_baseline_below_0.60": float((b < 0.60).mean()),
        "void_for_this_comparison": bool((b < 0.60).mean() > 0.90
                                         or (b >= 0.70).mean() < 0.05),
    }


def real_cells(atlas: Path) -> dict:
    """Prevalence of every published veterinary cell that carries a share."""
    payload = json.loads(atlas.read_text())
    rows = []
    for cohort_record in payload["cohorts"]:
        for agent, record in cohort_record["agents"].items():
            if "kappa" not in record:
                continue
            p = float(record["prevalence"])
            rows.append({"organism": cohort_record["organism"],
                         "host_group": cohort_record["host_group"],
                         "matrix": cohort_record["matrix"], "agent": agent,
                         "n": int(record["n"]), "prevalence": p,
                         "baseline_marginal": float(max(p, 1.0 - p)),
                         "kappa": float(record["kappa"]),
                         "estimable": bool(record["estimable"])})
    p = np.array([r["prevalence"] for r in rows])
    b = np.array([r["baseline_marginal"] for r in rows])
    return {"cells": rows, "leverage": _summary(b, p)}


def synthetic_cells(cells: list[dict]) -> dict:
    """Realised marginal prevalence of every synthetic cell, by generation.

    Only the generator runs here. No estimator is called, so this stays a
    description of the design and not a comparator number.
    """
    rows = []
    for cell in cells:
        rng = np.random.default_rng(SEED + 977 * cell["index"])
        kwargs = {"binary": True, "unbalance": cell["unbalance"]}
        if "prevalence" in cell:
            kwargs["prevalence"] = cell["prevalence"]
        draws = []
        for _ in range(cell["replicates"]):
            y, _lineage, _truth = cohort(rng, cell["n_groups"],
                                         cell["group_size"], cell["share"],
                                         **kwargs)
            draws.append(float(y.mean()))
        draws = np.asarray(draws)
        mean_p = float(draws.mean())
        rows.append(dict(cell, realised_prevalence=mean_p,
                         realised_prevalence_sd=float(draws.std(ddof=1)),
                         realised_prevalence_mc_se=float(
                             draws.std(ddof=1) / np.sqrt(draws.size)),
                         baseline_marginal=float(max(mean_p, 1.0 - mean_p))))
    p = np.array([r["realised_prevalence"] for r in rows])
    b = np.array([r["baseline_marginal"] for r in rows])
    return {"cells": rows, "leverage": _summary(b, p)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--atlas", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    shipped_binary = [c for c in design() if c["binary"]]
    payload = {
        "prereg_sha256": ("4d54749cf39ed17db74ad997fec31dd617"
                          "38a5500c77961b8809627780438df9"),
        "R1_real": real_cells(args.atlas),
        "S1_shipped_binary": synthetic_cells(shipped_binary),
        "S2_prevalence_sweep": synthetic_cells(sweep_design()),
    }
    (args.out / "leverage_check.json").write_text(
        json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    for name in ("R1_real", "S1_shipped_binary", "S2_prevalence_sweep"):
        print(name, json.dumps(payload[name]["leverage"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
