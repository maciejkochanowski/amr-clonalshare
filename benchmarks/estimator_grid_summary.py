#!/usr/bin/env python3
"""Write the package's ``validation_grid.json`` from an estimator benchmark.

    python benchmarks/estimator_grid_summary.py <dir>/estimator_benchmark.json \
        src/amr_clonalshare/validation_grid.json

Every report quotes this file: the coverage each interval reached on the
grid of ``estimator_benchmark.py``. Nothing here is computed beyond counting
the cells of that file.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

#: The estimators the report quotes, with the file's name for each.
QUOTED = ("clonal_share", "clonal_share_superpopulation", "realised_share")


def _block(cells, name):
    values = [c[name]["coverage"] for c in cells if c[name]["coverage"] is not None]
    return {"mean": round(sum(values) / len(values), 3), "cells": len(values),
            "runs": sum(c["replicates"] for c in cells if c[name]["coverage"] is not None),
            "min_cell": round(min(values), 3), "max_cell": round(max(values), 3)}


def summarise(payload: dict) -> dict:
    cells = payload["cells"]
    provenance = payload["provenance"]
    coverage = {}
    for name in QUOTED:
        coverage[name] = {kind: _block([c for c in cells if c["binary"] == (kind == "binary")], name)
                          for kind in ("binary", "continuous")}
    few = [c["clonal_share"]["coverage"] for c in cells
           if c["n_groups"] == 5 and c["share"] > 0 and c["clonal_share"]["coverage"] is not None]
    more = [c["clonal_share"]["coverage"] for c in cells
            if c["n_groups"] >= 10 and c["share"] > 0 and c["clonal_share"]["coverage"] is not None]
    coverage["clonal_share"]["few_lineages"] = {
        "lineages": 5, "min_cell": round(min(few), 3), "max_cell": round(max(few), 3),
        "what": ("coverage of the lineage-membership share by the lineage-bootstrap interval "
                 "in the five-lineage cells with a non-zero share; from ten lineages up the "
                 f"cells covered at {min(more):.2f} to {max(more):.2f}")}
    groups = sorted({c["n_groups"] for c in cells})
    sizes = sorted({c["group_size"] for c in cells})
    prevalences = sorted({c.get("prevalence", 0.25) for c in cells if c["binary"]})
    return {
        "what": ("Interval coverage measured on the release's validation grid: simulated cohorts "
                 "with a known share, one estimator interval per run, the fraction of runs whose "
                 "interval contained the truth. These values belong to the release, not to any run "
                 "that cites them. Each estimator is scored against the quantity it estimates: "
                 "clonal_share against the lineage-membership share of this collection, "
                 "clonal_share_superpopulation against the share a fresh draw of lineages would "
                 "show, realised_share against the design-corrected component ratio the exact "
                 f"interval is derived for. The grid spans {groups[0]} to {groups[-1]:,} lineages, "
                 f"{sizes[0]} to {sizes[-1]} isolates a lineage, normal and two-point lineage laws, "
                 "and binary prevalences of "
                 + " and ".join(f"{round(100 * p)}" for p in prevalences) + " per cent."),
        "source": (f"estimator_benchmark.json written by benchmarks/estimator_benchmark.py, "
                   f"seed {provenance['seed']}, nominal level {1 - provenance['alpha']:.2f}"),
        "coverage": coverage,
    }


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    payload = json.loads(Path(argv[0]).read_text(encoding="utf-8"))
    Path(argv[1]).write_text(json.dumps(summarise(payload), indent=1) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
