#!/usr/bin/env python3
"""Hold each module to its own coverage floor.

A single figure over the whole package is easy to satisfy with a well-tested
periphery around an estimator nobody exercises. The floors below are per
module, set from the measured coverage at the time they were written and
raised as the tests grow; a module that carries a reported number has a
floor a reviewer would accept for it, and a module that falls under its floor
fails the build, whatever the package total says.

    pytest --cov=src/amr_clonalshare --cov-report=json:coverage.json
    python scripts/check_coverage_floors.py coverage.json
"""
from __future__ import annotations

import json
import pathlib
import sys

# Percent of statements, per module. A module absent from this table has no
# floor of its own and only counts toward the package total.
FLOORS = {
    # Estimators and the arithmetic that carries a reported number.
    "attribution": 92,
    "censored": 92,
    "clonality": 92,
    "realised": 90,
    "evalues": 93,
    "inference": 84,
    "stats": 84,
    "lineage": 90,
    "archephy": 85,
    "small_n": 85,
    "baselines": 88,
    "influence": 88,
    "fusion": 86,
    "tva": 72,
    "phenotype": 85,
    # Pipeline, input and output.
    "core": 75,
    "cli": 82,
    "config": 84,
    "io": 70,
    "jsonio": 80,
    "qc": 95,
    "schemas": 60,
    "report": 95,
    "report_html": 75,
    "synthetic": 94,
}
PACKAGE_FLOOR = 80


def main(argv: list[str]) -> int:
    path = pathlib.Path(argv[1] if len(argv) > 1 else "coverage.json")
    data = json.loads(path.read_text(encoding="utf-8"))
    measured = {}
    for filename, entry in data["files"].items():
        name = pathlib.Path(filename).name
        if name.endswith(".py"):
            measured[name[:-3]] = float(entry["summary"]["percent_covered"])
    failures = []
    for module, floor in FLOORS.items():
        got = measured.get(module)
        if got is None:
            failures.append("%-14s not measured (floor %d %%)" % (module, floor))
        elif got < floor:
            failures.append("%-14s %5.1f %% is under its floor of %d %%"
                            % (module, got, floor))
    total = float(data["totals"]["percent_covered"])
    if total < PACKAGE_FLOOR:
        failures.append("package        %5.1f %% is under its floor of %d %%"
                        % (total, PACKAGE_FLOOR))
    width = max(len(m) for m in measured)
    for module in sorted(measured):
        floor = FLOORS.get(module)
        print("%-*s %5.1f %%%s" % (width, module, measured[module],
                                   "   floor %d" % floor if floor else ""))
    print("package %.1f %% (floor %d)" % (total, PACKAGE_FLOOR))
    for line in failures:
        print("FAIL", line)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
