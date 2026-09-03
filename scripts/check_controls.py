#!/usr/bin/env python3
"""Assert the shipped controls behave as controls.

The positive control must recover the planted truth, the negative control must
report no clusters, and the adversarial control must be flagged as
lineage-confounded. A pipeline that passes only the positive control is not
demonstrating anything.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import json

ROOT = Path(__file__).resolve().parents[1]


def run(cfg: str, out: str) -> dict:
    subprocess.run(["amr-clonalshare", "--config", str(ROOT / cfg),
                    "--results-dir", out, "--quiet"], check=False)
    return json.loads((Path(out) / "summary.json").read_text())


def main() -> int:
    fails = []

    planted = run("examples/synthetic/planted.yaml", "/tmp/ci_planted")
    if planted["selected_k"] != 3:
        fails.append(f"positive control: selected_k = {planted['selected_k']}, expected 3")
    if not planted["structure_detected"]:
        fails.append("positive control: structure not detected")

    null = run("examples/synthetic/null.yaml", "/tmp/ci_null")
    if null["selected_k"] != 1 or not null["no_structure"]:
        fails.append(f"negative control: selected_k = {null['selected_k']}, "
                     f"no_structure = {null['no_structure']}; expected 1 / True")

    clonal = run("examples/synthetic/clonal.yaml", "/tmp/ci_clonal")
    if not clonal.get("lineage_confounded"):
        fails.append("adversarial control: not flagged as lineage-confounded")

    if fails:
        print("control gate FAILED:", *fails, sep="\n  ", file=sys.stderr)
        return 1
    print("control gate OK: planted k=3, null k=1, clonal flagged")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
