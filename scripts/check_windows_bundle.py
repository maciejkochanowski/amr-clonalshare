#!/usr/bin/env python3
"""Check a Windows folder built by build_windows_bundle.py, on Windows.

    python scripts/check_windows_bundle.py dist/windows/AMR-ClonalShare-1.0.0-windows

Runs, with the folder's own interpreter, the shipped S. suis configuration
and compares the lineage shares of every call and of every MIC ordering with
the shipped record to 1e-8, imports the local form, and runs its tests with
the folder's interpreter (the test runner is part of the pinned stack the
folder carries). Exits non-zero on any difference.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOLERANCE = 1e-8


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("folder", type=Path)
    parser.add_argument("--out", type=Path, default=ROOT / "dist" / "bundle-check")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--skip-tests", action="store_true")
    args = parser.parse_args(argv)
    python = args.folder / "python" / "python.exe"
    if not python.exists():
        print(f"no interpreter at {python}", file=sys.stderr)
        return 2
    results = args.out / "ssuis"
    if results.exists():
        import shutil
        shutil.rmtree(results)
    print("running the shipped configuration with the folder's interpreter")
    subprocess.run([str(python), "-X", "utf8", "-m", "amr_clonalshare.cli", "--config",
                    str(ROOT / "examples" / "ssuis" / "config.yaml"), "--results-dir", str(results),
                    "--threads", str(args.threads), "--quiet"], check=True, cwd=ROOT)
    got = json.loads((results / "clonal_share_result.json").read_text(encoding="utf-8"))
    want = json.loads((ROOT / "examples" / "ssuis" / "expected" / "clonal_share_result.json").read_text(encoding="utf-8"))
    worst = 0.0
    for agent, block in want["metadata_diagnostics"]["clonal_share"].items():
        other = got["metadata_diagnostics"]["clonal_share"][agent]
        worst = max(worst, abs(other["kappa_adj"] - block["kappa_adj"]))
    for agent, block in want["metadata_diagnostics"]["censored_share"]["per_agent"].items():
        other = got["metadata_diagnostics"]["censored_share"]["per_agent"][agent]
        worst = max(worst, abs(other["order"]["kappa_adj"] - block["order"]["kappa_adj"]))
        worst = max(worst, abs(other["order"]["latent_order_lower"] - block["order"]["latent_order_lower"]))
    print(f"largest difference from the shipped record: {worst:.2e} (tolerance {TOLERANCE:.0e})")
    if worst > TOLERANCE:
        return 1
    print("the local form imports")
    subprocess.run([str(python), "-X", "utf8", "-c", "import amr_clonalshare.gui as g; print(g.__doc__.splitlines()[0])"],
                   check=True)
    print("the command line of the folder")
    subprocess.run([str(args.folder / "amr-clonalshare.cmd"), "version"], check=True)
    if not args.skip_tests:
        # The interpreter of the folder does not see the test runner; the
        # folder keeps it under python\\Lib\\check-helpers for this check.
        # The embeddable interpreter ignores PYTHONPATH, so the path goes
        # onto sys.path in the command itself.
        helpers = args.folder / "python" / "Lib" / "check-helpers"
        print("the tests of the form, with the folder's interpreter")
        subprocess.run([str(python), "-X", "utf8", "-c",
                        "import sys; sys.path.insert(0, sys.argv[1]); import pytest; "
                        "sys.exit(pytest.main(sys.argv[2:]))", str(helpers),
                        "-q", "-p", "no:cacheprovider", str(ROOT / "tests" / "test_gui.py")],
                       check=True, cwd=ROOT)
    print("windows folder OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
