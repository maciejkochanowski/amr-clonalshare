"""Rewrite the golden artefacts under tests/golden from the current code.

The golden tests compare whole artefacts, so an intended change of wording, of
a key or of a reported number shows up here as a diff of a committed file.
Run this after such a change, read the diff, and commit it with the change.
It writes nothing outside tests/golden.
"""
from __future__ import annotations

import importlib.util
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _load_case_module():
    path = ROOT / "tests" / "test_golden_artefacts.py"
    spec = importlib.util.spec_from_file_location("golden_cases", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main() -> int:
    sys.path.insert(0, str(ROOT / "src"))
    cases = _load_case_module()
    written = 0
    for case in cases.CASES:
        target = cases.GOLDEN / case
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out"
            cases.run_case(case, out)
            target.mkdir(parents=True, exist_ok=True)
            for name in cases.case_files(case):
                shutil.copyfile(out / name, target / name)
                written += 1
        print(f"{case}: {len(cases.case_files(case))} artefacts")
    print(f"wrote {written} files under {cases.GOLDEN}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
