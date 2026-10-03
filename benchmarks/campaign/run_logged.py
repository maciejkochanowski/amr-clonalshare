#!/usr/bin/env python3
"""Run one step of the campaign and write a receipt beside its output.

    python benchmarks/campaign/run_logged.py --receipt OUT/RUN_RECEIPT.json \
        [--outputs OUT] -- python benchmarks/null_uniformity.py --out OUT

The receipt records what a reader needs to tie a result to the code that made
it: the command, the commit of the source tree, a digest of every Python file
of the package, the interpreter and the numerical stack, the Slurm job, the
start, the end and the exit status, and the sha256 of every file the step
wrote under ``--outputs``. The step itself is run unchanged.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def source_digest() -> dict:
    """sha256 of every Python file of the package, and one digest over all."""
    package = ROOT / "src" / "amr_clonalshare"
    files = {p.relative_to(package).as_posix(): _sha256(p) for p in sorted(package.rglob("*.py"))
             if "__pycache__" not in p.parts}
    whole = hashlib.sha256("".join(f"{k} {v}\n" for k, v in files.items()).encode()).hexdigest()
    return {"package_sha256": whole, "files": files}


def commit() -> str | None:
    try:
        return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], check=True,
                              capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def stack() -> dict:
    versions = {"python": platform.python_version(), "platform": platform.platform()}
    for name in ("numpy", "scipy", "pandas", "yaml", "threadpoolctl"):
        try:
            versions[name] = __import__(name).__version__
        except ImportError:
            versions[name] = None
    return versions


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--" not in argv:
        print(__doc__, file=sys.stderr)
        return 2
    split = argv.index("--")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--outputs", type=Path, action="append", default=[])
    args = parser.parse_args(argv[:split])
    command = argv[split + 1:]
    started, clock = _now(), time.perf_counter()
    status = subprocess.run(command, cwd=ROOT).returncode
    outputs = {}
    for directory in args.outputs:
        for path in sorted(directory.rglob("*")):
            if path.is_file() and path.resolve() != args.receipt.resolve():
                outputs[path.relative_to(directory).as_posix()] = _sha256(path)
    receipt = {"command": command, "exit_status": status, "started_utc": started,
               "finished_utc": _now(), "wall_seconds": round(time.perf_counter() - clock, 1),
               "commit": commit(), "source": source_digest(), "stack": stack(),
               "cpus_available": len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count(),
               "slurm": {k: os.environ.get(k) for k in ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID",
                                                        "SLURM_ARRAY_TASK_ID", "SLURM_CPUS_PER_TASK")},
               "outputs_sha256": outputs}
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=1) + "\n", encoding="utf-8")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
