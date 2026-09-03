#!/usr/bin/env python3
"""Exercise repository-level contract tests from the extracted 1.0.0 sdist."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import tarfile
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CANDIDATE = ROOT / "release" / "v1.0.0"
#: What an extracted distribution must contain. The manuscript is not in it:
#: the article ships in its own bundle, and the two contract tests that read it
#: skip when it is absent, which is what this run checks.
REQUIRED = (
    "CITATION.cff",
    "LICENSE",
    "Licence.txt",
    "MANIFEST.in",
    "examples/klebsiella/expected/summary.json",
    "examples/klebsiella/expected/cluster_result.json",
    "tests/conftest.py",
    "tests/test_public_release_contract.py",
    "tests/test_manuscript_claims.py",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(command: list[str], *, cwd: Path | None = None) -> dict[str, object]:
    start_utc = utc_now()
    start = time.monotonic()
    proc = subprocess.run(
        command,
        cwd=cwd,
        env={**os.environ, "UV_CACHE_DIR": "/tmp/btc-uv-cache"},
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    record = {
        "command": command,
        "cwd": str(cwd) if cwd else None,
        "start_utc": start_utc,
        "end_utc": utc_now(),
        "elapsed_seconds": round(time.monotonic() - start, 3),
        "exit_code": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
    }
    if proc.returncode != 0:
        raise RuntimeError(json.dumps(record, indent=2))
    return record


def junit_counts(path: Path) -> dict[str, int]:
    root = ElementTree.parse(path).getroot()
    suites = [root] if root.tag.rsplit("}", 1)[-1] == "testsuite" else list(root)
    return {
        key: sum(int(item.attrib.get(key, 0)) for item in suites)
        for key in ("tests", "failures", "errors", "skipped")
    }


def verify(candidate: Path) -> dict[str, object]:
    uv = shutil.which("uv")
    if not uv:
        raise RuntimeError("uv executable not found")
    sdists = list((candidate / "dist").glob("amr_clonalshare-1.0.0.tar.gz"))
    if len(sdists) != 1:
        raise RuntimeError(f"expected one sdist, found {sdists}")
    sdist = sdists[0].resolve()
    temporary_root = Path(tempfile.mkdtemp(prefix="btc-1.0.0-sdist-cleanroom-"))
    extract_root = temporary_root / "extract"
    extract_root.mkdir()
    with tarfile.open(sdist, "r:gz") as archive:
        members = archive.getmembers()
        for member in members:
            destination = (extract_root / member.name).resolve()
            if not destination.is_relative_to(extract_root.resolve()):
                raise RuntimeError(f"unsafe archive member: {member.name}")
        archive.extractall(extract_root)
    roots = [item for item in extract_root.iterdir() if item.is_dir()]
    if len(roots) != 1:
        raise RuntimeError(f"expected one extracted root, found {roots}")
    source = roots[0]
    missing = [relative for relative in REQUIRED if not (source / relative).is_file()]
    if missing:
        raise RuntimeError(f"sdist contract evidence missing: {missing}")

    environment = temporary_root / "venv"
    create = run([uv, "venv", str(environment), "--python", "3.11"])
    python = environment / "bin" / "python"
    install = run([uv, "pip", "install", "--python", str(python), str(sdist), "pytest>=7"])
    junit_path = candidate / "sdist-cleanroom-junit.xml"
    pytest_record = run(
        [
            str(python), "-m", "pytest", "-q",
            "tests/test_public_release_contract.py",
            "tests/test_manuscript_claims.py",
            "--junitxml", str(junit_path.resolve()),
        ],
        cwd=source,
    )
    counts = junit_counts(junit_path)
    passed = counts["tests"] > 0 and counts["failures"] == 0 and counts["errors"] == 0
    return {
        "generated_utc": utc_now(),
        "status": "PASS" if passed else "FAIL",
                "sdist": {"path": str(sdist), "bytes": sdist.stat().st_size, "sha256": sha256(sdist)},
        "temporary_root": str(temporary_root),
        "required_files": {
            relative: sha256(source / relative) for relative in REQUIRED
        },
        "commands": {"create_environment": create, "install": install, "pytest": pytest_record},
        "junit": {**counts, "path": str(junit_path), "sha256": sha256(junit_path)},
        "scope": "Packaging-specific contract tests from an extracted sdist; the full scientific suite is recorded in the campaign receipt.",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    candidate = args.candidate.resolve()
    receipt = args.receipt or candidate / "SDIST_CLEANROOM_RECEIPT.json"
    payload = verify(candidate)
    receipt.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": payload["status"], "junit": payload["junit"]}, indent=2))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
