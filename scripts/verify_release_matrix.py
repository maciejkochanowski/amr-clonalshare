#!/usr/bin/env python3
"""Install the 1.0.0 sdist and wheel in clean Python 3.10-3.12 environments."""

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
import zipfile
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CANDIDATE = ROOT / "release" / "v1.0.0"
PYTHONS = ("3.10", "3.11", "3.12")
ARTIFACT_KINDS = ("wheel", "sdist")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(command: list[str], *, cwd: Path | None = None) -> dict[str, object]:
    started = utc_now()
    start = time.monotonic()
    completed = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env={**os.environ, "UV_CACHE_DIR": "/tmp/btc-uv-cache"},
        check=False,
    )
    record = {
        "command": command,
        "cwd": str(cwd) if cwd else None,
        "start_utc": started,
        "end_utc": utc_now(),
        "elapsed_seconds": round(time.monotonic() - start, 3),
        "exit_code": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }
    if completed.returncode != 0:
        raise RuntimeError(json.dumps(record, indent=2))
    return record


def archive_audit(wheel: Path, sdist: Path) -> dict[str, object]:
    with zipfile.ZipFile(wheel) as archive:
        wheel_files = sorted(archive.namelist())
    with tarfile.open(sdist, "r:gz") as archive:
        sdist_files = sorted(member.name for member in archive.getmembers() if member.isfile())

    wheel_licenses = [name for name in wheel_files if "/licenses/" in name]
    sdist_licenses = [name for name in sdist_files if name.endswith(("/LICENSE", "/Licence.txt"))]
    # The distribution carries no part of the article, so the manuscript is not
    # among the evidence a distribution must contain. The two contract tests
    # that read it skip when it is absent, which the cleanroom run exercises.
    required_sdist_evidence = [
        "CITATION.cff",
        "MANIFEST.in",
        "examples/klebsiella/expected/summary.json",
        "examples/klebsiella/expected/cluster_result.json",
        "tests/conftest.py",
    ]
    required_wheel_modules = [
        "amr_clonalshare/__init__.py",
        "amr_clonalshare/cli.py",
        "amr_clonalshare/core.py",
    ]
    return {
        "wheel": {
            "path": str(wheel),
            "sha256": sha256(wheel),
            "bytes": wheel.stat().st_size,
            "file_count": len(wheel_files),
            "licenses": wheel_licenses,
            "required_modules_present": {
                module: any(name.endswith(module) for name in wheel_files)
                for module in required_wheel_modules
            },
        },
        "sdist": {
            "path": str(sdist),
            "sha256": sha256(sdist),
            "bytes": sdist.stat().st_size,
            "file_count": len(sdist_files),
            "licenses": sdist_licenses,
            "pyproject_present": any(name.endswith("/pyproject.toml") for name in sdist_files),
            "readme_present": any(name.endswith("/README.md") for name in sdist_files),
            "contract_evidence_present": {
                item: any(name.endswith("/" + item) for name in sdist_files)
                for item in required_sdist_evidence
            },
        },
        "passed": (
            any(name.endswith("/LICENSE") for name in wheel_licenses)
            and any(name.endswith("/Licence.txt") for name in wheel_licenses)
            and any(name.endswith("/LICENSE") for name in sdist_licenses)
            and any(name.endswith("/Licence.txt") for name in sdist_licenses)
            and all(any(name.endswith(module) for name in wheel_files) for module in required_wheel_modules)
            and any(name.endswith("/pyproject.toml") for name in sdist_files)
            and any(name.endswith("/README.md") for name in sdist_files)
            and all(
                any(name.endswith("/" + item) for name in sdist_files)
                for item in required_sdist_evidence
            )
        ),
    }


SMOKE_PROGRAM = r"""
import importlib.metadata as md
import importlib.util
import inspect
import json
import platform
import sys
import amr_clonalshare
from amr_clonalshare import core

dependencies = {}
for dist_name in ["numpy", "pandas", "scipy", "scikit-learn", "PyYAML", "pandera"]:
    dependencies[dist_name] = md.version(dist_name)
required_modules = [
    "amr_clonalshare.cli",
    "amr_clonalshare.config",
    "amr_clonalshare.core",
    "amr_clonalshare.fusion",
    "amr_clonalshare.inference",
    "amr_clonalshare.jsonio",
    "amr_clonalshare.phenotype",
    "amr_clonalshare.stats",
]
payload = {
    "python": platform.python_version(),
    "executable": sys.executable,
    "distribution_version": md.version("amr-clonalshare"),
    "module_version": amr_clonalshare.__version__,
    "schema_version": "1.0" if '\"schema_version\": \"1.0\"' in inspect.getsource(core.run) else None,
    "dependencies": dependencies,
    "modules_present": {name: importlib.util.find_spec(name) is not None for name in required_modules},
    "installed_distribution_files": [str(path) for path in (md.files("amr-clonalshare") or [])],
}
assert payload["distribution_version"] == "1.0.0"
assert payload["module_version"] == "1.0.0"
assert payload["schema_version"] == "1.0"
assert all(payload["modules_present"].values())
print(json.dumps(payload, sort_keys=True))
"""


def verify(candidate: Path) -> dict[str, object]:
    uv = shutil.which("uv")
    if not uv:
        raise RuntimeError("uv executable not found")
    dist = candidate / "dist"
    wheels = list(dist.glob("amr_clonalshare-1.0.0-*.whl"))
    sdists = list(dist.glob("amr_clonalshare-1.0.0.tar.gz"))
    if len(wheels) != 1 or len(sdists) != 1:
        raise RuntimeError(f"Expected one wheel and one sdist, found {wheels} and {sdists}")
    artifacts = {"wheel": wheels[0].resolve(), "sdist": sdists[0].resolve()}
    archive = archive_audit(artifacts["wheel"], artifacts["sdist"])
    if not archive["passed"]:
        raise RuntimeError(f"Archive audit failed: {archive}")

    records: list[dict[str, object]] = []
    for python_version in PYTHONS:
        for kind in ARTIFACT_KINDS:
            temp_root = Path(
                tempfile.mkdtemp(prefix=f"btc-1.0.0-py{python_version.replace('.', '')}-{kind}-")
            )
            environment = temp_root / "venv"
            create = run([uv, "venv", str(environment), "--python", python_version])
            python = environment / "bin" / "python"
            install = run([uv, "pip", "install", "--python", str(python), str(artifacts[kind])])
            smoke = run([str(python), "-c", SMOKE_PROGRAM])
            smoke_payload = json.loads(str(smoke["stdout"]).strip())
            cli_version = run([str(environment / "bin" / "amr-clonalshare"), "--version"])
            cli_help = run([str(environment / "bin" / "amr-clonalshare"), "--help"])
            expected_prefix = python_version + "."
            passed = (
                str(smoke_payload["python"]).startswith(expected_prefix)
                and str(cli_version["stdout"]).strip() == "amr-clonalshare 1.0.0"
                and "--config" in str(cli_help["stdout"])
            )
            records.append(
                {
                    "python_requested": python_version,
                    "artifact_kind": kind,
                    "artifact_path": str(artifacts[kind]),
                    "artifact_sha256": sha256(artifacts[kind]),
                    "temporary_environment": str(environment),
                    "create_environment": create,
                    "install": install,
                    "smoke": smoke_payload,
                    "cli_version": str(cli_version["stdout"]).strip(),
                    "cli_help_contains_config": "--config" in str(cli_help["stdout"]),
                    "passed": passed,
                }
            )
    return {
        "generated_utc": utc_now(),
        "candidate": str(candidate),
        "product_version": "1.0.0",
                "archive_audit": archive,
        "matrix": records,
        "status": "PASS" if all(record["passed"] for record in records) else "FAIL",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    candidate = args.candidate.resolve()
    receipt = args.receipt or candidate / "INSTALL_MATRIX_RECEIPT.json"
    payload = verify(candidate)
    receipt.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": payload["status"], "receipt": str(receipt), "checks": len(payload["matrix"])}, indent=2))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
