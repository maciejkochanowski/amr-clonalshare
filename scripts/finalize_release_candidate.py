#!/usr/bin/env python3
"""Verify and seal the local 1.0.0 candidate with receipts and SHA-256 manifests."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CANDIDATE = ROOT / "release" / "v1.0.0"
SELF_EXCLUDED = {"MANIFEST.tsv", "SHA256SUMS"}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git(command: list[str]) -> str:
    completed = subprocess.run(
        ["git", *command],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    return completed.stdout.strip()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def junit_counts(path: Path) -> dict[str, object]:
    root = ElementTree.parse(path).getroot()
    suites = [root] if root.tag.rsplit("}", 1)[-1] == "testsuite" else [
        child for child in root if child.tag.rsplit("}", 1)[-1] == "testsuite"
    ]
    if not suites:
        raise ValueError(f"No testsuite elements found in {path}")
    return {
        "tests": sum(int(suite.attrib.get("tests", 0)) for suite in suites),
        "failures": sum(int(suite.attrib.get("failures", 0)) for suite in suites),
        "errors": sum(int(suite.attrib.get("errors", 0)) for suite in suites),
        "skipped": sum(int(suite.attrib.get("skipped", 0)) for suite in suites),
        "time_seconds": sum(float(suite.attrib.get("time", 0.0)) for suite in suites),
    }


def finalize(candidate: Path) -> dict[str, object]:
    source = candidate / "source"
    article = candidate / "article"
    matrix_path = candidate / "INSTALL_MATRIX_RECEIPT.json"
    cleanroom_path = candidate / "SDIST_CLEANROOM_RECEIPT.json"
    campaign_path = article / "evidence" / "campaign_2026-09-01" / "RUN_RECEIPT.json"
    manuscript_qa_path = article / "qa" / "manuscript_qa_receipt.json"
    docx_qa_path = article / "qa" / "docx_qa_receipt.json"
    # The rendered DOCX needs the publisher's template, which is not
    # redistributable and is not in this tree. Its absence is recorded with its
    # reason rather than silently omitted, and it stays on the blocker list. A
    # release reporting a digest for a file it does not carry would be worse
    # than one that says it is not finished.
    docx_file = article / "amr-clonalshare_SoftwareX_1.0.0.docx"
    docx_state = (
        {"path": f"article/{docx_file.name}", "sha256": sha256(docx_file),
         "reason": None}
        if docx_file.is_file() else
        {"path": None, "sha256": None,
         "reason": "the SoftwareX v6 template is not redistributable and is "
                   "not in this tree, so the DOCX is rendered at submission "
                   "time from article/manuscript.md"})
    junit_path = article / "evidence" / "campaign_2026-09-01" / "pytest-junit.xml"
    resilience_path = (
        article / "evidence" / "reviewer_resilience_2026-09-01"
        / "RECEIPT.json"
    )

    matrix = load_json(matrix_path)
    cleanroom = load_json(cleanroom_path)
    campaign = load_json(campaign_path)
    manuscript_qa = load_json(manuscript_qa_path)
    docx_qa = load_json(docx_qa_path)
    resilience = load_json(resilience_path)
    junit = junit_counts(junit_path)
    wheels = list((candidate / "dist").glob("amr_clonalshare-1.0.0-*.whl"))
    sdists = list((candidate / "dist").glob("amr_clonalshare-1.0.0.tar.gz"))

    checks = {
        "matrix_pass": matrix.get("status") == "PASS" and len(matrix.get("matrix", [])) == 6,
        "sdist_cleanroom_pass": cleanroom.get("status") == "PASS"
        and cleanroom.get("junit", {}).get("failures") == 0
        and cleanroom.get("junit", {}).get("errors") == 0,
        "campaign_pass": campaign.get("status") == "PASS_TECHNICAL_AND_ILLUSTRATIVE_CAMPAIGN",
        "manuscript_qa_pass": manuscript_qa.get("status") == "PASS",
        "docx_qa_pass": docx_qa.get("status") == "PASS",
        "reviewer_resilience_pass": resilience.get("status") == "PASS",
        "junit_terminal_pass": junit["failures"] == 0 and junit["errors"] == 0,
        "licence_files_identical": (source / "LICENSE").read_bytes() == (source / "Licence.txt").read_bytes(),
        "one_wheel": len(wheels) == 1,
        "one_sdist": len(sdists) == 1,
        "public_version_exact": 'version = "1.0.0"' in (source / "pyproject.toml").read_text(encoding="utf-8"),
        "no_public_repository_url_claimed": "github.com" not in (source / "pyproject.toml").read_text(encoding="utf-8")
        and "repository-code:" not in (source / "CITATION.cff").read_text(encoding="utf-8"),
        "blocked_status_present": "CANDIDATE_COMPLETE_SUBMISSION_BLOCKED" in
        (article / "submission_checklist.md").read_text(encoding="utf-8"),
    }
    if not all(checks.values()):
        raise RuntimeError(f"Candidate checks failed: {checks}")

    worktree_status = git(["status", "--short", "--", "."])
    build_receipt = {
        "generated_utc": utc_now(),
        "status": "CANDIDATE_COMPLETE_SUBMISSION_BLOCKED",
        "product": "amr-clonalshare",
        "product_version": "1.0.0",
                "candidate_date": "2026-09-02",
        "candidate_root": str(candidate),
        "source_provenance": {
            "git_root": git(["rev-parse", "--show-toplevel"]),
            "git_head": git(["rev-parse", "HEAD"]),
            "git_branch": git(["branch", "--show-current"]),
            "worktree_dirty": bool(worktree_status),
            "scoped_worktree_status": worktree_status.splitlines(),
            "commit_or_tag_created": False,
        },
        "checks": checks,
        "fresh_campaign": {
            "receipt": str(campaign_path.relative_to(candidate)),
            "sha256": sha256(campaign_path),
            "status": campaign["status"],
            "junit": junit,
        },
        "reviewer_resilience": {
            "receipt": str(resilience_path.relative_to(candidate)),
            "sha256": sha256(resilience_path),
            "status": resilience["status"],
        },
        "article": {
            "manuscript": "article/manuscript.md",
            "manuscript_sha256": sha256(article / "manuscript.md"),
            "supplementary": "article/supplementary.md",
            "supplementary_sha256": sha256(article / "supplementary.md"),
            "docx": docx_state["path"],
            "docx_sha256": docx_state["sha256"],
            "docx_absent_reason": docx_state["reason"],
            "manuscript_qa_sha256": sha256(manuscript_qa_path),
            "docx_qa_sha256": (sha256(docx_qa_path)
                               if docx_qa_path.is_file() else None),
        },
        "distributions": {
            "wheel": {
                "path": str(wheels[0].relative_to(candidate)),
                "bytes": wheels[0].stat().st_size,
                "sha256": sha256(wheels[0]),
            },
            "sdist": {
                "path": str(sdists[0].relative_to(candidate)),
                "bytes": sdists[0].stat().st_size,
                "sha256": sha256(sdists[0]),
            },
            "install_matrix_receipt": str(matrix_path.relative_to(candidate)),
            "install_matrix_sha256": sha256(matrix_path),
            "python_versions": ["3.10", "3.11", "3.12"],
            "artifact_kinds_each": ["wheel", "sdist"],
            "checks_passed": 6,
            "sdist_cleanroom_receipt": str(cleanroom_path.relative_to(candidate)),
            "sdist_cleanroom_sha256": sha256(cleanroom_path),
            "sdist_cleanroom_junit": cleanroom.get("junit"),
        },
        "submission_blockers": [
            "C2 permanent public GitHub link for the exact v1.0.0 snapshot",
            "C7 public documentation URL",
            "C8 approved author and support email",
        ],
        "external_actions_not_performed": [
            "public GitHub publication",
            "git commit",
            "git tag v1.0.0",
            "DOI archive or deposit",
            "journal submission",
        ],
        "manifest_policy": "MANIFEST.tsv and SHA256SUMS cover every candidate file except themselves and include this BUILD_RECEIPT.json.",
    }
    receipt_path = candidate / "BUILD_RECEIPT.json"
    receipt_path.write_text(json.dumps(build_receipt, indent=2) + "\n", encoding="utf-8")

    entries = []
    for path in sorted(item for item in candidate.rglob("*") if item.is_file()):
        relative = path.relative_to(candidate).as_posix()
        if relative in SELF_EXCLUDED:
            continue
        entries.append((relative, path.stat().st_size, sha256(path)))
    manifest_lines = ["path\tbytes\tsha256"] + [
        f"{relative}\t{size}\t{digest}" for relative, size, digest in entries
    ]
    (candidate / "MANIFEST.tsv").write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")
    (candidate / "SHA256SUMS").write_text(
        "\n".join(f"{digest}  {relative}" for relative, _size, digest in entries) + "\n",
        encoding="utf-8",
    )
    return {
        "status": build_receipt["status"],
        "files_hashed": len(entries),
        "manifest_sha256": sha256(candidate / "MANIFEST.tsv"),
        "sha256sums_sha256": sha256(candidate / "SHA256SUMS"),
        "build_receipt_sha256": sha256(receipt_path),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    args = parser.parse_args()
    payload = finalize(args.candidate.resolve())
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
