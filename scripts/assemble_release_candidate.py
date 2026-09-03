#!/usr/bin/env python3
"""Assemble the immutable input tree for the SoftwareX release 1.0.0."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TARGET = ROOT / "release" / "v1.0.0"

SOURCE_FILES = [
    "README.md",
    "NOVELTY_EVIDENCE.md",
    "CHANGELOG.md",
    "CITATION.cff",
    "CODE_OF_CONDUCT.md",
    "CONTRIBUTING.md",
    "LICENSE",
    "Licence.txt",
    "MANIFEST.in",
    "pyproject.toml",
    "uv.lock",
]
SOURCE_DIRS = [
    "src",
    "tests",
    "docs",
    "examples",
    "benchmarks",
    "audit_checks",
    "scripts",
]
ARTICLE_FILES = [
    "README.md",
    "manuscript.md",
    "supplementary.md",
    "highlights.txt",
    "cover_letter.md",
    "submission_checklist.md",
    "claim-evidence.tsv",
]
#: Copied when present. The rendered DOCX needs the publisher's template, which
#: is not redistributable, so a bundle built without it is complete apart from
#: the render and says so rather than failing.
OPTIONAL_ARTICLE_FILES = [
    "amr-clonalshare_SoftwareX_1.0.0.docx",
]
ARTICLE_QA_FILES = [
    "a11y_audit.json",
    "style_lint.json",
    "manuscript_qa_receipt.json",
    "docx_qa_receipt.json",
    "template_distillation.md",
]


#: Run scratch, not source. `audit_checks/_tmp/` and `_permuted/` hold YAML
#: configurations written by an audit script against the tree as it stood that
#: day, and they carry the absolute paths of the machine that wrote them. The
#: tree-wide `tools/check_paths.py` lists both as workbenches that are not
#: distributed; copying them here made that statement false and put twenty
#: `/home/<someone>/` paths into the release.
SCRATCH_DIRECTORIES = {"__pycache__", ".pytest_cache", ".hypothesis",
                       ".DS_Store", "_tmp", "_permuted"}


def ignore_generated(_directory: str, names: list[str]) -> set[str]:
    return {
        name
        for name in names
        if name in SCRATCH_DIRECTORIES or name.endswith((".pyc", ".pyo"))
    }


def copy_file(source: Path, target: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def copy_directory(source: Path, target: Path) -> None:
    if not source.is_dir():
        raise FileNotFoundError(source)
    shutil.copytree(source, target, ignore=ignore_generated)


def assemble(target: Path) -> None:
    if target.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing candidate directory: {target}"
        )
    source_root = target / "source"
    article_root = target / "article"
    target.mkdir(parents=True)

    for relative in SOURCE_FILES:
        copy_file(ROOT / relative, source_root / relative)
    for relative in SOURCE_DIRS:
        copy_directory(ROOT / relative, source_root / relative)
    # Nothing under paper/ is copied here. The software bundle and the article
    # bundle are separate products, and a manuscript sitting inside an
    # installation invites the manuscript to be read as a description of what
    # was installed. The two contract tests that read the retired methods
    # manuscript skip when it is absent and run in the article bundle, where
    # `qa_manuscript.py` pins the same digest.
    copy_file(ROOT / ".zenodo.json", source_root / ".zenodo.json")
    copy_file(
        ROOT / ".github" / "workflows" / "ci.yml",
        source_root / ".github" / "workflows" / "ci.yml",
    )

    softwarex = ROOT / "paper" / "softwarex"
    for relative in ARTICLE_FILES:
        copy_file(softwarex / relative, article_root / relative)
    for relative in OPTIONAL_ARTICLE_FILES:
        if (softwarex / relative).is_file():
            copy_file(softwarex / relative, article_root / relative)
    copy_file(ROOT / "paper" / "README.md", article_root / "PAPER_README.md")
    # The retired methods manuscript belongs with the article, not with the
    # software, and its digest is pinned by the article QA.
    copy_file(ROOT / "paper" / "manuscript.md",
              article_root / "retired_methods_manuscript.md")
    copy_directory(softwarex / "figures", article_root / "figures")
    copy_directory(
        softwarex / "evidence" / "campaign_2026-09-01",
        article_root / "evidence" / "campaign_2026-09-01",
    )
    copy_directory(
        softwarex / "evidence" / "reviewer_resilience_2026-09-01",
        article_root / "evidence" / "reviewer_resilience_2026-09-01",
    )
    censored = softwarex / "evidence" / "censored_2026-09-01"
    if censored.is_dir():
        copy_directory(censored, article_root / "evidence" / "censored_2026-09-01")
    # The cross-species atlas, the realised-share campaign and the
    # pre-registration written before the run it scores. Each directory is
    # named for the day it was produced, so a reader can tell which
    # evidence predates which.
    for name in ("atlas_2026-09-02", "realised_2026-09-02", "prereg",
                 "vet_atlas_2026-09-02", "benchmark_2026-09-02"):
        directory = softwarex / "evidence" / name
        if directory.is_dir():
            copy_directory(directory, article_root / "evidence" / name)
    copy_directory(softwarex / "review", article_root / "review")
    for relative in ARTICLE_QA_FILES:
        if (softwarex / "qa" / relative).is_file():
            copy_file(softwarex / "qa" / relative, article_root / "qa" / relative)
    if (softwarex / "qa" / "docx_render_verified").is_dir():
        copy_directory(
            softwarex / "qa" / "docx_render_verified",
            article_root / "qa" / "docx_render_verified",
        )
    copy_directory(softwarex / "scripts", article_root / "scripts")

    release_readme = target / "RELEASE_README.md"
    release_readme.write_text(
        "# amr-clonalshare 1.0.0\n\n"
        "Status: **CANDIDATE_COMPLETE_SUBMISSION_BLOCKED**\n\n"
        "Three bundles, each complete on its own and none containing another.\n\n"
        "- `source/` is the repository bundle: the package, its tests, the "
        "examples, the documentation, the benchmarks and the release tooling. "
        "It contains no part of the article.\n"
        "- `dist/` is the Python distribution bundle, the sdist and the wheel "
        "built from `source/`. It contains no part of the article either, and "
        "the two contract tests that read the retired methods manuscript skip "
        "when it is absent.\n"
        "- `article/` is the article bundle: the SoftwareX manuscript in "
        "Markdown and as the rendered DOCX, its supplement, figures, cover "
        "letter, checklist, claim register, the campaign, calibration and "
        "cross-species evidence, the QA receipts including the page-by-page "
        "render that was inspected, and the retired methods manuscript it "
        "supersedes.\n\n"
        "The product version is exactly `1.0.0`; there is no candidate "
        "revision, and a rebuild from the same commit produces the same "
        "version.\n\n"
        "Three fields remain author actions and are the reason the status is "
        "blocked rather than released: the permanent repository link, the "
        "public documentation URL and the approved contact email. No tag, DOI, "
        "public upload or journal submission is claimed here.\n\n"
        "`dist/`, the Python installation-matrix receipt, `BUILD_RECEIPT.json`, "
        "`MANIFEST.tsv` and `SHA256SUMS` are added by the subsequent verified "
        "build.\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
    args = parser.parse_args()
    assemble(args.target.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
