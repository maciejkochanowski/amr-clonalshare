#!/usr/bin/env python3
"""One entry point for the 1.0.0 release: check, assemble, build, seal.

    python scripts/prepare_release.py            # everything
    python scripts/prepare_release.py --dry-run  # report what would run

The product version is always `1.0.0`. There is no candidate revision: a
rebuild from the same commit produces the same version, and what identifies a
build is its commit and its manifest rather than a label incremented by hand.

Order matters and each step is a gate for the next.

1. **Preconditions.** The working tree must be clean, the version in
   `pyproject.toml` must be exactly 1.0.0, and it must agree with
   `CITATION.cff` and `.zenodo.json`. A release built from a dirty tree cannot
   be reproduced from its own commit.
2. **Tests and article QA.** The full suite, then `qa_manuscript.py`, which
   checks the claim register hashes against the files they name.
3. **Assemble.** The isolated source and article trees.
4. **Build.** Source distribution and wheel, into the release tree.
5. **Seal.** Receipts, `MANIFEST.tsv` and `SHA256SUMS`.

The three fields that remain author actions - the permanent repository link,
the public documentation URL and the contact email - are not filled by this
script and are reported at the end, because a release that invents them would
be worse than one that says it is not finished.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / "release" / "v1.0.0"
VERSION = "1.0.0"

AUTHOR_ACTIONS = [
    ("C2", "permanent repository link for the v1.0.0 tag"),
    ("C7", "public documentation URL"),
    ("C8", "approved contact and support email"),
]


def run(command: list[str], *, dry: bool, cwd: Path = ROOT) -> int:
    printable = " ".join(str(part) for part in command)
    print(f"  $ {printable}", flush=True)
    if dry:
        return 0
    return subprocess.run(command, cwd=cwd).returncode


def preconditions() -> list[str]:
    problems = []
    # Scoped to this package. The repository holds fourteen projects and the
    # others' work in progress is not this release's business; what has to be
    # committed is what this release is built from.
    # `release/` is this script's own output, so its state is a result rather
    # than an input and does not make the tree dirty.
    dirty = [line for line in subprocess.run(
        ["git", "status", "--porcelain", "--", ".", ":!release"], cwd=ROOT,
        capture_output=True, text=True).stdout.splitlines()]
    if dirty:
        problems.append(f"this package has {len(dirty)} uncommitted changes; "
                        f"commit before releasing. First: {dirty[0][3:]}")

    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    if f'version = "{VERSION}"' not in pyproject:
        problems.append(f"pyproject.toml does not declare version {VERSION}")

    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    if f'version: "{VERSION}"' not in citation:
        problems.append(f"CITATION.cff does not declare version {VERSION}")

    zenodo = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
    if zenodo.get("version") != VERSION:
        problems.append(f".zenodo.json declares {zenodo.get('version')!r}")

    if (ROOT / "LICENSE").read_bytes() != (ROOT / "Licence.txt").read_bytes():
        problems.append("LICENSE and Licence.txt differ")
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-tests", action="store_true",
                        help="only when the suite has just been run")
    args = parser.parse_args()
    python = sys.executable

    print("1. preconditions")
    problems = preconditions()
    for problem in problems:
        print(f"  FAIL {problem}")
    if problems and not args.dry_run:
        return 1
    if not problems:
        print("  all clear")

    print("2. tests and article QA")
    if not args.skip_tests:
        if run([python, "-m", "pytest", "-q"], dry=args.dry_run):
            return 1
    for stage in ("build_supplement.py", "build_claim_register.py",
                  "prose_lint.py"):
        if run([python, f"paper/softwarex/scripts/{stage}"], dry=args.dry_run):
            return 1
    if run([python, "paper/softwarex/scripts/qa_manuscript.py",
            "--receipt", "paper/softwarex/qa/manuscript_qa_receipt.json"],
           dry=args.dry_run):
        return 1

    print("3. assemble")
    if run([python, "scripts/assemble_release_candidate.py",
            "--target", str(RELEASE)], dry=args.dry_run):
        return 1
    if not args.dry_run:
        (RELEASE / ".zenodo.json").write_bytes(
            (ROOT / ".zenodo.json").read_bytes())

    print("4. build the distributions")
    dist = RELEASE / "dist"
    if run([python, "-m", "build", "--sdist", "--wheel",
            "--outdir", str(dist), str(RELEASE / "source")], dry=args.dry_run):
        print("  `python -m build` is required: pip install build")
        return 1

    # `python -m build` writes a .gitignore containing `*` into its output
    # directory. Left there it excludes the very artifacts the release manifest
    # lists, and the failure is silent: SHA256SUMS names files git refuses to
    # store, and nobody finds out until a release is verified from a clone.
    stray = dist / ".gitignore"
    if stray.is_file() and not args.dry_run:
        stray.unlink()
        print(f"  removed {stray.relative_to(ROOT)}, which would have excluded "
              f"the built distributions from version control")

    print("5. verify the built artifacts")
    # An installation matrix and a clean-room run of the packaging contract
    # tests from the extracted sdist. Both write receipts that the seal step
    # requires, so a release cannot be sealed without them.
    for script, label in (("verify_release_matrix.py",
                           "install on Python 3.10 to 3.12"),
                          ("verify_sdist_cleanroom.py",
                           "contract tests from the extracted sdist")):
        print(f"  {label}")
        if run([python, f"scripts/{script}", "--candidate", str(RELEASE)],
               dry=args.dry_run):
            return 1

    print("6. seal")
    if run([python, "scripts/finalize_release_candidate.py",
            "--candidate", str(RELEASE)], dry=args.dry_run):
        return 1

    print("\nremaining author actions, which this script does not invent:")
    for field, description in AUTHOR_ACTIONS:
        print(f"  {field}  {description}")
    print("\nupload order once those exist: GitHub tag v1.0.0, then Zenodo "
          "from the tag, then PyPI from release/v1.0.0/dist/.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
