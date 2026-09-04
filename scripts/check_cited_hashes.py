#!/usr/bin/env python3
"""Check that every sha256 quoted in a document still matches the file it names.

A document that quotes a digest is making a checkable claim. When the file is
regenerated and the digest is not, the claim becomes false silently, and the
reader has no way to tell. This script reads every markdown file given to it,
finds each `path/to/file` ... sha256 <prefix> citation, resolves the path
against the package root, and compares.

Exit status is 1 if any citation is wrong, so it can gate a release.

    python scripts/check_cited_hashes.py [document ...]
"""
from __future__ import annotations

import hashlib
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
CITATION = re.compile(
    r"`(?P<path>[A-Za-z0-9_./\-]+\.(?:json|tsv|csv|md))`"
    r"[\s\S]{0,90}?sha256\s+(?P<digest>[0-9a-f]{8,64})"
)


def digest(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 16), b""):
            h.update(block)
    return h.hexdigest()


def check(document: pathlib.Path) -> tuple[int, int, int]:
    text = document.read_text(encoding="utf-8")
    ok = wrong = unresolved = 0
    for match in CITATION.finditer(text):
        # A document may name a file relative to the package root or
        # relative to itself; a citation that resolves either way is a claim
        # and gets checked. Anything else is prose that happens to look like a
        # path, and is counted rather than failed.
        target = ROOT / match.group("path")
        if not target.is_file():
            target = document.parent / match.group("path")
        if not target.is_file():
            unresolved += 1
            continue
        actual = digest(target)
        if actual.startswith(match.group("digest")):
            ok += 1
        else:
            wrong += 1
            print(
                "%s: %s cited as %s, actual %s"
                % (document.name, match.group("path"),
                   match.group("digest"), actual[: len(match.group("digest"))])
            )
    return ok, wrong, unresolved


def main(argv: list[str]) -> int:
    documents = [pathlib.Path(a) for a in argv[1:]]
    if not documents:
        # The digests are quoted in the article material as well as in the
        # package documentation, and a gate that reads only the latter checks
        # nothing. Superseded drafts are left out: they cite the files as they
        # stood when they were written, which is what a superseded draft is
        # for.
        documents = (sorted(ROOT.glob("*.md"))
                     + sorted(ROOT.glob("docs/*.md"))
                     + sorted(ROOT.glob("paper/*.md"))
                     + sorted(p for p in ROOT.glob("paper/softwarex/*.md")
                              if "_prev" not in p.name))
    ok = wrong = unresolved = 0
    for document in documents:
        if not document.is_file():
            print("no such document: %s" % document)
            return 2
        a, b, c = check(document)
        ok, wrong, unresolved = ok + a, wrong + b, unresolved + c
    print(
        "%d citation(s) verified, %d wrong, %d naming a file outside the package"
        % (ok, wrong, unresolved)
    )
    return 1 if wrong else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
