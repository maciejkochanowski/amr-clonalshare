#!/usr/bin/env python3
"""Assert that no module imports a sibling tool package.

Uses the AST rather than a text search, so a mention in a docstring or a comment
cannot trip the gate — which is what happened to the previous grep-based check.
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

SIBLINGS = ("bact_mdr_profiler", "bact_assoc_net", "bact_phylotrait")
SRC = Path(__file__).resolve().parents[1] / "src" / "amr_clonalshare"


def imported_names(tree: ast.AST):
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                yield a.name
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                yield node.module


def main() -> int:
    bad = []
    for f in sorted(SRC.rglob("*.py")):
        tree = ast.parse(f.read_text(encoding="utf-8"), filename=str(f))
        for name in imported_names(tree):
            root = name.split(".")[0]
            if root in SIBLINGS:
                bad.append(f"{f.name} imports {name}")
    if bad:
        print("isolation gate FAILED:", *bad, sep="\n  ", file=sys.stderr)
        return 1
    print(f"isolation gate OK ({len(list(SRC.rglob('*.py')))} modules checked)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
