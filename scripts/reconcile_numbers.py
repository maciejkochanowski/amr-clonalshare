#!/usr/bin/env python3
"""Check every number the article quotes against the record it comes from.

    python scripts/reconcile_numbers.py --manuscript manuscript.md \
        --supplement supplement.md --supplement2 supplement2.md [--manifest scripts/reconcile_numbers.json]

The manifest lists claims. Each claim names the document it appears in, the
literal text the document must carry, and an expression that computes the
same quantity from the stored records (the example records, the empirical
analyses, the confirmatory verdicts and summaries, the campaign accounting,
the profile, the seed-stability check and the mutation record), with the
format that turns the value into the literal. The script evaluates every
expression, formats it, and reports a claim as failed when the formatted value
is not the literal, or when the literal does not appear in the document. It
exits non-zero on any failure, so a rebuilt article can be checked in one
step after the records are regenerated.

A claim may carry a context, a phrase with "{}" where the literal stands
(alternatives separated by "|"), so that a short number is looked for in its
sentence rather than anywhere in the document.

The documents may be given as the Markdown sources or as the text of the
built PDFs (`pdftotext -layout`); several files per document are read as one.
Unicode minus signs, thin spaces and non-breaking spaces are folded before
the search, so the literal in the manifest is written with an ASCII hyphen.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "scripts" / "reconcile_numbers.json"

WORDS = {0: "zero", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five",
         6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten", 11: "eleven",
         12: "twelve", 13: "thirteen", 14: "fourteen", 15: "fifteen",
         16: "sixteen", 17: "seventeen", 18: "eighteen", 19: "nineteen",
         20: "twenty", 30: "thirty", 40: "forty", 50: "fifty"}

_cache: dict[str, Any] = {}


def J(path: str):
    """A JSON record of the repository, read once."""
    key = "J:" + path
    if key not in _cache:
        _cache[key] = json.loads((ROOT / path).read_text(encoding="utf-8"))
    return _cache[key]


def C(path: str) -> list[dict[str, str]]:
    """The rows of a CSV file of the repository, as strings."""
    key = "C:" + path
    if key not in _cache:
        with (ROOT / path).open(encoding="utf-8", newline="") as handle:
            _cache[key] = list(csv.DictReader(handle))
    return _cache[key]


def T(path: str) -> str:
    key = "T:" + path
    if key not in _cache:
        _cache[key] = (ROOT / path).read_text(encoding="utf-8")
    return _cache[key]


CONSTANTS: dict[str, Any] = {}


def K(name: str):
    """A constant of the manifest (a commit the text names, for instance)."""
    return CONSTANTS[name]


def G(pattern: str) -> list[str]:
    """Files of the repository matching a glob, relative to its root."""
    return sorted(p.relative_to(ROOT).as_posix() for p in ROOT.glob(pattern))


def rows(path: str, **where) -> list[dict[str, str]]:
    """CSV rows whose named columns equal the given values."""
    out = []
    for row in C(path):
        if all(str(row.get(k)) == str(v) for k, v in where.items()):
            out.append(row)
    return out


def row(path: str, **where) -> dict[str, str]:
    found = rows(path, **where)
    if len(found) != 1:
        raise ValueError(f"{path}: {len(found)} rows match {where}")
    return found[0]


def f(x) -> float:
    return float(x)


def word(n) -> str:
    return WORDS[int(n)]


def pct(x, digits=1) -> str:
    return f"{100.0 * float(x):.{digits}f}%"


def signed(x, digits=3) -> str:
    return f"{float(x):+.{digits}f}"


def argmin(d: dict) -> str:
    return min(d, key=lambda k: d[k])


def argmax(d: dict) -> str:
    return max(d, key=lambda k: d[k])


def betainv(q: float, a: float, b: float) -> float:
    """Quantile of the beta distribution, for a Clopper-Pearson bound."""
    from scipy.stats import beta
    return float(beta.ppf(q, a, b))


def fold(text: str) -> str:
    """Fold the typographic characters of a built document to ASCII."""
    for old, new in (("\u2212", "-"), ("\u2013", "-"), ("\u2264", "<="), ("\u00a0", " "),
                     ("\u202f", " "), ("\u2009", " "), ("\u2019", "'"), ("\u00d7", "x")):
        text = text.replace(old, new)
    return text


def load_documents(paths: list[str]) -> str:
    parts = []
    for path in paths:
        parts.append(Path(path).read_text(encoding="utf-8", errors="replace"))
    text = fold("\n".join(parts))
    # Line breaks inside a sentence and the column gaps of a table laid out
    # as text are one space for the search, and so are the cell separators
    # of a Markdown table.
    return re.sub(r"\s+", " ", text.replace("|", " "))


def evaluate(expr: str):
    env = {"J": J, "C": C, "T": T, "G": G, "K": K, "ROOT": ROOT, "re": re, "rows": rows, "row": row, "f": f, "word": word, "betainv": betainv,
           "pct": pct, "signed": signed, "argmin": argmin, "argmax": argmax,
           "math": math, "statistics": statistics, "min": min, "max": max,
           "sum": sum, "len": len, "sorted": sorted, "abs": abs, "round": round,
           "int": int, "float": float, "str": str, "all": all, "any": any,
           "set": set, "list": list, "dict": dict, "enumerate": enumerate,
           "zip": zip, "range": range}
    env["__builtins__"] = {}
    # The helpers are globals, so that a lambda or a comprehension inside an
    # expression sees them too.
    return eval(expr, env)  # the manifest is part of the repository


def render(value, fmt: str | None) -> str:
    if fmt is None or fmt == "":
        return str(value)
    if fmt == "word":
        return word(value)
    if fmt == "Word":
        return word(value).capitalize()
    return fmt.format(value)


def check(claim: dict, docs: dict[str, str]) -> tuple[bool, str]:
    text = re.sub(r"\s+", " ", fold(str(claim["text"])))
    try:
        value = evaluate(claim["expr"])
        got = fold(render(value, claim.get("format")))
    except Exception as error:  # a broken expression is a failed claim
        return False, f"{claim['id']}: expression failed: {error!r}"
    if text == "<generated>":
        text = re.sub(r"\s+", " ", got)
    elif got != text:
        return False, f"{claim['id']}: the records give {got!r}, the article says {text!r}"
    if text == "":
        return True, f"{claim['id']}: no row"
    missing = []
    for name in claim["where"]:
        doc = docs.get(name)
        if doc is None:
            continue
        contexts = [re.sub(r"\s+", " ", fold(c)) for c in str(claim.get("context", "{}")).split("|")]
        count = sum(doc.count(c.replace("{}", text)) for c in contexts)
        if count == 0:
            missing.append(name)
        elif claim.get("at_least", 1) > count:
            missing.append(f"{name} ({count} of {claim['at_least']} occurrences)")
    if missing:
        return False, f"{claim['id']}: {text!r} not found in {', '.join(missing)}"
    return True, f"{claim['id']}: {text}"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--manuscript", nargs="*", default=[], help="text or Markdown files of the manuscript")
    parser.add_argument("--supplement", nargs="*", default=[], help="text or Markdown files of the supplement")
    parser.add_argument("--supplement2", nargs="*", default=[], help="text or Markdown files of the second supplement (the E. coli example)")
    parser.add_argument("--manifest", default=str(MANIFEST))
    parser.add_argument("--only", default=None, help="check the claims whose id starts with this prefix")
    parser.add_argument("--no-tables", action="store_true",
                        help="skip the table cells (ids starting with 'table'); a cell laid out over two lines "
                             "of a PDF's text cannot be found as one string, so check tables on the Markdown sources")
    parser.add_argument("--quiet", action="store_true", help="print failures only")
    args = parser.parse_args(argv)
    manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    CONSTANTS.update(manifest.get("constants", {}))
    docs = {}
    if args.manuscript:
        docs["manuscript"] = load_documents(args.manuscript)
    if args.supplement:
        docs["supplement"] = load_documents(args.supplement)
    if args.supplement2:
        docs["supplement2"] = load_documents(args.supplement2)
    if not docs:
        print("no document given; the expressions are evaluated and formatted only", file=sys.stderr)
    failures = 0
    checked = 0
    for claim in manifest["claims"]:
        if args.only and not claim["id"].startswith(args.only):
            continue
        if args.no_tables and claim["id"].startswith("table"):
            continue
        if docs and not any(name in docs for name in claim["where"]):
            continue
        checked += 1
        ok, message = check(claim, docs)
        if not ok:
            failures += 1
            print("FAIL " + message)
        elif not args.quiet:
            print("ok   " + message)
    print(f"{checked} claims checked, {failures} failed")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
