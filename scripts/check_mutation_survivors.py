"""Fail on a surviving mutant that nobody has accounted for.

`mutmut run` leaves a list of mutants the tests did not kill. Some are
equivalent to the original program -- a dtype that was already right, a bound
that the admissible inputs never reach -- and no test can kill them, because
there is nothing to observe; deciding which ones those are is not something a
program can do. Others are real gaps, and saying so is more useful than a
score. Both are written down, by class, with a verdict and a reason, in
`tests/mutation_survivors.json`, and this script compares a fresh run against
that record:

    mutmut results > survivors.txt
    python scripts/check_mutation_survivors.py survivors.txt

It fails when a mutant survives that the record does not name, which is a new
gap; when the record names a mutant that no longer survives, which is a stale
entry to delete; and when a class carries no verdict or no reason. That turns
a list of survivors from a note in a changelog into something that cannot rot
unnoticed.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RECORD = ROOT / "tests" / "mutation_survivors.json"
LINE = re.compile(r"^\s*(?P<name>[\w.]+):\s*(?P<status>\w+)\s*$")
VERDICTS = ("equivalent", "gap")


def read_results(path: Path) -> set[str]:
    """The mutants a run reports as survived; timeouts are not survivors."""
    survivors = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        found = LINE.match(line)
        if found and found.group("status") == "survived":
            survivors.add(found.group("name"))
    return survivors


def read_record(path: Path = RECORD) -> tuple[dict[str, str], list[str]]:
    """Map every named mutant to its class, and report a malformed class."""
    stored = json.loads(path.read_text(encoding="utf-8"))
    placed: dict[str, str] = {}
    problems = []
    for name, entry in stored["classes"].items():
        if entry.get("verdict") not in VERDICTS:
            problems.append(f"class {name}: verdict must be one of {VERDICTS}")
        if not entry.get("reason"):
            problems.append(f"class {name}: no reason given")
        if entry.get("verdict") == "gap" and not entry.get("what_would_kill_it"):
            problems.append(f"class {name}: a gap must say what would kill it")
        for mutant in entry["mutants"]:
            if mutant in placed:
                problems.append(f"{mutant}: in both {placed[mutant]} and {name}")
            placed[mutant] = name
    return placed, problems


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        print("usage: check_mutation_survivors.py <mutmut results output>",
              file=sys.stderr)
        return 2
    survived = read_results(Path(argv[0]))
    placed, problems = read_record()
    unaccounted = sorted(survived - set(placed))
    stale = sorted(set(placed) - survived)
    for name in unaccounted:
        problems.append(f"unaccounted survivor: {name}")
    for name in stale:
        problems.append(f"stale entry, this mutant no longer survives: {name}")
    for problem in problems:
        print(problem, file=sys.stderr)
    print(f"{len(survived)} survivors, {len(placed)} accounted for, "
          f"{len(unaccounted)} not, {len(stale)} stale")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
