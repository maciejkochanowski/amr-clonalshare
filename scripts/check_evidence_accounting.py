#!/usr/bin/env python3
"""Hold the shipped atlas records to the accounting they claim.

A table that reports fewer units than it held, or refuses one without saying
why, or prints a number that was never estimated, is read as selection. The
checks below are the ones that would have caught each of those in the records
this package ships, and they run on the files themselves rather than on the
code that wrote them, so a record edited by hand is checked too.

    python scripts/check_evidence_accounting.py [directory ...]

Exit status is 1 when any check fails, so it can gate a release.
"""
from __future__ import annotations

import json
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
REFUSALS = {
    "fewer than 50 isolates tested",
    "no variance",
    "fewer than 2 lineages in the analysed subset",
}


def finite(x) -> bool:
    return isinstance(x, (int, float)) and not math.isnan(float(x))


def check_cohort(where: str, present, agents: dict, out: list) -> None:
    """Every agent present leaves one record, refused ones say why, and no
    reported cell carries an interval that was never estimated."""
    if isinstance(present, int) and len(agents) != present:
        out.append("%s: %d agents present but %d recorded"
                   % (where, present, len(agents)))
    for name, cell in agents.items():
        if not isinstance(cell, dict):
            out.append("%s/%s: record is not an object" % (where, name))
            continue
        reason = cell.get("skipped")
        if reason:
            if reason not in REFUSALS:
                out.append("%s/%s: refusal reason outside the closed set: %r"
                           % (where, name, reason))
            continue
        lo, hi = cell.get("ci_low"), cell.get("ci_high")
        if finite(lo) and finite(hi) and hi - lo == 0:
            out.append("%s/%s: reported with an interval of zero width, which "
                       "is the shape of a cell that was never estimated"
                       % (where, name))
        n_groups = (cell.get("realised") or {}).get("n_groups")
        if isinstance(n_groups, int) and n_groups < 2:
            out.append("%s/%s: reported on %d lineage(s); a share needs two"
                       % (where, name, n_groups))
        realised = cell.get("realised")
        if isinstance(realised, dict) and not realised.get("estimable"):
            if not realised.get("reason"):
                out.append("%s/%s: the realised estimand was refused without a "
                           "stated reason" % (where, name))


def check_file(path: pathlib.Path, out: list) -> int:
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
    except (ValueError, OSError) as exc:
        out.append("%s: unreadable (%s)" % (path.name, exc))
        return 0

    seen = 0
    if isinstance(doc, dict) and isinstance(doc.get("cohorts"), list):
        for cohort in doc["cohorts"]:
            label = "%s/%s/%s" % (cohort.get("organism"),
                                  cohort.get("host_group"),
                                  cohort.get("matrix"))
            present = cohort.get("n_agents_present")
            recorded = cohort.get("n_agents_recorded")
            missing = cohort.get("n_agents_not_recorded")
            if isinstance(missing, int) and isinstance(recorded, int):
                # A record restated after the fact reconciles by count.
                if recorded + missing != present:
                    out.append("%s: %d recorded plus %d not recorded is not "
                               "%s present" % (label, recorded, missing,
                                               present))
                present = None
            check_cohort(label, present, cohort.get("agents") or {}, out)
            seen += 1
    elif isinstance(doc, dict) and "agents" in doc:
        present = doc.get("n_drugs_available")
        recorded = doc.get("n_agents_recorded")
        missing = doc.get("n_agents_not_recorded")
        if isinstance(missing, int) and isinstance(recorded, int):
            # A record restated after the fact reconciles by count, because
            # the raw release it was built from is not retained.
            if recorded + missing != present:
                out.append("%s: %d recorded plus %d not recorded is not %s "
                           "present" % (path.stem, recorded, missing, present))
            present = None
        check_cohort(doc.get("organism") or path.stem, present, doc["agents"],
                     out)
        seen += 1
    return seen


def main(argv: list[str]) -> int:
    roots = [pathlib.Path(a) for a in argv[1:]]
    if not roots:
        roots = [ROOT / "paper" / "softwarex" / "evidence",
                 ROOT / "benchmarks"]
    problems: list[str] = []
    cohorts = files = 0
    for root in roots:
        if not root.exists():
            continue
        for path in sorted(root.rglob("*.json")):
            n = check_file(path, problems)
            if n:
                files += 1
                cohorts += n
    for line in problems:
        print(line)
    print("%d file(s), %d cohort(s) checked, %d problem(s)"
          % (files, cohorts, len(problems)))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
