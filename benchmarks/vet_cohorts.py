#!/usr/bin/env python3
"""Size the veterinary cohorts before any of them is analysed.

    python vet_cohorts.py --raw <raw_dir> --out <dir>

The atlas the article already carries reads every isolate of a species
together, which mixes a hospital urine from a person with a supermarket
chicken breast. That is the right cohort for asking whether the estimator
behaves the same way across species, and the wrong one for asking anything
veterinary. This script cuts the same source by host and by where in the chain
the sample was taken, and reports how many isolates, lineages and analysable
agents each cut would have, so the analysis budget is spent on cells that can
answer something rather than on cells that will be refused.

Nothing here estimates anything. The thresholds are the ones the published
atlas already uses, so a cell counted here as analysable is analysable under
the same rule.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import pandas as pd

from vet_source_taxonomy import ANIMAL_GROUPS, classify

MIN_ISOLATES = 60
MIN_DRUG_N = 50
MIN_MINOR = 0.02


def parse_ast(field: str) -> dict[str, float]:
    out: dict[str, float] = {}
    if not isinstance(field, str):
        return out
    for part in field.strip('"').split(","):
        if "=" not in part:
            continue
        drug, _, call = part.rpartition("=")
        call = call.strip().upper()
        if call in ("S", "I", "R"):
            out[drug.strip().lower()] = 1.0 if call in ("I", "R") else 0.0
    return out


def load(raw: Path, organism: str) -> pd.DataFrame | None:
    ast_path = raw / f"{organism}.ast.tsv"
    cl_path = raw / f"{organism}.clusters.tsv"
    if not ast_path.is_file() or not cl_path.is_file():
        return None
    ast = pd.read_csv(ast_path, sep="\t", dtype=str, on_bad_lines="skip")
    clusters = pd.read_csv(cl_path, sep="\t", dtype=str, on_bad_lines="skip")
    if "PDS_acc" not in clusters.columns or "target_acc" not in ast.columns:
        return None
    joined = ast.merge(clusters[["target_acc", "PDS_acc"]], on="target_acc",
                       how="inner")
    joined = joined[joined["PDS_acc"].notna() & (joined["PDS_acc"] != "NULL")]
    if joined.empty:
        return None
    labels = [classify(h, s) for h, s in
              zip(joined.get("host", ""), joined.get("isolation_source", ""))]
    joined = joined.assign(
        host_group=[c["group"] for c in labels],
        host_name=[c["host"] for c in labels],
        matrix=[c["matrix"] for c in labels],
        rule=[c["rule"] for c in labels])
    return joined


def size_cohort(frame: pd.DataFrame) -> dict:
    rows = [parse_ast(s) for s in frame["AST_phenotypes"]]
    agents = sorted({d for r in rows for d in r})
    lineage = frame["PDS_acc"].tolist()
    sizes = Counter(lineage)
    analysable = []
    for agent in agents:
        vals = [r[agent] for r in rows if agent in r]
        if len(vals) < MIN_DRUG_N:
            continue
        p = sum(vals) / len(vals)
        if p < MIN_MINOR or p > 1 - MIN_MINOR:
            continue
        analysable.append({"agent": agent, "n": len(vals), "prevalence": p})
    return {
        "n_isolates": int(len(frame)),
        "n_lineages": int(len(sizes)),
        "support": float(sum(v for v in sizes.values() if v >= 2) / len(frame)),
        "n_agents_analysable": len(analysable),
        "agents": sorted(analysable, key=lambda a: -a["n"]),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--raw", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    organisms = sorted(p.name[:-len(".ast.tsv")]
                       for p in args.raw.glob("*.ast.tsv"))
    cohorts, skipped = {}, {}
    for organism in organisms:
        frame = load(args.raw, organism)
        if frame is None:
            skipped[organism] = "no joinable release"
            continue
        for group in ANIMAL_GROUPS + ("human",):
            sub = frame[frame["host_group"] == group]
            if len(sub) < MIN_ISOLATES:
                continue
            key = f"{organism}|{group}"
            cohorts[key] = size_cohort(sub) | {"organism": organism,
                                               "host_group": group,
                                               "matrix": "any"}
            for matrix in sorted(x for x in sub["matrix"].unique() if x):
                cut = sub[sub["matrix"] == matrix]
                if len(cut) < MIN_ISOLATES:
                    continue
                cohorts[f"{key}|{matrix}"] = size_cohort(cut) | {
                    "organism": organism, "host_group": group,
                    "matrix": matrix}

    payload = {"thresholds": {"min_isolates": MIN_ISOLATES,
                              "min_isolates_per_agent": MIN_DRUG_N,
                              "min_minor_class_share": MIN_MINOR},
               "organisms_read": organisms, "skipped": skipped,
               "cohorts": cohorts}
    (args.out / "vet_cohort_sizes.json").write_text(
        json.dumps(payload, indent=1) + "\n", encoding="utf-8")

    print(f"{'cohort':<58} {'isolates':>8} {'lineages':>8} {'support':>8} {'agents':>7}")
    total_cells = 0
    for key, c in sorted(cohorts.items(),
                         key=lambda kv: -kv[1]["n_agents_analysable"]):
        if c["n_agents_analysable"] == 0:
            continue
        if c["matrix"] == "any":
            total_cells += c["n_agents_analysable"]
        print(f"{key:<58} {c['n_isolates']:>8} {c['n_lineages']:>8} "
              f"{c['support']:>8.3f} {c['n_agents_analysable']:>7}")
    print(f"\ncells at host-group level: {total_cells}")
    print(f"cohorts written: {len(cohorts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
