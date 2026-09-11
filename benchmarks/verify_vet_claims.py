#!/usr/bin/env python3
"""Adversarial verification of the two veterinary claims in the illustrative example.

    python verify_vet_claims.py --raw <raw> --out <file.json>

CLAIM A (as written in the manuscript). The poultry-versus-swine contrast in
nalidixic acid, matrix animal_food, is carried by a lineage whose stated
mechanism is a gyrA substitution.

CLAIM B (as written). The contrast is a property of the cohort, not of a period.

WHAT WOULD MEAN THIS SCRIPT MEASURED THE WRONG THING. If the serovar column is
empty for most poultry animal_food isolates, or if fewer than 20 isolates carry
the agent in either arm, the breakdown below is not a test of either claim and
the run is void. Both conditions are checked and reported.

Nothing here is fitted. Every number is a count or a proportion read off the
same table the atlas used, at the same filters.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ORGANISM = "Salmonella"
SEROVAR_MARK = " serovar "
AGENT = "nalidixic acid"
PERIOD_CUT = 2016


def parse_ast(field):
    out = {}
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


def serovar(name):
    if not isinstance(name, str) or SEROVAR_MARK not in name:
        return ""
    return name.split(SEROVAR_MARK, 1)[1].strip()


def year(value):
    if not isinstance(value, str) or len(value) < 4:
        return np.nan
    head = value[:4]
    return float(head) if head.isdigit() else np.nan


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True, type=Path)
    ap.add_argument("--taxonomy", required=True, type=Path,
                    help="directory holding source_taxonomy.py")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    sys.path.insert(0, str(args.taxonomy))
    from source_taxonomy import classify

    ast = pd.read_csv(args.raw / f"{ORGANISM}.ast.tsv", sep="\t", dtype=str,
                      on_bad_lines="skip")
    clusters = pd.read_csv(args.raw / f"{ORGANISM}.clusters.tsv", sep="\t",
                           dtype=str, on_bad_lines="skip")
    df = ast.merge(clusters[["target_acc", "PDS_acc"]], on="target_acc",
                   how="inner")
    df = df[df["PDS_acc"].notna() & (df["PDS_acc"] != "NULL")]

    labels = [classify(h, s) for h, s in
              zip(df.get("host", ""), df.get("isolation_source", ""))]
    df = df.assign(host_group=[c["group"] for c in labels],
                   matrix=[c["matrix"] for c in labels],
                   serovar=[serovar(n) for n in df["scientific_name"]])

    date_col = next((c for c in ("collection_date", "target_creation_date",
                                 "Create date") if c in df.columns), None)
    df["year"] = ([year(v) for v in df[date_col]] if date_col
                  else [np.nan] * len(df))

    calls = [parse_ast(v) for v in df.get("AST_phenotypes", [""] * len(df))]
    df["ns"] = [c.get(AGENT, np.nan) for c in calls]

    report = {"agent": AGENT, "date_column": date_col,
              "n_joined": int(len(df)),
              "n_with_agent": int(df["ns"].notna().sum()),
              "columns": sorted(df.columns.tolist())}

    frame = df[df["ns"].notna() & (df["matrix"] == "animal_food")]
    report["void_checks"] = {}

    arms = {}
    for host in ("poultry", "swine"):
        arm = frame[frame["host_group"] == host]
        if arm.empty:
            arms[host] = {"n": 0}
            continue
        typed = arm[arm["serovar"] != ""]
        by_serovar = (typed.groupby("serovar")["ns"]
                      .agg(n="size", ns="sum").sort_values("n", ascending=False))
        by_serovar["prop_ns"] = by_serovar["ns"] / by_serovar["n"]
        top = by_serovar.head(12)

        # share of all non-susceptible isolates contributed by each serovar
        total_ns = float(typed["ns"].sum())
        top = top.assign(share_of_ns=(top["ns"] / total_ns) if total_ns else np.nan)

        with_year = arm[arm["year"].notna()]
        early = with_year[with_year["year"] < PERIOD_CUT]
        late = with_year[with_year["year"] >= PERIOD_CUT]

        arms[host] = {
            "n": int(len(arm)),
            "n_ns": int(arm["ns"].sum()),
            "prop_ns": float(arm["ns"].mean()),
            "n_serotyped": int(len(typed)),
            "serotyped_share": float(len(typed) / len(arm)),
            "total_ns_serotyped": total_ns,
            "by_serovar": json.loads(top.reset_index().to_json(orient="records")),
            "n_with_year": int(len(with_year)),
            "early": {"n": int(len(early)), "n_ns": int(early["ns"].sum()),
                      "prop_ns": float(early["ns"].mean()) if len(early) else None},
            "late": {"n": int(len(late)), "n_ns": int(late["ns"].sum()),
                     "prop_ns": float(late["ns"].mean()) if len(late) else None},
        }
        report["void_checks"][host] = {
            "serotyped_share_below_0.5": bool(len(typed) / len(arm) < 0.5),
            "fewer_than_20_ns": bool(arm["ns"].sum() < 20),
        }

    report["arms"] = arms
    args.out.write_text(json.dumps(report, indent=2))
    print(json.dumps({k: report[k] for k in
                      ("agent", "date_column", "n_joined", "n_with_agent",
                       "void_checks")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
