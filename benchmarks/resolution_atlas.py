#!/usr/bin/env python3
"""The clonal share of resistance at two typing resolutions.

    python resolution_atlas.py --raw <raw> --out <dir> --cohort <i>
    python resolution_atlas.py --out <dir> --aggregate

Registered in prereg_resolution.json before this script was run.

WHAT THIS ANSWERS. The published veterinary atlas defines a lineage as an NCBI
Pathogen Detection SNP cluster. At that resolution most animal isolates are
singletons, so within-lineage variation is scarce and the estimator refuses the
cell. Of 256 candidate host contrasts in that atlas none was identifiable. The
estimand, however, is the share carried by the lineages a cohort holds, and
lineage is a definition. This script runs the same cohorts, the same agents and
the same estimators at a second, coarser definition taken from the serovar, and
reports the pair.

WHY SEROVAR. A veterinary diagnostic laboratory serotypes routinely and rarely
holds SNP clusters, so a share at serovar resolution is one it can reproduce.

WHAT WOULD MEAN IT MEASURED THE WRONG THING. A coarser grouping can manufacture
between-group variance rather than measure it. The permutation control is
therefore run at BOTH resolutions: labels are permuted inside the analysed
subset, and the coarse arm is void if that control stops returning zero.
"""
from __future__ import annotations

import argparse
import json
import platform
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from amr_clonalshare.attribution import clonal_share
from amr_clonalshare.realised import realised_share

from vet_source_taxonomy import ANIMAL_GROUPS, classify

ORGANISM = "Salmonella"
MIN_ISOLATES = 60
MIN_DRUG_N = 50
MIN_MINOR = 0.02
FOLDS = 5
REPEATS = 10
N_BOOT = 300
N_PERM = 150
SEED = 42
PERMUTATION_SEED = 20260902
GROUPS = ANIMAL_GROUPS + ("human",)
SEROVAR_MARK = " serovar "


def parse_ast(field: str) -> dict[str, float]:
    """One AST_phenotypes cell as agent -> non-susceptible indicator."""
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


def serovar(name) -> str:
    """The serovar, or the empty string when the record stops above it.

    An isolate recorded only as Salmonella enterica, or as the subspecies, has
    not been serotyped. It is dropped from the serovar arm rather than pooled,
    because one pooled unserotyped group would act as a large artificial
    lineage and inflate the between-lineage term it is supposed to measure.
    """
    if not isinstance(name, str) or SEROVAR_MARK not in name:
        return ""
    tail = name.split(SEROVAR_MARK, 1)[1].strip()
    return tail if tail else ""


def load(raw: Path) -> pd.DataFrame:
    ast = pd.read_csv(raw / f"{ORGANISM}.ast.tsv", sep="\t", dtype=str,
                      on_bad_lines="skip")
    clusters = pd.read_csv(raw / f"{ORGANISM}.clusters.tsv", sep="\t",
                           dtype=str, on_bad_lines="skip")
    joined = ast.merge(clusters[["target_acc", "PDS_acc"]], on="target_acc",
                       how="inner")
    joined = joined[joined["PDS_acc"].notna() & (joined["PDS_acc"] != "NULL")]
    labels = [classify(h, s) for h, s in
              zip(joined.get("host", ""), joined.get("isolation_source", ""))]
    return joined.assign(host_group=[c["group"] for c in labels],
                         matrix=[c["matrix"] for c in labels],
                         serovar=[serovar(n) for n in joined["scientific_name"]])


def cohort_plan(raw: Path) -> list[dict]:
    frame = load(raw)
    plan = []
    for group in GROUPS:
        sub = frame[frame["host_group"] == group]
        if len(sub) < MIN_ISOLATES:
            continue
        plan.append({"host_group": group, "matrix": "any", "n": int(len(sub))})
        for matrix in sorted(x for x in sub["matrix"].unique() if x):
            cut = sub[sub["matrix"] == matrix]
            if len(cut) < MIN_ISOLATES or len(cut) == len(sub):
                continue
            plan.append({"host_group": group, "matrix": matrix,
                         "n": int(len(cut))})
    return plan


def estimate(y, labels, tag):
    """Both estimators plus the permutation control, on one label vector."""
    rng = np.random.default_rng(PERMUTATION_SEED)
    permuted = list(rng.permutation(list(labels)))
    real = clonal_share(y, list(labels), folds=FOLDS, repeats=REPEATS,
                        n_boot=N_BOOT, n_perm=N_PERM, seed=SEED).as_dict()
    null = clonal_share(y, permuted, folds=FOLDS, repeats=REPEATS,
                        n_boot=N_BOOT, n_perm=N_PERM, seed=SEED).as_dict()
    cond = realised_share(y, list(labels))
    cond_null = realised_share(y, permuted)
    sizes = pd.Series(list(labels)).value_counts()
    return {
        "resolution": tag,
        "n_groups": int(len(sizes)),
        "support": float(sizes[sizes >= 2].sum() / len(labels)),
        "largest_group_share": float(sizes.iloc[0] / len(labels)),
        "clonal": {"est": real["kappa_adj"], "lo": real["ci_low"],
                   "hi": real["ci_high"], "support": real["support"],
                   "estimable": bool(real["estimable"])},
        "clonal_permuted": {"est": null["kappa_adj"], "lo": null["ci_low"],
                            "hi": null["ci_high"],
                            "excludes_zero": bool(null["ci_low"] > 0)},
        "realised": {"est": cond.kappa, "lo": cond.ci_low, "hi": cond.ci_high,
                     "n_groups": cond.n_groups,
                     "estimable": bool(cond.estimable), "reason": cond.reason},
        "realised_permuted": {"est": cond_null.kappa, "lo": cond_null.ci_low,
                              "hi": cond_null.ci_high,
                              "estimable": bool(cond_null.estimable)},
    }


def run_cohort(raw: Path, cell: dict) -> dict:
    frame = load(raw)
    sub = frame[frame["host_group"] == cell["host_group"]]
    if cell["matrix"] != "any":
        sub = sub[sub["matrix"] == cell["matrix"]]
    rows = [parse_ast(s) for s in sub["AST_phenotypes"]]
    agents = sorted({d for r in rows for d in r})
    pds = sub["PDS_acc"].tolist()
    sero = sub["serovar"].tolist()
    has_sero = [bool(s) for s in sero]

    counts = pd.Series([s for s in sero if s]).value_counts()
    out = dict(cell)
    out |= {"n_isolates": int(len(sub)),
            "n_without_serovar": int(len(sub) - sum(has_sero)),
            "serovar_drop_fraction": float(1 - sum(has_sero) / len(sub)),
            "serovar_composition": {k: int(v) for k, v in
                                    counts.head(12).items()},
            "n_serovars": int(len(counts)),
            "agents": {}}

    for agent in agents:
        y = np.array([r.get(agent, np.nan) for r in rows], dtype=float)
        ok = np.isfinite(y)
        if ok.sum() < MIN_DRUG_N:
            continue
        yy = y[ok]
        prevalence = float(yy.mean())
        if prevalence < MIN_MINOR or prevalence > 1 - MIN_MINOR:
            out["agents"][agent] = {"n": int(ok.sum()),
                                    "prevalence": prevalence,
                                    "skipped": "no variance"}
            continue
        idx = np.where(ok)[0]
        rec = {"n": int(ok.sum()), "prevalence": prevalence, "arms": {}}
        rec["arms"]["PDS"] = estimate(yy, [pds[i] for i in idx], "PDS")

        keep = [i for i in idx if has_sero[i]]
        ys = np.array([y[i] for i in keep], dtype=float)
        rec["n_serovar_arm"] = int(len(keep))
        if len(keep) < MIN_DRUG_N:
            rec["arms"]["SEROVAR"] = {"resolution": "SEROVAR",
                                      "skipped": "below MIN_DRUG_N after drop"}
        elif not (MIN_MINOR <= float(ys.mean()) <= 1 - MIN_MINOR):
            rec["arms"]["SEROVAR"] = {"resolution": "SEROVAR",
                                      "skipped": "no variance after drop",
                                      "prevalence": float(ys.mean())}
        else:
            rec["prevalence_serovar_arm"] = float(ys.mean())
            rec["arms"]["SEROVAR"] = estimate(ys, [sero[i] for i in keep],
                                              "SEROVAR")
        out["agents"][agent] = rec
        p = rec["arms"]["PDS"]["clonal"]
        s = rec["arms"].get("SEROVAR", {}).get("clonal", {})
        print(f"  {agent:<30} n={ok.sum():>5} "
              f"PDS {p['est']:+.3f} est={p['estimable']} | "
              f"SERO {s.get('est', float('nan')):+.3f} "
              f"est={s.get('estimable')}", flush=True)
    return out


def _provenance() -> dict:
    return {"generated": str(date.today()), "seed": SEED,
            "permutation_seed": PERMUTATION_SEED,
            "python": platform.python_version(),
            "numpy": np.__version__, "pandas": pd.__version__,
            "organism": ORGANISM,
            "thresholds": {"min_isolates": MIN_ISOLATES,
                           "min_isolates_per_agent": MIN_DRUG_N,
                           "min_minor_class_share": MIN_MINOR},
            "estimator": {"folds": FOLDS, "repeats": REPEATS,
                          "bootstrap": N_BOOT, "permutations": N_PERM}}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--cohort", type=int)
    ap.add_argument("--plan", action="store_true")
    ap.add_argument("--aggregate", action="store_true")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    if args.plan:
        plan = cohort_plan(args.raw)
        (args.out / "cohort_plan.json").write_text(json.dumps(plan, indent=1))
        for i, c in enumerate(plan):
            print(i, c)
        return 0

    if args.aggregate:
        cells = sorted((args.out / "cells").glob("*.json"))
        out = {"provenance": _provenance(),
               "cohorts": [json.loads(p.read_text()) for p in cells]}
        (args.out / "resolution_atlas.json").write_text(json.dumps(out, indent=1))
        print(f"aggregated {len(cells)} cohorts")
        return 0

    plan = cohort_plan(args.raw)
    if args.cohort is None or not 0 <= args.cohort < len(plan):
        print(f"cohort index must be in 0..{len(plan) - 1}")
        return 2
    cell = plan[args.cohort]
    print(f"cohort {args.cohort}: {cell}", flush=True)
    result = run_cohort(args.raw, cell)
    d = args.out / "cells"
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{args.cohort:03d}.json").write_text(json.dumps(result, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
