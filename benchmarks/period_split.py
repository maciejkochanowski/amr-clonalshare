#!/usr/bin/env python3
"""The clonal share of nalidixic acid non-susceptibility in poultry, by period.

    python period_split.py --raw <raw> --taxonomy <dir> --out <file.json>

WHAT THIS ANSWERS. verify_vet_claims.py showed that non-susceptibility in
poultry animal_food rose from 0.0059 before 2016 to 0.1823 from 2016 onward,
and that one serovar carries five sixths of it. The question here is whether
the share of that resistance carried by lineage is a property of the whole
cohort or of the late period only, and whether the rise is one lineage entering
the cohort rather than resistance spreading across lineages.

DESIGN. One cell (poultry, animal_food, nalidixic acid), two typing resolutions
(NCBI PDS SNP cluster; serovar), three frames (all years; before 2016; 2016
onward). Same estimators, same settings as resolution_atlas.py. Permutation
control at every resolution and frame.

VOID CONDITIONS, stated before the run. A frame with fewer than 20
non-susceptible isolates is not a test and is reported as void. A permutation
control that leaves zero at any resolution voids that resolution. If the
pre-2016 frame holds fewer than 20 non-susceptible isolates the comparison
of shares across periods is void and only the counts are reported.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from amr_clonalshare.attribution import clonal_share
from amr_clonalshare.realised import realised_share

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_vet_claims import parse_ast, serovar, year, ORGANISM, AGENT

AGENTS = ("nalidixic acid", "ciprofloxacin")
FOLDS, REPEATS, N_BOOT, N_PERM = 5, 10, 300, 150
SEED, PERMUTATION_SEED = 42, 20260902
#: The registered cut is 2016, anchored on the 2014 first detection of the
#: emergent Infantis clone in United States retail meat with a lag to retail
#: sampling; it was fixed before any period share was computed. 2015 and 2017
#: are run as a sensitivity so that the reading does not rest on one year.
PERIOD_CUT = 2016
SENSITIVITY_CUTS = (2015, 2017)
MIN_NS = 20
#: "-:r:1,5" is the antigenic formula of Infantis (6,7,14:r:1,5) with the O
#: group unread. A curated arm merges it into Infantis so the cost of taking
#: serovar strings as written is measured rather than asserted.
SEROVAR_MERGE = {"-:r:1,5": "Infantis"}


def share_block(y, lab, tag):
    y = np.asarray(y, float)
    lab = np.asarray(lab, object)
    out = {"n": int(y.size), "n_ns": int(y.sum()), "prevalence": float(y.mean()),
           "n_groups": int(pd.unique(lab).size)}
    out["void"] = bool(y.sum() < MIN_NS)
    if out["void"]:
        return out
    r = clonal_share(y, lab, folds=FOLDS, repeats=REPEATS, n_boot=N_BOOT,
                     n_perm=N_PERM, seed=SEED)
    out.update(clonal=float(r.kappa_adj), lo=float(r.ci_low), hi=float(r.ci_high),
               support=float(r.support), estimable=bool(r.estimable),
               p_value=float(r.p_value))
    rng = np.random.default_rng(PERMUTATION_SEED)
    p = clonal_share(y, rng.permutation(lab), folds=FOLDS, repeats=REPEATS,
                     n_boot=N_BOOT, n_perm=N_PERM, seed=SEED)
    out.update(perm=float(p.kappa_adj), perm_lo=float(p.ci_low),
               perm_hi=float(p.ci_high))
    rs = realised_share(y, lab)
    out.update(realised=float(rs.kappa), realised_estimable=bool(rs.estimable),
               realised_reason=rs.reason)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True, type=Path)
    ap.add_argument("--taxonomy", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    sys.path.insert(0, str(a.taxonomy))
    from source_taxonomy import classify

    ast = pd.read_csv(a.raw / f"{ORGANISM}.ast.tsv", sep="\t", dtype=str,
                      on_bad_lines="skip")
    clusters = pd.read_csv(a.raw / f"{ORGANISM}.clusters.tsv", sep="\t",
                           dtype=str, on_bad_lines="skip")
    df = ast.merge(clusters[["target_acc", "PDS_acc"]], on="target_acc")
    df = df[df["PDS_acc"].notna() & (df["PDS_acc"] != "NULL")]
    lab = [classify(h, s) for h, s in zip(df["host"], df["isolation_source"])]
    df = df.assign(host_group=[c["group"] for c in lab],
                   matrix=[c["matrix"] for c in lab],
                   serovar=[serovar(n) for n in df["scientific_name"]],
                   year=[year(v) for v in df["collection_date"]],
                   ns=[parse_ast(v).get(AGENT, np.nan) for v in df["AST_phenotypes"]])
    base = df[(df["host_group"] == "poultry") & (df["matrix"] == "animal_food")
              & (df["serovar"] != "") & df["year"].notna()]
    res = {}
    for agent in AGENTS:
        calls = [parse_ast(v).get(agent, np.nan) for v in base["AST_phenotypes"]]
        cell = base.assign(ns=calls)
        cell = cell[cell["ns"].notna()]
        frames = {"all": cell}
        for cut in (PERIOD_CUT,) + SENSITIVITY_CUTS:
            frames[f"before_{cut}"] = cell[cell["year"] < cut]
            frames[f"from_{cut}"] = cell[cell["year"] >= cut]
        res[agent] = {}
        for fname, fr in frames.items():
            curated = fr["serovar"].replace(SEROVAR_MERGE)
            res[agent][fname] = {
                "PDS": share_block(fr["ns"], fr["PDS_acc"], "PDS"),
                "SEROVAR": share_block(fr["ns"], fr["serovar"], "SEROVAR"),
                "SEROVAR_CURATED": share_block(fr["ns"], curated, "SEROVAR_CURATED")}
            top = (fr.groupby("serovar")["ns"].agg(n="size", ns="sum")
                   .sort_values("ns", ascending=False).head(5))
            res[agent][fname]["top_serovars_by_ns"] = json.loads(
                top.reset_index().to_json(orient="records"))
            res[agent][fname]["country"] = {
                str(k): int(v) for k, v in
                fr["geo_loc_name"].fillna("(missing)").str.split(":").str[0]
                .value_counts().head(8).items()}
    cell = base

    report = {"provenance": {"date": date.today().isoformat(),
                             "python": platform.python_version(),
                             "numpy": np.__version__, "pandas": pd.__version__,
                             "folds": FOLDS, "repeats": REPEATS, "n_boot": N_BOOT,
                             "n_perm": N_PERM, "seed": SEED,
                             "permutation_seed": PERMUTATION_SEED,
                             "period_cut": PERIOD_CUT,
                             "sensitivity_cuts": SENSITIVITY_CUTS,
                             "serovar_merge": SEROVAR_MERGE, "min_ns": MIN_NS},
              "cell": {"host": "poultry", "matrix": "animal_food", "agents": AGENTS,
                       "n_serotyped_dated": int(len(cell))},
              "frames": res}
    a.out.write_text(json.dumps(report, indent=2))
    print(json.dumps(report["frames"], indent=1)[:3000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
