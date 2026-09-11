#!/usr/bin/env python3
"""The clonal share of resistance, cut by the host the isolates came from.

    python vet_atlas.py --raw <raw> --out <dir> --cohort <i>
    python vet_atlas.py --raw <raw> --out <dir> --decompose
    python vet_atlas.py --out <dir> --aggregate

WHAT THIS ANSWERS. The cross-species atlas the article already carries reads
every isolate of a species together, so a hospital urine and a supermarket
chicken breast enter the same cell. That is the right cohort for asking whether
the estimator behaves the same way across species and the wrong one for asking
anything veterinary. Here the same source, the same lineage variable and the
same estimator are used, and the only thing that changes is that the isolates
are first cut by host and by where in the chain the sample was taken.

WHY THE CUT IS THE POINT. A veterinary laboratory does not hold a species, it
holds a species in a host: broilers, or fattening pigs, or dogs. The question
it can act on is whether resistance in that host travels with the clone, which
calls for movement control and biosecurity, or arises within lineages, which
calls for stewardship. A share pooled over hosts answers neither.

REGISTERED BEFORE THE RUN. vet_expectations.json, written before this script
was run and hashed into the receipt, states four expectations and three
conditions that would mean the design measured the wrong thing.

WHAT WOULD MEAN IT MEASURED THE WRONG THING. Every cell is run again with the
lineage labels permuted inside the analysed subset; that arm must return a
share indistinguishable from zero. The permutation is inside the subset and not
over the whole species, because permuting first and subsetting afterwards would
draw labels from a larger pool and the lineage size distribution of the control
would not match the real one.
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
from amr_clonalshare.clonality import decompose_prevalence_difference
from amr_clonalshare.realised import realised_share

from agent_screen import MIN_DRUG_N, MIN_MINOR, reconcile, screen_agent
from vet_source_taxonomy import ANIMAL_GROUPS, classify

MIN_ISOLATES = 60
FOLDS = 5
REPEATS = 10
N_BOOT = 300
N_PERM = 150
SEED = 42
PERMUTATION_SEED = 20260902
DECOMPOSITION_BOOT = 2000
GROUPS = ANIMAL_GROUPS + ("human",)


def parse_ast(field: str) -> dict[str, float]:
    """One AST_phenotypes cell as agent -> non-susceptible indicator.

    Only S, I and R are calls. ND means the agent was not determined for that
    isolate and is dropped rather than read as susceptible: in one organism in
    this release almost every cell is ND, and treating it as a call would
    manufacture a cohort out of absent measurements.
    """
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
    ast_path, cl_path = raw / f"{organism}.ast.tsv", raw / f"{organism}.clusters.tsv"
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
    return joined.assign(host_group=[c["group"] for c in labels],
                         host_name=[c["host"] for c in labels],
                         matrix=[c["matrix"] for c in labels])


def cohort_plan(raw: Path) -> list[dict]:
    """Every cut worth running, in a fixed order, so a task index is stable."""
    plan = []
    for organism in sorted(p.name[:-len(".ast.tsv")]
                           for p in raw.glob("*.ast.tsv")):
        frame = load(raw, organism)
        if frame is None:
            continue
        for group in GROUPS:
            sub = frame[frame["host_group"] == group]
            if len(sub) < MIN_ISOLATES:
                continue
            plan.append({"organism": organism, "host_group": group,
                         "matrix": "any", "n": int(len(sub))})
            for matrix in sorted(x for x in sub["matrix"].unique() if x):
                cut = sub[sub["matrix"] == matrix]
                if len(cut) < MIN_ISOLATES or len(cut) == len(sub):
                    continue
                plan.append({"organism": organism, "host_group": group,
                             "matrix": matrix, "n": int(len(cut))})
    return plan


def select(frame: pd.DataFrame, cell: dict) -> pd.DataFrame:
    sub = frame[frame["host_group"] == cell["host_group"]]
    if cell["matrix"] != "any":
        sub = sub[sub["matrix"] == cell["matrix"]]
    return sub


def run_cohort(raw: Path, cell: dict) -> dict:
    frame = load(raw, cell["organism"])
    sub = select(frame, cell)
    rows = [parse_ast(s) for s in sub["AST_phenotypes"]]
    agents = sorted({d for r in rows for d in r})
    lineage = sub["PDS_acc"].tolist()
    sizes = pd.Series(lineage).value_counts()
    rng = np.random.default_rng(PERMUTATION_SEED)
    out = dict(cell)
    out |= {"n_isolates": int(len(sub)), "n_lineages": int(len(sizes)),
            "support": float(sizes[sizes >= 2].sum() / len(sub)),
            "n_agents_present": len(agents),
            "host_names": sorted({str(h) for h in sub["host_name"]
                                  if isinstance(h, str) and h}),
            "agents": {}}
    for agent in agents:
        y = np.array([r.get(agent, np.nan) for r in rows], dtype=float)
        screen = screen_agent(y, lineage)
        if not screen.analyse:
            out["agents"][agent] = screen.record
            continue
        yy, ll = screen.y, list(screen.lineage)
        prevalence = screen.record["prevalence"]
        permuted = list(rng.permutation(ll))
        conditional = realised_share(yy, ll)
        real = clonal_share(yy, ll, folds=FOLDS, repeats=REPEATS,
                            n_boot=N_BOOT, n_perm=N_PERM, seed=SEED).as_dict()
        null = clonal_share(yy, permuted, folds=FOLDS, repeats=REPEATS,
                            n_boot=N_BOOT, n_perm=N_PERM, seed=SEED).as_dict()
        out["agents"][agent] = {
            "n": screen.record["n"], "prevalence": prevalence,
            "kappa": real["kappa_adj"], "ci_low": real["ci_low"],
            "ci_high": real["ci_high"], "support": real["support"],
            "estimable": bool(real["estimable"]),
            "kappa_permuted": null["kappa_adj"],
            "permuted_ci": [null["ci_low"], null["ci_high"]],
            "realised": {"kappa": conditional.kappa,
                         "ci_low": conditional.ci_low,
                         "ci_high": conditional.ci_high,
                         "superpopulation_low": conditional.superpopulation_low,
                         "superpopulation_high": conditional.superpopulation_high,
                         "n_groups": conditional.n_groups,
                         "estimable": conditional.estimable,
                         "reason": conditional.reason}}
        print(f"  {agent:<30} n={screen.record['n']:>5} p={prevalence:.3f} "
              f"kappa={real['kappa_adj']:+.3f} perm={null['kappa_adj']:+.3f}",
              flush=True)
    reconcile(len(agents), out["agents"])
    return out


def run_decomposition(raw: Path) -> list[dict]:
    """Split a prevalence difference between two hosts of the same species.

    The question a veterinary reader asks of two hosts is not only which has
    more resistance but why: whether the hosts carry different lineages, or the
    same lineages resistant at different rates. The second is only identified
    on lineages both hosts hold, so the estimator states the shared support it
    rests on and withholds the component when that support is too thin. A
    refusal here is the expected outcome for hosts whose lineage pools barely
    overlap, and it is reported as the answer rather than worked around.
    """
    results = []
    for organism in sorted(p.name[:-len(".ast.tsv")]
                           for p in raw.glob("*.ast.tsv")):
        frame = load(raw, organism)
        if frame is None:
            continue
        present = [g for g in ANIMAL_GROUPS + ("human",)
                   if (frame["host_group"] == g).sum() >= MIN_ISOLATES]
        for i, a in enumerate(present):
            for b in present[i + 1:]:
                fa = frame[frame["host_group"] == a]
                fb = frame[frame["host_group"] == b]
                ra = [parse_ast(s) for s in fa["AST_phenotypes"]]
                rb = [parse_ast(s) for s in fb["AST_phenotypes"]]
                shared = ({d for r in ra for d in r}
                          & {d for r in rb for d in r})
                for agent in sorted(shared):
                    ya = np.array([r.get(agent, np.nan) for r in ra], float)
                    yb = np.array([r.get(agent, np.nan) for r in rb], float)
                    oa, ob = np.isfinite(ya), np.isfinite(yb)
                    if oa.sum() < MIN_DRUG_N or ob.sum() < MIN_DRUG_N:
                        continue
                    la = [fa["PDS_acc"].tolist()[i] for i in np.where(oa)[0]]
                    lb = [fb["PDS_acc"].tolist()[i] for i in np.where(ob)[0]]
                    pa, pb = float(ya[oa].mean()), float(yb[ob].mean())
                    if abs(pa - pb) < 0.02:
                        continue
                    res = decompose_prevalence_difference(
                        ya[oa], la, yb[ob], lb, n_boot=DECOMPOSITION_BOOT,
                        rng=np.random.default_rng(SEED))
                    results.append({"organism": organism, "host_a": a,
                                    "host_b": b, "agent": agent,
                                    "n_a": int(oa.sum()), "n_b": int(ob.sum()),
                                    "prevalence_a": pa, "prevalence_b": pb,
                                    "result": res})
                    print(f"  {organism} {a} vs {b} {agent}: "
                          f"{pa:.3f} vs {pb:.3f}", flush=True)
    return results


def _provenance() -> dict:
    return {"generated": date.today().isoformat(), "seed": SEED,
            "permutation_seed": PERMUTATION_SEED,
            "python": platform.python_version(), "numpy": np.__version__,
            "pandas": pd.__version__,
            "thresholds": {"min_isolates": MIN_ISOLATES,
                           "min_isolates_per_agent": MIN_DRUG_N,
                           "min_minor_class_share": MIN_MINOR},
            "estimator": {"folds": FOLDS, "repeats": REPEATS,
                          "bootstrap": N_BOOT, "permutations": N_PERM,
                          "decomposition_bootstrap": DECOMPOSITION_BOOT}}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--raw", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--cohort", type=int)
    ap.add_argument("--plan", action="store_true")
    ap.add_argument("--decompose", action="store_true")
    ap.add_argument("--aggregate", action="store_true")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    if args.plan:
        plan = cohort_plan(args.raw)
        (args.out / "cohort_plan.json").write_text(
            json.dumps(plan, indent=1) + "\n", encoding="utf-8")
        for i, c in enumerate(plan):
            print(f"{i:>3} {c['organism']:<32} {c['host_group']:<18} "
                  f"{c['matrix']:<20} n={c['n']}")
        print(f"cohorts: {len(plan)}")
        return 0

    if args.decompose:
        results = run_decomposition(args.raw)
        (args.out / "vet_decomposition.json").write_text(
            json.dumps({"provenance": _provenance(), "rows": results},
                       indent=1) + "\n", encoding="utf-8")
        print(f"decomposition rows: {len(results)}")
        return 0

    if args.aggregate:
        cells = [json.loads(p.read_text())
                 for p in sorted((args.out / "cohorts").glob("*.json"))]
        (args.out / "vet_atlas.json").write_text(
            json.dumps({"provenance": _provenance(), "cohorts": cells},
                       indent=1) + "\n", encoding="utf-8")
        print(f"cohorts aggregated: {len(cells)}")
        return 0

    plan = json.loads((args.out / "cohort_plan.json").read_text())
    cell = plan[args.cohort]
    print(f"[{args.cohort}] {cell['organism']} {cell['host_group']} "
          f"{cell['matrix']} n={cell['n']}", flush=True)
    result = run_cohort(args.raw, cell)
    (args.out / "cohorts").mkdir(exist_ok=True)
    (args.out / "cohorts" / f"{args.cohort:03d}.json").write_text(
        json.dumps(result, indent=1) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
