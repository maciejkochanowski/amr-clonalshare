#!/usr/bin/env python3
"""The clonal share across species, on public NCBI Pathogen Detection data.

    bash benchmarks/fetch_pathogen_detection.sh <raw_dir>
    python benchmarks/atlas_cross_species.py --raw <raw_dir> --out <dir>

WHAT THIS IS FOR. The worked examples in the article are one species, one
laboratory and one panel, so they cannot show whether the estimator behaves
the same way elsewhere. This script runs the same estimator on a source that
shares none of those choices: the lineage is the SNP cluster NCBI assigns
(PDS accession), the trait is non-susceptibility read from the AST_phenotypes
metadata field, and both come from the same public release, so nothing here
depends on a cut-off we derived.

WHAT WOULD MEAN IT MEASURED THE WRONG THING. Every species is also run with
the lineage labels permuted within the analysed subset. That arm has to return
a share indistinguishable from zero. If it does not, the estimator is reading
something other than lineage and the whole table is void.

WHAT THE NUMBERS ARE NOT. The AST phenotypes are contributed by many
laboratories, under different standards and over different years, and the SNP
clusters are recomputed at every release. A share here is therefore a
measurement of how the estimator behaves on heterogeneous public data, not an
epidemiological statement about the species.
"""
from __future__ import annotations

import argparse
import json
import statistics
from statistics import median
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from amr_clonalshare.attribution import clonal_share
from amr_clonalshare.realised import realised_share

MIN_ISOLATES = 60
MIN_DRUG_N = 50
MIN_MINOR = 0.02          # an agent needs both classes present to have variance
FOLDS = 5
REPEATS = 10
N_BOOT = 300
N_PERM = 150
SEED = 42
PERMUTATION_SEED = 20260901


def parse_ast(field: str) -> dict[str, float]:
    """Read one AST_phenotypes cell into agent -> non-susceptible indicator."""
    out: dict[str, float] = {}
    if not isinstance(field, str):
        return out
    for part in field.strip('"').split(","):
        if "=" not in part:
            continue
        drug, _, call = part.rpartition("=")
        call = call.strip().upper()
        drug = drug.strip().lower()
        if call in ("S", "I", "R"):
            out[drug] = 1.0 if call in ("I", "R") else 0.0
    return out


def run_species(raw: Path, organism: str) -> dict:
    ast = pd.read_csv(raw / f"{organism}.ast.tsv", sep="\t", dtype=str,
                      on_bad_lines="skip")
    clusters = pd.read_csv(raw / f"{organism}.clusters.tsv", sep="\t",
                           dtype=str, on_bad_lines="skip")
    if "PDS_acc" not in clusters.columns or "target_acc" not in ast.columns:
        return {"organism": organism, "error": "missing join columns"}
    joined = ast.merge(clusters[["target_acc", "PDS_acc"]], on="target_acc",
                       how="inner")
    joined = joined[joined["PDS_acc"].notna() & (joined["PDS_acc"] != "NULL")]
    if len(joined) < MIN_ISOLATES:
        return {"organism": organism, "n": int(len(joined)),
                "error": "too few joined isolates"}

    rows = [parse_ast(s) for s in joined["AST_phenotypes"]]
    agents = sorted({d for r in rows for d in r})
    lineage = joined["PDS_acc"].tolist()
    sizes = pd.Series(lineage).value_counts()
    support = float(sizes[sizes >= 2].sum() / len(lineage))

    rng = np.random.default_rng(PERMUTATION_SEED)
    out = {"organism": organism, "n_isolates": int(len(joined)),
           "n_lineages": int(len(set(lineage))), "support": support,
           "n_drugs_available": len(agents), "agents": {}}
    for agent in agents:
        y = np.array([r.get(agent, np.nan) for r in rows], dtype=float)
        ok = np.isfinite(y)
        if ok.sum() < MIN_DRUG_N:
            continue
        yy = y[ok]
        ll = [lineage[i] for i in np.where(ok)[0]]
        prevalence = float(yy.mean())
        if prevalence < MIN_MINOR or prevalence > 1 - MIN_MINOR:
            out["agents"][agent] = {"n": int(ok.sum()),
                                    "prevalence": prevalence,
                                    "skipped": "no variance"}
            continue
        # The labels of the analysed subset are permuted, not those of the
        # whole cohort: permuting globally and then subsetting draws labels
        # from a larger pool, so the lineage size distribution would differ
        # from the real one and the control would not be comparable.
        permuted = list(rng.permutation(ll))
        # The realised share answers the other question on the same data: how
        # much sits between the lineages this cohort holds, rather than between
        # lineages drawn afresh. It is what a cohort with few lineages can
        # actually support, so it is recorded for every cell.
        conditional = realised_share(yy, ll)
        real = clonal_share(yy, ll, folds=FOLDS, repeats=REPEATS,
                            n_boot=N_BOOT, n_perm=N_PERM, seed=SEED).as_dict()
        null = clonal_share(yy, permuted, folds=FOLDS, repeats=REPEATS,
                            n_boot=N_BOOT, n_perm=N_PERM,
                            seed=SEED).as_dict()
        # kappa_adj is the bias-corrected estimate and the one reported;
        # kappa is the raw skill, biased low by the cross-validation penalty.
        out["agents"][agent] = {
            "n": int(ok.sum()), "prevalence": prevalence,
            "kappa": real["kappa_adj"], "kappa_raw": real["kappa"],
            "null_mean": real["null_mean"],
            "ci_low": real["ci_low"], "ci_high": real["ci_high"],
            "support": real["support"], "estimable": bool(real["estimable"]),
            "kappa_permuted": null["kappa_adj"],
            "kappa_permuted_raw": null["kappa"],
            "permuted_ci": [null["ci_low"], null["ci_high"]],
            "realised": {
                "kappa": conditional.kappa,
                "ci_low": conditional.ci_low,
                "ci_high": conditional.ci_high,
                "superpopulation_low": conditional.superpopulation_low,
                "superpopulation_high": conditional.superpopulation_high,
                "n_groups": conditional.n_groups,
                "effective_group_size": conditional.effective_group_size,
                "residual_excess_kurtosis":
                    conditional.residual_excess_kurtosis,
                "estimable": conditional.estimable,
                "reason": conditional.reason}}
        print(f"  {agent:<30} n={ok.sum():>6} p={prevalence:.3f} "
              f"kappa={real['kappa_adj']:+.3f} "
              f"permuted={null['kappa_adj']:+.3f}", flush=True)
    return out


def build_table(results: list[dict]) -> list[dict]:
    table = []
    for record in results:
        for agent, cell in (record.get("agents") or {}).items():
            if "kappa" not in cell:
                continue
            row = {"species": record["organism"], "agent": agent,
                   "kappa": cell["kappa"],
                   "kappa_permuted": cell["kappa_permuted"],
                   "prevalence": cell["prevalence"],
                   "ci_low": cell["ci_low"],
                   "ci_high": cell["ci_high"]}
            conditional = cell.get("realised") or {}
            row["realised_kappa"] = conditional.get("kappa")
            row["realised_ci_low"] = conditional.get("ci_low")
            row["realised_ci_high"] = conditional.get("ci_high")
            row["realised_estimable"] = conditional.get("estimable")
            row["superpopulation_ci_low"] = conditional.get(
                "superpopulation_low")
            row["superpopulation_ci_high"] = conditional.get(
                "superpopulation_high")
            row["residual_excess_kurtosis"] = conditional.get(
                "residual_excess_kurtosis")
            row["n_lineages_analysed"] = conditional.get("n_groups")
            table.append(row)
    return table


def spearman(a: list[float], b: list[float]) -> float:
    """Rank correlation, written out to keep the dependency list unchanged."""
    def ranks(values: list[float]) -> list[float]:
        order = sorted(range(len(values)), key=lambda i: values[i])
        out = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while (j + 1 < len(order)
                   and values[order[j + 1]] == values[order[i]]):
                j += 1
            shared = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                out[order[k]] = shared
            i = j + 1
        return out

    x, y = ranks(a), ranks(b)
    n = len(x)
    if n < 2:
        return 0.0
    mx, my = sum(x) / n, sum(y) / n
    num = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    den = (sum((xi - mx) ** 2 for xi in x)
           * sum((yi - my) ** 2 for yi in y)) ** 0.5
    return num / den if den else 0.0


def summarise(table: list[dict], results: list[dict]) -> dict:
    observed = [r["kappa"] for r in table]
    permuted = [r["kappa_permuted"] for r in table]
    per_species = {}
    for record in results:
        cells = [r for r in table if r["species"] == record["organism"]]
        if not cells:
            continue
        per_species[record["organism"]] = {
            "n_isolates": record["n_isolates"],
            "n_lineages": record["n_lineages"],
            "support": record["support"],
            "n_agents": len(cells),
            "kappa_median": statistics.median(c["kappa"] for c in cells),
            "kappa_min": min(c["kappa"] for c in cells),
            "kappa_max": max(c["kappa"] for c in cells),
            "kappa_permuted_median": statistics.median(
                c["kappa_permuted"] for c in cells)}
    permuted_cells = [cell for record in results
                      for cell in (record.get("agents") or {}).values()
                      if "permuted_ci" in cell]
    sizes = [v["n_isolates"] for v in per_species.values()]
    medians = [v["kappa_median"] for v in per_species.values()]
    return {
        "generated": date.today().isoformat(),
        "cells": len(table),
        "species_reported": len(per_species),
        "isolates_reported": sum(sizes),
        "lineages_reported": sum(v["n_lineages"] for v in per_species.values()),
        "kappa_median": statistics.median(observed) if observed else None,
        "kappa_permuted_median": (statistics.median(permuted)
                                  if permuted else None),
        "kappa_permuted_mean": (statistics.fmean(permuted)
                                if permuted else None),
        "cells_ci_excludes_zero": sum(1 for r in table if r["ci_low"] > 0),
        "realised_cells_estimable": sum(1 for r in table
                                        if r.get("realised_estimable")),
        "realised_median_width": (median(
            r["realised_ci_high"] - r["realised_ci_low"] for r in table
            if r.get("realised_estimable")) if any(
                r.get("realised_estimable") for r in table) else None),
        # Both widths come from the same fit, so the comparison is of two
        # questions on one cohort and not of two estimators on two scales.
        "superpopulation_median_width": (median(
            r["superpopulation_ci_high"] - r["superpopulation_ci_low"]
            for r in table if r.get("realised_estimable")) if any(
                r.get("realised_estimable") for r in table) else None),
        "permuted_cells_ci_excludes_zero": sum(
            1 for c in permuted_cells
            if c["permuted_ci"][0] > 0.0 or c["permuted_ci"][1] < 0.0),
        "spearman_isolates_vs_median_share": spearman(sizes, medians),
        "per_species": per_species}


def write_outputs(out: Path, results: list[dict]) -> None:
    table = build_table(results)
    (out / "atlas_table.json").write_text(
        json.dumps(table, indent=1) + "\n", encoding="utf-8")
    summary = summarise(table, results)
    (out / "atlas_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"\ncells {summary['cells']} across "
          f"{summary['species_reported']} species; median kappa "
          f"{summary['kappa_median']:.3f}, median permuted "
          f"{summary['kappa_permuted_median']:+.4f}")
    print(f"written: {out}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, required=True,
                        help="directory holding <organism>.ast.tsv and "
                             "<organism>.clusters.tsv")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--organism", action="append",
                        help="restrict to these organisms; default is every "
                             "pair found under --raw")
    parser.add_argument("--aggregate-only", action="store_true",
                        help="skip the fitting and rebuild the table and the "
                             "summary from the per-organism files already "
                             "under --out, which is how the species are run "
                             "one to a task on a cluster")
    args = parser.parse_args()

    if args.aggregate_only:
        results = [json.loads(p.read_text(encoding="utf-8"))
                   for p in sorted((args.out / "results").glob("*.json"))]
        write_outputs(args.out, results)
        return 0

    organisms = args.organism or sorted(
        p.name[:-len(".ast.tsv")] for p in args.raw.glob("*.ast.tsv"))
    if not organisms:
        raise SystemExit(f"no *.ast.tsv under {args.raw}")

    (args.out / "results").mkdir(parents=True, exist_ok=True)
    results = []
    for organism in organisms:
        print(f"=== {organism}", flush=True)
        try:
            record = run_species(args.raw, organism)
        except Exception as exc:                      # noqa: BLE001
            record = {"organism": organism,
                      "error": f"{type(exc).__name__}: {exc}"}
            print(f"  failed: {record['error']}", flush=True)
        results.append(record)
        (args.out / "results" / f"{organism}.json").write_text(
            json.dumps(record, indent=1) + "\n", encoding="utf-8")

    write_outputs(args.out, results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
