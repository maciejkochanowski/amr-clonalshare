#!/usr/bin/env python3
"""Operating envelope of the prevalence-difference decomposition.

The decomposition is an exact algebraic identity, so nothing about it can be
wrong in the sense a formula can be wrong. What can be wrong is the interval:
the components are estimated from two finite samples, and a lineage that the
population holds in both collections can easily be sampled into only one. When
that happens the estimator charges the whole of that lineage to composition,
because the declared convention says a lineage that is not there is a fact
about the mix. A *sampling* zero is then read as a *population* zero, and the
within-lineage component loses part of its support.

This study measures how far that goes, over the grid a veterinary cohort
actually occupies: a few hundred isolates, tens to a hundred sequence types,
one dominant clone and a long singleton tail. It reports bias, the coverage of
the nominal 95 % percentile interval, the type I error when a component is
truly zero, and power when it is not.

The estimand is the population value of the same functional the estimator
computes, with the same convention applied at population level: a lineage with
zero share in one collection contributes to composition only. Any other choice
would be measuring the estimator against a target it does not aim at.

    python benchmarks/decomposition_calibration.py --replicates 1000 --out results/

The output is one JSON record per scenario cell plus a receipt naming the seed,
the grid and the environment. Nothing here reads or writes package state.
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from amr_clonalshare.clonality import (
    decompose_prevalence_difference,
)

# --- population shapes ------------------------------------------------------

SHARE_PROFILES = {
    # every lineage equally abundant: the friendliest case, and not the one a
    # submission-driven collection produces
    "even": lambda n_lineages: np.full(n_lineages, 1.0 / n_lineages),
    # one clone at 30 % and a Zipf tail: the shape of the shipped S. suis
    # cohort, where ST1 holds 31 % and 60 of 108 types are singletons
    "dominant": lambda n_lineages: _dominant(n_lineages),
}


def _dominant(n_lineages: int) -> np.ndarray:
    tail = 1.0 / np.arange(1, n_lineages)
    tail = 0.70 * tail / tail.sum()
    return np.r_[0.30, tail]


EFFECTS = {
    # (composition shift, within-lineage rate shift) in prevalence points
    "null": (0.00, 0.00),
    "composition_only": (0.10, 0.00),
    "within_only": (0.00, 0.10),
    "offsetting": (0.12, -0.12),
}


def build_population(n_lineages: int, profile: str, turnover: float,
                     base_prevalence: float, effect: str,
                     rng: np.random.Generator) -> dict:
    """Two populations over one lineage set, with a known decomposition.

    ``turnover`` is the fraction of lineages truly absent from one collection.
    The composition effect is created by moving share towards the lineages that
    carry the trait most; the within-lineage effect by shifting every shared
    lineage's rate. Both are then measured on the constructed populations
    rather than assumed, so the recorded truth is the truth of the population
    that was actually built.
    """
    shares = SHARE_PROFILES[profile](n_lineages)
    order = rng.permutation(n_lineages)
    shares = shares[order]

    # Within-lineage rates spread around the base prevalence on the logit
    # scale, so that no rate is pinned at a boundary.
    logit = np.log(base_prevalence / (1 - base_prevalence))
    rates_b = 1.0 / (1.0 + np.exp(-(logit + rng.normal(0.0, 0.8, n_lineages))))

    composition_shift, within_shift = EFFECTS[effect]

    # within-lineage: shift every rate by the requested amount on the
    # prevalence scale, clipped away from the boundaries
    rates_a = np.clip(rates_b + within_shift, 0.005, 0.995)

    # composition: tilt shares in A towards the higher-rate lineages, scaled so
    # that the composition component of the constructed population equals the
    # requested value rather than merely pointing the right way. The tilt is a
    # zero-sum vector, so the shares still sum to one before clipping.
    shares_b = shares.copy()
    if composition_shift != 0.0:
        weight = rates_b - rates_b.mean()
        mid = (rates_a + rates_b) / 2.0
        scale = float((weight * mid).sum())
        tilt = (composition_shift * weight / scale if abs(scale) > 1e-12
                else np.zeros_like(weight))
        shares_a = np.clip(shares + tilt, 1e-9, None)
        shares_a /= shares_a.sum()
    else:
        shares_a = shares.copy()

    # turnover: a block of lineages exists in one collection only
    n_absent = int(round(turnover * n_lineages))
    absent_from_b = np.zeros(n_lineages, dtype=bool)
    absent_from_a = np.zeros(n_lineages, dtype=bool)
    if n_absent:
        half = max(1, n_absent // 2)
        absent_from_b[:half] = True
        absent_from_a[half:n_absent] = True
        shares_a = np.where(absent_from_a, 0.0, shares_a)
        shares_b = np.where(absent_from_b, 0.0, shares_b)
        shares_a = shares_a / shares_a.sum()
        shares_b = shares_b / shares_b.sum()

    truth = population_components(shares_a, rates_a, shares_b, rates_b)
    return {"shares_a": shares_a, "rates_a": rates_a,
            "shares_b": shares_b, "rates_b": rates_b, **truth}


def population_components(wa, pa, wb, pb) -> dict:
    """The estimand: the same functional, evaluated on the population.

    The convention that defines the estimator is applied here too. A lineage
    with zero share in one collection has no rate there, and is given the
    other's, so its within-lineage term is zero and its whole contribution sits
    in composition.
    """
    present_a, present_b = wa > 0, wb > 0
    pa_f = np.where(present_a, pa, pb)
    pb_f = np.where(present_b, pb, pa)
    composition = float(((wa - wb) * (pa_f + pb_f) / 2.0).sum())
    within = float(((wa + wb) / 2.0 * (pa_f - pb_f)).sum())
    shared = present_a & present_b
    return {
        "true_composition": composition,
        "true_within_lineage": within,
        "true_difference": float((wa * pa_f).sum() - (wb * pb_f).sum()),
        "true_shared_support": float((wa[shared].sum() + wb[shared].sum()) / 2),
    }


def draw(population: dict, n_a: int, n_b: int,
         rng: np.random.Generator) -> tuple:
    """One cohort pair from the population."""
    out = []
    for side, n in (("a", n_a), ("b", n_b)):
        shares = population[f"shares_{side}"]
        rates = population[f"rates_{side}"]
        codes = rng.choice(shares.size, size=n, p=shares)
        y = (rng.random(n) < rates[codes]).astype(float)
        out.append((y, np.array([f"L{c}" for c in codes], dtype=object)))
    return out[0][0], out[0][1], out[1][0], out[1][1]


# --- one scenario cell ------------------------------------------------------

def run_cell(cell: dict) -> dict:
    rng = np.random.default_rng(cell["seed"])
    hits = {"composition": 0, "within_lineage": 0}
    excl = {"composition": 0, "within_lineage": 0}
    err = {"composition": [], "within_lineage": []}
    widths = {"composition": [], "within_lineage": []}
    supports, gated, truths = [], 0, []

    for _ in range(cell["replicates"]):
        population = build_population(
            cell["n_lineages"], cell["profile"], cell["turnover"],
            cell["base_prevalence"], cell["effect"], rng)
        ya, lina, yb, linb = draw(population, cell["n_a"], cell["n_b"], rng)
        result = decompose_prevalence_difference(
            ya, lina, yb, linb, n_boot=cell["n_boot"], rng=rng,
            min_shared_support=cell["min_shared_support"])
        if result["status"] != "ok":
            continue
        supports.append(result["shared_support_isolate_share"])
        gated += 0 if result["within_lineage_estimable"] else 1
        truths.append((population["true_composition"],
                       population["true_within_lineage"]))
        for key, true_value in (("composition", population["true_composition"]),
                                ("within_lineage",
                                 population["true_within_lineage"])):
            low, high = result[f"{key}_ci95"]
            hits[key] += int(low <= true_value <= high)
            excl[key] += int(low * high > 0)
            err[key].append(result[key] - true_value)
            widths[key].append(high - low)

    n = len(supports)
    record = {k: cell[k] for k in
              ("n_lineages", "profile", "turnover", "base_prevalence",
               "effect", "n_a", "n_b", "n_boot", "replicates",
               "min_shared_support", "seed")}
    record["n_usable"] = n
    record["mean_shared_support"] = float(np.mean(supports)) if n else None
    record["fraction_gated"] = float(gated / n) if n else None
    for key in ("composition", "within_lineage"):
        true_mean = float(np.mean([t[0 if key == "composition" else 1]
                                   for t in truths])) if n else float("nan")
        record[f"{key}_truth"] = true_mean
        record[f"{key}_coverage"] = float(hits[key] / n) if n else None
        record[f"{key}_bias"] = float(np.mean(err[key])) if n else None
        record[f"{key}_rmse"] = (float(np.sqrt(np.mean(np.square(err[key]))))
                                 if n else None)
        record[f"{key}_mean_ci_width"] = (float(np.mean(widths[key]))
                                          if n else None)
        # An interval excluding zero is a discovery. When the truth is zero it
        # is a false one, which is why the same count is reported under two
        # names and the caller reads the one the cell's truth licenses.
        record[f"{key}_reject_rate"] = float(excl[key] / n) if n else None
    return record


def build_grid(args) -> list:
    cells, seed = [], args.seed
    for n_lineages in args.lineages:
        for profile in args.profiles:
            for turnover in args.turnover:
                for base_prevalence in args.prevalence:
                    for effect in args.effects:
                        for n in args.cohort:
                            seed += 1
                            cells.append({
                                "n_lineages": n_lineages, "profile": profile,
                                "turnover": turnover,
                                "base_prevalence": base_prevalence,
                                "effect": effect, "n_a": n, "n_b": n,
                                "n_boot": args.n_boot,
                                "replicates": args.replicates,
                                "min_shared_support": args.min_shared_support,
                                "seed": seed})
    return cells


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--n-boot", type=int, default=999)
    parser.add_argument("--seed", type=int, default=20260831)
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--min-shared-support", type=float, default=0.5)
    parser.add_argument("--lineages", type=int, nargs="+", default=[10, 30, 100])
    parser.add_argument("--cohort", type=int, nargs="+",
                        default=[50, 100, 200, 500, 1000])
    parser.add_argument("--turnover", type=float, nargs="+",
                        default=[0.0, 0.2, 0.5])
    parser.add_argument("--prevalence", type=float, nargs="+",
                        default=[0.10, 0.40])
    parser.add_argument("--profiles", nargs="+",
                        default=["even", "dominant"])
    parser.add_argument("--effects", nargs="+", default=list(EFFECTS))
    parser.add_argument("--out", type=Path,
                        default=Path(__file__).resolve().parent
                        / "results_decomposition_calibration")
    args = parser.parse_args()

    cells = build_grid(args)
    args.out.mkdir(parents=True, exist_ok=True)
    started = time.time()
    print(f"{len(cells)} cells x {args.replicates} replicates "
          f"x {args.n_boot} bootstrap on {args.workers} workers", flush=True)

    records = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for done, record in enumerate(pool.map(run_cell, cells, chunksize=1), 1):
            records.append(record)
            if done % 25 == 0 or done == len(cells):
                print(f"  {done}/{len(cells)}  "
                      f"{time.time() - started:.0f}s", flush=True)

    (args.out / "cells.json").write_text(json.dumps(records, indent=1))
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat()
        .replace("+00:00", "Z"),
        "script": Path(__file__).name,
        "grid": {k: getattr(args, k) for k in
                 ("lineages", "cohort", "turnover", "prevalence", "profiles",
                  "effects", "replicates", "n_boot", "seed",
                  "min_shared_support")},
        "n_cells": len(cells),
        "n_decompositions": len(cells) * args.replicates,
        "elapsed_seconds": round(time.time() - started, 1),
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "platform": platform.platform(),
            "cpu_count": os.cpu_count(),
        },
    }
    (args.out / "RECEIPT.json").write_text(json.dumps(receipt, indent=1))
    print(f"written: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
