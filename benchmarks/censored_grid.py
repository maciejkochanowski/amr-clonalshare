#!/usr/bin/env python3
"""The three readings of a panel over the design grid.

    python benchmarks/censored_grid.py --out <dir>

One design is not an operating envelope. This driver runs the same readings as
``censored_calibration.py`` -- the exact value, the dilution a doubling panel
records, and a single cut at four positions -- over five cohort shapes, from 15
lineages of 50 isolates to 100 of 6, and four true shares including a null,
at 200 replicates a cell: 120 cells, 4,000 simulated cohorts read six ways,
24,000 readings in all. It writes one row per cell in the layout the supplement's design-grid table is read
from. The interval scored is the variance-ratio interval the package reports;
the cluster bootstrap and the likelihood-ratio width are scored only when
``--n-boot`` is positive, because they are checks and not the interval.

The latent model is the one-way random-effects model of the calibration
driver: a lineage effect with variance ``share / (1 - share)`` and a residual
of unit variance, so the true share is the design value exactly. The panel is
the ten-well doubling panel of the calibration driver, and the dilution
reading passes through :func:`amr_clonalshare.censored.intervals_from_mic`,
the function a run uses, so that the lattice the package infers from recorded
values is part of what is measured.
"""
from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path

import numpy as np

from amr_clonalshare.censored import (CENSORED_GROUP_LIMIT,
                                         SINGLE_CUT_PREVALENCE,
                                         censored_clonal_share,
                                         intervals_from_binary,
                                         intervals_from_mic)

SEED = 20260901
SHAPES = ((15, 50), (30, 8), (30, 25), (60, 12), (100, 6))
SHARES = (0.0, 0.2, 0.5, 0.8)
CUTS = (-1.0, 0.0, 1.0, 2.0)
WELLS = np.arange(-4.0, 6.0)


def cohort(rng, n_lineages: int, per_lineage: int, share: float):
    tau = np.sqrt(share / (1.0 - share)) if share < 1.0 else 0.0
    code = np.repeat(np.arange(n_lineages), per_lineage)
    mu = rng.normal(0.0, tau, n_lineages)
    z = mu[code] + rng.normal(0.0, 1.0, code.size)
    return z, np.array(["L%d" % g for g in code], dtype=object)


def recorded_mic(z):
    """The concentration a doubling panel records: the lowest tested well at
    or above the latent value, the end wells absorbing what lies beyond."""
    idx = np.clip(np.searchsorted(WELLS, z, side="left"), 0, WELLS.size - 1)
    return 2.0 ** WELLS[idx]


def one(lo, hi, lineage, truth: float, n_boot: int) -> dict:
    r = censored_clonal_share(lo, hi, lineage, n_boot=n_boot,
                              profile=n_boot > 0, seed=1)

    def covers(low, high):
        return bool(np.isfinite(low) and np.isfinite(high)
                    and low <= truth <= high)

    return {"kappa": r.kappa, "estimable": bool(r.estimable),
            "interval_width": r.ci_high - r.ci_low,
            "interval_covers": covers(r.ci_low, r.ci_high),
            "bootstrap_width": r.boot_high - r.boot_low,
            "bootstrap_covers": covers(r.boot_low, r.boot_high),
            "profile_width": r.profile_high - r.profile_low,
            "profile_covers": covers(r.profile_low, r.profile_high),
            "share_in_censored_groups": r.share_in_censored_groups}


def summarise(rows: list[dict], truth: float, prevalence: list[float]) -> dict:
    est = np.array([r["kappa"] for r in rows], dtype=float)
    ok = np.isfinite(est)
    out = {"finite": int(ok.sum()),
           "mean": float(est[ok].mean()) if ok.any() else float("nan"),
           "bias": float(est[ok].mean() - truth) if ok.any() else float("nan"),
           "sd": float(est[ok].std(ddof=1)) if ok.sum() > 1 else float("nan"),
           "rmse": float(np.sqrt(((est[ok] - truth) ** 2).mean()))
           if ok.any() else float("nan")}
    for key in ("interval", "bootstrap", "profile"):
        widths = np.array([r[key + "_width"] for r in rows], dtype=float)
        scored = np.isfinite(widths).any()
        out[key + "_width"] = float(np.nanmean(widths)) if scored else float("nan")
        out[key + "_coverage"] = (float(np.mean([r[key + "_covers"] for r in rows]))
                                  if scored else float("nan"))
    out["prevalence"] = float(np.mean(prevalence)) if prevalence else float("nan")
    out["share_in_censored_groups"] = float(np.mean(
        [r["share_in_censored_groups"] for r in rows]))
    out["estimable_fraction"] = float(np.mean([r["estimable"] for r in rows]))
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path,
                        default=Path("benchmarks/results_censored"))
    parser.add_argument("--reps", type=int, default=200)
    parser.add_argument("--n-boot", type=int, default=0)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    readings = ["point", "interval"] + ["cut%+.0f" % c for c in CUTS]
    rows = []
    cell = 0
    for (n_lineages, per_lineage) in SHAPES:
        for share in SHARES:
            rng = np.random.default_rng(SEED + cell)
            arms: dict[str, list] = {name: [] for name in readings}
            prevalence: dict[str, list] = {name: [] for name in readings}
            for _ in range(args.reps):
                z, lineage = cohort(rng, n_lineages, per_lineage, share)
                arms["point"].append(one(z, z, lineage, share, args.n_boot))
                lo, hi = intervals_from_mic(recorded_mic(z))
                arms["interval"].append(one(lo, hi, lineage, share, args.n_boot))
                for cut in CUTS:
                    name = "cut%+.0f" % cut
                    y = (z > cut).astype(float)
                    prevalence[name].append(float(y.mean()))
                    blo, bhi = intervals_from_binary(y, cutoff_log2=cut)
                    arms[name].append(one(blo, bhi, lineage, share, args.n_boot))
            for name in readings:
                entry = {"cell": cell, "lineages": n_lineages,
                         "isolates_per_lineage": per_lineage,
                         "true_share": share, "reading": name,
                         "replicates": args.reps,
                         "estimability_limit": CENSORED_GROUP_LIMIT,
                         "single_cut_prevalence_window": list(SINGLE_CUT_PREVALENCE)}
                entry.update(summarise(arms[name], share, prevalence[name]))
                rows.append(entry)
            cell += 1
            print(f"cell {cell}/{len(SHAPES) * len(SHARES)} done", flush=True)

    payload = {"generated": date.today().isoformat(), "seed": SEED,
               "cells_expected": len(rows), "cells_present": len(rows),
               "cells_missing": 0, "replicates_per_cell": args.reps,
               "simulated_cohorts": args.reps * len(SHAPES) * len(SHARES),
               "bootstrap_draws": args.n_boot, "rows": rows}
    target = args.out / "censored_grid.json"
    target.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    print("written", target)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
