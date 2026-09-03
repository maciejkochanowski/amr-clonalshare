#!/usr/bin/env python3
"""Calibration of the interval-censored clonal share.

    python benchmarks/censored_calibration.py --out <dir>

Four questions, each answered by simulation against a known latent share, and
each producing a number the article quotes.

**What does coarsening cost?** One latent log2 concentration is read three
ways: exactly, as the interval a doubling panel records, and as a single
non-wild-type call. The estimand is the same in all three, so the difference
between them is the price of the reading and nothing else.

**Where does that price fall?** The call is swept across the panel so that the
prevalence above the cut-off runs from a tail to the middle. A cut in a tail
separates few isolates and is where the call carries least.

**Which interval answers which question?** The cluster bootstrap resamples
whole lineages and its width is governed by how many lineages the cohort holds;
the profile likelihood reads the information each isolate carries. Coverage of
both is measured, and so is how each width responds to the reading.

**How much of a cohort may lie beyond the panel?** A lineage every one of whose
readings is one-sided has no identified latent mean. The censoring is swept
from none to nearly total, and the estimability limit in ``censored.py`` is
read off the point where the bias starts to move.
"""
from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path

import numpy as np

from amr_clonalshare.censored import (CENSORED_GROUP_LIMIT,
                                         censored_clonal_share,
                                         intervals_from_binary)

TAU2 = 1.0
SIG2 = 1.0
TRUE = TAU2 / (TAU2 + SIG2)
WELLS = np.arange(-4.0, 6.0)
SEED = 20260901


def cohort(rng, n_lineages: int, per_lineage: int):
    code = np.repeat(np.arange(n_lineages), per_lineage)
    mu = rng.normal(0.0, np.sqrt(TAU2), n_lineages)
    z = mu[code] + rng.normal(0.0, np.sqrt(SIG2), code.size)
    return z, code, np.array(["L%d" % g for g in code], dtype=object)


def dilution(z):
    """The interval a doubling panel records for each latent value."""
    idx = np.searchsorted(WELLS, z, side="left")
    lo = np.where(idx > 0, WELLS[np.clip(idx - 1, 0, WELLS.size - 1)], -np.inf)
    hi = np.where(idx < WELLS.size, WELLS[np.clip(idx, 0, WELLS.size - 1)],
                  np.inf)
    return lo, hi


def summarise(rows: list[dict]) -> dict:
    est = np.array([r["kappa"] for r in rows], dtype=float)
    ok = np.isfinite(est)
    out = {
        "n": int(ok.sum()),
        "mean": float(est[ok].mean()),
        "bias": float(est[ok].mean() - TRUE),
        "sd": float(est[ok].std(ddof=1)),
        "rmse": float(np.sqrt(((est[ok] - TRUE) ** 2).mean())),
    }
    for key in ("interval", "bootstrap", "profile"):
        out[key + "_width"] = float(np.nanmean(
            [r[key + "_width"] for r in rows]))
        out[key + "_coverage"] = float(np.mean(
            [r[key + "_covers"] for r in rows]))
    return out


def one(lo, hi, lineage, n_boot: int) -> dict:
    r = censored_clonal_share(lo, hi, lineage, n_boot=n_boot, seed=1)

    def covers(low, high):
        return bool(np.isfinite(low) and np.isfinite(high)
                    and low <= TRUE <= high)

    return {
        "kappa": r.kappa,
        "interval_width": r.ci_high - r.ci_low,
        "interval_covers": covers(r.ci_low, r.ci_high),
        "bootstrap_width": r.boot_high - r.boot_low,
        "bootstrap_covers": covers(r.boot_low, r.boot_high),
        "profile_width": r.profile_high - r.profile_low,
        "profile_covers": covers(r.profile_low, r.profile_high),
    }


def reading_arms(reps: int, n_lineages: int, per_lineage: int,
                 n_boot: int) -> dict:
    cuts = (-1.0, 0.0, 1.0, 2.0)
    arms: dict[str, list] = {"point": [], "interval": []}
    for cut in cuts:
        arms["binary_cut_%+.0f" % cut] = []
    prevalence = {("binary_cut_%+.0f" % c): [] for c in cuts}
    rng = np.random.default_rng(SEED)
    for _ in range(reps):
        z, _code, lineage = cohort(rng, n_lineages, per_lineage)
        arms["point"].append(one(z, z, lineage, n_boot))
        lo, hi = dilution(z)
        arms["interval"].append(one(lo, hi, lineage, n_boot))
        for cut in cuts:
            key = "binary_cut_%+.0f" % cut
            y = (z > cut).astype(float)
            prevalence[key].append(float(y.mean()))
            blo, bhi = intervals_from_binary(y, cutoff_log2=cut)
            arms[key].append(one(blo, bhi, lineage, n_boot))
    out = {}
    for name, rows in arms.items():
        entry = summarise(rows)
        if name in prevalence:
            entry["prevalence"] = float(np.mean(prevalence[name]))
        out[name] = entry
    return out


def censoring_sweep(reps: int, n_lineages: int, per_lineage: int) -> list[dict]:
    rows = []
    for cut in np.arange(-1.0, 3.25, 0.25):
        rng = np.random.default_rng(SEED + 1)
        errors, shares, counts = [], [], []
        for _ in range(reps):
            z, code, lineage = cohort(rng, n_lineages, per_lineage)
            below = z <= cut
            lo = np.where(below, -np.inf, cut)
            hi = np.where(below, cut, np.inf)
            r = censored_clonal_share(lo, hi, lineage, n_boot=0,
                                      profile=False, seed=1)
            if np.isfinite(r.kappa):
                errors.append(r.kappa - TRUE)
                shares.append(r.share_in_censored_groups)
                counts.append(r.n_groups_fully_censored)
        if not errors:
            continue
        rows.append({
            "cut_log2": float(cut),
            "isolate_share_in_censored_lineages": float(np.mean(shares)),
            "lineages_wholly_censored": float(np.mean(counts)),
            "bias": float(np.mean(errors)),
            "sd": float(np.std(errors, ddof=1)),
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path,
                        default=Path("benchmarks/results_censored"))
    parser.add_argument("--reps", type=int, default=200)
    parser.add_argument("--lineages", type=int, default=30)
    parser.add_argument("--per-lineage", type=int, default=25)
    parser.add_argument("--n-boot", type=int, default=200)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    payload = {
        "generated": date.today().isoformat(),
        "design": {
            "true_share": TRUE,
            "between_lineage_variance": TAU2,
            "within_lineage_variance": SIG2,
            "replicates": args.reps,
            "lineages": args.lineages,
            "isolates_per_lineage": args.per_lineage,
            "bootstrap_draws": args.n_boot,
            "panel_wells": int(WELLS.size),
            "seed": SEED,
        },
        "readings": reading_arms(args.reps, args.lineages, args.per_lineage,
                                 args.n_boot),
        "censoring_sweep": censoring_sweep(args.reps, args.lineages,
                                           args.per_lineage),
        "estimability_limit": CENSORED_GROUP_LIMIT,
    }
    target = args.out / "censored_calibration.json"
    target.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["readings"], indent=2))
    print("written", target)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
