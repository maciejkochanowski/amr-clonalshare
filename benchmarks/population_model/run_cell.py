#!/usr/bin/env python3
"""Run replicates of one cell of the population-model study.

    python -m benchmarks.population_model.run_cell --cell 17 --start 0 \
        --replicates 250 --method lr --output OUT/lr/cell_017_0000.csv
    python -m benchmarks.population_model.run_cell --cell 230 --start 0 \
        --replicates 25 --method general --critical-value 6.5 --output ...

``lr`` fits the profile likelihood and records the likelihood-ratio statistic
at the generating rho; with ``--critical-value`` it also records the interval
the fixed cut-off gives. ``general`` records the interval of the general
method. ``benefit`` records the profile estimate beside the observed-scale
lineage-membership share of the same isolates. One row per replicate; a
failed fit is a row with its reason, never a missing row.
"""
from __future__ import annotations

import argparse
import csv
import math
import time
from pathlib import Path

import numpy as np

from amr_clonalshare import _population_probit_numerics as numerics
from benchmarks.population_model.design import cells, draw, rng_for

FIELDS = ("cell", "family", "replicate", "method", "rho", "prevalence", "status",
          "rho_hat", "ci_low", "ci_high", "lr_at_truth", "covered", "share_kappa_adj",
          "share_latent", "failure", "seconds")


def _lr(k, m, rho, critical):
    out = numerics.fit_profile(k, m, critical_value=critical or numerics.LR95,
                               truth_rho=rho, compute_interval=critical is not None)
    row = dict(status=out["status"], rho_hat=out["rho_hat"], lr_at_truth=out["lr_at_truth"],
               ci_low=out["ci_low"] if critical is not None else None,
               ci_high=out["ci_high"] if critical is not None else None)
    if out["status"] != "ok":
        # An unidentified set is [0, 1] and carries no statistic.
        row["lr_at_truth"] = 0.0
    return row


def _general(k, m, rho, seed, key, session):
    from amr_clonalshare import general_probit_icc
    out = general_probit_icc(k, m, seed=seed, case_key=key, compute=session).as_dict()
    return dict(status=out["status"], rho_hat=out["rho_hat"], ci_low=out["ci_low"],
                ci_high=out["ci_high"], failure=out.get("failure_reason"))


def _benefit(k, m, rho, seed):
    from amr_clonalshare.attribution import clonal_share
    row = _lr(k, m, rho, None)
    y = np.concatenate([np.r_[np.ones(c), np.zeros(n - c)] for c, n in zip(k, m)])
    lineage = np.repeat(np.arange(len(m)), m)
    share = clonal_share(y, lineage, folds=5, repeats=5, n_perm=49, n_boot=0, seed=seed)
    row.update(share_kappa_adj=share.kappa_adj, share_latent=share.latent_share)
    return row


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cell", type=int, required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--replicates", type=int, required=True)
    parser.add_argument("--method", choices=("lr", "general", "benefit"), required=True)
    parser.add_argument("--critical-value", type=float)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    cell = cells()[args.cell]
    if args.output.exists():
        raise FileExistsError(args.output)
    session = None
    if args.method == "general":
        from amr_clonalshare import GeneralComputeSession
        session = GeneralComputeSession()
    rows = []
    for replicate in range(args.start, args.start + args.replicates):
        rng = rng_for(cell, replicate)
        k, m = draw(cell, rng)
        seed = int(rng.integers(0, 2**31 - 1))
        started = time.perf_counter()
        row = dict(cell=cell["cell"], family=cell["family"], replicate=replicate, method=args.method,
                   rho=cell["rho"], prevalence=cell["prevalence"])
        try:
            if args.method == "lr":
                row.update(_lr(k, m, cell["rho"], args.critical_value))
            elif args.method == "general":
                row.update(_general(k, m, cell["rho"], seed, f"cell{cell['cell']}-rep{replicate}", session))
            else:
                row.update(_benefit(k, m, cell["rho"], seed))
        except (RuntimeError, ValueError, FloatingPointError, OverflowError, ArithmeticError, MemoryError) as error:
            row.update(status="numerical_failure", failure=str(error)[:200])
        low, high = row.get("ci_low"), row.get("ci_high")
        if low is not None and high is not None and math.isfinite(low) and math.isfinite(high):
            row["covered"] = int(low - 1e-12 <= cell["rho"] <= high + 1e-12)
        row["seconds"] = round(time.perf_counter() - started, 3)
        rows.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output.with_suffix(".partial")
    with partial.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
