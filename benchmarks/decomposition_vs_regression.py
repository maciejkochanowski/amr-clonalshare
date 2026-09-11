#!/usr/bin/env python3
"""The decomposition against the regression a reader will propose instead.

The first question any statistician asks about the Kitagawa split is why it is
not a regression. Adjust for lineage and read the coefficient on collection:
that is the within-lineage rate change, and it comes with a standard error, a
model, and a century of familiarity. The question deserves a measurement rather
than an argument, so this script runs three estimators of the same quantity on
the same simulated cohorts:

* **unadjusted** - the difference in prevalence, which is what surveillance
  reports and which cannot separate the two mechanisms at all;
* **fixed-effect logistic** - one intercept per lineage plus a collection
  indicator, converted to the prevalence scale by marginal standardisation
  (g-computation) so that it estimates the same functional as the Kitagawa
  within-lineage component rather than a log-odds ratio;
* **random-intercept logistic** - the same model with the lineage effects drawn
  from a normal, fitted by Gauss-Hermite quadrature, standardised the same way.

With equal cohort sizes the pooled empirical lineage distribution is exactly
``(w_A + w_B) / 2``, so all three prevalence-scale estimators target the same
estimand and the comparison is fair rather than rhetorical.

Two things are measured that decide the question.

**Separation.** A sequence type seen a handful of times, all wild type, drives
its own intercept to minus infinity. The count of lineages in that state is
reported, because it is the reason the maximum-likelihood fit is not available
on a cohort of this shape and why a penalty has to be introduced - and a
penalty is a choice that moves the estimate.

**Extrapolation.** A lineage observed in only one collection carries no
information about a within-lineage change, and the Kitagawa convention gives
it a within-lineage term of exactly zero. The regression does not: it fits that
lineage an intercept and then applies the shared collection effect to it,
predicting a rate change for a lineage that was never seen under both
conditions. The size of that extrapolation is reported as the gap between
standardising over all lineages and over the shared ones only.

    python benchmarks/decomposition_vs_regression.py --replicates 300
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
from scipy.optimize import minimize
from scipy.special import expit

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from amr_clonalshare.clonality import (
    decompose_prevalence_difference,
)
from decomposition_calibration import build_population, draw

_GH_NODES = 20


def _neg_loglik_fe(theta, codes, y, side, n_lineages, ridge):
    alpha, beta = theta[:n_lineages], theta[n_lineages]
    eta = alpha[codes] + beta * side
    # log(1 + exp(eta)) computed stably
    log1pexp = np.logaddexp(0.0, eta)
    nll = float(np.sum(log1pexp - y * eta))
    grad_eta = expit(eta) - y
    grad = np.zeros_like(theta)
    np.add.at(grad, codes, grad_eta)
    grad[n_lineages] = float(np.sum(grad_eta * side))
    if ridge > 0:
        nll += 0.5 * ridge * float(np.sum(alpha ** 2))
        grad[:n_lineages] += ridge * alpha
    return nll, grad


def fit_fixed_effect(codes, y, side, n_lineages, ridge=1e-3):
    """Logistic regression with a free intercept per lineage."""
    theta0 = np.zeros(n_lineages + 1)
    out = minimize(_neg_loglik_fe, theta0, jac=True, method="L-BFGS-B",
                   args=(codes, y, side, n_lineages, ridge),
                   options={"maxiter": 500})
    return out.x[:n_lineages], float(out.x[n_lineages]), bool(out.success)


def _neg_loglik_ri(theta, codes, y, side, n_lineages):
    mu, beta, log_sigma = theta
    sigma = np.exp(log_sigma)
    nodes, weights = np.polynomial.hermite_e.hermegauss(_GH_NODES)
    weights = weights / np.sqrt(2.0 * np.pi)
    fixed = mu + beta * side
    # per lineage: log sum_q w_q exp( sum_i loglik_i(node_q) )
    total = 0.0
    for lineage in range(n_lineages):
        mask = codes == lineage
        if not mask.any():
            continue
        eta = fixed[mask][:, None] + sigma * nodes[None, :]
        ll = np.sum(y[mask][:, None] * eta - np.logaddexp(0.0, eta), axis=0)
        total += float(np.log(np.sum(weights * np.exp(ll - ll.max()))) + ll.max())
    return -total


def fit_random_intercept(codes, y, side, n_lineages):
    """Random-intercept logistic fitted by Gauss-Hermite marginal likelihood."""
    theta0 = np.array([0.0, 0.0, np.log(0.5)])
    out = minimize(_neg_loglik_ri, theta0, method="Nelder-Mead",
                   args=(codes, y, side, n_lineages),
                   options={"maxiter": 2000, "xatol": 1e-4, "fatol": 1e-4})
    mu, beta, log_sigma = out.x
    sigma = float(np.exp(log_sigma))
    nodes, weights = np.polynomial.hermite_e.hermegauss(_GH_NODES)
    weights = weights / np.sqrt(2.0 * np.pi)
    fixed = mu + beta * side
    posterior = np.zeros(n_lineages)
    for lineage in range(n_lineages):
        mask = codes == lineage
        if not mask.any():
            continue
        eta = fixed[mask][:, None] + sigma * nodes[None, :]
        ll = np.sum(y[mask][:, None] * eta - np.logaddexp(0.0, eta), axis=0)
        post = weights * np.exp(ll - ll.max())
        posterior[lineage] = float(np.sum(post * sigma * nodes) / np.sum(post))
    return float(mu), float(beta), posterior, bool(out.success)


def standardise(intercepts, beta, codes, keep_mask):
    """Marginal standardisation of a fitted model to the pooled lineage mix.

    The contrast is the average predicted prevalence under the later
    collection minus the average under the earlier one, taken over the isolates
    in ``keep_mask``. With equal cohort sizes that average is over the pooled
    lineage distribution, which is the weight the Kitagawa within-lineage
    component uses.
    """
    if not keep_mask.any():
        return float("nan")
    base = intercepts[codes[keep_mask]]
    return float(np.mean(expit(base + beta) - expit(base)))


def run_cell(cell: dict) -> dict:
    rng = np.random.default_rng(cell["seed"])
    rows = {name: [] for name in
            ("kitagawa", "fe_all", "fe_shared", "ri_all", "ri_shared",
             "unadjusted")}
    truths, separated, fe_failed, ri_failed, supports = [], [], 0, 0, []

    for _ in range(cell["replicates"]):
        population = build_population(
            cell["n_lineages"], cell["profile"], cell["turnover"],
            cell["base_prevalence"], cell["effect"], rng)
        ya, lina, yb, linb = draw(population, cell["n"], cell["n"], rng)

        labels = sorted(set(lina.tolist()) | set(linb.tolist()))
        index = {v: i for i, v in enumerate(labels)}
        n_lineages = len(labels)
        codes = np.fromiter((index[v] for v in np.r_[lina, linb]),
                            dtype=np.intp, count=lina.size + linb.size)
        y = np.r_[ya, yb]
        side = np.r_[np.ones(ya.size), np.zeros(yb.size)]

        present_a = np.zeros(n_lineages, dtype=bool)
        present_b = np.zeros(n_lineages, dtype=bool)
        present_a[np.unique(codes[:ya.size])] = True
        present_b[np.unique(codes[ya.size:])] = True
        shared = present_a & present_b
        shared_mask = shared[codes]
        supports.append(float(shared_mask.mean()))

        # a lineage whose outcomes are all equal separates its own intercept
        totals = np.bincount(codes, minlength=n_lineages)
        ones = np.bincount(codes, weights=y, minlength=n_lineages)
        separated.append(int(np.sum((ones == 0) | (ones == totals))))

        result = decompose_prevalence_difference(
            ya, lina, yb, linb, n_boot=0, rng=rng, min_shared_support=0.0)
        rows["kitagawa"].append(result["within_lineage"])
        rows["unadjusted"].append(result["difference"])

        alpha, beta, ok = fit_fixed_effect(codes, y, side, n_lineages,
                                           ridge=cell["ridge"])
        fe_failed += int(not ok)
        rows["fe_all"].append(standardise(alpha, beta, codes,
                                          np.ones_like(codes, dtype=bool)))
        rows["fe_shared"].append(standardise(alpha, beta, codes, shared_mask))

        mu, beta_r, posterior, ok_r = fit_random_intercept(
            codes, y, side, n_lineages)
        ri_failed += int(not ok_r)
        rows["ri_all"].append(standardise(mu + posterior, beta_r, codes,
                                          np.ones_like(codes, dtype=bool)))
        rows["ri_shared"].append(standardise(mu + posterior, beta_r, codes,
                                             shared_mask))
        truths.append(population["true_within_lineage"])

    truth = float(np.mean(truths))
    record = {k: cell[k] for k in
              ("n_lineages", "profile", "turnover", "base_prevalence",
               "effect", "n", "replicates", "ridge", "seed")}
    record.update({
        "true_within_lineage": truth,
        "mean_shared_support": float(np.mean(supports)),
        "mean_separated_lineages": float(np.mean(separated)),
        "mean_lineages": float(cell["n_lineages"]),
        "fixed_effect_fit_failures": fe_failed,
        "random_intercept_fit_failures": ri_failed,
    })
    for name, values in rows.items():
        arr = np.asarray(values, dtype=float)
        good = arr[np.isfinite(arr)]
        record[f"{name}_bias"] = (float(np.mean(good) - truth) if good.size
                                  else None)
        record[f"{name}_rmse"] = (float(np.sqrt(np.mean((good - truth) ** 2)))
                                  if good.size else None)
        record[f"{name}_mean"] = float(np.mean(good)) if good.size else None
    record["extrapolation_gap_fixed_effect"] = (
        record["fe_all_mean"] - record["fe_shared_mean"]
        if record["fe_all_mean"] is not None else None)
    record["extrapolation_gap_random_intercept"] = (
        record["ri_all_mean"] - record["ri_shared_mean"]
        if record["ri_all_mean"] is not None else None)
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--replicates", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20260901)
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--ridge", type=float, default=1e-3)
    parser.add_argument("--lineages", type=int, nargs="+", default=[10, 30, 100])
    parser.add_argument("--cohort", type=int, nargs="+", default=[100, 300])
    parser.add_argument("--turnover", type=float, nargs="+", default=[0.0, 0.4])
    parser.add_argument("--prevalence", type=float, nargs="+", default=[0.40])
    parser.add_argument("--profiles", nargs="+", default=["dominant"])
    parser.add_argument("--effects", nargs="+",
                        default=["null", "within_only", "offsetting"])
    parser.add_argument("--out", type=Path,
                        default=Path(__file__).resolve().parent
                        / "results_decomposition_vs_regression")
    args = parser.parse_args()

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
                                "effect": effect, "n": n,
                                "replicates": args.replicates,
                                "ridge": args.ridge, "seed": seed})

    args.out.mkdir(parents=True, exist_ok=True)
    started = time.time()
    print(f"{len(cells)} cells x {args.replicates} replicates "
          f"on {args.workers} workers", flush=True)
    records = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for done, record in enumerate(pool.map(run_cell, cells, chunksize=1), 1):
            records.append(record)
            if done % 10 == 0 or done == len(cells):
                print(f"  {done}/{len(cells)}  {time.time() - started:.0f}s",
                      flush=True)

    (args.out / "cells.json").write_text(json.dumps(records, indent=1))
    (args.out / "RECEIPT.json").write_text(json.dumps({
        "generated_utc": datetime.now(timezone.utc).isoformat()
        .replace("+00:00", "Z"),
        "script": Path(__file__).name,
        "grid": {k: getattr(args, k) for k in
                 ("lineages", "cohort", "turnover", "prevalence", "profiles",
                  "effects", "replicates", "ridge", "seed")},
        "n_cells": len(cells),
        "elapsed_seconds": round(time.time() - started, 1),
        "gauss_hermite_nodes": _GH_NODES,
        "environment": {"python": platform.python_version(),
                        "numpy": np.__version__,
                        "platform": platform.platform()},
    }, indent=1))
    print(f"written: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
