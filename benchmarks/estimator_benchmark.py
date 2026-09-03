#!/usr/bin/env python3
"""The clonal share against the estimators a reader would reach for instead.

    python estimator_benchmark.py --cell <i> --out <dir>
    python estimator_benchmark.py --aggregate --out <dir>
    python estimator_benchmark.py --real <cohort.json> --out <dir>

WHY THIS EXISTS. "How much of a resistance phenotype sits between lineages" is
not a new question, and a reader who has not met this package will ask why the
answer should not come from a variance component. Four alternatives are run
here on identical data, with the truth known, so the comparison is a
measurement rather than an argument.

THE ESTIMATORS, and what each targets.

  in-sample R2 of the lineage means. The quantity a reader computes first. It
  is the fraction of the sum of squares that falls between lineages, and it is
  biased upwards by exactly the number of lineages divided by the number of
  isolates, so on a cohort typed at high resolution it approaches one whatever
  the trait is. It is included because it is what a spreadsheet produces, and
  because the size of its bias is the reason the rest of this exists.

  the analysis-of-variance variance component, the classical one-way estimator
  of tau2 / (tau2 + sigma2) formed from the two mean squares (Searle, Casella
  and McCulloch). Unbiased for the ratio under the model, and the interval
  comes from inverting the central F.

  restricted maximum likelihood in the same model. This is what a linear mixed
  model fits when a lineage random intercept is the only random term, and it
  is the estimator the bacterial heritability literature uses under the name
  LMM. Mai and colleagues report that the LMM over-estimates heritability
  markedly in the middle of the range (Bioinformatics Advances 2023;3:vbad027,
  doi:10.1093/bioadv/vbad027), and that claim is testable here: the same grid
  is run and the bias is measured rather than assumed.

  clonal_share, this package: out-of-sample skill of the lineage-mean
  predictor, debiased on a permuted-label run.

  realised_share, this package: the share carried by the lineages the cohort
  holds, from inverting the noncentral F.

WHAT WOULD MEAN THIS MEASURED THE WRONG THING. The truth is computed from the
lineage effects drawn in each replicate, not from the sample that produced the
estimate. If every estimator returned the same value the design would have no
leverage and the comparison would be empty; if the in-sample statistic were
not the most biased of the five at high lineage counts, the generator would not
be producing the regime the comparison is about. Both are checked and reported
in the verdict.
"""
from __future__ import annotations

import argparse
import json
import platform
import time
from datetime import date
from pathlib import Path

import numpy as np
import scipy
from scipy import optimize, stats

from amr_clonalshare.attribution import clonal_share
from amr_clonalshare.realised import realised_share

SEED = 20260902
ALPHA = 0.05


# ---------------------------------------------------------------- generators
def cohort(rng, n_groups, size, share, binary=False, prevalence=0.25,
           unbalance="balanced"):
    """A cohort with a known between-lineage share, and that share returned."""
    if unbalance == "balanced":
        sizes = np.full(n_groups, size, dtype=int)
    else:
        sizes = 1 + rng.poisson(max(size - 1, 1), n_groups)
    tau = np.sqrt(share / (1.0 - share)) if share < 1.0 else 0.0
    effects = rng.normal(0.0, tau, n_groups)
    lineage = np.repeat(np.arange(n_groups), sizes)
    if binary:
        # On the observed scale the truth is the intraclass correlation of a
        # Bernoulli whose probability varies by lineage, which is a different
        # number from the latent share and is computed here rather than assumed.
        p = np.clip(stats.norm.cdf(stats.norm.ppf(prevalence) + effects), 1e-6,
                    1 - 1e-6)
        y = rng.binomial(1, np.repeat(p, sizes)).astype(float)
        w = sizes / sizes.sum()
        pbar = float((w * p).sum())
        between = float((w * (p - pbar) ** 2).sum())
        within = float((w * p * (1 - p)).sum())
        truth = between / (between + within)
        # On the observed scale the two questions coincide only in the limit;
        # the superpopulation target for a binary trait is not the latent
        # share either, so it is reported as the same observed-scale quantity
        # and the binary arm is read as a robustness check rather than as a
        # coverage statement about two estimands.
        return y, lineage, {"realised": truth, "superpopulation": truth}
    y = np.repeat(effects, sizes) + rng.normal(0.0, 1.0, sizes.sum())
    realised = ((effects - effects.mean()) ** 2).sum() / (n_groups - 1)
    return y, lineage, {"realised": realised / (realised + 1.0),
                        "superpopulation": share}


# ---------------------------------------------------------------- estimators
def _mean_squares(y, lineage):
    codes, code = np.unique(lineage, return_inverse=True)
    g = codes.size
    n = y.size
    counts = np.bincount(code, minlength=g).astype(float)
    sums = np.bincount(code, weights=y, minlength=g)
    means = sums / np.maximum(counts, 1)
    grand = y.mean()
    ssb = float((counts * (means - grand) ** 2).sum())
    ssw = float(((y - means[code]) ** 2).sum())
    n0 = (n - (counts ** 2).sum() / n) / (g - 1)
    return ssb, ssw, g, n, float(n0)


def in_sample_r2(y, lineage):
    ssb, ssw, g, n, n0 = _mean_squares(y, lineage)
    total = ssb + ssw
    return {"point": ssb / total if total > 0 else float("nan"),
            "low": float("nan"), "high": float("nan")}


def anova_component(y, lineage):
    """The classical one-way variance-component share, with an F interval."""
    ssb, ssw, g, n, n0 = _mean_squares(y, lineage)
    df1, df2 = g - 1, n - g
    if df1 < 1 or df2 < 1 or ssw <= 0:
        return {"point": float("nan"), "low": float("nan"),
                "high": float("nan")}
    msb, msw = ssb / df1, ssw / df2
    tau2 = max((msb - msw) / n0, 0.0)
    point = tau2 / (tau2 + msw) if (tau2 + msw) > 0 else float("nan")
    f = msb / msw
    lo_r = (f / stats.f.ppf(1 - ALPHA / 2, df1, df2) - 1.0) / n0
    hi_r = (f / stats.f.ppf(ALPHA / 2, df1, df2) - 1.0) / n0
    to_share = lambda r: max(r, 0.0) / (max(r, 0.0) + 1.0)
    return {"point": point, "low": to_share(lo_r), "high": to_share(hi_r)}


def reml_component(y, lineage):
    """REML in the one-way random-effects model, profiled on the ratio.

    A linear mixed model with a lineage random intercept and nothing else has a
    restricted likelihood that depends on the two variances only through their
    ratio, so the fit is a one-dimensional maximisation and needs no optimiser
    library. The interval is the profile-likelihood one at a chi-square
    cut-off on one degree of freedom, which is what a mixed-model package
    reports.
    """
    codes, code = np.unique(lineage, return_inverse=True)
    g, n = codes.size, y.size
    counts = np.bincount(code, minlength=g).astype(float)
    sums = np.bincount(code, weights=y, minlength=g)
    if g < 2 or n <= g:
        return {"point": float("nan"), "low": float("nan"),
                "high": float("nan")}

    def neg_restricted(log_ratio):
        ratio = float(np.exp(log_ratio))
        w = 1.0 / (1.0 + ratio * counts)          # weight of each lineage mean
        mu = float((w * sums / counts * counts).sum() / (w * counts).sum())
        # Residual sum of squares under this ratio, profiled over sigma2.
        within = float(((y - (sums / counts)[code]) ** 2).sum())
        between = float((counts * w * (sums / counts - mu) ** 2).sum())
        quad = within + between
        if quad <= 0:
            return np.inf
        logdet = float(np.sum(np.log1p(ratio * counts)))
        logdet_x = float(np.log((w * counts).sum()))
        return 0.5 * ((n - 1) * np.log(quad) + logdet + logdet_x)

    grid = np.linspace(-12.0, 6.0, 181)
    values = np.array([neg_restricted(v) for v in grid])
    start = grid[int(np.argmin(values))]
    fit = optimize.minimize_scalar(neg_restricted, bracket=None,
                                   bounds=(start - 0.5, start + 0.5),
                                   method="bounded")
    best, best_value = float(fit.x), float(fit.fun)
    ratio = float(np.exp(best))
    point = ratio / (1.0 + ratio)
    cut = best_value + 0.5 * stats.chi2.ppf(1 - ALPHA, 1)
    low_r, high_r = 0.0, np.inf
    below = grid[values <= cut]
    if below.size:
        low_r = float(np.exp(below.min()))
        high_r = float(np.exp(below.max()))
        if below.min() <= grid[0] + 1e-9:
            low_r = 0.0
    to_share = lambda r: r / (1.0 + r) if np.isfinite(r) else 1.0
    return {"point": point, "low": to_share(low_r), "high": to_share(high_r)}


def package_clonal(y, lineage):
    r = clonal_share(y, lineage, folds=5, repeats=10, n_boot=200, n_perm=100,
                     seed=SEED).as_dict()
    return {"point": r["kappa_adj"], "low": r["ci_low"], "high": r["ci_high"]}


def package_realised(y, lineage):
    r = realised_share(y, lineage)
    return {"point": r.kappa, "low": r.ci_low, "high": r.ci_high,
            "estimable": r.estimable}


ESTIMATORS = {
    "in_sample_r2": in_sample_r2,
    "anova_component": anova_component,
    "reml_lmm": reml_component,
    "clonal_share": package_clonal,
    "realised_share": package_realised,
}

#: Which quantity each estimator is aiming at. Scoring an estimator against a
#: target it does not aim at is the mistake this table exists to prevent: the
#: variance component and the mixed model answer what a fresh draw of lineages
#: would show, and only realised_share answers what the lineages in hand do.
#: Both scores are recorded for every estimator so the difference is visible.
TARGET = {
    "in_sample_r2": "superpopulation",
    "anova_component": "superpopulation",
    "reml_lmm": "superpopulation",
    "clonal_share": "superpopulation",
    "realised_share": "realised",
}


# --------------------------------------------------------------- the design
def design() -> list[dict]:
    cells, i = [], 0
    for binary in (False, True):
        for n_groups in (10, 30, 100):
            for size in (5, 20):
                for share in (0.0, 0.1, 0.3, 0.5, 0.7):
                    cells.append({"index": i, "binary": binary,
                                  "n_groups": n_groups, "group_size": size,
                                  "share": share, "unbalance": "balanced",
                                  "replicates": 400})
                    i += 1
    for share in (0.1, 0.3, 0.5):
        cells.append({"index": i, "binary": False, "n_groups": 30,
                      "group_size": 20, "share": share,
                      "unbalance": "poisson", "replicates": 400})
        i += 1
    return cells


def run_cell(cell: dict) -> dict:
    rng = np.random.default_rng(SEED + 977 * cell["index"])
    stats_by = {k: {"err": [], "cover": 0, "cover_other": 0, "width": [],
                    "seconds": 0.0, "estimable": 0}
                for k in ESTIMATORS}
    for _ in range(cell["replicates"]):
        y, lineage, truth = cohort(rng, cell["n_groups"], cell["group_size"],
                                   cell["share"], binary=cell["binary"],
                                   unbalance=cell["unbalance"])
        for name, fn in ESTIMATORS.items():
            t0 = time.perf_counter()
            out = fn(y, lineage)
            stats_by[name]["seconds"] += time.perf_counter() - t0
            own = truth[TARGET[name]]
            other = truth["realised" if TARGET[name] == "superpopulation"
                          else "superpopulation"]
            point, lo, hi = out["point"], out["low"], out["high"]
            if np.isfinite(point):
                stats_by[name]["err"].append(point - own)
            if np.isfinite(lo) and np.isfinite(hi):
                stats_by[name]["cover"] += int(lo <= own <= hi)
                stats_by[name]["cover_other"] += int(lo <= other <= hi)
                stats_by[name]["width"].append(hi - lo)
            stats_by[name]["estimable"] += int(out.get("estimable", True))
    result = dict(cell)
    reps = cell["replicates"]
    for name, s in stats_by.items():
        err = np.array(s["err"], dtype=float)
        width = np.array(s["width"], dtype=float)
        result[name] = {
            "bias": float(err.mean()) if err.size else None,
            "rmse": float(np.sqrt((err ** 2).mean())) if err.size else None,
            "target": TARGET[name],
            "coverage": s["cover"] / reps if width.size else None,
            "coverage_of_the_other_target": (s["cover_other"] / reps
                                             if width.size else None),
            "median_width": float(np.median(width)) if width.size else None,
            "seconds_per_call": s["seconds"] / reps,
            "estimable_share": s["estimable"] / reps,
            "n_scored": int(err.size),
        }
    return result


def aggregate(out: Path) -> dict:
    rows = [json.loads(p.read_text())
            for p in sorted((out / "cells").glob("*.json"))]
    verdict = {}
    gaussian = [r for r in rows if not r["binary"] and r["share"] > 0]
    if gaussian:
        naive = float(np.mean([r["in_sample_r2"]["bias"] for r in gaussian]))
        reml = float(np.mean([r["reml_lmm"]["bias"] for r in gaussian]))
        ours = float(np.mean([r["clonal_share"]["bias"] for r in gaussian]))
        verdict["mean_bias_in_sample_r2"] = naive
        verdict["mean_bias_reml_lmm"] = reml
        verdict["mean_bias_clonal_share"] = ours
        verdict["in_sample_is_the_most_biased"] = naive > max(abs(reml),
                                                              abs(ours))
        spread = [abs(r["in_sample_r2"]["bias"] - r["clonal_share"]["bias"])
                  for r in gaussian]
        verdict["design_has_leverage"] = bool(max(spread) > 0.05)
    payload = {"provenance": {"generated": date.today().isoformat(),
                              "seed": SEED, "alpha": ALPHA,
                              "python": platform.python_version(),
                              "numpy": np.__version__,
                              "scipy": scipy.__version__},
               "verdict": verdict, "cells": rows}
    (out / "estimator_benchmark.json").write_text(
        json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(verdict, indent=2))
    print(f"cells: {len(rows)}")
    return payload


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--cell", type=int)
    ap.add_argument("--plan", action="store_true")
    ap.add_argument("--aggregate", action="store_true")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.plan:
        cells = design()
        print(f"cells: {len(cells)}")
        return 0
    if args.aggregate:
        aggregate(args.out)
        return 0
    cells = design()
    cell = cells[args.cell]
    print(f"[{args.cell}] binary={cell['binary']} G={cell['n_groups']} "
          f"m={cell['group_size']} share={cell['share']} "
          f"unbalance={cell['unbalance']}", flush=True)
    result = run_cell(cell)
    (args.out / "cells").mkdir(exist_ok=True)
    (args.out / "cells" / f"{args.cell:03d}.json").write_text(
        json.dumps(result, indent=1) + "\n", encoding="utf-8")
    for name in ESTIMATORS:
        r = result[name]
        print(f"  {name:<18} bias={r['bias']:+.4f} rmse={r['rmse']:.4f} "
              f"cover={r['coverage'] if r['coverage'] is None else round(r['coverage'],4)} "
              f"width={r['median_width'] if r['median_width'] is None else round(r['median_width'],4)} "
              f"{r['seconds_per_call']*1000:.2f} ms")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
