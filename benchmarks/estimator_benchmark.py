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
  predictor, debiased on a permuted-label run, with its lineage-bootstrap
  interval. Its target is eta, the share of the collection's variance that
  lineage membership carries, B_w / (B_w + sigma^2) with B_w the
  isolate-weighted variance of the lineage means; and the same estimate with
  its superpopulation interval, a within-lineage bootstrap under a chi-square
  layer for the draw of the lineages, scored against the fresh-draw target.

  realised_share, this package: the design-corrected component ratio rho_c of
  equation (1) of the article, S_a^2 / (S_a^2 + sigma^2) with S_a^2 = B_w /
  (1 - sum w_g^2), from inverting the noncentral F. eta and rho_c are two
  functionals of the same lineage means and differ by (G - 1) / G on a
  balanced design; each estimator is scored against its own.

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
from amr_clonalshare.latent import observed_share
from amr_clonalshare.realised import realised_share

SEED = 20260902
ALPHA = 0.05


# ---------------------------------------------------------------- generators
#: Share of lineages that carry the effect under the two-point law: one
#: lineage in ten is a carrier, the rest sit at a common low level, the
#: pattern a resistant clone makes inside a collection.
TWO_POINT_CARRIERS = 0.1


def _binary_intercept(prevalence: float, tau: float, law: str) -> float:
    """The liability intercept at which the prevalence marginal over the
    lineage effects equals ``prevalence``."""
    if law == "normal":
        return float(stats.norm.ppf(prevalence) * np.sqrt(1.0 + tau * tau))
    pi = TWO_POINT_CARRIERS
    delta = tau * np.sqrt((1.0 - pi) / pi)
    levels = np.array([delta, -delta * pi / (1.0 - pi)])
    weights = np.array([pi, 1.0 - pi])

    def marginal(c):
        return float((weights * stats.norm.cdf(c + levels)).sum())

    lo, hi = -20.0, 20.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if marginal(mid) < prevalence:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def cohort(rng, n_groups, size, share, binary=False, prevalence=0.25,
           unbalance="balanced", law="normal"):
    """A cohort with a known between-lineage share, and that share returned.

    ``law`` is the distribution of the lineage effects. ``normal`` is the law
    the chi-square layer of the species interval is derived under;
    ``two_point`` is a lineage effect that is either high, in a tenth of the
    lineages on average, or at one common low level, with the same variance
    ``tau^2`` as the normal law. Under it the realised share of a draw can sit
    far from the superpopulation share, because a cohort of ten lineages
    holds one carrier on average and sometimes none, which is the case the
    envelope on the species interval exists for.
    """
    if unbalance == "balanced":
        sizes = np.full(n_groups, size, dtype=int)
    else:
        sizes = 1 + rng.poisson(max(size - 1, 1), n_groups)
    tau = np.sqrt(share / (1.0 - share)) if share < 1.0 else 0.0
    if law == "normal":
        effects = rng.normal(0.0, tau, n_groups)
    elif law == "two_point":
        pi = TWO_POINT_CARRIERS
        delta = tau * np.sqrt((1.0 - pi) / pi)
        carrier = rng.random(n_groups) < pi
        effects = np.where(carrier, delta, -delta * pi / (1.0 - pi))
    else:
        raise ValueError(f"unknown lineage law {law!r}")
    lineage = np.repeat(np.arange(n_groups), sizes)
    w = sizes / sizes.sum()
    design = 1.0 - float((w ** 2).sum())   # n0 (G - 1) / n
    if binary:
        # The intercept is set so that the prevalence marginal over the lineage
        # effects equals ``prevalence``: under the normal law the liability has
        # variance 1 + tau^2, so the intercept is Phi^-1(p) sqrt(1 + tau^2);
        # under the two-point law it solves the mixture equation numerically.
        # Placing Phi^-1(p) on the liability directly, as the first release
        # did, gave a marginal prevalence that rose with the share, 0.32
        # rather than 0.25 at a share of 0.5.
        c = _binary_intercept(prevalence, tau, law)
        p = np.clip(stats.norm.cdf(c + effects), 1e-6, 1 - 1e-6)
        y = rng.binomial(1, np.repeat(p, sizes)).astype(float)
        pbar = float((w * p).sum())
        between = float((w * (p - pbar) ** 2).sum())
        within = float((w * p * (1 - p)).sum())
        # The superpopulation share on the observed scale: the correlation
        # between two isolates of one lineage under the threshold model,
        # integrated over the law of the lineage effects. Under the normal
        # law that is the tetrachoric quantity ``observed_share`` computes;
        # under the two-point law it is the mixture of the two levels.
        if law == "two_point":
            pi = TWO_POINT_CARRIERS
            p_hi = stats.norm.cdf(c + delta)
            p_lo = stats.norm.cdf(c - delta * pi / (1.0 - pi))
            var_p = pi * (1 - pi) * (p_hi - p_lo) ** 2
            within_p = pi * p_hi * (1 - p_hi) + (1 - pi) * p_lo * (1 - p_lo)
            superpopulation = float(var_p / (var_p + within_p))
        else:
            superpopulation = float(observed_share(share, prevalence))
    else:
        y = np.repeat(effects, sizes) + rng.normal(0.0, 1.0, sizes.sum())
        a_w = float((w * effects).sum())
        between = float((w * (effects - a_w) ** 2).sum())
        within = 1.0
        superpopulation = share
    # Three quantities, three targets (Section 2.2 of the article). ``eta`` is
    # the share of this collection's variance that lineage membership carries,
    # B_w / (B_w + sigma^2) with B_w the isolate-weighted variance of the
    # lineage means: the quantity a lineage-mean predictor scored out of
    # sample converges to, and the target of ``clonal_share``. ``rho`` is the
    # design-corrected component ratio of equation (1), S_a^2 / (S_a^2 +
    # sigma^2) with S_a^2 = B_w / (1 - sum w_g^2), the noncentrality the
    # noncentral-F inversion of ``realised_share`` is exact for. The two
    # differ by the factor 1 - sum w_g^2, i.e. (G - 1) / G when balanced, so
    # they part company where lineages are few; scoring one estimator against
    # the other's target is a definitional error, not a property of either.
    s2 = between / design if design > 0 else float("nan")
    return y, lineage, {"eta": between / (between + within),
                        "rho": s2 / (s2 + within),
                        "superpopulation": superpopulation}


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


_last_clonal = {}


def _clonal(y, lineage):
    """One call of clonal_share per cohort, shared by its two interval rows."""
    key = (y.tobytes(), np.asarray(lineage).tobytes())
    if _last_clonal.get("key") != key:
        _last_clonal["key"] = key
        _last_clonal["r"] = clonal_share(y, lineage, folds=5, repeats=10,
                                         n_boot=200, n_perm=100,
                                         seed=SEED).as_dict()
    return _last_clonal["r"]


def package_clonal(y, lineage):
    r = _clonal(y, lineage)
    return {"point": r["kappa_adj"], "low": r["ci_low"], "high": r["ci_high"]}


def package_clonal_superpopulation(y, lineage):
    """The same estimate with the interval for a fresh draw of lineages."""
    r = _clonal(y, lineage)
    return {"point": r["kappa_adj"], "low": r["superpopulation_low"],
            "high": r["superpopulation_high"]}


def package_realised(y, lineage):
    r = realised_share(y, lineage)
    return {"point": r.kappa, "low": r.ci_low, "high": r.ci_high,
            "estimable": r.estimable}


ESTIMATORS = {
    "in_sample_r2": in_sample_r2,
    "anova_component": anova_component,
    "reml_lmm": reml_component,
    "clonal_share": package_clonal,
    "clonal_share_superpopulation": package_clonal_superpopulation,
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
    "clonal_share": "eta",
    "clonal_share_superpopulation": "superpopulation",
    "realised_share": "rho",
}
#: The "other" score of each estimator: the fresh-draw target for the two
#: collection-in-hand estimators, and the design-corrected collection target
#: for the three that estimate the fresh-draw share.
OTHER = {name: ("superpopulation" if t != "superpopulation" else "rho")
         for name, t in TARGET.items()}
#: clonal_share scores the lineage means of the lineages in hand on held-out
#: isolates, so the quantity it estimates is the share those lineages carry;
#: its bootstrap draws those lineages whole and does not carry the extra
#: uncertainty of a fresh draw of lineages, which is why on the first run of
#: this grid it covered the superpopulation target at 0.708 to 0.833 with ten
#: lineages and the realised target at 0.90 to 0.99. Both scores are still
#: recorded, so the reader can see the difference rather than take it on trust.


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
    # The ends of the lineage-count range: a collection of five lineages, and
    # the many-lineages-few-isolates regime a genome-cluster label produces
    # (the cross-species atlas holds 22,671 isolates in 5,904 clusters).
    for binary in (False, True):
        for n_groups, size in ((5, 5), (5, 20), (300, 3), (1000, 3)):
            for share in (0.0, 0.1, 0.3, 0.5, 0.7):
                cells.append({"index": i, "binary": binary,
                              "n_groups": n_groups, "group_size": size,
                              "share": share, "unbalance": "balanced",
                              "replicates": 400})
                i += 1
    # Lineage effects that are not normal: a carrier law, the shape a
    # resistant clone makes, against which the chi-square layer of the
    # species interval was not derived.
    for binary in (False, True):
        for n_groups in (10, 30, 100):
            for share in (0.1, 0.3, 0.5, 0.7):
                cells.append({"index": i, "binary": binary,
                              "n_groups": n_groups, "group_size": 20,
                              "share": share, "unbalance": "balanced",
                              "law": "two_point", "replicates": 400})
                i += 1
    # A rare call: eight per cent prevalence, the level of the poultry
    # quinolone cell of the article, where the exact realised interval
    # refuses by construction and clonal_share is the estimate.
    for n_groups in (30, 100):
        for share in (0.1, 0.3, 0.5, 0.7):
            cells.append({"index": i, "binary": True, "n_groups": n_groups,
                          "group_size": 20, "share": share,
                          "unbalance": "balanced", "prevalence": 0.08,
                          "replicates": 400})
            i += 1
    for c in cells:
        c.setdefault("law", "normal")
        c.setdefault("prevalence", 0.25)
    return cells


def run_cell(cell: dict) -> dict:
    rng = np.random.default_rng(SEED + 977 * cell["index"])
    stats_by = {k: {"err": [], "err_other": [], "cover": 0, "cover_other": 0, "width": [],
                    "seconds": 0.0, "estimable": 0}
                for k in ESTIMATORS}
    for _ in range(cell["replicates"]):
        y, lineage, truth = cohort(rng, cell["n_groups"], cell["group_size"],
                                   cell["share"], binary=cell["binary"],
                                   prevalence=cell["prevalence"],
                                   unbalance=cell["unbalance"],
                                   law=cell["law"])
        for name, fn in ESTIMATORS.items():
            t0 = time.perf_counter()
            out = fn(y, lineage)
            stats_by[name]["seconds"] += time.perf_counter() - t0
            own = truth[TARGET[name]]
            other = truth[OTHER[name]]
            point, lo, hi = out["point"], out["low"], out["high"]
            if np.isfinite(point):
                stats_by[name]["err"].append(point - own)
                stats_by[name]["err_other"].append(point - other)
            if np.isfinite(lo) and np.isfinite(hi):
                stats_by[name]["cover"] += int(lo <= own <= hi)
                stats_by[name]["cover_other"] += int(lo <= other <= hi)
                stats_by[name]["width"].append(hi - lo)
            stats_by[name]["estimable"] += int(out.get("estimable", True))
    result = dict(cell)
    reps = cell["replicates"]
    for name, s in stats_by.items():
        err = np.array(s["err"], dtype=float)
        err_other = np.array(s["err_other"], dtype=float)
        width = np.array(s["width"], dtype=float)
        result[name] = {
            "bias": float(err.mean()) if err.size else None,
            "rmse": float(np.sqrt((err ** 2).mean())) if err.size else None,
            "bias_of_the_other_target": (float(err_other.mean())
                                         if err_other.size else None),
            "rmse_of_the_other_target": (float(np.sqrt((err_other ** 2).mean()))
                                         if err_other.size else None),
            "target": TARGET[name],
            "other_target": OTHER[name],
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
    gaussian = [r for r in rows if not r["binary"] and r["share"] > 0
                and r.get("law", "normal") == "normal"]
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
          f"unbalance={cell['unbalance']} law={cell['law']} "
          f"prevalence={cell['prevalence']}", flush=True)
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
