#!/usr/bin/env python3
"""Operating characteristics of the realised share, and of the question beside it.

    python benchmarks/realised_calibration.py --cell 0 --out <dir>
    python benchmarks/realised_calibration.py --aggregate --out <dir>

ESTIMAND. Two questions share a point estimate and must not share an interval.
The superpopulation share asks what would sit between lineages drawn afresh from
the species; the realised share asks what sits between the lineages this cohort
holds. Each interval is scored against its own truth, computed from the effects
actually drawn in that replicate, never from the sample that built the interval.

DESIGN. A calibration family the method is derived against, four validation
families each breaking one assumption at matched theoretical moments, and an
anchor arm that re-derives the reference law from the model statement rather
than from the implementation: given the realised effects, the between-lineage
sum of squares is a noncentral chi-square and the within-lineage sum of squares
an independent central one, so the anchor draws those two laws and never forms a
trait value at all.

VOID CONDITION, fixed before the run. If the anchor arm does not reproduce the
calibration arm to within four Monte-Carlo standard errors, the table is void
and nothing in it is reported. If the cross-scoring arm does not break coverage
when each interval is scored against the other estimand's truth, the two
questions are not distinguishable and the realised interval has no reason to
exist.

DETECTION. Coverage says the interval contains the truth; it does not say
whether the interval is any use as a decision. A laboratory reads one thing off
it: does the lower endpoint stand above zero, that is, do these lineages differ
in this trait at all. Sensitivity is the share of replicates where it does when
the lineages really do differ, and specificity the share where it does not when
they do not. The detection family measures both over a share grid at two
lineage counts, on the same cells and the same draws that give the coverage.

SECOND VOID CONDITION, fixed before the run. At a true realised share of
exactly zero every lineage effect is exactly zero, the noncentrality is zero,
and the F ratio is central, so the lower endpoint must stand above zero at the
nominal one-sided rate alpha / 2 = 0.025 and no oftener. A measured rate of
zero would mean the arm cannot fire and measures nothing; a rate near alpha
would mean the interval is not the inversion it is claimed to be. Either voids
the detection table.

WHAT IS RECORDED THOUGH IT FAILED. Two repairs for heavy tails were tried and
are kept as cells rather than deleted: a Box-type deflation of the degrees of
freedom by the estimated residual kurtosis, and a within-lineage percentile
bootstrap intended as an assumption-free fallback.
"""
from __future__ import annotations

import argparse
import json
import platform
from datetime import date
from pathlib import Path

import numpy as np
import scipy
from scipy import stats

from amr_clonalshare.realised import (realised_interval, superpopulation_interval)

ALPHA = 0.05
SEED = 20260902
BOOTSTRAP_DRAWS = 300


# ------------------------------------------------------------------ generators
def _effects(rng, n_groups, share, shape="gaussian"):
    tau = np.sqrt(share / (1.0 - share)) if share < 1.0 else 0.0
    if shape == "two-point":
        return tau * rng.choice([-1.0, 1.0], n_groups)
    return rng.normal(0.0, tau, n_groups)


def _noise(rng, size, kind):
    """Residuals at unit theoretical variance and zero theoretical mean.

    Standardising a draw by its own sample moments would condition the sample
    and remove the randomness the interval exists to quantify, so every family
    is matched on the moments of its law.
    """
    if kind == "gaussian":
        return rng.normal(0.0, 1.0, size)
    if kind.startswith("t"):
        df = float(kind[1:])
        return rng.standard_t(df, size) / np.sqrt(df / (df - 2.0))
    if kind == "skew":
        alpha = 4.0
        delta = alpha / np.sqrt(1.0 + alpha ** 2)
        raw = stats.skewnorm.rvs(alpha, size=size, random_state=rng)
        return ((raw - delta * np.sqrt(2.0 / np.pi))
                / np.sqrt(1.0 - 2.0 * delta ** 2 / np.pi))
    if kind == "contaminated":
        eps, shift = 0.10, 3.0
        raw = rng.normal(0.0, 1.0, size) + (rng.random(size) < eps) * shift
        return ((raw - eps * shift)
                / np.sqrt(1.0 + eps * (1.0 - eps) * shift ** 2))
    raise ValueError(f"unknown noise family {kind!r}")


def _bernoulli_cohort(rng, n_groups, sizes, share, prevalence):
    """A binary trait with a known realised share on the observed scale.

    The atlas reads non-susceptibility as a 0/1 indicator, so the share it
    reports is the intraclass correlation of that indicator and not a
    liability-scale quantity. The residuals of a Bernoulli are not Gaussian by
    construction: at a middling prevalence they are platykurtic, which the
    F approximation tolerates, and at an extreme one they are heavy, which the
    kurtosis gate is there to catch. This arm is what licenses quoting the
    realised interval on a binary panel at all.
    """
    latent = _effects(rng, n_groups, share, "gaussian")
    cut = stats.norm.ppf(prevalence)
    probability = stats.norm.cdf(latent + cut)
    y = np.concatenate([(rng.random(int(sizes[i])) < probability[i]).astype(float)
                        for i in range(n_groups)])
    lineage = np.repeat(np.arange(n_groups), sizes)
    total = sizes.sum()
    n0 = (total - (sizes ** 2).sum() / total) / (n_groups - 1)
    weighted = np.average(probability, weights=sizes)
    between = ((sizes * (probability - weighted) ** 2).sum()
               / ((n_groups - 1) * n0))
    within = float((sizes * probability * (1.0 - probability)).sum() / total)
    return y, lineage, between / (between + within)


def _cohort(rng, n_groups, sizes, share, noise, effect_shape):
    effects = _effects(rng, n_groups, share, effect_shape)
    y = np.concatenate([effects[i] + _noise(rng, int(sizes[i]), noise)
                        for i in range(n_groups)])
    lineage = np.repeat(np.arange(n_groups), sizes)
    weighted = np.average(effects, weights=sizes)
    n0 = (sizes.sum() - (sizes ** 2).sum() / sizes.sum()) / (n_groups - 1)
    realised = (sizes * (effects - weighted) ** 2).sum() / ((n_groups - 1) * n0)
    return y, lineage, realised / (realised + 1.0)


def _mean_squares(y, lineage):
    _, code = np.unique(lineage, return_inverse=True)
    g = int(code.max()) + 1
    n = y.size
    counts = np.bincount(code, minlength=g).astype(float)
    means = np.bincount(code, weights=y, minlength=g) / counts
    ssb = float((counts * (means - y.mean()) ** 2).sum())
    residual = y - means[code]
    ssw = float((residual ** 2).sum())
    n0 = (n - (counts ** 2).sum() / n) / (g - 1)
    return (ssb / (g - 1)) / (ssw / (n - g)), g, n, n0, residual


# ----------------------------------------------------------------- refuted arms
def _box_corrected(f_ratio, df1, df2, n0, g, residual):
    """Deflate the degrees of freedom by the estimated residual kurtosis.

    Recorded because it was tried and refuted, not because it is used.
    """
    excess = float(stats.kurtosis(residual, fisher=True, bias=False))
    inflate = max(1.0 + excess / 2.0, 1.0)
    return realised_interval(f_ratio, df1 / inflate, df2 / inflate, n0, g,
                             ALPHA)


def _within_lineage_bootstrap(y, lineage, rng, draws=BOOTSTRAP_DRAWS):
    """Resample isolates inside each lineage, holding the lineage set fixed.

    Recorded because it was tried and refuted: it under-covers the very target
    it was built for, under the process the method is derived against.
    """
    _, code = np.unique(lineage, return_inverse=True)
    g = int(code.max()) + 1
    index = [np.flatnonzero(code == i) for i in range(g)]
    out = np.empty(draws)
    for b in range(draws):
        drawn = np.concatenate([idx[rng.integers(0, idx.size, idx.size)]
                                for idx in index])
        f, gg, n, n0, _ = _mean_squares(y[drawn], code[drawn])
        msw = 1.0
        between = max((f - 1.0) / n0, 0.0)
        out[b] = between / (between + msw / 1.0) if between >= 0 else 0.0
    return tuple(np.percentile(out, [100 * ALPHA / 2, 100 * (1 - ALPHA / 2)]))


# ------------------------------------------------------------------- one cell
def _sizes(rng, n_groups, size, unbalance):
    if unbalance == "balanced":
        return np.full(n_groups, size, dtype=int)
    if unbalance == "poisson":
        return 2 + rng.poisson(max(size - 2, 1), n_groups)
    if unbalance == "long-tailed":
        base = np.full(n_groups, 2, dtype=int)
        heavy = rng.random(n_groups) < 0.25
        base[heavy] = size * 4
        return base
    raise ValueError(unbalance)


def run_cell(cell: dict) -> dict:
    rng = np.random.default_rng(SEED + cell["index"])
    g, m, share = cell["n_groups"], cell["group_size"], cell["share"]
    reps = cell["replicates"]
    counters = {k: 0 for k in ("realised", "superpopulation",
                               "realised_wrong_target",
                               "superpopulation_wrong_target",
                               "box", "bootstrap", "gate_open",
                               "detect_realised", "detect_superpopulation")}
    widths = {k: [] for k in ("realised", "superpopulation", "box",
                              "bootstrap")}
    kurtoses = []
    for _ in range(reps):
        sizes = _sizes(rng, g, m, cell["unbalance"])
        if cell["arm"] == "anchor":
            effects = _effects(rng, g, share, cell["effects"])
            realised_truth = ((effects - effects.mean()) ** 2).sum() / (g - 1)
            realised_truth = realised_truth / (realised_truth + 1.0)
            lam = m * ((effects - effects.mean()) ** 2).sum()
            ssb = stats.ncx2.rvs(g - 1, lam, random_state=rng)
            ssw = stats.chi2.rvs(g * m - g, random_state=rng)
            f_ratio = (ssb / (g - 1)) / (ssw / (g * m - g))
            n, n0, residual = g * m, float(m), None
        elif cell["prevalence"] is not None:
            y, lineage, realised_truth = _bernoulli_cohort(
                rng, g, sizes, share, cell["prevalence"])
            f_ratio, g, n, n0, residual = _mean_squares(y, lineage)
        else:
            y, lineage, realised_truth = _cohort(rng, g, sizes, share,
                                                 cell["noise"],
                                                 cell["effects"])
            f_ratio, g, n, n0, residual = _mean_squares(y, lineage)
        df1, df2 = float(g - 1), float(n - g)
        low, high = realised_interval(f_ratio, df1, df2, n0, g, ALPHA)
        slow, shigh = superpopulation_interval(f_ratio, df1, df2, n0, ALPHA)
        counters["realised"] += low <= realised_truth <= high
        if cell["prevalence"] is None:
            counters["superpopulation"] += slow <= share <= shigh
            counters["realised_wrong_target"] += low <= share <= high
        counters["superpopulation_wrong_target"] += (slow <= realised_truth
                                                     <= shigh)
        # The decision a laboratory takes off the interval: does the lower
        # endpoint stand above zero. Counted on every cell, so sensitivity and
        # specificity come from the same draws as the coverage.
        counters["detect_realised"] += low > 0.0
        counters["detect_superpopulation"] += slow > 0.0
        widths["realised"].append(high - low)
        widths["superpopulation"].append(shigh - slow)
        if residual is not None:
            excess = float(stats.kurtosis(residual, fisher=True, bias=False))
            kurtoses.append(excess)
            if cell["refuted_arms"]:
                blow, bhigh = _box_corrected(f_ratio, df1, df2, n0, g, residual)
                counters["box"] += blow <= realised_truth <= bhigh
                widths["box"].append(bhigh - blow)
                clow, chigh = _within_lineage_bootstrap(y, lineage, rng)
                counters["bootstrap"] += clow <= realised_truth <= chigh
                widths["bootstrap"].append(chigh - clow)
    out = dict(cell)
    out["monte_carlo_se"] = float(np.sqrt(0.05 * 0.95 / reps))
    # On a binary trait the superpopulation share is not defined against the
    # latent variance the cohort was drawn with, so those two counters are not
    # scored at all. They are reported as absent rather than as zero: a
    # coverage of exactly zero would read as a failure instead of a column that
    # was never filled.
    unscored = ({"superpopulation", "realised_wrong_target"}
                if cell["prevalence"] is not None else set())
    if not cell["refuted_arms"]:
        unscored |= {"box", "bootstrap"}
    for key, hits in counters.items():
        if key.startswith("detect_"):
            out[key] = hits / reps
            continue
        out[f"coverage_{key}"] = None if key in unscored else hits / reps
    for key, values in widths.items():
        out[f"median_width_{key}"] = (float(np.median(values)) if values
                                      else None)
    out["median_excess_kurtosis"] = (float(np.median(kurtoses)) if kurtoses
                                     else None)
    return out


# ---------------------------------------------------------------- the design
def _cell(index, label, family, **kw):
    base = dict(index=index, label=label, family=family, arm="observed",
                n_groups=30, group_size=20, share=0.30, noise="gaussian",
                effects="gaussian", unbalance="balanced", replicates=4000,
                refuted_arms=False, prevalence=None)
    base.update(kw)
    return base


def design() -> list[dict]:
    cells, i = [], 0

    def add(label, family, **kw):
        nonlocal i
        cells.append(_cell(i, label, family, **kw))
        i += 1

    add("gaussian", "calibration")
    add("model law drawn directly", "anchor", arm="anchor")
    for kind, name in (("t4", "t(4) residuals"), ("t6", "t(6) residuals"),
                       ("skew", "skew-normal residuals"),
                       ("contaminated", "contaminated residuals")):
        add(name, "validation", noise=kind)
    add("two-point lineage effects", "validation", effects="two-point")
    add("cross-scoring", "falsifiability")

    for g in (5, 10, 15, 20, 30, 50, 80, 120):
        add(f"lineages = {g}", "grid: lineages", n_groups=g, replicates=2000)
    for m in (3, 5, 10, 20, 40, 80):
        add(f"isolates per lineage = {m}", "grid: isolates", n_groups=20,
            group_size=m, replicates=2000)
    for share in (0.0, 0.05, 0.10, 0.30, 0.50, 0.70, 0.90):
        add(f"true share = {share:.2f}", "grid: share", share=share,
            replicates=2000)
    for shape in ("poisson", "long-tailed"):
        add(f"lineage sizes {shape}", "grid: unbalance", unbalance=shape,
            replicates=2000)
    for df in (4, 5, 6, 8, 10, 15, 30):
        add(f"t({df}) residuals", "gate: kurtosis", noise=f"t{df}",
            replicates=3000)
    add("gaussian residuals", "gate: kurtosis", replicates=3000)
    for prevalence in (0.05, 0.10, 0.25, 0.50):
        add(f"binary trait, prevalence {prevalence:.2f}", "binary",
            prevalence=prevalence, replicates=3000)
    for g in (10, 20, 30, 60):
        add(f"binary trait, {g} lineages", "binary: lineages",
            prevalence=0.25, n_groups=g, replicates=2000)
    add("refuted repairs, gaussian", "refuted", replicates=1000,
        refuted_arms=True)
    add("refuted repairs, t(4) residuals", "refuted", replicates=1000,
        refuted_arms=True, noise="t4")
    # Appended last on purpose: each cell is seeded from its index, so adding a
    # family at the end leaves every cell already reported bit-for-bit as it
    # was, and the campaign can be re-run as one piece rather than as two.
    for g in (10, 30):
        for share in (0.00, 0.02, 0.05, 0.10, 0.20):
            add(f"{g} lineages, true share {share:.2f}", "detection",
                n_groups=g, share=share, replicates=4000)
    return cells


def _provenance() -> dict:
    return {"generated": date.today().isoformat(), "seed": SEED,
            "alpha": ALPHA, "python": platform.python_version(),
            "numpy": np.__version__, "scipy": scipy.__version__,
            "bootstrap_draws": BOOTSTRAP_DRAWS}


def aggregate(out: Path) -> dict:
    rows = [json.loads(p.read_text()) for p in sorted((out / "cells").glob("*.json"))]
    by_label = {r["label"]: r for r in rows if r["family"] in ("calibration",
                                                              "anchor")}
    verdict = {}
    if "gaussian" in by_label and "model law drawn directly" in by_label:
        cal = by_label["gaussian"]["coverage_realised"]
        anc = by_label["model law drawn directly"]["coverage_realised"]
        se = by_label["gaussian"]["monte_carlo_se"]
        verdict["anchor_reproduces_calibration"] = abs(anc - cal) <= 4 * se
        verdict["anchor_difference"] = anc - cal
    cross = [r for r in rows if r["family"] == "falsifiability"]
    if cross:
        c = cross[0]
        verdict["estimands_distinguishable"] = (
            c["coverage_realised_wrong_target"] < 0.90
            and c["coverage_superpopulation_wrong_target"] > 0.98)
    gate = sorted((r for r in rows if r["family"] == "gate: kurtosis"),
                  key=lambda r: r["median_excess_kurtosis"] or 0.0)
    limit = None
    for row in gate:
        if row["coverage_realised"] >= 0.93:
            limit = row["median_excess_kurtosis"]
    verdict["kurtosis_limit_from_the_curve"] = limit
    # The second void condition, scored here rather than read off the table by
    # eye. At a true share of zero every effect is zero, so the F ratio is
    # central and the lower endpoint stands above zero at exactly the nominal
    # one-sided rate. A rate of zero would mean the arm cannot fire.
    detect = [r for r in rows if r["family"] == "detection"]
    if detect:
        nulls = [r for r in detect if r["share"] == 0.0]
        fired = [r for r in detect if r["share"] > 0.0]
        reps = min(r["replicates"] for r in nulls)
        tol = 4.0 * float(np.sqrt(ALPHA / 2 * (1 - ALPHA / 2) / reps))
        verdict["nominal_one_sided_rate"] = ALPHA / 2
        verdict["false_positive_rate_at_zero"] = [r["detect_realised"]
                                                  for r in nulls]
        verdict["false_positive_rate_matches_nominal"] = all(
            abs(r["detect_realised"] - ALPHA / 2) <= tol for r in nulls)
        verdict["detection_arm_has_leverage"] = (
            max(r["detect_realised"] for r in fired) > 0.5)
    payload = {"provenance": _provenance(), "verdict": verdict, "cells": rows}
    (out / "realised_calibration.json").write_text(
        json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(verdict, indent=2))
    print(f"cells: {len(rows)}")
    print(f"written: {out / 'realised_calibration.json'}")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cell", type=int)
    parser.add_argument("--aggregate", action="store_true")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()
    cells = design()
    if args.list:
        for c in cells:
            print(f"{c['index']:>3} {c['family']:<20} {c['label']}")
        print(f"total {len(cells)}")
        return 0
    args.out.mkdir(parents=True, exist_ok=True)
    if args.aggregate:
        aggregate(args.out)
        return 0
    if args.cell is None:
        raise SystemExit("give --cell, --aggregate or --list")
    (args.out / "cells").mkdir(exist_ok=True)
    result = run_cell(cells[args.cell])
    path = args.out / "cells" / f"{args.cell:03d}.json"
    path.write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")
    show = lambda v: "     -" if v is None else f"{v:.4f}"
    print(f"{result['family']:<20}{result['label']:<32}"
          f"realised {show(result['coverage_realised'])}  "
          f"superpopulation {show(result['coverage_superpopulation'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
