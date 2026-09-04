#!/usr/bin/env python3
"""Is the clonal-share permutation p-value uniform when there is nothing to
find, and how fast does it find something when there is?

A p-value that is valid is uniform under the null, or stochastically larger.
Reporting only the rejection rate at 0.05 checks one point of that
distribution; this script checks the whole of it. It is the frequentist
counterpart of simulation-based calibration: draw cohorts from the null the
test claims to control, run the test, and compare the p-values against the
distribution they must have.

The permutation p-value is discrete. With ``n_perm`` permutations it takes the
values ``(b + 1) / (n_perm + 1)`` and under exchangeability it is uniform on
that grid, not on the unit interval. The check is therefore against the grid
uniform: the largest gap between the empirical distribution and the grid
distribution, with its null reference drawn by Monte Carlo from the grid. A
test against the continuous uniform would reject a correct p-value for being
discrete.

    python benchmarks/null_uniformity.py --out benchmarks/results_null_uniformity
    python benchmarks/null_uniformity.py --quick
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from datetime import date
from pathlib import Path

import numpy as np
import scipy

from amr_clonalshare import attribution
from amr_clonalshare.clonality import _clopper_pearson

ALPHAS = (0.01, 0.05, 0.10)
N_PERM = 199
SETTINGS = dict(folds=5, repeats=5, n_boot=0, n_perm=N_PERM)


def graded_lineages(rng, n: int, n_lineages: int) -> np.ndarray:
    """Lineage sizes that fall off geometrically, as surveillance cohorts do:
    a few large clones and a tail of small ones."""
    w = 0.7 ** np.arange(n_lineages)
    return rng.choice(n_lineages, size=n, p=w / w.sum())


def draw_cohort(rng, n: int, n_lineages: int, spread: float,
                prevalence: float = 0.35):
    """Trait calls whose per-lineage rate is ``prevalence`` shifted by a
    lineage effect of standard deviation ``spread`` on the probability scale.
    ``spread = 0`` is the null: the label carries nothing."""
    lineage = graded_lineages(rng, n, n_lineages)
    effect = rng.normal(0.0, spread, size=n_lineages) if spread > 0 else np.zeros(n_lineages)
    rate = np.clip(prevalence + effect[lineage], 0.02, 0.98)
    y = (rng.random(n) < rate).astype(float)
    return y, lineage


def grid_gap(p: np.ndarray, n_perm: int) -> float:
    """Largest distance between the empirical distribution of ``p`` and the
    uniform distribution on the permutation grid, evaluated at the grid."""
    grid = np.arange(1, n_perm + 2) / (n_perm + 1)
    ecdf = np.searchsorted(np.sort(p), grid, side="right") / p.size
    return float(np.max(np.abs(ecdf - grid)))


def grid_gap_reference(rng, n_rep: int, n_perm: int, n_draws: int) -> np.ndarray:
    grid = np.arange(1, n_perm + 2) / (n_perm + 1)
    out = np.empty(n_draws)
    for i in range(n_draws):
        out[i] = grid_gap(rng.choice(grid, size=n_rep), n_perm)
    return out


def run(n_rep: int, n_power: int, n: int, n_lineages: int, seed: int) -> dict:
    root = np.random.SeedSequence(seed)
    null_ss, power_ss, ref_ss = root.spawn(3)

    # The null: one child stream per replicate, so a replicate can be rerun
    # on its own and the set does not depend on how many were asked for.
    p_null, kappa_null = [], []
    for child in null_ss.spawn(n_rep):
        rng = np.random.default_rng(child)
        y, lineage = draw_cohort(rng, n, n_lineages, spread=0.0)
        r = attribution.clonal_share(y, lineage, seed=rng, **SETTINGS)
        p_null.append(r.p_value)
        kappa_null.append(r.kappa_adj)
    p_null = np.asarray(p_null)
    kappa_null = np.asarray(kappa_null)

    gap = grid_gap(p_null, N_PERM)
    ref = grid_gap_reference(np.random.default_rng(ref_ss), n_rep, N_PERM, 4000)
    gap_p = float((np.sum(ref >= gap) + 1) / (ref.size + 1))

    level = {}
    for a in ALPHAS:
        x = int(np.sum(p_null <= a))
        lo, hi = _clopper_pearson(x, n_rep)
        level["%.2f" % a] = {"rejections": x, "rate": x / n_rep,
                             "cp95_low": lo, "cp95_high": hi,
                             "holds": hi <= a + 0.02 or x / n_rep <= a}

    # Power against a graded lineage effect, at 0.05.
    spreads = (0.0, 0.05, 0.10, 0.15, 0.20, 0.30)
    power = {}
    for spread, ss in zip(spreads, power_ss.spawn(len(spreads))):
        rej, kap = 0, []
        for child in ss.spawn(n_power):
            rng = np.random.default_rng(child)
            y, lineage = draw_cohort(rng, n, n_lineages, spread=spread)
            r = attribution.clonal_share(y, lineage, seed=rng, **SETTINGS)
            rej += int(r.p_value <= 0.05)
            kap.append(r.kappa_adj)
        lo, hi = _clopper_pearson(rej, n_power)
        power["%.2f" % spread] = {"rejections": rej, "power": rej / n_power,
                                  "cp95_low": lo, "cp95_high": hi,
                                  "median_kappa_adj": float(np.nanmedian(kap))}

    return {
        "what": "uniformity of the clonal-share permutation p-value under a "
                "null with no lineage effect, and power against a graded "
                "lineage effect",
        "date": date.today().isoformat(),
        "design": {"n_isolates": n, "n_lineages": n_lineages,
                   "lineage_sizes": "geometric, ratio 0.7",
                   "prevalence": 0.35, "n_null_replicates": n_rep,
                   "n_power_replicates_per_spread": n_power,
                   "estimator_settings": SETTINGS, "seed": seed,
                   "seeding": "numpy SeedSequence, one spawned child per "
                              "replicate"},
        "null": {
            "grid_gap": gap,
            "grid_gap_reference_draws": int(ref.size),
            "grid_gap_p_value": gap_p,
            "uniform_on_the_grid": gap_p > 0.05,
            "level": level,
            "kappa_adj_mean": float(np.nanmean(kappa_null)),
            "kappa_adj_sd": float(np.nanstd(kappa_null)),
            "p_value_deciles": [float(v) for v in
                                np.quantile(p_null, np.linspace(0.1, 0.9, 9))],
        },
        "power_at_0.05": power,
        "reading_rule": "The p-value is called uniform when the grid-gap "
                        "reference p-value exceeds 0.05, and the level holds at "
                        "alpha when the observed rejection rate is at or below "
                        "alpha or its exact 95 % upper limit is within 0.02 of "
                        "it. Both rules were fixed before the run.",
        "environment": {"python": platform.python_version(),
                        "numpy": np.__version__, "scipy": scipy.__version__},
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="benchmarks/results_null_uniformity")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=20260904)
    args = ap.parse_args()
    q = args.quick
    result = run(n_rep=200 if q else 2000, n_power=50 if q else 400,
                 n=240, n_lineages=40, seed=args.seed)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / ("null_uniformity_quick.json" if q else "null_uniformity.json")
    text = json.dumps(result, indent=2)
    path.write_text(text + "\n", encoding="utf-8")
    print("wrote %s sha256 %s" % (path, hashlib.sha256(text.encode()).hexdigest()[:12]))
    print(json.dumps(result["null"], indent=2)[:1500])
    print(json.dumps(result["power_at_0.05"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
