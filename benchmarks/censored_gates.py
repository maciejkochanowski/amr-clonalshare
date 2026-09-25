#!/usr/bin/env python3
"""The evidence for the two estimability gates of the moment estimator.

    python benchmarks/censored_gates.py --out benchmarks/results_censored_gates

``censored.censored_clonal_share`` refuses a share in two situations, and each
limit in ``censored.py`` is read off one sweep here, against a known latent
share, on cohorts of 30 lineages of 25 isolates with Gaussian lineage effects
and residuals.

**Lineages wholly beyond the panel** (``CENSORED_GROUP_LIMIT``). A single cut
point is moved across the latent distribution at a true share of 0.5, so that
more and more lineages hold only one-sided readings. The sweep records the
share of isolates in such lineages and the bias of the estimate.

**A single cut point in a tail** (``SINGLE_CUT_PREVALENCE``). A binary call is
made at a prevalence from 4 to 50 per cent at true shares of 0, 0.3 and 0.5;
the sweep records the bias of the estimate and the coverage of its
approximate F interval.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import norm

from amr_clonalshare.censored import censored_clonal_share, intervals_from_binary

SEED = 20260925
LINEAGES, PER_LINEAGE = 30, 25


def cohort(rng, share):
    code = np.repeat(np.arange(LINEAGES), PER_LINEAGE)
    effects = rng.normal(0.0, np.sqrt(share), LINEAGES)
    z = effects[code] + rng.normal(0.0, np.sqrt(1.0 - share), code.size)
    return z, np.array([f"L{g}" for g in code], dtype=object)


def beyond_panel(reps):
    rows = []
    for level, cut in enumerate(np.linspace(-1.0, 3.0, 21)):
        rng = np.random.default_rng([SEED, 1, level])
        errors, shares = [], []
        for _ in range(reps):
            z, lineage = cohort(rng, 0.5)
            below = z <= cut
            r = censored_clonal_share(np.where(below, -np.inf, cut), np.where(below, cut, np.inf), lineage)
            shares.append(r.share_in_censored_groups)
            if np.isfinite(r.kappa):
                errors.append(r.kappa - 0.5)
        rows.append(dict(cut=float(cut), replicates=reps, scored=len(errors),
                         isolate_share_in_censored_lineages=float(np.mean(shares)),
                         bias=float(np.mean(errors)) if errors else None,
                         sd=float(np.std(errors, ddof=1)) if len(errors) > 1 else None))
    return rows


def single_cut(reps):
    rows = []
    for i, share in enumerate((0.0, 0.3, 0.5)):
        for j, prevalence in enumerate((0.04, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20, 0.24, 0.30, 0.40, 0.50)):
            rng = np.random.default_rng([SEED, 2, i, j])
            cut = float(norm.ppf(1.0 - prevalence))
            errors, covered, finite = [], 0, 0
            for _ in range(reps):
                z, lineage = cohort(rng, share)
                lo, hi = intervals_from_binary((z > cut).astype(float), cutoff_log2=cut)
                r = censored_clonal_share(lo, hi, lineage)
                if np.isfinite(r.kappa):
                    errors.append(r.kappa - share)
                if np.isfinite(r.ci_low) and np.isfinite(r.ci_high):
                    finite += 1
                    covered += int(r.ci_low <= share <= r.ci_high)
            rows.append(dict(true_share=share, prevalence=prevalence, replicates=reps,
                             scored=len(errors), bias=float(np.mean(errors)) if errors else None,
                             intervals=finite, coverage=covered / reps))
    return rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--reps", type=int, default=200)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    payload = {"design": {"seed": SEED, "lineages": LINEAGES, "isolates_per_lineage": PER_LINEAGE,
                          "replicates": args.reps},
               "beyond_panel": beyond_panel(args.reps), "single_cut": single_cut(args.reps)}
    (args.out / "censored_gates.json").write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
