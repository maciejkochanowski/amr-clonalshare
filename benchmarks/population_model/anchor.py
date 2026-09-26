#!/usr/bin/env python3
"""Recompute the population-model likelihood by adaptive quadrature.

    python -m benchmarks.population_model.anchor --output OUT/anchor.json

For datasets drawn from the validation cells the log-likelihood the package
maximises is recomputed at its fitted point and at the generating rho, by
``scipy.integrate.quad`` over the lineage effect, from the model statement:
an isolate is positive with probability Phi((a + u) / sqrt(1 - rho)) given
its lineage effect u ~ N(0, rho), where Phi(a) is the marginal prevalence.
The study is void if any difference exceeds 1e-6 (PROTOCOL.md).
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy import integrate
from scipy.special import gammaln
from scipy.stats import norm

from amr_clonalshare import _population_probit_numerics as numerics
from benchmarks.population_model.design import SEEDS, cells, draw

TOLERANCE = 1e-6


def _log_integral(kk, mm, a, rho):
    """log of the integral over z ~ N(0, 1) of the binomial probability of kk."""
    scale = math.sqrt(1.0 - rho)
    root = math.sqrt(rho)

    def log_kernel(z):
        eta = (a + root * z) / scale
        return kk * norm.logcdf(eta) + (mm - kk) * norm.logcdf(-eta) - 0.5 * z * z
    grid = np.linspace(-12, 12, 4801)
    values = log_kernel(grid)
    peak, mode = float(values.max()), float(grid[values.argmax()])
    value, _ = integrate.quad(lambda z: math.exp(float(log_kernel(z)) - peak), -12, 12,
                              epsabs=0, epsrel=1e-12, limit=500, points=[mode])
    return peak + math.log(value) - 0.5 * math.log(2 * math.pi)


def loglik(k, m, a, rho):
    total = 0.0
    pairs, freq = np.unique(np.c_[k, m], axis=0, return_counts=True)
    for (kk, mm), f in zip(pairs, freq):
        log_choose = gammaln(mm + 1) - gammaln(kk + 1) - gammaln(mm - kk + 1)
        if rho == 0.0:
            p = norm.cdf(a)
            term = kk * math.log(p) + (mm - kk) * math.log1p(-p)
        else:
            term = _log_integral(int(kk), int(mm), a, rho)
        total += f * (log_choose + term)
    return total


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--datasets", type=int, default=300)
    args = parser.parse_args(argv)
    pool = [c for c in cells() if c["family"] == "validation" and max(c["sizes"]) <= 200]
    rng = np.random.default_rng(SEEDS["anchor"])
    rows = []
    for i in range(args.datasets):
        cell = pool[int(rng.integers(len(pool)))]
        k, m = draw(cell, np.random.default_rng([SEEDS["anchor"], i]))
        if k.sum() in (0, m.sum()) or (m >= 2).sum() < 2:
            continue
        like = numerics.ProfileLikelihood(k, m)
        rho_hat, a_hat, best = like.mle()
        rho_hat = min(max(rho_hat, 0.0), 1 - 1e-9)
        nll_truth, a_truth = like.profile(cell["rho"])
        package = (-like.nll(a_hat, rho_hat), -nll_truth)
        anchor = (loglik(k, m, a_hat, rho_hat), loglik(k, m, a_truth, cell["rho"]))
        rows.append(dict(cell=cell["cell"], dataset=i, rho=cell["rho"], rho_hat=rho_hat,
                         package=package, anchor=anchor,
                         difference=max(abs(x - y) for x, y in zip(package, anchor))))
    worst = max(r["difference"] for r in rows)
    payload = dict(tolerance=TOLERANCE, datasets=len(rows), worst_difference=worst,
                   void=bool(worst > TOLERANCE), rows=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: payload[k] for k in ("datasets", "worst_difference", "void")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
