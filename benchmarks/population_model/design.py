"""Prespecified cells of the population-model study; see PROTOCOL.md.

Every cell names its family, the lineage sizes, the marginal prevalence p, the
liability ICC rho and the law of the lineage effects. Nothing here reads a
result.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.stats import norm

#: Root seeds, one per family, fixed before any run.
SEEDS = {"calibration": 20260926101, "validation": 20260926202,
         "robustness": 20260926303, "benefit": 20260926404, "anchor": 20260926505}

#: Lineage sizes of the S. suis example collection (hierBAPS clusters).
SSUIS_SIZES = (161, 60, 57, 50, 44, 42, 40, 36, 31, 22, 20, 17, 13, 11, 11, 10,
               8, 7, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 1)


def skewed(groups: int) -> list[int]:
    """A few large lineages and a long tail of small ones and singletons."""
    return [max(1, round(60 * 0.85 ** i)) for i in range(groups)]


def _sizes(pattern: str, groups: int) -> list[int]:
    if pattern.startswith("balanced"):
        return [int(pattern[len("balanced"):])] * groups
    if pattern == "unequal":
        cycle = (1, 2, 5, 10, 20, 50)
        return [cycle[i % len(cycle)] for i in range(groups)]
    if pattern == "skewed":
        return skewed(groups)
    if pattern == "ssuis":
        return list(SSUIS_SIZES)
    raise ValueError(pattern)


def cells() -> list[dict]:
    out = []

    def add(family, pattern, groups, p, rho, law="gaussian", replicates=2000):
        out.append(dict(cell=len(out), family=family, pattern=pattern, groups=groups,
                        sizes=_sizes(pattern, groups), prevalence=p, rho=rho, law=law,
                        replicates=replicates))

    # Calibration of the fixed cut-off: the Gaussian model on a grid.
    for groups in (10, 30, 100):
        for pattern in ("balanced5", "balanced20", "unequal"):
            for p in (0.05, 0.15, 0.30, 0.50):
                for rho in (0.0, 0.05, 0.2, 0.4, 0.6, 0.8):
                    add("calibration", pattern, groups, p, rho)
    # Validation: the Gaussian model off the calibration grid, fresh seed.
    for groups in (15, 50, 200):
        for pattern in ("balanced10", "skewed"):
            for p in (0.02, 0.10, 0.25, 0.85):
                for rho in (0.1, 0.3, 0.5, 0.7, 0.9):
                    add("validation", pattern, groups, p, rho)
    for p in (0.04, 0.24, 0.54, 0.85):
        for rho in (0.0, 0.1, 0.3, 0.5, 0.7, 0.9):
            add("validation", "ssuis", len(SSUIS_SIZES), p, rho)
    # Robustness: one assumption broken at a time.
    for law in ("t4", "two_point", "informative_sizes"):
        for p in (0.10, 0.25):
            for rho in (0.3, 0.7):
                add("robustness", "skewed" if law == "informative_sizes" else "balanced10",
                    50, p, rho, law=law)
    # Benefit: one liability ICC read at different prevalences.
    for rho in (0.3, 0.6):
        for p in (0.02, 0.05, 0.10, 0.20, 0.35, 0.50):
            add("benefit", "balanced20", 30, p, rho, replicates=500)
    return out


def lineage_effects(cell: dict, rng: np.random.Generator) -> np.ndarray:
    """Effects with variance rho on the liability scale, by theoretical moments."""
    groups, rho = cell["groups"], cell["rho"]
    law = cell["law"]
    if law in ("gaussian", "informative_sizes"):
        z = rng.standard_normal(groups)
    elif law == "t4":
        z = rng.standard_t(4, groups) / math.sqrt(2.0)
    elif law == "two_point":
        q = 0.2
        z = (rng.random(groups) < q).astype(float)
        z = (z - q) / math.sqrt(q * (1 - q))
    else:
        raise ValueError(law)
    return math.sqrt(rho) * z


def draw(cell: dict, rng: np.random.Generator):
    """Counts and sizes from a threshold model on normal liabilities.

    The liability of an isolate is its lineage effect plus an independent
    normal residual of variance 1 - rho, and the isolate is positive when the
    liability exceeds the (1 - p) quantile of a standard normal. This is the
    model's own statement, written without the package's likelihood.
    """
    sizes = np.asarray(cell["sizes"], dtype=int)
    effects = lineage_effects(cell, rng)
    if cell["law"] == "informative_sizes":
        # Larger lineages carry larger effects.
        sizes = np.sort(sizes)[np.argsort(np.argsort(effects))]
    cut = norm.ppf(1.0 - cell["prevalence"])
    residual_sd = math.sqrt(1.0 - cell["rho"])
    counts = np.array([int((effects[g] + residual_sd * rng.standard_normal(m) > cut).sum())
                       for g, m in enumerate(sizes)])
    return counts, sizes


def rng_for(cell: dict, replicate: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence(
        [SEEDS[cell["family"]], cell["cell"], replicate]))
