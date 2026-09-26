"""Check of the Gaussian MIC model against the dilution readings.

The calibrated interval assumes that the latent log2 MIC is Gaussian within
and between lineages. Two statistics compare the recorded readings with what
the fitted model predicts:

``marginal``    the deviance of the readings of every panel and covariate level
                over that panel's dilution intervals, against the fitted normal
                law of a reading; it sees a second mode or a heavy tail in the
                pooled distribution.
``conditional`` the Pearson statistic of every lineage's readings over the
                dilution intervals, against their prediction given the rest of
                the lineage (the posterior of its effect); it sees readings
                that the lineage's own centre and the residual scale do not
                explain.

Both are calibrated by the parametric bootstrap that the null-wise test runs
at the maximum-likelihood estimate, so they cost no extra fits; the check
rejects when twice the smaller p-value is at most 0.05. Exact readings take
no part. The check says whether the readings look like the model; it does
not say which other model holds.
"""
from __future__ import annotations

import numpy as np
from scipy.special import ndtr

from ._mic_panels import _record_offset

#: Points of the grid that carries the posterior of a lineage effect.
GRID = 201
#: Level of the check: rejected when min(p_marginal, p_conditional) <= LEVEL / 2.
LEVEL = .05


def _categories(lo, hi, index, edges):
    """Dilution interval of every reading within its own panel (0..len(edges))."""
    cat = np.empty(len(lo), dtype=int)
    for k, e in enumerate(edges):
        m = index == k
        cat[m] = np.searchsorted(e, np.where(np.isfinite(hi[m]), hi[m], np.inf), side='left')
    return cat


def statistics(problem, fit, panels, exact):
    """The two statistics of one dataset under one fit; None when no reading is an interval."""
    edges, index = panels
    use = ~exact
    if not use.any() or not np.isfinite(fit.total_sd) or fit.total_sd <= 0:
        return None
    lo, hi = problem.lo, problem.hi
    mu = fit.mean + np.asarray(_record_offset(problem, fit)) * np.ones(problem.n)
    sd, rho = fit.total_sd, min(max(fit.rho, 0.), 1.)
    cat = _categories(lo, hi, index, edges)
    level = problem.record_level
    key = np.column_stack([index, level]) if level.size else index[:, None]
    marginal = 0.
    groups, inverse = np.unique(key[use], axis=0, return_inverse=True)
    inverse = inverse.ravel()
    rows = np.flatnonzero(use)
    for g in range(len(groups)):
        m = rows[inverse == g]
        e = edges[groups[g][0]]
        p = np.diff(np.r_[0., ndtr((e - mu[m[0]]) / sd), 1.])
        obs = np.bincount(cat[m], minlength=len(p)).astype(float)
        exp = len(m) * np.maximum(p, 1e-300)
        nz = obs > 0
        marginal += float(2 * np.sum(obs[nz] * np.log(obs[nz] / exp[nz])))
    tau, sig = sd * np.sqrt(rho), sd * np.sqrt(max(1 - rho, 0.))
    if sig <= 0:
        return marginal, 0.
    u = np.linspace(-7 * tau, 7 * tau, GRID) if tau > 1e-8 * sd else np.zeros(1)
    log_prior = -.5 * (u / tau) ** 2 if tau > 1e-8 * sd else np.zeros(1)
    conditional = 0.
    code = problem.code
    for g in range(problem.groups):
        m = np.flatnonzero((code == g) & use)
        if len(m) < 2:
            continue
        eta = mu[m][:, None] + u[None, :]
        mass = ndtr((hi[m][:, None] - eta) / sig) - ndtr((lo[m][:, None] - eta) / sig)
        log_w = np.log(np.clip(mass, 1e-300, None)).sum(0) + log_prior
        w = np.exp(log_w - log_w.max())
        w /= w.sum()
        cell = key[m]
        cells, cinv = np.unique(cell, axis=0, return_inverse=True)
        cinv = cinv.ravel()
        for c in range(len(cells)):
            members = m[cinv == c]
            e = edges[cells[c][0]]
            cdf = ndtr((e[None, :, None] - mu[members][:, None, None] - u[None, None, :]) / sig)
            pk = np.diff(np.concatenate([np.zeros((len(members), 1, len(u))), cdf,
                                         np.ones((len(members), 1, len(u)))], axis=1), axis=1)
            expected = np.maximum((pk @ w).sum(0), 1e-12)
            observed = np.bincount(cat[members], minlength=len(e) + 1).astype(float)
            conditional += float(np.sum((observed - expected) ** 2 / expected))
    return marginal, conditional


def check(observed, simulated):
    """p-values of the two statistics and the decision, from the simulated datasets."""
    if observed is None:
        return None
    sims = np.array([s for s in simulated if s is not None], dtype=float).reshape(-1, 2)
    n = len(sims)
    if n == 0:
        return None
    p = [(1 + int(np.sum(sims[:, j] >= observed[j] - 1e-12))) / (n + 1) for j in (0, 1)]
    combined = min(1., 2 * min(p))
    return {"marginal_deviance": float(observed[0]), "conditional_pearson": float(observed[1]),
            "p_marginal": float(p[0]), "p_conditional": float(p[1]), "p_value": float(combined),
            "rejected": bool(combined <= LEVEL), "simulated_datasets": int(n),
            "rule": "twice the smaller of the two bootstrap p-values, at most 0.05"}
