"""Check of the Gaussian law of lineage effects behind the population model.

Kept apart from ``_population_probit_numerics`` so that the numerical source
the population model's calibration and validation record stays unchanged.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.special import ndtr, ndtri

from ._population_probit_numerics import ProfileLikelihood, validate_counts


#: Points of the grid on which the unrestricted law of lineage probabilities
#: is fitted by the mixing check.
MIXING_GRID = 200
#: Simulated datasets of the mixing check.
MIXING_DRAWS = 99


def _npmle_loglik(k, m, iterations=400):
    """Log likelihood at the nonparametric maximum-likelihood law of the
    lineage probabilities, fitted by EM on a fixed grid (Laird 1978)."""
    from scipy.stats import binom
    grid = np.linspace(.5/MIXING_GRID, 1-.5/MIXING_GRID, MIXING_GRID)
    like = binom.pmf(k[:, None], m[:, None], grid[None, :])
    w = np.full(MIXING_GRID, 1/MIXING_GRID)
    for _ in range(iterations):
        w = w*(like/(like@w)[:, None]).mean(0)
    return float(np.log(like@w).sum())


def _local_nll(k, m, a, rho):
    """Minimum of the negative log likelihood over (intercept, rho), searched
    from a start near the maximum. The simulated datasets of the check are
    drawn from the fit, so their maxima lie near it; a local search stands in
    for the global scan of :func:`fit_profile`, and a maximum it missed can
    only enlarge the simulated statistic, which makes the check conservative."""
    from scipy.optimize import minimize
    ll = ProfileLikelihood(k, m)
    def f(x):
        r = 1/(1+math.exp(-x[1]))
        v = ll.nll(x[0], min(max(r, 1e-9), 1-1e-9))
        return v if np.isfinite(v) else 1e300
    r0 = min(max(rho, .01), .99)
    best = minimize(f, [a, math.log(r0/(1-r0))], method="Nelder-Mead",
                    options={"xatol": 1e-6, "fatol": 1e-9, "maxiter": 600})
    zero = ll.nll(float(ndtri(k.sum()/m.sum())), 0.)
    return float(min(best.fun, zero))


def mixing_check(counts, sizes, native, *, seed, draws=MIXING_DRAWS):
    """Check of the Gaussian law of lineage effects against an unrestricted one.

    The statistic is twice the log-likelihood gap between the nonparametric
    maximum-likelihood law of the lineage probabilities and the fitted probit
    model; it is calibrated by datasets simulated from the fit, with the
    lineage sizes held. Lineage effects that fall into a few classes, a
    resistant clone among susceptible lineages, give a large gap; there the
    fixed-cut-off interval loses its coverage, which is why the check is run.
    """
    k, m = validate_counts(counts, sizes)
    if not native.get("identified") or native.get("nll") is None:
        return None
    rho, prevalence = float(native["rho_hat"]), float(native["prevalence_hat"])
    a = float(ndtri(prevalence))
    observed = 2*(_npmle_loglik(k, m) + float(native["nll"]))
    rng = np.random.default_rng(seed)
    simulated = []
    for _ in range(int(draws)):
        u = rng.normal(0., np.sqrt(rho), len(m))
        p = ndtr((a+u)/np.sqrt(max(1-rho, 1e-12)))
        ks = rng.binomial(m, p)
        if ks.sum() in (0, m.sum()):
            continue
        try:
            nll = _local_nll(ks, m, a, rho)
        except (RuntimeError, ValueError, FloatingPointError, OverflowError):
            continue
        simulated.append(2*(_npmle_loglik(ks, m) + nll))
    if not simulated:
        return None
    simulated = np.asarray(simulated)
    pvalue = (1+int(np.sum(simulated >= observed-1e-9)))/(len(simulated)+1)
    return {"statistic": float(observed), "p_value": float(pvalue), "rejected": bool(pvalue <= .05),
            "simulated_datasets": int(len(simulated)),
            "rule": "nonparametric against Gaussian law of lineage effects, bootstrap p at most 0.05"}
