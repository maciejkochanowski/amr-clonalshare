"""Adaptive quadrature for the one-dimensional Gaussian lineage effect.

Nodes are centred and scaled at the conditional mode. Successive node counts
must agree; difficult groups fall back to independently adaptive integration.
This evaluates a likelihood, not a confidence-interval calibration.
"""
from functools import lru_cache
import numpy as np
from scipy.integrate import quad
from scipy.special import logsumexp


@lru_cache(maxsize=4)
def _nodes(order):
    x, w = np.polynomial.hermite_e.hermegauss(order)
    return x, np.log(w) - 0.5*np.log(2*np.pi)


def adaptive_loglik(a, b, code, groups, grand, tau2, sigma, logmass, moments):
    """Evaluate grouped interval probabilities with adaptive Gaussian nodes."""
    counts = np.bincount(code, minlength=groups).astype(float)
    present = counts > 0
    exact = np.isfinite(a) & (a == b)
    variance = sigma*sigma
    def observation_loglik(effect):
        mu = grand + effect
        cell = logmass((a[:, None]-mu)/sigma, (b[:, None]-mu)/sigma)
        density = -0.5*((b[:, None]-mu)/sigma)**2 - np.log(sigma) - 0.5*np.log(2*np.pi)
        return np.where(exact[:, None], density, cell)
    def grouped(cell):
        return np.stack([np.bincount(code, weights=cell[:, k], minlength=groups)
                         for k in range(cell.shape[1])], axis=1)
    def joint(effect):
        return grouped(observation_loglik(effect[code, None]))[:, 0] - effect**2/(2*tau2)
    proxy = np.where(np.isfinite(a) & np.isfinite(b),
                     0.5*(np.where(np.isfinite(a), a, 0)+np.where(np.isfinite(b), b, 0)),
                     np.where(np.isfinite(a), a, np.where(np.isfinite(b), b, grand)))
    mode = np.bincount(code, weights=proxy-grand, minlength=groups)/np.maximum(counts, 1)
    mode *= tau2/(tau2 + variance/np.maximum(counts, 1))
    for _ in range(60):
        mean, var = moments(a, b, grand+mode[code], sigma)
        score = np.bincount(code, weights=(mean-grand-mode[code])/variance,
                            minlength=groups) - mode/tau2
        precision = 1/tau2 + np.bincount(code, weights=(1-np.clip(var/variance, 0, 1))/variance,
                                        minlength=groups)
        step = score/precision
        if np.max(np.abs(step[present])/(1+np.abs(mode[present]))) < 1e-10:
            break
        before = joint(mode)
        for _ in range(30):
            candidate = mode+step
            bad = joint(candidate) < before-1e-10
            if not bad.any():
                break
            step[bad] *= 0.5
        mode = candidate
    else:
        # Rounding can hold the Newton step above tolerance at the maximum of a
        # very narrow conditional density. The mode only centres the nodes, so
        # accept it when the last applied step is negligible on the scale of the
        # conditional spread; the node-count and error checks below still apply.
        if np.max(np.abs(step[present])*np.sqrt(precision[present])) > 1e-6:
            raise ArithmeticError("MIC quadrature: conditional mode did not converge")
    scale = 1/np.sqrt(precision)
    previous = None
    converged = ~present
    for order in (20, 40, 80, 160):
        x, logw = _nodes(order)
        u = mode[:, None] + scale[:, None]*x[None, :]
        terms = grouped(observation_loglik(u[code]))
        terms += logw[None, :] - u*u/(2*tau2) + 0.5*x[None, :]**2
        value = logsumexp(terms, axis=1) + np.log(scale/np.sqrt(tau2))
        if previous is not None:
            converged = (~present) | (np.abs(value-previous) <= 1e-10*(1+np.abs(value)))
            if converged.all():
                return float(value[present].sum())
        previous = value
    for g in np.flatnonzero(present & ~converged):
        mask = code == g
        aa, bb, ee = a[mask], b[mask], exact[mask]
        def log_integrand(z):
            u = mode[g] + scale[g]*z
            cell = logmass((aa-grand-u)/sigma, (bb-grand-u)/sigma)
            density = -0.5*((bb-grand-u)/sigma)**2 - np.log(sigma) - 0.5*np.log(2*np.pi)
            return float(np.where(ee, density, cell).sum() - u*u/(2*tau2))
        peak = log_integrand(0.0)
        integral, error = quad(lambda z: np.exp(log_integrand(z)-peak), -np.inf, np.inf,
                               epsabs=1e-11, epsrel=1e-10, limit=200)
        if integral <= 0 or error > 1e-7*integral:
            raise ArithmeticError("MIC quadrature failed its integration error check")
        value[g] = np.log(integral)+peak+np.log(scale[g])-0.5*np.log(2*np.pi*tau2)
    return float(value[present].sum())
