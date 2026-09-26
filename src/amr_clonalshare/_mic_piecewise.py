"""Full-real-line composite quadrature for sharply skewed lineage likelihoods.

An arctangent map covers the entire real line. Panel edges, transition scales
and the conditional mode define integration subintervals. Successive Legendre
orders check the normalising integral and analytic score expectations. This
is a numerical alternative, not a different statistical likelihood.
"""
from __future__ import annotations
from functools import lru_cache
import numpy as np
from scipy.special import roots_legendre, logsumexp
from ._normal_numerics import normal_interval_logmass, normal_truncated_moments


@lru_cache(maxsize=16)
def _rule(order):
    return roots_legendre(order)


def _group(a, b, n, exact, mode, scale, mean, sigma, tau2, strict):
    tau = np.sqrt(tau2)
    edges = np.unique(np.r_[a[np.isfinite(a)], b[np.isfinite(b)]])
    transitions = (edges[:, None]-mean+sigma*np.array([-8., -3., 0., 3., 8.])).ravel()
    anchors = np.r_[transitions, mode+scale*np.array([-4., -1., 0., 1., 4.]),
                    tau*np.array([-10., -4., -1., 0., 1., 4., 10.])]
    cuts = np.unique(np.r_[0., np.arctan((anchors-mode)/scale)/np.pi+.5, 1.])
    left, width = cuts[:-1], np.diff(cuts)
    previous = None
    tolerance = 2e-11 if strict else 2e-9
    for order in (12, 24, 48, 96, 192):
        nodes, weights = _rule(order)
        t = (left[:, None]+width[:, None]*(nodes+1)/2).ravel()
        angle = np.pi*(t-.5)
        u = mode+scale*np.tan(angle)
        eta = mean+u
        mass = normal_interval_logmass((a[:, None]-eta)/sigma, (b[:, None]-eta)/sigma)
        if exact.any():
            mass[exact] = (-.5*((a[exact, None]-eta)/sigma)**2
                           - np.log(sigma)-.5*np.log(2*np.pi))
        logweights = (np.log(width[:, None]*weights/2).ravel()+np.log(scale*np.pi)
                      -2*np.log(np.cos(angle))-.5*np.log(2*np.pi*tau2))
        terms = n@mass-u*u/(2*tau2)+logweights
        ll = float(logsumexp(terms))
        if not np.isfinite(ll):
            raise ArithmeticError('nonfinite composite lineage likelihood')
        # Avoid zero times overflowing score terms outside floating-point mass.
        live = terms-ll > -700.
        posterior = np.exp(terms[live]-ll)
        ey, vy = normal_truncated_moments(a[:, None], b[:, None], eta[None, live], sigma)
        delta = ey-eta[None, live]
        per_bin = n*(delta@posterior)/sigma**2
        value = np.array([ll, np.dot(posterior, n@delta/sigma**2),
            np.dot(posterior, n@((vy+delta*delta)/sigma**2-1)),
            np.dot(posterior, u[live]**2/tau2-1)])
        if not np.isfinite(value).all():
            raise ArithmeticError('nonfinite composite lineage score')
        if previous is not None:
            errors = np.abs(value-previous)
            relative_score = np.max(errors[1:]/np.maximum(1., np.abs(value[1:])))
            if errors[0] <= tolerance and relative_score <= tolerance:
                return value, float(max(errors[0], relative_score)), per_bin
        previous = value
    raise ArithmeticError('composite lineage likelihood or score failed precision check')


def integrate_lineages(problem, counts, modes, scales, mean, sigma, tau2, *, strict=False):
    """Return log likelihood, mean/log-sigma/log-tau scores, error and the
    mean score of every bin (for the fixed effects), with error check."""
    total = np.zeros(4)
    error = 0.
    per_bin = np.zeros(len(problem.a))
    for nn, mode, scale in zip(counts, modes, scales):
        present = nn > 0
        value, local_error, local_bin = _group(problem.a[present], problem.b[present], nn[present],
            problem.exact[present], mode, scale, mean, sigma, tau2, strict)
        total += value
        error += local_error
        per_bin[present] += local_bin
    return (*map(float, total), float(error), per_bin)
