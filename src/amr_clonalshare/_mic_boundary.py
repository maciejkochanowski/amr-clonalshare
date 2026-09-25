"""Likelihood suprema when every lineage occupies one observation interval.

For fixed mean and total scale, the event that all members of a lineage
fall in the same interval is contained in the event for one member. The
upper bound is attained as residual variance tends to zero (rho tends to
one). Thus the unrestricted supremum reduces to independent, interval-
censored lineage values, with each lineage counted once, not n_g times.
Degenerate nuisance limits are reported explicitly and are never simulated
as though they were finite fitted Gaussian distributions.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize
from ._normal_numerics import normal_interval_logmass, normal_truncated_moments


def monomorphic_boundary_fit(problem):
    """Return the rho=1 limiting maximum, or None when the shortcut is invalid."""
    from ._mic_likelihood import MICFit
    if problem.all_exact or problem.exact.any():
        return None
    occupied = problem.counts > 0
    if not np.all(occupied.sum(axis=1) == 1):
        return None
    counts = occupied.sum(axis=0).astype(float)
    used = counts > 0
    a, b, counts = problem.a[used], problem.b[used], counts[used]
    groups = float(counts.sum())
    if len(counts) == 1:
        mean = float((a[0]+b[0])/2) if np.isfinite(a[0]+b[0]) else problem.center
        return MICFit(mean, 0., 1., 0., True, 0, 0, 0.,
            'rho=1 degenerate-scale supremum; no finite nuisance maximum')
    if len(counts) == 2:
        ll = float(np.dot(counts, np.log(counts/groups)))
        if np.isneginf(a[0]) and np.isposinf(b[1]) and b[0] < a[1]:
            return MICFit(float('nan'), float('inf'), 1., ll, True, 0, 0, 0.,
                'rho=1 infinite-scale binomial supremum; no finite nuisance maximum')
        if b[0] == a[1]:
            return MICFit(float(b[0]), 0., 1., ll, True, 0, 0, 0.,
                'rho=1 adjacent-bin zero-scale supremum; no finite nuisance maximum')
    unit, center = problem.unit, problem.center
    aa, bb = (a-center)/unit, (b-center)/unit
    evaluations = 0
    def objective(theta):
        nonlocal evaluations
        evaluations += 1
        mean, sd = float(theta[0]), float(np.exp(theta[1]))
        logmass = normal_interval_logmass((aa-mean)/sd, (bb-mean)/sd)
        ey, vy = normal_truncated_moments(aa, bb, mean, sd)
        delta = ey-mean
        gradient = np.array([np.dot(counts, delta)/sd**2,
            np.dot(counts, (vy+delta**2)/sd**2-1)])
        return -float(np.dot(counts, logmass)), -gradient
    fits = [minimize(objective, start, jac=True, method='L-BFGS-B',
        bounds=[(-100., 100.), (-14., 14.)],
        options={'ftol': 1e-13, 'gtol': 1e-8, 'maxiter': 300, 'maxls': 40})
        for start in ([0., 0.], [0., np.log(2.)])]
    best = min(fits, key=lambda fit: fit.fun)
    edge = abs(best.x[0]) > 99.999 or abs(best.x[1]) > 13.999
    ok = bool(best.success and np.isfinite(best.fun) and not edge)
    return MICFit(float(center+unit*best.x[0]), float(unit*np.exp(best.x[1])),
        1., float(-best.fun), ok, evaluations, 0, 0.,
        'rho=1 limiting grouped-normal maximum' if ok else
        'unresolved limiting grouped-normal nuisance maximum')
