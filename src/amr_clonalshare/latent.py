"""The clonal share of a binary trait restated on the latent scale.

A susceptibility call is read on the observed scale: the share is the
intraclass correlation of the 0/1 indicator, the scale a prevalence table uses,
and every interval in this package is measured on it. A mixed model with a
binomial link reports a different number for the same cohort, the intraclass
correlation of an unobserved continuous liability whose threshold crossing is
the call, and a reader who fits ``rptR`` or ``lme4`` will want the two side by
side.

The two are one-to-one under the threshold model. Let the liability of isolate
``i`` in lineage ``g`` be ``a_g + e_ig`` with ``a_g ~ N(0, rho_L)`` and
``e_ig ~ N(0, 1 - rho_L)``, and let the call be ``1`` where the liability
exceeds ``t = Phi^{-1}(1 - p)`` for prevalence ``p``. Two isolates of one
lineage then share a bivariate normal liability with correlation ``rho_L``, and
the observed-scale correlation of their calls is

    rho_obs = (Phi_2(-t, -t; rho_L) - p^2) / (p (1 - p)),

where ``Phi_2`` is the standard bivariate normal distribution function. The map
is increasing in ``rho_L`` with ``rho_obs(0) = 0`` and ``rho_obs(1) = 1``, so it
is inverted by bisection, and the endpoints of an interval carry across at a
fixed prevalence, their coverage unchanged by a monotone map. The package
passes the observed prevalence and treats it as fixed, so the latent interval
is conditional on it; ``latent_interval`` states how far that assumption
reaches.

``Phi_2(h, h; rho)`` is computed as ``Phi(h)^2`` plus the integral of the
bivariate density over the correlation from zero to ``rho``, which after the
substitution ``r = sin(theta)`` is a bounded smooth integrand,

    Phi_2(h, h; rho) = Phi(h)^2 + (1 / 2 pi) int_0^{arcsin rho} exp(-h^2 / (1 + sin theta)) d theta,

so a Gauss-Legendre rule of modest order gives it to machine precision. At
``h = 0`` the integral is ``arcsin(rho) / 2 pi`` exactly, which is the check
the tests use.

This is the transformation Dempster and Lerner (1950, Genetics 35:212) and
Robertson (in the same paper's appendix) derived for heritability on the
liability scale, computed exactly rather than by their first-order
approximation, which overshoots one at large shares.
"""
from __future__ import annotations

from math import asin, erf, sqrt

import numpy as np
from scipy.special import ndtri as _ndtri

__all__ = ["bivariate_normal_cdf_equal", "observed_share", "latent_share",
           "latent_interval"]

_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(48)


def _phi(h: float) -> float:
    return 0.5 * (1.0 + erf(h / sqrt(2.0)))


def bivariate_normal_cdf_equal(h: float, rho: float) -> float:
    """``P(Z_1 <= h, Z_2 <= h)`` for standard normals with correlation ``rho``."""
    rho = float(np.clip(rho, -1.0, 1.0))
    base = _phi(h) ** 2
    if rho == 0.0:
        return base
    upper = asin(rho)
    theta = 0.5 * upper * (_NODES + 1.0)
    integrand = np.exp(-h * h / (1.0 + np.sin(theta)))
    integral = 0.5 * upper * float(np.sum(_WEIGHTS * integrand))
    return float(base + integral / (2.0 * np.pi))


def observed_share(latent: float, prevalence: float) -> float:
    """The observed-scale share a latent share ``latent`` produces at
    ``prevalence`` under the threshold model."""
    p = float(prevalence)
    if not 0.0 < p < 1.0:
        return float("nan")
    if latent <= 0.0:
        return 0.0
    if latent >= 1.0:
        return 1.0
    t = float(_ndtri(1.0 - p))
    joint = bivariate_normal_cdf_equal(-t, latent)
    return float(np.clip((joint - p * p) / (p * (1.0 - p)), 0.0, 1.0))


def latent_share(observed: float, prevalence: float, *, tol: float = 1e-10) -> float:
    """The latent-scale share whose threshold-model image is ``observed``.

    A share at or below zero has no latent counterpart and maps to zero; a
    share at or above one maps to one. Outside ``0 < prevalence < 1`` the
    threshold is not defined and NaN is returned.
    """
    p = float(prevalence)
    if not np.isfinite(observed) or not 0.0 < p < 1.0:
        return float("nan")
    if observed <= 0.0:
        return 0.0
    if observed >= 1.0:
        return 1.0
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if observed_share(mid, p) < observed:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return float(0.5 * (lo + hi))


def latent_interval(low: float, high: float, prevalence: float):
    """Both endpoints carried to the latent scale at the given prevalence.

    The map is monotone, so at a fixed prevalence the coverage of the interval
    is unchanged by the transformation. The prevalence the package passes is
    the observed one, treated as fixed: its own sampling error is not carried
    into the latent endpoints. Moving the prevalence by two standard errors
    moves a latent share of 0.77 by 0.001 on a cohort of 6,915 isolates at
    8.5 % prevalence and a latent share of 0.52 by 0.008 on 677 isolates at
    23 %; on a cohort of a hundred isolates the shift reaches 0.03, and the
    latent interval is then to be read as conditional on the observed
    prevalence.
    """
    return latent_share(low, prevalence), latent_share(high, prevalence)
