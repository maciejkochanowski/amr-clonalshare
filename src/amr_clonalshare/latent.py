"""The clonal share of a binary trait restated on the latent scale.

A susceptibility call is read on the observed scale: the share is the
intraclass correlation of the 0/1 indicator, the scale a prevalence table uses,
and the binary-interval validation is on that scale. A Gaussian probit
threshold model also defines the intraclass correlation of an unobserved
continuous liability whose threshold crossing is the call. The map below
uses that model; it is not the latent ICC of every binomial-link model.

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
known, fixed prevalence for the matching population parameter. The package
passes the estimated prevalence; treating that estimate as fixed does not
establish conditional or unconditional coverage. Moreover, the image of a
finite collection's lineage-membership share is a Gaussian-model-equivalent
summary, not that collection's actual latent variance share. The collection
and population intervals must retain their distinct targets.

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

    A monotone map preserves coverage for a matching target at a known,
    fixed prevalence. When prevalence is estimated from the same clustered
    sample, this plug-in transformation does not propagate its uncertainty
    or covariance with the observed-scale share, and is not a calibrated
    conditional confidence procedure. The returned endpoints are therefore
    an exploratory model-equivalent restatement. A population latent ICC
    additionally requires the Gaussian probit model and a population-scale
    observed interval; transforming a collection interval does not change
    its target into a population parameter.
    """
    return latent_share(low, prevalence), latent_share(high, prevalence)
