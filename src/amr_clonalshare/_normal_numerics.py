"""Normal interval probabilities and moments without tail subtraction loss.

Ordinary intervals use reflected log-CDF differences. Short intervals are
integrated on their own scale. Far-tail moments are evaluated as moments of
the offset from the finite endpoint, so their variance is not the difference
of two numbers of order endpoint**2. These are numerical identities, not a
statistical calibration of any estimator using them.
"""
from __future__ import annotations

import numpy as np
from scipy.special import log_ndtr, logsumexp

_LOG2PI = float(np.log(2*np.pi))
_GL_X, _GL_W = np.polynomial.legendre.leggauss(48)
_GL_LOGW = np.log(_GL_W)
_LAG_X, _LAG_W = np.polynomial.laguerre.laggauss(64)


def normal_interval_logmass(a, b):
    """Return log P(a < Z <= b), including zero-width and unbounded intervals."""
    a, b = np.broadcast_arrays(np.asarray(a, dtype=float), np.asarray(b, dtype=float))
    reflected = a > 0
    left = np.where(reflected, -b, a)
    right = np.where(reflected, -a, b)
    la, lb = log_ndtr(left), log_ndtr(right)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        result = np.asarray(lb + np.log(-np.expm1(np.minimum(la-lb, 0.))))
        width = b-a
        middle = a + width/2
        short = (np.isfinite(a) & np.isfinite(b) & (width > 0)
                 & (width*(1+np.abs(middle)) < 1e-3))
    if np.any(short):
        half = width[short]/2
        z = middle[short, None] + half[:, None]*_GL_X
        result[short] = (np.log(half) + logsumexp(_GL_LOGW - .5*z*z, axis=1)
                         - .5*_LOG2PI)
    return result


def _bounded_offsets(lower, width):
    """Mean and variance of an offset on [0,width] with normal density."""
    delta = width[:, None]*(1+_GL_X)/2
    logw = _GL_LOGW - lower[:, None]*delta - .5*delta*delta
    weight = np.exp(logw-logsumexp(logw, axis=1)[:, None])
    mean = np.sum(weight*delta, axis=1)
    variance = np.sum(weight*(delta-mean[:, None])**2, axis=1)
    return mean, variance


def normal_truncated_moments(a, b, mean, sd):
    """Mean and variance of N(mean,sd**2) restricted to (a,b].

    Equal finite endpoints denote an exact observation. Infinite endpoints
    are allowed only on the corresponding open side. Public callers validate
    interval ordering; no positive-width interval is silently made exact.
    """
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError('normal standard deviation must be positive and finite')
    a, b, mean = np.broadcast_arrays(np.asarray(a, dtype=float), np.asarray(b, dtype=float),
                                     np.asarray(mean, dtype=float))
    shape = a.shape
    a, b, mean = a.ravel(), b.ravel(), mean.ravel()
    al, be = (a-mean)/sd, (b-mean)/sd
    exact = np.isfinite(a) & (a == b)
    logmass = normal_interval_logmass(al, be)
    finite_a, finite_b = np.isfinite(al), np.isfinite(be)
    aa, bb = np.where(finite_a, al, 0.), np.where(finite_b, be, 0.)
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        pa = np.where(finite_a & ~exact, np.exp(-.5*aa*aa-.5*_LOG2PI-logmass), 0.)
        pb = np.where(finite_b & ~exact, np.exp(-.5*bb*bb-.5*_LOG2PI-logmass), 0.)
        standard_mean = pa-pb
        standard_var = 1+aa*pa-bb*pb-standard_mean**2
        width = be-al
        short = (finite_a & finite_b & ~exact & (width > 0)
                 & (width*(1+np.abs(al+width/2)) < .1))
    # Direct centred moments avoid cancellation in a very short interval.
    if np.any(short):
        offset, var = _bounded_offsets(al[short], width[short])
        standard_mean[short], standard_var[short] = al[short]+offset, var
    # In a far tail the distribution is concentrated at its finite endpoint.
    for negative in (False, True):
        lower = -be if negative else al
        upper = -al if negative else be
        tail = (lower >= 8) & ~exact & (upper > lower) & ~short
        if not np.any(tail):
            continue
        indices = np.flatnonzero(tail)
        low, high = lower[tail], upper[tail]
        offset, variance = np.empty(low.size), np.empty(low.size)
        unbounded = np.isposinf(high)
        if np.any(unbounded):
            scale = low[unbounded]
            # x = lower*(Z-lower); the leading density is exp(-x).
            weight = _LAG_W * np.exp(-.5*(_LAG_X/scale[:, None])**2)
            weight /= weight.sum(axis=1)[:, None]
            ex = np.sum(weight*_LAG_X, axis=1)
            offset[unbounded] = ex/scale
            variance[unbounded] = np.sum(weight*(_LAG_X-ex[:, None])**2, axis=1)/scale**2
        bounded = ~unbounded
        if np.any(bounded):
            # Omitting offsets with lower*offset > 50 loses < exp(-50)
            # of the leading exponential tail; it prevents a broad quadrature
            # interval from missing the boundary peak.
            span = np.minimum(high[bounded]-low[bounded], 50/low[bounded])
            offset[bounded], variance[bounded] = _bounded_offsets(low[bounded], span)
        standard_mean[indices] = (-1 if negative else 1)*(low+offset)
        standard_var[indices] = variance
    result_mean = mean + sd*standard_mean
    result_var = sd*sd*np.clip(standard_var, 0., 1.)
    result_mean[exact], result_var[exact] = a[exact], 0.
    return result_mean.reshape(shape), result_var.reshape(shape)
