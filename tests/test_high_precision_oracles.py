"""The same quantities, recomputed independently to thirty digits.

A test that compares an implementation with another library checks that the
two agree, which is weaker than it looks when both use the same method. The
oracles here are written for the test alone: each one states the definition of
the quantity and evaluates it with `mpmath` at a precision where the
floating-point result cannot be the reference. Where a closed form in exact
rational arithmetic exists, that is used instead and the comparison is exact.

The formula an oracle uses is deliberately not the formula the package uses.
The bivariate normal distribution function, for instance, is computed here by
conditioning on the first coordinate, while the package integrates Plackett's
identity over the correlation: agreement then says the value is right, not
that the same series was summed twice.

Each tolerance is argued where it is written. None of them is a number chosen
until the test passed.
"""
from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
import pytest

from amr_clonalshare._normal_numerics import (normal_interval_logmass,
                                              normal_truncated_moments)
from amr_clonalshare.censored import marginal_loglik
from amr_clonalshare.clonality import _clopper_pearson
from amr_clonalshare.latent import bivariate_normal_cdf_equal
from amr_clonalshare.stats import fisher_exact_p

# A development dependency: the package itself never needs it, and a checkout
# without it skips this file rather than failing it.
mpmath = pytest.importorskip("mpmath")

DPS = 40


def _mp(x):
    """A float as an mpmath number, exactly.

    Not through `repr`: the shortest decimal that round-trips is not the same
    number as the float, and on an interval of width 1e-8 an endpoint moved by
    one unit in the last place changes the width by one part in 10^8. An
    oracle built that way reports an error the implementation does not have.
    `mpmath.mpf` converts a binary float without loss.
    """
    return mpmath.mpf(float(x))


# --------------------------------------------------------------------------- #
# Oracles
# --------------------------------------------------------------------------- #

def _mp_interval_mass(a, b):
    """P(a < Z <= b) at DPS digits."""
    with mpmath.workdps(DPS):
        lo = mpmath.mpf("-inf") if a == -math.inf else _mp(a)
        hi = mpmath.mpf("inf") if b == math.inf else _mp(b)
        return mpmath.ncdf(hi) - mpmath.ncdf(lo)


def _mp_truncated_moments(a, b, mean, sd):
    """Mean and variance of N(mean, sd^2) restricted to (a, b], by quadrature."""
    with mpmath.workdps(DPS):
        m, s = _mp(mean), _mp(sd)
        lo = mpmath.mpf("-inf") if a == -math.inf else _mp(a)
        hi = mpmath.mpf("inf") if b == math.inf else _mp(b)
        density = lambda x: mpmath.npdf(x, m, s)
        mass = mpmath.quad(density, [lo, hi])
        first = mpmath.quad(lambda x: x * density(x), [lo, hi]) / mass
        second = mpmath.quad(lambda x: (x - first) ** 2 * density(x), [lo, hi]) / mass
        return float(first), float(second)


def _mp_bivariate_equal(h, rho):
    """P(Z1 <= h, Z2 <= h) by conditioning on Z1, not by Plackett's identity."""
    with mpmath.workdps(DPS):
        hh, r = _mp(h), _mp(rho)
        if abs(r) == 1:
            return float(mpmath.ncdf(hh)) if r > 0 else float(
                max(mpmath.mpf(0), 2 * mpmath.ncdf(hh) - 1))
        root = mpmath.sqrt(1 - r * r)
        inner = lambda x: mpmath.npdf(x) * mpmath.ncdf((hh - r * x) / root)
        return float(mpmath.quad(inner, [mpmath.mpf("-inf"), hh]))


def _mp_marginal_loglik(lo, hi, code, G, grand, tau2, sigma):
    """The definition: one integral over the random intercept of each lineage.

        sum_g log int phi(u; 0, tau2) prod_i [Phi((hi-grand-u)/s)
                                              - Phi((lo-grand-u)/s)] du
    """
    with mpmath.workdps(DPS):
        s = _mp(sigma)
        mu = _mp(grand)
        t2 = _mp(tau2)
        total = mpmath.mpf(0)
        for g in range(G):
            members = [i for i, c in enumerate(code) if c == g]
            if not members:
                continue

            def product(u, members=members):
                out = mpmath.mpf(1)
                for i in members:
                    a, b = lo[i], hi[i]
                    upper = (mpmath.mpf("inf") if b == math.inf
                             else (_mp(b) - mu - u) / s)
                    lower = (mpmath.mpf("-inf") if a == -math.inf
                             else (_mp(a) - mu - u) / s)
                    out *= mpmath.ncdf(upper) - mpmath.ncdf(lower)
                return out

            if t2 == 0:
                total += mpmath.log(product(mpmath.mpf(0)))
                continue
            sd = mpmath.sqrt(t2)
            integrand = lambda u: mpmath.npdf(u, 0, sd) * product(u)
            # Ten standard deviations of the intercept hold the mass to well
            # beyond the working precision.
            value = mpmath.quad(integrand, [-10 * sd, 0, 10 * sd])
            total += mpmath.log(value)
        return float(total)


def _exact_fisher(n11, n10, n01, n00):
    """The hypergeometric two-sided p-value in exact rational arithmetic."""
    row1, row2 = n11 + n10, n01 + n00
    col1, total = n11 + n01, n11 + n10 + n01 + n00

    def weight(k):
        return Fraction(math.comb(row1, k) * math.comb(row2, col1 - k),
                        math.comb(total, col1))

    support = [k for k in range(max(0, col1 - row2), min(row1, col1) + 1)]
    observed = weight(n11)
    # The conventional two-sided rule: every table no more probable than the
    # observed one, with a rational comparison so no rounding decides it.
    return float(sum(weight(k) for k in support if weight(k) <= observed))


# --------------------------------------------------------------------------- #
# The comparisons
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("a, b", [
    (-1.0, 1.0), (0.0, 0.5), (-math.inf, -2.0), (3.0, math.inf),
    (7.0, 9.0),          # a far tail, where a difference of two Phi cancels
    (-9.0, -8.0),
    (1.0, 1.0 + 1e-7),   # a very short interval, the other cancellation
    (-40.0, -39.0),      # beyond the range where Phi itself is representable
])
def test_the_log_mass_of_an_interval_matches_the_definition(a, b):
    """log(Phi(b) - Phi(a)), recomputed at forty digits.

    Tolerance: the package computes this in float64 through log_ndtr and
    expm1, whose relative error is a few units in the last place; 1e-13
    relative on the logarithm is two orders of magnitude above that and far
    below anything that would hide a wrong branch.
    """
    got = float(normal_interval_logmass(a, b))
    with mpmath.workdps(DPS):
        want = float(mpmath.log(_mp_interval_mass(a, b)))
    assert got == pytest.approx(want, rel=1e-13, abs=1e-13)


@pytest.mark.parametrize("a, b, mean, sd", [
    (-1.0, 1.0, 0.0, 1.0),
    (-math.inf, 0.0, 0.5, 2.0),
    (2.0, math.inf, 0.0, 1.0),
    (4.0, 6.0, 0.0, 1.0),
    (-3.0, -2.5, 1.0, 0.5),
])
def test_the_truncated_moments_match_the_definition(a, b, mean, sd):
    """Tolerance: quadrature on both sides, so 1e-10 relative; the package's
    own accuracy target for these moments is the sixth digit."""
    got_mean, got_var = normal_truncated_moments(a, b, mean, sd)
    want_mean, want_var = _mp_truncated_moments(a, b, mean, sd)
    assert float(got_mean) == pytest.approx(want_mean, rel=1e-10, abs=1e-12)
    assert float(got_var) == pytest.approx(want_var, rel=1e-10, abs=1e-12)


@pytest.mark.parametrize("h", [-2.0, -0.5, 0.0, 0.75, 2.5])
@pytest.mark.parametrize("rho", [0.0, 0.05, 0.4, 0.9, 0.999])
def test_the_bivariate_normal_distribution_function_matches_the_definition(h, rho):
    """The package integrates Plackett's identity; the oracle conditions on
    the first coordinate. Tolerance: 1e-12 absolute, the accuracy the
    package's own docstring claims for its Gauss-Legendre rule."""
    assert bivariate_normal_cdf_equal(h, rho) == pytest.approx(
        _mp_bivariate_equal(h, rho), rel=0, abs=1e-12)


def _censored_fixture():
    rng = np.random.default_rng(4242)
    code = np.repeat(np.arange(4), [3, 4, 2, 3])
    y = rng.normal(size=4)[code] * 0.7 + rng.normal(size=code.size) * 0.9
    lo = np.floor(y)
    hi = lo + 1.0
    lo[y < -0.8] = -np.inf          # left censored at the lowest well
    hi[y < -0.8] = -0.8
    lo[y > 1.2] = 1.2               # right censored at the highest well
    hi[y > 1.2] = np.inf
    return lo, hi, code.astype(np.int64), 4


@pytest.mark.parametrize("grand, tau2, sigma", [
    (0.0, 0.5, 1.0),
    (0.3, 0.05, 0.8),
    (-0.4, 2.0, 1.5),
    (0.0, 0.0, 1.0),     # no between-lineage variance: the integral collapses
    (0.0, 1e-6, 0.6),    # a nearly degenerate intercept, the hard case for
                         # Gauss-Hermite nodes
])
def test_the_censored_likelihood_matches_its_own_definition(grand, tau2, sigma):
    """The headline estimator rests on this number.

    Tolerance: 1e-9 relative on a log likelihood of order ten, that is about
    1e-8 absolute. The package integrates with adaptive Gauss-Hermite and
    claims convergence to a tolerance of 1e-8; a disagreement above this is
    the integration failing, not the arithmetic.
    """
    lo, hi, code, G = _censored_fixture()
    got = marginal_loglik(lo, hi, code, G, grand, tau2, sigma)
    want = _mp_marginal_loglik(list(lo), list(hi), list(code), G,
                               grand, tau2, sigma)
    assert got == pytest.approx(want, rel=1e-9, abs=1e-8)


def test_the_censored_likelihood_matches_the_definition_under_heavy_censoring():
    """Every reading on an end well: the case a dilution panel produces and
    the one where a quadrature rule is most easily wrong."""
    code = np.array([0, 0, 0, 1, 1, 2, 2, 2], dtype=np.int64)
    lo = np.array([-np.inf, -np.inf, -np.inf, 2.0, 2.0, -np.inf, 2.0, 2.0])
    hi = np.array([-1.0, -1.0, -1.0, np.inf, np.inf, -1.0, np.inf, np.inf])
    got = marginal_loglik(lo, hi, code, 3, 0.2, 0.7, 1.1)
    want = _mp_marginal_loglik(list(lo), list(hi), list(code), 3, 0.2, 0.7, 1.1)
    assert got == pytest.approx(want, rel=1e-9, abs=1e-8)


@pytest.mark.parametrize("table", [
    (3, 1, 1, 3), (10, 2, 3, 15), (0, 5, 5, 0), (1, 0, 0, 1), (7, 7, 7, 7),
    (20, 1, 1, 20),
])
def test_the_exact_test_is_the_rational_hypergeometric_sum(table):
    """Exact on both sides: the comparison is to the last bit that float64
    can carry, so the tolerance is the representation error alone."""
    assert fisher_exact_p(*table) == pytest.approx(_exact_fisher(*table),
                                                   rel=1e-12, abs=1e-15)


@pytest.mark.parametrize("successes, trials", [(0, 20), (1, 20), (7, 20),
                                               (19, 20), (20, 20), (3, 1000)])
def test_the_exact_binomial_interval_solves_its_defining_equation(successes, trials):
    """A Clopper-Pearson limit is defined by the tail probability it makes
    equal to alpha/2, not by the Beta quantile that computes it. The oracle
    evaluates that tail at forty digits.

    Tolerance: 1e-10 on a probability, three orders below the width of any
    interval this test looks at.
    """
    alpha = 0.05
    low, high = _clopper_pearson(successes, trials, 1.0 - alpha)
    with mpmath.workdps(DPS):
        def upper_tail(p, k):
            return sum(mpmath.binomial(trials, j) * p ** j * (1 - p) ** (trials - j)
                       for j in range(k, trials + 1))

        def lower_tail(p, k):
            return sum(mpmath.binomial(trials, j) * p ** j * (1 - p) ** (trials - j)
                       for j in range(0, k + 1))

        if successes == 0:
            assert low == 0.0
        else:
            assert float(upper_tail(_mp(low), successes)) == \
                pytest.approx(alpha / 2, rel=0, abs=1e-10)
        if successes == trials:
            assert high == 1.0
        else:
            assert float(lower_tail(_mp(high), successes)) == \
                pytest.approx(alpha / 2, rel=0, abs=1e-10)
