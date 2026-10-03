"""The same quantities, recomputed independently to thirty digits.

A test that compares an implementation with another library checks that the
two agree, which is weaker than it looks when both use the same method. The
oracles here are written for the test alone: each one states the definition of
the quantity and evaluates it with `mpmath` at a precision where the
floating-point result cannot be the reference. Where a closed form in exact
rational arithmetic exists, that is used instead and the comparison is exact.

The formula an oracle uses is deliberately not the formula the package uses:
agreement then says the value is right, not that the same series was summed
twice.

Each tolerance is argued where it is written. None of them is a number chosen
until the test passed.
"""
from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
import pytest

from amr_clonalshare.clonality import _clopper_pearson
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


# --------------------------------------------------------------------------- #
# The share of the lineages in hand and its bootstrap
# --------------------------------------------------------------------------- #

def _held_out_numerator(values, code, fold):
    """The corrected numerator of one fold draw from its definition, in exact
    rational arithmetic: every held-out isolate predicted by the mean of its
    lineage in training, less the variance of that mean, the lineage's
    unbiased variance over its training count."""
    total = Fraction(0)
    lineages = sorted(set(code))
    variance = {}
    for g in lineages:
        xs = [values[i] for i in range(len(code)) if code[i] == g]
        mean = sum(xs, Fraction(0)) / len(xs)
        variance[g] = sum(((x - mean) ** 2 for x in xs), Fraction(0)) / (len(xs) - 1)
    for f in set(fold):
        for i in range(len(code)):
            if fold[i] != f:
                continue
            train = [values[j] for j in range(len(code)) if fold[j] != f and code[j] == code[i]]
            mean = sum(train, Fraction(0)) / len(train)
            total += (values[i] - mean) ** 2 - variance[code[i]] / len(train)
    return total


def test_over_every_dealing_the_corrected_numerator_is_the_within_lineage_sum_of_squares():
    """The identity behind the fold average: over every dealing of the folds,
    equally likely, the corrected held-out error of a lineage of m isolates
    averages to m times its unbiased variance. Every priority order of six
    isolates and every first fold is enumerated, in rational arithmetic, and
    the package's closed form is compared with the result."""
    import itertools
    from amr_clonalshare.attribution import _dealt_folds, _fold_average

    values = [Fraction(v) for v in (0, 1, 3, 2, 2, 7)]
    code = np.array([0, 0, 0, 1, 1, 1])
    k = 2

    class Fixed:
        def __init__(self, priority, start):
            self.priority, self.start = priority, start

        def permutation(self, n):
            return np.asarray(self.priority)

        def integers(self, k):
            return self.start

    total, count = Fraction(0), 0
    for priority in itertools.permutations(range(6)):
        for start in range(k):
            fold = _dealt_folds(code, k, Fixed(priority, start))
            total += _held_out_numerator(values, list(code), list(fold))
            count += 1
    average = total / count
    expected = Fraction(0)
    for g in (0, 1):
        xs = [values[i] for i in range(6) if code[i] == g]
        mean = sum(xs, Fraction(0)) / 3
        expected += 3 * sum(((x - mean) ** 2 for x in xs), Fraction(0)) / 2   # m s^2, s^2 = SS / (m - 1)
    assert average == expected
    # the closed form uses the same numerator: with a ratio of one it is
    # 1 - numerator / SST
    y = np.array([float(v) for v in values])
    sst = float(sum((v - sum(values, Fraction(0)) / 6) ** 2 for v in values))
    assert _fold_average(y, code, 1.) == pytest.approx(1 - float(expected) / sst, rel=1e-14)


def test_the_share_of_a_law_is_its_between_lineage_share_of_the_variance():
    """eta of a smoothed law from its definition, in rational arithmetic:
    every isolate of a cell draws value v_k with probability law[c, k]; the
    share is the variance of the lineage means over the total variance."""
    from amr_clonalshare.attribution import _share_of
    law = [[Fraction(1, 2), Fraction(1, 3), Fraction(1, 6)],
           [Fraction(1, 5), Fraction(0), Fraction(4, 5)],
           [Fraction(1, 4), Fraction(1, 4), Fraction(1, 2)]]
    values = [Fraction(-1), Fraction(1, 2), Fraction(3)]
    cell_of = [0, 0, 0, 1, 1, 2, 2, 2, 2]          # cells 0, 1 in lineage 0; cell 2 is lineage 1
    cell_lineage = [0, 0, 1]
    n = len(cell_of)
    mean = [sum((p * v for p, v in zip(row, values)), Fraction(0)) for row in law]
    second = [sum((p * v * v for p, v in zip(row, values)), Fraction(0)) for row in law]
    size = [sum(1 for c in cell_of if cell_lineage[c] == g) for g in (0, 1)]
    mu = [sum((mean[c] for c in cell_of if cell_lineage[c] == g), Fraction(0)) / size[g] for g in (0, 1)]
    grand = sum((mean[c] for c in cell_of), Fraction(0)) / n
    total = sum((second[c] for c in cell_of), Fraction(0)) / n - grand ** 2
    between = sum((Fraction(size[g], n) * (mu[g] - grand) ** 2 for g in (0, 1)), Fraction(0))
    got = _share_of(np.array([float(v) for v in values]), np.array(size, dtype=float),
                    np.array([[float(p) for p in row] for row in law]), np.array(cell_of),
                    np.array(cell_lineage))
    assert got == pytest.approx(float(between / total), rel=1e-14)


# --------------------------------------------------------------------------- #
# The bounds on the unobserved ordering
# --------------------------------------------------------------------------- #

def _rational_pieces(counts):
    """p, the reading ranges and the vertex of an order, in fractions, for
    readings of one panel (the ranges are the cumulative reading shares)."""
    G, B = len(counts), len(counts[0])
    n = sum(sum(r) for r in counts)
    P = [Fraction(sum(counts[g][r] for g in range(G)), n) for r in range(B)]
    H = [sum(P[:r + 1], Fraction(0)) for r in range(B)]
    L = [H[r] - P[r] for r in range(B)]
    w = [Fraction(sum(counts[g]), n) for g in range(G)]

    def vertex(order):
        x = [Fraction(0)] * G
        for r in range(B):
            if P[r] == 0:
                continue
            top = H[r]
            for g in order:                                      # top slice first
                share = Fraction(counts[g][r], n) / P[r]         # pi_gr
                x[g] += P[r] * share * (top - (H[r] - L[r]) * share / 2)
                top -= (H[r] - L[r]) * share
        return x
    return w, vertex, L, H, P


@pytest.mark.parametrize("counts", [
    [[3, 1, 0], [0, 2, 2], [1, 1, 3]],
    [[2, 0, 1, 1], [1, 2, 0, 1], [0, 1, 3, 0], [1, 1, 1, 1]],
    [[4, 1], [1, 4], [2, 2], [0, 3], [3, 0]],
])
def test_the_exact_upper_end_is_the_largest_share_over_every_order_in_rational_arithmetic(counts):
    import itertools
    from amr_clonalshare.latent_order import _Polytope, _upper_exact
    w, vertex, L, H, _ = _rational_pieces(counts)
    best = max(sum((x * x / wg for x, wg in zip(vertex(o), w)), Fraction(0))
               for o in itertools.permutations(range(len(counts))))
    poly = _Polytope(np.array(counts, dtype=float), np.array([float(v) for v in L]),
                     np.array([float(v) for v in H]))
    assert _upper_exact(poly) == pytest.approx(float(best), rel=1e-13)


@pytest.mark.parametrize("counts", [
    [[3, 1, 0], [0, 2, 2]],
    [[5, 0], [2, 3]],
    [[1, 1, 1, 1], [0, 0, 3, 1]],
    [[2, 2], [2, 2]],
])
def test_for_two_lineages_the_dual_certificate_meets_the_guaranteed_share_exactly(counts):
    """Two lineages on one panel: rho_min in closed form (x_1 is the weighted
    share of F(V) clipped to its range [l_1, u_1]) and the largest dual
    certificate, 12 G(e)_+^2 / Q(e) over the two directions, both in rational
    arithmetic, are equal (Theorem 2)."""
    w, vertex, L, H, P = _rational_pieces(counts)
    n = sum(sum(r) for r in counts)
    top, bottom = vertex([0, 1])[0], vertex([1, 0])[0]            # u_1 and l_1
    total = sum(vertex([0, 1]), Fraction(0))
    x1 = min(max(w[0] * total, bottom), top)
    rho_min = 12 * (x1 ** 2 / w[0] + (total - x1) ** 2 / w[1] - Fraction(1, 4))

    def dual(g, h):
        # P(Y_h wholly below Y_g) from the reading counts, one panel
        below = sum((Fraction(counts[h][a] * counts[g][b], sum(counts[h]) * sum(counts[g]))
                     for a in range(len(P)) for b in range(len(P)) if a < b), Fraction(0))
        omega = Fraction(sum(counts[g]) * sum(counts[h]), n * n)
        value = omega * (below - Fraction(1, 2))
        return 12 * max(value, Fraction(0)) ** 2 / (w[0] * w[1])   # Q(e) for z = e_g, centred
    assert max(dual(0, 1), dual(1, 0)) == rho_min
