"""Properties the arithmetic must satisfy for every input, not for the inputs
someone thought of.

Each test states an identity or an order relation and lets Hypothesis search
for a counterexample. Where the property is an equation, the right-hand side
is computed in exact rational arithmetic with :mod:`fractions`, so the test is
an oracle rather than a second copy of the floating-point code. Where the
property is a bound or a monotonicity, the check is exact by construction.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import comb

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from amr_clonalshare import attribution
from amr_clonalshare.archephy import _exact_upper_tail
from amr_clonalshare.clonality import _clopper_pearson, _kitagawa, _shares_and_rates
from amr_clonalshare.inference import merge_pvalues
from amr_clonalshare.stats import (benjamini_hochberg, binary_distance_matrix,
                                   effective_dimension, jaccard_distance_matrix,
                                   permutation_pvalue, uninformative_rows)

# --------------------------------------------------------------------------- #
# Prevalence decomposition: the two components sum to the difference exactly
# --------------------------------------------------------------------------- #

# A collection is a list of (isolates, positives) per lineage; a lineage with
# zero isolates is absent from that collection, which is the case the
# rate-filling convention exists for.
_collection = st.lists(
    st.tuples(st.integers(0, 12), st.integers(0, 12)).map(
        lambda t: (t[0], min(t[1], t[0]))),
    min_size=1, max_size=6)


def _expand(counts):
    codes, y = [], []
    for code, (n, k) in enumerate(counts):
        codes.extend([code] * n)
        y.extend([1] * k + [0] * (n - k))
    return np.asarray(codes, dtype=np.intp), np.asarray(y, dtype=float)


@given(_collection, _collection)
def test_kitagawa_components_sum_to_the_prevalence_difference(a, b):
    """composition + within == prevalence(A) - prevalence(B), as rationals.

    Each share and each rate is a ratio of integers, so the identity can be
    checked exactly. Absent lineages take the other collection's rate, which
    is the convention the documentation states, and the identity must survive
    it.
    """
    L = max(len(a), len(b))
    a = a + [(0, 0)] * (L - len(a))
    b = b + [(0, 0)] * (L - len(b))
    na, nb = sum(n for n, _ in a), sum(n for n, _ in b)
    if na == 0 or nb == 0:
        return
    ca, ya = _expand(a)
    cb, yb = _expand(b)
    wa, pa = _shares_and_rates(ca, ya, L)
    wb, pb = _shares_and_rates(cb, yb, L)
    composition, within = _kitagawa(wa, pa, wb, pb)

    exact = (Fraction(sum(k for _, k in a), na)
             - Fraction(sum(k for _, k in b), nb))
    assert abs(float(composition + within) - float(exact)) < 1e-12

    # The same identity, term by term, from the integer counts alone.
    comp_x = within_x = Fraction(0)
    for (n_a, k_a), (n_b, k_b) in zip(a, b):
        w_a, w_b = Fraction(n_a, na), Fraction(n_b, nb)
        p_a = Fraction(k_a, n_a) if n_a else None
        p_b = Fraction(k_b, n_b) if n_b else None
        if p_a is None and p_b is None:
            continue
        p_a = p_b if p_a is None else p_a
        p_b = p_a if p_b is None else p_b
        comp_x += (w_a - w_b) * (p_a + p_b) / 2
        within_x += (w_a + w_b) / 2 * (p_a - p_b)
    assert abs(float(composition) - float(comp_x)) < 1e-12
    assert abs(float(within) - float(within_x)) < 1e-12


# --------------------------------------------------------------------------- #
# Exact binomial interval
# --------------------------------------------------------------------------- #

@given(st.integers(1, 2000).flatmap(
    lambda n: st.tuples(st.just(n), st.integers(0, n))))
def test_clopper_pearson_brackets_the_proportion_and_respects_the_edges(nx):
    n, x = nx
    low, high = _clopper_pearson(x, n)
    assert 0.0 <= low <= x / n <= high <= 1.0
    if x == 0:
        assert low == 0.0
    if x == n:
        assert high == 1.0
    if 0 < x < n:
        assert low < x / n < high


@given(st.integers(2, 500).flatmap(
    lambda n: st.tuples(st.just(n), st.integers(0, n - 1))))
def test_clopper_pearson_moves_with_the_count(nx):
    """One more success can only move both limits up."""
    n, x = nx
    l0, h0 = _clopper_pearson(x, n)
    l1, h1 = _clopper_pearson(x + 1, n)
    assert l1 >= l0 and h1 >= h0


# --------------------------------------------------------------------------- #
# Permutation p-value
# --------------------------------------------------------------------------- #

@given(st.lists(st.integers(-50, 50), min_size=0, max_size=60),
       st.integers(-60, 60))
def test_permutation_pvalue_is_the_phipson_smyth_ratio(null, observed):
    p = permutation_pvalue(null, float(observed))
    B = len(null)
    b = sum(1 for s in null if s >= observed)
    assert p == float(Fraction(b + 1, B + 1))
    assert 0.0 < p <= 1.0


@given(st.lists(st.integers(-50, 50), min_size=1, max_size=60),
       st.integers(-60, 59))
def test_permutation_pvalue_falls_as_the_observed_statistic_rises(null, obs):
    assert (permutation_pvalue(null, float(obs + 1))
            <= permutation_pvalue(null, float(obs)))


# --------------------------------------------------------------------------- #
# Benjamini-Hochberg
# --------------------------------------------------------------------------- #

_pfamily = st.lists(st.integers(0, 1000).map(lambda k: Fraction(k, 1000)),
                    min_size=1, max_size=40)


@given(_pfamily, st.integers(1, 20).map(lambda k: Fraction(k, 40)))
def test_benjamini_hochberg_matches_the_step_up_definition(family, q):
    """adjusted_(i) = min_{j >= i} m p_(j) / j, computed exactly, and the
    rejection set is exactly the set with adjusted p at or below q."""
    raw = np.array([float(p) for p in family])
    adjusted, reject = benjamini_hochberg(raw, q=float(q))
    m = len(family)
    order = sorted(range(m), key=lambda i: family[i])
    exact = [Fraction(0)] * m
    running = Fraction(1)
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, family[i] * m / rank)
        exact[i] = running
    for i in range(m):
        assert abs(adjusted[i] - float(exact[i])) < 1e-12
        assert bool(reject[i]) == (float(exact[i]) <= float(q))


@given(_pfamily, st.integers(1, 19).map(lambda k: Fraction(k, 40)))
def test_benjamini_hochberg_rejections_grow_with_q(family, q):
    raw = np.array([float(p) for p in family])
    _, r_small = benjamini_hochberg(raw, q=float(q))
    _, r_large = benjamini_hochberg(raw, q=float(q + Fraction(1, 40)))
    assert not np.any(r_small & ~r_large)


# --------------------------------------------------------------------------- #
# Merging exchangeable p-values
# --------------------------------------------------------------------------- #

@given(st.lists(st.floats(0.0, 1.0), min_size=1, max_size=25),
       st.sampled_from(["ruger", "twice_mean", "mmb"]))
def test_symmetric_merges_are_probabilities_and_ignore_the_order(ps, how):
    merged = merge_pvalues(ps, how=how)
    assert 0.0 <= merged <= 1.0
    # Reversing the order changes the order of a floating-point sum, so the
    # comparison allows the rounding that a sum in another order carries.
    assert abs(merge_pvalues(list(reversed(ps)), how=how) - merged) < 1e-12


@given(st.lists(st.floats(0.0, 1.0), min_size=1, max_size=25))
def test_exchangeable_ruger_never_exceeds_the_plain_ruger_merge(ps):
    """The exchangeable rule takes the best prefix of the given order and the
    full family is one of those prefixes, so it can only come out at or below
    the plain rule. It is order-dependent on purpose, which is why the package
    reports its order sensitivity separately rather than hiding it."""
    merged = merge_pvalues(ps, how="exchangeable_ruger")
    assert 0.0 <= merged <= merge_pvalues(ps, how="ruger") + 1e-12


@given(st.lists(st.integers(0, 1000).map(lambda k: Fraction(k, 1000)),
                min_size=1, max_size=25))
def test_twice_the_mean_is_twice_the_mean(ps):
    merged = merge_pvalues([float(p) for p in ps], how="twice_mean")
    exact = min(Fraction(1), 2 * sum(ps) / len(ps))
    assert abs(merged - float(exact)) < 1e-12


# --------------------------------------------------------------------------- #
# Joint homoplasy tail: inclusion-exclusion against brute-force enumeration
# --------------------------------------------------------------------------- #

def _tail_by_enumeration(E, c, m):
    """P(M >= m) where trait t places c_t changes on distinct edges chosen
    uniformly among the E, independently across traits, and M counts the
    edges every trait hit. Every placement is enumerated, so the answer is
    an exact rational."""
    edges = range(E)
    hits = Fraction(0)
    total = 1
    for ct in c:
        total *= comb(E, ct)

    def rec(t, common):
        nonlocal hits
        if t == len(c):
            if len(common) >= m:
                hits += 1
            return
        for s in combinations(edges, c[t]):
            rec(t + 1, common & set(s))
    rec(0, set(edges))
    return hits / total


@given(st.integers(1, 5).flatmap(
    lambda E: st.tuples(st.just(E),
                        st.lists(st.integers(0, E), min_size=2, max_size=3),
                        st.integers(1, E))))
@settings(max_examples=60, deadline=None)
def test_exact_upper_tail_agrees_with_enumeration(Ecm):
    E, c, m = Ecm
    expected = _tail_by_enumeration(E, c, m)
    assert abs(_exact_upper_tail(E, c, m) - float(expected)) < 1e-12


@given(st.integers(1, 30).flatmap(
    lambda E: st.tuples(st.just(E), st.lists(st.integers(0, E),
                                             min_size=2, max_size=4))))
def test_exact_upper_tail_is_a_survival_function(Ec):
    E, c = Ec
    assert _exact_upper_tail(E, c, 0) == 1.0
    prev = 1.0
    for m in range(1, E + 2):
        cur = _exact_upper_tail(E, c, m)
        assert 0.0 <= cur <= prev
        prev = cur


# --------------------------------------------------------------------------- #
# Clonal share: a lineage label is a label
# --------------------------------------------------------------------------- #

@given(st.integers(0, 2 ** 31 - 1),
       st.lists(st.integers(0, 5), min_size=24, max_size=60))
@settings(max_examples=8, deadline=None)
def test_clonal_share_does_not_depend_on_what_the_lineages_are_called(seed, lin):
    """Renaming every lineage with a bijection that reverses the sort order
    must leave every number in the result unchanged. The estimator is told it
    receives a label and nothing more; this is the test of that promise."""
    rng = np.random.default_rng(seed)
    lineage = np.asarray(lin)
    y = (rng.random(lineage.size) < 0.3 + 0.1 * lineage).astype(float)
    if y.min() == y.max():
        return
    names_a = np.array(["L%d" % v for v in lineage])
    names_b = np.array(["Z%d" % (9 - v) for v in lineage])
    kw = dict(folds=3, repeats=2, n_boot=15, n_perm=15, seed=seed)
    ra = attribution.clonal_share(y, names_a, **kw).as_dict()
    rb = attribution.clonal_share(y, names_b, **kw).as_dict()
    assert ra.keys() == rb.keys()
    # The bootstrap draws lineages by their code, so under a renaming that
    # reorders the codes the same stream picks different lineages: the
    # interval is the same in distribution but not bit for bit. Everything
    # else is a function of positions and must not move.
    for key in ra:
        if key in ("ci_low", "ci_high"):
            continue
        a, b = ra[key], rb[key]
        if isinstance(a, float) and np.isnan(a):
            assert isinstance(b, float) and np.isnan(b), key
        else:
            assert a == b, key


@given(st.integers(0, 2 ** 31 - 1), st.integers(1, 4))
@settings(max_examples=6, deadline=None)
def test_the_estimate_does_not_depend_on_how_the_null_is_averaged(seed, extra):
    """``null_repeats`` changes what the p-value compares against and nothing
    else. The estimate, its interval and the penalty are read before the extra
    fold draws are taken, from the same stream in the same order, so they
    must be the same numbers to the last bit."""
    rng = np.random.default_rng(seed)
    lineage = rng.integers(0, 8, size=60)
    y = (rng.random(60) < 0.2 + 0.05 * lineage).astype(float)
    if y.min() == y.max():
        return
    kw = dict(folds=3, repeats=4, n_boot=25, n_perm=10, seed=seed)
    one = attribution.clonal_share(y, lineage, null_repeats=1, **kw)
    more = attribution.clonal_share(y, lineage, null_repeats=extra, **kw)
    for name in ("kappa", "kappa_adj", "ci_low", "ci_high", "null_mean",
                 "cv_sd", "support"):
        a, b = getattr(one, name), getattr(more, name)
        assert (a == b) or (np.isnan(a) and np.isnan(b)), name


# --------------------------------------------------------------------------- #
# The branches a first pass of mutation testing found unexercised
# --------------------------------------------------------------------------- #

@given(st.lists(st.integers(-50, 50), min_size=0, max_size=60),
       st.integers(-60, 60))
def test_the_lower_tail_counts_the_other_side(null, observed):
    p = permutation_pvalue(null, float(observed), tail="less")
    b = sum(1 for s in null if s <= observed)
    assert p == float(Fraction(b + 1, len(null) + 1))


def test_an_unknown_tail_is_refused():
    with pytest.raises(ValueError):
        permutation_pvalue([0.0, 1.0], 0.5, tail="both")


@given(_pfamily, st.integers(1, 20).map(lambda k: Fraction(k, 40)))
def test_arbitrary_dependence_is_the_step_up_at_q_over_the_harmonic_sum(family, q):
    """Benjamini-Yekutieli: the same step-up with every adjusted value scaled
    by H_m = sum 1/i, capped at one. Computed exactly beside the float."""
    raw = np.array([float(p) for p in family])
    adjusted, reject = benjamini_hochberg(raw, q=float(q), dependence="arbitrary")
    m = len(family)
    h = sum(Fraction(1, i) for i in range(1, m + 1))
    order = sorted(range(m), key=lambda i: family[i])
    exact = [Fraction(0)] * m
    running = Fraction(1)
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, family[i] * h * m / rank)
        exact[i] = running
    for i in range(m):
        assert abs(adjusted[i] - float(exact[i])) < 1e-12
        assert bool(reject[i]) == (float(exact[i]) <= float(q))
    plain, _ = benjamini_hochberg(raw, q=float(q))
    assert np.all(adjusted >= plain - 1e-12)


def test_a_non_finite_pvalue_stops_the_family_unless_omission_is_asked_for():
    raw = np.array([0.01, np.nan, 0.5])
    with pytest.raises(ValueError):
        benjamini_hochberg(raw)
    adjusted, reject = benjamini_hochberg(raw, nan_policy="omit")
    assert np.isnan(adjusted[1]) and not reject[1]
    kept, _ = benjamini_hochberg(np.array([0.01, 0.5]))
    assert adjusted[0] == kept[0] and adjusted[2] == kept[1]
    with pytest.raises(ValueError):
        benjamini_hochberg(raw, nan_policy="drop")
    with pytest.raises(ValueError):
        benjamini_hochberg(np.array([0.5]), dependence="positive")


_binary = st.integers(2, 12).flatmap(
    lambda n: st.integers(1, 10).flatmap(
        lambda p: st.lists(st.lists(st.integers(0, 1), min_size=p, max_size=p),
                           min_size=n, max_size=n)))


@given(_binary, st.sampled_from(["jaccard", "dice", "simple_matching", "hamming"]))
def test_binary_distances_match_their_definitions_pair_by_pair(rows, metric):
    """a, b, c, d counted by hand for every pair, as exact rationals."""
    X = np.array(rows, dtype=float)
    D = binary_distance_matrix(X, metric=metric, undefined_pair="nan")
    n, p = X.shape
    assert D.shape == (n, n) and np.all(np.diag(D) == 0.0)
    for i in range(n):
        for j in range(i + 1, n):
            a = int(np.sum((X[i] == 1) & (X[j] == 1)))
            b = int(np.sum((X[i] == 1) & (X[j] == 0)))
            c = int(np.sum((X[i] == 0) & (X[j] == 1)))
            if metric == "jaccard":
                want = None if a + b + c == 0 else 1 - Fraction(a, a + b + c)
            elif metric == "dice":
                want = None if 2 * a + b + c == 0 else 1 - Fraction(2 * a, 2 * a + b + c)
            else:
                want = Fraction(b + c, p)
            if want is None:
                assert np.isnan(D[i, j]) and np.isnan(D[j, i])
            else:
                assert abs(D[i, j] - float(want)) < 1e-12
                assert D[i, j] == D[j, i]


def test_an_empty_union_is_scored_by_the_stated_convention():
    X = np.array([[0, 0, 0], [0, 0, 0], [1, 0, 0]], dtype=float)
    assert binary_distance_matrix(X, undefined_pair="identical")[0, 1] == 0.0
    assert binary_distance_matrix(X, undefined_pair="distinct")[0, 1] == 1.0
    assert np.isnan(binary_distance_matrix(X, undefined_pair="nan")[0, 1])
    # A row that is all zero is uninformative; a row with any call is not.
    assert uninformative_rows(X).tolist() == [True, True, False]
    with pytest.raises(ValueError):
        binary_distance_matrix(X, metric="cosine")
    with pytest.raises(ValueError):
        binary_distance_matrix(X, undefined_pair="zero")
    assert np.array_equal(jaccard_distance_matrix(X, undefined_pair="distinct"),
                          binary_distance_matrix(X, metric="jaccard",
                                                 undefined_pair="distinct"))


@given(st.integers(2, 30), st.integers(1, 8))
def test_effective_dimension_lies_between_one_and_the_rank(n, p):
    rng = np.random.default_rng(n * 31 + p)
    X = (rng.random((n, p)) < 0.5).astype(float)
    if np.allclose(X, X[0]):
        return
    d = effective_dimension(X)
    assert 1.0 - 1e-9 <= d <= min(n, p) + 1e-9
    # The columns of a full factorial design are pairwise uncorrelated, so the
    # spectrum is flat and the effective dimension is the number of columns.
    F = np.array([[(i >> k) & 1 for k in range(3)] for i in range(8)], dtype=float)
    assert abs(effective_dimension(F) - 3.0) < 1e-9
