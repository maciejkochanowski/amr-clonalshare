"""The permutation test of the lineage share against references it does not use.

The statistic is the between-lineage sum of squares. For one binary trait it is
the Pearson chi-square statistic of the lineage-by-outcome table times
p(1 - p), and for mid-rank scores of distinct values the Kruskal-Wallis
statistic times (n + 1) / (12 n); both identities are checked against SciPy.
The Monte Carlo p-value is checked against the exact permutation p-value,
enumerated over every distinct labelling, without strata and within strata.
"""
from __future__ import annotations

import itertools

import numpy as np
import pytest
from scipy.stats import chi2_contingency, kruskal

from amr_clonalshare import attribution
from amr_clonalshare.attribution import _between_lineage_ss, _codes


def _centred(x, code):
    x = np.asarray(x, dtype=float).reshape(-1, 1)
    return _between_lineage_ss(x, code) - x.sum() ** 2 / x.size


def test_for_a_binary_trait_the_statistic_is_the_chi_square_of_the_table():
    rng = np.random.default_rng(4)
    code = np.repeat(np.arange(6), [3, 7, 4, 9, 2, 5])
    y = (rng.random(code.size) < 0.2 + 0.1 * code).astype(float)
    table = np.array([[np.sum((code == g) & (y == v)) for v in (0, 1)] for g in range(6)])
    chi2 = chi2_contingency(table, correction=False)[0]
    p = y.mean()
    assert _centred(y, code) == pytest.approx(chi2 * p * (1 - p), rel=1e-12)


def test_for_mid_rank_scores_the_statistic_is_the_kruskal_wallis_statistic():
    rng = np.random.default_rng(5)
    code = np.repeat(np.arange(5), [4, 6, 3, 8, 5])
    values = rng.normal(size=code.size) + 0.5 * code
    n = values.size
    scores = (np.argsort(np.argsort(values)) + 0.5) / n
    h = kruskal(*[values[code == g] for g in range(5)]).statistic
    assert _centred(scores, code) == pytest.approx(h * (n + 1) / (12 * n), rel=1e-12)


def _labellings(code):
    """Every distinct arrangement of the labels, each once."""
    return sorted(set(itertools.permutations(code.tolist())))


def _exact_p(y, code, strata=None):
    """The exact permutation p-value of the between-lineage sum of squares
    of the trait scored within its strata (centred on the stratum mean;
    without strata the trait itself), over the labellings permuted within
    strata."""
    y = np.asarray(y, dtype=float)
    if strata is not None:
        y = y - np.array([y[strata == s].mean() for s in strata])
    x = y.reshape(-1, 1)
    observed = _between_lineage_ss(x, code)
    if strata is None:
        arrangements = [np.array(a) for a in _labellings(code)]
    else:
        parts = []
        for s in np.unique(strata):
            idx = np.flatnonzero(strata == s)
            parts.append((idx, _labellings(code[idx])))
        arrangements = []
        for combo in itertools.product(*[p[1] for p in parts]):
            a = code.copy()
            for (idx, _), labels in zip(parts, combo):
                a[idx] = labels
            arrangements.append(a)
    stats = np.array([_between_lineage_ss(x, a) for a in arrangements])
    return float(np.mean(stats >= observed - 1e-11 * abs(observed))), len(arrangements)


def _monte_carlo_p(y, lineage, strata=None, n_perm=20000):
    r = attribution.layer_clonal_share(y.reshape(-1, 1), lineage, folds=2, repeats=1, n_boot=0,
                                       n_perm=n_perm, seed=11, strata=strata)
    return r.p_value


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_the_monte_carlo_p_value_is_the_exact_permutation_p_value(seed):
    rng = np.random.default_rng(seed)
    code = np.repeat(np.arange(3), [3, 3, 2])
    y = rng.normal(size=code.size) + np.array([0., .8, 1.6])[code]
    exact, count = _exact_p(y, code)
    assert count == 560
    mc = _monte_carlo_p(y, np.array([f"L{c}" for c in code], dtype=object))
    assert abs(mc - exact) <= 4 * np.sqrt(exact * (1 - exact) / 20000) + 2 / 20000


@pytest.mark.parametrize("seed", [3, 4])
def test_within_strata_the_monte_carlo_p_value_is_the_exact_stratified_one(seed):
    rng = np.random.default_rng(seed)
    code = np.array([0, 0, 1, 1, 2, 0, 2, 2, 1])
    strata = np.array(["A"] * 5 + ["B"] * 4, dtype=object)
    # a stratum offset that a permutation across strata would read as lineage
    y = rng.normal(size=code.size) + 3.0 * (strata == "B") + 0.7 * code
    exact, count = _exact_p(y, code, _codes(strata))
    assert count == 30 * 12
    mc = _monte_carlo_p(y, np.array([f"L{c}" for c in code], dtype=object), strata=strata)
    assert abs(mc - exact) <= 4 * np.sqrt(exact * (1 - exact) / 20000) + 2 / 20000


def test_labellings_with_the_same_values_in_the_same_lineages_tie_exactly():
    """Swapping two equal readings between positions of one lineage must give
    the same statistic to the last bit, so that such a permutation counts as
    reaching the observed value."""
    rng = np.random.default_rng(8)
    x = rng.choice([0.1, 0.35, 0.6, 0.85], size=40).reshape(-1, 1)
    code = np.repeat(np.arange(8), 5)
    order = rng.permutation(40)
    within = np.concatenate([np.flatnonzero(code == g)[rng.permutation(5)] for g in range(8)])
    assert _between_lineage_ss(x[within], code) == _between_lineage_ss(x, code)
    assert _between_lineage_ss(x[order], code[order]) == _between_lineage_ss(x, code)


def test_for_a_binary_trait_ties_with_the_observed_value_count_as_reaching_it():
    """A call takes few values, so many labellings tie with the observed sum
    of squares; they must count towards the p-value (here 0.49 with the ties
    and 0.23 without them)."""
    code = np.repeat(np.arange(3), [3, 3, 2])
    y = np.array([1., 0., 0., 1., 0., 0., 1., 1.])
    exact, count = _exact_p(y, code)
    assert count == 560 and exact == pytest.approx(272 / 560)
    x = y.reshape(-1, 1)
    observed = _between_lineage_ss(x, code)
    ties = np.mean([np.isclose(_between_lineage_ss(x, np.array(a)), observed) for a in _labellings(code)])
    assert ties > 0.2
    mc = _monte_carlo_p(y, np.array([f"L{c}" for c in code], dtype=object))
    assert abs(mc - exact) <= 4 * np.sqrt(exact * (1 - exact) / 20000) + 2 / 20000


def test_the_statistic_of_several_traits_is_the_sum_over_the_traits():
    rng = np.random.default_rng(9)
    code = np.repeat(np.arange(4), [5, 3, 6, 2])
    X = rng.normal(size=(code.size, 3))
    whole = _between_lineage_ss(X, code)
    parts = [_between_lineage_ss(X[:, [j]], code) for j in range(3)]
    assert whole == pytest.approx(sum(parts), rel=1e-13)
    sizes = np.bincount(code)
    direct = [np.sum(np.bincount(code, weights=X[:, j]) ** 2 / sizes) for j in range(3)]
    assert parts == pytest.approx(direct, rel=1e-13)


def test_the_sum_within_a_lineage_does_not_depend_on_the_order_of_its_readings():
    """(0.1 + 0.2) + 0.3 and 0.3 + 0.2 + 0.1 differ in the last bit, and summed
    in the order of the rows the statistic takes three values over the
    orders of these readings; summed in sorted order it takes one, for one
    trait and for several."""
    assert (0.1 + 0.2) + 0.3 != (0.3 + 0.2) + 0.1
    code = np.array([0, 0, 0, 1, 1, 1])
    X = np.array([[0.1, 0.3], [0.2, 0.2], [0.3, 0.1], [0.3, 0.1], [0.2, 0.2], [0.1, 0.3]])
    values, single = set(), set()
    for first in itertools.permutations(range(3)):
        for second in itertools.permutations(range(3, 6)):
            rows = list(first) + list(second)
            values.add(_between_lineage_ss(X[rows], code))
            single.add(_between_lineage_ss(X[rows, :1], code))
    assert len(values) == 1 and len(single) == 1
