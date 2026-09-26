"""The shapes at the edge of every estimator's domain.

A cohort of one. Two lineages, which is the fewest that admits a between
comparison. A trait that never varies. A bootstrap of one draw and of two. A
support exactly at the threshold rather than near it. Each of these sits on a
guard in the code, and a guard no test stands on is a guard that can be moved
by one without anything noticing.

Every assertion here is the documented behaviour: a refusal by name, a result
that declares itself not estimable with a reason, or a number that is forced
by the shape of the input rather than by the draw.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from amr_clonalshare import clonality, evalues, missingness, qc
from amr_clonalshare.attribution import clonal_share
from amr_clonalshare.realised import realised_share
from amr_clonalshare.stats import (benjamini_hochberg, effective_dimension,
                                   fisher_exact_p, permutation_pvalue)

BUDGET = dict(folds=2, repeats=1, n_boot=20, n_perm=9, seed=0)


# --------------------------------------------------------------------------- #
# The smallest cohorts
# --------------------------------------------------------------------------- #

def test_one_isolate_is_not_a_collection():
    result = clonal_share([1.0], ["A"], **BUDGET)
    assert result.estimable is False
    assert result.n == 1 and result.n_groups == 1


def test_two_lineages_are_the_fewest_that_admit_a_comparison():
    """One lineage cannot be compared with anything; two can."""
    one = clonal_share([0.0, 1.0, 0.0, 1.0], ["A"] * 4, **BUDGET)
    two = clonal_share([0.0, 0.0, 1.0, 1.0], ["A", "A", "B", "B"], **BUDGET)
    assert one.n_groups == 1 and one.estimable is False
    assert two.n_groups == 2 and np.isfinite(two.kappa)


def test_group_adequacy_admits_exactly_two_lineages():
    """The gate is two, not three: a two-lineage cohort is estimable."""
    two = qc.group_adequacy(pd.Series(["A", "A", "B", "B"]))
    one = qc.group_adequacy(pd.Series(["A", "A", "A", "A"]))
    assert two["n_groups"] == 2 and two["estimable"] is True
    assert one["n_groups"] == 1 and one["estimable"] is False


def test_two_repeated_lineages_are_enough_whatever_the_support():
    """Support is reported, not gated: two lineages with two or more
    isolates admit the share however many singletons sit beside them, and
    one repeated lineage does not."""
    labels = ["A"] * 2 + ["B"] * 2 + [f"S{i}" for i in range(16)]
    adequacy = qc.group_adequacy(pd.Series(labels))
    assert adequacy["support"] == pytest.approx(0.2)
    assert adequacy["estimable"] is True
    labels = ["A"] * 5 + [f"S{i}" for i in range(5)]
    below = qc.group_adequacy(pd.Series(labels))
    assert below["support"] == pytest.approx(0.5)
    assert below["estimable"] is False


def test_a_trait_that_never_varies_has_nothing_to_attribute():
    for value in (0.0, 1.0):
        result = clonal_share([value] * 8, ["A"] * 4 + ["B"] * 4, **BUDGET)
        assert result.prevalence == value
        assert not np.isfinite(result.kappa) or result.kappa == 0.0 or not result.estimable


def test_a_constant_outcome_leaves_the_variance_ratio_undefined():
    result = realised_share([1.0] * 12, ["A"] * 6 + ["B"] * 6)
    assert result.estimable is False
    assert result.reason


def test_one_populated_lineage_is_refused_by_name():
    result = realised_share([0.0, 1.0, 0.5, 0.25], ["A"] * 4)
    assert result.estimable is False
    assert "1 populated lineage" in result.reason


# --------------------------------------------------------------------------- #
# The smallest budgets
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("n_boot", [0, 1, 2])
def test_the_smallest_bootstrap_budgets_are_honoured(n_boot):
    y = np.array([0.0, 1.0] * 6)
    lineage = np.repeat(["A", "B", "C"], 4)
    result = clonal_share(y, lineage, folds=2, repeats=1, n_boot=n_boot,
                          n_perm=9, seed=0)
    if n_boot == 0:
        assert not np.isfinite(result.ci_low) and not np.isfinite(result.ci_high)
        assert result.reason == ""
    else:
        assert np.isfinite(result.kappa_adj)


@pytest.mark.parametrize("n_perm", [0, 1, 2])
def test_the_smallest_permutation_budgets_are_honoured(n_perm):
    y = np.array([0.0, 1.0] * 6)
    lineage = np.repeat(["A", "B", "C"], 4)
    result = clonal_share(y, lineage, folds=2, repeats=1, n_boot=20,
                          n_perm=n_perm, seed=0)
    if n_perm == 0:
        assert not np.isfinite(result.p_floor)
    else:
        assert result.p_floor == pytest.approx(1.0 / (n_perm + 1))


def test_a_decomposition_with_one_bootstrap_draw_reports_no_spread():
    """With fewer than two draws there is no standard deviation to report."""
    ya, la = [1, 0, 1, 0], ["A", "A", "B", "B"]
    yb, lb = [0, 0, 1, 1], ["A", "A", "B", "B"]
    one = clonality.decompose_prevalence_difference(ya, la, yb, lb, n_boot=1)
    two = clonality.decompose_prevalence_difference(ya, la, yb, lb, n_boot=2)
    assert one["composition_se"] == 0.0 or np.isnan(one["composition_se"])
    assert np.isfinite(two["composition_se"])


# --------------------------------------------------------------------------- #
# Empty and single-element arrays
# --------------------------------------------------------------------------- #

def test_an_empty_family_of_p_values_is_an_empty_answer():
    adjusted, rejected = benjamini_hochberg([])
    assert list(adjusted) == [] and list(rejected) == []


def test_a_family_of_one_is_its_own_p_value():
    adjusted, rejected = benjamini_hochberg([0.01], 0.05)
    assert list(rejected) == [True] and adjusted[0] == pytest.approx(0.01)


def test_effective_dimension_of_one_column_is_one():
    assert effective_dimension(np.zeros((5, 1)) + np.arange(5)[:, None]) == \
        pytest.approx(1.0)


def test_effective_dimension_of_a_constant_block_is_undefined():
    """No column varies, so there is no correlation matrix to take."""
    assert np.isnan(effective_dimension(np.ones((6, 3))))
    assert effective_dimension(np.zeros((6, 0))) == 0.0


def test_an_empty_permutation_set_gives_the_only_p_value_it_can():
    assert permutation_pvalue([], 1.0) == pytest.approx(1.0)


def test_an_exact_test_on_an_empty_table_is_one():
    assert fisher_exact_p(0, 0, 0, 0) == pytest.approx(1.0)


def test_bounds_of_an_empty_collection_are_empty():
    bounds = missingness.finite_collection_bounds([])
    assert bounds["status"] == "empty"
    assert bounds["lower_bound"] is None and bounds["upper_bound"] is None


def test_bounds_of_a_collection_with_no_observed_outcome_span_everything():
    bounds = missingness.finite_collection_bounds([None, None, None])
    assert bounds["lower_bound"] == 0.0 and bounds["upper_bound"] == 1.0
    assert bounds["observed_count"] == 0


def test_a_difference_of_two_empty_records_is_empty():
    empty = missingness.finite_collection_bounds([])
    difference = missingness.difference_bounds(empty, empty)
    assert difference["status"] == "empty"
    assert difference["lower_bound"] is None


def test_evidence_needs_two_lineages_and_says_so_when_it_has_one():
    result = evalues.e_process([0.0, 1.0, 0.0, 1.0], ["A"] * 4, folds=2,
                               repeats=1, seed=0)
    assert result.n_groups == 1
    assert result.n_splits == 0 or result.e_value == pytest.approx(1.0, abs=0.5)


def test_the_smallest_e_value_products_are_the_identities():
    assert evalues.combine_independent([]) == pytest.approx(1.0)
    assert evalues.combine_independent([2.0]) == pytest.approx(2.0)
    assert evalues.combine_within_cohort([3.0]) == pytest.approx(3.0)


def test_an_e_bh_family_of_one_uses_the_whole_budget():
    decision = evalues.e_bh([40.0], alpha=0.05)
    assert decision["n_rejected"] == 1
    assert decision["threshold"] == pytest.approx(1.0 / 0.05)
