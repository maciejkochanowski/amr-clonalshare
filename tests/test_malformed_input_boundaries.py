"""Public-boundary checks: malformed data must not become evidence."""
import numpy as np
import pandas as pd
import pytest
from amr_clonalshare.stats import benjamini_hochberg, permutation_pvalue
from amr_clonalshare.attribution import clonal_share
from amr_clonalshare.evalues import e_process
from amr_clonalshare.realised import realised_share
from amr_clonalshare.censored import censored_clonal_share, marginal_loglik


@pytest.mark.parametrize('p', [[-.01, .5], [1.01], [[.1, .2]]])
def test_step_up_rejects_invalid_probability_vectors(p):
    with pytest.raises(ValueError, match='p-value|one-dimensional'):
        benjamini_hochberg(p)


@pytest.mark.parametrize('q', [0., 1., -1., np.nan, np.inf])
def test_step_up_rejects_invalid_error_rates_even_for_empty_family(q):
    with pytest.raises(ValueError, match='q'):
        benjamini_hochberg([], q=q)


@pytest.mark.parametrize('null,observed', [([0., np.nan], 1.), ([0., 1.], np.nan), ([np.inf], 1.)])
def test_permutation_nonfinite_values_do_not_become_minimum_p(null, observed):
    with pytest.raises(ValueError, match='finite'):
        permutation_pvalue(null, observed)


@pytest.mark.parametrize('setting,value', [('folds', 2.5), ('folds', True), ('repeats', 0), ('repeats', -1), ('n_boot', -1), ('n_perm', -1), ('null_repeats', 0)])
def test_attribution_rejects_invalid_execution_budgets(setting, value):
    kw = dict(folds=2, repeats=2, n_boot=0, n_perm=3, seed=1)
    kw[setting] = value
    with pytest.raises(ValueError, match=setting):
        clonal_share(np.tile([0., 1.], 12), np.repeat(np.arange(4), 6), **kw)


@pytest.mark.parametrize('setting,value', [('folds', 1), ('folds', 2.5), ('repeats', 0)])
def test_split_evidence_rejects_invalid_budgets(setting, value):
    kw = dict(folds=2, repeats=2, seed=1)
    kw[setting] = value
    with pytest.raises(ValueError, match=setting):
        e_process(np.tile([0., 1.], 12), np.repeat(np.arange(4), 6), **kw)


@pytest.mark.parametrize('missing', [pd.NA, pd.NaT, ' NULL ', ' N/A '])
def test_realised_route_uses_the_same_missing_label_definition(missing):
    labels = np.array(['A'] * 10 + ['B'] * 10 + [missing] * 4, dtype=object)
    values = np.r_[np.linspace(-1., 0., 10), np.linspace(1., 2., 10), [7.] * 4]
    result = realised_share(values, labels)
    assert result.n == 20 and result.n_groups == 2
    assert result.n_dropped_untyped == 4


def test_mic_singletons_do_not_identify_two_variance_components():
    y = np.linspace(-1., 1., 12)
    result = censored_clonal_share(y, y, np.arange(12))
    assert not result.estimable
    assert 'within-lineage' in result.reason


def test_uninformative_intervals_are_not_counted_as_readings():
    y = np.r_[np.linspace(-1., 0., 10), np.linspace(1., 2., 10)]
    labels = np.repeat(['A', 'B'], 10)
    lo, hi = y.copy(), y.copy()
    lo[:3], hi[:3] = -np.inf, np.inf
    result = censored_clonal_share(lo, hi, labels)
    assert result.n == 17
    assert result.n_dropped_uninformative == 3


@pytest.mark.parametrize('tau2,sigma', [(-1., 1.), (1., 0.), (1., -1.), (np.nan, 1.)])
def test_public_marginal_likelihood_rejects_invalid_variances(tau2, sigma):
    with pytest.raises(ValueError, match='sigma|tau2|finite'):
        marginal_loglik([0., 1.], [0., 1.], [0, 0], 1, 0., tau2, sigma)
