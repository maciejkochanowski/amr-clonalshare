"""Numerical support and explicit no-resampling contracts."""
import numpy as np
import pytest
from amr_clonalshare import censored as c
from amr_clonalshare.config import from_dict


def test_zero_surveillance_bootstrap_budget_retains_point_only_route():
    cfg = from_dict({'dataset': {'name': 'point-only', 'metadata': 'm.csv',
        'lineage_column': 'lineage', 'phenotype': 'p.csv'},
        'surveillance': {'n_boot': 0}})
    assert cfg.surveillance.n_boot == 0


@pytest.mark.parametrize('lower,width', [(12., 1e-8), (40., 1e-8), (-40., 1e-8), (0., 1e-10)])
def test_tiny_truncation_interval_has_supported_mean_and_variance(lower, width):
    upper = lower + width
    mean, var = c._trunc_moments(np.array([lower]), np.array([upper]), 0., 1.)
    assert lower <= mean[0] <= upper
    assert 0 <= var[0] <= (upper-lower)**2/4
    assert var[0] == pytest.approx((upper-lower)**2/12, rel=1e-5, abs=1e-35)
    # The local-midpoint density error is O((1 + mid**2)*width**2).
    middle = (lower+upper)/2
    expected = np.log(upper-lower)-0.5*middle**2-0.5*np.log(2*np.pi)
    actual = c.marginal_loglik([lower], [upper], np.array([0]), 1, 0., 0., 1.)
    assert actual == pytest.approx(expected, abs=2e-11)


def test_extreme_one_sided_moments_match_mills_ratio_expansion():
    a = 1000.
    mean, var = c._trunc_moments(np.array([a]), np.array([np.inf]), 0., 1.)
    expected_mean = a + 1/a - 2/a**3 + 10/a**5
    expected_var = 1/a**2 - 6/a**4 + 50/a**6
    assert mean[0] == pytest.approx(expected_mean, abs=1e-10)
    assert var[0] == pytest.approx(expected_var, rel=1e-7)


@pytest.mark.parametrize('low,high', [(np.inf,np.inf),(-np.inf,-np.inf)])
def test_censored_fit_rejects_infinite_exact_observations(low, high):
    with pytest.raises(ValueError, match='endpoint|infinite|interval'):
        c.censored_clonal_share(np.full(12,low), np.full(12,high), np.repeat(['a','b'],6))
