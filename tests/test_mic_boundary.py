"""A monomorphic interval in every group has a valid rho=1 likelihood limit."""
import importlib
import numpy as np
import pytest
from amr_clonalshare._mic_likelihood import GaussianMICLikelihood


def api():
    spec = importlib.util.find_spec('amr_clonalshare._mic_boundary')
    assert spec is not None, 'The boundary likelihood is not implemented'
    return importlib.import_module('amr_clonalshare._mic_boundary')


def test_four_equiprobable_categories_have_analytic_boundary_maximum():
    code = np.repeat(np.arange(4), 10)
    lo = np.repeat([-np.inf, -1., 0., 1.], 10)
    hi = np.repeat([-1., 0., 1., np.inf], 10)
    problem = GaussianMICLikelihood(lo, hi, code)
    fit = api().monomorphic_boundary_fit(problem)
    assert fit is not None and fit.converged
    assert fit.rho == 1.
    assert fit.loglik == pytest.approx(4*np.log(.25), abs=1e-7)
    assert fit.mean == pytest.approx(0., abs=1e-5)


def test_tail_only_boundary_has_binomial_supremum():
    code = np.repeat(np.arange(8), 5)
    lo = np.repeat([-np.inf]*3+[1.]*5, 5)
    hi = np.repeat([-1.]*3+[np.inf]*5, 5)
    fit = api().monomorphic_boundary_fit(GaussianMICLikelihood(lo, hi, code))
    assert fit.loglik == pytest.approx(3*np.log(3/8)+5*np.log(5/8))


def test_mixed_bins_within_a_group_do_not_use_boundary_shortcut():
    code = np.array([0, 0, 1, 1])
    lo = np.array([-1., 0., -1., 0.])
    hi = lo+1.
    assert api().monomorphic_boundary_fit(GaussianMICLikelihood(lo, hi, code)) is None


def test_unrestricted_fit_uses_certified_boundary_and_restricted_fit_still_works():
    code = np.repeat(np.arange(4), 10)
    lo = np.repeat([-np.inf, -1., 0., 1.], 10)
    hi = np.repeat([-1., 0., 1., np.inf], 10)
    problem = GaussianMICLikelihood(lo, hi, code)
    full = problem.fit()
    assert full.converged and full.rho == 1.
    fixed = problem.fit(rho=.8, start=full)
    assert fixed.converged and fixed.loglik <= full.loglik + 1e-8


def test_infinite_scale_supremum_can_warm_start_finite_null():
    code = np.repeat(np.arange(8), 5)
    lo = np.repeat([-np.inf]*3+[1.]*5, 5)
    hi = np.repeat([-1.]*3+[np.inf]*5, 5)
    problem = GaussianMICLikelihood(lo, hi, code)
    full = problem.fit()
    assert full.converged and full.rho == 1.
    assert np.isinf(full.total_sd)
    fixed = problem.fit(rho=.7, start=full)
    assert np.isfinite(fixed.loglik)
    assert fixed.loglik <= full.loglik+1e-7


def test_boundary_likelihood_ignores_within_lineage_repetition():
    a = [-np.inf, -1., 0., 1.]
    b = [-1., 0., 1., np.inf]
    two = GaussianMICLikelihood(np.repeat(a,2), np.repeat(b,2), np.repeat(np.arange(4),2))
    many = GaussianMICLikelihood(np.repeat(a,50), np.repeat(b,50), np.repeat(np.arange(4),50))
    assert api().monomorphic_boundary_fit(two).loglik == pytest.approx(api().monomorphic_boundary_fit(many).loglik)
