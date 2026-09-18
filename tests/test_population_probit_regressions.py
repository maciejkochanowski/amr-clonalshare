"""Regression cases identified by failed independent V2 validation, now development."""
import numpy as np
import pytest
from amr_clonalshare._population_probit_numerics import ProfileLikelihood,fit_profile

CUT=6.9268385202558775

def test_high_endpoint_excludes_impossible_perfect_correlation():
    k=[3,20,0,20,20,20,0,0,0,20];m=[20]*10
    fit=fit_profile(k,m,critical_value=CUT)
    assert .999<fit['ci_high']<1
    ll=ProfileLikelihood(k,m)
    residual=2*(ll.profile(fit['ci_high'])[0]-fit['nll'])-CUT
    assert abs(residual)<1e-5

def test_mle_finds_interior_dip_between_old_coarse_grid_points():
    ll=ProfileLikelihood([1,1,0,0,2,0,0,1,0,1],[1,1,1,1,20]*2)
    rho,a,best=ll.mle()
    assert rho==pytest.approx(.22661644549703952,abs=1e-5)
    assert best==pytest.approx(9.834580801526776,abs=1e-7)

def test_large_groups_have_finite_nuisance_likelihood_and_interval():
    k=[0,0,4,6,511,1,1,2,10,449,1,0,2,9,420,1,1,5,4,399]
    m=[1,2,5,12,1000]*4
    ll=ProfileLikelihood(k,m)
    assert np.isfinite(ll.nll(8.,.1))
    fit=fit_profile(k,m,critical_value=CUT)
    assert fit['identified'] and fit['ci_low']<=fit['rho_hat']<=fit['ci_high']
