"""Scientific oracles for grouped-binomial probit likelihood and inversion."""
import math
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import ndtr
from scipy.stats import binom
from amr_clonalshare._population_probit_numerics import count_probabilities, fit_profile, selected_count_probabilities


@pytest.mark.parametrize('a,rho', [(-1.372,.1),(-1.372,.6),(-1.372,.95),
                                  (-.674,.95),(0.,.9999),(-2.5,.9999)])
def test_integrated_count_probabilities_against_adaptive_oracle(a,rho):
    m=20
    probs=count_probabilities(m,a,rho)
    def oracle(k):
        def pmf(z):
            eta=(a+math.sqrt(rho)*z)/math.sqrt(1-rho)
            return math.comb(m,k)*ndtr(eta)**k*ndtr(-eta)**(m-k)
        scale=math.sqrt((1-rho)/rho)
        points=-a/math.sqrt(rho)+scale*np.array([-10,-5,-3,-1,0,1,3,5,10])
        return quad(lambda z: pmf(z)
                    *math.exp(-z*z/2)/math.sqrt(2*math.pi),-12,12,
                    points=points[(points>-12)&(points<12)],epsabs=2e-12,epsrel=2e-10,limit=300)[0]
    expected=np.array([oracle(k) for k in range(m+1)])
    assert np.max(np.abs(probs-expected))<2e-9
    assert abs(probs.sum()-1)<1e-10
    assert abs(np.arange(m+1)@probs/m-ndtr(a))<1e-10


def test_null_is_exact_binomial_and_success_failure_symmetry():
    for a in (-2.,0.,1.2):
        assert np.allclose(count_probabilities(20,a,0),binom.pmf(np.arange(21),20,ndtr(a)),rtol=1e-11,atol=1e-14)
        for rho in (.1,.95,.9999,1.):
            assert np.allclose(count_probabilities(20,a,rho),count_probabilities(20,-a,rho)[::-1],atol=1e-12)


def test_constant_data_are_unidentified_and_not_a_precise_zero():
    for k in (np.zeros(10,dtype=int),np.full(10,20)):
        fit=fit_profile(k,np.full(10,20))
        assert fit['rho_hat'] is None and fit['identified'] is False
        assert fit['ci_low']==0 and fit['ci_high']==1
        assert fit['status']=='constant_outcome_unidentified'


def test_boundary_high_icc_and_complementary_counts():
    k=np.array([0,0,0,0,0,0,0,0,20,20])
    fit=fit_profile(k,np.full(10,20))
    other=fit_profile(20-k,np.full(10,20))
    assert fit['rho_hat']==1. and fit['ci_high']==1.
    assert 0<fit['ci_low']<1
    assert abs(fit['ci_low']-other['ci_low'])<1e-6
    assert abs(fit['prevalence_hat']-.2)<1e-8


def test_invalid_counts_refused():
    for k,m in [([21],[20]),([-1],[20]),([.5],[20]),([0],[0]),([0,1],[20])]:
        with pytest.raises(ValueError):
            fit_profile(k,m)


@pytest.mark.parametrize('m,a,rho',[(1,-1.372,.95),(70,-1.372,.2499),
                                  (70,-1.372,.2501),(70,-1.372,.95),
                                  (70,-.674,.9999)])
def test_unequal_size_quadrature_and_hybrid_switch(m,a,rho):
    probs=count_probabilities(m,a,rho)
    scale=math.sqrt((1-rho)/rho)
    points=-a/math.sqrt(rho)+scale*np.array([-10,-5,-3,-1,0,1,3,5,10])
    for k in sorted(set([0,1,m//4,m//2,m-1,m])):
        def integrand(z):
            eta=(a+math.sqrt(rho)*z)/math.sqrt(1-rho)
            return math.comb(m,k)*ndtr(eta)**k*ndtr(-eta)**(m-k)*math.exp(-z*z/2)/math.sqrt(2*math.pi)
        expected=quad(integrand,-12,12,points=points[(points>-12)&(points<12)],
                      epsabs=2e-12,epsrel=2e-10,limit=300)[0]
        assert abs(probs[k]-expected)<2e-9
    assert abs(probs.sum()-1)<1e-10
    assert abs(np.arange(m+1)@probs/m-ndtr(a))<1e-10


def test_singletons_have_no_icc_information():
    fit=fit_profile([0,0,0,1,1],[1,1,1,1,1])
    assert fit['ci_low']==0 and fit['ci_high']==1 and fit['rho_hat'] is None
    assert fit['status']=='insufficient_repeated_groups'


@pytest.mark.parametrize('m',[20,100,500,1000,5000])
def test_exact_uniform_count_law_for_half_correlation(m):
    # a=0,rho=.5 gives q=Phi(Z) uniform, hence Beta-binomial(1,1).
    k=tuple(sorted(set([0,1,m//20,m//4,m//2,m-1,m])))
    p=selected_count_probabilities(m,0.,.5,k)
    assert np.max(np.abs(p*(m+1)-1.))<2e-8


def test_invalid_truth_and_unrepresentable_counts_rejected_before_early_returns():
    for truth in (-.5,2.,np.nan,np.inf):
        with pytest.raises(ValueError):
            fit_profile([0,1],[1,1],truth_rho=truth)
        with pytest.raises(ValueError):
            fit_profile([0,0],[3,3],truth_rho=truth)
    with pytest.raises(ValueError):
        fit_profile([0],[1e20])
