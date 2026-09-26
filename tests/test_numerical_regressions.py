# Numerical and public-API regression tests.
import numpy as np
import pytest
from scipy.special import log_ndtr
from scipy.stats import multivariate_normal, truncnorm
from amr_clonalshare import censored as c
from amr_clonalshare.config import ConfigError, from_dict

@pytest.mark.parametrize("a,b", [(8.,9.), (10.,11.), (40.,41.), (8.,8.001)])
def test_upper_tail_probability_matches_survival_reference(a,b):
    expected = log_ndtr(-a) + np.log(-np.expm1(log_ndtr(-b)-log_ndtr(-a)))
    actual = c.marginal_loglik(np.array([a]),np.array([b]),np.array([0]),1,0.,0.,1.)
    assert actual == pytest.approx(expected, abs=1e-9)

@pytest.mark.parametrize("a,b", [(40.,np.inf), (40.,41.), (-41.,-40.)])
def test_truncated_moments_respect_support_in_extreme_tails(a,b):
    mean,var = c._trunc_moments(np.array([a]),np.array([b]),0.,1.)
    expected_mean, expected_var = truncnorm.stats(a,b,moments="mv")
    assert a <= mean[0] <= b
    assert mean[0] == pytest.approx(expected_mean, rel=1e-9)
    assert var[0] == pytest.approx(expected_var, abs=2e-7)

@pytest.mark.parametrize("sizes,tau2,sigma", [([30,30,30],4.,.2),([1,5,30],2.,.3),([3,7,11],.01,2.)])
def test_exact_random_intercept_likelihood_matches_multivariate_normal(sizes,tau2,sigma):
    labels=np.repeat(np.arange(len(sizes)),sizes)
    y=np.concatenate([np.linspace(.3,.7,n) + .13*g for g,n in enumerate(sizes)])
    expected=sum(multivariate_normal.logpdf(y[labels==g],mean=np.zeros(n),cov=sigma**2*np.eye(n)+tau2*np.ones((n,n))) for g,n in enumerate(sizes))
    actual=c.marginal_loglik(y,y,labels,len(sizes),0.,tau2,sigma)
    assert actual == pytest.approx(expected,abs=1e-8)

def test_explicit_wells_never_snap_api_values():
    with pytest.raises(ValueError,match="well|panel|concentration"):
        c.intervals_from_mic([1.3],wells=[.5,1.,2.],operators=["="])

@pytest.mark.parametrize("bad", [-1.,.25,2.,np.inf])
def test_binary_interval_api_rejects_nonbinary_values(bad):
    with pytest.raises(ValueError,match="binary|zero|one|0/1"):
        c.intervals_from_binary([0.,1.,bad])

def test_binary_interval_api_preserves_missing():
    lo,hi=c.intervals_from_binary([0.,1.,np.nan])
    assert np.isnan(lo[2]) and np.isnan(hi[2])

def test_from_dict_rejects_invalid_budget_without_separate_validation():
    with pytest.raises(ConfigError,match="n_boot"):
        from_dict({"dataset":{"name":"audit","metadata":"metadata.csv","lineage_column":"lineage","phenotype":"calls.csv"},"attribution":{"n_boot":-1}})

@pytest.mark.parametrize("center,width,sigma,tau2", [(.5,.1,.2,4.), (1.5,.2,.15,9.), (-.7,.8,.3,2.)])
def test_interval_random_intercept_matches_adaptive_integration(center,width,sigma,tau2):
    from scipy.integrate import quad
    from scipy.optimize import minimize_scalar
    from scipy.stats import norm
    n=30
    a=center-width/2
    b=center+width/2
    def log_integrand(u):
        aa=(a-u)/sigma
        bb=(b-u)/sigma
        if aa>0:
            log_right,log_left=log_ndtr(-aa),log_ndtr(-bb)
        else:
            log_right,log_left=log_ndtr(bb),log_ndtr(aa)
        mass=log_right+np.log(-np.expm1(min(log_left-log_right,0.)))
        return n*mass+norm.logpdf(u,scale=np.sqrt(tau2))
    mode=minimize_scalar(lambda u:-log_integrand(u),bounds=(center-5,center+5),method="bounded").x
    peak=log_integrand(mode)
    integral=quad(lambda u:np.exp(log_integrand(u)-peak),mode-20,mode,epsabs=1e-11)[0]
    integral+=quad(lambda u:np.exp(log_integrand(u)-peak),mode,mode+20,epsabs=1e-11)[0]
    expected=np.log(integral)+peak
    actual=c.marginal_loglik(np.full(n,a),np.full(n,b),np.zeros(n,dtype=int),1,0.,tau2,sigma)
    assert actual == pytest.approx(expected,abs=1e-7)


def test_narrow_conditional_density_integrates_despite_rounded_newton_steps():
    from pathlib import Path
    import pandas as pd
    from scipy.special import logsumexp
    from amr_clonalshare.attribution import _codes
    data = Path(__file__).resolve().parents[1]/'examples'/'ssuis'/'data'
    mic = pd.read_csv(data/'mic_panel.csv', dtype={'genome_id': str})
    meta = pd.read_csv(data/'metadata.csv', dtype={'genome_id': str}).set_index('genome_id')
    block = mic[mic.antibiotic.eq('cefquinome')].set_index('genome_id')
    ids = [i for i in block.index if i in meta.index]
    lo, hi = c.intervals_from_mic(block.loc[ids, 'measurement'].to_numpy(float),
                                  treat_end_wells_as_censored=True)
    code = _codes(meta.loc[ids, 'baps_cluster'].to_numpy(object))
    groups = int(code.max())+1
    grand, tau2, sigma = -5.677722891698857, 0.15114066487503733, 0.012300079545753732
    actual = c.marginal_loglik(lo, hi, code, groups, grand, tau2, sigma)
    u = np.linspace(-6., 6., 60001)
    expected = 0.
    for g in range(groups):
        m = code == g
        ll = sum(c._normal_interval_logmass((a-grand-u)/sigma, (b-grand-u)/sigma)
                 for a, b in zip(lo[m], hi[m]))
        expected += logsumexp(ll - u*u/(2*tau2)) + np.log(u[1]-u[0]) - 0.5*np.log(2*np.pi*tau2)
    assert actual == pytest.approx(expected, abs=1e-6)
