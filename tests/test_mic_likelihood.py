import importlib
import numpy as np
import pytest
from amr_clonalshare.censored import marginal_loglik


def api():
    spec=importlib.util.find_spec('amr_clonalshare._mic_likelihood')
    assert spec is not None, 'Full MIC likelihood fitting is not implemented'
    return importlib.import_module('amr_clonalshare._mic_likelihood')


def data():
    rng=np.random.default_rng(81733); c=np.repeat(np.arange(5),[3,6,2,10,4])
    y=rng.normal(size=5)[c]+rng.normal(size=len(c))*.8
    lo=np.floor(y); hi=lo+1
    lo[y<=-1]=-np.inf; hi[y<=-1]=-1
    lo[y>1]=1; hi[y>1]=np.inf
    return lo,hi,c


@pytest.mark.parametrize('rho',[0.,.1,.5,.9])
def test_histogram_agrees_with_independently_adaptive_likelihood(rho):
    m=api(); lo,hi,c=data(); problem=m.GaussianMICLikelihood(lo,hi,c)
    ll,grad=problem.evaluate(.17,1.3,rho,order=48)
    expected=marginal_loglik(lo,hi,c,5,.17,rho*1.3**2,1.3*np.sqrt(1-rho))
    assert ll == pytest.approx(expected,abs=2e-7)
    assert np.isfinite(grad).all()


def test_histogram_gradient_matches_finite_difference():
    m=api(); lo,hi,c=data(); p=m.GaussianMICLikelihood(lo,hi,c)
    x=np.array([.1,np.log(1.2),.45]); value,gradient=p.objective(x,order=48)
    numeric=[]
    for k in range(3):
        step=np.zeros(3); step[k]=1e-5
        numeric.append((p.objective(x+step,order=48)[0]-p.objective(x-step,order=48)[0])/(2e-5))
    assert gradient == pytest.approx(numeric,rel=3e-5,abs=3e-5)


def test_fixed_rho_profiles_mean_and_scale():
    m=api(); lo,hi,c=data(); p=m.GaussianMICLikelihood(lo,hi,c)
    fit=p.fit(rho=.45)
    assert fit.converged
    initial=p.evaluate(0,1,.45,order=48)[0]
    assert fit.loglik >= initial-1e-6
    assert fit.rho == .45


def test_unrestricted_fit_dominates_fixed_candidates():
    m=api(); lo,hi,c=data(); p=m.GaussianMICLikelihood(lo,hi,c)
    fit=p.fit()
    assert fit.converged
    for r in (0,.25,.5,.85):
        assert fit.loglik >= p.fit(rho=r).loglik-1e-5


def test_singleton_rich_optimizer_excursion_is_resolved():
    from benchmarks.mic_inference.calibrate import generate,design_grid
    a,b,c,_,_=generate(design_grid()[32],0,'pilot')
    fit=api().GaussianMICLikelihood(a,b,c).fit()
    assert fit.converged


@pytest.mark.parametrize('rho',[.01,.5,.95])
def test_singleton_likelihood_is_integrated_analytically(rho):
    m=api(); a,b,c=data()
    c=c.copy(); c[:3]=np.arange(5,8)
    _,c=np.unique(c,return_inverse=True)
    p=m.GaussianMICLikelihood(a,b,c)
    ll,_=p.evaluate(.3,1.2,rho,order=48)
    expected=marginal_loglik(a,b,c,len(np.unique(c)),.3,1.44*rho,1.2*np.sqrt(1-rho))
    assert ll==pytest.approx(expected,abs=3e-7)


def test_high_rho_profile_has_verified_quadrature():
    from benchmarks.mic_inference.calibrate import generate,design_grid
    a,b,c,_,_=generate(design_grid()[50],0,'pilot')
    fit=api().GaussianMICLikelihood(a,b,c).fit(rho=.995)
    assert fit.converged,fit.as_dict()
    assert fit.quadrature_error<1e-5


@pytest.mark.parametrize('cell,rep',[(56,2),(11,1)])
def test_asymmetric_panel_profiles_match_independent_integrator(cell,rep):
    from benchmarks.mic_inference.calibrate import generate,design_grid
    a,b,c,_,_=generate(design_grid()[cell],rep,'pilot')
    p=api().GaussianMICLikelihood(a,b,c)
    fit=p.fit(rho=.995)
    assert fit.converged,fit.as_dict()
    independent=marginal_loglik(a,b,c,len(np.unique(c)),fit.mean,fit.rho*fit.total_sd**2,fit.total_sd*np.sqrt(1-fit.rho))
    assert fit.loglik==pytest.approx(independent,abs=1e-5)


def test_covariate_gradient_matches_finite_difference():
    m=api(); lo,hi,c=data()
    laboratory=np.array(['A','B','A','B','A','B','C','A','B','A','B','C','A','B','A','B','A','B','C','A','B','A','B','A','B'])[:len(c)]
    p=m.GaussianMICLikelihood(lo,hi,c,covariate=laboratory)
    assert p.levels == (('A','B','C'),) and p.n_coefficients == 2
    x=np.array([.1,np.log(1.2),.45,.3,-.4]); value,gradient=p.objective(x,order=48)
    numeric=[]
    for k in range(5):
        step=np.zeros(5); step[k]=1e-5
        numeric.append((p.objective(x+step,order=48)[0]-p.objective(x-step,order=48)[0])/(2e-5))
    assert gradient == pytest.approx(numeric,rel=3e-5,abs=3e-5)
    # a covariate of one level is the model without one
    same=m.GaussianMICLikelihood(lo,hi,c,covariate=np.full(len(c),'A'))
    assert same.n_coefficients == 0
    assert same.evaluate(.1,1.2,.45)[0] == pytest.approx(m.GaussianMICLikelihood(lo,hi,c).evaluate(.1,1.2,.45)[0])


def test_covariate_shift_is_removed_from_the_variance_components():
    m=api(); rng=np.random.default_rng(7)
    c=np.repeat(np.arange(24),12); n=len(c)
    laboratory=np.where(rng.random(n)<.85,c%2,1-c%2)
    y=rng.normal(0,np.sqrt(.1),24)[c]+rng.normal(0,np.sqrt(.9),n)+1.5*laboratory
    lo=np.floor(y); hi=lo+1
    plain=m.GaussianMICLikelihood(lo,hi,c).fit()
    adjusted=m.GaussianMICLikelihood(lo,hi,c,covariate=laboratory).fit()
    assert adjusted.converged and len(adjusted.coefficients) == 1
    assert adjusted.coefficients[0] == pytest.approx(1.5,abs=.3)
    # the shift is aligned with lineage, so the plain fit attributes part of it to lineage
    assert adjusted.rho < plain.rho
    assert abs(adjusted.rho-.1) < abs(plain.rho-.1)
    restricted=m.GaussianMICLikelihood(lo,hi,c,covariate=laboratory).fit(rho=.1)
    assert restricted.converged and restricted.loglik <= adjusted.loglik+1e-6


def test_exact_readings_with_a_covariate_use_the_general_likelihood():
    m=api(); rng=np.random.default_rng(11)
    c=np.repeat(np.arange(8),6); y=rng.normal(size=8)[c]+rng.normal(size=48)+np.where(np.arange(48)%2==0,0.,.8)
    p=m.GaussianMICLikelihood(y,y,c,covariate=np.arange(48)%2)
    assert p.all_exact and not p.closed_form
    fit=p.fit()
    assert fit.converged and 0 <= fit.rho < 1 and len(fit.coefficients) == 1


def test_two_covariates_have_additive_effects_and_matching_gradients():
    m=api(); rng=np.random.default_rng(19)
    c=np.repeat(np.arange(16),8); n=len(c)
    laboratory=np.where(rng.random(n)<.85,c%2,1-c%2); year=rng.integers(0,3,n)
    y=rng.normal(0,.6,16)[c]+rng.normal(0,.8,n)+.7*laboratory+np.array([0.,.4,.9])[year]
    lo=np.floor(y); hi=lo+1
    p=m.GaussianMICLikelihood(lo,hi,c,covariate=[laboratory,year])
    assert p.n_levels == (2,3) and p.n_coefficients == 3
    x=np.array([.1,np.log(1.1),.3,.5,.3,.8]); value,gradient=p.objective(x,order=48)
    numeric=[]
    for k in range(6):
        step=np.zeros(6); step[k]=1e-5
        numeric.append((p.objective(x+step,order=48)[0]-p.objective(x-step,order=48)[0])/(2e-5))
    assert gradient == pytest.approx(numeric,rel=3e-5,abs=3e-5)
    fit=p.fit()
    assert fit.converged and len(fit.coefficients) == 3
    assert fit.coefficients[0] == pytest.approx(.7,abs=.35)
    assert fit.coefficients[2] > fit.coefficients[1]


def test_a_covariate_level_censored_on_one_side_everywhere_is_refused():
    lo,hi,c=data()
    level=np.where(np.arange(lo.size)%3==0,'B','A')
    lo[level=='B']=2.; hi[level=='B']=np.inf
    with pytest.raises(ValueError,match='censored on one side'):
        api().GaussianMICLikelihood(lo,hi,c,covariate=level)
