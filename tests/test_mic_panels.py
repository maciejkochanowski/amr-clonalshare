import numpy as np
import pytest

from amr_clonalshare import _mic_panels as m
from amr_clonalshare._mic_likelihood import GaussianMICLikelihood


def data():
    c=np.repeat(np.arange(6),8)
    rng=np.random.default_rng(2092)
    y=rng.normal(size=6)[c]+rng.normal(size=len(c))
    edges=np.arange(-3.,4.)
    k=np.searchsorted(edges,y,side='left')
    return np.r_[-np.inf,edges][k],np.r_[edges,np.inf][k],c,edges


def test_panel_boundary_convention():
    a,b=m.observe_panel(np.array([-2.,-1.,0.,1.,2.]),np.array([-1.,0.,1.]))
    np.testing.assert_array_equal(a,[-np.inf,-np.inf,-1.,0.,1.])
    np.testing.assert_array_equal(b,[-1.,-1.,0.,1.,np.inf])


def test_panel_is_explicit_and_validated():
    a,b,c,e=data()
    with pytest.raises(ValueError,match='panel'):
        m._validated_problem(a,b,c,None)
    with pytest.raises(ValueError,match='panel'):
        m._validated_problem(a,b,c,e+.1)


def test_adaptive_modes_survive_optimizer_boundary_excursions():
    a,b,c,e=data(); p=GaussianMICLikelihood(a,b,c); fit=p.fit()
    rng=np.random.default_rng(73)
    for _ in range(8):
        effects=rng.normal(0,fit.total_sd*np.sqrt(fit.rho),p.groups)
        latent=fit.mean+effects[p.code]+rng.normal(0,fit.total_sd*np.sqrt(1-fit.rho),p.n)
    aa,bb=m.observe_panel(latent,e)
    result=GaussianMICLikelihood(aa,bb,c).fit(start=fit,multistart=False)
    assert result.converged


def test_simulated_datasets_carry_the_covariate_levels_by_name():
    lo,hi,c,edges=data()
    level=np.array(['L%d'%(i%12) for i in range(lo.size)],dtype=object)
    problem=m._validated_problem(lo,hi,c,edges,covariate=[level])[0]
    names=m._record_covariate(problem)
    assert len(names)==1 and names[0].tolist()==level.tolist()


def test_one_bare_edge_array_is_read_as_one_panel():
    lo,hi,c,edges=data()
    a,b=m.observe_panels(np.array([-2.,0.5,9.]),edges)
    np.testing.assert_array_equal(a,[-3.,0.,3.])
    np.testing.assert_array_equal(b,[-2.,1.,np.inf])
