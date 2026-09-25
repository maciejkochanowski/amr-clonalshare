import importlib
import numpy as np
import pytest
from scipy.stats import f


def api():
    spec = importlib.util.find_spec('amr_clonalshare.mic_inference')
    assert spec is not None, 'The calibrated MIC inference module is not implemented'
    return importlib.import_module('amr_clonalshare.mic_inference')


def fixture_data():
    rng = np.random.default_rng(97241)
    c = np.repeat(np.arange(6), 9)
    y = rng.normal(0, 1.4, 6)[c] + rng.normal(size=c.size)
    return y, c


def test_balanced_pivot_matches_classical_anova():
    m = api()
    y,c = fixture_data()
    r = m.exact_gaussian_interval(y,c)
    means = y.reshape(6,9).mean(axis=1)
    msb = 9*np.sum((means-means.mean())**2)/5
    msw = np.sum((y-means[c])**2)/48
    q = (msb/msw)/f.ppf([.975,.025],5,48)
    expected = np.clip((q-1)/(q+8),0,1)
    assert [r.low,r.high] == pytest.approx(expected, abs=1e-9)


def test_exact_location_scale_label_invariance():
    m=api(); y,c=fixture_data(); order=np.random.default_rng(44).permutation(len(y))
    a=m.exact_gaussian_interval(y,c)
    b=m.exact_gaussian_interval((7+2.3*y)[order],np.array(['L'+str(v) for v in c])[order])
    assert [a.low,a.high] == pytest.approx([b.low,b.high],abs=1e-9)


def test_unbalanced_pivot_matches_dense_weighted_projection():
    m=api(); rng=np.random.default_rng(15); sizes=np.array([1,2,7,3,9])
    c=np.repeat(np.arange(5),sizes); y=rng.normal(size=len(c))+rng.normal(size=5)[c]
    means=np.bincount(c,weights=y)/sizes; rho=.43; lam=rho/(1-rho)
    v=np.diag(lam+1/sizes); inv=np.linalg.inv(v); one=np.ones(5)
    mu=(one@inv@means)/(one@inv@one)
    q=(means-mu)@inv@(means-mu)
    sse=np.sum((y-means[c])**2)
    assert m.exact_f_statistic(y,c,rho) == pytest.approx(q/4/(sse/(len(y)-5)),rel=1e-11)


@pytest.mark.parametrize('y,c',[(np.arange(8.),np.arange(8)),(np.ones(12),np.repeat([0,1],6))])
def test_exact_degenerate_refused(y,c):
    with pytest.raises(ValueError): api().exact_gaussian_interval(y,c)


@pytest.mark.parametrize('alpha',[0,1,-.1,float('nan')])
def test_invalid_alpha(alpha):
    y,c=fixture_data()
    with pytest.raises(ValueError): api().exact_gaussian_interval(y,c,alpha=alpha)
