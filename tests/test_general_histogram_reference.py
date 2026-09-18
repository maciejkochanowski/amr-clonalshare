"""Deterministic numerical references; these tests do not validate an inferential rule."""
import importlib
import numpy as np
import pytest
from amr_clonalshare import _population_probit_numerics as frozen


def api():
    return importlib.import_module('amr_clonalshare._general_probit.histogram_table_research.histogram_table')


def test_boundary_entries_and_uniform_large_group_reference():
    module=api();theta=[(0.,0.),(0.,.5),(0.,1.),(-1.,.25)]
    table=module.build_table([1,2,683],theta)
    for m in (1,2,683):
        start=table.offsets[m]
        for j,(a,rho) in enumerate(theta):
            expected=frozen.selected_count_log_probabilities(m,a,rho,tuple(range(m+1)))
            np.testing.assert_allclose(table.logp[start:start+m+1,j],expected,atol=0,rtol=0)
    np.testing.assert_allclose(table.logp[table.offsets[683]:table.offsets[683]+684,1],-np.log(684),atol=1e-9)
    assert table.max_normalization_error<1e-8


@pytest.mark.parametrize('backend',['dense','csr'])
def test_masked_batch_likelihood_matches_frozen(backend):
    module=api();sizes=[1,2,2,7];theta=[(-1.,0.),(0.,.25),(0.,.5),(0.,1.),(1.,1.)]
    table=module.build_table(sizes,theta)
    counts=np.array([[0,0,0,0],[1,2,2,7],[0,1,2,3],[1,0,2,0]])
    h=module.encode_histograms(counts,sizes,table)
    result=module.score_histograms(h,table,backend=backend,batch_block=2,parameter_block=2)
    expected=np.array([[-frozen.ProfileLikelihood(k,sizes).nll(a,rho) for a,rho in theta] for k in counts])
    assert not np.isnan(result).any()
    assert np.array_equal(np.isneginf(result),np.isneginf(expected))
    np.testing.assert_allclose(result,expected,rtol=1e-12,atol=1e-12)
    assert np.isneginf(result[2,-1]) and np.isfinite(result[0,-1])


def test_histogram_permutation_and_wrong_geometry():
    module=api();sizes=[7,1,2,2];counts=np.array([[3,0,1,2]])
    table=module.build_table(sizes,[(0.,0.),(0.,.5)])
    h=module.encode_histograms(counts,sizes,table)
    reverse=module.encode_histograms(counts[:,::-1],sizes[::-1],table)
    np.testing.assert_array_equal(h.toarray(),reverse.toarray())
    assert h.sum()==4 and table.geometry==tuple(sorted(sizes))
    with pytest.raises(ValueError,match='geometry'):
        module.encode_histograms([[0,0,0,0]],[7,1,2,3],table)


def test_grid_score_handles_impossible_null_without_nan():
    module=api();theta=[(0.,0.),(0.,.5),(0.,1.)]
    ll=np.array([[-4.,-3.,-np.inf],[-4.,-3.,-2.]])
    result=module.grid_lr_scores(ll,theta,1.)
    assert np.isposinf(result[0]) and result[1]==0
    with pytest.raises(ValueError):module.grid_lr_scores(ll,theta,.25)
    with pytest.raises(ValueError):module.grid_lr_scores(np.full((1,3),-np.inf),theta,0.)


@pytest.mark.parametrize('sizes,theta',[([],[(0.,.5)]),([0],[(0.,.5)]),([2],[(0.,1.1)]),([2],[(np.nan,.5)]),([2],[])])
def test_invalid_table_input(sizes,theta):
    with pytest.raises(ValueError):api().build_table(sizes,theta)


def test_memory_budget_and_count_validation():
    module=api()
    with pytest.raises(MemoryError):module.build_table([683],[(0.,.5)],memory_budget_bytes=1)
    table=module.build_table([2],[(0.,.5)])
    for counts in ([[3]],[[.5]],[[np.nan]],[[0,1]]):
        with pytest.raises(ValueError):module.encode_histograms(counts,[2],table)
