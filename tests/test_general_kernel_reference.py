import numpy as np
import pytest
from scipy.special import ndtri

from amr_clonalshare._general_probit.profile_bootstrap_research.bootstrap import restricted_fit, BootstrapEngine, inclusive_rank, node_components, infer_profile_bootstrap, reference_counts


def test_restricted_boundary_zero_is_weighted_prevalence():
    r=restricted_fit([1,8],[2,10],0.)
    assert r['status']=='ok'
    assert r['a']==pytest.approx(ndtri(.75))


def test_restricted_boundary_one_counts_groups_not_observations():
    r=restricted_fit([2,0],[2,10],1.)
    assert r['a']==0.
    assert restricted_fit([1,0],[2,10],1.)['status']=='impossible_null'


def test_one_repeat_is_not_refused_and_complement_is_symmetric():
    r=restricted_fit([1],[2],.5)
    assert r['status']=='ok' and abs(r['a'])<1e-5
    a=restricted_fit([0,1,1],[1,1,5],.6)
    b=restricted_fit([1,0,4],[1,1,5],.6)
    assert a['a']==pytest.approx(-b['a'],abs=2e-6)
    assert a['nll']==pytest.approx(b['nll'],abs=1e-9)


def test_rank_inclusive_ties_and_infinities():
    assert inclusive_rank(1.,[0.,1.,2.])==.75
    assert inclusive_rank(float('inf'),[0.,1.,2.])==.25
    with pytest.raises(ValueError):inclusive_rank(1.,[np.nan])


def test_components_keep_disconnected_and_full_cells():
    assert node_components([0.,.5,1.],[True,False,True])==[[0.,.25],[.75,1.]]
    assert node_components([0.,.5,1.],[True,True,True])==[[0.,1.]]
    assert node_components([0.,.5,1.],[False,False,False])==[]


def test_node_is_reproducible_and_permutation_invariant():
    x=BootstrapEngine([2,5],a_grid=[-1.,0.,1.],alternative_rhos=[0.,.5,1.])
    y=BootstrapEngine([5,2],a_grid=[-1.,0.,1.],alternative_rhos=[0.,.5,1.])
    a=x.evaluate_node([1,4],.5,B=19,seed=713,case_key='fixture')
    b=y.evaluate_node([4,1],.5,B=19,seed=713,case_key='fixture')
    assert a['exceedances']==b['exceedances']
    assert a['p_value']==b['p_value']
    assert a['observed_statistic']==b['observed_statistic']
    assert a['reference_counts_sha256']==b['reference_counts_sha256']


def test_rho_one_impossible_null_has_minimal_rank():
    x=BootstrapEngine([2,2],a_grid=[-1.,0.,1.],alternative_rhos=[0.,.5,1.])
    r=x.evaluate_node([1,0],1.,B=19,seed=713,case_key='fixture')
    assert r['status']=='impossible_null' and r['p_value']==.05


def test_grid_score_uses_added_null_and_nonnegative_statistic():
    x=BootstrapEngine([2,2],a_grid=[-1.,0.,1.],alternative_rhos=[0.,1.])
    scores=x.statistics([[1,1],[0,2]],.5)
    assert np.isfinite(scores).all() and (scores>=0).all()


def test_constant_and_singleton_are_explicit_conservative_branches():
    r=infer_profile_bootstrap([0,1],[1,1],B=19,rho_grid=[0.,.5,1.])
    assert r['status']=='structurally_unidentified' and r['components']==[[0.,1.]]
    r=infer_profile_bootstrap([0,0],[2,5],B=19,rho_grid=[0.,.5,1.])
    assert r['status']=='constant_outcome' and r['components']==[[0.,1.]]


def test_invalid_input_and_simulation_budget_are_rejected():
    with pytest.raises(ValueError):restricted_fit([3],[2],.5)
    with pytest.raises(ValueError):restricted_fit([1],[2],-1.)
    x=BootstrapEngine([2],a_grid=[-1.,0.,1.],alternative_rhos=[0.,1.])
    with pytest.raises(ValueError):x.evaluate_node([1],.5,B=0,seed=713,case_key='fixture')


@pytest.mark.parametrize('rho',[0.,.4,1.])
def test_reference_extension_is_a_literal_count_prefix(rho):
    a=reference_counts([1,2,5],-.7,rho,B=19,seed=713,case_key='prefix')
    b=reference_counts([1,2,5],-.7,rho,B=59,seed=713,case_key='prefix')
    assert np.array_equal(a,b[:19])


def test_empty_is_unavailable_and_constant_arguments_still_validated():
    assert infer_profile_bootstrap([],[])['status']=='unavailable'
    with pytest.raises(ValueError):infer_profile_bootstrap([0],[2],B=0)
    with pytest.raises(ValueError):infer_profile_bootstrap([0],[2],alphas=[1.])


def test_numerical_failure_retains_unresolved_cells_but_is_not_usable(monkeypatch):
    engine=BootstrapEngine([2,2],a_grid=[-1.,0.,1.],alternative_rhos=[0.,.5,1.])
    def fail(*args,**kwargs):raise RuntimeError('deliberate numerical fixture')
    monkeypatch.setattr(engine,'evaluate_node',fail)
    r=infer_profile_bootstrap([1,0],[2,2],B=19,rho_grid=[0.,.5,1.],engine=engine,refinement_levels=0)
    assert r['status']=='numerical_incomplete' and not r['usable']
    assert len(r['numerical_errors'])==3
    assert all(s['components']==[[0.,1.]] for s in r['sets'])


def test_constant_reference_statistics_use_analytic_prevalence_limits():
    e=BootstrapEngine([2,5],a_grid=[-1.,0.,1.],alternative_rhos=[0.,.5,1.])
    for r in (0.,.33,1.):assert np.array_equal(e.statistics([[0,0],[2,5]],r),[0.,0.])


def test_supplied_engine_cannot_reinterpret_paired_sizes():
    e=BootstrapEngine([2,5],a_grid=[-1.,0.,1.],alternative_rhos=[0.,.5,1.])
    for sizes in ([5,2],[2,6]):
        with pytest.raises(ValueError):infer_profile_bootstrap([1,1],sizes,engine=e,B=19)


def test_canonical_mesh_nodes_have_identical_adaptive_global_streams():
    from amr_clonalshare._general_probit.profile_bootstrap_research.bootstrap import canonical_node
    adaptive=canonical_node((canonical_node(.1)+canonical_node(.125))/2)
    global_node=canonical_node(36/320)
    assert adaptive.hex()==global_node.hex()
    assert np.array_equal(reference_counts([2,5],-.7,adaptive,B=19,seed=713,case_key='mesh'),
                          reference_counts([2,5],-.7,global_node,B=19,seed=713,case_key='mesh'))


def test_byte_budget_retains_more_than_eight_small_null_tables():
    e=BootstrapEngine([2,5],a_grid=[-1.,0.,1.],alternative_rhos=[0.,.5,1.])
    for i in range(10):e._null(i/10)
    assert len(e.null_tables)==10
    assert sum(t.nbytes for t in e.null_tables.values())+e.alternative.nbytes<=e.budget

