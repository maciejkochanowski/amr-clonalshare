import math
import numpy as np
import pytest
from amr_clonalshare._general_probit_upper import mc_upper_bound, invert_grid_pvalues, simulate_scores, all_pairs_reference


def test_outward_grid_hull_retains_late_acceptances():
    grid=[0.,.25,.5,.75,1.]
    assert invert_grid_pvalues(grid,[.8,.01,.8,.01,.01],.05)==.75
    assert invert_grid_pvalues(grid,[.01]*5,.05)==0.
    assert invert_grid_pvalues(grid,[.01,.01,.01,.8,.01],.05)==1.


def test_exact_equality_rejects_inclusive_rank_tail():
    assert invert_grid_pvalues([0.,.5,1.],[.05,.05,.05],.05)==0.


def test_zero_centrality_and_singletons_no_restriction():
    r=mc_upper_bound([0,4,1],[3,4,1],seed=7,B=19,grid=[0.,.5,1.])
    assert r['upper']==1 and r['estimate'] is None and r['lower']==0
    assert r['simulated_replicates']==0


def test_seeded_reproducible_permutation_and_no_global_rng():
    state=np.random.get_state()
    a=mc_upper_bound([1,2,0],[2,5,3],seed=789,B=99,grid=[0.,.5,1.])
    b=mc_upper_bound([0,1,2],[3,2,5],seed=789,B=99,grid=[0.,.5,1.])
    for key in ['upper','pvalues','exceedances','grid','seed','B','score']:
        assert a[key]==b[key]
    after=np.random.get_state()
    assert all(np.array_equal(x,y) for x,y in zip(state,after))
    assert all(p>0 for p in a['pvalues'])
    assert a['pvalues'][-1]==1/100


def test_rank_values_from_seeded_reference():
    sizes=[2,3]; B=31; grid=[0.,.5,1.]; seed=99
    streams=np.random.SeedSequence(seed).spawn(len(grid))
    expected=[]
    for rho,stream in zip(grid,streams):
        samples=simulate_scores(sizes,rho,B,np.random.default_rng(stream),score='centrality')
        expected.append((1+int(np.count_nonzero(samples>=2)))/(B+1))
    assert mc_upper_bound([1,1],sizes,seed=seed,B=B,grid=grid)['pvalues']==expected


def test_simulation_boundaries_and_centrality():
    rng=np.random.default_rng(100)
    assert np.all(simulate_scores([2,5,8],1.,30,rng)==0)
    values=simulate_scores([2,3,8],.5,100,rng)
    assert np.all((values>=0)&(values<=6))
    mixed=simulate_scores([2,3,8],.5,100,np.random.default_rng(100),score='mixed')
    assert np.all((mixed>=0)&(mixed<=3))


@pytest.mark.parametrize('alpha',[.05,.025])
def test_one_pair_analytic_reference(alpha):
    assert all_pairs_reference([1],[2],alpha=alpha)['upper']==pytest.approx(math.cos(math.pi*alpha),abs=1e-14)


def test_pair_simulation_matches_analytic_mean():
    rho=.7; n=2000
    values=simulate_scores([2]*3,rho,n,np.random.default_rng(20260911))
    assert abs(values.mean()-3*math.acos(rho)/math.pi)<.05


@pytest.mark.parametrize('kwargs',[{'B':0},{'B':2.5},{'seed':-1},{'grid':[.1,1.]},{'grid':[0.,.5,.5,1.]},{'alpha':.1},{'score':'adaptive'}])
def test_invalid_protocol(kwargs):
    params=dict(seed=1,B=19,grid=[0.,.5,1.]); params.update(kwargs)
    with pytest.raises(ValueError): mc_upper_bound([1],[2],**params)

def test_rank_superuniformity_by_exact_bernoulli_enumeration():
    # Exact finite-B distribution for one mixed-pair score under a dominated law.
    from scipy.stats import binom
    B=19; null=.3
    for actual in [0.,.1,.3]:
        for alpha in [.05,.1,.25,.5]:
            rejection=0.
            for observed in [0,1]:
                observed_prob=actual if observed else 1-actual
                for exceed in range(B+1):
                    # If score observed=0, every simulation ties or exceeds.
                    prob=float(binom.pmf(exceed,B,null)) if observed else float(exceed==B)
                    p=(1+exceed)/(B+1)
                    rejection+=observed_prob*prob*(p<=alpha)
            assert rejection<=alpha+1e-13


def test_hull_failure_implies_floor_grid_rejection_exhaustively():
    import itertools
    grid=[0.,.25,.5,.75,1.]
    for accepted in itertools.product([False,True],repeat=len(grid)):
        p=[1. if x else .01 for x in accepted]
        upper=invert_grid_pvalues(grid,p,.05)
        for truth in [.01,.25,.38,.5,.99,1.]:
            floor=max(i for i,r in enumerate(grid) if r<=truth)
            if upper<truth:
                assert not accepted[floor]
