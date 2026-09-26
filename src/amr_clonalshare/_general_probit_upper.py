"""Finite-Monte-Carlo Gaussian-liability upper bound on a prespecified grid.

Probability theory covers the Monte Carlo rank construction
averaged over its randomness under the stated model. No ICC point estimate.
"""
from __future__ import annotations
import math
import platform
import time
import numpy as np
import scipy
from scipy.special import betaincinv, ndtr
from scipy.stats import binom

DEFAULT_GRID=tuple(float(x) for x in np.linspace(0.,1.,101))


def _design(counts,sizes):
    try:
        k=np.asarray(counts,dtype=float); m=np.asarray(sizes,dtype=float)
    except (TypeError,ValueError,OverflowError) as exc: raise ValueError('Require integer count/size vectors') from exc
    if k.ndim!=1 or m.ndim!=1 or len(k)!=len(m) or not len(m): raise ValueError('Require matching nonempty one-dimensional vectors')
    if not np.isfinite(k).all() or not np.isfinite(m).all() or not np.equal(k,np.floor(k)).all() or not np.equal(m,np.floor(m)).all():
        raise ValueError('Counts and sizes must be finite integers')
    if np.any(m<1) or np.any(m>10000000) or np.any(k<0) or np.any(k>m): raise ValueError('Require 1<=size<=10000000 and 0<=count<=size')
    return k.astype(np.int64),m.astype(np.int64)


def _grid(values):
    try: grid=np.asarray(values,dtype=float)
    except (TypeError,ValueError,OverflowError) as exc: raise ValueError('Invalid grid') from exc
    if grid.ndim!=1 or len(grid)<2 or not np.isfinite(grid).all() or grid[0]!=0. or grid[-1]!=1. or np.any(np.diff(grid)<=0):
        raise ValueError('Prespecified grid must increase strictly from 0 to 1')
    return grid


def invert_grid_pvalues(grid,pvalues,alpha):
    """Upper hull of every accepted cell; preserves late/nonmonotone acceptances."""
    grid=_grid(grid); p=np.asarray(pvalues,dtype=float)
    if p.shape!=grid.shape or not np.isfinite(p).all() or np.any(p<0) or np.any(p>1): raise ValueError('Invalid grid p-values')
    accepted=np.flatnonzero(p>alpha)
    return 0. if len(accepted)==0 else float(grid[min(int(accepted[-1])+1,len(grid)-1)])


def simulate_scores(sizes,rho,B,rng,*,score='centrality',batch_size=None):
    """IID least-favourable a=0 counts, with fixed sizes; no outcome filtering."""
    _,sizes=_design(np.zeros(len(sizes)),sizes)
    sizes=np.sort(sizes[sizes>=2])
    if score not in ('centrality','mixed'): raise ValueError('Unknown predeclared score')
    if not np.isfinite(rho) or not 0<=rho<=1: raise ValueError('rho must lie in [0,1]')
    if isinstance(B,bool) or not isinstance(B,(int,np.integer)) or B<1: raise ValueError('B must be positive integer')
    out=np.zeros(B,dtype=np.int64)
    if rho==1. or not len(sizes): return out
    if batch_size is None: batch_size=max(1,min(256,1000000//len(sizes)))
    if not isinstance(batch_size,(int,np.integer)) or batch_size<1: raise ValueError('Invalid batch size')
    for begin in range(0,B,batch_size):
        end=min(B,begin+batch_size)
        if rho==0.: p=.5
        else:
            z=rng.standard_normal((end-begin,len(sizes)))
            p=ndtr(math.sqrt(rho/(1-rho))*z)
        counts=rng.binomial(sizes,p,size=(end-begin,len(sizes)))
        if score=='centrality': out[begin:end]=np.minimum(counts,sizes-counts).sum(axis=1)
        else: out[begin:end]=((counts>0)&(counts<sizes)).sum(axis=1)
    return out


def mc_upper_bound(counts,sizes,*,seed,B=9999,grid=DEFAULT_GRID,alpha=.05,score='centrality'):
    """Return only [0,U]; grid/B/score/seed protocol must not be outcome-selected.

    The inclusive rank p=(1+#null scores >= observed)/(B+1) is used without
    mid-p or adaptive stopping. Each grid point has a separate spawned stream.
    """
    started=time.perf_counter(); counts,sizes=_design(counts,sizes); grid=_grid(grid)
    if isinstance(seed,bool) or not isinstance(seed,(int,np.integer)) or seed<0: raise ValueError('seed must be nonnegative integer')
    if isinstance(B,bool) or not isinstance(B,(int,np.integer)) or B<1: raise ValueError('B must be positive integer')
    if alpha not in (.05,.025): raise ValueError('Predeclared alpha must be .05 or .025')
    if score not in ('centrality','mixed'): raise ValueError('Unknown predeclared score')
    observed=int(np.minimum(counts,sizes-counts).sum()) if score=='centrality' else int(((counts>0)&(counts<sizes)).sum())
    repeated=int((sizes>=2).sum()); pvalues=[]; exceedances=[]; simulated=0
    streams=np.random.SeedSequence(int(seed)).spawn(len(grid))
    for rho,stream in zip(grid,streams):
        if observed==0:
            exceed=B
        elif rho==1.:
            # Exactly zero centrality at the limiting all-or-none model.
            exceed=0
        else:
            values=simulate_scores(sizes,float(rho),B,np.random.default_rng(stream),score=score)
            exceed=int(np.count_nonzero(values>=observed)); simulated+=B
        exceedances.append(exceed); pvalues.append((1+exceed)/(B+1))
    upper=invert_grid_pvalues(grid,pvalues,alpha)
    return dict(method='gaussian_liability_'+score+'_finite_mc_upper',lower=0.,upper=upper,estimate=None,
                score=score,observed_score=observed,n=int(sizes.sum()),n_groups=len(sizes),n_repeated_groups=repeated,
                alpha=alpha,confidence_level=1-alpha,seed=int(seed),B=int(B),grid=grid.tolist(),
                max_grid_cell_width=float(np.diff(grid).max()),pvalues=pvalues,exceedances=exceedances,
                grid_acceptance=[p>alpha for p in pvalues],simulated_replicates=simulated,
                estimated_binomial_variates=simulated*repeated,
                status='no_upper_restriction' if observed==0 else ('empty_grid_acceptance_reported_as_zero' if not any(p>alpha for p in pvalues) else 'ok'),
                monte_carlo_rejects_at_zero=pvalues[0]<=alpha,
                tail_convention='inclusive >= tail, plus one numerator/denominator; reject p<=alpha',
                inversion='upper hull of next grid endpoints for every accepted grid point; zero if none',
                coverage_scope='finite-Monte-Carlo coverage averaged over algorithm randomness under the stated model; not conditional on a fixed deterministic seed',
                numerical_scope='Gaussian/binomial pseudorandom simulation in floating point; no quadrature or plug-in prevalence',
                assumptions='independent Gaussian random group effects and binomial counts; fixed or noninformative sizes; predeclared score, grid, B and RNG protocol; no outcome-selected groups',
                joint_coverage='standalone one-sided bound; no intersection or joint confidence claim with another procedure',
                reproducibility=dict(rng='PCG64 via SeedSequence.spawn per fixed grid index',python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sorted_sizes=sorted(map(int,sizes))),
                seconds=time.perf_counter()-started)


def all_pairs_reference(counts,sizes,*,alpha=.05):
    """Numerical analytic reference: mixed pairs are Bernoulli(acos(rho)/pi)."""
    counts,sizes=_design(counts,sizes)
    if alpha not in (.05,.025) or np.any((sizes!=1)&(sizes!=2)): raise ValueError('Reference requires sizes 1/2 and alpha .05/.025')
    groups=int((sizes==2).sum()); observed=int(((sizes==2)&(counts==1)).sum())
    if observed==0: upper=1.; pzero=1.
    else:
        pzero=float(binom.sf(observed-1,groups,.5))
        q=float(betaincinv(observed,groups-observed+1,alpha))
        upper=0. if pzero<=alpha else max(0.,math.cos(math.pi*q))
    return dict(method='analytic_all_pairs_reference',lower=0.,upper=upper,estimate=None,alpha=alpha,n_pairs=groups,observed_mixed=observed,p_at_zero=pzero)
