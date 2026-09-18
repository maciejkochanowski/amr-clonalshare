"""Optional Gaussian population-ICC profile likelihood for grouped binary data.

Standalone research prototype. a is the marginal probit intercept, p=Phi(a).
The population liability variance is one and rho is its random-intercept share.
No finite-collection membership interpretation or universal coverage guarantee.
"""
from __future__ import annotations
from functools import lru_cache
import math
import time
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import gammaln, log_ndtr, logsumexp, ndtr, ndtri, roots_hermitenorm, roots_legendre

LR95 = 3.841458820694124
GH_ORDER = 96
GL_ORDER = 256
MIN_REPEATED_GROUPS = 2
MAX_GROUP_SIZE = 1000000

@lru_cache(maxsize=64)
def _conditional_rule(m):
    order=max(GL_ORDER,int(16*math.ceil(32*math.sqrt(m+1)/16)))
    l,lw=roots_legendre(order)
    return 10*l,10*lw

@lru_cache(maxsize=32)
def _normal_rule(m):
    z,hw=roots_hermitenorm(max(GH_ORDER,4*m))
    logweights=np.full_like(hw,-np.inf)
    np.log(hw,out=logweights,where=hw>0)
    logweights-=.5*math.log(2*math.pi)
    return z,logweights

@lru_cache(maxsize=128)
def _conditional_basis(m,indices):
    k=np.asarray(indices)
    lc=gammaln(m+1)-gammaln(k+1)-gammaln(m-k+1)
    w,ww=_conditional_rule(m)
    basis=lc[:,None]+k[:,None]*log_ndtr(w)[None,:]+(m-k)[:,None]*log_ndtr(-w)[None,:]
    return k,lc,basis,w,ww

def count_probabilities(m, a, rho):
    """Full integrated count PMF for diagnostics (m<=10000).

    The likelihood uses selected_count_probabilities and does not materialize
    probabilities for every possible count in a large group.
    """
    if not (np.isfinite(m) and int(m)==m and 1<=m<=10000):
        raise ValueError('Full PMF requires integer 1<=m<=10000; use selected counts for larger groups')
    return selected_count_probabilities(int(m),a,rho,tuple(range(int(m)+1)))

def selected_count_probabilities(m,a,rho,indices):
    return np.exp(selected_count_log_probabilities(m,a,rho,indices))

def selected_count_log_probabilities(m,a,rho,indices):
    if not (np.isfinite(m) and int(m)==m and 1<=m<=MAX_GROUP_SIZE and np.isfinite(a) and 0<=rho<=1):
        raise ValueError('Require valid integer group size, finite marginal probit intercept, and 0<=rho<=1')
    m=int(m)
    raw=np.asarray(indices,dtype=float)
    if raw.ndim!=1 or not np.isfinite(raw).all() or not np.equal(raw,np.floor(raw)).all() or np.any(raw<0) or np.any(raw>m):
        raise ValueError('Selected counts must lie in 0..m')
    k=raw.astype(int)
    lc=gammaln(m+1)-gammaln(k+1)-gammaln(m-k+1)
    if m==1:
        return np.where(k==0,log_ndtr(-a),log_ndtr(a))
    if rho==0:
        return lc+k*log_ndtr(a)+(m-k)*log_ndtr(-a)
    if rho==1:
        return np.where(k==0,log_ndtr(-a),np.where(k==m,log_ndtr(a),-np.inf))
    tau=math.sqrt(rho/(1-rho))
    mu=a/math.sqrt(1-rho)
    if rho<=.25:
        z,log_hw=_normal_rule(m)
        eta=mu+tau*z
        return lc+logsumexp(k[:,None]*log_ndtr(eta)[None,:]
                       +(m-k)[:,None]*log_ndtr(-eta)[None,:]+log_hw[None,:],axis=1)
    # Integrate in conditional-probit coordinates to resolve rho near one.
    k,lc,basis,w,ww=_conditional_basis(m,tuple(map(int,k)))
    normal_weight=-.5*((w-mu)/tau)**2+np.log(ww)-math.log(tau*math.sqrt(2*math.pi))
    out=logsumexp(basis+normal_weight[None,:],axis=1)
    # Beyond +/-10 the conditional binomial is effectively constant.
    out[k==0]=np.logaddexp(out[k==0],log_ndtr((-10-mu)/tau))
    out[k==m]=np.logaddexp(out[k==m],log_ndtr((mu-10)/tau))
    return out

def validate_counts(counts, sizes):
    k=np.asarray(counts,dtype=float)
    m=np.asarray(sizes,dtype=float)
    if k.ndim!=1 or m.ndim!=1 or k.size!=m.size or k.size<1:
        raise ValueError('counts and sizes must be equal nonempty one-dimensional vectors')
    if not (np.isfinite(k).all() and np.isfinite(m).all() and
            np.equal(k,np.floor(k)).all() and np.equal(m,np.floor(m)).all() and
            (m>=1).all() and (m<=MAX_GROUP_SIZE).all() and (k>=0).all() and (k<=m).all()):
        raise ValueError(f'Require integer 1<=sizes<={MAX_GROUP_SIZE} and counts between zero and size')
    return k.astype(int),m.astype(int)

class ProfileLikelihood:
    """Compressed sufficient counts and memoized nuisance profiles."""
    def __init__(self,counts,sizes):
        self.k,self.m=validate_counts(counts,sizes)
        self.parts=[]
        for m in np.unique(self.m):
            counts,freq=np.unique(self.k[self.m==m],return_counts=True)
            self.parts.append((int(m),counts,freq))
        self.total=int(self.m.sum()); self.successes=int(self.k.sum())
        self.cache={}; self.n_eval=0

    def nll(self,a,rho):
        self.n_eval+=1
        value=0.
        for m,counts,freq in self.parts:
            logp=selected_count_log_probabilities(m,float(a),float(rho),counts)
            if np.any(~np.isfinite(logp)):
                return float('inf')
            value-=float(freq@logp)
        return value

    def profile(self,rho):
        rho=float(rho)
        if rho in self.cache:
            return self.cache[rho]
        if rho==0:
            a=float(ndtri(self.successes/self.total))
            result=(self.nll(a,0),a)
        elif rho==1:
            if np.any((self.k>0)&(self.k<self.m)):
                result=(float('inf'),float('nan'))
            else:
                a=float(ndtri(float((self.k==self.m).mean())))
                result=(self.nll(a,1),a)
        else:
            opt=minimize_scalar(lambda a:self.nll(a,rho),bounds=(-9.,9.),method='bounded',
                                options={'xatol':2e-7,'maxiter':100})
            if not opt.success or not np.isfinite(opt.fun) or abs(opt.x)>8.999:
                raise RuntimeError('Nuisance optimization failed')
            result=(float(opt.fun),float(opt.x))
        self.cache[rho]=result
        return result

    def mle(self):
        # Refine every interval, including those whose endpoints do not
        # bracket a grid-local minimum. Retain exact boundary candidates.
        # Finite bracketing is a numerical diagnostic, not a global theorem.
        grid=np.array([0.,.0625,.125,.1875,.25,.35,.5,.65,.8,.9,.98,1.])
        vals=np.array([self.profile(r)[0] for r in grid])
        candidates=[(float(v),float(r)) for v,r in zip(vals,grid)]
        for lower,upper in zip(grid[:-1],grid[1:]):
            opt=minimize_scalar(lambda r:self.profile(r)[0],bounds=(lower,upper),
                                method='bounded',options={'xatol':2e-9,'maxiter':150})
            if not opt.success:
                raise RuntimeError('ICC interval optimization failed')
            candidates.append((float(opt.fun),float(opt.x)))
        best,rho=min(candidates)
        return rho,self.profile(rho)[1],best

def fit_profile(counts,sizes,*,critical_value=LR95,truth_rho=None,compute_interval=True):
    """Profile-LR set for population rho; critical value must be prespecified.

    All-singleton data contain no rho information. A minimum of two groups
    with repeat observations is a declared reporting policy, not a theorem.
    Constant overall outcomes retain the uninformative identified set [0,1].
    """
    start=time.perf_counter()
    k,m=validate_counts(counts,sizes)
    if not np.isfinite(critical_value) or critical_value<=0:
        raise ValueError('critical_value must be positive and finite')
    if truth_rho is not None and (not np.isfinite(truth_rho) or not 0<=truth_rho<=1):
        raise ValueError('truth_rho must be finite and in [0,1]')
    repeated=int((m>=2).sum())
    total=int(m.sum()); successes=int(k.sum())
    base=dict(n=total,n_groups=len(k),n_repeated_groups=repeated,
              repeat_isolate_support=float(m[m>=2].sum()/total),critical_value=float(critical_value))
    if successes in (0,total) or repeated<MIN_REPEATED_GROUPS:
        status='constant_outcome_unidentified' if successes in (0,total) else 'insufficient_repeated_groups'
        return dict(base,rho_hat=None,prevalence_hat=successes/total,ci_low=0.,ci_high=1.,
                    identified=False,informative=False,status=status,lr_at_truth=0. if truth_rho is not None else None,
                    nll=None,n_eval=0,seconds=time.perf_counter()-start)
    ll=ProfileLikelihood(k,m)
    rho,a,best=ll.mle()
    def score(r):
        return max(0.,2*(ll.profile(r)[0]-best))-critical_value
    low=high=None
    topology='not_checked'
    if compute_interval:
        low=0. if score(0.)<=0 else brentq(score,0.,rho,xtol=5e-15,rtol=9e-16,maxiter=150)
        high=1. if score(1.)<=0 else brentq(score,rho,1.,xtol=5e-15,rtol=9e-16,maxiter=150)
        for endpoint in (low,high):
            if 0.<endpoint<1. and abs(score(endpoint))>1e-5:
                raise RuntimeError('Profile endpoint likelihood-ratio residual exceeds tolerance')
        # Guard against obvious disconnected LR regions and an overlooked
        # profile minimum. A finite diagnostic grid is not a proof of global
        # topology; return a numerical failure instead of hiding a detected one.
        diagnostic=np.unique(np.r_[np.linspace(0.,1.,17),rho,
                                   (low+rho)/2,(high+rho)/2])
        for r in diagnostic:
            value=score(float(r))
            expected_inside=low<=r<=high
            if (value<=-1e-5 and not expected_inside) or (value>1e-5 and expected_inside):
                raise RuntimeError('Disconnected or numerically inconsistent profile likelihood set')
            if ll.profile(float(r))[0]<best-1e-5:
                raise RuntimeError('Profile diagnostic found a better maximum')
        topology='17_point_profile_scan_passed'
    lr=None if truth_rho is None else max(0.,2*(ll.profile(float(truth_rho))[0]-best))
    return dict(base,rho_hat=float(rho),prevalence_hat=float(ndtr(a)),ci_low=low,ci_high=high,
                identified=True,informative=bool(low is not None and high-low<.999999),status='ok',
                lr_at_truth=lr,nll=float(best),n_eval=ll.n_eval,profile_topology_check=topology,
                seconds=time.perf_counter()-start)
