"""Parametric bootstrap from the unrestricted fit, with a profile inversion.

The comparison method of the exact-readings campaign (`calibrate.py`): the
critical value is taken from datasets simulated at the fitted model, and the
profile is inverted numerically. It undercovers near the boundary of the
parameter space, which is why the package calibrates its interval null-wise
(`amr_clonalshare.mic_null_bootstrap`) instead; it stays here so that the
campaign can be rerun as reported.
"""
from __future__ import annotations

from numbers import Integral

import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi2

from amr_clonalshare._mic_likelihood import GaussianMICLikelihood
from amr_clonalshare._mic_panels import _as_panels, _record_covariate, _record_offset, observe_panels
from amr_clonalshare.mic_inference import _alpha


def bootstrap_critical(statistics,alpha=.05):
    alpha=_alpha(alpha)
    t=np.asarray(statistics,dtype=float)
    if t.ndim!=1 or len(t)==0 or np.isnan(t).any() or np.any(t<0):
        raise ValueError('bootstrap statistics must be nonnegative, with +infinity for failures')
    k=int(np.ceil((len(t)+1)*(1-alpha)))
    cutoff=float(np.sort(t)[k-1]) if k<=len(t) else float('inf')
    return max(float(chi2.ppf(1-alpha,1)),cutoff)


def calibrate_profile(problem, panels, exact, *, n_boot, seed, alpha=.05, order=24):
    """Return full fit, critical value and every bootstrap LR statistic.

    ``panels`` is the pair ``_panels`` returns (one edge array per panel and the
    panel index of every record) or, for records read on one panel, its edge
    array."""
    panels=_as_panels(panels,problem.n)
    if isinstance(n_boot,bool) or not isinstance(n_boot,Integral) or n_boot<19:
        raise ValueError('n_boot must be an integer at least 19 (199 or more is preferable)')
    if seed is None or isinstance(seed,bool) or not isinstance(seed,Integral) or seed<0:
        raise ValueError('an explicit nonnegative integer seed is required')
    _alpha(alpha)
    fit=problem.fit(order=order)
    if (not fit.converged or fit.rho >= 1. or not np.isfinite(fit.mean)
            or not np.isfinite(fit.total_sd) or fit.total_sd <= 0.):
        # A limiting nuisance fit cannot be simulated.
        return fit,float('inf'),np.full(n_boot,np.inf)
    rng=np.random.default_rng(seed)
    statistics=np.full(n_boot,np.inf)
    for k in range(n_boot):
        effects=rng.normal(0,fit.total_sd*np.sqrt(fit.rho),problem.groups)
        latent=fit.mean+_record_offset(problem,fit)+effects[problem.code]+rng.normal(0,fit.total_sd*np.sqrt(1-fit.rho),problem.n)
        a,b=observe_panels(latent,panels)
        a[exact]=latent[exact]; b[exact]=latent[exact]
        try:
            boot=GaussianMICLikelihood(a,b,problem.code,covariate=_record_covariate(problem))
            unrestricted=boot.fit(start=fit,order=order,multistart=False)
            restricted=boot.fit(rho=fit.rho,start=unrestricted,order=order,multistart=False)
            if not unrestricted.converged or unrestricted.loglik<restricted.loglik-1e-5:
                unrestricted=boot.fit(start=restricted,order=order,multistart=True)
            if (not unrestricted.converged or not restricted.converged
                    or unrestricted.loglik<restricted.loglik-1e-5):
                continue
            statistics[k]=max(0.,2*(unrestricted.loglik-restricted.loglik))
        except (ValueError,ArithmeticError,FloatingPointError):
            continue
    return fit,bootstrap_critical(statistics,alpha),statistics


def profile_bounds(problem, fit, critical, *, order=24):
    """Numerically invert the full profile and close the parameter boundaries.

    Scan a fixed grid plus the fitted value before rooting the outer crossing
    points. Retain the hull if the accepted set has more than one component.
    Any unresolved profile evaluation returns [0,1], not a spuriously narrow
    interval. Root tolerances are expanded outward.
    """
    if not fit.converged or not np.isfinite(critical):
        return 0.,1.,1
    cache={float(fit.rho):(0.,fit)}
    failed=0
    def lr(rho):
        nonlocal failed
        key=float(rho)
        if key in cache:
            return cache[key][0]
        nearest=min(cache,key=lambda p:abs(p-key))
        restricted=problem.fit(rho=key,start=cache[nearest][1],order=order,multistart=False)
        if not restricted.converged or restricted.loglik>fit.loglik+1e-4:
            restricted=problem.fit(rho=key,order=order,multistart=True)
        if not restricted.converged or restricted.loglik>fit.loglik+1e-4:
            failed+=1
            raise ArithmeticError('profile fit unresolved or exceeds the unrestricted maximum')
        value=max(0.,2*(fit.loglik-restricted.loglik))
        cache[key]=(value,restricted)
        return value
    try:
        grid=np.unique(np.r_[0.,.1,.25,.5,.75,.9,.97,.995,fit.rho])
        values=np.array([lr(r) for r in grid])
        inside=np.flatnonzero(values<=critical)
        if not len(inside):
            return 0.,1.,failed+1
        first,last=int(inside[0]),int(inside[-1])
        low=0. if first==0 else brentq(lambda r:lr(r)-critical,grid[first-1],grid[first],xtol=1e-6)
        if last==len(grid)-1:
            high=1.
        else:
            high=brentq(lambda r:lr(r)-critical,grid[last],grid[last+1],xtol=1e-6)
        return max(0.,float(low)-2e-6),min(1.,float(high)+2e-6),failed
    except (ValueError,ArithmeticError,FloatingPointError):
        return 0.,1.,failed+1
