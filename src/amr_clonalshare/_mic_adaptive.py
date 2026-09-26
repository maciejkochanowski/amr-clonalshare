"""Checked adaptive integration when Gauss-Hermite cannot resolve a tail.

Normalisation and analytic scores are integrated together after log-density
centring, over the full real line. No prior-tail truncation is introduced.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import quad_vec
from ._normal_numerics import normal_interval_logmass, normal_truncated_moments


def integrate_lineages_reference(problem, counts, modes, scales, mean, sigma, tau2, *, strict=False):
    ll=dm=ds=dt=0.
    total_error=0.
    sig2=sigma*sigma
    for nn,mode,scale in zip(counts,modes,scales):
        present=nn>0
        a,b,n,exact=problem.a[present],problem.b[present],nn[present],problem.exact[present]
        def pieces(u, moments):
            eta=mean+u
            masses=normal_interval_logmass((a-eta)/sigma,(b-eta)/sigma)
            if exact.any():
                masses=np.where(exact,-.5*((a-eta)/sigma)**2-np.log(sigma)-.5*np.log(2*np.pi),masses)
            logweight=float(np.dot(n,masses)-u*u/(2*tau2))
            if not moments:
                return logweight
            ey,vy=normal_truncated_moments(a,b,eta,sigma)
            delta=ey-eta
            return logweight,np.array([1.,np.dot(n,delta)/sig2,
                np.dot(n,(vy+delta*delta)/sig2-1.),u*u/tau2-1.])
        peak=pieces(mode,False)
        def integrand(z):
            u=mode+scale*z
            logweight=pieces(u,False)-peak
            if logweight < -740. or not np.isfinite(logweight):
                return np.zeros(4)
            _,scores=pieces(u,True)
            return np.exp(logweight)*scores
        eps=2e-11 if strict else 2e-9
        integral,error,info=quad_vec(integrand,-np.inf,np.inf,epsabs=eps,epsrel=eps,
            norm='max',limit=400,full_output=True)
        if (not info.success or not np.isfinite(integral).all() or integral[0]<=0
                or error>1e-6*max(1.,float(np.max(np.abs(integral))))):
            raise ArithmeticError('adaptive lineage likelihood or score failed precision check')
        ll+=peak+np.log(scale)-.5*np.log(2*np.pi*tau2)+np.log(integral[0])
        dm+=integral[1]/integral[0]
        ds+=integral[2]/integral[0]
        dt+=integral[3]/integral[0]
        total_error+=float(error/integral[0])
    return float(ll),float(dm),float(ds),float(dt),total_error


def integrate_lineages(problem, counts, modes, scales, mean, sigma, tau2, *, strict=False):
    """Use checked composite integration; retain independent adaptive fallback.

    The composite rule also returns the mean score of every bin, from which the
    fixed-effect scores follow; the fallback does not, and returns None there.
    """
    from ._mic_piecewise import integrate_lineages as composite
    arguments = (problem, counts, modes, scales, mean, sigma, tau2)
    try:
        return composite(*arguments, strict=strict)
    except ArithmeticError:
        return (*integrate_lineages_reference(*arguments, strict=strict), None)
