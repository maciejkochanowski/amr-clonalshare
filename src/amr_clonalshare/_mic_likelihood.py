"""Fully profiled Gaussian MIC likelihood using grouped observation counts.

Adaptive Gauss-Hermite nodes follow each lineage's conditional mode. Scores
are posterior expectations of complete-data derivatives. Final fits check a
higher quadrature order. This numerical procedure does not itself establish
coverage of likelihood-ratio confidence limits.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache

import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar
from scipy.special import logsumexp, roots_hermitenorm

from .attribution import _codes, _is_untyped
from ._normal_numerics import normal_interval_logmass, normal_truncated_moments

_LOG2PI = float(np.log(2*np.pi))


@lru_cache(maxsize=8)
def _nodes(order):
    x, w = roots_hermitenorm(order)
    with np.errstate(divide='ignore'):
        return x, np.log(w)-.5*_LOG2PI


@dataclass(frozen=True)
class MICFit:
    mean: float
    total_sd: float
    rho: float
    loglik: float
    converged: bool
    #: likelihood evaluations spent on this fit
    evaluations: int
    #: Gauss-Hermite order of the final evaluation; 0 when no quadrature was
    #: needed (closed form) or when the adaptive integrator took over, which
    #: ``message`` tells apart
    quadrature_order: int
    quadrature_error: float
    message: str
    #: fixed effects of the covariate levels after the first, in level order;
    #: empty without a covariate
    coefficients: tuple = ()

    def as_dict(self):
        return asdict(self)


def _covariate_columns(covariate):
    """Normalise ``covariate`` to a list of per-reading level vectors."""
    if covariate is None:
        return []
    if isinstance(covariate, np.ndarray) and covariate.ndim == 2:
        return [covariate[:, k] for k in range(covariate.shape[1])]
    if isinstance(covariate, (list, tuple)) and len(covariate) and all(
            isinstance(c, (list, tuple, np.ndarray)) for c in covariate):
        return list(covariate)
    return [covariate]


class _ShiftedBins:
    """The bins of a likelihood shifted by the fixed effect of their level."""
    def __init__(self, a, b, exact):
        self.a, self.b, self.exact = a, b, exact


class GaussianMICLikelihood:
    """Gaussian lineage model for interval readings, with optional categorical
    fixed effects.

    ``covariate`` names, per reading, the level of one categorical covariate
    (a laboratory, a country), or is a sequence of such vectors, one per
    covariate. The model is then y = mean + sum of beta[level] over the
    covariates + lineage effect + residual, with beta of the first level of
    each covariate, in sorted order, fixed at zero. The fixed effects are
    estimated jointly with the mean, the scale and rho, so that a shift
    between levels is removed from both variance components rather than
    absorbed into whichever of them it is aligned with.
    """
    def __init__(self, lo, hi, lineage, covariate=None):
        a, b = np.asarray(lo,dtype=float), np.asarray(hi,dtype=float)
        labels = np.asarray(lineage,dtype=object)
        if a.ndim != 1 or a.shape != b.shape or a.shape != labels.shape:
            raise ValueError('matching one-dimensional interval and label arrays are required')
        columns = _covariate_columns(covariate)
        self.levels = []          # level names of every covariate
        level_columns = []
        for raw in columns:
            raw = np.asarray(raw, dtype=object).ravel()
            if raw.shape != a.shape:
                raise ValueError('every covariate must name a level for every reading')
            if any(_is_untyped(v) for v in raw):
                raise ValueError('missing covariate levels must be handled explicitly')
            names = sorted({str(v) for v in raw})
            level_columns.append(np.searchsorted(names, np.asarray([str(v) for v in raw])))
            self.levels.append(tuple(names))
        self.levels = tuple(self.levels)
        self.n_levels = tuple(len(names) for names in self.levels)
        self.n_coefficients = int(sum(max(n-1, 0) for n in self.n_levels))
        # per-reading level of every covariate (readings x covariates)
        self.record_level = (np.column_stack(level_columns) if level_columns
                             else np.zeros((a.size, 0), dtype=int))
        if a.size < 3 or np.isnan(a).any() or np.isnan(b).any() or np.any(a>b):
            raise ValueError('at least three ordered intervals without NaN are required')
        if np.isposinf(a).any() or np.isneginf(b).any() or np.any(np.isneginf(a)&np.isposinf(b)):
            raise ValueError('infinite exact or wholly uninformative intervals are not supported')
        if any(_is_untyped(v) for v in labels):
            raise ValueError('missing lineage labels must be handled explicitly')
        self.lo, self.hi = a.copy(), b.copy()
        self.code = _codes(labels)
        self.sizes = np.bincount(self.code).astype(float)
        self.groups, self.n = len(self.sizes), len(a)
        if self.groups<2 or self.n<=self.groups:
            raise ValueError('at least two lineages and within-lineage replication are required')
        bins, inverse = np.unique(np.column_stack([a,b,self.record_level]),axis=0,return_inverse=True)
        inverse = inverse.ravel()   # NumPy 2.0.0 returns it as a column
        self.a, self.b = bins[:,0], bins[:,1]
        self.level = bins[:,2:].astype(int)   # bins x covariates
        self.offset = np.zeros(len(bins))
        self.exact = np.isfinite(self.a)&(self.a==self.b)
        self.counts = np.zeros((self.groups,len(bins)))
        np.add.at(self.counts,(self.code,inverse),1)
        self.all_exact = bool(np.all(a==b))
        # the closed form of the one-way model needs one mean for all readings
        self.closed_form = self.all_exact and self.n_coefficients == 0
        for column, n in enumerate(self.n_levels):
            for k in range(n):
                rows = self.record_level[:,column] == k
                if np.isneginf(a[rows]).all() or np.isposinf(b[rows]).all():
                    raise ValueError('covariate level %r is censored on one side in every reading; '
                                     'its effect is not identified' % (self.levels[column][k],))
        finite = np.r_[a[np.isfinite(a)],b[np.isfinite(b)]]
        if np.ptp(finite)==0:
            raise ValueError('a single cut point or constant observations do not identify the full scale model')
        proxy = np.where(np.isfinite(a)&np.isfinite(b),
                         (np.where(np.isfinite(a),a,0)+np.where(np.isfinite(b),b,0))/2,
                         np.where(np.isfinite(a),a,b))
        width = b-a; width = width[np.isfinite(width)&(width>0)]
        self.unit = max(float(np.std(proxy)),float(np.median(width)) if len(width) else 0.,np.finfo(float).eps)
        self.center = float(np.mean(proxy))
        self.group_proxy = np.bincount(self.code,weights=proxy)/self.sizes
        self.means = np.bincount(self.code,weights=a)/self.sizes if self.all_exact else None
        self.sse = float(np.sum((a-self.means[self.code])**2)) if self.all_exact else None
        self.evaluations = 0
        self._cache = None

    def _cells(self, mean, sd, with_moments=True):
        """Log masses and truncated moments of every (lineage, bin, node) cell.

        ``mean`` is the lineage-level latent mean (groups, nodes); the fixed
        effect of each bin's covariate level is added here, so every caller
        receives ``eta`` as the full cell mean.
        """
        aa, bb = self.a[None,:,None], self.b[None,:,None]
        eta = mean[:,None,:]+self.offset[None,:,None]
        lm = normal_interval_logmass((aa-eta)/sd,(bb-eta)/sd)
        exact = self.exact[None,:,None]
        density = -.5*((np.where(exact,bb,0)-eta)/sd)**2-np.log(sd)-.5*_LOG2PI
        lm = np.where(exact,density,lm)
        if not with_moments:
            return lm
        ey, vy = normal_truncated_moments(aa,bb,eta,sd)
        return lm, ey, vy, eta

    def _by_level(self, per_bin):
        """Sum a per-bin mean score over the bins of each non-reference level,
        covariate by covariate, in the order of the coefficient vector."""
        if not self.n_coefficients:
            return np.zeros(0)
        return np.concatenate([np.bincount(self.level[:,column], weights=per_bin, minlength=n)[1:]
                               for column, n in enumerate(self.n_levels)])

    def _split_coefficients(self, coefficients):
        """The coefficient vector as one array per covariate, reference level included as zero."""
        values = np.asarray(coefficients, dtype=float).ravel()
        if values.size != self.n_coefficients or not np.isfinite(values).all():
            raise ValueError('one finite coefficient per covariate level after the first is required')
        betas, start = [], 0
        for n in self.n_levels:
            betas.append(np.r_[0., values[start:start+n-1]])
            start += n-1
        return betas

    def set_coefficients(self, coefficients):
        if not self.n_coefficients:
            self.offset = np.zeros(len(self.level))
            return
        betas = self._split_coefficients(coefficients)
        self.offset = np.sum([beta[self.level[:,column]] for column, beta in enumerate(betas)], axis=0)

    def evaluate(self, mean, total_sd, rho, *, order=24, coefficients=(), beta_gradient=True):
        """Log likelihood and its gradient in (mean, log total_sd, rho, coefficients...).

        ``order`` <= 0 selects the checked adaptive integration; with a
        covariate its gradient in the coefficients is taken by central
        differences of the adaptive log likelihood, which ``beta_gradient``
        switches off inside those differences.
        """
        if not np.isfinite(mean) or not np.isfinite(total_sd) or total_sd<=0 or not 0<=rho<1:
            raise ValueError('finite mean, positive total_sd and 0 <= rho < 1 required')
        self.set_coefficients(coefficients)
        self.evaluations += 1
        v=total_sd**2; sig2=(1-rho)*v; tau2=rho*v; sd=np.sqrt(sig2)
        if self.closed_form:
            d=self.means-mean; vg=sig2+self.sizes*tau2
            ll=-.5*(self.n*_LOG2PI+(self.n-self.groups)*np.log(sig2)
                      +np.sum(np.log(vg))+self.sse/sig2+np.sum(self.sizes*d*d/vg))
            dm=np.sum(self.sizes*d/vg)
            ds=-(self.n-self.groups)+self.sse/sig2-np.sum(sig2/vg)+np.sum(self.sizes*d*d*sig2/vg**2)
            dt=-np.sum(self.sizes*tau2/vg)+np.sum(self.sizes**2*d*d*tau2/vg**2)
            if rho>0:
                dr=.5*(dt/rho-ds/(1-rho))
            else:
                d_tau=.5*np.sum(self.sizes**2*d*d/sig2**2-self.sizes/sig2)
                dr=v*d_tau-.5*ds
            return float(ll),np.array([dm,ds+dt,dr])
        if rho==0:
            lm,ey,vy,eta=self._cells(np.full((self.groups,1),mean),sd)
            diff=ey[:,:,0]-eta[:,:,0]
            per_bin=np.sum(self.counts*diff/sig2,axis=0)
            score=np.sum(self.counts*diff/sig2,axis=1)
            hessian=np.sum(self.counts*(vy[:,:,0]/sig2-1)/sig2,axis=1)
            ds=np.sum(self.counts*((vy[:,:,0]+diff**2)/sig2-1))
            dr=.5*v*np.sum(score**2+hessian)-.5*ds
            return float(np.sum(self.counts*lm[:,:,0])),np.r_[score.sum(),ds,dr,self._by_level(per_bin)]
        multiple=self.sizes>1
        counts=self.counts[multiple]
        sizes=self.sizes[multiple]
        group_proxy=self.group_proxy[multiple]
        singleton_counts=self.counts[~multiple].sum(axis=0)
        singleton_ll=singleton_dm=singleton_ds=0.
        singleton_bin=np.zeros(len(self.level))
        if singleton_counts.sum():
            slm,sey,svy,seta=self._cells(np.array([[mean]]),total_sd)
            delta=sey[0,:,0]-seta[0,:,0]
            singleton_ll=float(np.dot(singleton_counts,slm[0,:,0]))
            singleton_bin=singleton_counts*delta/v
            singleton_dm=float(singleton_bin.sum())
            singleton_ds=float(np.dot(singleton_counts,(svy[0,:,0]+delta**2)/v-1))
        mode=(group_proxy-mean)*tau2/(tau2+sig2/sizes)
        for _ in range(80):
            mode_lm,ey,vy,eta=self._cells((mean+mode)[:,None],sd)
            score=np.sum(counts*(ey[:,:,0]-eta[:,:,0])/sig2,axis=1)-mode/tau2
            precision=1/tau2+np.sum(counts*(1-vy[:,:,0]/sig2)/sig2,axis=1)
            step=score/precision
            step=np.clip(step,-2*total_sd,2*total_sd)
            if np.max(np.abs(step)/(1+np.abs(mode)))<1e-9:
                break
            before=np.sum(counts*mode_lm[:,:,0],axis=1)-mode*mode/(2*tau2)
            for _backtrack in range(40):
                candidate=mode+step
                lm_candidate=self._cells((mean+candidate)[:,None],sd,with_moments=False)
                after=np.sum(counts*lm_candidate[:,:,0],axis=1)-candidate*candidate/(2*tau2)
                bad=after<before-1e-11
                if not bad.any():
                    break
                step[bad]*=.5
            else:
                raise ArithmeticError('conditional mode line search failed')
            mode=candidate
        else:
            # Each conditional log density is strictly concave. A bracketed
            # score root is slower but avoids declaring a valid parameter
            # point impossible when a Newton step approaches a panel edge.
            for g in range(len(sizes)):
                present=counts[g]>0
                aa,bb,nn=self.a[present],self.b[present],counts[g,present]
                shift=self.offset[present]
                def score_at(u):
                    ey,_=normal_truncated_moments(aa,bb,mean+shift+u,sd)
                    return float(np.dot(nn,ey-mean-shift-u)/sig2-u/tau2)
                bound=max(total_sd,abs(group_proxy[g]-mean)+total_sd)
                for _expand in range(60):
                    if score_at(-bound)>=0 and score_at(bound)<=0:
                        break
                    bound*=2
                else:
                    raise ArithmeticError('conditional mode score could not be bracketed')
                mode[g]=brentq(score_at,-bound,bound,xtol=1e-12,rtol=1e-14)
            _,_,vy,_=self._cells((mean+mode)[:,None],sd)
            precision=1/tau2+np.sum(counts*(1-vy[:,:,0]/sig2)/sig2,axis=1)
        scale=1/np.sqrt(precision)
        if order <= 0:
            from ._mic_adaptive import integrate_lineages
            # A fixed effect shifts the bins of its level: integrating the shifted
            # bins at the lineage mean is the same integral, so the adaptive
            # integrators need no covariate of their own.
            shifted=_ShiftedBins(self.a-self.offset,self.b-self.offset,self.exact)
            ll,dm,ds,dt,error=integrate_lineages(shifted,counts,mode,scale,mean,sd,tau2,strict=order<0)
            dr=.5*(dt/rho-ds/(1-rho))
            gradient=[dm+singleton_dm,ds+dt+singleton_ds,dr]
            if self.n_coefficients and beta_gradient:
                beta=np.asarray(coefficients,dtype=float).ravel()
                for k in range(self.n_coefficients):
                    h=1e-5*max(1.,abs(beta[k]))
                    plus,minus=beta.copy(),beta.copy(); plus[k]+=h; minus[k]-=h
                    up=self.evaluate(mean,total_sd,rho,order=order,coefficients=plus,beta_gradient=False)[0]
                    down=self.evaluate(mean,total_sd,rho,order=order,coefficients=minus,beta_gradient=False)[0]
                    gradient.append((up-down)/(2*h))
                self.set_coefficients(coefficients)
            elif self.n_coefficients:
                gradient.extend([np.nan]*self.n_coefficients)
            return ll+singleton_ll,np.array(gradient)
        x,lw=_nodes(int(order))
        u=mode[:,None]+scale[:,None]*x
        lm,ey,vy,eta=self._cells(mean+u,sd)
        terms=np.sum(counts[:,:,None]*lm,axis=1)+lw-u*u/(2*tau2)+.5*x*x
        lz=logsumexp(terms,axis=1)
        weights=np.exp(terms-lz[:,None])
        ll=np.sum(lz+np.log(scale/np.sqrt(tau2)))
        diff=ey-eta
        per_bin=np.sum(weights[:,None,:]*counts[:,:,None]*diff/sig2,axis=(0,2))
        dm=per_bin.sum()
        ds=np.sum(weights*np.sum(counts[:,:,None]*((vy+diff**2)/sig2-1),axis=1))
        dt=np.sum(weights*(u*u/tau2-1))
        dr=.5*(dt/rho-ds/(1-rho))
        return float(ll+singleton_ll),np.r_[dm+singleton_dm,ds+dt+singleton_ds,dr,
                                            self._by_level(per_bin+singleton_bin)]

    def objective(self,x,*,order=24):
        key=(tuple(np.asarray(x,dtype=float)),int(order))
        if self._cache is not None and self._cache[0]==key:
            return self._cache[1]
        value,gradient=self.evaluate(float(x[0]),float(np.exp(x[1])),float(x[2]),order=order,
                                     coefficients=x[3:])
        result=(-value,-gradient)
        self._cache=(key,result)
        return result

    def _exact_fit(self,rho):
        before=self.evaluations
        def at(r):
            lam=r/(1-r); w=self.sizes/(1+self.sizes*lam)
            mu=float(np.dot(w,self.means)/w.sum())
            ss=self.sse+np.dot(w,(self.means-mu)**2)
            sig2=ss/self.n
            ll=-.5*(self.n*(_LOG2PI+1+np.log(sig2))+np.log1p(self.sizes*lam).sum())
            return float(ll),mu,float(np.sqrt(sig2/(1-r)))
        if self.sse<=0:
            raise ValueError('positive within-lineage variation is required')
        if rho is None:
            candidates=[(0.,at(0.))]
            # Bracket local optima on a deterministic partition, retaining the boundary.
            grid=np.array([0,.1,.3,.6,.85,.97,.9999999])
            for left,right in zip(grid[:-1],grid[1:]):
                fit=minimize_scalar(lambda r:-at(r)[0],bounds=(left,right),method='bounded',options={'xatol':1e-10})
                candidates.append((float(fit.x),at(float(fit.x))))
            rho,parts=max(candidates,key=lambda t:t[1][0])
        else:
            parts=at(rho)
        ll,mu,sd=parts
        return MICFit(mu,sd,float(rho),ll,True,self.evaluations-before,0,0.,'analytic Gaussian likelihood')

    def fit(self,*,rho=None,start=None,order=24,verify=True,multistart=True,box_maximum=False):
        if rho is not None and (not np.isfinite(rho) or not 0<=rho<1):
            raise ValueError('fixed rho must lie in [0,1)')
        if self.closed_form:
            return self._exact_fit(rho)
        if rho is None and not self.n_coefficients:
            from ._mic_boundary import monomorphic_boundary_fit
            boundary = monomorphic_boundary_fit(self)
            if boundary is not None:
                return boundary
        p=self.n_coefficients
        lower=[self.center-20*self.unit,np.log(self.unit)-6,1e-7]+[-20*self.unit]*p
        upper=[self.center+20*self.unit,np.log(self.unit)+6,1-1e-7]+[20*self.unit]*p
        fixed=rho is not None
        free=[0,1]+list(range(3,3+p))   # indices optimised when rho is fixed
        starts=[]
        if (start is not None and np.isfinite(start.mean)
                and np.isfinite(start.total_sd) and start.total_sd > 0):
            warm_beta=(list(start.coefficients) if len(start.coefficients)==p else [0.]*p)
            starts.append([start.mean,np.log(start.total_sd),start.rho]+warm_beta)
        if not starts or multistart:
            starts.extend([[self.center,np.log(self.unit),r]+[0.]*p for r in ([rho] if fixed else [.15,.75])])
        candidates=[]
        before=self.evaluations
        for initial in starts:
            initial=np.clip(initial,lower,upper)
            if fixed:
                def fun(x):
                    v,g=self.objective(np.r_[x[:2],rho,x[2:]],order=order)
                    return v,g[free]
                result=minimize(fun,initial[free],jac=True,method='L-BFGS-B',
                                bounds=[(lower[k],upper[k]) for k in free],
                                options={'ftol':1e-10,'gtol':2e-6,'maxiter':180,'maxls':30})
                theta=np.r_[result.x[:2],rho,result.x[2:]]
            else:
                result=minimize(lambda x:self.objective(x,order=order),initial,jac=True,method='L-BFGS-B',
                                bounds=list(zip(lower,upper)),options={'ftol':1e-10,'gtol':2e-6,'maxiter':180,'maxls':30})
                theta=result.x
            if np.isfinite(result.fun):
                success=bool(result.success)
                if not success and 'ABNORMAL' in str(result.message) and np.max(np.abs(result.jac))<1e-3:
                    # The line search can fail on rounding noise when the start
                    # is already the maximum; a vanishing gradient is convergence.
                    success=True
                candidates.append((float(-result.fun),theta,success,str(result.message)))
        if not fixed:
            zero=self.fit(rho=0.,start=start,order=order,verify=False,multistart=False)
            candidates.append((zero.loglik,np.r_[zero.mean,np.log(zero.total_sd),0.,zero.coefficients],
                               zero.converged,zero.message))
        if not candidates:
            raise ArithmeticError('no finite likelihood fit')
        ll,theta,success,message=max(candidates,key=lambda v:v[0])
        error=0.
        if verify:
            checked=self.evaluate(theta[0],np.exp(theta[1]),theta[2],order=-1 if order<=0 else 2*order,
                                  coefficients=theta[3:],beta_gradient=False)[0]
            error=abs(checked-ll)
            if error>1e-6 and 0<order<384:
                warm=MICFit(theta[0],np.exp(theta[1]),theta[2],ll,success,0,order,error,message,
                            tuple(float(c) for c in theta[3:]))
                return self.fit(rho=rho,start=warm,order=2*order,verify=True,multistart=False,
                                box_maximum=box_maximum)
            if error>1e-5 and order>=384:
                warm=MICFit(theta[0],np.exp(theta[1]),theta[2],ll,success,0,order,error,message,
                            tuple(float(c) for c in theta[3:]))
                return self.fit(rho=rho,start=warm,order=0,verify=True,multistart=False,
                                box_maximum=box_maximum)
            success=success and error<=1e-5
            ll=checked
        edge=any(abs(theta[k]-lower[k])<1e-5 or abs(theta[k]-upper[k])<1e-5 for k in [0,1]+list(range(3,3+p)))
        if not fixed and theta[2]>1-1e-5:
            edge=True
        if edge and box_maximum:
            # Maximum over the bounded nuisance box; used only where both
            # likelihoods of a statistic are maximised over the same box.
            message='maximum on the nuisance parameter box; '+message
        elif edge:
            success=False; message='optimizer reached a numerical parameter boundary'
        return MICFit(float(theta[0]),float(np.exp(theta[1])),float(theta[2]),float(ll),bool(success),
                      self.evaluations-before,int(order),float(error),message,
                      tuple(float(c) for c in theta[3:]))
