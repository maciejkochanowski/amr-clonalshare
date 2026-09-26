"""Restricted-nuisance profile bootstrap; approximate inference."""
from collections import OrderedDict
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import ndtr,ndtri
from amr_clonalshare import _population_probit_numerics as numerical

_TABLE_PATH=Path(__file__).resolve().parents[1]/'histogram_table/histogram_table.py'
_SPEC=importlib.util.spec_from_file_location('histogram_table_api',_TABLE_PATH)
table_api=importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name]=table_api
_SPEC.loader.exec_module(table_api)
HISTOGRAM_SHA256='e4c2be7ac3a468c749ca61d3a73905238aed7c33b877c46c3819388de4a0a595'
NUMERICAL_SHA256='5bb56b972497aa951f0d70b8bf14f56e83a915a23a03f7c854336666f82594e8'
A_GRID=tuple([-8.]+[i/10 for i in range(-60,61)]+[8.])
ALT_RHOS=(0.,.01,.025,.05,.1,.2,.35,.5,.65,.8,.9,.95,.975,.99,.995,.999,1.)
RHO_GRID=tuple(sorted(set([i/40 for i in range(41)]+[.001,.005,.01,.99,.995,.999])))


def _rho(rho):
    rho=float(rho)
    if not math.isfinite(rho) or not 0<=rho<=1:raise ValueError('rho outside [0,1]')
    return rho


def canonical_node(rho):
    """Internal mesh convention: nearest decimal at twelve fractional places."""
    return float(round(_rho(rho),12))


def restricted_fit(counts,sizes,rho):
    """Continuous numerical nuisance fit with rho fixed; no repeat-count refusal.

    Exact boundary formulae, expanding finite bracket, bounded scalar solve and
    fixed local diagnostic probes. This is not a certified likelihood bound.
    """
    rho=_rho(rho);k,m=numerical.validate_counts(counts,sizes)
    n=int(m.sum());successes=int(k.sum());start=time.perf_counter()
    if successes in (0,n):
        return dict(status='constant_outcome',a=None,prevalence=successes/n,nll=0.,n_eval=0,seconds=time.perf_counter()-start)
    likelihood=numerical.ProfileLikelihood(k,m)
    if rho==0:
        p=successes/n;a=float(ndtri(p));value=likelihood.nll(a,0.)
        bracket=None
    elif rho==1:
        if np.any((k>0)&(k<m)):
            return dict(status='impossible_null',a=None,prevalence=None,nll=None,n_eval=0,seconds=time.perf_counter()-start)
        p=float(np.mean(k==m));a=float(ndtri(p));value=likelihood.nll(a,1.)
        bracket=None
    else:
        center=float(ndtri(successes/n));center=min(8.5,max(-8.5,center))
        center_value=likelihood.nll(center,rho);bracket=None
        if not math.isfinite(center_value):raise RuntimeError('No finite central nuisance value')
        for width in (1.,2.,4.,8.,18.):
            low=max(-9.,center-width);high=min(9.,center+width)
            lv=likelihood.nll(low,rho);hv=likelihood.nll(high,rho)
            if lv>=center_value and hv>=center_value:
                bracket=(low,high);break
        if bracket is None:raise RuntimeError('Nuisance minimum not bracketed inside [-9,9]')
        opt=minimize_scalar(lambda a:likelihood.nll(a,rho),bounds=bracket,method='bounded',
                            options=dict(xatol=2e-7,maxiter=150))
        if not opt.success or not math.isfinite(opt.fun) or abs(opt.x)>=8.999:
            raise RuntimeError('Restricted nuisance optimization failed')
        a=float(opt.x);value=float(opt.fun);p=float(ndtr(a))
        probes=[center,bracket[0],bracket[1],max(-9.,a-1e-4),min(9.,a+1e-4)]
        if min(likelihood.nll(x,rho) for x in probes)<value-1e-7:
            raise RuntimeError('Restricted nuisance diagnostic found a better value')
    if not math.isfinite(value):raise RuntimeError('Nonfinite restricted likelihood')
    return dict(status='ok',a=a,prevalence=p,nll=float(value),n_eval=likelihood.n_eval,
                bracket=list(bracket) if bracket else None,seconds=time.perf_counter()-start)


def inclusive_rank(observed,reference):
    observed=float(observed);reference=np.asarray(reference,dtype=float)
    if reference.ndim!=1 or not len(reference) or np.isnan(reference).any() or math.isnan(observed):
        raise ValueError('Invalid rank values')
    return (1+int(np.count_nonzero(reference>=observed)))/(len(reference)+1)


def node_components(nodes,accepted):
    nodes=np.asarray(nodes,dtype=float);accepted=np.asarray(accepted,dtype=bool)
    if nodes.ndim!=1 or accepted.shape!=nodes.shape or len(nodes)<2 or nodes[0]!=0 or nodes[-1]!=1 or np.any(np.diff(nodes)<=0):
        raise ValueError('Ordered complete rho grid required')
    edges=np.concatenate(([0.],(nodes[:-1]+nodes[1:])/2,[1.]))
    result=[]
    for i,keep in enumerate(accepted):
        if keep:
            lo,hi=float(edges[i]),float(edges[i+1])
            if result and result[-1][1]==lo:result[-1][1]=hi
            else:result.append([lo,hi])
    return result


def _rng(seed,case_key,rho):
    encoded=json.dumps([int(seed),str(case_key),float(rho).hex(),'restricted_bootstrap_v1'],separators=(',',':')).encode()
    words=np.frombuffer(hashlib.sha256(encoded).digest(),dtype='<u4').tolist()
    return np.random.default_rng(words)


def reference_counts(sizes,a,rho,*,B,seed,case_key):
    if isinstance(B,bool) or int(B)!=B or B<1:raise ValueError('B must be a positive integer')
    rho=_rho(rho);a=float(a)
    if not math.isfinite(a):raise ValueError('Finite generating intercept required')
    _,sizes=numerical.validate_counts(np.zeros(len(sizes)),sizes)
    rng=_rng(seed,case_key,rho)
    generated=np.empty((int(B),len(sizes)),dtype=np.int64)
    for b in range(int(B)):
        z=rng.normal(size=len(sizes))
        q=(a+z>=0).astype(float) if rho==1 else ndtr((a+math.sqrt(rho)*z)/math.sqrt(1-rho))
        generated[b]=rng.binomial(sizes,q)
    return generated


class BootstrapEngine:
    """One exact geometry, a fixed statistic grid and bounded null-table cache."""
    def __init__(self,sizes,*,a_grid=A_GRID,alternative_rhos=ALT_RHOS,memory_budget_bytes=256*2**20):
        if hashlib.sha256(_TABLE_PATH.read_bytes()).hexdigest()!=HISTOGRAM_SHA256:
            raise RuntimeError('Histogram API source identity changed')
        if hashlib.sha256(Path(numerical.__file__).read_bytes()).hexdigest()!=NUMERICAL_SHA256:
            raise RuntimeError('Numerical source identity changed')
        _,self.input_sizes=numerical.validate_counts(np.zeros(len(sizes)),sizes)
        self.order=np.argsort(self.input_sizes,kind='stable');self.sizes=self.input_sizes[self.order]
        self.a_grid=tuple(map(float,a_grid));self.alternative_rhos=tuple(map(_rho,alternative_rhos))
        self.budget=int(memory_budget_bytes)
        self.alternative=table_api.build_table(self.sizes,[(a,r) for r in self.alternative_rhos for a in self.a_grid],memory_budget_bytes=self.budget)
        self.null_tables=OrderedDict()

    def _null(self,rho):
        if rho in self.null_tables:
            self.null_tables.move_to_end(rho);return self.null_tables[rho]
        remaining=self.budget-self.alternative.nbytes
        required=self.alternative.logp.shape[0]*len(self.a_grid)*17
        while self.null_tables and sum(x.nbytes for x in self.null_tables.values())+required>remaining:
            self.null_tables.popitem(last=False)
        t=table_api.build_table(self.sizes,[(a,rho) for a in self.a_grid],memory_budget_bytes=remaining)
        self.null_tables[rho]=t
        return t

    def statistics(self,canonical_counts,rho):
        rho=_rho(rho)
        h=table_api.encode_histograms(canonical_counts,self.sizes,self.alternative)
        alt=table_api.score_histograms(h,self.alternative,backend='csr',batch_block=256,parameter_block=256).max(axis=1)
        null=table_api.score_histograms(h,self._null(rho),backend='csr',batch_block=256,parameter_block=256).max(axis=1)
        if np.any(~np.isfinite(alt)):raise RuntimeError('No finite alternative statistic value')
        result=np.maximum(0.,2*(alt-null))
        counts=np.asarray(canonical_counts)
        if counts.ndim==1:counts=counts[None,:]
        totals=counts.sum(axis=1)
        # Both unrestricted and restricted maxima are one at p=0/1.
        result[(totals==0)|(totals==self.sizes.sum())]=0.
        return result

    def evaluate_node(self,counts,rho,*,B,seed,case_key):
        if isinstance(B,bool) or int(B)!=B or B<1:raise ValueError('B must be a positive integer')
        B=int(B);rho=_rho(rho);k,_=numerical.validate_counts(counts,self.input_sizes);k=k[self.order]
        fit=restricted_fit(k,self.sizes,rho)
        if fit['status']=='constant_outcome':
            return dict(rho=rho,status=fit['status'],p_value=1.,B=B,exceedances=B,fit=fit)
        if fit['status']=='impossible_null':
            return dict(rho=rho,status=fit['status'],p_value=1/(B+1),B=B,exceedances=0,fit=fit,
                        observed_statistic=None,statistic_is_positive_infinity=True)
        # Row interleaving keeps every reference-count prefix invariant in B.
        generated=reference_counts(self.sizes,fit['a'],rho,B=B,seed=seed,case_key=case_key)
        observed=float(self.statistics(k,rho)[0]);reference=self.statistics(generated,rho)
        exceedances=int(np.count_nonzero(reference>=observed))
        return dict(rho=rho,status='ok',p_value=inclusive_rank(observed,reference),B=B,exceedances=exceedances,
                    fit=fit,observed_statistic=observed,reference_counts_sha256=hashlib.sha256(generated.astype('<i8').tobytes()).hexdigest(),
                    reference_statistics_sha256=hashlib.sha256(reference.astype('<f8').tobytes()).hexdigest())


def infer_profile_bootstrap(counts,sizes,*,B=999,seed=2026091202,case_key='query',rho_grid=RHO_GRID,
                            alphas=(.05,.04),refinement_levels=3,engine=None):
    if isinstance(B,bool) or int(B)!=B or B<1:raise ValueError('B must be a positive integer')
    alphas=tuple(float(a) for a in alphas)
    if not alphas or any(not 0<a<1 for a in alphas):raise ValueError('Invalid alpha')
    if isinstance(refinement_levels,bool) or int(refinement_levels)!=refinement_levels or refinement_levels<0:
        raise ValueError('Invalid refinement depth')
    nodes=sorted(set(map(canonical_node,rho_grid)))
    node_components(nodes,[True]*len(nodes))
    if np.asarray(counts).shape==(0,) and np.asarray(sizes).shape==(0,):
        return dict(status='unavailable',usable=False,sets=[],numerical_errors=[])
    k,m=numerical.validate_counts(counts,sizes)
    if engine is not None and not np.array_equal(engine.input_sizes,m):
        raise ValueError('Supplied engine must use the same ordered paired sizes')
    if np.all(m==1) or k.sum() in (0,m.sum()):
        return dict(status='structurally_unidentified' if np.all(m==1) else 'constant_outcome',
                    components=[[0.,1.]],ci_low=0.,ci_high=1.,numerical_errors=[],usable=True,
                    inference='approximate_restricted_nuisance_bootstrap',
                    sets=[dict(alpha=a,components=[[0.,1.]],ci_low=0.,ci_high=1.,empty=False) for a in alphas])
    engine=engine or BootstrapEngine(m)
    records={};errors=[]
    def evaluate(r):
        try:records[r]=engine.evaluate_node(k,r,B=B,seed=seed,case_key=case_key)
        except (RuntimeError,FloatingPointError,MemoryError) as exc:
            errors.append(dict(rho=r,error=type(exc).__name__+': '+str(exc)))
            records[r]=dict(rho=r,status='unresolved',p_value=None)
    for r in nodes:evaluate(r)
    for level in range(refinement_levels):
        extra=[]
        for lo,hi in zip(nodes[:-1],nodes[1:]):
            if hi-lo<=.003125+1e-15:continue
            p=records[lo]['p_value'];q=records[hi]['p_value']
            need=p is None or q is None or any((p>a)!=(q>a) for a in alphas) or (.02<=p<=.08) or (.02<=q<=.08)
            if need:extra.append(canonical_node((lo+hi)/2))
        for r in extra:evaluate(r)
        nodes=sorted(records)
    sets=[]
    for alpha in alphas:
        components=node_components(nodes,[records[r]['p_value'] is None or records[r]['p_value']>alpha for r in nodes])
        sets.append(dict(alpha=alpha,components=components,ci_low=components[0][0] if components else None,
                         ci_high=components[-1][1] if components else None,empty=not components))
    return dict(status='numerical_incomplete' if errors else 'ok',usable=not errors,numerical_errors=errors,sets=sets,
                nodes=[records[r] for r in nodes],maximum_grid_gap=max(np.diff(nodes)),B=B,p_value_resolution=1/(B+1),
                inference='approximate_restricted_nuisance_bootstrap',n_repeated_groups=int(np.count_nonzero(m>=2)),
                memory_budget_scope='Retained likelihood-table arrays only; excludes score matrices, references and integration temporaries',
                internal_mesh_identity='rho rounded to twelve decimal places before evaluation and stream assignment',
                alternative_table_identity=engine.alternative.identity())
