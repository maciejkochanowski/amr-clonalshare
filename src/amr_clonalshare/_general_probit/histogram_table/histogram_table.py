"""Histogram table of grouped-count probabilities; it computes no interval by itself.

A finite grid defines a different statistic from adaptive continuous profiling.
The quadrature of the profile likelihood supplies the table entries. Impossible bins are masked
explicitly so zero histogram entries never multiply negative infinity.
"""
from dataclasses import dataclass
from collections import Counter
import hashlib
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.special import logsumexp

from amr_clonalshare import _population_probit_numerics as numerics

NUMERICAL_SHA256='5bb56b972497aa951f0d70b8bf14f56e83a915a23a03f7c854336666f82594e8'


@dataclass(frozen=True)
class HistogramTable:
    geometry: tuple
    theta: tuple
    offsets: dict
    logp: np.ndarray
    finite_values: np.ndarray
    impossible: np.ndarray
    max_normalization_error: float
    numerical_source_sha256: str

    @property
    def nbytes(self):
        return self.logp.nbytes+self.finite_values.nbytes+self.impossible.nbytes

    def identity(self):
        return dict(geometry=list(self.geometry),theta=[[float(a).hex(),float(r).hex()] for a,r in self.theta],
            logp_shape=list(self.logp.shape),logp_sha256=hashlib.sha256(self.logp.tobytes(order='C')).hexdigest(),
            numerical_source_sha256=self.numerical_source_sha256,normalization='profile_quadrature',
            max_normalization_error=self.max_normalization_error,dtype=str(self.logp.dtype),nbytes=self.nbytes)


def build_table(sizes,theta,*,memory_budget_bytes=256*2**20,count_block=128):
    _,m=numerics.validate_counts(np.zeros(len(sizes)),sizes)
    geometry=tuple(sorted(m.tolist()))
    points=np.asarray(theta,dtype=float)
    if points.ndim!=2 or points.shape[1]!=2 or len(points)==0 or not np.isfinite(points).all():
        raise ValueError('Require nonempty finite (a,rho) parameter pairs')
    if np.any((points[:,1]<0)|(points[:,1]>1)):
        raise ValueError('rho must lie in zero to one')
    if count_block<1:raise ValueError('Positive count block required')
    source=hashlib.sha256(Path(numerics.__file__).read_bytes()).hexdigest()
    if source!=NUMERICAL_SHA256:raise RuntimeError('the numerical source differs from the recorded one')
    offsets={};bins=0
    for size in sorted(set(geometry)):
        offsets[size]=bins;bins+=size+1
    # Three retained arrays: float64 raw, float64 finite, bool impossible.
    if bins*len(points)*17>memory_budget_bytes:
        raise MemoryError('Table exceeds the declared retained-array memory budget')
    logp=np.empty((bins,len(points)),dtype=np.float64)
    normalization=[]
    for size,offset in offsets.items():
        for j,(a,rho) in enumerate(points):
            for begin in range(0,size+1,count_block):
                end=min(size+1,begin+count_block)
                logp[offset+begin:offset+end,j]=numerics.selected_count_log_probabilities(
                    size,float(a),float(rho),tuple(range(begin,end)))
        normalization.extend(np.abs(np.expm1(logsumexp(logp[offset:offset+size+1],axis=0))).tolist())
    if np.isnan(logp).any() or np.isposinf(logp).any():raise RuntimeError('Invalid table values')
    impossible=np.isneginf(logp);finite=np.where(impossible,0.,logp)
    for array in (logp,finite,impossible):array.flags.writeable=False
    return HistogramTable(geometry,tuple(map(tuple,points)),offsets,logp,finite,impossible,max(normalization),source)


def encode_histograms(counts,sizes,table):
    _,m=numerics.validate_counts(np.zeros(len(sizes)),sizes)
    if tuple(sorted(m.tolist()))!=table.geometry:raise ValueError('Wrong exact geometry')
    k=np.asarray(counts,dtype=float)
    if k.ndim==1:k=k[None,:]
    if k.ndim!=2 or k.shape[1]!=len(m) or len(k)==0 or not np.isfinite(k).all():
        raise ValueError('Require nonempty dataset-by-group count matrix')
    if np.any(k!=np.floor(k)) or np.any(k<0) or np.any(k>m[None,:]):
        raise ValueError('Counts must be integer and between zero and group size')
    offsets=np.array([table.offsets[int(size)] for size in m],dtype=np.int64)
    columns=(k.astype(np.int64)+offsets[None,:]).ravel()
    rows=np.repeat(np.arange(len(k)),len(m))
    h=sparse.csr_matrix((np.ones(len(columns)),(rows,columns)),shape=(len(k),table.logp.shape[0]))
    h.sum_duplicates();h.sort_indices()
    return h


def score_histograms(histograms,table,*,backend='csr',batch_block=256,parameter_block=256):
    if backend not in ('csr','dense') or batch_block<1 or parameter_block<1:
        raise ValueError('Invalid scoring backend or block sizes')
    h=sparse.csr_matrix(histograms,dtype=np.float64,copy=True)
    if h.shape[1]!=table.logp.shape[0] or not np.isfinite(h.data).all() or np.any(h.data<0):
        raise ValueError('Invalid histogram matrix')
    if np.any(h.data!=np.floor(h.data)):
        raise ValueError('Histogram entries must be nonnegative integer frequencies')
    h.sum_duplicates();h.sort_indices();h.eliminate_zeros()
    if np.any(np.asarray(h.sum(axis=1)).ravel()!=len(table.geometry)):
        raise ValueError('Histogram rows must encode this complete geometry')
    for size,multiplicity in Counter(table.geometry).items():
        offset=table.offsets[size]
        if np.any(np.asarray(h[:,offset:offset+size+1].sum(axis=1)).ravel()!=multiplicity):
            raise ValueError('Histogram size-block multiplicity disagrees with the geometry')
    output=np.empty((h.shape[0],table.logp.shape[1]),dtype=np.float64)
    for begin in range(0,h.shape[0],batch_block):
        end=min(h.shape[0],begin+batch_block);part=h[begin:end]
        if backend=='dense':part=part.toarray()
        for left in range(0,table.logp.shape[1],parameter_block):
            right=min(table.logp.shape[1],left+parameter_block)
            values=np.asarray(part@table.finite_values[:,left:right])
            bad=np.asarray(part@table.impossible[:,left:right])>0
            values[bad]=-np.inf;output[begin:end,left:right]=values
    return output


def grid_lr_scores(loglikelihoods,theta,null_rho):
    """Finite-grid LR score; not a confidence set by itself.

    No equivalence to the adaptive continuous likelihood optimum is assumed.
    Calibration, nuisance treatment and a continuum contract are absent here.
    """
    values=np.asarray(loglikelihoods,dtype=float);points=np.asarray(theta,dtype=float)
    if values.ndim!=2 or points.shape!=(values.shape[1],2) or np.isnan(values).any() or np.isposinf(values).any():
        raise ValueError('Invalid grid likelihood matrix')
    null=points[:,1]==null_rho
    if not null.any():raise ValueError('Null rho absent from the finite grid')
    maximum=values.max(axis=1)
    if not np.isfinite(maximum).all():raise ValueError('No finite alternative grid point')
    return np.maximum(0.,2*(maximum-values[:,null].max(axis=1)))
