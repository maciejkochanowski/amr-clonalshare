"""Prespecified MIC stress grid; one independent output per cell and implementation.

Usage: python mic_calibration.py --cell 0 --replicates 2000 --output out/cell.csv
Model-based population rho is the target. Failure/refusal remain in denominators.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from amr_clonalshare.censored import censored_clonal_share


def design_grid():
    cells=[]
    for groups in (8,30):
        for size_pattern in ('balanced','uneven','singleton_rich'):
            for rho in (0.,.1,.5,.9):
                for reading in ('exact','wide_panel','narrow_panel'):
                    cells.append(dict(groups=groups,size_pattern=size_pattern,rho=rho,
                                      reading=reading,residual='normal',domain='core'))
    for residual in ('t4','contaminated_normal'):
        for rho in (.1,.5,.9):
            for reading in ('exact','wide_panel','narrow_panel'):
                cells.append(dict(groups=30,size_pattern='uneven',rho=rho,reading=reading,
                                  residual=residual,domain='distribution_stress'))
    for rho in (0.,.1,.5,.9):
        for reading in ('exact','wide_panel','narrow_panel'):
            cells.append(dict(groups=30,size_pattern='all_singletons',rho=rho,reading=reading,
                              residual='normal',domain='nonidentifiable_control'))
    for i,c in enumerate(cells): c['cell']=i
    return cells


def generate(cell,replicate):
    rng=np.random.default_rng(np.random.SeedSequence([20260915,cell['cell'],replicate]))
    g=cell['groups']; pattern=cell['size_pattern']; rho=cell['rho']
    if pattern=='balanced': sizes=np.full(g,20)
    elif pattern=='uneven': sizes=np.maximum(2,np.rint(np.exp(rng.normal(2.2,1.,g))).astype(int)).clip(2,100)
    elif pattern=='singleton_rich': sizes=np.r_[np.ones(g//2,dtype=int),np.full(g-g//2,20)]
    else: sizes=np.ones(g,dtype=int)
    labels=np.repeat(np.arange(g),sizes)
    effects=rng.normal(0.,math.sqrt(rho),g)
    if cell['residual']=='t4': residual=rng.standard_t(4,len(labels))/math.sqrt(2.)
    elif cell['residual']=='contaminated_normal':
        scales=np.where(rng.random(len(labels))<.05,5.,1.)
        residual=rng.normal(size=len(labels))*scales/math.sqrt(.95+.05*25.)
    else: residual=rng.normal(size=len(labels))
    y=effects[labels]+math.sqrt(1.-rho)*residual
    if cell['reading']=='exact': lo,hi=y.copy(),y.copy()
    else:
        edge=4. if cell['reading']=='wide_panel' else 1.
        lo=np.floor(y); hi=lo+1.
        lo[y<=-edge]=-np.inf; hi[y<=-edge]=-edge
        lo[y>edge]=edge; hi[y>edge]=np.inf
    return lo,hi,labels


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cell',type=int,required=True)
    ap.add_argument('--replicates',type=int,default=2000)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--implementation',default='candidate')
    args=ap.parse_args()
    if args.replicates<1: ap.error('invalid replicate budget')
    cells=design_grid()
    if not 0<=args.cell<len(cells): ap.error('unknown cell')
    cell=cells[args.cell]; args.output.parent.mkdir(parents=True,exist_ok=True)
    if args.output.exists(): raise FileExistsError(args.output)
    fields=['implementation','cell','replicate','rho','n','groups','end_fraction','kappa','low','high',
            'estimable','covered','width','absolute_error',
            'fit_seconds','fit_converged','fit_iterations','fit_max_scaled_change','failure']
    started=time.perf_counter(); temporary=args.output.with_suffix('.partial')
    with temporary.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields); writer.writeheader()
        for rep in range(args.replicates):
            lo,hi,labels=generate(cell,rep)
            row=dict(implementation=args.implementation,cell=args.cell,replicate=rep,rho=cell['rho'],
                     n=len(labels),groups=cell['groups'],end_fraction=float(np.mean(~np.isfinite(lo)|~np.isfinite(hi))),
                     covered=False,estimable=False,width=1.,failure='')
            t=time.perf_counter()
            try:
                r=censored_clonal_share(lo,hi,labels)
                row.update(kappa=r.kappa,low=r.ci_low,high=r.ci_high,estimable=bool(r.estimable),
                           covered=bool(r.estimable and r.ci_low<=cell['rho']<=r.ci_high),
                           width=float(r.ci_high-r.ci_low) if r.estimable and np.isfinite(r.ci_high-r.ci_low) else 1.,
                           absolute_error=abs(r.kappa-cell['rho']) if r.estimable else math.nan,
                           failure='' if r.estimable else r.reason)
            except Exception as exc:
                row['failure']=f'{type(exc).__name__}: {exc}'
            row['fit_seconds']=time.perf_counter()-t
            if 'r' in locals():
                row.update(fit_converged=getattr(r,'fit_converged',None), fit_iterations=getattr(r,'fit_iterations',None), fit_max_scaled_change=getattr(r,'fit_max_scaled_change',None))
                del r
            writer.writerow(row)
        handle.flush(); os.fsync(handle.fileno())
    temporary.rename(args.output)
    import amr_clonalshare.censored as module
    receipt=dict(cell=cell,replicates=args.replicates,
                 implementation=args.implementation,elapsed_seconds=time.perf_counter()-started,
                 python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                 module_path=module.__file__,source_sha256=hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest(),
                 data_sha256=hashlib.sha256(args.output.read_bytes()).hexdigest(),
                 slurm_job_id=os.getenv('SLURM_JOB_ID'),slurm_array_task_id=os.getenv('SLURM_ARRAY_TASK_ID'))
    args.output.with_suffix('.receipt.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt,indent=2))


if __name__=='__main__': main()
