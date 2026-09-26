"""Independent confirmation of MIC population-rho intervals; one atomic chunk."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy.stats import chi2

from amr_clonalshare.censored import censored_clonal_share
from amr_clonalshare.mic_inference import exact_gaussian_interval
from amr_clonalshare._mic_likelihood import GaussianMICLikelihood
from amr_clonalshare._mic_panels import observe_panel
from benchmarks.mic_inference.profile_bootstrap import calibrate_profile, profile_bounds
from benchmarks.empirical.mic_calibration import design_grid

ROOT_SEED=20260915071
#: Seed stream of the calibration; stream 1 is reserved for unit-test datasets.
CALIBRATION_STREAM=2


def generate(cell,replicate,stream=CALIBRATION_STREAM):
    rng=np.random.default_rng(np.random.SeedSequence([ROOT_SEED,stream,cell['cell'],replicate]))
    g=cell['groups']; pattern=cell['size_pattern']; rho=cell['rho']
    if pattern=='balanced': sizes=np.full(g,20)
    elif pattern=='uneven': sizes=np.maximum(2,np.rint(np.exp(rng.normal(2.2,1.,g))).astype(int)).clip(2,100)
    elif pattern=='singleton_rich': sizes=np.r_[np.ones(g//2,dtype=int),np.full(g-g//2,20)]
    else: sizes=np.ones(g,dtype=int)
    code=np.repeat(np.arange(g),sizes)
    effects=rng.normal(0,np.sqrt(rho),g)
    if cell['residual']=='t4': residual=rng.standard_t(4,len(code))/np.sqrt(2.)
    elif cell['residual']=='contaminated_normal':
        scale=np.where(rng.random(len(code))<.05,5.,1.)
        residual=rng.normal(size=len(code))*scale/np.sqrt(2.2)
    else: residual=rng.normal(size=len(code))
    y=effects[code]+np.sqrt(1-rho)*residual
    edges=np.arange(-4.,5.) if cell['reading']=='wide_panel' else np.arange(-1.,2.)
    if cell['reading']=='exact': a,b=y.copy(),y.copy()
    else: a,b=observe_panel(y,edges)
    inner_seed=int(np.random.SeedSequence([ROOT_SEED,stream,cell['cell'],replicate,991]).generate_state(1)[0])
    return a,b,code,edges,inner_seed


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cell',type=int,required=True)
    ap.add_argument('--start',type=int,default=0)
    ap.add_argument('--replicates',type=int,default=20)
    ap.add_argument('--bootstrap',type=int,default=199)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--bounds-every',type=int,default=1)
    args=ap.parse_args()
    if args.start<0 or args.replicates<1 or args.bootstrap<19 or args.bounds_every<1:
        ap.error('invalid simulation arguments')
    cells=design_grid()
    if not 0<=args.cell<len(cells): ap.error('unknown cell')
    cell=cells[args.cell]; args.output.parent.mkdir(parents=True,exist_ok=True)
    if args.output.exists(): raise FileExistsError(args.output)
    rows=[]; stored_statistics=[]; started=time.perf_counter()
    for rep in range(args.start,args.start+args.replicates):
        t=time.perf_counter(); a,b,c,edges,seed=generate(cell,rep)
        row=dict(cell=args.cell,replicate=rep,rho=cell['rho'],n=len(c),groups=cell['groups'],
                 reading=cell['reading'],size_pattern=cell['size_pattern'],residual=cell['residual'],
                 domain=cell['domain'],end_fraction=float(np.mean(~np.isfinite(a)|~np.isfinite(b))),
                 moment_low=np.nan,moment_high=np.nan,moment_reportable=False,moment_covered=False,
                 moment_estimate=np.nan,interval_low=0.,interval_high=1.,interval_reportable=False,interval_covered=True,
                 interval_test_covered=True,ml_estimate=np.nan,chi_square_covered=False,lr_truth=np.nan,
                 critical=np.inf,bootstrap_attempted=0,bootstrap_failed=0,profile_failures=0,
                 raw_set_empty=False,failure='',seconds=0.,quadrature_error=np.nan)
        try:
            moment=censored_clonal_share(a,b,c)
            row.update(moment_low=moment.ci_low,moment_high=moment.ci_high,moment_reportable=moment.estimable,
                       moment_covered=bool(moment.estimable and moment.ci_low<=cell['rho']<=moment.ci_high),
                       moment_estimate=moment.kappa)
        except Exception as exc:
            row['moment_failure']=f'{type(exc).__name__}: {exc}'
        try:
            if cell['reading']=='exact':
                result=exact_gaussian_interval(a,c)
                row.update(interval_low=result.low,interval_high=result.high,interval_reportable=True,
                           interval_covered=result.low<=cell['rho']<=result.high,
                           interval_test_covered=result.low<=cell['rho']<=result.high,
                           raw_set_empty=result.raw_set_empty,critical=np.nan)
                stored_statistics.append(np.full(args.bootstrap,np.nan))
            else:
                problem=GaussianMICLikelihood(a,b,c)
                fit,critical,statistics=calibrate_profile(problem,edges,a==b,n_boot=args.bootstrap,seed=seed)
                stored_statistics.append(statistics)
                row.update(ml_estimate=fit.rho,critical=critical,bootstrap_attempted=args.bootstrap,
                           bootstrap_failed=int(np.isinf(statistics).sum()),quadrature_error=fit.quadrature_error)
                if not fit.converged:
                    raise ArithmeticError('unrestricted fit did not converge: '+fit.message)
                truth=problem.fit(rho=cell['rho'],start=fit,multistart=False)
                if not truth.converged or truth.loglik>fit.loglik+1e-4:
                    raise ArithmeticError('profile at truth unresolved')
                lr=max(0.,2*(fit.loglik-truth.loglik))
                row.update(lr_truth=lr,chi_square_covered=lr<=chi2.ppf(.95,1),interval_test_covered=lr<=critical)
                if rep%args.bounds_every==0:
                    low,high,nfailed=profile_bounds(problem,fit,critical)
                    row.update(interval_low=low,interval_high=high,profile_failures=nfailed,
                               interval_covered=low<=cell['rho']<=high,
                               interval_reportable=bool(np.isfinite(critical) and nfailed==0))
                else:
                    row.update(interval_low=np.nan,interval_high=np.nan,interval_covered=lr<=critical,
                               interval_reportable=bool(np.isfinite(critical)))
                if not np.isfinite(critical): row['failure']='bootstrap critical value unresolved'
                if row['profile_failures']: row['failure']='profile endpoints unresolved; full range returned'
        except Exception as exc:
            row['failure']=f'{type(exc).__name__}: {exc}'
            if len(stored_statistics)<len(rows)+1:
                stored_statistics.append(np.full(args.bootstrap,np.inf))
        row['seconds']=time.perf_counter()-t
        rows.append(row)
    keys=sorted(set().union(*(r.keys() for r in rows)))
    temporary=args.output.with_suffix('.partial')
    with temporary.open('w',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=keys); writer.writeheader(); writer.writerows(rows)
        out.flush();os.fsync(out.fileno())
    temporary.rename(args.output)
    np.savez_compressed(args.output.with_suffix('.bootstrap.npz'),statistics=np.stack(stored_statistics))
    sources={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
             [Path(__file__),*Path('src/amr_clonalshare').glob('*.py')]}
    receipt=dict(arguments=vars(args)|{'output':str(args.output)},cell=cell,source_hashes=sources,
                 root_seed=ROOT_SEED,elapsed_seconds=time.perf_counter()-started,
                 python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                 csv_sha256=hashlib.sha256(args.output.read_bytes()).hexdigest(),
                 slurm_job_id=os.getenv('SLURM_JOB_ID'),slurm_array_task_id=os.getenv('SLURM_ARRAY_TASK_ID'))
    args.output.with_suffix('.receipt.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps({'cell':args.cell,'replicates':len(rows),'seconds':receipt['elapsed_seconds'],
                      'reportable':sum(r['interval_reportable'] for r in rows),
                      'covered':sum(r['interval_covered'] for r in rows)}))


if __name__=='__main__':main()
