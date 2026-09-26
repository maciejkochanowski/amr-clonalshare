"""Reanalyse one agent with the null-wise MIC interval and preserve all tests."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import time
import numpy as np
from amr_clonalshare.censored import censored_clonal_share
from amr_clonalshare.mic_null_bootstrap import null_bootstrap_interval
from amr_clonalshare.mic_inference_cli import _strict

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--manifest',type=Path,required=True)
p.add_argument('--agent-index',type=int,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--workers',type=int,default=1)
p.add_argument('--bootstrap',type=int,default=199,help='simulated datasets per tested value of rho')
p.add_argument('--covariate-column',help='column of the input table whose levels enter the model as fixed effects')
a=p.parse_args()
records=json.loads(a.manifest.read_text())['records']
record=records[a.agent_index]
source=a.manifest.parent/record['file']
assert hashlib.sha256(source.read_bytes()).hexdigest()==record['sha256']
with source.open(newline='') as f: rows=list(csv.DictReader(f))
lo=np.array([float(r['lo']) for r in rows]);hi=np.array([float(r['hi']) for r in rows])
labels=np.array([r['lineage'] for r in rows]);panel=np.array([r['panel'] for r in rows])
covariate=np.array([r[a.covariate_column] for r in rows]) if a.covariate_column else None
a.output.mkdir(parents=True,exist_ok=True)
output=a.output/(record['agent']+'.json')
if output.exists(): raise FileExistsError(output)
started=time.perf_counter()
moment=censored_clonal_share(lo,hi,labels)
interval=null_bootstrap_interval(lo,hi,labels,panel_edges=record['panel_edges'],panel=panel,covariate=covariate,n_boot=a.bootstrap,seed=record['seed'],workers=a.workers)
import amr_clonalshare.mic_null_bootstrap as module
package=Path(module.__file__).parent
sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in package.glob('*.py')}
result=dict(observation=record,moment=moment.as_dict(),interval=interval.as_dict(),
    elapsed_seconds=time.perf_counter()-started,source_sha256=sources,
    slurm_job_id=os.getenv('SLURM_JOB_ID'),slurm_array_task_id=os.getenv('SLURM_ARRAY_TASK_ID'),
    bootstrap=a.bootstrap,covariate_column=a.covariate_column,seed=record['seed'],driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
temporary=output.with_suffix('.partial');temporary.write_text(json.dumps(_strict(result),indent=2,allow_nan=False));temporary.replace(output)
print(json.dumps(dict(agent=record['agent'],estimate=interval.fit.rho,low=interval.low,high=interval.high,
    reportable=interval.reportable,nulls=len(interval.tests),seconds=result['elapsed_seconds'])))
