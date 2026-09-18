"""Prepare the retained S. suis example under an explicit inferred-panel assumption.

Observed extrema are not independent evidence of tested laboratory endpoints.
This program keeps the inferred-panel working assumption solely for a matched-method
illustration. It does not label the resulting intervals as verified panels. The
panel is inferred within each testing laboratory, as the source study names
them, so that the lowest well of one laboratory is not read as an interior
dilution because the other laboratory tested lower. The country of isolation
is written beside every reading so that ``run_empirical.py`` can enter it as
a fixed effect (``--covariate-column country``); it distinguishes the two
laboratories as well as the two countries tested in one of them.
"""
import argparse
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from amr_clonalshare.censored import intervals_by_panel
from amr_clonalshare._mic_panels import _validated_problem

ROOT=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--country',help='keep only isolates from this country (one laboratory in the S. suis source); default: all isolates')
parser.add_argument('--output',type=Path,default=ROOT/'benchmarks/results_mic_release/mic_empirical_inputs')
parser.add_argument('--seed-base',type=int,default=20260916100)
args=parser.parse_args()
OUT=args.output
OUT.mkdir(exist_ok=True,parents=True)
DATA=ROOT/'examples/ssuis/data'
metadata=pd.read_csv(DATA/'metadata.csv',dtype=str,keep_default_na=False).set_index('genome_id')
if args.country:
    metadata=metadata[metadata.isolation_country.eq(args.country)]
    if metadata.empty: raise SystemExit('no isolates from '+args.country)
measurements=pd.read_csv(DATA/'mic_panel.csv',dtype=str,keep_default_na=False)
measurements=measurements[measurements.genome_id.isin(metadata.index)]
assert not measurements.duplicated(['genome_id','antibiotic']).any()
records=[]
for index,(agent,frame) in enumerate(measurements.groupby('antibiotic',sort=True)):
    frame=frame.set_index('genome_id').reindex(metadata.index)
    if frame.measurement.isna().any():
        raise ValueError('Incomplete aligned measurement frame: '+agent)
    panel=frame.testing_laboratory.to_numpy()
    lo,hi=intervals_by_panel(frame.measurement.astype(float),panel,operators=frame.measurement_sign)
    edges={}
    for laboratory in np.unique(panel):
        here=panel==laboratory
        finite=np.unique(np.r_[lo[here][np.isfinite(lo[here])],hi[here][np.isfinite(hi[here])]])
        edges[str(laboratory)]=(np.arange(finite.min(),finite.max()+.5,1.) if np.allclose(finite,np.round(finite)) else finite).tolist()
    _validated_problem(lo,hi,metadata.baps_cluster,edges,panel)
    output=OUT/(agent+'.csv')
    pd.DataFrame(dict(lo=lo,hi=hi,lineage=metadata.baps_cluster.to_numpy(),panel=panel,
                      country=metadata.isolation_country.to_numpy(),
                      isolate_id=metadata.index)).to_csv(output,index=False)
    records.append(dict(index=index,agent=agent,file=output.name,panel_edges=edges,
        n=len(lo),groups=metadata.baps_cluster.nunique(),
        panel_provenance='assumed_from_recorded_lattice_and_extrema_within_testing_laboratory; not_verified_laboratory_panel',
        endpoint_assumption='rc1_inferred_end_wells; same_intervals_for_both_methods',
        seed=args.seed_base+index,sha256=hashlib.sha256(output.read_bytes()).hexdigest()))
(OUT/'manifest.json').write_text(json.dumps(dict(records=records,subset=dict(country=args.country,isolates=len(metadata)),
    scientific_scope='Method comparison under an explicit observation assumption; not clinical inference',
    source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in [DATA/'metadata.csv',DATA/'mic_panel.csv',Path(__file__)]}),indent=2))
print('Prepared',len(records),'agents;',sum(x['n'] for x in records),'measurements')
