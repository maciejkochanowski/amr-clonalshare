"""Verify receipts and summarize coverage without hiding refusals."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm


def wilson(success,total,alpha=.05):
    if total==0: return (float('nan'),float('nan'))
    z=norm.ppf(1-alpha/2); p=success/total
    center=(p+z*z/(2*total))/(1+z*z/total)
    half=z*np.sqrt(p*(1-p)/total+z*z/(4*total*total))/(1+z*z/total)
    return float(center-half),float(center+half)


def summarize(directory,output):
    paths=sorted(Path(directory).glob('cell_*.csv'))
    if not paths: raise ValueError('no complete simulation chunks found')
    tables=[]; hashes={}; receipts=[]
    for path in paths:
        r=json.loads(path.with_suffix('.receipt.json').read_text())
        if hashlib.sha256(path.read_bytes()).hexdigest()!=r['csv_sha256']:
            raise ValueError(f'CSV checksum mismatch: {path}')
        frame=pd.read_csv(path)
        if len(frame)!=r['arguments']['replicates']:
            raise ValueError(f'replication count mismatch: {path}')
        expected=np.arange(r['arguments']['start'],r['arguments']['start']+len(frame))
        if not np.array_equal(frame.replicate,expected): raise ValueError(f'wrong replicate IDs: {path}')
        for name,digest in r['source_hashes'].items():
            key=name[name.index('benchmarks/'):] if 'benchmarks/' in name else name
            if key in hashes and hashes[key]!=digest: raise ValueError(f'mixed sources: {key}')
            hashes[key]=digest
        stats=np.load(path.with_suffix('.bootstrap.npz'))['statistics']
        if stats.shape!=(len(frame),r['arguments']['bootstrap']):
            raise ValueError(f'bootstrap matrix shape mismatch: {path}')
        tables.append(frame); receipts.append(r)
    data=pd.concat(tables,ignore_index=True)
    if data.duplicated(['cell','replicate']).any(): raise ValueError('duplicate datasets')
    rows=[]
    for cell,group in data.groupby('cell'):
        row={k:group.iloc[0][k] for k in ['cell','rho','groups','reading','size_pattern','residual','domain']}
        row['attempts']=len(group)
        for side in ['moment','interval']:
            ok=group[side+'_reportable'].astype(bool)
            covered=group[side+'_covered'].astype(bool)
            count=int(ok.sum()); success=int((ok&covered).sum())
            low,high=wilson(success,count)
            row.update({side+'_reportable':count,side+'_refused':len(group)-count,
                        side+'_covered':success,side+'_conditional_coverage':success/count if count else np.nan,
                        side+'_coverage_wilson_low':low,side+'_coverage_wilson_high':high,
                        side+'_reportability':count/len(group),
                        side+'_median_width':float((group.loc[ok,side+'_high']-group.loc[ok,side+'_low']).median()) if count else np.nan})
        both=group.moment_reportable&group.interval_reportable
        diff=group.loc[both,'interval_covered'].astype(float)-group.loc[both,'moment_covered'].astype(float)
        row['paired_n']=int(both.sum()); row['paired_coverage_difference']=float(diff.mean()) if len(diff) else np.nan
        row['paired_difference_mcse']=float(diff.std(ddof=1)/np.sqrt(len(diff))) if len(diff)>1 else np.nan
        row['interval_full_range_returns']=int((~group.interval_reportable).sum())
        row['interval_coverage_including_full_ranges']=float(group.interval_covered.mean())
        row['chi_square_coverage']=float(group.loc[group.interval_reportable,'chi_square_covered'].mean()) if row['reading']!='exact' else np.nan
        row['bootstrap_failed']=int(group.bootstrap_failed.sum())
        row['bootstrap_slots']=int(group.bootstrap_attempted.sum())
        row['profile_failures']=int(group.profile_failures.sum())
        row['raw_empty_sets']=int(group.raw_set_empty.sum())
        row['max_quadrature_error']=float(group.quadrature_error.max())
        row['median_seconds']=float(group.seconds.median())
        row['coverage_membership_mismatches']=int(((group.interval_covered!=group.interval_test_covered)&group.interval_reportable).sum())
        rows.append(row)
    summary=pd.DataFrame(rows)
    output=Path(output); output.parent.mkdir(parents=True,exist_ok=True)
    summary.to_csv(output,index=False)
    counts=dict(chunks=len(paths),datasets=len(data),cells=int(data.cell.nunique()),
                sum_fit_seconds=float(data.seconds.sum()),
                source_hashes=hashes,failures=data.failure.fillna('').value_counts().to_dict())
    output.with_suffix('.verification.json').write_text(json.dumps(counts,indent=2,allow_nan=False))
    data.to_csv(output.with_name(output.stem+'_replicates.csv.gz'),index=False,compression='gzip')
    print(summary[['cell','reading','size_pattern','rho','attempts','moment_conditional_coverage','interval_conditional_coverage','interval_reportability','interval_median_width']].to_string(index=False))
    print(json.dumps({k:v for k,v in counts.items() if k not in ('source_hashes','failures')},indent=2))
    return summary,data


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory');p.add_argument('output')
    a=p.parse_args();summarize(a.directory,a.output)
