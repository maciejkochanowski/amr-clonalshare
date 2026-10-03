"""The S. suis analyses beyond the configured run, on the recorded data.

The comparison of the population cluster with the sequence type (MLST) per
agent, with the paired intervals of its five terms; the decomposition of the
change in the ceftiofur call between the earlier and later UK isolates; and
the lineage share of the ordering of every agent's MIC readings within
testing laboratory x country of isolation, with the Benjamini-Yekutieli
selection across agents. Every analysis runs with the package's default
budgets.

    python benchmarks/empirical/ssuis_analysis.py [output dir]
"""
import sys
from pathlib import Path
import hashlib
import json
import time
import numpy as np
import pandas as pd
from amr_clonalshare.comparison import compare_lineage_definitions
from amr_clonalshare.clonality import decompose_prevalence_difference
from amr_clonalshare.jsonio import write_json

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT/'benchmarks/results_empirical/ssuis'
OUT.mkdir(parents=True, exist_ok=True)
DATA = ROOT/'examples/ssuis/data'
meta = pd.read_csv(DATA/'metadata.csv', dtype=str, keep_default_na=False).set_index('genome_id')
calls = pd.read_csv(DATA/'calls_long.csv', dtype=str, keep_default_na=False)
assert not calls.duplicated(['genome_id','antibiotic']).any()
calls['value'] = calls['call'].map({'1': 1., '0': 0.})
assert calls['value'].notna().all()
panel = calls.pivot(index='genome_id', columns='antibiotic', values='value').reindex(meta.index)
assert panel.shape == (677,13) and panel.notna().all().all()
rows, comparisons = [], {}
started = time.perf_counter()
for agent in panel:
    print('comparison', agent, flush=True)
    result = compare_lineage_definitions(panel[agent].to_numpy(), meta.baps_cluster, meta.mlst,
        ids=meta.index.to_numpy(), name_a='Population cluster', name_b='MLST',
        seed=20260915, folds=5, repeats=20, n_perm=999, n_boot=999)
    comparisons[agent] = result
    write_json(result, OUT/(agent+'_comparison.json'))
    for arm, item in result['arms'].items():
        r = item['share']
        rows.append(dict(agent=agent, arm=arm, n=item['n'], n_positive=item['n_positive'],
            groups=r['n_groups'], support=r['support'], estimate=r['kappa_adj'],
            lower=r['observed_low'], upper=r['observed_high'], estimable=r['estimable']))
pd.DataFrame(rows).to_csv(OUT/'lineage_comparisons.csv', index=False)
uk = meta.isolation_country.eq('United Kingdom') & meta.collection_year.ne('')
early = uk & pd.to_numeric(meta.collection_year, errors='coerce').lt(2013)
late = uk & pd.to_numeric(meta.collection_year, errors='coerce').ge(2013)
y = panel.ceftiofur
counts = []
for name, mask in [('Earlier',early), ('Later',late)]:
    for scope, selected in [('All records',mask), ('MLST available',mask & meta.mlst.ne(''))]:
        counts.append(dict(period=name, scope=scope, n=int(selected.sum()),
            positive=int(y[selected].sum()), prevalence=float(y[selected].mean())))
contrast = decompose_prevalence_difference(y[late], meta.baps_cluster[late],
    y[early], meta.baps_cluster[early], n_boot=1999, rng=np.random.default_rng(20260915))
write_json(contrast, OUT/'ceftiofur_decomposition.json')
pd.DataFrame(counts).to_csv(OUT/'ceftiofur_counts.csv', index=False)
files = [DATA/'metadata.csv', DATA/'calls_long.csv', Path(__file__)]
files += sorted((ROOT/'src/amr_clonalshare').glob('*.py'))
receipt = dict(seed=20260915, elapsed_seconds=time.perf_counter()-started,
    comparison_settings=dict(folds=5, repeats=20, n_perm=999, n_boot=999),
    decomposition_bootstraps=1999,
    phenotype='Retained threshold-defined binary outcomes; not reclassified clinical S/I/R',
    source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
write_json(receipt, OUT/'receipt.json')
print(json.dumps({'counts': counts, 'seconds': receipt['elapsed_seconds']}, indent=2))
