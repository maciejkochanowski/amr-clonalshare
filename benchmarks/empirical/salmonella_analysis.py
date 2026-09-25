"""Compare two lineage definitions of the poultry-meat Salmonella example on matched records."""
import os
from pathlib import Path
import hashlib
import time
import numpy as np
import pandas as pd
from amr_clonalshare.comparison import compare_lineage_definitions
from amr_clonalshare.jsonio import write_json

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT/'examples/salmonella_poultry/data'
OUT = Path(os.environ.get('AMR_RESULTS_OUT', ROOT/'benchmarks/results_mic_release'))/'external'
OUT.mkdir(parents=True, exist_ok=True)
meta = pd.read_csv(DATA/'metadata.csv', dtype=str, keep_default_na=False).set_index('genome_id')
calls = pd.read_csv(DATA/'calls_long.csv', dtype=str, keep_default_na=False)
selected = calls[calls.antibiotic.eq('nalidixic acid')].copy()
if selected.empty:
    selected = calls[calls.antibiotic.str.replace('_',' ').eq('nalidixic acid')].copy()
assert not selected.empty and not selected.genome_id.duplicated().any()
selected['value'] = selected.call.map({'1': 1., '0': 0.})
assert selected.value.notna().all()
y = selected.set_index('genome_id').value.reindex(meta.index)
started = time.perf_counter()
result = compare_lineage_definitions(y.to_numpy(), meta.pds_cluster, meta.serovar,
    ids=meta.index.to_numpy(), name_a='SNP cluster', name_b='Serovar',
    seed=20260915, folds=5, repeats=20, n_perm=199, n_boot=399)
result['source_definition'] = 'Archived I/R-combined non-susceptibility; not harmonized R-only clinical resistance'
write_json(result, OUT/'nalidixic_matched_comparison.json')
rows = []
for arm, item in result['arms'].items():
    r = item['share']
    rows.append(dict(arm=arm, n=item['n'], n_positive=item['n_positive'],
        groups=r['n_groups'], support=r['support'], estimate=r['kappa_adj'],
        lower=r['ci_low'], upper=r['ci_high'], estimable=r['estimable']))
pd.DataFrame(rows).to_csv(OUT/'nalidixic_arms.csv', index=False)
inputs = [DATA/'metadata.csv', DATA/'calls_long.csv', Path(__file__)]
inputs += list((ROOT/'src/amr_clonalshare').glob('*.py'))
write_json(dict(elapsed_seconds=time.perf_counter()-started,
    input_frame=len(meta), observed=int(np.isfinite(y).sum()), positive=int(y.sum()),
    source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}), OUT/'receipt.json')
print(pd.DataFrame(rows).to_string(index=False))
print(result['difference_decomposition'])
