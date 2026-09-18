"""Verify null-bootstrap evidence and distinguish test coverage from interval coverage."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from benchmarks.mic_inference.summarize import summarize, wilson


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--tasks', type=Path)
    args = parser.parse_args()
    files = sorted(args.directory.glob('cell_*.csv'))
    for path in files:
        receipt = json.loads(path.with_suffix('.receipt.json').read_text())
        matrix = path.with_suffix('.bootstrap.npz')
        if hashlib.sha256(matrix.read_bytes()).hexdigest() != receipt['bootstrap_sha256']:
            raise ValueError(f'bootstrap checksum mismatch: {path}')
        for name, digest in receipt['bounds_artifacts'].items():
            if hashlib.sha256((path.parent/name).read_bytes()).hexdigest() != digest:
                raise ValueError(f'interval inversion checksum mismatch: {name}')
    summary, data = summarize(args.directory, args.output)
    if args.tasks is not None:
        expected = set()
        for line in args.tasks.read_text().splitlines():
            cell, start, reps, boot, bounds = map(int, line.split())
            expected.update((cell, rep) for rep in range(start, start+4*reps))
        actual = set(zip(data.cell.astype(int), data.replicate.astype(int)))
        if actual != expected:
            raise ValueError(f'incomplete campaign: {len(expected-actual)} missing, {len(actual-expected)} extra')
    details = []
    for (phase, cell), group in data.groupby(['phase', 'cell']):
        ok = group.get('test_reportable', group.new_reportable).astype(bool)
        accepted = group.get('pointwise_accepted', group.new_test_covered).astype(bool)
        success = int((ok & accepted).sum())
        count = int(ok.sum())
        low, high = wilson(success, count)
        computed = group.bounds_computed.astype(bool)
        resolved = computed & group.new_reportable.astype(bool)
        anti = group.anti_conservative_inversion_mismatch.astype(bool)
        rows = group.loc[resolved]
        details.append(dict(phase=phase, cell=cell,
            test_reportable=count, test_covered=success,
            test_conditional_coverage=success/count if count else np.nan,
            test_wilson_low=low, test_wilson_high=high,
            test_coverage_including_no_rejections=float(accepted.mean()),
            inversion_attempted=int(computed.sum()), inversion_resolved=int(resolved.sum()),
            inversion_covered=int((resolved & group.new_covered).sum()),
            inversion_median_width=float((rows.new_high-rows.new_low).median()) if len(rows) else np.nan,
            anti_conservative_inversion_mismatches=int(anti.sum())))
    summary = summary.merge(pd.DataFrame(details), on=['phase', 'cell'], validate='one_to_one')
    summary.to_csv(args.output, index=False)
    print(summary[['cell', 'test_conditional_coverage', 'test_reportable',
        'inversion_attempted', 'inversion_resolved', 'anti_conservative_inversion_mismatches']].to_string(index=False))
    verification_path = args.output.with_suffix('.verification.json')
    verification = json.loads(verification_path.read_text())
    verification.update(bootstrap_and_inversion_hashes_verified=True,
        task_table_complete=args.tasks is not None,
        pointwise_test_is_primary=True,
        inversion_attempted=int(summary.inversion_attempted.sum()),
        inversion_resolved=int(summary.inversion_resolved.sum()),
        anti_conservative_inversion_mismatches=int(summary.anti_conservative_inversion_mismatches.sum()))
    verification_path.write_text(json.dumps(verification, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
