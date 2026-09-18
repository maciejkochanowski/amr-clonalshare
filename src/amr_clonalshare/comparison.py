"""Compare two recorded lineage definitions on a common, explicit isolate frame.

This is a descriptive sensitivity analysis, not a causal decomposition. The
four arms separate a change of labels from a change of included records. The
reported differences have no paired confidence interval.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from numbers import Integral
from pathlib import Path
from typing import Any

import numpy as np

from .attribution import _is_untyped, clonal_share
from .jsonio import to_jsonable, write_json

__all__ = ['compare_lineage_definitions', 'main']


def compare_lineage_definitions(y, lineage_a, lineage_b, *, ids=None,
                                name_a: str = 'A', name_b: str = 'B',
                                seed: int = 42, folds: int = 5,
                                repeats: int = 20, n_perm: int = 200,
                                n_boot: int = 0) -> dict[str, Any]:
    """Compare available and common-record adjusted Brier shares.

    Outcomes must be 0/1/NaN. Supplied IDs must be unique; their string order
    fixes record order before random splits so input row order cannot change
    the result. With omitted IDs, positions identify records. Common-arm
    outcomes, IDs and random split seeds are identical, but group labels vary.
    A nonestimable arm retains diagnostics; it is not promoted to a valid
    estimate by this comparison. Individual arm intervals are not intervals
    for the paired differences. No sampling or causal interpretation is added.
    """
    values = np.asarray(y, dtype=float)
    a, b = np.asarray(lineage_a, dtype=object), np.asarray(lineage_b, dtype=object)
    if values.ndim != 1 or a.ndim != 1 or b.ndim != 1 or not (len(values) == len(a) == len(b)):
        raise ValueError('outcome and both lineage vectors must be one-dimensional and equally long')
    if np.any(np.isinf(values)) or np.any(np.isfinite(values) & ~np.isin(values, [0., 1.])):
        raise ValueError('outcomes must be binary 0/1 or NaN')
    if isinstance(seed, bool) or not isinstance(seed, Integral) or seed < 0:
        raise ValueError('seed must be a nonnegative integer shared by all four arms')
    if not name_a or not name_b or name_a == name_b:
        raise ValueError('lineage names must be nonempty and distinct')
    raw_ids = np.arange(values.size).astype(str) if ids is None else np.asarray(ids, dtype=object)
    if raw_ids.ndim != 1 or raw_ids.size != values.size:
        raise ValueError('ids must have one entry per outcome')
    if any(x is None or str(x).strip() == '' or str(x) in ('<NA>', 'nan', 'NaT') for x in raw_ids):
        raise ValueError('ids must be nonmissing and unique')
    names = np.asarray([str(x) for x in raw_ids])
    if np.unique(names).size != names.size:
        raise ValueError('ids must be unique after string conversion')
    order = np.argsort(names, kind='stable')
    names, values, a, b = names[order], values[order], a[order], b[order]
    observed = np.isfinite(values)
    available_a = observed & np.array([not _is_untyped(x) for x in a], dtype=bool)
    available_b = observed & np.array([not _is_untyped(x) for x in b], dtype=bool)
    common = available_a & available_b
    settings = dict(seed=int(seed), folds=folds, repeats=repeats, n_perm=n_perm, n_boot=n_boot)

    def arm(mask, labels, label_name):
        result = clonal_share(values[mask], labels[mask], **settings)
        digest = hashlib.sha256(json.dumps(names[mask].tolist(), ensure_ascii=False,
                                          separators=(',', ':')).encode()).hexdigest()
        return {'lineage_definition': label_name, 'n': int(mask.sum()),
                'n_positive': int(values[mask].sum()), 'id_sha256': digest,
                'share': result.as_dict()}

    arms = {'a_available': arm(available_a, a, name_a),
            'a_common': arm(common, a, name_a),
            'b_common': arm(common, b, name_b),
            'b_available': arm(available_b, b, name_b)}
    points = [arms[k]['share']['kappa_adj'] for k in arms]
    difference = dict.fromkeys(['selection_from_a', 'label_difference_on_common',
                               'selection_to_b', 'total_difference', 'identity_residual'])
    if np.isfinite(points).all():
        x0, x1, x2, x3 = points
        difference.update(selection_from_a=x1-x0, label_difference_on_common=x2-x1,
                          selection_to_b=x3-x2, total_difference=x3-x0,
                          identity_residual=(x1-x0)+(x2-x1)+(x3-x2)-(x3-x0))
    all_eligible = all(x['share']['estimable'] for x in arms.values())
    status = ('computed' if all_eligible else 'diagnostic_only: one or more arms not estimable')
    return to_jsonable({
        'schema_version': 'lineage_comparison_1.0', 'status': status,
        'counts': dict(supplied=int(values.size), observed=int(observed.sum()),
                       available_a=int(available_a.sum()), available_b=int(available_b.sum()),
                       common=int(common.sum())),
        'settings': settings, 'names': {'a': name_a, 'b': name_b}, 'arms': arms,
        'difference_decomposition': difference,
        'paired_uncertainty': 'not computed; individual intervals are not intervals for differences',
        'interpretation': ('Descriptive, not causal. The middle difference changes labels on the same '
                           'observed records. The outer differences change record inclusion. Common '
                           'records need not represent either original collection or a population.'),
    })


def main(argv=None) -> int:
    """CLI for an aligned, one-row-per-isolate CSV; output publication is atomic."""
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ('input', 'id-column', 'outcome', 'lineage-a', 'lineage-b', 'output'):
        parser.add_argument('--' + option, required=True)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--folds', type=int, default=5)
    parser.add_argument('--repeats', type=int, default=20)
    parser.add_argument('--permutations', type=int, default=200)
    parser.add_argument('--bootstraps', type=int, default=0)
    args = parser.parse_args(argv)
    try:
        import pandas as pd
        from .outputs import _publish, _manifest, validate_destination
        validate_destination(args.output)
        data = pd.read_csv(args.input, dtype=str, keep_default_na=False)
        required = [args.id_column, args.outcome, args.lineage_a, args.lineage_b]
        missing = set(required) - set(data.columns)
        if missing:
            raise ValueError(f'CSV columns missing: {sorted(missing)}')
        raw = data[args.outcome].str.strip()
        raw = raw.mask(raw.str.lower().isin(['', 'nan', 'na', 'n/a', 'null']))
        outcome = pd.to_numeric(raw, errors='raise').to_numpy(dtype=float)
        result = compare_lineage_definitions(
            outcome, data[args.lineage_a], data[args.lineage_b], ids=data[args.id_column],
            name_a=args.lineage_a, name_b=args.lineage_b, seed=args.seed,
            folds=args.folds, repeats=args.repeats, n_perm=args.permutations, n_boot=args.bootstraps)
        result['input_sha256'] = hashlib.sha256(Path(args.input).read_bytes()).hexdigest()

        def writer(directory):
            write_json(result, directory / 'comparison.json')
            with (directory / 'arms.csv').open('w', newline='', encoding='utf-8') as handle:
                fields = ['arm', 'definition', 'n', 'n_positive', 'n_groups', 'n_groups_repeated',
                          'effective_groups', 'largest_group_share', 'support',
                          'kappa_adj', 'estimable', 'id_sha256']
                table = csv.DictWriter(handle, fieldnames=fields)
                table.writeheader()
                for key, item in result['arms'].items():
                    share = item['share']
                    table.writerow(dict(arm=key, definition=item['lineage_definition'], n=item['n'],
                                        n_positive=item['n_positive'], n_groups=share['n_groups'],
                                        n_groups_repeated=share.get('n_groups_repeated'),
                                        effective_groups=share.get('effective_groups'),
                                        largest_group_share=share.get('largest_group_share'),
                                        support=share['support'], kappa_adj=share['kappa_adj'],
                                        estimable=share['estimable'], id_sha256=item['id_sha256']))
            text = ('# Matched-record lineage comparison\n\n' + result['interpretation'] + '\n\n'
                    + 'Status: ' + result['status'] + '\n\n' + result['paired_uncertainty'] + '\n\n'
                    + 'Common observed records: ' + str(result['counts']['common']) + '.\n\n'
                    + 'Read `arms.csv` with the eligibility flags; `comparison.json` records '
                    + 'the four arms, split settings and descriptive differences.\n')
            (directory / 'report.md').write_text(text, encoding='utf-8')
            _manifest(directory, 'matched_lineage_comparison')
        _publish(args.output, writer)
    except (OSError, ValueError, TypeError) as exc:
        print(f'Comparison not completed: {exc}', file=sys.stderr)
        return 2
    print(f'Comparison saved to {args.output}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
