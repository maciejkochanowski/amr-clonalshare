"""Compare two recorded lineage definitions on a common, explicit isolate frame.

This is a descriptive sensitivity analysis, not a causal decomposition. Six
analyses separate a change of the included records from a change of labels.
The share is scored on the isolates of repeated lineages, so a definition that
leaves more isolates in singleton lineages also scores fewer isolates; the
records shared by both definitions are therefore narrowed to those that sit in
a repeated lineage under both (the common scorable records, found by removing
isolates until every remaining one has a partner under each definition), and
the difference splits into five terms that sum to it exactly: record selection
into the shared frame, isolates set aside as singletons under the first
definition, relabelling on the common scorable records, isolates set aside as
singletons under the second definition, and record selection out of the
shared frame. Only the middle term changes the labels on the same scored
isolates.

The terms have paired intervals. The records and both labellings are held
fixed, and the outcome of every record is drawn again from the smoothed distribution
of its cell, the records that share both labels (Jeffreys' posterior mean
rate, (positives + 1/2) / (records + 1)); the six analyses are recomputed on
the same drawn outcomes, each set against its own permuted control, and every
term is studentized by its own standard error, the infinitesimal jackknife of
the shares over the records with the variance of every cell taken under its
smoothed distribution, and centred on the terms of the drawn scenario, computed exactly
with the distribution of every lineage that of a record drawn from it, so that a
lineage whose records sit in cells of different rates carries that spread
within it: the bootstrap-t interval, as for the share itself. Arms that score
the same records under the same labels are one analysis and share their fold
noise, so a term between them is zero in every draw. An interval is therefore
for the change the terms describe on these records and labels, with the
outcomes as the random part.
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

from .attribution import _codes, _debias, _fold_average, _is_untyped, _jumped, _share_derivatives, clonal_share
from .jsonio import to_jsonable, write_json

__all__ = ['compare_lineage_definitions', 'main']


TERMS = ['selection_from_a', 'singletons_under_a', 'label_difference_on_scorable',
         'singletons_under_b', 'selection_to_b']


def _arm_truth(p, labels, mask):
    """Lineage share of the represented lineages of one arm when record i is positive with
    probability p_i: the repeated lineages of the arm, records weighted alike.
    The distribution of a lineage is that of a record drawn from it, positive with the
    mean rate mu of its records, so its variance is mu (1 - mu): a lineage
    whose records sit in cells of different rates carries that spread within
    it, as the share of Appendix A.2.1 of the article defines it."""
    code = _codes(labels[mask])
    sizes = np.bincount(code).astype(float)
    keep = sizes[code] >= 2
    code, q = _codes(labels[mask][keep]), p[mask][keep]
    sizes = np.bincount(code).astype(float)
    if sizes.size < 2:
        return float('nan')
    w = sizes / sizes.sum()
    mu = np.bincount(code, weights=q) / sizes
    within = float(w @ (mu * (1. - mu)))
    between = float(w @ (mu - w @ mu) ** 2)
    return between / (between + within) if between + within > 0 else float('nan')


def _scored(labels, mask):
    """The records of ``mask`` the arm scores: those in a repeated lineage."""
    out = np.zeros(mask.size, dtype=bool)
    code = _codes(labels[mask])
    out[np.flatnonzero(mask)] = np.bincount(code)[code] >= 2
    return out


def _arm_share(y, labels, mask, share):
    """The share of one arm on drawn outcomes: its average over the fold
    assignments, set against the arm's permuted control from the data."""
    code = _codes(labels[mask])
    sizes = np.bincount(code)
    keep = sizes[code] >= 2
    values = y[mask][keep]
    if values.size < 2 or np.ptp(values) <= 0:
        return float('nan')
    code = _codes(labels[mask][keep])
    return _debias(_fold_average(values, code, share['fold_ratio']), share['null_mean'])


#: The contrasts of the six arm shares: the five terms, the change on the
#: common records and the total.
_CONTRASTS = np.array([[-1, 1, 0, 0, 0, 0], [0, -1, 1, 0, 0, 0], [0, 0, -1, 1, 0, 0],
                       [0, 0, 0, -1, 1, 0], [0, 0, 0, 0, -1, 1], [0, -1, 0, 0, 1, 0],
                       [-1, 0, 0, 0, 0, 1]], dtype=float)


def _arm_derivatives(y, labels, mask, share):
    """The derivative of one arm's share with respect to the weight of every
    record at outcome 0 and at outcome 1 (:func:`_share_derivatives`), zero
    for records the arm does not score: (records, 2)."""
    out = np.zeros((y.size, 2))
    idx = np.flatnonzero(mask)
    code = _codes(labels[mask])
    keep = np.bincount(code)[code] >= 2
    idx = idx[keep]
    values = y[idx]
    if values.size < 2 or np.ptp(values) <= 0:
        return None
    code = _codes(labels[idx])
    counts = np.zeros((int(code.max()) + 1, 2))
    np.add.at(counts, (code, values.astype(int)), 1.)
    D = _share_derivatives(np.array([0., 1.]), counts, share['fold_ratio'], share['null_mean'])
    if D is None:
        return None
    out[idx] = D[code]
    return out


def _contrast_se(y, labels, masks, shares, cell, spread, keys):
    """Standard errors of the seven contrasts: the jackknife over the records,
    the variance of every cell (records sharing both labels) under its
    Jeffreys-smoothed rate, and the fold spread of the distinct arms."""
    arms = [_arm_derivatives(y, lab, m, share) for lab, m, share in zip(labels, masks, shares)]
    if any(arm is None for arm in arms):
        return None
    D = np.einsum('ca,iak->cik', _CONTRASTS, np.stack(arms, axis=1))     # contrast, record, outcome
    seen = np.isfinite(y)
    C = int(cell.max()) + 1
    m = np.bincount(cell[seen], minlength=C).astype(float)
    rate = (np.bincount(cell[seen], weights=y[seen], minlength=C) + .5) / (m + 1.)
    first = np.full(C, -1)
    first[cell[seen][::-1]] = np.flatnonzero(seen)[::-1]
    variance = np.zeros(_CONTRASTS.shape[0])
    present = first >= 0
    # every record of a cell has the same two derivatives
    gap = D[:, first[present], 1] - D[:, first[present], 0]
    q = rate[present]
    size = m[present]
    variance += np.sum(size * size / np.maximum(size - 1., 1.) * q * (1. - q) * gap * gap, axis=1)
    # the fold noise of the distinct analyses
    distinct = {}
    for j, key in enumerate(keys):
        distinct.setdefault(key, []).append(j)
    for members in distinct.values():
        coefficient = _CONTRASTS[:, members].sum(1)
        variance += (coefficient * spread[members[0]]) ** 2
    return np.sqrt(variance)


def _paired_intervals(values, a, b, masks, points, shares, *, repeats, n_boot, rng):
    """Bootstrap-t intervals of the five terms, the change on the common
    records and the total (module docstring), and the draws that gave all
    six shares."""
    observed = np.isfinite(values)
    cell = _codes(np.asarray([f'{x!r}|{z!r}' for x, z in zip(a, b)], dtype=object))
    positives = np.bincount(cell, weights=np.where(observed, values, 0.))
    records = np.bincount(cell, weights=observed.astype(float))
    p = ((positives + .5) / (records + 1.))[cell]
    labels = [a, a, a, b, b, b]
    truth = _CONTRASTS @ np.array([_arm_truth(p, lab, m) for lab, m in zip(labels, masks)])
    estimate = _CONTRASTS @ np.asarray(points, dtype=float)
    # arms that score the same records under the same labels are one
    # analysis, whatever singletons they set aside: the same fold assignments
    # gave them the same estimate, and they share their fold noise
    keys = [(id(lab), _scored(lab, m).tobytes()) for lab, m in zip(labels, masks)]
    first = {key: keys.index(key) for key in keys}
    spread = np.array([share['cv_sd'] / np.sqrt(repeats) / (1. - share['null_mean'])
                       if np.isfinite(share['cv_sd']) else 0. for share in shares])
    names = TERMS + ['definition_difference_on_common', 'total_difference']
    se_data = _contrast_se(values, labels, masks, shares, cell, spread, keys)
    if se_data is None or not np.isfinite(truth).all():
        return {name: None for name in names}, 0
    draws = []
    for _ in range(n_boot):
        y = np.where(observed, (rng.random(values.size) < p).astype(float), np.nan)
        noise = rng.normal(0., 1., len(keys)) * spread
        values_b = np.array([_arm_share(y, lab, m, share) + noise[first[key]]
                             for lab, m, share, key in zip(labels, masks, shares, keys)])
        if not np.isfinite(values_b).all():
            continue
        se = _contrast_se(y, labels, masks, shares, cell, spread, keys)
        if se is None:
            continue
        deviation = _CONTRASTS @ values_b - truth
        # a term that is zero by construction has no spread and no deviation
        draws.append(np.where(se > 0, deviation / np.where(se > 0, se, 1.), 0.))
    if len(draws) < .9 * n_boot or len(draws) < 2:
        return {name: None for name in names}, len(draws)
    t = np.asarray(draws)
    low = estimate - se_data * np.quantile(t, .975, axis=0)
    high = estimate - se_data * np.quantile(t, .025, axis=0)
    return {name: [float(lo), float(hi)] for name, lo, hi in zip(names, low, high)}, len(draws)


def compare_lineage_definitions(y, lineage_a, lineage_b, *, ids=None,
                                name_a: str = 'A', name_b: str = 'B',
                                seed: int = 42, folds: int = 5,
                                repeats: int = 20, n_perm: int = 999,
                                n_boot: int = 999) -> dict[str, Any]:
    """Compare the lineage share of a call under two lineage definitions.

    Outcomes must be 0/1/NaN. Supplied IDs must be unique; their string order
    fixes record order before random splits so input row order cannot change
    the result. With omitted IDs, positions identify records. Common-arm
    outcomes, IDs and random split seeds are identical, but group labels vary.
    A nonestimable arm retains diagnostics; it is not promoted to a valid
    estimate by this comparison. ``n_boot`` draws give the interval of every
    arm and the paired intervals of the terms (module docstring); 0 gives
    points only. No causal interpretation is added.
    """
    values = np.asarray(y, dtype=float)
    a, b = np.asarray(lineage_a, dtype=object), np.asarray(lineage_b, dtype=object)
    if values.ndim != 1 or a.ndim != 1 or b.ndim != 1 or not (len(values) == len(a) == len(b)):
        raise ValueError('outcome and both lineage vectors must be one-dimensional and equally long')
    if np.any(np.isinf(values)) or np.any(np.isfinite(values) & ~np.isin(values, [0., 1.])):
        raise ValueError('outcomes must be binary 0/1 or NaN')
    if isinstance(seed, bool) or not isinstance(seed, Integral) or seed < 0:
        raise ValueError('seed must be a nonnegative integer shared by all six analyses')
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
    scorable = common.copy()
    while True:
        # an isolate stays while its lineage has another member left under
        # both definitions; removals can strand a partner, so repeat
        kept = scorable.copy()
        for labels in (a, b):
            text = np.asarray([str(x) for x in labels], dtype=object)
            _, code, counts = np.unique(text[kept], return_inverse=True, return_counts=True)
            ok = np.zeros(values.size, dtype=bool)
            ok[np.flatnonzero(kept)] = counts[code.ravel()] >= 2
            kept &= ok
        if (kept == scorable).all():
            break
        scorable = kept
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
            'a_scorable': arm(scorable, a, name_a),
            'b_scorable': arm(scorable, b, name_b),
            'b_common': arm(common, b, name_b),
            'b_available': arm(available_b, b, name_b)}
    points = [arms[k]['share']['kappa_adj'] for k in arms]
    terms = TERMS
    difference = dict.fromkeys(terms + ['definition_difference_on_common', 'total_difference',
                                        'identity_residual'])
    if np.isfinite(points).all():
        steps = np.diff(points)
        difference.update(dict(zip(terms, steps.tolist())))
        difference.update(definition_difference_on_common=points[4] - points[1],
                          total_difference=points[-1] - points[0],
                          identity_residual=float(steps.sum() - (points[-1] - points[0])))
    all_eligible = all(x['share']['estimable'] for x in arms.values())
    status = ('computed' if all_eligible else 'diagnostic_only: one or more arms not estimable')
    paired: dict[str, Any] = {name: None for name in TERMS + ['definition_difference_on_common',
                                                            'total_difference']}
    used = 0
    if n_boot and all_eligible and np.isfinite(points).all():
        masks = [available_a, common, scorable, scorable, common, available_b]
        rng = np.random.Generator(_jumped(np.random.default_rng([int(seed), 2])))
        paired, used = _paired_intervals(values, a, b, masks, points,
                                         [arms[k]['share'] for k in arms],
                                         repeats=repeats, n_boot=int(n_boot), rng=rng)
    return to_jsonable({
        'schema_version': 'lineage_comparison_1.1', 'status': status,
        'counts': dict(supplied=int(values.size), observed=int(observed.sum()),
                       available_a=int(available_a.sum()), available_b=int(available_b.sum()),
                       common=int(common.sum()), common_scorable=int(scorable.sum())),
        'settings': settings, 'names': {'a': name_a, 'b': name_b}, 'arms': arms,
        'difference_decomposition': difference,
        'paired_intervals': paired,
        'paired_n_boot_used': int(used),
        'paired_uncertainty': ('95% bootstrap-t intervals, records and labels fixed, outcomes drawn '
                               'again from the smoothed distribution of the records sharing both labels'
                               if any(v is not None for v in paired.values()) else
                               'not computed'),
        'interpretation': ('Descriptive, not causal. The middle difference changes the labels on the '
                           'same scored isolates, the records that sit in a repeated lineage under both '
                           'definitions. The terms beside it are the isolates each definition sets aside '
                           'as singletons, and the outer terms change record inclusion. The difference '
                           'on all common records is the sum of the three middle terms. Common records '
                           'need not represent either original collection or a population.'),
    })


def _parser() -> argparse.ArgumentParser:
    from . import __version__
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version=f'amr-clonalshare {__version__}')
    for option, text in (('input', 'CSV with one row per isolate'),
                         ('id-column', 'column naming the isolate'),
                         ('outcome', 'column holding the 0/1 outcome'),
                         ('lineage-a', 'column holding the first lineage definition'),
                         ('lineage-b', 'column holding the second lineage definition'),
                         ('output', 'new directory for the results')):
        parser.add_argument('--' + option, required=True, help=text)
    parser.add_argument('--seed', type=int, default=42, help='random seed (default 42)')
    parser.add_argument('--folds', type=int, default=5, help='cross-validation folds (default 5)')
    parser.add_argument('--repeats', type=int, default=20, help='fold assignments averaged (default 20)')
    parser.add_argument('--permutations', type=int, default=999,
                        help='permuted-label runs for the control and p-value, at least 1 (default 999)')
    parser.add_argument('--bootstraps', '--bootstrap', type=int, default=999,
                        help='bootstrap draws for the interval of every arm and the paired intervals '
                             'of the terms; 0 reports points only (default 999)')
    return parser


def run_comparison(argv) -> Path:
    """Run the comparison the command line describes and publish it; raises on a refusal.

    ``main`` wraps this with the messages of the command line; the local form
    (``gui.py``) calls it directly and shows the exception it raises."""
    args = _parser().parse_args(argv)
    if args.permutations < 1 or args.folds < 2 or args.repeats < 1 or args.bootstraps < 0:
        raise ValueError('--permutations and --repeats must be at least 1, --folds at least 2 '
                         'and --bootstraps not negative')
    import pandas as pd
    from .io import _read_table
    from .outputs import _publish, _manifest, validate_destination
    validate_destination(args.output)
    # read as the analyses read their tables: comma, semicolon or tab, or a workbook
    data = _read_table(args.input, 'comparison table', dtype=str)
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
                      'n_singletons_set_aside', 'effective_groups', 'largest_group_share', 'support',
                      'kappa_adj', 'observed_low', 'observed_high', 'p_value', 'null_mean',
                      'estimable', 'id_sha256']
            table = csv.DictWriter(handle, fieldnames=fields)
            table.writeheader()
            for key, item in result['arms'].items():
                share = item['share']
                table.writerow(dict(arm=key, definition=item['lineage_definition'], n=item['n'],
                                    n_positive=item['n_positive'], n_groups=share['n_groups'],
                                    n_groups_repeated=share.get('n_groups_repeated'),
                                    n_singletons_set_aside=share.get('n_singletons_set_aside'),
                                    effective_groups=share.get('effective_groups'),
                                    largest_group_share=share.get('largest_group_share'),
                                    support=share['support'], kappa_adj=share['kappa_adj'],
                                    observed_low=share.get('observed_low'),
                                    observed_high=share.get('observed_high'),
                                    p_value=share.get('p_value'), null_mean=share.get('null_mean'),
                                    estimable=share['estimable'], id_sha256=item['id_sha256']))
        contrasts = (('selection_from_a', 'Record selection into the shared frame'),
                     ('singletons_under_a', 'Isolates set aside as singletons, first definition'),
                     ('label_difference_on_scorable', 'Relabelling on the common scored isolates'),
                     ('singletons_under_b', 'Isolates set aside as singletons, second definition'),
                     ('selection_to_b', 'Record selection out of the shared frame'),
                     ('definition_difference_on_common', 'Change of definition on the shared records'),
                     ('total_difference', 'Total difference'),
                     ('identity_residual', 'Residual of the identity'))
        def limits(key):
            pair = (result.get('paired_intervals') or {}).get(key)
            return '' if not pair else format(pair[0], '+.4f') + ' to ' + format(pair[1], '+.4f')
        rows = ''.join(
            '| ' + title + ' | '
            + ('not available' if result['difference_decomposition'][key] is None
               else format(result['difference_decomposition'][key], '+.4f')) + ' | '
            + limits(key) + ' |\n'
            for key, title in contrasts)
        arm_titles = (('a_available', 'all records the first definition labels'),
                      ('a_common', 'records both definitions label'),
                      ('a_scorable', 'records scored under both definitions'),
                      ('b_scorable', 'records scored under both definitions'),
                      ('b_common', 'records both definitions label'),
                      ('b_available', 'all records the second definition labels'))

        def figure(value, spec='.3f'):
            return '' if value is None else format(value, spec)

        arm_rows = ''.join(
            '| ' + result['arms'][key]['lineage_definition'] + ' | ' + title + ' | '
            + str(result['arms'][key]['n']) + ' | '
            + figure(result['arms'][key]['share'].get('kappa_adj')) + ' | '
            + ('' if result['arms'][key]['share'].get('observed_low') is None else
               figure(result['arms'][key]['share']['observed_low']) + ' to '
               + figure(result['arms'][key]['share'].get('observed_high'))) + ' | '
            + figure(result['arms'][key]['share'].get('p_value')) + ' |\n'
            for key, title in arm_titles if key in result['arms'])
        settings = result.get('settings') or {}
        text = ('# Matched-record lineage comparison\n\n' + result['interpretation'] + '\n\n'
                + 'Status: ' + result['status'] + '\n\n'
                + 'Records supplied: ' + str(result['counts'].get('supplied'))
                + '; labelled by the first definition: ' + str(result['counts'].get('available_a'))
                + ', by the second: ' + str(result['counts'].get('available_b'))
                + '; common observed records: ' + str(result['counts']['common'])
                + ', of which scored under both definitions: '
                + str(result['counts']['common_scorable']) + '. Seed ' + str(settings.get('seed'))
                + ', ' + str(settings.get('n_perm')) + ' permutations, ' + str(settings.get('n_boot'))
                + ' bootstrap draws.\n\n'
                + '| Definition | Records | n | Lineage share | 95% interval | p |\n'
                + '| --- | --- | ---: | ---: | ---: | ---: |\n' + arm_rows + '\n'
                + '| Contrast | Difference in the lineage share | 95% interval |\n| --- | ---: | ---: |\n'
                + rows + '\n'
                + 'The first five contrasts sum to the total difference, and the three '
                + 'middle ones to the change of definition on the shared records; the residual '
                + 'records that identity. Interval for a contrast: '
                + result['paired_uncertainty'] + '.\n\n'
                + 'Read `arms.csv` with the eligibility flags; `comparison.json` records '
                + 'the six analyses, split settings and descriptive differences.\n')
        (directory / 'report.md').write_text(text, encoding='utf-8')
        _manifest(directory, 'matched_lineage_comparison')
    _publish(args.output, writer)
    return Path(args.output)


def main(argv=None) -> int:
    """CLI for an aligned, one-row-per-isolate CSV; output publication is atomic."""
    try:
        output = run_comparison(argv)
    except (OSError, ValueError, TypeError) as exc:
        print(f'Comparison not completed: {exc}', file=sys.stderr)
        return 2
    print(f'Comparison saved to {output}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
