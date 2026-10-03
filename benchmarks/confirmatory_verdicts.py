#!/usr/bin/env python3
"""Apply the rules of CONFIRMATORY_PROTOCOL.md to the outputs of the campaign.

    python benchmarks/confirmatory_verdicts.py --grid OUT/estimator_grid/estimator_benchmark.json \
        --order OUT/order --conditional OUT/conditional --output benchmarks/results_confirmatory

Writes ``verdicts.json`` (every design, its rate, its Clopper-Pearson bound and
its verdict under every rule that applies) and ``verdicts.csv`` beside it, and
prints the designs that fail a rule. Nothing is decided anywhere else.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

LIBERAL_COVERAGE, STRINGENT_COVERAGE = .925, .945
LIBERAL_RATE, STRINGENT_RATE = .075, .055
ALPHA = .05


def cp_lower(k: int, n: int, level: float) -> float:
    return 0. if k == 0 else float(stats.beta.ppf(level, k, n - k + 1))


def cp_upper(k: int, n: int, level: float) -> float:
    return 1. if k == n else float(stats.beta.ppf(1. - level, k + 1, n - k))


def coverage_verdicts(rows, study, procedure):
    """rows: (design, covered, reported, replicates) -> one verdict each."""
    out = []
    D = max(len(rows), 1)
    for design, covered, reported, replicates in rows:
        if reported == 0:
            out.append(dict(study=study, procedure=procedure, design=design, replicates=replicates,
                            reported=0, rate=None, bound=None, verdict='no interval reported'))
            continue
        low = cp_lower(covered, reported, ALPHA / D)
        high = cp_upper(covered, reported, ALPHA / D)
        verdict = ('meets' if low >= LIBERAL_COVERAGE else
                   'falls short' if high < LIBERAL_COVERAGE else 'undecided')
        out.append(dict(study=study, procedure=procedure, design=design, replicates=replicates,
                        reported=reported, rate=covered / reported, bound=low, upper=high,
                        mc_se=float(np.sqrt(covered / reported * (1 - covered / reported) / reported)),
                        verdict=verdict, stringent=bool(low >= STRINGENT_COVERAGE)))
    return out


def rate_verdicts(rows, study, procedure):
    """rows: (design, events, trials): the rate of rejection at 5 % or of a
    limit above its target, which must stay at most 0.075."""
    out = []
    D = max(len(rows), 1)
    for design, events, trials in rows:
        if trials == 0:
            continue
        high = cp_upper(events, trials, ALPHA / D)
        out.append(dict(study=study, procedure=procedure, design=design, replicates=trials,
                        rate=events / trials, bound=high,
                        mc_se=float(np.sqrt(events / trials * (1 - events / trials) / trials)),
                        verdict='meets' if high <= LIBERAL_RATE else 'fails',
                        stringent=bool(high <= STRINGENT_RATE)))
    return out


def theorem_verdicts(rows, study):
    """rows: (design, within, trials): every dataset within the bounds."""
    return [dict(study=study, procedure='theorem 1', design=d, replicates=n, exceptions=n - k,
                 verdict='holds' if k == n else 'refuted') for d, k, n in rows]


def study1(path: Path, rerun: Path | None = None):
    cells = json.loads(path.read_text())['cells']
    name = lambda c: (f"{'call' if c['binary'] else 'trait'} G={c['n_groups']} m={c['group_size']} "
                      f"share={c['share']} {c['unbalance']} {c['law']} p={c['prevalence']}")
    cover = [(name(c), c['clonal_share']['n_covered'], c['clonal_share']['n_intervals'], c['n_replicates'])
             for c in cells]
    # A cell whose coverage is undecided on its 2,000 datasets is decided on
    # 10,000 of the same stream (CONFIRMATORY_PROTOCOL.md, Section 3); the
    # rerun decides nothing else.
    if rerun is not None:
        again = {c['index']: c for c in (json.loads(p.read_text()) for p in sorted(rerun.glob('[0-9][0-9][0-9].json')))}
        first = {v['design']: v for v in coverage_verdicts(cover, 1, 'interval for the lineages in hand')}
        for i, c in enumerate(cells):
            if first[name(c)]['verdict'] == 'undecided':
                r = again.get(c['index'])
                if r is None:
                    raise SystemExit(f'study 1: cell {c["index"]} is undecided and has no rerun')
                if r['n_replicates'] < 10000:
                    raise SystemExit(f'study 1: the rerun of cell {c["index"]} has {r["n_replicates"]} datasets')
                cover[i] = (name(c), r['clonal_share']['n_covered'], r['clonal_share']['n_intervals'],
                            r['n_replicates'])
    level = [(name(c), c['n_rejections'], c['n_replicates']) for c in cells if c['share'] == 0]
    calls = [c for c in cells if c.get('call_latent_bounds')]
    theorem = [(name(c), c['call_latent_bounds']['sample_share_within_bounds'],
                c['call_latent_bounds']['n_bounded']) for c in calls]
    limit = [(name(c), c['call_latent_bounds']['n_bounded'] - c['call_latent_bounds']['limit_below_truth'],
              c['call_latent_bounds']['n_bounded']) for c in calls]
    return (coverage_verdicts(cover, 1, 'interval for the lineages in hand')
            + rate_verdicts(level, 1, 'level of the test')
            + theorem_verdicts(theorem, 1)
            + rate_verdicts(limit, 1, 'lower limit above rho'))


def _rows(directory: Path):
    by = defaultdict(list)
    for p in sorted(directory.rglob('*.json')):
        for row in json.loads(p.read_text()):
            by[row['design']].append(row)
    for design, rows in by.items():
        seen = [r['replicate'] for r in rows]
        if len(seen) != len(set(seen)):
            raise SystemExit(f'{directory}: design {design} has a replicate twice')
    return by


def study2(directory: Path):
    by = _rows(directory)
    cover, level, level_pooled, theorem, limit = [], [], [], [], []
    for design, rows in by.items():
        rep = [r for r in rows if r['observed_covered'] is not None]
        if rows[0]['rho'] > 0:
            cover.append((design, sum(r['observed_covered'] for r in rep), len(rep), len(rows)))
        else:
            level.append((design, sum(r['p_value'] <= ALPHA for r in rows), len(rows)))
            if 'pooled_p_value' in rows[0]:
                level_pooled.append((design, sum(r['pooled_p_value'] <= ALPHA for r in rows), len(rows)))
        bounded = [r for r in rows if np.isfinite(r['latent_lower'])]
        theorem.append((design, sum(r['sample_within_bounds'] for r in bounded), len(bounded)))
        limit.append((design, sum(not r['latent_limit_below'] for r in bounded), len(bounded)))
    pooled = rate_verdicts(level_pooled, 2, 'level of the test ignoring the laboratory')
    for v in pooled:
        v['verdict'] = 'reported: ' + v['verdict']
    return (coverage_verdicts(cover, 2, 'interval for the lineages in hand')
            + rate_verdicts(level, 2, 'level of the test within laboratories') + pooled
            + theorem_verdicts(theorem, 2) + rate_verdicts(limit, 2, 'lower limit above rho'))


def study3(directory: Path):
    out = []
    obs = _rows(directory / 'observed') if (directory / 'observed').exists() else {}
    cover = []
    for design, rows in obs.items():
        rep = [r for r in rows if r['reported']]
        cover.append((design, sum(bool(r['covered']) for r in rep), len(rep), len(rows)))
    out += coverage_verdicts(cover, 3, 'interval for the lineages in hand')
    lat = _rows(directory / 'latent') if (directory / 'latent').exists() else {}
    theorem, above_min, above_rho = [], [], []
    for design, rows in lat.items():
        bounded = [r for r in rows if np.isfinite(r['lower'])]
        theorem.append((design, sum(r['within'] for r in bounded), len(bounded)))
        above_rho.append((design, sum(not r['limit_below_rho'] for r in bounded), len(bounded)))
        with_min = [r for r in bounded if r['limit_below_rho_min'] is not None]
        if with_min:
            above_min.append((design, sum(not r['limit_below_rho_min'] for r in with_min), len(with_min)))
    out += theorem_verdicts(theorem, 3)
    out += rate_verdicts(above_min, 3, 'lower limit above rho_min')
    out += rate_verdicts(above_rho, 3, 'lower limit above rho')
    units = _rows(directory / 'units') if (directory / 'units').exists() else {}
    within = [(d, sum(r['p_within_units'] <= ALPHA for r in rows), len(rows)) for d, rows in units.items()]
    ignoring = [(d, sum(r['p_ignoring_units'] <= ALPHA for r in rows), len(rows)) for d, rows in units.items()]
    out += rate_verdicts(within, 3, 'level of the test within units')
    ignored = rate_verdicts(ignoring, 3, 'level of the test ignoring units')
    for v in ignored:
        v['verdict'] = 'reported: ' + v['verdict']
    out += ignored
    comp = _rows(directory / 'comparison') if (directory / 'comparison').exists() else {}
    cover = []
    for design, rows in comp.items():
        for term in rows[0]['truth']:
            # a term zero by construction in the design, zero in its truth
            # and in every estimate, is not scored
            if all(r['truth'][term] == 0. and r['estimate'][term] == 0. for r in rows):
                continue
            rep = [r for r in rows if r['covered'][term] is not None]
            cover.append((f'{design}: {term}', sum(r['covered'][term] for r in rep), len(rep), len(rows)))
    out += coverage_verdicts(cover, 3, 'paired intervals of the comparison')
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--grid', type=Path)
    ap.add_argument('--grid-rerun', type=Path,
                    help='cells/ of the 10,000-dataset rerun of the undecided cells of Study 1')
    ap.add_argument('--order', type=Path)
    ap.add_argument('--conditional', type=Path)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    verdicts = []
    if args.grid:
        verdicts += study1(args.grid, args.grid_rerun)
    if args.order:
        verdicts += study2(args.order)
    if args.conditional:
        verdicts += study3(args.conditional)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'verdicts.json').write_text(json.dumps(verdicts, indent=1) + '\n')
    keys = sorted({k for v in verdicts for k in v})
    with (args.output / 'verdicts.csv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(verdicts)
    bad = [v for v in verdicts if v['verdict'] in ('falls short', 'fails', 'refuted', 'undecided')]
    for v in bad:
        print(f"study {v['study']} | {v['procedure']} | {v['design']}: {v['verdict']} "
              f"(rate {v.get('rate')}, bound {v.get('bound')})")
    print(f'{len(verdicts)} verdicts, {len(bad)} not met')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
