#!/usr/bin/env python3
"""Per-design operating characteristics of the confirmatory campaign.

    python benchmarks/confirmatory_summaries.py --grid OUT/estimator_grid/estimator_benchmark.json \
        [--grid-rerun OUT/estimator_grid_rerun/cells] --order OUT/order \
        --conditional OUT/conditional --output benchmarks/results_confirmatory

The verdicts of the protocol are those of ``confirmatory_verdicts.py``; this
script writes what Section 4 of the protocol reports beside them and decides
nothing: for every design the number of datasets, the mean estimate and its
target, the median width of an interval and the share of datasets without
one, and for the lower limit the rate at which it is positive and the median
ratio of the limit to its target. Writes ``design_summaries.json`` and the
digest of every raw output file it read (``raw_outputs_sha256.json``), so that
the raw outputs, which are too large to ship with the package, can be checked
against the summaries.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def _rows(directory: Path, digests: dict):
    by = defaultdict(list)
    for p in sorted(directory.rglob('*.json')):
        digests[str(p.relative_to(directory.parent))] = hashlib.sha256(p.read_bytes()).hexdigest()
        for row in json.loads(p.read_text()):
            by[row['design']].append(row)
    return by


def _median(values):
    values = np.asarray([v for v in values if v is not None and np.isfinite(v)], dtype=float)
    return float(np.median(values)) if values.size else None


def _mean(values):
    values = np.asarray([v for v in values if v is not None and np.isfinite(v)], dtype=float)
    return float(values.mean()) if values.size else None


def study1(path: Path, rerun: Path | None):
    cells = json.loads(path.read_text())['cells']
    if rerun is not None:
        again = {c['index']: c for c in (json.loads(p.read_text()) for p in sorted(rerun.glob('[0-9][0-9][0-9].json')))}
    else:
        again = {}
    out = []
    for c in cells:
        s = c['clonal_share']
        row = dict(study=1, index=c['index'], binary=c['binary'], n_groups=c['n_groups'], group_size=c['group_size'],
                   share=c['share'], unbalance=c['unbalance'], law=c['law'], prevalence=c['prevalence'],
                   replicates=c['n_replicates'], coverage=s['coverage'], median_width=s['median_width'],
                   bias=s['bias'], rmse=s['rmse'], estimable_share=s['estimable_share'],
                   rejection=c['n_rejections'] / c['n_replicates'],
                   other_estimators={k: {m: c[k][m] for m in ('bias', 'rmse', 'coverage', 'median_width')}
                                     for k in c if isinstance(c[k], dict) and 'rmse' in c[k] and k != 'clonal_share'})
        if c['index'] in again:
            r = again[c['index']]
            row['rerun'] = dict(replicates=r['n_replicates'], coverage=r['clonal_share']['coverage'],
                                median_width=r['clonal_share']['median_width'],
                                rejection=r['n_rejections'] / r['n_replicates'])
        b = c.get('call_latent_bounds')
        if b:
            row['bounds'] = dict(n=b['n_bounded'], exceptions=b['n_bounded'] - b['sample_share_within_bounds'],
                                 limit_positive=b['limit_positive'] / b['n_bounded'],
                                 mean_limit=b['mean_limit'], mean_truth=b['mean_truth'],
                                 mean_lower=b['mean_lower'], upper_exact=b['upper_exact'] / b['n_bounded'])
        out.append(row)
    return out


def study2(directory: Path, digests: dict):
    out = []
    for design, rows in sorted(_rows(directory, digests).items()):
        rep = [r for r in rows if r['observed_covered'] is not None]
        bounded = [r for r in rows if np.isfinite(r['latent_lower'])]
        out.append(dict(
            study=2, design=design, rho=rows[0]['rho'], replicates=len(rows),
            mean_truth=_mean(r['truth'] for r in rows), mean_estimate=_mean(r['estimate'] for r in rows),
            coverage=(sum(r['observed_covered'] for r in rep) / len(rep)) if rep else None,
            not_reported=1 - len(rep) / len(rows),
            median_width=_median(r['observed_high'] - r['observed_low'] for r in rep),
            rejection=sum(r['p_value'] <= .05 for r in rows) / len(rows),
            rejection_ignoring_laboratory=(sum(r['pooled_p_value'] <= .05 for r in rows) / len(rows)
                                           if 'pooled_p_value' in rows[0] else None),
            mean_end_wells=_mean(r['end_wells'] for r in rows),
            theorem_exceptions=sum(not r['sample_within_bounds'] for r in bounded),
            upper_exact=sum(r['latent_upper_exact'] for r in bounded) / len(bounded) if bounded else None,
            mean_latent_truth=_mean(r['latent_truth'] for r in bounded),
            mean_lower=_mean(r['latent_lower'] for r in bounded),
            limit_above_rho=sum(not r['latent_limit_below'] for r in bounded) / len(bounded) if bounded else None,
            limit_positive=sum(r['latent_limit'] > 0 for r in bounded) / len(bounded) if bounded else None,
            median_limit_over_rho=(_median(r['latent_limit'] / r['latent_truth'] for r in bounded
                                           if r['latent_truth'] > 0) if bounded else None)))
    return out


def study3(directory: Path, digests: dict):
    out = []
    if (directory / 'observed').exists():
        for design, rows in sorted(_rows(directory / 'observed', digests).items()):
            rep = [r for r in rows if r['reported']]
            out.append(dict(study=3, family='observed', design=design, replicates=len(rows), truth=rows[0]['truth'],
                            mean_estimate=_mean(r['estimate'] for r in rows),
                            coverage=sum(bool(r['covered']) for r in rep) / len(rep) if rep else None,
                            not_reported=1 - len(rep) / len(rows),
                            median_width=_median(r['high'] - r['low'] for r in rep),
                            rejection=sum(r['p_value'] <= .05 for r in rows) / len(rows)))
    if (directory / 'latent').exists():
        for design, rows in sorted(_rows(directory / 'latent', digests).items()):
            bounded = [r for r in rows if np.isfinite(r['lower'])]
            with_min = [r for r in bounded if r['limit_below_rho_min'] is not None]
            out.append(dict(study=3, family='latent', design=design, replicates=len(rows), rho=rows[0]['rho'],
                            rho_min=rows[0]['rho_min'], theorem_exceptions=sum(not r['within'] for r in bounded),
                            mean_lower=_mean(r['lower'] for r in bounded), mean_sample=_mean(r['sample'] for r in bounded),
                            upper_exact=sum(r['exact'] for r in bounded) / len(bounded),
                            limit_above_rho=sum(not r['limit_below_rho'] for r in bounded) / len(bounded),
                            limit_above_rho_min=(sum(not r['limit_below_rho_min'] for r in with_min) / len(with_min)
                                                 if with_min else None),
                            limit_positive=sum(r['limit'] > 0 for r in bounded) / len(bounded),
                            median_limit_over_rho=(_median(r['limit'] / r['rho'] for r in bounded)
                                                   if rows[0]['rho'] > 0 else None)))
    if (directory / 'units').exists():
        for design, rows in sorted(_rows(directory / 'units', digests).items()):
            out.append(dict(study=3, family='units', design=design, replicates=len(rows),
                            rejection_within_units=sum(r['p_within_units'] <= .05 for r in rows) / len(rows),
                            rejection_ignoring_units=sum(r['p_ignoring_units'] <= .05 for r in rows) / len(rows)))
    if (directory / 'comparison').exists():
        for design, rows in sorted(_rows(directory / 'comparison', digests).items()):
            for term in rows[0]['truth']:
                rep = [r for r in rows if r['covered'][term] is not None]
                out.append(dict(study=3, family='comparison', design=design, term=term, replicates=len(rows),
                                truth=rows[0]['truth'][term], mean_estimate=_mean(r['estimate'][term] for r in rows),
                                coverage=sum(r['covered'][term] for r in rep) / len(rep) if rep else None,
                                not_reported=1 - len(rep) / len(rows),
                                median_width=_median(r['intervals'][term][1] - r['intervals'][term][0]
                                                     for r in rep)))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--grid', type=Path)
    ap.add_argument('--grid-rerun', type=Path)
    ap.add_argument('--order', type=Path)
    ap.add_argument('--conditional', type=Path, action='append', default=[],
                    help='a directory with observed/, latent/, units/ or comparison/; may be given more than once')
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    digests: dict = {}
    rows = []
    if args.grid:
        rows += study1(args.grid, args.grid_rerun)
    if args.order:
        rows += study2(args.order, digests)
    for directory in args.conditional:
        rows += study3(directory, digests)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'design_summaries.json').write_text(json.dumps(rows, indent=1) + '\n')
    (args.output / 'raw_outputs_sha256.json').write_text(json.dumps(digests, indent=1, sort_keys=True) + '\n')
    print(f'{len(rows)} designs summarised, {len(digests)} raw files digested')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
