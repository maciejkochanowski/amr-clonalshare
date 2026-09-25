"""Independent calibration of null-wise MIC bootstrap tests and selected intervals."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import time
import numpy as np
import scipy
from amr_clonalshare.censored import censored_clonal_share
from amr_clonalshare._mic_panels import observe_panel
from amr_clonalshare.mic_null_bootstrap import null_bootstrap_test, null_bootstrap_interval
from benchmarks.empirical.mic_calibration import design_grid

ROOT_SEED = 20260916301
#: Seed stream of the calibration; stream 1 is reserved for unit-test datasets.
CALIBRATION_STREAM = 2


def null_design_grid():
    cells = design_grid()
    for groups in (8, 30):
        for pattern in ('uneven', 'singleton_rich'):
            for rho in (.05, .35, .8, .97):
                cells.append(dict(cell=len(cells), groups=groups, size_pattern=pattern,
                    rho=rho, reading='narrow_panel', residual='normal',
                    domain='off_grid_nuisance', mean=.6, total_sd=1.7))
    return cells


def generate(cell, replicate, *, root_seed=ROOT_SEED, stream=CALIBRATION_STREAM):
    rng = np.random.default_rng(np.random.SeedSequence(
        [root_seed, stream, cell['cell'], replicate]))
    g, rho = cell['groups'], cell['rho']
    if cell['size_pattern'] == 'balanced':
        sizes = np.full(g, 20)
    elif cell['size_pattern'] == 'uneven':
        sizes = np.rint(np.exp(rng.normal(2.2, 1., g))).astype(int).clip(2, 100)
    elif cell['size_pattern'] == 'singleton_rich':
        sizes = np.r_[np.ones(g//2, dtype=int), np.full(g-g//2, 20)]
    else:
        sizes = np.ones(g, dtype=int)
    labels = np.repeat(np.arange(g), sizes)
    effects = rng.normal(0, np.sqrt(rho), g)
    if cell['residual'] == 't4':
        residual = rng.standard_t(4, len(labels))/np.sqrt(2.)
    elif cell['residual'] == 'contaminated_normal':
        scales = np.where(rng.random(len(labels)) < .05, 5., 1.)
        residual = rng.normal(size=len(labels))*scales/np.sqrt(2.2)
    else:
        residual = rng.normal(size=len(labels))
    y = cell.get('mean', 0.) + cell.get('total_sd', 1.)*(effects[labels]+np.sqrt(1-rho)*residual)
    edges = np.arange(-4., 5.) if cell['reading'] == 'wide_panel' else np.arange(-1., 2.)
    a, b = observe_panel(y, edges)
    seed = int(np.random.SeedSequence([root_seed, stream, cell['cell'], replicate, 991]).generate_state(1)[0])
    return a, b, labels, edges, seed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cell', type=int, required=True)
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--replicates', type=int, default=25)
    parser.add_argument('--bootstrap', type=int, default=199)
    parser.add_argument('--bounds-every', type=int, default=0)
    parser.add_argument('--root-seed', type=int, default=ROOT_SEED)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.start < 0 or args.replicates < 1 or args.bootstrap < 19 or args.bounds_every < 0 or args.root_seed < 0:
        parser.error('invalid simulation arguments')
    cells = null_design_grid()
    if not 0 <= args.cell < len(cells) or cells[args.cell]['reading'] == 'exact':
        parser.error('choose an interval-reading design cell')
    cell = cells[args.cell]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists():
        raise FileExistsError(args.output)
    rows, stored = [], []
    started = time.perf_counter()
    for rep in range(args.start, args.start+args.replicates):
        tick = time.perf_counter()
        a, b, labels, edges, seed = generate(cell, rep, root_seed=args.root_seed)
        row = dict(cell=cell['cell'], replicate=rep, rho=cell['rho'], root_seed=args.root_seed,
            n=len(labels), groups=cell['groups'], reading=cell['reading'],
            size_pattern=cell['size_pattern'], residual=cell['residual'], domain=cell['domain'])
        row.update(generation_mean=cell.get('mean', 0.), generation_total_sd=cell.get('total_sd', 1.),
            end_fraction=float(np.mean(~np.isfinite(a) | ~np.isfinite(b))),
            moment_low=np.nan, moment_high=np.nan, moment_reportable=False,
            moment_covered=False, moment_estimate=np.nan, interval_low=0., interval_high=1.,
            interval_reportable=False, interval_covered=True, interval_test_covered=True,
            test_reportable=False, pointwise_accepted=True,
            ml_estimate=np.nan, chi_square_covered=False, lr_truth=np.nan,
            critical=np.inf, bootstrap_attempted=0, bootstrap_failed=0,
            profile_failures=0, raw_set_empty=False, failure='', seconds=0.,
            quadrature_error=np.nan, bounds_computed=False,
            bounds_requested=bool(args.bounds_every and rep % args.bounds_every == 0), null_mean=np.nan,
            null_total_sd=np.nan, pvalue=np.nan, anti_conservative_inversion_mismatch=False)
        stats = np.full(args.bootstrap, np.inf)
        try:
            moment = censored_clonal_share(a, b, labels)
            row.update(moment_low=moment.ci_low, moment_high=moment.ci_high,
                moment_reportable=moment.estimable, moment_estimate=moment.kappa,
                moment_covered=bool(moment.estimable and moment.ci_low <= cell['rho'] <= moment.ci_high))
        except Exception as exc:
            row['moment_failure'] = f'{type(exc).__name__}: {exc}'
        try:
            test = null_bootstrap_test(a, b, labels, rho=cell['rho'], panel_edges=edges,
                                      n_boot=args.bootstrap, seed=seed)
            stats = np.array(test.bootstrap_statistics)
            row.update(ml_estimate=test.unrestricted_fit.rho, lr_truth=test.statistic,
                critical=test.critical, pvalue=test.pvalue, interval_reportable=test.reportable,
                interval_covered=test.accepted, interval_test_covered=test.accepted,
                test_reportable=test.reportable, pointwise_accepted=test.accepted,
                bootstrap_attempted=test.bootstrap_attempted, bootstrap_failed=test.bootstrap_failed,
                null_mean=test.null_fit.mean, null_total_sd=test.null_fit.total_sd,
                quadrature_error=max(test.unrestricted_fit.quadrature_error, test.null_fit.quadrature_error))
            row['chi_square_covered'] = bool(test.statistic <= scipy.stats.chi2.ppf(.95, 1))
            if test.reportable:
                row.update(interval_low=np.nan, interval_high=np.nan)
            else:
                row['failure'] = test.status
            if args.bounds_every and rep % args.bounds_every == 0:
                row['bounds_computed'] = True
                interval = null_bootstrap_interval(a, b, labels, panel_edges=edges,
                    n_boot=args.bootstrap, seed=seed)
                covered = interval.low <= cell['rho'] <= interval.high
                row.update(interval_low=interval.low, interval_high=interval.high,
                    interval_covered=covered, interval_reportable=interval.reportable,
                    bounds_computed=True, profile_failures=int(not interval.reportable),
                    anti_conservative_inversion_mismatch=bool(interval.reportable and test.accepted and not covered))
                artifact = args.output.with_name(args.output.stem+f'.bounds_{rep}.npz')
                np.savez_compressed(artifact,
                    rho=np.array([t.rho for t in interval.tests]),
                    pvalue=np.array([t.pvalue for t in interval.tests]),
                    statistic=np.array([t.statistic for t in interval.tests]),
                    critical=np.array([t.critical for t in interval.tests]),
                    statistics=np.array([t.bootstrap_statistics for t in interval.tests]))
                if not interval.reportable:
                    row['failure'] = interval.status
        except Exception as exc:
            row['failure'] = f'{type(exc).__name__}: {exc}'
            row.update(interval_reportable=False, interval_low=0., interval_high=1.,
                       interval_covered=True, interval_test_covered=True)
        row['seconds'] = time.perf_counter()-tick
        rows.append(row)
        stored.append(stats)
    temporary = args.output.with_suffix('.partial')
    keys = sorted(set().union(*(r.keys() for r in rows)))
    with temporary.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    matrix = args.output.with_suffix('.bootstrap.npz')
    matrix_temporary = matrix.with_suffix('.partial.npz')
    np.savez_compressed(matrix_temporary, statistics=np.stack(stored))
    sources = [Path(__file__), Path('benchmarks/empirical/mic_calibration.py'),
               *Path('src/amr_clonalshare').glob('*.py')]
    receipt = dict(arguments=vars(args) | {'output': str(args.output)}, cell=cell,
        source_hashes={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        root_seed=args.root_seed, elapsed_seconds=time.perf_counter()-started,
        python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
        csv_sha256=hashlib.sha256(temporary.read_bytes()).hexdigest(),
        bootstrap_sha256=hashlib.sha256(matrix_temporary.read_bytes()).hexdigest(),
        bounds_artifacts={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in args.output.parent.glob(args.output.stem+'.bounds_*.npz')},
        slurm_job_id=os.getenv('SLURM_JOB_ID'), slurm_array_task_id=os.getenv('SLURM_ARRAY_TASK_ID'))
    receipt_path = args.output.with_suffix('.receipt.json')
    receipt_path.write_text(json.dumps(receipt, indent=2, allow_nan=False))
    matrix_temporary.replace(matrix)
    temporary.replace(args.output)
    print(json.dumps({'cell': args.cell, 'datasets': len(rows),
        'reportable': sum(r['interval_reportable'] for r in rows),
        'covered': sum(r['interval_covered'] for r in rows), 'seconds': receipt['elapsed_seconds']}))


if __name__ == '__main__':
    main()
