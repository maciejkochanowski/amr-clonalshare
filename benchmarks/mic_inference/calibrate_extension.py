"""Extension of the null-wise calibration to heavy censoring, within-lineage mixtures and laboratory effects.

The designs below were added after the main campaign to cover three features of
real MIC collections that the main designs do not: dilution panels on which most
readings fall on an end well, a two-component residual distribution within
lineages, and a laboratory effect that is partly aligned with lineage. The driver
reuses the loop, the outputs and the receipts of ``calibrate_null.py``; only the
design grid and the data-generating function differ. Each design is a separate
campaign root with its own root seed, so the main results are untouched.

The laboratory designs come in three analyses: pooled without a laboratory
term, restricted to one laboratory, and pooled with the laboratory as a fixed
effect of the MIC model (``covariate``); the third was added after the first
two had shown that a pooled analysis without the term is centred on the
between-lineage share of the mixed collection rather than on the generating
value. The last three designs add a second covariate, the year of isolation,
and adjust for both.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np
from amr_clonalshare._mic_panels import observe_panel
import benchmarks.mic_inference.calibrate_null as base
import amr_clonalshare.mic_null_bootstrap as _null

ROOT_SEED = 20260917401
END_SHIFT = {'end_0.6': -1.3, 'end_0.75': -1.7}   # mean of y for a unit total SD and cut points -1, 0, 1
MIXTURE_SHARE, MIXTURE_SHIFT = .3, 2.5
LAB_EFFECT, LAB_ALIGNMENT = .5, .85
YEAR_DRIFT = .3   # shift per year of isolation, three years, independent of lineage and laboratory


def extension_design_grid():
    cells = []
    for level in ('end_0.6', 'end_0.75'):
        for groups in (8, 30):
            for pattern in ('uneven', 'singleton_rich'):
                for rho in (.1, .5, .9):
                    cells.append(dict(cell=len(cells), groups=groups, size_pattern=pattern, rho=rho,
                                      reading='narrow_panel', residual='normal', domain='heavy_censoring',
                                      mean=END_SHIFT[level], total_sd=1., level=level))
    for reading in ('wide_panel', 'narrow_panel'):
        for rho in (.1, .5, .9):
            cells.append(dict(cell=len(cells), groups=30, size_pattern='uneven', rho=rho, reading=reading,
                              residual='two_component', domain='bimodal_within_lineage'))
    for analysis in ('pooled', 'stratified', 'adjusted'):
        for reading in ('wide_panel', 'narrow_panel'):
            for rho in (.1, .5, .9):
                cells.append(dict(cell=len(cells), groups=30, size_pattern='uneven', rho=rho, reading=reading,
                                  residual='normal', domain='laboratory_' + analysis, laboratory=analysis))
    for rho in (.1, .5, .9):
        cells.append(dict(cell=len(cells), groups=30, size_pattern='uneven', rho=rho, reading='narrow_panel',
                          residual='normal', domain='country_year_adjusted', laboratory='adjusted', year=True))
    return cells


# The laboratory of every isolate of the dataset generated last; the loop of
# calibrate_null.py passes intervals and lineages only, so the wrappers below
# add it as the covariate of the adjusted analysis.
_CURRENT = {'covariate': None}


def _test_with_covariate(*args, **kwargs):
    return _null.null_bootstrap_test(*args, covariate=_CURRENT['covariate'], **kwargs)


def _interval_with_covariate(*args, **kwargs):
    return _null.null_bootstrap_interval(*args, covariate=_CURRENT['covariate'], **kwargs)


def generate(cell, replicate, *, root_seed=ROOT_SEED, stream=base.CALIBRATION_STREAM):
    rng = np.random.default_rng(np.random.SeedSequence([root_seed, stream, cell['cell'], replicate]))
    g, rho = cell['groups'], cell['rho']
    if cell['size_pattern'] == 'uneven':
        sizes = np.rint(np.exp(rng.normal(2.2, 1., g))).astype(int).clip(2, 100)
    elif cell['size_pattern'] == 'singleton_rich':
        sizes = np.r_[np.ones(g//2, dtype=int), np.full(g-g//2, 20)]
    else:
        sizes = np.full(g, 20)
    labels = np.repeat(np.arange(g), sizes)
    effects = rng.normal(0, np.sqrt(rho), g)
    if cell['residual'] == 'two_component':
        shifted = rng.random(len(labels)) < MIXTURE_SHARE
        raw = rng.normal(size=len(labels)) + MIXTURE_SHIFT*shifted
        residual = (raw - MIXTURE_SHARE*MIXTURE_SHIFT)/np.sqrt(1. + MIXTURE_SHARE*(1-MIXTURE_SHARE)*MIXTURE_SHIFT**2)
    else:
        residual = rng.normal(size=len(labels))
    y = cell.get('mean', 0.) + cell.get('total_sd', 1.)*(effects[labels] + np.sqrt(1-rho)*residual)
    if 'laboratory' in cell:
        # Isolates of even-numbered lineages are tested mostly in laboratory 1,
        # odd-numbered mostly in laboratory 0, and the laboratories differ by a
        # fixed shift that the lineage-only model cannot see.
        aligned = rng.random(len(labels)) < LAB_ALIGNMENT
        laboratory = np.where(labels % 2 == 0, aligned, ~aligned).astype(int)
        y = y + LAB_EFFECT*(2*laboratory - 1)
        if cell['laboratory'] == 'stratified':
            keep = laboratory == 1
            y, labels = y[keep], labels[keep]
        if cell.get('year'):
            # A second covariate, the year of isolation, drawn independently of
            # lineage and laboratory, with a drift of YEAR_DRIFT per year; the
            # adjusted analysis carries both covariates as fixed effects.
            year = rng.integers(0, 3, len(labels))
            y = y + YEAR_DRIFT*year
            laboratory = np.column_stack([laboratory, year])
        _CURRENT['covariate'] = laboratory if cell['laboratory'] == 'adjusted' else None
    else:
        _CURRENT['covariate'] = None
    edges = np.arange(-4., 5.) if cell['reading'] == 'wide_panel' else np.arange(-1., 2.)
    a, b = observe_panel(y, edges)
    seed = int(np.random.SeedSequence([root_seed, stream, cell['cell'], replicate, 991]).generate_state(1)[0])
    return a, b, labels, edges, seed


def main():
    base.null_design_grid = extension_design_grid
    base.generate = generate
    base.ROOT_SEED = ROOT_SEED
    base.null_bootstrap_test = _test_with_covariate
    base.null_bootstrap_interval = _interval_with_covariate
    base.main()
    # The receipt written by the shared loop names calibrate_null.py; record this
    # driver beside it so the extension results are traceable to both files.
    import sys
    output = Path(sys.argv[sys.argv.index('--output') + 1])
    receipt_path = output.with_suffix('.receipt.json')
    receipt = json.loads(receipt_path.read_text())
    receipt['extension_driver_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    receipt['extension_root_seed'] = ROOT_SEED
    receipt_path.write_text(json.dumps(receipt, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
