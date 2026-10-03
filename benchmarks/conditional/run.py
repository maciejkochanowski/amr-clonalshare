#!/usr/bin/env python3
"""One chunk of replicates of one design of the conditional study.

    python benchmarks/conditional/run.py --family observed --design two_labs \
        --start 0 --replicates 500 --output out/observed/two_labs_00000.json

Every replicate draws the isolates of the design again from its fixed distributions,
from the stream SeedSequence([ROOT_SEED, family, design, replicate]), and
records what the package returns beside the truth written in designs.py. The
package runs with its defaults (5 folds, 20 repeats, 999 bootstrap draws and
999 permutations) unless the family says otherwise.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import designs as D

from amr_clonalshare.comparison import compare_lineage_definitions
from amr_clonalshare.latent_order import order_bounds
from amr_clonalshare.mic_order import mic_order_share

ROOT_SEED = 20261003311
FAMILIES = {'observed': D.observed_designs, 'latent': D.latent_designs,
            'units': D.unit_designs, 'comparison': D.comparison_designs}
TERMS = ['selection_from_a', 'singletons_under_a', 'label_difference_on_scorable',
         'singletons_under_b', 'selection_to_b', 'definition_difference_on_common', 'total_difference']


def _stream(family, design, replicate, part=0):
    f = list(FAMILIES).index(family)
    return np.random.default_rng(np.random.SeedSequence([ROOT_SEED, f, design, replicate, part]))


def _seed(family, design, replicate):
    return int(_stream(family, design, replicate, 1).integers(2 ** 62))


def observed(d, k, rep):
    lo, hi, lin, lab = d.draw(_stream('observed', k, rep))
    r = mic_order_share(lo, hi, lin, strata=lab if d.S > 1 else None, seed=_seed('observed', k, rep))
    ok = bool(r['estimable']) and np.isfinite(r['observed_low'])
    return dict(truth=d.eta, estimate=r['kappa_adj'], low=r['observed_low'], high=r['observed_high'],
                se=r['observed_se'], p_value=r['p_value'], estimable=bool(r['estimable']),
                reported=ok, covered=bool(r['observed_low'] <= d.eta <= r['observed_high']) if ok else None)


def latent(d, k, rep):
    lo, hi, lin, lab, x = d.draw(_stream('latent', k, rep))
    r = order_bounds(lo, hi, lin, strata=lab if d.S > 1 else None, seed=_seed('latent', k, rep))
    sample = D.sample_latent_share(x, lin, lab)
    slack = 1e-9
    return dict(rho=d.rho, rho_min=d.rho_min, sample=sample, lower=r['latent_order_lower'],
                upper=r['latent_order_upper'], upper_bound=r['latent_order_upper_bound'],
                exact=bool(r['latent_order_upper_exact']), sharp=bool(r['latent_order_sharp']),
                limit=r['latent_order_lower_limit'],
                within=bool(r['latent_order_lower'] - slack <= sample <= r['latent_order_upper_bound'] + slack),
                limit_below_rho=bool(r['latent_order_lower_limit'] <= d.rho),
                limit_below_rho_min=(bool(r['latent_order_lower_limit'] <= d.rho_min)
                                     if np.isfinite(d.rho_min) else None))


def units(d, k, rep):
    lo, hi, lin, farm, lab = d.draw(_stream('units', k, rep))
    strata = lab if d.labs > 1 else None
    seed = _seed('units', k, rep)
    within = mic_order_share(lo, hi, lin, strata=strata, units=farm, n_boot=0, seed=seed)
    pooled = mic_order_share(lo, hi, lin, strata=strata, n_boot=0, seed=seed)
    return dict(p_within_units=within['p_value'], p_ignoring_units=pooled['p_value'],
                estimate=within['kappa_adj'], n_strata=within['n_strata'])


def _repeated_both(a, b, mask):
    """Records of mask in a lineage of two or more under both labellings,
    found by removing records until every remaining one has a partner."""
    keep = mask.copy()
    while True:
        drop = np.zeros_like(keep)
        for labels in (a, b):
            names = np.asarray([str(x) for x in labels], dtype=object)
            values, counts = np.unique(names[keep], return_counts=True)
            single = set(values[counts < 2])
            drop |= keep & np.isin(names, list(single))
        if not drop.any():
            return keep
        keep &= ~drop


def _arm(p, labels, mask):
    names = np.asarray([str(x) for x in labels], dtype=object)
    values, code, counts = np.unique(names[mask], return_inverse=True, return_counts=True)
    keep = counts[code] >= 2
    q = p[mask][keep]
    code = np.unique(code[keep], return_inverse=True)[1]
    size = np.bincount(code).astype(float)
    if size.size < 2:
        return float('nan')
    w = size / size.sum()
    mu = np.bincount(code, weights=q) / size
    between = float(w @ (mu - w @ mu) ** 2)
    # the distribution of a lineage is that of a record drawn from it: a merged
    # lineage carries the spread of the rates of its parts within it
    within = float(w @ (mu * (1. - mu)))
    return between / (between + within)


def comparison_truth(d):
    typed_a = np.array([x is not None for x in d.a])
    typed_b = np.array([x is not None for x in d.b])
    common = typed_a & typed_b
    scorable = _repeated_both(d.a, d.b, common)
    arms = np.array([_arm(d.p, d.a, typed_a), _arm(d.p, d.a, common), _arm(d.p, d.a, scorable),
                     _arm(d.p, d.b, scorable), _arm(d.p, d.b, common), _arm(d.p, d.b, typed_b)])
    steps = np.diff(arms)
    return dict(zip(TERMS, np.r_[steps, arms[4] - arms[1], arms[5] - arms[0]].tolist()))


def comparison(d, k, rep):
    y = d.draw(_stream('comparison', k, rep))
    r = compare_lineage_definitions(y, d.a, d.b, ids=d.ids, seed=_seed('comparison', k, rep),
                                    n_perm=199, n_boot=999)
    truth = comparison_truth(d)
    out = dict(truth=truth, estimate={t: r['difference_decomposition'][t] for t in TERMS},
               intervals=r['paired_intervals'])
    out['covered'] = {t: (None if r['paired_intervals'][t] is None else
                          bool(r['paired_intervals'][t][0] - 1e-12 <= truth[t] <= r['paired_intervals'][t][1] + 1e-12))
                      for t in TERMS}
    return out


RUN = {'observed': observed, 'latent': latent, 'units': units, 'comparison': comparison}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--family', choices=list(FAMILIES), required=True)
    ap.add_argument('--design', required=True)
    ap.add_argument('--start', type=int, default=0)
    ap.add_argument('--replicates', type=int, default=10)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    designs = FAMILIES[args.family]()
    names = [d.name for d in designs]
    if args.design not in names:
        ap.error(f'unknown design {args.design!r}; one of {names}')
    if args.output.exists():
        raise FileExistsError(args.output)
    k = names.index(args.design)
    d = designs[k]
    rows = []
    for rep in range(args.start, args.start + args.replicates):
        t0 = time.perf_counter()
        row = dict(family=args.family, design=d.name, replicate=rep, **RUN[args.family](d, k, rep))
        row['seconds'] = time.perf_counter() - t0
        rows.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output.with_suffix('.partial')
    partial.write_text(json.dumps(rows, allow_nan=True, default=float))
    partial.replace(args.output)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
