#!/usr/bin/env python3
"""Experiment G again, with the thresholds compared on the same data.

The shipped experiment G threads one random stream through every threshold in
turn, so replicate *i* at threshold 0.5 and replicate *i* at threshold 0.7
are different cohorts. The cells are still each valid, but their difference
carries sampling noise that has nothing to do with the threshold, and two
runs of different size disagreed by more than that noise should allow.

Here every replicate is one cohort, drawn once from its own seed, and every
threshold is applied to that cohort with the same split stream. A difference
between two thresholds is then a difference between two rules on identical
data, and the paired counts say how often the rules disagree.

Seeding follows NumPy's SeedSequence: the root seed is spawned into one child
per replicate, and each child is spawned into a data stream and a split
stream. A replicate can be recomputed alone from its index, and the set of
replicates does not depend on how the run was cut into array tasks.

    python benchmarks/block_threshold_crn.py --start 0 --count 2000 --out DIR
    python benchmarks/block_threshold_crn.py --aggregate DIR
"""
from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import platform
import sys
from datetime import date
from pathlib import Path

import numpy as np
import scipy
import sklearn

sys.path.insert(0, str(Path(__file__).resolve().parent))
from calibration_study import _binary_null, _km

from amr_clonalshare.clonality import _clopper_pearson
from amr_clonalshare.inference import correlation_blocks, feature_split_test

ROOT_SEED = 20260904
THRESHOLDS = (0.0, 0.3, 0.5, 0.7, 0.85, 0.9)
N, P, BLOCK, N_SPLITS = 200, 30, 5, 11
ALPHA = 0.05


def one_replicate(index: int) -> dict:
    child = np.random.SeedSequence(ROOT_SEED).spawn(index + 1)[index]
    data_ss, split_ss = child.spawn(2)
    X = _binary_null(np.random.default_rng(data_ss), N, P, block_size=BLOCK)
    out = {"index": index}
    for t in THRESHOLDS:
        g = correlation_blocks(X, threshold=t) if t > 0 else None
        # A fresh generator from the same seed for every threshold: the split
        # draws are common to the thresholds, so what differs is the rule.
        r = feature_split_test(X, 2, cluster_fn=_km,
                               rng=np.random.default_rng(split_ss),
                               n_splits=N_SPLITS, groups=g)
        out["%.2f" % t] = {
            "status": r.get("status"),
            "p": r.get("p_value"),
            "units": r.get("n_split_units"),
            "cross": r.get("max_abs_corr_between_units"),
        }
    return out


def run_chunk(start: int, count: int, out: Path, procs: int = 1) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    indices = range(start, start + count)
    if procs > 1:
        from multiprocessing import Pool
        with Pool(procs) as pool:
            rows = pool.map(one_replicate, indices, chunksize=4)
    else:
        rows = [one_replicate(i) for i in indices]
    path = out / ("part_%06d_%06d.json" % (start, start + count))
    path.write_text(json.dumps(rows, separators=(",", ":")) + "\n")
    return path


def aggregate(out: Path) -> dict:
    rows = []
    for part in sorted(out.glob("part_*.json")):
        rows.extend(json.loads(part.read_text()))
    rows.sort(key=lambda r: r["index"])
    indices = [r["index"] for r in rows]
    assert indices == list(range(len(rows))), "replicates missing or duplicated"
    n = len(rows)
    keys = ["%.2f" % t for t in THRESHOLDS]

    cells = {}
    reject = {}
    for k in keys:
        ok = [r for r in rows if r[k]["status"] == "ok"]
        rej = np.array([r[k]["p"] <= ALPHA for r in ok])
        reject[k] = {r["index"]: bool(r[k]["p"] <= ALPHA) for r in ok}
        lo, hi = _clopper_pearson(int(rej.sum()), rej.size)
        cells[k] = {
            "n_usable": int(rej.size), "n_not_ok": n - int(rej.size),
            "rejections": int(rej.sum()), "type_I": float(rej.mean()),
            "cp95_low": lo, "cp95_high": hi,
            "holds_level": hi <= ALPHA,
            "mean_units": float(np.mean([r[k]["units"] for r in ok])),
            "mean_cross": float(np.mean([r[k]["cross"] for r in ok])),
        }

    # Paired contrasts: the same cohorts under two rules.
    pairs = {}
    for a, b in (("0.50", "0.70"), ("0.70", "0.85"), ("0.30", "0.50"),
                 ("0.85", "0.90")):
        common = sorted(set(reject[a]) & set(reject[b]))
        ra = np.array([reject[a][i] for i in common])
        rb = np.array([reject[b][i] for i in common])
        both, only_a, only_b = int((ra & rb).sum()), int((ra & ~rb).sum()), int((~ra & rb).sum())
        diff = float(rb.mean() - ra.mean())
        # Exact interval for the paired difference from the discordant
        # counts: with d = only_a + only_b discordant pairs the difference is
        # (only_b - only_a) / m, and only_b | d is binomial.
        d = only_a + only_b
        if d:
            lo_b, hi_b = _clopper_pearson(only_b, d)
            dlo, dhi = (2 * lo_b - 1) * d / len(common), (2 * hi_b - 1) * d / len(common)
        else:
            dlo = dhi = 0.0
        pairs["%s_minus_%s" % (b, a)] = {
            "n_pairs": len(common), "both_reject": both,
            "only_%s" % a: only_a, "only_%s" % b: only_b,
            "difference": diff, "ci95_low": dlo, "ci95_high": dhi,
        }

    src = inspect.getsource(_binary_null).encode()
    return {
        "what": "type-I error of the feature-split test by automatic blocking "
                "threshold, every threshold on the same cohorts",
        "date": date.today().isoformat(),
        "design": {"n_replicates": n, "n": N, "p": P, "true_block_size": BLOCK,
                   "n_splits": N_SPLITS, "thresholds": list(THRESHOLDS),
                   "alpha": ALPHA, "root_seed": ROOT_SEED,
                   "seeding": "SeedSequence(root).spawn(n)[i] -> (data, split); "
                              "the split stream is re-seeded identically for "
                              "each threshold",
                   "generator": "benchmarks/calibration_study.py::_binary_null",
                   "generator_source_sha256": hashlib.sha256(src).hexdigest()},
        "cells": cells,
        "paired": pairs,
        "environment": {"python": platform.python_version(),
                        "numpy": np.__version__, "scipy": scipy.__version__,
                        "scikit_learn": sklearn.__version__},
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--start", type=int)
    ap.add_argument("--count", type=int)
    ap.add_argument("--aggregate", action="store_true")
    ap.add_argument("--procs", type=int, default=1,
                    help="worker processes; each replicate is one task")
    args = ap.parse_args()
    out = Path(args.out)
    if args.aggregate:
        result = aggregate(out)
        path = out / "block_threshold_crn.json"
        text = json.dumps(result, indent=2)
        path.write_text(text + "\n")
        print("wrote %s sha256 %s" % (path, hashlib.sha256(text.encode()).hexdigest()[:12]))
        for k, v in result["cells"].items():
            print("  threshold %s  type-I %.4f  [%.4f, %.4f]  units %.2f" % (
                k, v["type_I"], v["cp95_low"], v["cp95_high"], v["mean_units"]))
        for k, v in result["paired"].items():
            print("  %s  %+.4f  [%+.4f, %+.4f]" % (k, v["difference"], v["ci95_low"], v["ci95_high"]))
        return 0
    if args.start is None or args.count is None:
        ap.error("--start and --count are required unless --aggregate")
    path = run_chunk(args.start, args.count, out, procs=args.procs)
    print("wrote", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
