#!/usr/bin/env python3
"""Size of the ArchePhy-CS null on traits that genuinely evolved independently.

The exact inclusion-exclusion null is exact -- for a null in which each trait's
Fitch change edges are a uniform random ``c_t``-subset of the ``E`` edges. Real
change edges are not uniform: they concentrate on short and terminal branches,
so two traits that evolved **independently** on the same tree share change edges
more often than uniform subsets do. The exact null is therefore anticonservative,
and "no resampling and no asymptotics means calibration at small n by
construction" does not follow.

This script measures it. Traits are simulated under a symmetric two-state CTMC
on simulated birth-death trees -- i.e. the null hypothesis is *true by
construction* -- and the rejection rate of each method is the measured size.

    python benchmarks/archephy_calibration.py --quick
    python benchmarks/archephy_calibration.py
"""
from __future__ import annotations

import argparse
import json
import pathlib
import random
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from amr_clonalshare.archephy import (  # noqa: E402
    _ctmc_leaf_states,
    archephy_cs_test,
    load_tree,
)

METHODS = ("ctmc", "exact", "poisson", "weighted")


def _tree(n_taxa: int, seed: int):
    import dendropy
    t = dendropy.simulate.treesim.birth_death_tree(
        birth_rate=1.0, death_rate=0.4, num_extant_tips=n_taxa,
        rng=random.Random(seed))
    return t.as_string(schema="newick").strip()


def _independent_traits(newick: str, rate: float, k: int, rng,
                        clock: str = "strict"):
    """k traits evolved independently on the same tree, plus their leaf labels.

    ``clock="relaxed"`` multiplies every branch length by a shared per-edge
    Gamma(0.3) draw before evolving the traits. This is the honest adversarial
    case for the CTMC null: the traits are still independent, but the *rate*
    varies along the tree in a way the null does not model, so both traits are
    pushed towards changing on the same fast edges. A null calibrated only
    against data from its own generative model would miss it -- which is the
    criticism this arm exists to answer, and it is the same criticism that broke
    the one-factor continuum null in round 4.
    """
    import dendropy
    tree, leaves = load_tree(newick)
    if clock == "relaxed":
        t2 = dendropy.Tree.get(data=newick, schema="newick")
        mult = {}
        for nd in t2.preorder_node_iter():
            if nd.edge is not None and nd.edge.length:
                mult[id(nd)] = float(rng.gamma(0.3, 1 / 0.3))
        for nd, tgt in zip(t2.preorder_node_iter(), tree.preorder_node_iter()):
            if tgt.edge is not None and tgt.edge.length and id(nd) in mult:
                tgt.edge.length = max(tgt.edge.length * mult[id(nd)], 1e-6)
    cols = [_ctmc_leaf_states(tree, leaves, rate, rng) for _ in range(k)]
    return np.column_stack(cols), leaves


def size_study(n_sim: int, n_taxa: int, rate: float, k: int, n_mc: int,
               alpha: float = 0.05, n_sim_cheap: int = 800,
               clock: str = "strict") -> dict:
    """Size of each null. The closed-form methods are free, so they get many more
    replicates: distinguishing 0.077 from 0.050 needs a standard error well
    under 0.01, which 60 replicates (SE 0.028) cannot deliver."""
    # Separate streams: the data must not depend on how much randomness the
    # methods happen to consume, or the cells are not comparable across runs.
    data_rng = np.random.default_rng(11)
    hits = {m: 0 for m in METHODS}
    usable = {m: 0 for m in METHODS}
    m_obs_all = []
    for b in range(max(n_sim, n_sim_cheap)):
        newick = _tree(n_taxa, seed=1000 + b)
        X, leaves = _independent_traits(newick, rate, k, data_rng, clock=clock)
        if X.sum(axis=0).min() == 0 or X.sum(axis=0).max() == len(leaves):
            continue                      # a constant trait has no change edges
        for m in METHODS:
            cheap = m in ("exact", "poisson")
            if b >= (n_sim_cheap if cheap else n_sim):
                continue
            try:
                r = archephy_cs_test(newick, X, leaves, list(range(k)),
                                     method=m, n_mc=n_mc,
                                     rng=np.random.default_rng(7000 + b),
                                     alpha=alpha)
            except Exception:
                continue
            usable[m] += 1
            if r["p_value"] <= alpha:
                hits[m] += 1
            if m == "ctmc":
                m_obs_all.append(r["m_obs"])
    return {
        "n_taxa": n_taxa, "rate": rate, "k": k, "n_mc": n_mc, "clock": clock,
        "n_usable": usable,
        "size_at_0.05": {m: (hits[m] / usable[m]) if usable[m] else None
                         for m in METHODS},
        "standard_error": {m: (float(np.sqrt((hits[m] / usable[m]) *
                                             (1 - hits[m] / usable[m]) / usable[m]))
                               if usable[m] else None) for m in METHODS},
        "mean_m_obs": float(np.mean(m_obs_all)) if m_obs_all else None,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    n_sim = 60 if args.quick else 250
    n_mc = 200 if args.quick else 500
    n_cheap = 600 if args.quick else 3000

    # The last cell is the misspecification arm: traits still independent, but a
    # shared per-edge rate multiplier the CTMC null does not model.
    cells = [(60, 3.0, 2, "strict"), (60, 1.0, 2, "strict"),
             (30, 3.0, 2, "strict"), (60, 3.0, 2, "relaxed")]
    out = {"description": "size of each ArchePhy-CS null on traits simulated "
                          "under independent two-state CTMCs on birth-death "
                          "trees, i.e. with the null true by construction",
           "n_sim": n_sim, "alpha": 0.05, "cells": {}}
    for n_taxa, rate, k, clock in cells:
        key = f"n_taxa={n_taxa}/rate={rate}/k={k}/clock={clock}"
        print(key, "...", flush=True)
        out["cells"][key] = size_study(n_sim, n_taxa, rate, k, n_mc,
                                       n_sim_cheap=n_cheap, clock=clock)
        print("  ", out["cells"][key]["size_at_0.05"], flush=True)

    path = pathlib.Path(args.out or (pathlib.Path(__file__).parent / "results"
                                     / "archephy_size.json"))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {path}")
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
