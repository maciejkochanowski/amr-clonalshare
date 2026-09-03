#!/usr/bin/env python3
"""Regenerate the deterministic synthetic demo used by CI.

Three cohorts, each a pair of binary layers plus a metadata file:

  planted/   n=240, k_true=3, overlap=0.10  -- development positive control
  planted_confirmation/
              n=240, same frozen design, independent data seed -- confirmation
  null/      n=180, k_true=1                -- no clusters at all
  clonal/    n=200, 25 lineages             -- all structure is clonal

The point of shipping all three is that a demo which only ever shows the
easy positive case cannot demonstrate that the tool declines to report
structure when there is none.

Usage:  python examples/synthetic/make_data.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

from amr_clonalshare.synthetic import synth_cluster_archetypes, synth_lineage_cohort

HERE = Path(__file__).resolve().parent


def _write(sub, amr, vir, meta):
    d = HERE / "data" / sub
    d.mkdir(parents=True, exist_ok=True)
    amr.rename_axis("isolate").to_csv(d / "amr.csv")
    vir.rename_axis("isolate").to_csv(d / "vir.csv")
    meta.rename_axis("isolate").to_csv(d / "metadata.csv")
    print(f"  {sub}: amr{amr.shape} vir{vir.shape} meta{meta.shape}")


def main():
    print("writing synthetic demo cohorts:")
    amr, vir, labels = synth_cluster_archetypes(
        n=240, p_amr=24, p_vir=36, k_true=3, overlap=0.10, seed=17)
    lineages = np.repeat(np.arange(60), 4)
    np.random.default_rng(1701).shuffle(lineages)
    meta = pd.DataFrame({"planted_cluster": labels,
                         "lineage": [f"L{x:03d}" for x in lineages]},
                        index=amr.index)
    _write("planted", amr, vir, meta)

    # The design above was chosen on seed 17. Seed 117 is a locked independent
    # realization: it must pass unchanged or the positive-control claim fails.
    amr, vir, labels = synth_cluster_archetypes(
        n=240, p_amr=24, p_vir=36, k_true=3, overlap=0.10, seed=117)
    lineages = np.repeat(np.arange(60), 4)
    np.random.default_rng(11701).shuffle(lineages)
    meta = pd.DataFrame({"planted_cluster": labels,
                         "lineage": [f"L{x:03d}" for x in lineages]},
                        index=amr.index)
    _write("planted_confirmation", amr, vir, meta)

    amr, vir, labels = synth_cluster_archetypes(
        n=180, p_amr=18, p_vir=30, k_true=1, overlap=0.25, seed=8)
    meta = pd.DataFrame({"planted_cluster": labels,
                         "lineage": [f"L{i % 45:02d}" for i in range(len(labels))]},
                        index=amr.index)
    _write("null", amr, vir, meta)

    amr, vir, meta = synth_lineage_cohort(
        n=200, n_lineages=25, p_amr=15, p_vir=20, seed=9)
    _write("clonal", amr, vir, meta)


if __name__ == "__main__":
    main()
