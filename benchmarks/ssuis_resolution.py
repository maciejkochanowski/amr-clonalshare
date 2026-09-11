"""The S. suis reading at three definitions of a lineage.

    python benchmarks/ssuis_resolution.py --out benchmarks/results_ssuis_resolution

WHY THIS EXISTS. The article reads the S. suis cohort with the lineage defined
as the source study's hierBAPS population cluster, which is present for every
isolate. A reader may ask whether the ordering of the agents by mechanism is a
property of that label. The cohort also carries a serotype for 541 isolates
and a multilocus sequence type for 458, so the same estimator is run with each
of the three as the lineage, on the whole cohort and on the isolates that hold
all three labels, and the reading the article quotes -- how the agents whose
determinants sit on the chromosome order against those whose determinants
sit on mobile elements, per agent pair and per determinant set -- is
recomputed at each.

The mechanism classes are read from the genetic location the source study
reports per determinant (``ssuis_mechanism.py``) and are used to score the
table, never to fit anything.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from amr_clonalshare.attribution import SUPPORT_THRESHOLD, clonal_share

SEED = 20260901
sys.path.insert(0, str(ROOT / "benchmarks"))
from ssuis_mechanism import DETERMINANTS, load_ssuis, mechanism_class, ordering

#: The class of each agent's determinants in this cohort, read once from the
#: source table: chromosomal, mobile, or mixed where the cohort holds both.
MECHANISM = {agent: mechanism_class(load_ssuis()[0], keys)[0]
             for agent, keys in DETERMINANTS.items()}
LABELS = {"hierbaps_cluster": "baps_cluster", "serotype": "serovar",
          "sequence_type": "mlst"}


def load():
    D = ROOT / "examples/ssuis/data"
    ribo = pd.read_csv(D / "ribo.csv")
    cell = pd.read_csv(D / "cell.csv")
    meta = pd.read_csv(D / "metadata.csv", dtype=str)
    df = (ribo.merge(cell, on="genome_id")
              .merge(meta[["genome_id"] + list(LABELS.values())], on="genome_id"))
    traits = [c for c in ribo.columns if c != "genome_id"] + \
             [c for c in cell.columns if c != "genome_id"]
    return df, traits


def run(df, traits, label, n_boot, n_perm):
    lineage = df[label].to_numpy()
    rows = []
    for j, agent in enumerate(traits):
        y = df[agent].to_numpy(float)
        r = clonal_share(y, lineage, seed=SEED + j, n_boot=n_boot, n_perm=n_perm)
        d = r.as_dict()
        rows.append(dict(agent=agent,
                         mechanism=MECHANISM[agent],
                         kappa_adj=d["kappa_adj"], ci_low=d["ci_low"],
                         ci_high=d["ci_high"], null_mean=d["null_mean"],
                         p_value=d["p_value"], p_floor=d["p_floor"],
                         support=d["support"], estimable=d["estimable"],
                         n=d["n"], n_groups=d["n_groups"],
                         prevalence=d["prevalence"]))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="benchmarks/results_ssuis_resolution")
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    n_boot, n_perm = (100, 100) if a.quick else (400, 200)
    t0 = time.time()
    df, traits = load()
    complete = df.dropna(subset=list(LABELS.values()))
    res = {"provenance": {"seed": SEED, "n_boot": n_boot, "n_perm": n_perm,
                          "support_threshold": SUPPORT_THRESHOLD,
                          "n_isolates": int(len(df)),
                          "n_with_all_three_labels": int(len(complete)),
                          "python": sys.version.split()[0],
                          "numpy": np.__version__, "pandas": pd.__version__},
           "labels": {}}
    ref_col = LABELS["hierbaps_cluster"]
    for name, col in LABELS.items():
        sub = df.dropna(subset=[col])
        whole = run(sub, traits, col, n_boot, n_perm)
        both = run(complete, traits, col, n_boot, n_perm)
        res["labels"][name] = {
            "column": col, "n_isolates": int(len(sub)),
            "n_levels": int(sub[col].nunique()),
            "per_agent": whole, "ordering": ordering(whole),
            "per_agent_common_isolates": both,
            "ordering_common_isolates": ordering(both),
        }
        # A label is present on a different set of isolates from the next, so
        # comparing two labels on their own subsets confounds the label with
        # the isolates it was available for. The reference label is therefore
        # run again on exactly the isolates this label covers: the difference
        # between that arm and this label's own is the label, and the
        # difference between that arm and the reference's full cohort is the
        # loss of isolates.
        if col != ref_col:
            ref_here = run(sub, traits, ref_col, n_boot, n_perm)
            res["labels"][name]["per_agent_reference_same_isolates"] = ref_here
            res["labels"][name]["ordering_reference_same_isolates"] = \
                ordering(ref_here)
            res["labels"][name]["reference_column"] = ref_col
    res["provenance"]["runtime_s"] = round(time.time() - t0, 1)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "ssuis_resolution.json").write_text(
        json.dumps(res, indent=1, default=float) + "\n", encoding="utf-8")
    for name, v in res["labels"].items():
        o = v["ordering"]
        oc = v["ordering_common_isolates"]
        print(f"{name:<18} n={v['n_isolates']:>4} levels={v['n_levels']:>4} "
              f"support={np.median([r['support'] for r in v['per_agent']]):.3f} "
              f"chrom={o['mean_chromosomal']:.3f} non-mobile={o['mean_non_mobile']:.3f} "
              f"mobile={o['mean_mobile']:.3f} "
              f"pairs={o['pairs_non_mobile_higher']}/{o['pairs_non_mobile_vs_mobile']} "
              f"sets={o['set_pairs_non_mobile_higher']}/{o['set_pairs_total']} | "
              f"common n={len(complete)}: pairs={oc['pairs_non_mobile_higher']}/"
              f"{oc['pairs_non_mobile_vs_mobile']}")
        rr = v.get("ordering_reference_same_isolates")
        if rr:
            print(f"{'':<18} reference label on the same {v['n_isolates']} "
                  f"isolates: pairs={rr['pairs_non_mobile_higher']}/"
                  f"{rr['pairs_non_mobile_vs_mobile']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
