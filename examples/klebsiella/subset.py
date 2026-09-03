#!/usr/bin/env python3
"""Build the case-study layers from a full Kleborate output table.

This script is the specification of how `data/` should be derived. It does not
replay the history of the committed files: the selection that produced them was
not recorded (see DATA_PROVENANCE.md). Run it against a Kleborate table to
regenerate an equivalent, fully documented cohort.

    python examples/klebsiella/subset.py --kleborate FULL.tsv \
        --n 1500 --seed 42 --out examples/klebsiella/data/

Selection is a **stratified de-replication**: at most `--max-per-st` isolates
per sequence type, sampled without replacement under `--seed`, then a random
sample of the survivors down to `--n`. De-replicating by lineage before
sampling is deliberate. A surveillance collection samples outbreak clones far
more heavily than the population they came from, and an unweighted sample from
it makes the number of discovered archetypes a function of submission
behaviour: on this cohort, running the pipeline before and after one-per-ST
de-replication changes the selected k. Capping per-ST representation does not
remove that dependence, so the case study still reports the lineage
diagnostics; it makes the starting point defensible rather than silent.

Column conventions produced:
  amr.csv    Strain_ID + one column per acquired AMR class / beta-lactamase subclass
  vir.csv    Strain_ID + one column per virulence gene
  kloc.csv   Strain_ID + one one-hot column per K locus  (K_<locus>)
  otype.csv  Strain_ID + one one-hot column per O antigen (O_<type>)
  metadata.csv  Strain_ID + ST, Country, Region, Year, Source, K_locus, O_type,
                virulence_score, resistance_score, Convergence_event
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

AMR_CLASSES = ["AGly", "Col", "Fcyn", "Flq", "Gly", "MLS", "Phe", "Rif", "Sul",
               "Tet", "Tgc", "Tmt", "Bla", "Bla_inhR", "Bla_ESBL",
               "Bla_ESBL_inhR", "Bla_Carb"]
VIR_LOCI = {
    "ybt": ["ybtS", "ybtX", "ybtQ", "ybtP", "ybtA", "irp2", "irp1", "ybtU",
            "ybtT", "ybtE", "fyuA"],
    "clb": ["clbA", "clbB", "clbC", "clbD", "clbE", "clbF", "clbG", "clbH",
            "clbI", "clbL", "clbM", "clbN", "clbO", "clbP", "clbQ"],
    "iuc": ["iucA", "iucB", "iucC", "iucD", "iutA"],
    "iro": ["iroB", "iroC", "iroD", "iroN"],
    "rmp": ["rmpA", "rmpD", "rmpC"],
}
META_COLS = ["ST", "Country", "Region", "Year", "Source", "K_locus", "O_type",
             "virulence_score", "resistance_score", "Convergence_event"]


def stratified_dereplicate(df: pd.DataFrame, st_col: str, max_per_st: int,
                           n: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    keep = []
    for _, grp in df.groupby(st_col, sort=True):
        idx = grp.index.to_numpy()
        if idx.size > max_per_st:
            idx = rng.choice(idx, size=max_per_st, replace=False)
        keep.extend(idx.tolist())
    keep = np.sort(np.asarray(keep))
    if keep.size > n:
        keep = np.sort(rng.choice(keep, size=n, replace=False))
    return df.loc[keep]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kleborate", required=True,
                    help="Kleborate output table (TSV) for the full cohort")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--n", type=int, default=1500)
    ap.add_argument("--max-per-st", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--strain-column", default="strain")
    args = ap.parse_args()

    src = pd.read_csv(args.kleborate, sep="\t")
    if args.strain_column not in src.columns:
        raise SystemExit(f"column {args.strain_column!r} not in {args.kleborate}")
    src = src.rename(columns={args.strain_column: "Strain_ID"}).reset_index(drop=True)
    st_col = "ST" if "ST" in src.columns else None
    if st_col is None:
        raise SystemExit("Kleborate table must carry an 'ST' column")

    sel = stratified_dereplicate(src, st_col, args.max_per_st, args.n, args.seed)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    amr = pd.DataFrame({"Strain_ID": sel["Strain_ID"].to_numpy()})
    for c in AMR_CLASSES:
        col = sel[c] if c in sel.columns else pd.Series("-", index=sel.index)
        amr[c] = (~col.astype(str).isin(["-", "", "nan", "NA"])).astype(int)
    amr.to_csv(out / "amr.csv", index=False)

    vir = pd.DataFrame({"Strain_ID": sel["Strain_ID"].to_numpy()})
    for locus, genes in VIR_LOCI.items():
        for g in genes:
            col = sel[g] if g in sel.columns else pd.Series("-", index=sel.index)
            vir[g] = (~col.astype(str).isin(["-", "", "nan", "NA"])).astype(int)
    vir.to_csv(out / "vir.csv", index=False)

    for src_col, prefix, fname in (("K_locus", "K_", "kloc.csv"),
                                   ("O_type", "O_", "otype.csv")):
        if src_col not in sel.columns:
            continue
        dummies = pd.get_dummies(sel[src_col].astype(str), prefix=prefix.rstrip("_"),
                                 prefix_sep="_").astype(int)
        dummies.insert(0, "Strain_ID", sel["Strain_ID"].to_numpy())
        dummies.to_csv(out / fname, index=False)

    meta = pd.DataFrame({"Strain_ID": sel["Strain_ID"].to_numpy()})
    for c in META_COLS:
        meta[c] = sel[c].to_numpy() if c in sel.columns else ""
    meta.to_csv(out / "metadata.csv", index=False)

    print(f"wrote {len(sel)} isolates to {out} "
          f"({sel[st_col].nunique()} sequence types, "
          f"max {sel[st_col].value_counts().max()} per ST)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
