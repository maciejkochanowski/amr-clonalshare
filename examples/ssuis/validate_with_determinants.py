#!/usr/bin/env python3
"""Check the decomposition's verdict against an independent data layer.

    python examples/ssuis/validate_with_determinants.py
    python examples/ssuis/validate_with_determinants.py --json out.json

The decomposition is computed on laboratory phenotype: an MIC compared with a
cut-off. The source study also reports, for the same isolates, the presence or
absence of 43 resistance determinants called from the genome assemblies. Those
calls come from a different instrument and a different pipeline, and they are
not used anywhere in the phenotype analysis.

That makes them a test rather than an illustration. If the decomposition is
measuring what it claims, then for an agent whose resistance has a single known
determinant the two decompositions should agree: an increase attributed to a
change in lineage composition on the phenotype should be attributed to a change
in lineage composition on the gene, and a within-lineage rate change should
appear as a within-lineage rate change in gene carriage. The components are
estimated from different measurements of the same isolates, so agreement is
evidence and disagreement would be a problem.

Three agents in this panel have a single dominant determinant and are testable:
the MLS_B group through *ermB*, the tetracyclines through *tetO*, and tiamulin
through *lsaE*. The beta-lactams are not testable here, and the reason is
biological rather than a gap in the data: the source study found beta-lactam
resistance in this species to involve many core-genome variants of small
effect appearing in a characteristic order, so there is no single gene column
to decompose.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from amr_clonalshare.clonality import decompose_prevalence_difference

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

#: (phenotype agent, determinant column, printable determinant name). Each
#: pairing is the mechanism the source study reports for that agent.
PAIRS = [
    ("erythromycin", "G__MLSB_ermB_M_O", "ermB"),
    ("tylosin", "G__MLSB_ermB_M_O", "ermB"),
    ("tilmicosin", "G__MLSB_ermB_M_O", "ermB"),
    ("lincomycin", "G__MLSB_ermB_M_O", "ermB"),
    ("tetracycline", "G__Tetra_tetO_M_O", "tetO"),
    ("doxycycline", "G__Tetra_tetO_M_O", "tetO"),
    ("tiamulin", "G__Pleuro_lsaE_M_O", "lsaE"),
]
CONTRASTS = {
    "period": ("UK 2013-2014", "UK 2009-2011"),
    "country": ("Canada", "United Kingdom"),
}


def load() -> pd.DataFrame:
    meta = pd.read_csv(DATA / "metadata.csv")
    source = pd.read_csv(DATA / "source_table_s1.csv", dtype=str)
    determinants = [c for c in source.columns if c.startswith("G__")]
    calls = source.set_index("Strain")[determinants].map(
        lambda v: 1.0 if str(v).strip().lower() in ("yes", "1", "true") else 0.0)
    frame = meta.merge(pd.read_csv(DATA / "ribo.csv"), on="genome_id") \
                .merge(pd.read_csv(DATA / "cell.csv"), on="genome_id") \
                .merge(calls.reset_index().rename(
                    columns={"Strain": "source_strain"}), on="source_strain")
    return frame


def split(frame: pd.DataFrame, contrast: str) -> tuple:
    if contrast == "period":
        uk = frame[(frame["isolation_country"] == "United Kingdom")
                   & frame["collection_year"].notna()]
        return uk[uk["collection_year"] >= 2013], uk[uk["collection_year"] < 2013]
    return (frame[frame["isolation_country"] == "Canada"],
            frame[frame["isolation_country"] == "United Kingdom"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--contrast", choices=sorted(CONTRASTS),
                        default="country")
    parser.add_argument("--lineage", default="baps_cluster")
    parser.add_argument("--n-boot", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--json", type=Path, default=None)
    args = parser.parse_args()

    frame = load()
    side_a, side_b = split(frame, args.contrast)
    name_a, name_b = CONTRASTS[args.contrast]
    print(f"{name_a} (n={len(side_a)}) versus {name_b} (n={len(side_b)}), "
          f"lineage {args.lineage}")
    print(f"{'agent':14s} {'layer':10s} {'variable':10s} "
          f"{'diff':>7s} {'mix':>8s} {'rate':>8s}   {'[rate 95% CI]':>18s}")

    records, rows = {}, []
    for agent, column, gene in PAIRS:
        pair = {}
        for layer, variable in (("phenotype", agent), ("genotype", column)):
            result = decompose_prevalence_difference(
                side_a[variable], side_a[args.lineage],
                side_b[variable], side_b[args.lineage],
                n_boot=args.n_boot, rng=np.random.default_rng(args.seed))
            pair[layer] = result
            label = agent if layer == "phenotype" else gene
            print(f"{agent:14s} {layer:10s} {label:10s} "
                  f"{result['difference'] * 100:+7.1f} "
                  f"{result['composition'] * 100:+8.1f} "
                  f"{result['within_lineage'] * 100:+8.1f}   "
                  f"[{result['within_lineage_ci95'][0] * 100:+7.1f},"
                  f"{result['within_lineage_ci95'][1] * 100:+7.1f}]")
        records[agent] = pair
        rows.append((pair["phenotype"]["within_lineage"],
                     pair["genotype"]["within_lineage"],
                     pair["phenotype"]["composition"],
                     pair["genotype"]["composition"]))

    array = np.asarray(rows)
    within = float(np.corrcoef(array[:, 0], array[:, 1])[0, 1])
    composition = float(np.corrcoef(array[:, 2], array[:, 3])[0, 1])
    largest = float(np.max(np.abs(array[:, 0] - array[:, 1])))
    print(f"\nacross {len(rows)} agent-determinant pairs, phenotype against "
          f"genotype:")
    print(f"  within-lineage components correlate {within:+.3f}")
    print(f"  composition components correlate    {composition:+.3f}")
    print(f"  largest within-lineage disagreement {largest * 100:.1f} points")
    print("  the two layers are different measurements of the same isolates; "
          "agreement is\n  evidence that the split is not an artefact of "
          "either measurement")

    if args.json:
        args.json.write_text(json.dumps(
            {"contrast": args.contrast, "lineage": args.lineage,
             "n_a": int(len(side_a)), "n_b": int(len(side_b)),
             "agreement": {"within_lineage_correlation": within,
                           "composition_correlation": composition,
                           "largest_within_lineage_disagreement": largest},
             "per_agent": records}, indent=2, default=float))
        print(f"\nwritten: {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
