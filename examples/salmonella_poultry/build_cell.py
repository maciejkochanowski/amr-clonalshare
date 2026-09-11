#!/usr/bin/env python3
"""Cut the poultry-meat *Salmonella* cell the article reads out of one NCBI
Pathogen Detection release, and write it in the two tables the run takes.

    python examples/salmonella_poultry/build_cell.py --raw <dir>

``<dir>`` holds ``Salmonella.ast.tsv`` and ``Salmonella.clusters.tsv`` as
``benchmarks/fetch_pathogen_detection.sh`` writes them for the release
accession pinned in ``benchmarks/pathogen_detection_releases.tsv``
(PDG000000002.4210). The cell is every isolate the release places in a SNP
cluster, that the host and sampling-matrix rules of
``benchmarks/vet_source_taxonomy.py`` assign to poultry and to a food matrix,
and that carries a serovar and a collection date: 7,049 isolates. The panel is
every susceptibility call the release records for them, one row per isolate
and antimicrobial, with S written as susceptible and I or R as
non-susceptible, the reading NCBI applies to its own calls. Nothing is
imputed and no isolate is dropped by a decision made here; an agent an
isolate was not tested against has no row.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
sys.path.insert(0, str(HERE.parent.parent / "benchmarks"))
from verify_vet_claims import parse_ast, serovar, year  # noqa: E402
from vet_source_taxonomy import classify  # noqa: E402

RELEASE = "PDG000000002.4210"
CALL = {0.0: "susceptible", 1.0: "non-susceptible"}


def build(raw: Path):
    ast = pd.read_csv(raw / "Salmonella.ast.tsv", sep="\t", dtype=str,
                      on_bad_lines="skip")
    clusters = pd.read_csv(raw / "Salmonella.clusters.tsv", sep="\t",
                           dtype=str, on_bad_lines="skip")
    df = ast.merge(clusters[["target_acc", "PDS_acc"]], on="target_acc")
    df = df[df["PDS_acc"].notna() & (df["PDS_acc"] != "NULL")]
    lab = [classify(h, s) for h, s in zip(df["host"], df["isolation_source"])]
    df = df.assign(host_group=[c["group"] for c in lab],
                   matrix=[c["matrix"] for c in lab],
                   serovar=[serovar(n) for n in df["scientific_name"]],
                   collection_year=[year(v) for v in df["collection_date"]])
    cell = df[(df["host_group"] == "poultry") & (df["matrix"] == "animal_food")
              & (df["serovar"] != "") & df["collection_year"].notna()]
    cell = cell.sort_values("target_acc").reset_index(drop=True)
    metadata = pd.DataFrame({
        "genome_id": cell["target_acc"],
        "serovar": cell["serovar"],
        "pds_cluster": cell["PDS_acc"],
        "collection_year": cell["collection_year"].astype(int),
        "period": ["from_2016" if y >= 2016 else "before_2016"
                   for y in cell["collection_year"]],
        "country": cell["geo_loc_name"].fillna("").str.split(":").str[0],
        "isolation_source": cell["isolation_source"].fillna(""),
    })
    rows = []
    for acc, field in zip(cell["target_acc"], cell["AST_phenotypes"]):
        for agent, call in sorted(parse_ast(field).items()):
            rows.append((acc, agent, CALL[call]))
    calls = pd.DataFrame(rows, columns=["genome_id", "antibiotic", "call"])
    return metadata, calls


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--raw", required=True, type=Path)
    a = ap.parse_args()
    metadata, calls = build(a.raw)
    DATA.mkdir(exist_ok=True)
    metadata.to_csv(DATA / "metadata.csv", index=False)
    calls.to_csv(DATA / "calls_long.csv", index=False)
    receipt = {
        "release": RELEASE,
        "n_isolates": int(len(metadata)),
        "n_serovars": int(metadata["serovar"].nunique()),
        "n_clusters": int(metadata["pds_cluster"].nunique()),
        "years": [int(metadata["collection_year"].min()),
                  int(metadata["collection_year"].max())],
        "n_calls": int(len(calls)),
        "n_antimicrobials": int(calls["antibiotic"].nunique()),
        "cell": "host poultry, matrix animal_food, serovar and collection "
                "date present, SNP cluster assigned",
    }
    (DATA / "cell_receipt.json").write_text(json.dumps(receipt, indent=1) + "\n")
    print(json.dumps(receipt, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
