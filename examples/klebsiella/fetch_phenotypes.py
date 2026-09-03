#!/usr/bin/env python3
"""Fetch measured AST phenotypes for the case-study isolates from BV-BRC.

Why this script exists
----------------------
Every "external" criterion the case study previously scored its partition
against -- Kleborate's ``virulence_score`` and ``resistance_score`` -- is a
deterministic function of the very matrices being clustered (reconstructible at
99.7 % and 97.8 %; see DATA_PROVENANCE.md). Scoring a partition against them
measures internal consistency, not external agreement. The cohort had **no
phenotypic anchor at all**, which is a larger validation gap than any missing
benchmark against another integration method.

It does not have to stay that way. 1074 of the 1500 ``Strain_ID`` values are
assembly accessions, and BV-BRC (the successor to PATRIC) holds
laboratory-measured antimicrobial susceptibility results keyed to assemblies:

> Olson RD, Assaf R, Brettin T, et al. (2023). *Introducing the Bacterial and
> Viral Bioinformatics Resource Center (BV-BRC): a resource combining PATRIC,
> IRD and ViPR.* Nucleic Acids Research 51(D1):D678-D689.
> doi:10.1093/nar/gkac1003

Two things this script is careful about.

* **Laboratory evidence only.** ``genome_amr`` mixes measured results with
  AdaBoost predictions (``evidence = "Computational Method"``). A predicted
  phenotype is a function of the genome, so scoring against it would repeat
  exactly the circularity we are trying to escape. Only
  ``evidence = "Laboratory Method"`` records are kept, and the filter is
  asserted rather than assumed.
* **RefSeq/GenBank pairing.** BV-BRC keys on GenBank (``GCA_``) accessions
  while the cohort mostly carries RefSeq (``GCF_``) ones. Paired assemblies
  share the numeric part, so ``GCF_003830175.1`` is looked up as
  ``GCA_003830175.1``. Accessions that do not resolve are reported, not
  silently dropped.

Usage::

    python examples/klebsiella/fetch_phenotypes.py \
        --metadata examples/klebsiella/data/metadata.csv \
        --out examples/klebsiella/data/ast_phenotypes.csv
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

API = "https://www.bv-brc.org/api"
ASSEMBLY_RE = re.compile(r"^GC[FA]_\d+\.\d+$")


def _get(path: str, query: str, retries: int = 3, timeout: int = 90):
    url = f"{API}/{path}/?{query}&http_accept=application/json"
    last = None
    for attempt in range(retries):
        try:
            # BV-BRC's front end returns HTTP 403 for urllib's default
            # User-Agent. Identifying the client is enough; nothing is spoofed.
            req = urllib.request.Request(url, headers={
                "Accept": "application/json",
                "User-Agent": "amr-clonalshare/1.0.0 (research use)"})
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                return json.loads(fh.read().decode())
        except Exception as exc:                      # network flakiness only
            last = exc
            time.sleep(2 * (attempt + 1))
    print(f"  ! giving up on {path}: {last}", file=sys.stderr)
    return []


def _batched(seq, size):
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


def resolve_genome_ids(accessions, batch: int = 80) -> dict:
    """assembly accession (as given) -> BV-BRC genome_id."""
    out = {}
    gca = {a.replace("GCF_", "GCA_"): a for a in accessions}
    keys = sorted(gca)
    for n, chunk in enumerate(_batched(keys, batch), 1):
        q = "in(assembly_accession,({}))&select(genome_id,assembly_accession)&limit(1000)".format(
            ",".join(chunk))
        for rec in _get("genome", q):
            acc = rec.get("assembly_accession")
            if acc in gca:
                out[gca[acc]] = rec["genome_id"]
        print(f"  genome batch {n}: {len(out)} resolved so far", flush=True)
    return out


def fetch_ast(genome_ids, batch: int = 60) -> pd.DataFrame:
    """Laboratory-measured AST results for the given BV-BRC genome ids."""
    rows = []
    ids = sorted(set(genome_ids))
    for n, chunk in enumerate(_batched(ids, batch), 1):
        q = ("and(in(genome_id,({})),eq(evidence,%22Laboratory%20Method%22))"
             "&select(genome_id,antibiotic,resistant_phenotype,measurement,"
             "measurement_unit,laboratory_typing_method,evidence)&limit(20000)"
             ).format(",".join(chunk))
        recs = _get("genome_amr", q)
        rows.extend(recs)
        print(f"  amr batch {n}/{(len(ids) + batch - 1) // batch}: "
              f"{len(rows)} records", flush=True)
    df = pd.DataFrame(rows)
    if not df.empty:
        # The filter is the whole point of the script; verify it held.
        assert set(df["evidence"].unique()) == {"Laboratory Method"}, \
            f"non-laboratory evidence leaked in: {df['evidence'].unique()}"
    return df


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = Path(__file__).resolve().parent
    ap.add_argument("--metadata", default=str(here / "data" / "metadata.csv"))
    ap.add_argument("--out", default=str(here / "data" / "ast_phenotypes.csv"))
    ap.add_argument("--id-column", default="Strain_ID")
    args = ap.parse_args()

    meta = pd.read_csv(args.metadata)
    ids = meta[args.id_column].astype(str).tolist()
    acc = [i for i in ids if ASSEMBLY_RE.match(i)]
    print(f"{len(ids)} isolates, {len(acc)} with assembly accessions", flush=True)

    gmap = resolve_genome_ids(acc)
    print(f"resolved {len(gmap)}/{len(acc)} accessions to BV-BRC genomes",
          flush=True)

    ast = fetch_ast(list(gmap.values()))
    if ast.empty:
        print("no laboratory-method AST records found", file=sys.stderr)
        return 1

    rev = {v: k for k, v in gmap.items()}
    ast["Strain_ID"] = ast["genome_id"].map(rev)
    ast = ast.dropna(subset=["Strain_ID"])
    ast["antibiotic"] = ast["antibiotic"].astype(str).str.lower().str.strip()

    keep = ["Strain_ID", "genome_id", "antibiotic", "resistant_phenotype",
            "measurement", "measurement_unit", "laboratory_typing_method"]
    ast = ast[[c for c in keep if c in ast.columns]]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    ast.to_csv(args.out, index=False)

    n_iso = ast["Strain_ID"].nunique()
    print(f"\nwrote {args.out}")
    print(f"  {len(ast)} laboratory AST results")
    print(f"  {n_iso} isolates ({100 * n_iso / len(ids):.1f} % of the cohort)")
    print(f"  {ast['antibiotic'].nunique()} antibiotics")
    print(ast["antibiotic"].value_counts().head(15).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
