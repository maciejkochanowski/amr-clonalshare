#!/usr/bin/env python3
"""Attach the source study's population clusters and determinants to this cohort.

    python examples/ssuis/link_source_lineages.py            # write and validate
    python examples/ssuis/link_source_lineages.py --check    # validate only

**Why a second lineage definition is needed.** The MLST sequence type served by
BV-BRC is missing for 219 of 677 isolates, and it is not missing at random.
Among United Kingdom isolates collected from 2013 only 40 % carry a sequence
type, against 97 % of those collected before 2013, and the untyped isolates
carry about two more non-wild-type results out of thirteen than the typed ones.
A contrast between those periods computed on sequence-typed isolates is
therefore computed on subsets that are not comparable across its own arms: for
ceftiofur the typed subset falls by 11 points while the collection it is drawn
from rises by 9.

The source publication solves this. It reports a hierBAPS cluster for every one
of its 678 isolates, inferred from the core-genome alignment (Cheng, Connor,
Siren, Aanensen & Corander 2013, *Mol Biol Evol* 30:1224-1228). A hierBAPS
cluster is also the better lineage for this question than a seven-locus
sequence type, because it is what the population structure is rather than a
proxy for it.

**The join is on the strain identifier, not on a statistical match.** BV-BRC
carries the source study's own strain name inside `genome_name`
("Streptococcus suis SS1038"), and the longest source identifier contained in
that string identifies the isolate. All 677 isolates match, each to a distinct
strain, with no collisions.

**Two variables validate the join and neither takes part in it.** Collection
year is recorded independently in both tables and must agree on every isolate
where both hold it. The sixteen-agent MIC vector is recorded independently in
both and must agree cell for cell. Both are enforced here; the script writes
nothing if either fails.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
SOURCE = DATA / "source_table_s1.csv"
METADATA = DATA / "metadata.csv"
MIC_LONG = DATA / "mic_long.csv"
RECEIPT = DATA / "lineage_link_receipt.json"

#: Dilutions that submitters report at two precisions. Snapping them is a
#: presentation fix, not a change of value: no pair spans a doubling step.
DILUTION_ALIASES = {0.031: 0.03, 0.063: 0.06, 0.125: 0.12, 0.016: 0.015}

REQUIRED_MATCH_RATE = 1.0
REQUIRED_YEAR_AGREEMENT = 1.0
REQUIRED_MIC_AGREEMENT = 1.0


def canonical_mic(value) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "NA"
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return "NA"
    for reported, canonical in DILUTION_ALIASES.items():
        if abs(number - reported) < 1e-9:
            number = canonical
    return f"{number:.4g}"


def join_on_strain_name(metadata: pd.DataFrame,
                        source: pd.DataFrame) -> pd.Series:
    """Longest source strain identifier contained in the BV-BRC genome name.

    Longest wins so that ``SS103`` cannot claim a genome named ``SS1038``.
    """
    identifiers = sorted((s for s in source["Strain"].astype(str) if s),
                         key=len, reverse=True)

    def find(name: str) -> str | None:
        text = str(name)
        for identifier in identifiers:
            if identifier in text:
                return identifier
        return None

    return metadata["genome_name"].map(find)


def validate(metadata: pd.DataFrame, source: pd.DataFrame,
             strain: pd.Series, mic: pd.DataFrame) -> dict:
    matched = strain.notna()
    indexed = source.set_index("Strain")

    served_year = pd.to_numeric(metadata["collection_year"], errors="coerce")
    published_year = pd.to_numeric(
        indexed["Year"].reindex(strain[matched]).to_numpy(), errors="coerce")
    comparable = (np.isfinite(served_year[matched].to_numpy())
                  & np.isfinite(published_year))
    year_agreement = (
        float((served_year[matched].to_numpy()[comparable]
               == published_year[comparable]).mean())
        if comparable.any() else float("nan"))

    renamed = {c: ("florfenicol" if c.endswith("Florifenicol")
                   else c.replace("P__", "").split("_", 1)[1].lower())
               for c in source.columns if c.startswith("P__")}
    published_mic = indexed.rename(columns=renamed)
    served_mic = mic.pivot(index="genome_id", columns="antibiotic",
                           values="measurement")
    served_mic.index = "g" + served_mic.index.astype(str)
    agents = sorted(set(served_mic.columns) & set(renamed.values()))
    left = served_mic.reindex(metadata["genome_id"][matched])[agents] \
        .map(canonical_mic).to_numpy()
    right = published_mic.reindex(strain[matched])[agents] \
        .map(canonical_mic).to_numpy()
    mic_agreement = float((left == right).mean())

    clusters = indexed["NewBaps"].astype(str).reindex(strain[matched])
    clusters.index = metadata["genome_id"][matched].to_numpy()
    sequence_type = metadata.set_index("genome_id")["mlst"].reindex(
        clusters.index)
    both = sequence_type.notna().to_numpy()
    information = float(adjusted_mutual_info_score(
        sequence_type[both], clusters[both])) if both.sum() > 2 else float("nan")
    generator = np.random.default_rng(0)
    null = [adjusted_mutual_info_score(
        sequence_type[both], generator.permutation(clusters[both].to_numpy()))
        for _ in range(200)]
    return {
        "n_isolates": int(len(metadata)),
        "n_matched": int(matched.sum()),
        "match_rate": float(matched.mean()),
        "n_distinct_strains_claimed": int(strain.nunique()),
        "n_agents_compared": len(agents),
        "n_mic_cells_compared": int(left.size),
        "mic_agreement": mic_agreement,
        "n_year_comparable": int(comparable.sum()),
        "year_agreement": year_agreement,
        "n_with_sequence_type": int(both.sum()),
        "adjusted_mutual_information_with_mlst": information,
        "permutation_null_mean": float(np.mean(null)),
        "permutation_null_sd": float(np.std(null)),
        "z_against_permutation_null": float(
            (information - np.mean(null)) / np.std(null)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    source = pd.read_csv(SOURCE, dtype=str)
    metadata = pd.read_csv(METADATA, dtype=str)
    mic = pd.read_csv(MIC_LONG, dtype={"genome_id": str})

    strain = join_on_strain_name(metadata, source)
    checks = validate(metadata, source, strain, mic)
    for name, value in checks.items():
        print(f"  {name}: {value}")

    failures = []
    if checks["match_rate"] < REQUIRED_MATCH_RATE:
        failures.append(f"only {checks['n_matched']} of {checks['n_isolates']} "
                        f"isolates match a source strain")
    if checks["n_distinct_strains_claimed"] != checks["n_matched"]:
        failures.append("two genomes claim the same source strain")
    if checks["year_agreement"] < REQUIRED_YEAR_AGREEMENT:
        failures.append(f"collection year agrees on only "
                        f"{checks['year_agreement']:.4f} of matched isolates")
    if checks["mic_agreement"] < REQUIRED_MIC_AGREEMENT:
        failures.append(f"the MIC tables agree on only "
                        f"{checks['mic_agreement']:.4f} of compared cells")
    if failures:
        for failure in failures:
            print(f"FAILED: {failure}")
        return 1

    indexed = source.set_index("Strain")
    clusters = ("BAPS" + indexed["NewBaps"].astype(str)).reindex(
        strain).to_numpy()

    if args.check:
        current = metadata.get("baps_cluster")
        if current is None:
            print("state: baps_cluster is not in metadata.csv")
            return 1
        same = current.fillna("").to_numpy() == pd.Series(clusters).fillna("").to_numpy()
        print(f"state: {int(same.sum())}/{len(clusters)} labels already match")
        return 0 if bool(same.all()) else 1

    metadata["source_strain"] = strain.to_numpy()
    metadata["baps_cluster"] = clusters
    metadata.to_csv(METADATA, index=False)
    RECEIPT.write_text(json.dumps({
        "generated_utc": datetime.now(timezone.utc).isoformat()
        .replace("+00:00", "Z"),
        "source": "Hadjirin et al. 2021 BMC Biology 19:191, Additional file 1",
        "join": "longest source Strain identifier contained in BV-BRC "
                "genome_name",
        "validated_by": ["collection year", "sixteen-agent MIC vector"],
        "n_clusters": int(pd.Series(clusters).nunique()),
        "validation": checks,
        "thresholds": {
            "required_match_rate": REQUIRED_MATCH_RATE,
            "required_year_agreement": REQUIRED_YEAR_AGREEMENT,
            "required_mic_agreement": REQUIRED_MIC_AGREEMENT,
        },
    }, indent=1))
    print(f"written: {METADATA} and {RECEIPT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
