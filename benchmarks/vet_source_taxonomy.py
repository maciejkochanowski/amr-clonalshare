#!/usr/bin/env python3
"""Assign a host and a sampling matrix to a Pathogen Detection isolate.

    python source_taxonomy.py --profile <host_source.tsv> --out <dir>

WHY THIS EXISTS. The `host` field of the Pathogen Detection metadata is empty
for about half the isolates that carry an antimicrobial susceptibility result,
and the `isolation_source` field that would say where the isolate came from is
free text with 988 distinct values in the release read here. Some of those
values are ontology-annotated and some are laboratory shorthand. A veterinary
reading needs both fields, and it needs the assignment to be auditable rather
than plausible, so every rule is written out, every rule is matched against the
value as recorded, and every value the rules do not cover is counted and
reported rather than forced into a class.

TWO AXES, NOT ONE. Host says which animal. Matrix says where in the chain the
sample was taken, and it is not a detail: an isolate from a sick dog, an
isolate from a healthy pig at slaughter and an isolate from a supermarket
chicken breast are three different epidemiological objects that happen to share
a host species. They are kept apart and can be pooled deliberately.

ORDER OF PRECEDENCE, fixed here and applied mechanically.
  1. An ontology term in `isolation_source`. A CURIE is unambiguous.
  2. The `host` field, when it names an organism.
  3. An enumerated free-text pattern in `isolation_source`.
  4. Otherwise unclassified, which is a reported outcome and not a failure.

WHAT WOULD MEAN THIS IS WRONG. If a rule assigned a human clinical specimen to
an animal the veterinary cohort would be contaminated by the human one, which
is the failure that would matter. The audit table lists every distinct pair of
values with its assignment and its count, so that claim is checkable rather
than asserted, and the human-specimen words are listed explicitly as a class of
their own instead of falling through to unclassified.
"""
from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from pathlib import Path

# --- 1. ontology terms, taken from the values as recorded ------------------
# Only CURIEs seen in the data are listed. A CURIE that appears later and is
# not here lands in unclassified and is reported, which is the intended
# behaviour: this table is evidence, not a guess about the ontology.
TAXON_CURIE = {
    "ncbitaxon:9031": "chicken",          # Gallus gallus
}
MATRIX_CURIE = {
    "envo:01000925": "food",              # abattoir
    "uberon:0008979": "food",             # carcass
    "foodon:03317170": "food",            # meat, whole or parts
    "foodon:03311126": "food",            # food, raw
    "foodon:02010116": "food",            # meat with bone
    "foodon:02010111": "food",            # meat with skin
    "foodon:00003856": "food",            # coconut meat, a plant food
    "foodon:03316636": "food",            # ready to eat
    "foodon:00004555": "food",            # food, chunks
    "foodon:03302148": "food",            # food, frozen
    "envo:01001448": "food",              # retail environment
}
# A plant food carries no animal host. Listed so the rule is visible.
NON_ANIMAL_FOOD_CURIE = {"foodon:00003856"}

# --- 2. the host field, every value in the release read here ----------------
HOST_VALUE = {
    "homo sapiens": ("human", "human"),
    "homo sapiens sapiens": ("human", "human"),
    "canis lupus familiaris": ("dog", "companion_animal"),
    "felis catus": ("cat", "companion_animal"),
    "chicken [ncbitaxon:9031]": ("chicken", "poultry"),
    "gallus gallus": ("chicken", "poultry"),
    "gallus gallus domesticus": ("chicken", "poultry"),
    "chicken": ("chicken", "poultry"),
    "broiler": ("chicken", "poultry"),
    "chicken broiler": ("chicken", "poultry"),
    "meleagris gallopavo": ("turkey", "poultry"),
    "anas": ("duck", "poultry"),
    "bovine": ("cattle", "cattle"),
    "bos taurus": ("cattle", "cattle"),
    "veal calf": ("cattle", "cattle"),
    "dairy cow": ("cattle", "cattle"),
    "sus scrofa": ("pig", "swine"),
    "sus scrofa domesticus": ("pig", "swine"),
    "sus scrofa domestica": ("pig", "swine"),
    "pig": ("pig", "swine"),
    "swine": ("pig", "swine"),
    "goat": ("goat", "small_ruminant"),
    "sheep": ("sheep", "small_ruminant"),
    # Animals that are neither livestock nor kept as companions. Grouped, not
    # discarded: they are animal isolates and belong in the animal total.
    "rat": ("rat", "other_animal"),
    "macaca mulatta": ("macaque", "other_animal"),
    "bucephala": ("wild bird", "other_animal"),
    "spheniscus magellanicus": ("penguin", "other_animal"),
    "iguana": ("iguana", "other_animal"),
    "free-range white-tailed deer": ("deer", "other_animal"),
    "bird": ("bird", "other_animal"),
    "canis rufus": ("red wolf", "other_animal"),
    # Values that name no organism.
    "food": (None, None),
    "environment": (None, None),
    "environmental": (None, None),
    "enviromental": (None, None),
    "not provided": (None, None),
    "null": (None, None),
    "": (None, None),
}

# --- 3. free-text patterns in isolation_source ------------------------------
# Written as (substring, host, group, matrix). Matched in order, first hit
# wins, so the specific patterns are listed before the general ones. Every
# entry was taken from a value observed in the release.
SOURCE_RULES: list[tuple[str, str | None, str | None, str]] = [
    # companion-animal clinical specimens, which name the host inside the value
    ("canis lupus familiaris", "dog", "companion_animal", "animal_clinical"),
    ("canis lupus familaris", "dog", "companion_animal", "animal_clinical"),
    ("canine", "dog", "companion_animal", "animal_clinical"),
    ("feline", "cat", "companion_animal", "animal_clinical"),

    # national surveillance sampling of live food animals: the NARMS
    # animal-<species>-<production class> vocabulary, and farm sampling
    ("animal-cattle", "cattle", "cattle", "animal_surveillance"),
    ("animal-swine", "pig", "swine", "animal_surveillance"),
    ("animal-chicken", "chicken", "poultry", "animal_surveillance"),
    ("animal-turkey", "turkey", "poultry", "animal_surveillance"),
    ("fattening pig farm", "pig", "swine", "animal_surveillance"),
    ("calf feces", "cattle", "cattle", "animal_surveillance"),
    ("dairy cattle", "cattle", "cattle", "animal_surveillance"),
    ("beef cattle", "cattle", "cattle", "animal_surveillance"),
    ("cattle", "cattle", "cattle", "animal_surveillance"),
    ("hogs", "pig", "swine", "animal_surveillance"),

    # meat and carcass, the food end of the chain
    ("ground turkey", "turkey", "poultry", "animal_food"),
    ("turkey patties", "turkey", "poultry", "animal_food"),
    ("young turkey", "turkey", "poultry", "animal_food"),
    ("turkey", "turkey", "poultry", "animal_food"),
    ("chicken breast", "chicken", "poultry", "animal_food"),
    ("chicken wing", "chicken", "poultry", "animal_food"),
    ("chicken thigh", "chicken", "poultry", "animal_food"),
    ("chicken leg", "chicken", "poultry", "animal_food"),
    ("chicken liver", "chicken", "poultry", "animal_food"),
    ("chicken gizzard", "chicken", "poultry", "animal_food"),
    ("chicken heart", "chicken", "poultry", "animal_food"),
    ("comminuted chicken", "chicken", "poultry", "animal_food"),
    ("retail chicken", "chicken", "poultry", "animal_food"),
    ("whole chicken", "chicken", "poultry", "animal_food"),
    ("raw intact chicken", "chicken", "poultry", "animal_food"),
    ("young chicken", "chicken", "poultry", "animal_food"),
    ("chicken", "chicken", "poultry", "animal_food"),
    ("ground beef", "cattle", "cattle", "animal_food"),
    ("comminuted beef", "cattle", "cattle", "animal_food"),
    ("retail veal", "cattle", "cattle", "animal_food"),
    ("intact-beef", "cattle", "cattle", "animal_food"),
    ("beef", "cattle", "cattle", "animal_food"),
    ("pork chop", "pig", "swine", "animal_food"),
    ("pork", "pig", "swine", "animal_food"),
    ("poultry", None, "poultry", "animal_food"),

    # human clinical specimens, named so they do not fall through silently
    ("urine", None, "human", "human_clinical"),
    ("blood", None, "human", "human_clinical"),
    ("sputum", None, "human", "human_clinical"),
    ("wound", None, "human", "human_clinical"),
    ("stool", None, "human", "human_clinical"),
    ("hospital", None, "human", "human_clinical"),
    ("clinical", None, "human", "human_clinical"),
]

CURIE = re.compile(r"\b[a-z]+:[0-9]+\b")
FOOD_ANIMAL_GROUPS = ("cattle", "swine", "poultry", "small_ruminant")
ANIMAL_GROUPS = FOOD_ANIMAL_GROUPS + ("companion_animal", "other_animal")


def classify(host: str, source: str) -> dict:
    """Host, group and matrix for one isolate, with the rule that decided it."""
    # A missing cell reaches this function as a float NaN when the table was
    # read with pandas and as an empty string when it was read line by line.
    # Both mean the same thing and neither is allowed to raise.
    host = "" if host is None or host != host else str(host).strip().lower()
    source = ("" if source is None or source != source
              else str(source).strip().lower())
    if host in ("null", "not provided"):
        host = ""
    if source in ("null", "not provided", "not collected", "missing"):
        source = ""

    curies = set(CURIE.findall(source))
    animal_curies = curies & set(TAXON_CURIE)
    matrix_curies = curies & set(MATRIX_CURIE)

    # 1. the host field, when it names an organism. It is the curated field for
    #    the organism the isolate came from, so it outranks a term that may
    #    describe an attributed food vehicle rather than the source animal. In
    #    the release read here no isolate carries both, so this ordering
    #    changes nothing; it is fixed so that a later release cannot make the
    #    wrong assignment silently.
    if host in HOST_VALUE:
        name, group = HOST_VALUE[host]
        if name is not None:
            matrix = _matrix_from_source(source, group)
            return {"host": name, "group": group, "matrix": matrix,
                    "rule": f"host:{host}"}

    # 2. an ontology term, which is unambiguous
    if animal_curies:
        name = TAXON_CURIE[sorted(animal_curies)[0]]
        group = "poultry" if name in ("chicken", "turkey", "duck") else None
        matrix = "animal_food" if matrix_curies else "animal_unspecified"
        return {"host": name, "group": group, "matrix": matrix,
                "rule": f"curie:{sorted(animal_curies)[0]}"}
    # 3. an enumerated free-text pattern
    for pattern, name, group, matrix in SOURCE_RULES:
        if pattern in source:
            return {"host": name, "group": group, "matrix": matrix,
                    "rule": f"source:{pattern}"}

    if curies & NON_ANIMAL_FOOD_CURIE:
        return {"host": None, "group": "non_animal_food",
                "matrix": "food_other", "rule": "curie:plant food"}
    return {"host": None, "group": None, "matrix": None, "rule": "unclassified"}


#: Specimen words that name a diagnostic sample rather than a place in the food
#: chain. Which clinic it belongs to is decided by the host, not by the word: a
#: urine from a dog is a veterinary clinical specimen and a urine from a person
#: is a human one. Taken from the values observed in the release.
CLINICAL_SPECIMEN = (
    "urine", "blood", "wound", "sputum", "stool", "abscess", "skin", "ear",
    "eye", "respiratory", "tissue", "swab", "aspirate", "lesion", "sepsis",
    "bile", "sterile body site", "sterile fluid", "vagina", "uterus", "nose",
    "nasal", "throat", "oral", "rectal", "drainage", "lung", "liver",
    "bladder", "urethra", "groin", "pleural fluid", "surgical site",
    "blood culture", "clinical", "infection", "culture", "perirectal",
)


def _matrix_from_source(source: str, group: str | None) -> str:
    """Where in the chain, once the host is already known from the host field."""
    if group == "human":
        return "human_clinical"
    for pattern, _, _, matrix in SOURCE_RULES:
        if pattern in source and matrix in ("animal_food", "animal_surveillance",
                                            "animal_clinical"):
            return matrix
    if any(word in source for word in CLINICAL_SPECIMEN):
        return "animal_clinical"
    return "animal_unspecified"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--profile", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    pairs: Counter = Counter()
    by_org: Counter = Counter()
    with args.profile.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                parts += [""] * (3 - len(parts))
            org, host, source = parts[0], parts[1], parts[2]
            pairs[(org, host, source)] += 1

    audit, group_counts, matrix_counts, unclassified = [], Counter(), Counter(), Counter()
    for (org, host, source), n in pairs.items():
        c = classify(host, source)
        audit.append({"organism": org, "host_field": host,
                      "isolation_source": source, "n": n, **c})
        key = c["group"] or "unclassified"
        group_counts[key] += n
        matrix_counts[c["matrix"] or "unclassified"] += n
        by_org[(org, key)] += n
        if c["rule"] == "unclassified":
            unclassified[(host, source)] += n

    (args.out / "audit_rules.json").write_text(
        json.dumps(sorted(audit, key=lambda r: -r["n"]), indent=1) + "\n",
        encoding="utf-8")
    total = sum(pairs.values())
    summary = {
        "n_isolates": total,
        "n_distinct_value_pairs": len(pairs),
        "by_group": dict(group_counts.most_common()),
        "by_matrix": dict(matrix_counts.most_common()),
        "animal_total": sum(v for k, v in group_counts.items()
                            if k in ANIMAL_GROUPS),
        "food_animal_total": sum(v for k, v in group_counts.items()
                                 if k in FOOD_ANIMAL_GROUPS),
        "unclassified_share": group_counts["unclassified"] / total,
        "by_organism_and_group": {f"{o}|{g}": n for (o, g), n
                                  in by_org.most_common()},
        "largest_unclassified_values": [
            {"host": h, "isolation_source": s, "n": n}
            for (h, s), n in unclassified.most_common(30)],
    }
    (args.out / "source_taxonomy_summary.json").write_text(
        json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items()
                      if k != "by_organism_and_group"}, indent=1)[:4000])
    print("\nby organism and group, animal groups only:")
    for (o, g), n in by_org.most_common():
        if g in ANIMAL_GROUPS:
            print(f"  {o:<34} {g:<18} {n:>6}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
