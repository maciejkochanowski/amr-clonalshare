"""The S. suis clonal share per agent, per block, and against the determinants.

    python benchmarks/ssuis_mechanism.py --out benchmarks/results_ssuis_mechanism

Writes the numbers the article quotes for the shipped *Streptococcus suis*
cohort, as JSON, with a receipt: the out-of-sample share of each agent's
non-wild-type variance that the hierBAPS population cluster explains, the same
share for the two mechanism blocks and for the whole panel, and the check the
article rests on, which uses a layer the estimate never sees. The source
study calls 43 resistance determinants from the same genomes; carriage of the
determinants it reports for each agent's mechanism is scored as "any of
them", its clonal share is computed identically, and the two shares are
compared agent by agent. The mechanism classes, chromosomal, mobile or mixed,
are read from the genetic location the source study reports per determinant
and counted per carrier in this cohort; they are used only to score the table.
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
from amr_clonalshare.attribution import clonal_share, layer_clonal_share

SEED = 20260901


def per_agent(X, cols, lineage, layer_of, n_boot, n_perm):
    rows = []
    for j, c in enumerate(cols):
        r = clonal_share(X[:, j], lineage, seed=SEED + j,
                         n_boot=n_boot, n_perm=n_perm)
        d = r.as_dict()
        d.update(agent=c, layer=layer_of(c))
        rows.append(d)
    return rows


def load_ssuis():
    D = ROOT / "examples/ssuis/data"
    ribo = pd.read_csv(D / "ribo.csv")
    cell = pd.read_csv(D / "cell.csv")
    meta = pd.read_csv(D / "metadata.csv", dtype=str)
    src = pd.read_csv(D / "source_table_s1.csv", dtype=str)
    df = (ribo.merge(cell, on="genome_id")
              .merge(meta[["genome_id", "baps_cluster", "mlst", "source_strain"]],
                     on="genome_id")
              .merge(src[["Strain"] + [c for c in src.columns if c.startswith("G__")]],
                     left_on="source_strain", right_on="Strain", how="left"))
    RIBO = [c for c in ribo.columns if c != "genome_id"]
    CELL = [c for c in cell.columns if c != "genome_id"]
    return df, RIBO, CELL


# Agent -> the determinants the source study calls for that agent's mechanism.
# Carriage is scored as "any of these", because a lineage may reach the same
# phenotype by a different member of the same family and the question here is
# the mechanism class, not the allele. The keys match the column names of the
# source table, whose suffix codes the reported genetic location: ``_M_`` a
# mobile element (plasmid, ICE or transposon), ``_C_`` the chromosome, ``_V_``
# a chromosomal haplotype or variant.
DETERMINANTS = {
    "erythromycin": ["ermB", "ermG", "ermT"],
    "tylosin": ["ermB", "ermG", "ermT"],
    "tilmicosin": ["ermB", "ermG", "ermT"],
    "lincomycin": ["ermB", "linB", "lnuB", "lnuC"],
    "tetracycline": ["tetM__", "tetO_", "tetW__", "tet44", "tetL"],
    "doxycycline": ["tetM__", "tetO_", "tetW__", "tet44", "tetL"],
    "spectinomycin": ["ant_6_Ia", "ant_9_Ia", "aph_3_IIIa", "ant_4_Ib"],
    "tiamulin": ["lsaE", "vgaF"],
    "trimethoprim": ["dfrF", "dfrK", "DHFR_102", "DHFRPromoter"],
    "penicillin": ["PBP2B.hap", "PBP2X.hap", "PBP2X_551", "MraY.hap"],
    "amoxicillin": ["PBP2B.hap", "PBP2X.hap", "PBP2X_551", "MraY.hap"],
    "ceftiofur": ["PBP2B.hap", "PBP2X.hap", "PBP2X_551", "MraY.hap"],
    "cefquinome": ["PBP2B.hap", "PBP2X.hap", "PBP2X_551", "MraY.hap"],
}
#: The determinant sets, so that agents scored on one set are one unit.
DETERMINANT_SET = {
    "erythromycin": "erm", "tylosin": "erm", "tilmicosin": "erm",
    "lincomycin": "lnu/lin", "tetracycline": "tet", "doxycycline": "tet",
    "spectinomycin": "ant/aph", "tiamulin": "lsaE/vgaF",
    "trimethoprim": "dfr/DHFR", "penicillin": "PBP", "amoxicillin": "PBP",
    "ceftiofur": "PBP", "cefquinome": "PBP",
}


def _location(column: str) -> str:
    """The reported genetic location coded in a source-table column name."""
    code = column.split("_")[-2] if column.count("_") >= 2 else ""
    return {"M": "mobile", "C": "chromosomal", "V": "chromosomal"}.get(code, "unknown")


def mechanism_class(df, keys):
    """How an agent's determinants are carried in this cohort.

    ``chromosomal`` when every carrier's determinant is on the chromosome,
    ``mobile`` when every one is on a mobile element, and ``mixed`` when the
    cohort holds both, with the carrier counts of each so the reader can see
    which dominates. The class is the source study's location column read
    per determinant and counted per isolate; nothing here is fitted.
    """
    counts = {"mobile": 0, "chromosomal": 0}
    for c in df.columns:
        if c.startswith("G__") and any(k in c for k in keys):
            loc = _location(c)
            if loc in counts:
                counts[loc] += int((df[c].astype(str) == "yes").sum())
    if counts["mobile"] and counts["chromosomal"]:
        cls = "mixed"
    elif counts["mobile"]:
        cls = "mobile"
    else:
        cls = "chromosomal"
    return cls, counts


def determinant_carriage(df, keys):
    cols = [c for c in df.columns
            if c.startswith("G__") and any(k in c for k in keys)]
    if not cols:
        return None, []
    v = np.zeros(len(df), dtype=float)
    for c in cols:
        v = np.maximum(v, (df[c].astype(str) == "yes").astype(float).to_numpy())
    return v, cols


def ordering(rows):
    """How the shares sit against the way the determinants are carried.

    Two readings, both descriptive. Per agent: of the pairs of one agent
    whose determinants are chromosomal (or chromosomal for most carriers,
    the ``mixed`` class) against one whose determinants are mobile, how many
    place the former higher. Per determinant set: the agents scored on one
    set are one unit, so the sets are compared by the mean share of their
    agents; with three non-mobile sets against four mobile ones there are
    twelve such pairs. No test is attached, because the units are few and
    the agents within a set are not independent observations.
    """
    by = {r["agent"]: r for r in rows}
    cls = {a: r["mechanism"] for a, r in by.items()}
    share = {a: r["kappa_adj"] for a, r in by.items()}
    chrom = [a for a in by if cls[a] == "chromosomal"]
    mixed = [a for a in by if cls[a] == "mixed"]
    mob = [a for a in by if cls[a] == "mobile"]
    non_mobile = chrom + mixed
    sets = {}
    for a in by:
        sets.setdefault(DETERMINANT_SET.get(a, a), []).append(a)
    set_mean = {k: float(np.mean([share[a] for a in v])) for k, v in sets.items()}
    set_class = {k: cls[v[0]] for k, v in sets.items()}
    nm_sets = [k for k in sets if set_class[k] != "mobile"]
    m_sets = [k for k in sets if set_class[k] == "mobile"]
    return dict(
        n_chromosomal=len(chrom), n_mixed=len(mixed), n_mobile=len(mob),
        mean_chromosomal=float(np.mean([share[a] for a in chrom])) if chrom else None,
        mean_non_mobile=float(np.mean([share[a] for a in non_mobile])) if non_mobile else None,
        mean_mobile=float(np.mean([share[a] for a in mob])) if mob else None,
        pairs_chromosomal_vs_mobile=len(chrom) * len(mob),
        pairs_chromosomal_higher=int(sum(share[c] > share[m] for c in chrom for m in mob)),
        pairs_non_mobile_vs_mobile=len(non_mobile) * len(mob),
        pairs_non_mobile_higher=int(sum(share[c] > share[m] for c in non_mobile for m in mob)),
        determinant_sets={k: dict(agents=v, mechanism=set_class[k], mean_share=set_mean[k])
                          for k, v in sets.items()},
        set_pairs_total=len(nm_sets) * len(m_sets),
        set_pairs_non_mobile_higher=int(sum(set_mean[c] > set_mean[m]
                                            for c in nm_sets for m in m_sets)),
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="benchmarks/results_ssuis_mechanism")
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    n_boot, n_perm = (100, 100) if a.quick else (400, 200)
    t0 = time.time()
    res: dict = {}

    df, RIBO, CELL = load_ssuis()
    traits = RIBO + CELL
    Xs = df[traits].to_numpy(float)
    baps = df["baps_cluster"].to_numpy()

    res["ssuis_per_agent"] = per_agent(
        Xs, traits, baps, lambda c: "ribosomal" if c in RIBO else "cell_folate",
        n_boot, n_perm)
    res["ssuis_layers"] = {}
    for nm, cols in (("ribosomal", RIBO), ("cell_folate", CELL), ("all", traits)):
        r = layer_clonal_share(df[cols].to_numpy(float), baps, seed=SEED,
                               n_boot=n_boot, n_perm=n_perm)
        res["ssuis_layers"][nm] = r.as_dict()
    for r in res["ssuis_per_agent"]:
        cls, counts = mechanism_class(df, DETERMINANTS[r["agent"]])
        r["mechanism"] = cls
        r["determinant_carriers"] = counts
        r["determinant_set"] = DETERMINANT_SET[r["agent"]]
    res["ssuis_ordering"] = ordering(res["ssuis_per_agent"])

    # external layer: the determinants, which built nothing above
    ext = []
    for agent, keys in DETERMINANTS.items():
        v, cols = determinant_carriage(df, keys)
        if v is None or not (0.02 <= v.mean() <= 0.98):
            continue
        rg = clonal_share(v, baps, seed=SEED, n_boot=n_boot, n_perm=n_perm)
        ph = [r for r in res["ssuis_per_agent"] if r["agent"] == agent][0]
        ext.append(dict(agent=agent, mechanism=mechanism_class(df, keys)[0],
                        determinants=cols,
                        carriage=float(v.mean()),
                        kappa_genotype=rg.kappa_adj, kappa_phenotype=ph["kappa_adj"],
                        ci_genotype=[rg.ci_low, rg.ci_high]))
    res["ssuis_genotype_check"] = ext

    g = np.array([e["kappa_genotype"] for e in ext])
    p = np.array([e["kappa_phenotype"] for e in ext])
    res["ssuis_genotype_agreement"] = dict(
        n_pairs=len(ext),
        pearson=float(np.corrcoef(g, p)[0, 1]) if len(ext) > 2 else float("nan"),
        spearman=float(pd.Series(g).corr(pd.Series(p), method="spearman")),
        mean_abs_diff=float(np.abs(g - p).mean()),
        mobile_mean_phenotype=float(p[[e["mechanism"] == "mobile" for e in ext]].mean()),
        non_mobile_mean_phenotype=float(
            p[[e["mechanism"] != "mobile" for e in ext]].mean()),
        mobile_mean_genotype=float(g[[e["mechanism"] == "mobile" for e in ext]].mean()),
        non_mobile_mean_genotype=float(
            g[[e["mechanism"] != "mobile" for e in ext]].mean()),
    )

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "cells.json").write_text(json.dumps(res, indent=1, default=float))
    (out / "RECEIPT.json").write_text(json.dumps({
        "seed": SEED, "quick": a.quick, "n_boot": n_boot, "n_perm": n_perm,
        "n_isolates": int(len(df)), "n_agents": len(traits),
        "runtime_s": round(time.time() - t0, 1),
        "numpy": np.__version__, "pandas": pd.__version__,
        "python": sys.version.split()[0],
    }, indent=1))
    o = res["ssuis_ordering"]
    print(f"chromosomal {o['n_chromosomal']}, mixed {o['n_mixed']}, mobile {o['n_mobile']} agents; "
          f"chromosomal mean {o['mean_chromosomal']:.3f}, non-mobile mean {o['mean_non_mobile']:.3f}, "
          f"mobile mean {o['mean_mobile']:.3f}; pairs chromosomal>mobile "
          f"{o['pairs_chromosomal_higher']}/{o['pairs_chromosomal_vs_mobile']}, non-mobile>mobile "
          f"{o['pairs_non_mobile_higher']}/{o['pairs_non_mobile_vs_mobile']}; determinant sets "
          f"{o['set_pairs_non_mobile_higher']}/{o['set_pairs_total']}")
    print(f"wrote {out}   ({time.time() - t0:.0f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
