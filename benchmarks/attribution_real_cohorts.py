"""Clonal attribution on the two shipped cohorts, with an external check.

Writes the numbers the manuscript quotes, as JSON, with a receipt. Nothing here
is copied from a historical run: every figure in the article's attribution
section comes out of this script.

THE PARTITION IS THE PIPELINE'S, NOT THIS SCRIPT'S. An earlier version built a
Ward partition of the same calls and reported it under the same heading as the
article's own figures. That was wrong twice over: the article is about the
partition the software returns, and two partitions of one cohort reported side
by side as though they were one made the supplement contradict the main text,
in one case reversing which of the lineage and the partition explained more.
The partition is now read from the campaign the article cites and joined on the
isolate identifier, with the join asserted complete.

  python benchmarks/attribution_real_cohorts.py --out benchmarks/results_attribution_real
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from amr_clonalshare.attribution import (  # noqa: E402
    SUPPORT_THRESHOLD, attribute_partition, clonal_share, concordance_z,
    layer_clonal_share,
)

SEED = 20260901


def ward(X, k):
    return fcluster(linkage(pdist(X, "euclidean"), "ward"), k, "maxclust") - 1


def per_agent(X, cols, lineage, layer_of, n_boot, n_perm):
    rows = []
    for j, c in enumerate(cols):
        r = clonal_share(X[:, j], lineage, seed=SEED + j,
                         n_boot=n_boot, n_perm=n_perm)
        d = r.as_dict()
        d.update(agent=c, layer=layer_of(c))
        rows.append(d)
    return rows



CAMPAIGN = (ROOT / "paper" / "softwarex" / "evidence" / "campaign_2026-09-01"
            / "results")
ID_COLUMNS = ("isolate_id", "isolate", "genome_id", "Strain_ID",
              "id", "sample_id")


def campaign_partition(name: str, frame: pd.DataFrame) -> np.ndarray:
    """The cluster label the pipeline assigned, read from the campaign.

    Joined on the isolate identifier as text on both sides. An identifier that
    looks numeric is parsed as a float on one side and indexed as a string on
    the other, and a silent partial join here would put the wrong labels
    against the right traits, so an incomplete join raises instead.
    """
    table = pd.read_csv(CAMPAIGN / name / "assignment.tsv", sep="\t", dtype=str)
    column = next((c for c in ID_COLUMNS if c in frame.columns), None)
    if column is None:
        raise RuntimeError(f"{name}: no isolate identifier column in the cohort "
                           f"frame; looked for {ID_COLUMNS}")
    lookup = dict(zip(table["isolate_id"].astype(str), table["cluster"]))
    labels = [lookup.get(str(v)) for v in frame[column]]
    missing = sum(1 for v in labels if v is None)
    if missing:
        raise RuntimeError(f"{name}: {missing} of {len(labels)} isolates carry "
                           f"no campaign cluster; the join is not complete")
    codes = {v: i for i, v in enumerate(sorted(set(labels)))}
    return np.asarray([codes[v] for v in labels], dtype=int)


def partition_row(X, lineage, k, n_boot, labels=None):
    # ``labels`` is the partition the pipeline returned. The Ward fallback is
    # kept for the sensitivity row that deliberately uses a different
    # clustering route, and that row says so in its own key.
    lab = ward(X, k) if labels is None else np.asarray(labels, dtype=int)
    k = int(lab.max()) + 1
    a = attribute_partition(X, lab, lineage, seed=SEED, n_boot=n_boot)
    obs, mu, sd, z, p = concordance_z(lab, lineage, n_perm=500, seed=SEED)
    d = a.as_dict()
    d.update(k=k, partition_source="pipeline" if labels is not None else "ward",
             gate_concordance=obs, gate_null_mean=mu, gate_null_sd=sd,
             gate_z=z, gate_p=p, cluster_sizes=np.bincount(lab).tolist())
    return d, lab


# ------------------------------------------------------------------- S. suis
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
# the mechanism class, not the allele.
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
#: Which of those mechanisms the source study reports as mobile-element borne
#: rather than chromosomal. Used only to label the table, never to fit anything.
MOBILE = {"erythromycin", "tylosin", "tilmicosin", "lincomycin", "tetracycline",
          "doxycycline", "spectinomycin", "tiamulin"}


def determinant_carriage(df, keys):
    cols = [c for c in df.columns
            if c.startswith("G__") and any(k in c for k in keys)]
    if not cols:
        return None, []
    v = np.zeros(len(df), dtype=float)
    for c in cols:
        v = np.maximum(v, (df[c].astype(str) == "yes").astype(float).to_numpy())
    return v, cols


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="benchmarks/results_attribution_real")
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
    res["ssuis_partition"], lab_s = partition_row(
        Xs, baps, 5, n_boot, labels=campaign_partition("ssuis", df))

    # sensitivity: the same partition judged against the other lineage variable
    typed = df["mlst"].notna().to_numpy()
    res["ssuis_partition_mlst_typed_only"], _ = partition_row(
        Xs[typed], df["mlst"].to_numpy()[typed], 5, n_boot,
        labels=campaign_partition("ssuis", df)[typed])

    # external layer: the determinants, which built nothing above
    ext = []
    for agent, keys in DETERMINANTS.items():
        v, cols = determinant_carriage(df, keys)
        if v is None or not (0.02 <= v.mean() <= 0.98):
            continue
        rg = clonal_share(v, baps, seed=SEED, n_boot=n_boot, n_perm=n_perm)
        ph = [r for r in res["ssuis_per_agent"] if r["agent"] == agent][0]
        ext.append(dict(agent=agent, mechanism="mobile" if agent in MOBILE
                        else "chromosomal", determinants=cols,
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
        chromosomal_mean_phenotype=float(
            p[[e["mechanism"] == "chromosomal" for e in ext]].mean()),
        mobile_mean_genotype=float(g[[e["mechanism"] == "mobile" for e in ext]].mean()),
        chromosomal_mean_genotype=float(
            g[[e["mechanism"] == "chromosomal" for e in ext]].mean()),
    )

    # how the gate statistic scales while the attribution does not
    scaling = []
    rng = np.random.default_rng(SEED)
    for frac in (0.1, 0.2, 0.4, 0.7, 1.0):
        m = max(int(frac * len(df)), 60)
        sub = rng.choice(len(df), size=m, replace=False)
        d, lb = partition_row(Xs[sub], baps[sub], 5, max(n_boot // 4, 50),
                              labels=campaign_partition("ssuis", df)[sub])
        scaling.append(dict(n=m, lam=d["lam"], gate_z=d["gate_z"],
                            r2_lineage=d["r2_lineage"],
                            r2_partition=d["r2_partition"]))
    res["gate_scaling"] = scaling

    # --------------------------------------------------------- Klebsiella
    K = ROOT / "examples/klebsiella/data"
    amr = pd.read_csv(K / "amr.csv")
    kmeta = pd.read_csv(K / "metadata.csv", dtype=str)
    kk = amr.merge(kmeta[["Strain_ID", "ST"]], on="Strain_ID")
    classes = [c for c in amr.columns if c != "Strain_ID"]
    keep = [c for c in classes if 0.02 <= kk[c].mean() <= 0.98]
    Xk = kk[keep].to_numpy(float)
    st = kk["ST"].to_numpy()
    res["klebsiella_dropped_classes"] = sorted(set(classes) - set(keep))
    res["klebsiella_per_agent"] = per_agent(Xk, keep, st, lambda c: "amr_class",
                                            n_boot, n_perm)
    kleb_labels = campaign_partition("klebsiella", kk)
    res["klebsiella_partition"], _ = partition_row(Xk, st, 2, n_boot,
                                                   labels=kleb_labels)

    vc = pd.Series(st).value_counts()
    m = np.array([s in set(vc[vc >= 2].index) for s in st])
    res["klebsiella_nonsingleton_per_agent"] = per_agent(
        Xk[m], keep, st[m], lambda c: "amr_class", n_boot, n_perm)
    res["klebsiella_nonsingleton_partition"], _ = partition_row(
        Xk[m], st[m], 2, n_boot, labels=kleb_labels[m])
    res["klebsiella_nonsingleton_n"] = int(m.sum())

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "cells.json").write_text(json.dumps(res, indent=1, default=float))
    (out / "RECEIPT.json").write_text(json.dumps({
        "seed": SEED, "quick": a.quick, "n_boot": n_boot, "n_perm": n_perm,
        "support_threshold": SUPPORT_THRESHOLD,
        "runtime_s": round(time.time() - t0, 1),
        "numpy": np.__version__, "pandas": pd.__version__,
        "python": sys.version.split()[0],
    }, indent=1))
    print(json.dumps({k: (v if not isinstance(v, list) else f"[{len(v)} rows]")
                      for k, v in res.items()}, indent=1, default=float)[:400])
    print(f"wrote {out}  ({time.time() - t0:.0f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
