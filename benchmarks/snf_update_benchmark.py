"""Which cross-diffusion update is the reference one, and which recovers truth?

Two questions, answered separately.

**Provenance.** ``fusion.snf_fuse`` re-applies the row normalisation after every
diffusion step. An earlier release of this package documented that as a
deviation from "the reference implementation". That was wrong, and the error is
worth stating precisely because it changes the default: there are *two*
reference lineages, and they disagree with each other.

* Wang's MATLAB code, ``SNFtool`` <= 2.2 (CRAN, 2014) and ``snfpy`` (a port of
  the MATLAB) add a constant to the diagonal and let the off-diagonal row mass
  float: ``Wall[[j]] = nextW[[j]] + diag(n)`` / ``_B0_normalized(W, alpha)``.
  Normalisation happens once at the start and once at the end.
* ``SNFtool`` >= 2.2.1 (CRAN, 2017-11-24) replaced that line with
  ``Wall[[j]] <- normalize(nextW[[j]])`` and simultaneously replaced
  ``normalize(X) = X / rowSums(X)`` with
  ``X / (2 * (rowSums(X) - diag(X))); diag(X) <- 0.5`` -- i.e. per-iteration
  renormalisation with self-weight 1/2, which is what this package does and
  what equation (2) of the paper describes. That is still the current CRAN
  release (2.3.1, 2021-06-11).

So the update implemented here is the maintained reference update; the update
this package called "the reference" is the 2014 one that its own authors
replaced. This script verifies both claims by reading the sources
(``--show-provenance``) and by reproducing ``snfpy`` numerically to 1e-12.

**Recovery.** Provenance does not settle which update is *better*, so the second
half plants known structure and measures adjusted Rand index against it, over
replicates, at several signal strengths, plus the null and clonal controls and
the Klebsiella cohort.

Usage::

    python benchmarks/snf_update_benchmark.py --show-provenance
    python benchmarks/snf_update_benchmark.py --quick
    python benchmarks/snf_update_benchmark.py
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from amr_clonalshare.core import spectral_from_similarity  # noqa: E402
from amr_clonalshare.fusion import (  # noqa: E402
    _cross_diffuse,
    _rowstochastic,
    default_K,
    snf_fuse,
    snf_kernel,
)
from amr_clonalshare.influence import layer_influence  # noqa: E402
from amr_clonalshare.stats import binary_distance_matrix  # noqa: E402
from amr_clonalshare.synthetic import (  # noqa: E402
    synth_cluster_archetypes,
    synth_lineage_cohort,
)

UPDATES = ("renormalise", "add_identity")


# --------------------------------------------------------------------------
# 1. numerical agreement with snfpy
# --------------------------------------------------------------------------
def snfpy_agreement(n: int = 60, m: int = 3, K: int = 12, T: int = 20,
                    n_seeds: int = 5) -> dict:
    """Reproduce ``snf.snf`` with ``_cross_diffuse(update="add_identity")``.

    ``snfpy`` sparsifies with a percentile cutoff rather than an exact top-K, so
    we hand our loop *its* ``S`` and its ``P0``. What is being tested is the
    update rule, which is the thing in dispute; the sparsifier is a separate
    documented choice (see ``fusion.sparsify``).
    """
    try:
        import snf.compute as snfc
        from snf.compute import _find_dominate_set, snf as snfpy_snf
    except ImportError:
        return {"available": False,
                "note": "pip install snfpy to run this check"}

    # snfpy 0.2.2 calls sklearn's check_array with force_all_finite=, renamed to
    # ensure_all_finite= in scikit-learn 1.6. Shim it rather than pin sklearn.
    _orig_check_array = snfc.check_array

    def _check_array(a, **kw):
        if "force_all_finite" in kw:
            try:
                return _orig_check_array(a, **kw)
            except TypeError:
                kw["ensure_all_finite"] = kw.pop("force_all_finite")
        return _orig_check_array(a, **kw)

    snfc.check_array = _check_array

    diffs = []
    for seed in range(n_seeds):
        rng = np.random.default_rng(seed)
        affs = []
        for _ in range(m):
            A = rng.random((n, n))
            A = (A + A.T) / 2.0
            np.fill_diagonal(A, 0.0)
            affs.append(A)

        theirs = snfpy_snf(*[a.copy() for a in affs], K=K, t=T, alpha=1.0)

        idx = np.arange(n)
        P = [_rowstochastic(a) for a in affs]
        P = [(p + p.T) / 2.0 for p in P]
        S = [_find_dominate_set(p, K) for p in P]
        P = _cross_diffuse(P, S, T=T, alpha=1.0, update="add_identity", idx=idx)
        ours = sum(P) / m
        ours = _rowstochastic(ours)
        ours = (ours + ours.T + np.eye(n)) / 2.0

        denom = max(np.abs(theirs).max(), 1e-300)
        diffs.append(float(np.abs(ours - theirs).max() / denom))

    return {"available": True, "n": n, "m": m, "K": K, "T": T,
            "n_seeds": n_seeds,
            "max_relative_difference": max(diffs),
            "agrees": max(diffs) < 1e-12}


# --------------------------------------------------------------------------
# 2. recovery of planted structure
# --------------------------------------------------------------------------
def _binary(df):
    return np.array(df.to_numpy(), dtype=bool, copy=True)


def _kernel(X, K):
    D = binary_distance_matrix(np.asarray(X), metric="jaccard",
                               undefined_pair="identical")
    return snf_kernel(D, mu=0.5, K=K)


def _fuse_fn(update, K):
    alpha = 0.5 if update == "renormalise" else 1.0
    return lambda Ws: snf_fuse(Ws, K=K, T=20, alpha=alpha, update=update)


def _profile(mats, names, labels_true, k, update, K):
    """Fuse under one update rule and score the partition it produces."""
    kern = lambda X: _kernel(X, K)          # noqa: E731
    fuse = _fuse_fn(update, K)
    Ws = [kern(X) for X in mats]
    labels = spectral_from_similarity(fuse(Ws), k)
    solo = [spectral_from_similarity(Wl, k) for Wl in Ws]
    inf = layer_influence(mats, names, kernel_fn=kern, fuse_fn=fuse,
                          cluster_fn=spectral_from_similarity, k=k,
                          fused_labels=labels)
    return {
        "labels": labels,
        "ari_vs_truth": (float(adjusted_rand_score(labels_true, labels))
                         if labels_true is not None else float("nan")),
        "max_ari_vs_single_layer": max(
            float(adjusted_rand_score(labels, s)) for s in solo),
        "best_single_layer_ari_vs_truth": (
            max(float(adjusted_rand_score(labels_true, s)) for s in solo)
            if labels_true is not None else float("nan")),
        "n_eff_layers": float(inf["n_eff"]),
        "regime": str(inf["regime"]),
        "collapse": bool(inf["collapse"]),
        "n_clusters_found": int(len(np.unique(labels))),
    }


def planted_recovery(n_rep: int, overlaps=(0.05, 0.15, 0.25, 0.35),
                     n: int = 150, k_true: int = 3) -> dict:
    """Planted archetypes at four levels of between-cluster overlap."""
    out = {}
    for ov in overlaps:
        per_update = {u: {"ari": [], "n_eff": [], "collapsed": [],
                          "ari_vs_layer": []} for u in UPDATES}
        margin = []
        for rep in range(n_rep):
            amr, vir, truth = synth_cluster_archetypes(
                n=n, k_true=k_true, overlap=ov, seed=1000 + rep)
            K = default_K(n)
            mats = [_binary(amr), _binary(vir)]
            scored = {}
            for u in UPDATES:
                r = _profile(mats, ["amr", "vir"], truth, k_true, u, K)
                scored[u] = r
                per_update[u]["ari"].append(r["ari_vs_truth"])
                per_update[u]["n_eff"].append(r["n_eff_layers"])
                per_update[u]["collapsed"].append(r["collapse"])
                per_update[u]["ari_vs_layer"].append(r["max_ari_vs_single_layer"])
            margin.append(scored["renormalise"]["ari_vs_truth"]
                          - scored["add_identity"]["ari_vs_truth"])
        row = {}
        for u in UPDATES:
            d = per_update[u]
            row[u] = {
                "mean_ari_vs_truth": float(np.mean(d["ari"])),
                "sd_ari_vs_truth": float(np.std(d["ari"], ddof=1)) if n_rep > 1 else 0.0,
                "mean_n_eff_layers": float(np.mean(d["n_eff"])),
                "frac_layer_collapse_flagged": float(np.mean(d["collapsed"])),
                "mean_ari_fused_vs_best_single_layer": float(np.mean(d["ari_vs_layer"])),
            }
        margin = np.asarray(margin)
        row["paired_margin_renormalise_minus_add_identity"] = {
            "mean": float(margin.mean()),
            "sd": float(margin.std(ddof=1)) if n_rep > 1 else 0.0,
            "frac_renormalise_better": float((margin > 1e-12).mean()),
            "frac_tied": float((np.abs(margin) <= 1e-12).mean()),
        }
        out[f"overlap={ov}"] = row
    return {"n_rep": n_rep, "n": n, "k_true": k_true, "by_overlap": out}


def null_control(n_rep: int, n: int = 150) -> dict:
    """No structure at all: the fusion should not manufacture agreement."""
    res = {u: {"ari_between_reps": [], "n_eff": [], "regime": []} for u in UPDATES}
    for rep in range(n_rep):
        rng = np.random.default_rng(5000 + rep)
        mats = [rng.binomial(1, 0.3, size=(n, p)).astype(bool) for p in (20, 40)]
        K = default_K(n)
        for u in UPDATES:
            r = _profile(mats, ["amr", "vir"], None, 3, u, K)
            res[u]["ari_between_reps"].append(r["max_ari_vs_single_layer"])
            res[u]["n_eff"].append(r["n_eff_layers"])
            res[u]["regime"].append(r["regime"])
    return {
        "n_rep": n_rep, "n": n,
        "description": "unstructured binary layers; reported quantity is how "
                       "close the fused partition is to a single layer's, and "
                       "whether the influence diagnostic calls it out",
        **{u: {
            "mean_ari_fused_vs_best_single_layer":
                float(np.mean(res[u]["ari_between_reps"])),
            "mean_n_eff_layers": float(np.mean(res[u]["n_eff"])),
            "regimes": {r: res[u]["regime"].count(r)
                        for r in sorted(set(res[u]["regime"]))},
        } for u in UPDATES},
    }


def clonal_control(n_rep: int) -> dict:
    """All structure is clonal: both updates should track lineage, not archetype."""
    from sklearn.metrics import adjusted_mutual_info_score
    res = {u: {"ami_vs_lineage": [], "n_eff": []} for u in UPDATES}
    for rep in range(n_rep):
        amr, vir, meta = synth_lineage_cohort(n=300, n_lineages=30,
                                              seed=7000 + rep)
        lineage = meta["lineage"].to_numpy()
        n = len(lineage)
        K = default_K(n)
        mats = [_binary(amr), _binary(vir)]
        for u in UPDATES:
            r = _profile(mats, ["amr", "vir"], None, 3, u, K)
            res[u]["ami_vs_lineage"].append(
                float(adjusted_mutual_info_score(lineage, r["labels"])))
            res[u]["n_eff"].append(r["n_eff_layers"])
    return {"n_rep": n_rep,
            **{u: {"mean_ami_vs_lineage": float(np.mean(res[u]["ami_vs_lineage"])),
                   "mean_n_eff_layers": float(np.mean(res[u]["n_eff"]))}
               for u in UPDATES}}


def klebsiella_case() -> dict:
    """The real cohort: does either update collapse onto one layer?"""
    import pandas as pd
    root = pathlib.Path(__file__).resolve().parents[1]
    ddir = root / "examples" / "klebsiella" / "data"
    if not (ddir / "amr.csv").exists():
        return {"available": False}

    def load(name):
        df = pd.read_csv(ddir / name)
        return df.drop(columns=[df.columns[0]])

    amr, vir = load("amr.csv"), load("vir.csv")
    meta = pd.read_csv(ddir / "metadata.csv")
    n = len(amr)
    K = default_K(n)
    mats = [_binary(amr), _binary(vir)]
    Ws = [_kernel(X, K) for X in mats]
    solo = {nm: spectral_from_similarity(W, 2)
            for nm, W in zip(("amr", "vir"), Ws)}
    empty_vir = (vir.to_numpy().sum(axis=1) == 0).astype(int)
    empty_amr = (amr.to_numpy().sum(axis=1) == 0).astype(int)

    out = {"available": True, "n": n, "K": K}
    for u in UPDATES:
        r = _profile(mats, ["amr", "vir"], None, 2, u, K)
        lab = r["labels"]
        out[u] = {
            "ari_vs_amr_alone": float(adjusted_rand_score(lab, solo["amr"])),
            "ari_vs_vir_alone": float(adjusted_rand_score(lab, solo["vir"])),
            "ari_vs_empty_virulence_indicator": float(
                adjusted_rand_score(lab, empty_vir)),
            "ari_vs_empty_amr_indicator": float(
                adjusted_rand_score(lab, empty_amr)),
            "n_eff_layers": r["n_eff_layers"],
            "regime": r["regime"],
            "collapse": r["collapse"],
            "cluster_sizes": sorted(np.bincount(lab).tolist(), reverse=True),
            "mean_virulence_score_by_cluster": [
                float(meta["virulence_score"].to_numpy()[lab == c].mean())
                for c in sorted(np.unique(lab))],
        }
    return out


def show_provenance() -> None:
    print(__doc__.split("Usage::")[0].strip())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--show-provenance", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if args.show_provenance:
        show_provenance()
        return

    n_rep = 5 if args.quick else 20
    res = {}
    print("1. numerical agreement with snfpy ...", flush=True)
    res["snfpy_agreement"] = snfpy_agreement(n_seeds=2 if args.quick else 5)
    print("2. planted recovery ...", flush=True)
    res["planted_recovery"] = planted_recovery(n_rep)
    print("3. null control ...", flush=True)
    res["null_control"] = null_control(max(3, n_rep // 3))
    print("4. clonal control ...", flush=True)
    res["clonal_control"] = clonal_control(max(3, n_rep // 6))
    print("5. Klebsiella cohort ...", flush=True)
    res["klebsiella"] = klebsiella_case()

    out = pathlib.Path(args.out or (pathlib.Path(__file__).parent / "results"
                                    / "snf_update.json"))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=2))
    print(f"\nwrote {out}\n")
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
