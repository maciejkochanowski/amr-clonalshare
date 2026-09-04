#!/usr/bin/env python3
"""MOFA against the graph fusion, on the case study and on planted controls.

Why this benchmark is not optional
----------------------------------
The comparative multi-view literature is unanimous that the mandatory question
is whether a multi-view solution beats the best single view (Rappoport & Shamir
2018, *Nucleic Acids Res* 46:10546-10562; Tini et al. 2019, *Brief Bioinform*
20:1269-1279). Answering it with single-layer spectral partitions and a
Bernoulli mixture alone would leave MOFA and iClusterBayes unbenchmarked, and
they are the natural competitors.

They are more than competitors here. MOFA (Argelaguet et al. 2018, *Mol Syst
Biol* 14:e8124, doi:10.15252/msb.20178124) with a Bernoulli likelihood is the
multi-view generalisation of the very model this package uses as its continuum
null: ``inference.fit_latent_trait`` is a single-view, ``q``-factor logistic
latent-trait model. So MOFA answers two questions at once:

1. Is the structure in this cohort better described as a small number of
   continuous factors than as discrete groups? (The continuum null already says
   BIC wants at least three latent dimensions.)
2. How much of each view does each factor explain? Per-view variance
   decomposition is the layer-contribution diagnostic that does **not**
   degenerate at two layers, where ``layer_influence``'s ``delta_loo`` is
   algebraically ``1 - ari_solo`` of the other layer.

iClusterBayes (Mo et al. 2018, *Biostatistics* 19:71-86,
doi:10.1093/biostatistics/kxx017) is R/Bioconductor only and is **not** run
here; that is a stated gap, not a silent one.

    python benchmarks/latent_factor_benchmark.py --quick
    python benchmarks/latent_factor_benchmark.py
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from amr_clonalshare.core import spectral_from_similarity
from amr_clonalshare.fusion import default_K, snf_fuse, snf_kernel
from amr_clonalshare.inference import select_latent_dimension
from amr_clonalshare.stats import binary_distance_matrix
from amr_clonalshare.synthetic import synth_cluster_archetypes


def run_mofa(views: dict, n_factors: int = 10, seed: int = 1,
             iterations: int = 300) -> dict:
    """Fit MOFA with a Bernoulli likelihood; return factors and variance explained."""
    from mofapy2.run.entry_point import entry_point

    names = list(views)
    # mofapy2 wants [view][group] -> (samples, features)
    data = [[np.asarray(views[v], dtype=float)] for v in names]
    ent = entry_point()
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        ent.set_data_options(scale_views=False)
        ent.set_data_matrix(data, likelihoods=["bernoulli"] * len(names),
                            views_names=names, groups_names=["g"])
        ent.set_model_options(factors=n_factors, spikeslab_weights=True,
                              ard_weights=True)
        ent.set_train_options(iter=iterations, convergence_mode="fast",
                              seed=seed, verbose=False, gpu_mode=False)
        ent.build()
        ent.run()
    model = ent.model
    Z = np.asarray(model.nodes["Z"].getExpectation())[:, :]
    W = [np.asarray(w) for w in model.nodes["W"].getExpectation()]

    # Per-view variance explained by each factor. The likelihood is Bernoulli,
    # so `Z @ W.T` is a linear predictor, not a reconstruction of the 0/1 matrix:
    # squaring `Y - Z W^T` directly (which an earlier version of this function
    # did) gives identically zero R^2. The reconstruction is
    # `sigmoid(logit(mean_j) + Z[:, f] * W[j, f])`, scored against the
    # intercept-only model on the observed matrix.
    def _sig(a):
        return 1.0 / (1.0 + np.exp(-np.clip(a, -30, 30)))

    r2 = {}
    for vi, name in enumerate(names):
        Y = np.asarray(views[name], dtype=float)
        mean_j = np.clip(Y.mean(axis=0), 1e-4, 1 - 1e-4)
        base = np.log(mean_j / (1 - mean_j))
        denom = float(((Y - mean_j[None, :]) ** 2).sum())
        per = []
        for f in range(Z.shape[1]):
            pred = _sig(base[None, :] + np.outer(Z[:, f], W[vi][:, f]))
            per.append(float(max(0.0, 1.0 - ((Y - pred) ** 2).sum() / denom))
                       if denom > 0 else 0.0)
        r2[name] = per
    return {"factors": Z, "r2_per_view": r2, "n_factors": int(Z.shape[1])}


def _layers_klebsiella(root: pathlib.Path):
    d = root / "examples" / "klebsiella" / "data"
    amr = pd.read_csv(d / "amr.csv").set_index("Strain_ID")
    vir = pd.read_csv(d / "vir.csv").set_index("Strain_ID")
    groups = {"ybt": ["ybtS", "ybtX", "ybtQ", "ybtP", "ybtA", "irp2", "irp1",
                      "ybtU", "ybtT", "ybtE", "fyuA"],
              "clb": [f"clb{c}" for c in "ABCDEFGHILMNOPQ"],
              "iuc": ["iucA", "iucB", "iucC", "iucD", "iutA"],
              "iro": ["iroB", "iroC", "iroD", "iroN"],
              "rmp": ["rmpA", "rmpD", "rmpC"]}
    V = pd.DataFrame({g: vir[[c for c in m if c in vir.columns]].max(axis=1)
                      for g, m in groups.items()}, index=vir.index)
    keep = [c for c in amr.columns if amr[c].sum() >= 2]
    return {"amr": amr[keep], "vir": V}, list(amr.index)


def _fuse(views, k):
    Ws = []
    for M in views.values():
        X = np.array(M.to_numpy(), dtype=bool, copy=True)
        D = binary_distance_matrix(X, metric="jaccard", undefined_pair="identical")
        Ws.append(snf_kernel(D, mu=0.5, K=default_K(X.shape[0])))
    K = default_K(next(iter(views.values())).shape[0])
    return spectral_from_similarity(snf_fuse(Ws, K=K, T=20), k)


def klebsiella_arm(root: pathlib.Path, n_factors: int, iterations: int) -> dict:
    views, ids = _layers_klebsiella(root)
    fused = _fuse(views, 2)
    mofa = run_mofa(views, n_factors=n_factors, iterations=iterations)
    Z = mofa["factors"]

    # How many factors carry non-trivial variance in at least one view?
    r2 = mofa["r2_per_view"]
    active = [f for f in range(Z.shape[1])
              if max(r2[v][f] for v in r2) > 0.01]
    km = KMeans(n_clusters=2, n_init=10, random_state=0).fit_predict(
        Z[:, active] if active else Z)

    X_all = np.hstack([np.asarray(v, dtype=float) for v in views.values()])
    dim = select_latent_dimension(X_all, q_max=4)

    out = {
        "n_isolates": int(Z.shape[0]),
        "n_factors_requested": n_factors,
        "n_factors_with_r2_above_1pct": len(active),
        "r2_per_view_per_factor": {v: [round(x, 4) for x in r2[v]] for v in r2},
        "total_r2_per_view": {v: round(float(sum(r2[v])), 4) for v in r2},
        "ari_mofa_kmeans_vs_snf_partition": float(adjusted_rand_score(km, fused)),
        "latent_dimension_by_bic": dim,
    }
    # Which view does each active factor belong to? This is the layer-attribution
    # question that delta_loo cannot answer at two layers.
    dominant = []
    for f in active:
        vs = {v: r2[v][f] for v in r2}
        top = max(vs, key=vs.get)
        share = vs[top] / max(sum(vs.values()), 1e-12)
        dominant.append({"factor": f, "dominant_view": top,
                         "share_of_factor_r2": round(float(share), 3),
                         "r2": {v: round(float(x), 4) for v, x in vs.items()}})
    out["factor_attribution"] = dominant
    return out


def planted_arm(n_rep: int, n_factors: int, iterations: int) -> dict:
    """Can MOFA recover planted discrete archetypes, and does SNF beat it?"""
    rows = []
    for rep in range(n_rep):
        amr, vir, truth = synth_cluster_archetypes(n=150, k_true=3,
                                                   overlap=0.05, seed=500 + rep)
        views = {"amr": amr, "vir": vir}
        fused = _fuse(views, 3)
        mofa = run_mofa(views, n_factors=n_factors, iterations=iterations)
        km = KMeans(n_clusters=3, n_init=10, random_state=0).fit_predict(
            mofa["factors"])
        rows.append({"snf": float(adjusted_rand_score(truth, fused)),
                     "mofa_kmeans": float(adjusted_rand_score(truth, km))})
    df = pd.DataFrame(rows)
    return {"n_rep": n_rep,
            "mean_ari_vs_truth": {c: float(df[c].mean()) for c in df.columns},
            "sd_ari_vs_truth": {c: float(df[c].std(ddof=1)) if n_rep > 1 else 0.0
                                for c in df.columns}}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    root = pathlib.Path(__file__).resolve().parents[1]
    iters = 100 if args.quick else 300
    nf = 6 if args.quick else 10

    res = {"note": "iClusterBayes (Mo et al. 2018) is R/Bioconductor only and "
                   "is not run here; this is a stated gap"}
    print("1. planted control ...", flush=True)
    res["planted"] = planted_arm(3 if args.quick else 10, nf, iters)
    print("2. Klebsiella cohort ...", flush=True)
    res["klebsiella"] = klebsiella_arm(root, nf, iters)

    path = pathlib.Path(args.out or (root / "benchmarks" / "results"
                                     / "latent_factor.json"))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(res, indent=2, default=str))
    print(f"\nwrote {path}\n")
    print(json.dumps(res, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
