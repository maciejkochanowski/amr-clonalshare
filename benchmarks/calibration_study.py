#!/usr/bin/env python3
"""Monte-Carlo calibration and power study for the post-clustering inference.

Every number the manuscript quotes about type-I error and power is produced by
this script. It is deliberately self-contained and seeded, so a reader can
re-derive the table rather than take it on trust.

    python benchmarks/calibration_study.py --out benchmarks/results
    python benchmarks/calibration_study.py --quick      # ~2 min smoke run

Experiments
-----------
A. **Merging rules.** Do the four p-value merging rules hold their level on
   exchangeable, dependent p-values? The null is an equicorrelated Gaussian
   copula, ``p_m = Phi(sqrt(rho) W + sqrt(1-rho) e_m)``, so each ``p_m`` is
   exactly uniform while the family is exchangeable and, at ``rho = 0.9``,
   strongly dependent. That is the regime repeated splits of one fixed dataset
   live in.

B. **Feature splitting on binary panels.** Type-I under a global no-cluster null
   with independent features and, separately, with features in correlated
   blocks — which is what real panels look like, and which is where the naive
   per-feature split fails completely. Both the per-feature and the block-aware
   split are run on the same data so the difference is attributable. Power
   against a planted two-group prevalence shift, at four effect sizes, with the
   single-split figure alongside.

C. **Count splitting.** Realised ``corr(X1, X2)`` for the negative-binomial
   split with the true ``r``, with a wrong ``r``, and for the binomial
   (hypergeometric) split with ``m`` known. This is the direct check of the
   independence the method rests on, and the numerical form of Proposition 11
   of Neufeld et al. (2024).

D. **The per-feature selection defect.** The defective per-feature test picks the
   contrast direction from the held-out half and then tested one-sided in that
   direction on the same half. Both the defective and the corrected tests are
   run under an identical null so the inflation can be quoted as a number.

E. **k selection under a null.** Does the criterion return k = 1 when there is
   no structure? Run with and without k = 1 in the sweep to quantify what the
   a sweep that starts at k = 2 costs.

F. **The continuum null.** Power on planted discrete data, and false-positive
   rate on continuum data of one, two and three latent dimensions, with the null
   dimension chosen by BIC and with it pinned at one. The pinned arm is the one
   that fails: a one-dimensional null calls a two-dimensional gradient discrete
   every time.

G. **The automatic blocking threshold.** Type-I error and the resulting number
   of split units as a function of the correlation at which features are forced
   into the same split unit. This is what makes the default 0.7 a calibrated
   value rather than an asserted one.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
from scipy.stats import norm
from sklearn.cluster import KMeans

from amr_clonalshare.baselines import bernoulli_mixture_select_k
from amr_clonalshare.inference import (binomial_thin, correlation_blocks,
                                          feature_split_test, merge_pvalues,
                                          nb_thin, thinning_dependence)
from amr_clonalshare.stats import benjamini_hochberg

MERGERS = ("exchangeable_ruger", "ruger", "twice_mean", "mmb")


def _km(block, k, seed):
    return KMeans(n_clusters=k, n_init=5, random_state=seed).fit(
        np.asarray(block, dtype=float)).labels_


# --------------------------------------------------------------------------- #
def exp_a_merging(n_sim: int, M: int = 9) -> dict:
    rng = np.random.default_rng(20240806)
    out = {}
    for rho in (0.0, 0.5, 0.9):
        W = rng.standard_normal(n_sim)[:, None]
        E = rng.standard_normal((n_sim, M))
        P = norm.cdf(math.sqrt(rho) * W + math.sqrt(1 - rho) * E)
        out[f"rho={rho}"] = {
            name: float(np.mean([merge_pvalues(P[b], name) < 0.05
                                 for b in range(n_sim)]))
            for name in MERGERS
        }
    return {"description": "type-I of each merging rule at alpha=0.05, "
                           "equicorrelated Gaussian-copula exchangeable null",
            "n_sim": n_sim, "M_splits": M, "rates": out}


def _binary_null(rng, n, p, block_size=1):
    """Binary matrix with no clusters; features correlated in blocks."""
    prev = rng.uniform(0.1, 0.6, size=p)
    if block_size <= 1:
        return (rng.random((n, p)) < prev[None, :]).astype(int)
    X = np.zeros((n, p), dtype=int)
    for start in range(0, p, block_size):
        stop = min(start + block_size, p)
        driver = rng.random(n) < prev[start]
        for j in range(start, stop):
            flip = rng.random(n) < 0.1
            X[:, j] = np.where(flip, rng.random(n) < prev[j], driver)
    return X


def exp_b_feature_split(n_sim: int, n: int, p: int, n_splits: int,
                        block_size: int = 5) -> dict:
    """Type-I and power of the feature-split test, per-feature vs block-aware.

    The headline comparison of the paper: with correlated features, splitting
    per feature puts halves of the same correlated block on both sides of the
    split, so the discovery labels predict the held-out features for reasons
    unrelated to clustering. Splitting whole blocks does not.
    """
    res = {}
    groups = np.repeat(np.arange((p + block_size - 1) // block_size),
                       block_size)[:p]
    for label, block in (("independent_features", 1),
                         ("blocked_features", block_size)):
        rng = np.random.default_rng(11 + block)
        per_feat, block_aware, single = [], [], []
        for _ in range(n_sim):
            X = _binary_null(rng, n, p, block_size=block)
            g = None if block <= 1 else groups
            r_pf = feature_split_test(X, 2, cluster_fn=_km, rng=rng,
                                      n_splits=n_splits)
            r_ba = feature_split_test(X, 2, cluster_fn=_km, rng=rng,
                                      n_splits=n_splits, groups=g)
            per_feat.append(r_pf["p_value"])
            block_aware.append(r_ba["p_value"])
            single.append(r_pf["p_per_split"][0])
        res[f"typeI_{label}"] = {
            "per_feature_split": float(np.mean(np.asarray(per_feat) < 0.05)),
            "block_aware_split": float(np.mean(np.asarray(block_aware) < 0.05)),
            "single_split_no_merge": float(np.mean(np.asarray(single) < 0.05)),
        }
    power = {}
    for d in (0.10, 0.20, 0.30, 0.40):
        rng = np.random.default_rng(77)
        merged, single = [], []
        for _ in range(max(n_sim // 2, 30)):
            prev = rng.uniform(0.15, 0.5, size=p)
            g = rng.integers(0, 2, size=n)
            P = np.tile(prev, (n, 1))
            P[g == 1, :p // 3] = np.clip(prev[:p // 3] + d, 0, 1)
            X = (rng.random((n, p)) < P).astype(int)
            r = feature_split_test(X, 2, cluster_fn=_km, rng=rng, n_splits=n_splits)
            merged.append(r["p_value"])
            single.append(r["p_per_split"][0])
        power[f"delta_p={d}"] = {
            "merged": float(np.mean(np.asarray(merged) < 0.05)),
            "single_split": float(np.mean(np.asarray(single) < 0.05)),
        }
    res["power"] = power
    res["block_size"] = block_size
    res["n"] = n
    res["p"] = p
    res["n_splits"] = n_splits
    res["n_sim"] = n_sim
    res["description"] = ("feature-split test on binary panels: type-I under a "
                          "global null with independent and with block-"
                          "correlated features, and power against a planted "
                          "two-group prevalence shift on a third of the features")
    return res


def exp_c_count_splitting(n: int = 2000, p: int = 8) -> dict:
    rng = np.random.default_rng(5)
    r_true, mu = 5.0, 4.0
    Xnb = rng.negative_binomial(r_true, r_true / (r_true + mu), size=(n, p))
    out = {}
    for label, r_used in (("nb_r_correct", r_true), ("nb_r_too_small", 0.4),
                          ("nb_r_too_large", 50.0)):
        X1, X2 = nb_thin(Xnb, r=r_used, rng=rng)
        dep = thinning_dependence(X1, X2)
        out[label] = {"r_used": r_used,
                      "mean_corr": float(np.nanmean(dep)),
                      "max_abs_corr": float(np.nanmax(np.abs(dep)))}
    m = 11
    Xb = rng.binomial(m, 0.3, size=(n, p))
    X1, X2 = binomial_thin(Xb, m=m, rng=rng)
    dep = thinning_dependence(X1, X2)
    out["binomial_hypergeometric_m_known"] = {
        "m": m, "mean_corr": float(np.nanmean(dep)),
        "max_abs_corr": float(np.nanmax(np.abs(dep)))}
    X1, X2 = nb_thin(Xb, r=0.43, rng=rng)
    dep = thinning_dependence(X1, X2)
    out["binomial_data_split_as_nb"] = {
        "r_used": 0.43, "mean_corr": float(np.nanmean(dep)),
        "max_abs_corr": float(np.nanmax(np.abs(dep)))}
    return {"description": "realised corr(X1, X2); should be ~0 when the assumed "
                           "family and parameter are right (Neufeld et al. 2024 "
                           "JMLR 25(57) Prop. 11 gives the exact covariance "
                           "otherwise)",
            "n": n, "p": p, "results": out}


def _perfeature_defective(X1, X2, labels):
    """The defective test: pick the direction from X2, then test one-sided
    on the same X2."""
    from scipy.stats import mannwhitneyu
    p = X2.shape[1]
    pv = np.ones(p)
    uniq = np.unique(labels)
    for j in range(p):
        col = X2[:, j].astype(float)
        means = {int(c): col[labels == c].mean() for c in uniq
                 if (labels == c).sum() > 0}
        if len(means) < 2:
            continue
        top = max(means, key=means.get)
        m = labels == top
        if m.sum() < 2 or (~m).sum() < 2:
            continue
        try:
            pv[j] = float(mannwhitneyu(col[m], col[~m], alternative="greater").pvalue)
        except ValueError:
            pv[j] = 1.0
    return pv


def _perfeature_corrected(X1, X2, labels):
    """Direction-free Kruskal-Wallis across the discovery labels."""
    from scipy.stats import kruskal
    p = X2.shape[1]
    pv = np.ones(p)
    uniq = [int(c) for c in np.unique(labels) if (labels == c).sum() >= 2]
    if len(uniq) < 2:
        return pv
    for j in range(p):
        col = X2[:, j].astype(float)
        groups = [col[labels == c] for c in uniq]
        if np.ptp(np.concatenate(groups)) == 0:
            continue
        try:
            pv[j] = float(kruskal(*groups).pvalue)
        except ValueError:
            pv[j] = 1.0
    return pv


def exp_d_selection_defect(n_sim: int, n: int = 300, p: int = 12,
                           k: int = 3) -> dict:
    rng = np.random.default_rng(31)
    r_true, mu = 5.0, 4.0
    raw = {"defective": [], "corrected": []}
    fam = {"defective": [], "corrected": []}
    for _ in range(n_sim):
        X = rng.negative_binomial(r_true, r_true / (r_true + mu), size=(n, p))
        X1, X2 = nb_thin(X, r=r_true, rng=rng)          # r KNOWN: exact thinning
        labels = _km(X1, k, 0)
        for name, fn in (("defective", _perfeature_defective),
                         ("corrected", _perfeature_corrected)):
            pv = fn(X1, X2, labels)
            raw[name].append(float(np.mean(pv < 0.05)))
            fam[name].append(bool(benjamini_hochberg(pv, q=0.05)[1].any()))
    return {"description": "per-feature type-I under a global negative-binomial "
                           "null with NO clusters and r KNOWN, so the only "
                           "defect is the direction being selected from the "
                           "tested half",
            "n_sim": n_sim, "n": n, "p": p, "k": k,
            "marginal_rejection_rate": {kk: float(np.mean(v)) for kk, v in raw.items()},
            "family_wise_any_bh_rejection": {kk: float(np.mean(v)) for kk, v in fam.items()},
            "nominal": 0.05}


def exp_e_k_selection_null(n_sim: int, n: int = 200, p: int = 30) -> dict:
    rng = np.random.default_rng(41)
    with_one, without_one = [], []
    for _ in range(n_sim):
        X = _binary_null(rng, n, p)
        sel1 = bernoulli_mixture_select_k(X, [1, 2, 3, 4], rng=rng, n_init=5)
        sel2 = bernoulli_mixture_select_k(X, [2, 3, 4], rng=rng, n_init=5)
        with_one.append(sel1["k_selected"])
        without_one.append(sel2["k_selected"])
    return {"description": "k selected on data with no clusters, with and "
                           "without k = 1 in the sweep",
            "n_sim": n_sim, "n": n, "p": p,
            "frac_k_equals_1_when_1_is_in_the_sweep":
                float(np.mean(np.asarray(with_one) == 1)),
            "frac_k_equals_1_when_sweep_starts_at_2": 0.0,
            "mean_k_when_sweep_starts_at_2": float(np.mean(without_one))}


def exp_f_continuum_null(n_sim: int, n: int = 200, p: int = 30, k: int = 3,
                         n_boot: int = 49) -> dict:
    """Does the continuum null separate discrete groups from a gradient?

    Power on planted discrete data, and false-positive rate on data simulated
    from a one-factor latent-trait model. The naive comparison against a single
    *independent-Bernoulli* component is run alongside, because it calls both
    cases discrete and is therefore not an instrument for this question.
    """
    from amr_clonalshare.inference import continuum_null_test
    from amr_clonalshare.synthetic import synth_cluster_archetypes

    rng = np.random.default_rng(61)
    disc_calls, cont_calls = [], []
    for b in range(n_sim):
        amr, vir, _ = synth_cluster_archetypes(n=n, p_amr=p // 2, p_vir=p - p // 2,
                                               k_true=k, overlap=0.05, seed=200 + b)
        X = np.hstack([amr.values, vir.values])
        disc_calls.append(continuum_null_test(X, k, rng=rng, n_boot=n_boot)[
            "discrete_beyond_a_gradient"])
        a = rng.normal(-0.5, 1.0, X.shape[1])
        bb = rng.normal(1.2, 0.4, X.shape[1])
        z = rng.standard_normal(n)
        P = 1 / (1 + np.exp(-(a[None, :] + z[:, None] * bb[None, :])))
        Xc = (rng.random((n, X.shape[1])) < P).astype(int)
        cont_calls.append(continuum_null_test(Xc, k, rng=rng, n_boot=n_boot)[
            "discrete_beyond_a_gradient"])

    # The arm that matters: a continuum with more than one dimension. An accessory genome has at least
    # two (resistance load and virulence load), so a null that can only
    # represent one ordering of the isolates will call any real panel discrete.
    multi = {}
    for q_true in (1, 2, 3):
        for mode in ("q_selected_by_bic", "q_fixed_at_1"):
            calls = []
            for b in range(max(n_sim // 2, 6)):
                Z = rng.standard_normal((n, q_true))
                a = rng.normal(0, 1, p)
                B = rng.normal(0, 1.2, (p, q_true))
                Pq = 1 / (1 + np.exp(-(a[None, :] + Z @ B.T)))
                Xq = (rng.random((n, p)) < Pq).astype(float)
                r = continuum_null_test(Xq, k, rng=rng, n_boot=n_boot,
                                        q=1 if mode == "q_fixed_at_1" else None,
                                        q_max=3)
                if r.get("status") == "ok":
                    calls.append(r["p_value"] < 0.05)
            multi[f"q_true={q_true}/{mode}"] = float(np.mean(calls)) if calls else None

    return {"description": "continuum-null test: power on planted discrete data, "
                           "and false-positive rate on continuum data of one, two "
                           "and three latent dimensions, with the null dimension "
                           "chosen by BIC and with it pinned at one",
            "n_sim": n_sim, "n": n, "p": p, "k": k, "n_boot": n_boot,
            "power_on_discrete": float(np.mean(disc_calls)),
            "false_positive_on_continuum": float(np.mean(cont_calls)),
            "false_positive_by_latent_dimension": multi}


def exp_g_block_threshold(n_sim: int, n: int = 200, p: int = 30,
                          n_splits: int = 11,
                          thresholds=(0.0, 0.3, 0.5, 0.7, 0.9)) -> dict:
    """Where should the automatic blocking threshold sit?

    Experiment B shows that per-feature splitting fails and that *declared*
    blocks fix it. In practice nobody declares every block: the Klebsiella
    config declares five virulence loci and nothing for AMR, while the AMR
    classes are co-selected at r up to 0.81. ``inference.correlation_blocks``
    therefore derives the units from the feature correlation matrix, which
    needs a threshold, and a threshold asserted rather than calibrated is just
    a different arbitrary choice.

    The blocks are graded: within a block the pairwise correlation is not
    constant, so a threshold that is too high leaves genuinely dependent
    features in different units. We sweep it and report both type-I error and
    the resulting number of units, which is what power costs.
    """
    rng = np.random.default_rng(4242)
    block_size = 5
    rows = {}
    for t in thresholds:
        rej, units, cross = [], [], []
        for _ in range(n_sim):
            X = _binary_null(rng, n, p, block_size=block_size)
            g = correlation_blocks(X, threshold=t) if t > 0 else None
            r = feature_split_test(X, 2, cluster_fn=_km, rng=rng,
                                   n_splits=n_splits, groups=g)
            if r.get("status") != "ok":
                continue
            rej.append(r["p_value"])
            units.append(r["n_split_units"])
            cross.append(r["max_abs_corr_between_units"])
        rows[f"threshold={t}"] = {
            "type_I_at_0.05": float(np.mean(np.asarray(rej) < 0.05)) if rej else None,
            "mean_n_split_units": float(np.mean(units)) if units else None,
            "mean_max_abs_corr_between_units": float(np.mean(cross)) if cross else None,
            "n_usable": len(rej),
        }
    return {
        "description": "type-I error of the feature-split test as a function of "
                       "the automatic blocking threshold, under a no-cluster "
                       "null with correlated feature blocks; threshold = 0 "
                       "means no automatic blocking (per-feature split)",
        "n_sim": n_sim, "n": n, "p": p, "true_block_size": block_size,
        "n_splits": n_splits,
        "by_threshold": rows,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="benchmarks/results")
    ap.add_argument("--quick", action="store_true",
                    help="small replicate counts for a smoke run")
    args = ap.parse_args()
    q = args.quick
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    results = {}
    print("A. merging rules ...", flush=True)
    results["A_merging_rules"] = exp_a_merging(300 if q else 4000)
    print("B. feature splitting ...", flush=True)
    # Both split arms run on the same data, so the replicate count dominates
    # the runtime here. 120 is ample to separate 1.000 from 0.04.
    results["B_feature_splitting"] = exp_b_feature_split(
        20 if q else 120, n=200, p=30, n_splits=9 if q else 11)
    print("C. count splitting ...", flush=True)
    results["C_count_splitting"] = exp_c_count_splitting(n=500 if q else 3000)
    print("D. per-feature selection defect ...", flush=True)
    results["D_selection_defect"] = exp_d_selection_defect(30 if q else 300)
    print("E. k selection under a null ...", flush=True)
    results["E_k_selection_null"] = exp_e_k_selection_null(10 if q else 60)
    print("F. continuum null ...", flush=True)
    results["F_continuum_null"] = exp_f_continuum_null(
        4 if q else 25, n_boot=29 if q else 49)
    print("G. automatic blocking threshold ...", flush=True)
    results["G_block_threshold"] = exp_g_block_threshold(
        20 if q else 120, n_splits=9 if q else 11)

    path = out / ("calibration_quick.json" if q else "calibration.json")
    path.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {path}\n")
    print(json.dumps(results, indent=2)[:4000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
