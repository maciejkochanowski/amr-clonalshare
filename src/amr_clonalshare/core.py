"""core.py — the amr-clonalshare pipeline.

The analysis, in order:

1. load the configured layers and align them on a shared strain index;
2. gate features (prevalence and/or absolute count, with a protected list) and
   optionally collapse co-inherited feature groups to locus level;
3. per-layer distances - Jaccard family for ``wide_binary`` layers, a matching
   coefficient for ``one_hot`` layers - with the empty-union convention stated
   explicitly;
4. self-tuning affinities and cross-diffusion fusion (:mod:`.fusion`);
5. choose k (MDL gain, prediction strength, gap statistic, or Bernoulli-mixture
   BIC);
6. spectral clustering of the fused affinity, with a Monti consensus matrix for
   per-isolate confidence;
7. **diagnostics that decide whether the result means anything**: layer
   influence and fusion collapse (:mod:`.influence`), the uninformative-row
   stratum, feature-group collinearity, K-nearest-neighbour tie degeneracy, and
   lineage confounding (:mod:`.lineage`);
8. **baselines** the result must beat to be worth its complexity
   (:mod:`.baselines`);
9. descriptive archetype profiles (effect sizes; these are computed on the data
   the partition came from and are labelled as descriptive, not inferential);
10. **valid post-clustering inference** by multi-split data thinning
    (:mod:`.tva`), with :mod:`.archephy` as the small-cohort, phylogeny-aware
    alternative.

Every stochastic stage draws from its own generator spawned off the master
seed, so adding or removing a diagnostic cannot perturb the clustering.
"""
from __future__ import annotations

import math
from itertools import combinations
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.sparse.csgraph import minimum_spanning_tree
from sklearn.cluster import SpectralClustering
from sklearn.metrics import (
    adjusted_mutual_info_score,
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score,
)

from . import attribution as _attribution
from . import baselines as _baselines
from . import censored as _censored
from . import evalues as _evalues
from . import influence as _influence
from . import clonality as _clonality
from . import lineage as _lineage
from .config import Config
from .fusion import default_K, knn_tie_diagnostics, snf_fuse, snf_kernel
from .io import load_dataset
from .small_n import gap_statistic_k
from .stats import (
    benjamini_hochberg,
    binary_distance_matrix,
    effective_dimension,
    fisher_exact_p,
    permutation_pvalue,
    uninformative_rows,
)
from .synthetic import synth_cluster_archetypes

__all__ = ["run", "validate", "spectral_from_similarity"]

# Number of RNG streams spawned from the master seed, in fixed order.
_STREAMS = ("consensus", "mdl_null", "bootstrap", "tva", "influence",
            "recoverability", "lineage", "baselines")


def _spawn_rngs(seed: int) -> Dict[str, np.random.Generator]:
    children = np.random.SeedSequence(seed).spawn(len(_STREAMS))
    return {name: np.random.default_rng(ss) for name, ss in zip(_STREAMS, children)}


# =============================================================================
# Clustering primitives
# =============================================================================
def spectral_from_similarity(W: np.ndarray, k: int, *,
                             random_state: int = 0) -> np.ndarray:
    """Spectral clustering of a symmetric affinity matrix.

    ``k = 1`` returns a single cluster, so that "no archetypes" is a
    representable answer rather than something the code cannot express.
    """
    n = W.shape[0]
    if k <= 1:
        return np.zeros(n, dtype=int)
    W_sym = np.maximum(W, W.T)
    np.fill_diagonal(W_sym, 0)
    k = int(max(2, min(k, n - 1)))
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        sc = SpectralClustering(n_clusters=k, affinity="precomputed",
                                random_state=random_state,
                                assign_labels="discretize")
        return sc.fit_predict(W_sym)


def consensus_matrix(W: np.ndarray, k: int, rng: np.random.Generator,
                     B: int = 200, rho: float = 0.8) -> np.ndarray:
    """Monti consensus matrix ``M(i,j) = sum_b I_b(i,j) / sum_b C_b(i,j)``.

    Monti et al. (2003), Machine Learning 52:91-118,
    doi:10.1023/A:1023949509487. ``B`` subsamples at fraction ``rho``, spectral
    base clusterer.
    """
    n = W.shape[0]
    I = np.zeros((n, n))
    C = np.zeros((n, n))
    if k <= 1:
        return np.ones((n, n))
    n_sub = max(int(rho * n), k + 1)
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        for b in range(B):
            sub = rng.choice(n, size=n_sub, replace=False)
            W_sub = W[np.ix_(sub, sub)]
            try:
                labels = spectral_from_similarity(W_sub, k, random_state=b)
            except Exception:
                continue
            co = (labels[:, None] == labels[None, :]).astype(np.float64)
            C[np.ix_(sub, sub)] += 1.0
            I[np.ix_(sub, sub)] += co
    with np.errstate(invalid="ignore", divide="ignore"):
        M = np.where(C > 0, I / C, 0.0)
    np.fill_diagonal(M, 1.0)
    return (M + M.T) / 2


# =============================================================================
# Minimum description length
# =============================================================================
def xlog2x(p: float) -> float:
    """Binary entropy, with ``0 log 0 = 0``."""
    return 0.0 if p <= 0 or p >= 1 else -(p * math.log2(p) + (1 - p) * math.log2(1 - p))


def mdl_null(X: np.ndarray, *, complexity: str = "rissanen") -> float:
    """Code length of the i.i.d.-per-feature null model.

    ``sum_f n H(p_f)`` plus, unless ``complexity="legacy"``, the null model's
    own parameter cost ``(p/2) log2 n``. Charging the clustered model for its
    parameters while charging the null nothing is not a code-length comparison,
    which is the trap this function exists to avoid.
    """
    n, p = X.shape
    base = float(sum(n * xlog2x(float(pv)) for pv in X.mean(axis=0)))
    if complexity == "legacy":
        return base
    return base + (p / 2.0) * math.log2(max(n, 2))


def mdl_archetype(X: np.ndarray, labels: np.ndarray, *,
                  complexity: str = "rissanen") -> dict:
    """Two-part MDL decomposition of a partition of a binary matrix.

    ``L_assign = n log2 k`` encodes the labels; ``L_profiles`` encodes the
    per-cluster Bernoulli parameters; ``L_residuals`` is the Bernoulli code
    length of the data under the per-cluster Krichevsky-Trofimov estimate
    ``theta = (sum x + 0.5) / (n_c + 1)``.

    ``complexity="rissanen"`` (default) uses the standard parametric complexity
    ``(d/2) log2 n_c`` for ``d = p`` free parameters per cluster (Rissanen 1983,
    Ann Statist 11:416-431, doi:10.1214/aos/1176346150; Grunwald 2007, *The
    Minimum Description Length Principle*, MIT Press). ``complexity="legacy"``
    reproduces the unpenalised term ``p log2(n_c + 1)``, which is twice the
    standard one and therefore biases selection towards smaller k; it is
    available so that the effect on reported k can be measured.
    """
    n, p = X.shape
    uniq = np.unique(labels)
    k = len(uniq)
    L_assign = n * math.log2(max(k, 1))
    sizes = [int((labels == c).sum()) for c in uniq]
    if complexity == "legacy":
        L_profiles = sum(p * math.log2(nc + 1) for nc in sizes)
    else:
        L_profiles = sum((p / 2.0) * math.log2(max(nc, 2)) for nc in sizes)
    L_res = 0.0
    for c, nc in zip(uniq, sizes):
        mask = labels == c
        theta = (X[mask].sum(axis=0) + 0.5) / (nc + 1.0)
        theta = np.clip(theta, 1e-12, 1 - 1e-12)
        Xc = X[mask].astype(float)
        L_res += float(-(Xc * np.log2(theta) + (1 - Xc) * np.log2(1 - theta)).sum())
    L_null = mdl_null(X, complexity=complexity)
    return {"L_null": L_null, "L_assign": L_assign, "L_profiles": L_profiles,
            "L_residuals": L_res, "L_archetype": L_assign + L_profiles + L_res,
            "gain": L_null - (L_assign + L_profiles + L_res), "k": k,
            "complexity": complexity}


def mdl_null_calibration(X: np.ndarray, k: int, W_fused: np.ndarray,
                         rng: np.random.Generator,
                         refuse: Callable[[np.ndarray], np.ndarray],
                         n_perm: int = 200,
                         complexity: str = "rissanen") -> dict:
    """Permutation test for the MDL gain.

    The null permutes each feature column independently, destroying row
    structure while preserving marginals, and re-runs the whole
    distance -> fusion -> spectral -> MDL chain. The p-value is Phipson-Smyth
    corrected, so it can never be reported as zero, and the attainable minimum
    ``1 / (n_perm + 1)`` is reported alongside it: with a small ``n_perm`` the
    test simply has no resolution, and saying so is better than printing 0.0.
    """
    labels_obs = spectral_from_similarity(W_fused, k)
    gain_obs = mdl_archetype(X, labels_obs, complexity=complexity)["gain"]
    null_gains = []
    for _ in range(int(n_perm)):
        X_perm = X.copy()
        for f in range(X.shape[1]):
            rng.shuffle(X_perm[:, f])
        labels_perm = spectral_from_similarity(refuse(X_perm), k)
        null_gains.append(mdl_archetype(X_perm, labels_perm,
                                        complexity=complexity)["gain"])
    null_gains = np.asarray(null_gains, dtype=float)
    p_val = permutation_pvalue(null_gains, gain_obs, tail="greater")
    return {"k": int(k), "gain_obs": float(gain_obs),
            "null_mean": float(null_gains.mean()) if null_gains.size else float("nan"),
            "null_std": float(null_gains.std()) if null_gains.size else float("nan"),
            "n_perm": int(n_perm),
            "p_MDL": p_val,
            "p_MDL_minimum_attainable": float(1.0 / (int(n_perm) + 1)),
            "reject_H0": bool(p_val < 0.05),
            "note": "Phipson & Smyth (2010) corrected permutation p-value; "
                    "cannot be zero, and is bounded below by 1/(n_perm+1)"}


def mdl_per_archetype(X: np.ndarray, labels: np.ndarray) -> pd.DataFrame:
    """Per-archetype share of the description-length gain."""
    rows = []
    for c in np.unique(labels):
        mask = labels == c
        n_c = int(mask.sum())
        Xc = X[mask].astype(float)
        theta = np.clip((X[mask].sum(axis=0) + 0.5) / (n_c + 1.0), 1e-9, 1 - 1e-9)
        l_res = float(-(Xc * np.log2(theta) + (1 - Xc) * np.log2(1 - theta)).sum())
        l_null_c = float(sum(n_c * xlog2x(float(pv)) for pv in Xc.mean(axis=0)))
        rows.append({"cluster": int(c), "n_c": n_c,
                     "L_null_restricted": l_null_c, "L_residuals": l_res,
                     "MDL_contribution": l_null_c - l_res})
    return pd.DataFrame(rows)


# =============================================================================
# Consensus-matrix summaries
# =============================================================================
def pac_metric(M: np.ndarray, lo: float = 0.1, hi: float = 0.9) -> float:
    """Proportion of ambiguous clustering: consensus entries strictly in (lo, hi)."""
    iu = np.triu_indices_from(M, k=1)
    e = M[iu]
    return float(((e > lo) & (e < hi)).mean()) if e.size else 0.0


def cdf_area(M: np.ndarray) -> float:
    """Area under the empirical CDF of the consensus entries."""
    iu = np.triu_indices_from(M, k=1)
    entries = np.sort(M[iu])
    if entries.size == 0:
        return 0.0
    x = np.linspace(0, 1, entries.size)
    trapz = getattr(np, "trapezoid", None) or np.trapz
    return float(trapz(x, entries))


def delta_cdf_area(M_k: np.ndarray, M_km1: np.ndarray) -> float:
    """Relative gain in consensus CDF area from k-1 to k."""
    a_km1 = cdf_area(M_km1)
    return 0.0 if a_km1 == 0 else (cdf_area(M_k) - a_km1) / a_km1


def assignment_confidence(M: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """``conf(i)`` = mean in-cluster consensus minus best out-of-cluster mean."""
    n = M.shape[0]
    conf = np.zeros(n)
    uniq = np.unique(labels)
    idx = np.arange(n)
    for i in range(n):
        c_i = labels[i]
        in_mask = (labels == c_i) & (idx != i)
        in_avg = M[i, in_mask].mean() if in_mask.any() else 0.0
        best_out = 0.0
        for c in uniq:
            if c == c_i:
                continue
            mask = labels == c
            if mask.any():
                best_out = max(best_out, M[i, mask].mean())
        conf[i] = in_avg - best_out
    return conf


def label_stability(conf: np.ndarray, tau_robust: float = 0.4,
                    tau_fragile: float = 0.1) -> List[str]:
    """Map confidence to {robust, intermediate, fragile}."""
    return ["robust" if c >= tau_robust else "fragile" if c <= tau_fragile
            else "intermediate" for c in conf]


# =============================================================================
# Descriptive archetype profiles (NOT inference - see .tva for that)
# =============================================================================
def cohens_h(p1: float, p2: float) -> float:
    """``2 asin sqrt(p1) - 2 asin sqrt(p2)`` (Cohen 1988)."""
    return (2 * math.asin(math.sqrt(max(0.0, min(1.0, p1))))
            - 2 * math.asin(math.sqrt(max(0.0, min(1.0, p2)))))


def archetype_profiles(X: pd.DataFrame, labels: np.ndarray,
                       rng: np.random.Generator, *,
                       q_fdr: float = 0.05, n_boot: int = 500,
                       min_abs_h: float = 0.5, min_support: int = 3
                       ) -> pd.DataFrame:
    """Per-(cluster, feature) effect sizes with bootstrap intervals.

    **These are descriptive, not inferential.** The partition was estimated from
    the same matrix, so the Fisher p-values and bootstrap intervals below do not
    have their nominal frequentist meaning under the null of no clusters (Gao,
    Bien & Witten 2024, JASA 119:332-342, doi:10.1080/01621459.2022.2116331).
    They rank features by effect size; :func:`amr_clonalshare.tva.tva_report`
    supplies the valid test. The column names carry the ``descriptive_`` prefix
    in the output JSON for that reason.
    """
    rows = []
    if len(np.unique(labels)) < 2:
        # One cluster: there is no "rest" to contrast against.
        return pd.DataFrame(columns=[
            "cluster", "feature", "p_c", "p_rest", "delta_p", "delta_p_ci_lo",
            "delta_p_ci_hi", "cohens_h", "log_OR", "fisher_p", "support",
            "fisher_p_bh", "bh_reject", "is_defining_descriptive"])
    for c in np.unique(labels):
        mask = labels == c
        n_c = int(mask.sum())
        X_in, X_out = X[mask], X[~mask]
        n_out = len(X_out)
        for feat in X.columns:
            x_in = X_in[feat].to_numpy()
            x_out = X_out[feat].to_numpy()
            n11 = int((x_in == 1).sum())
            n10 = int((x_in == 0).sum())
            n01 = int((x_out == 1).sum())
            n00 = int((x_out == 0).sum())
            p_c = n11 / max(n_c, 1)
            p_r = n01 / max(n_out, 1)
            a, b, cc, dd = n11 + 0.5, n10 + 0.5, n01 + 0.5, n00 + 0.5
            boot = np.empty(n_boot)
            for t in range(n_boot):
                bi = rng.integers(0, max(len(x_in), 1), len(x_in))
                bo = rng.integers(0, max(len(x_out), 1), len(x_out))
                boot[t] = (x_in[bi].mean() if len(bi) else 0.0) - \
                          (x_out[bo].mean() if len(bo) else 0.0)
            ci_lo, ci_hi = np.percentile(boot, [2.5, 97.5])
            rows.append({
                "cluster": int(c), "feature": feat, "p_c": p_c, "p_rest": p_r,
                "delta_p": p_c - p_r, "delta_p_ci_lo": float(ci_lo),
                "delta_p_ci_hi": float(ci_hi), "cohens_h": cohens_h(p_c, p_r),
                "log_OR": math.log(a * dd / (b * cc)),
                "fisher_p": fisher_exact_p(n11, n10, n01, n00), "support": n11,
            })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["fisher_p_bh"] = float("nan")
    df["bh_reject"] = False
    for c in df["cluster"].unique():
        mask = df["cluster"] == c
        adj, rej = benjamini_hochberg(df.loc[mask, "fisher_p"].to_numpy(), q=q_fdr)
        df.loc[mask, "fisher_p_bh"] = adj
        df.loc[mask, "bh_reject"] = rej
    df["is_defining_descriptive"] = ((df["cohens_h"].abs() >= min_abs_h)
                                     & df["bh_reject"]
                                     & (df["support"] >= min_support))
    return df


# =============================================================================
# Fragility, trajectory, recoverability
# =============================================================================
def _bernoulli_profiles(X: np.ndarray, labels: np.ndarray
                        ) -> Tuple[np.ndarray, np.ndarray]:
    uniq = np.unique(labels)
    theta = np.zeros((len(uniq), X.shape[1]))
    for i, c in enumerate(uniq):
        mask = labels == c
        theta[i] = (X[mask].sum(axis=0) + 0.5) / (mask.sum() + 1.0)
    return np.clip(theta, 1e-9, 1 - 1e-9), uniq


def fragility_analysis(X: np.ndarray, labels: np.ndarray) -> dict:
    """Sensitivity of each assignment to a single feature call being wrong.

    For each isolate and feature, flip that feature and re-assign by maximum
    Bernoulli log-likelihood under the per-cluster profiles; ``F_isolate`` is
    the fraction of single flips that move the isolate.
    """
    n, p = X.shape
    theta, uniq = _bernoulli_profiles(X, labels)
    log_theta, log_not = np.log(theta), np.log(1 - theta)
    Xf = X.astype(np.float64)
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        LL = Xf @ log_theta.T + (1 - Xf) @ log_not.T
    baseline = LL.argmax(axis=1).astype(int)
    delta = log_theta - log_not
    sign_if = 1.0 - 2.0 * Xf
    best_c = np.full((n, p), -1, dtype=int)
    best_ll = np.full((n, p), -np.inf)
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        for c in range(len(uniq)):
            ll_if = LL[:, c:c + 1] + sign_if * delta[c, :]
            improve = ll_if > best_ll
            best_c = np.where(improve, c, best_c)
            best_ll = np.where(improve, ll_if, best_ll)
    boundary = (best_c != baseline[:, None]).astype(int)
    k = len(uniq)
    B = np.zeros((k, k), dtype=int)
    for c_from in range(k):
        m = baseline == c_from
        if not m.any():
            continue
        flips = best_c[m]
        for c_to in range(k):
            if c_to != c_from:
                B[c_from, c_to] = int((flips == c_to).sum())
    return {"F_isolate": boundary.sum(axis=1) / p,
            "F_feature": boundary.sum(axis=0) / n,
            "boundary_matrix": B,
            "baseline_assignment": baseline, "theta": theta}


def tati_graph(X: np.ndarray, labels: np.ndarray,
               feature_names: Optional[list] = None) -> dict:
    """Trait-acquisition trajectory: MST over binarised cluster centroids.

    Edges are oriented from lower to higher trait count. This is a descriptive
    summary of how the discovered profiles differ, **not** an evolutionary
    inference: it uses no phylogeny and no ancestral-state reconstruction. For a
    tree-aware statement about whether a trait combination arose repeatedly, use
    :mod:`amr_clonalshare.archephy`.
    """
    uniq = np.unique(labels)
    k, p = len(uniq), X.shape[1]
    centroids = np.zeros((k, p), dtype=int)
    for i, c in enumerate(uniq):
        centroids[i] = (X[labels == c].mean(axis=0) > 0.5).astype(int)
    D = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            D[i, j] = (centroids[i] != centroids[j]).sum()
    mst = minimum_spanning_tree(D).toarray()
    edges = []
    for i in range(k):
        for j in range(k):
            if mst[i, j] > 0:
                src, dst = (i, j) if centroids[i].sum() <= centroids[j].sum() else (j, i)
                gained = [feature_names[f] if feature_names else f for f in range(p)
                          if centroids[src][f] == 0 and centroids[dst][f] == 1]
                lost = [feature_names[f] if feature_names else f for f in range(p)
                        if centroids[src][f] == 1 and centroids[dst][f] == 0]
                edges.append({"edge_from": int(uniq[src]), "edge_to": int(uniq[dst]),
                              "hamming_distance": int(mst[i, j]),
                              "n_gained": len(gained), "n_lost": len(lost),
                              "traits_gained": gained[:10], "traits_lost": lost[:10]})
    return {"centroids": centroids, "distance_matrix": D, "edges": edges}


def label_recoverability(X: pd.DataFrame, labels: np.ndarray,
                         rng: np.random.Generator, *,
                         frac_observed: float = 0.5, n_demo: int = 20) -> dict:
    """How well a partial feature vector recovers the label the pipeline assigned.

    **This is not predictive validation.** The target is the pipeline's own
    output, so a high score says the clusters are internally coherent under a
    naive-Bayes read-out - not that they predict anything external. A value near
    1.0 is the expected, uninformative result. Use
    :func:`amr_clonalshare.baselines.external_agreement` against a real
    external labelling for the question this metric is often mistaken for.
    """
    n, p = X.shape
    n_obs = max(1, int(frac_observed * p))
    ids = rng.choice(n, size=min(n_demo, n), replace=False)
    rows = []
    for i in ids:
        X_train = np.delete(X.values, i, axis=0)
        labels_train = np.delete(labels, i)
        theta, uniq = _bernoulli_profiles(X_train, labels_train)
        pi = np.array([(labels_train == c).mean() for c in uniq])
        obs = np.zeros(p, dtype=bool)
        obs[rng.choice(p, size=n_obs, replace=False)] = True
        x = X.values[i]
        ll = ((x[obs] * np.log(theta[:, obs])
               + (1 - x[obs]) * np.log(1 - theta[:, obs])).sum(axis=1)
              + np.log(pi + 1e-12))
        rows.append({"isolate_id": X.index[i],
                     "assigned_cluster": int(labels[i]),
                     "recovered_cluster": int(uniq[int(np.argmax(ll))]),
                     "correct": int(uniq[int(np.argmax(ll))] == labels[i])})
    df = pd.DataFrame(rows)
    return {
        "records": df.to_dict("records"),
        "internal_label_recoverability": float(df["correct"].mean()) if len(df) else float("nan"),
        "n_evaluated": int(len(df)),
        "note": "recovery of the pipeline's own labels from half the features; "
                "internal coherence, not external predictive validity",
    }


# =============================================================================
# Cross-layer agreement
# =============================================================================
def hungarian_align(labels_a: np.ndarray, labels_b: np.ndarray) -> np.ndarray:
    """Relabel ``labels_a`` to best match ``labels_b``."""
    ua, ub = np.unique(labels_a), np.unique(labels_b)
    K = max(len(ua), len(ub))
    cost = np.zeros((K, K))
    for i, ca in enumerate(ua):
        for j, cb in enumerate(ub):
            cost[i, j] = -int(((labels_a == ca) & (labels_b == cb)).sum())
    row, col = linear_sum_assignment(cost)
    mapping = {int(ua[r]): int(ub[c]) if c < len(ub) else -1
               for r, c in zip(row, col) if r < len(ua)}
    return np.array([mapping.get(int(l), -1) for l in labels_a])


def layer_agreement(layer_labels: List[np.ndarray], layer_names: List[str]
                    ) -> dict:
    """Pairwise Hungarian-aligned agreement between per-layer partitions.

    ``conflict_per_isolate`` averages over **layer-layer** pairs only; the fused
    partition, when supplied as the first element, is reported separately as
    ``fused_vs_layer`` rather than being folded into the mean.
    """
    m = len(layer_labels)
    n = len(layer_labels[0])
    aligned = [layer_labels[0]] + [hungarian_align(l, layer_labels[0])
                                   for l in layer_labels[1:]]
    pair = {}
    for l1, l2 in combinations(range(m), 2):
        pair[(l1, l2)] = float((aligned[l1] == aligned[l2]).mean())
    layer_idx = list(range(1, m))
    conflict = np.zeros(n)
    n_pairs = 0
    for l1, l2 in combinations(layer_idx, 2):
        conflict += (aligned[l1] != aligned[l2]).astype(float)
        n_pairs += 1
    conflict /= max(n_pairs, 1)
    delta = np.zeros((m, n))
    for l in range(m):
        contrib = [(aligned[l] != aligned[lp]).astype(float)
                   for lp in range(m) if lp != l]
        delta[l] = np.mean(contrib, axis=0) if contrib else np.zeros(n)
    return {
        "aligned": aligned,
        "pair_agreement": {f"{layer_names[a]}|{layer_names[b]}": v
                           for (a, b), v in pair.items()},
        "fused_vs_layer": {layer_names[b]: pair[(0, b)] for b in layer_idx},
        "conflict_per_isolate": conflict,
        "conflict_mean_layer_pairs_only": float(conflict.mean()) if n_pairs else float("nan"),
        "n_layer_pairs": int(n_pairs),
        "delta": delta,
        "discordant_layer_per_isolate": np.argmax(delta, axis=0),
    }


def prediction_strength(W_fused: np.ndarray, k: int, *, train_frac: float = 0.5,
                        n_reps: int = 30, seed: int = 42) -> dict:
    """Prediction strength (Tibshirani & Walther 2005, JCGS 14:511-528).

    ``PS(k)`` is the minimum over test clusters of the proportion of test-cluster
    pairs that a train-based nearest-centroid rule co-assigns. The pair counting
    is vectorised rather than looped.
    """
    n = W_fused.shape[0]
    rng = np.random.default_rng(seed)
    scores = []
    for _ in range(n_reps):
        perm = rng.permutation(n)
        n_train = int(train_frac * n)
        tr, te = perm[:n_train], perm[n_train:]
        if len(te) < k + 1 or len(tr) < k + 1:
            continue
        labels_tr = spectral_from_similarity(W_fused[np.ix_(tr, tr)], k)
        labels_te = spectral_from_similarity(W_fused[np.ix_(te, te)], k)
        Wtt = W_fused[np.ix_(te, tr)]
        uniq_tr = np.unique(labels_tr)
        means = np.column_stack([Wtt[:, labels_tr == c].mean(axis=1) for c in uniq_tr])
        clf = uniq_tr[means.argmax(axis=1)]
        per_cluster = []
        for c in np.unique(labels_te):
            idx = np.where(labels_te == c)[0]
            if idx.size < 2:
                continue
            lab = clf[idx]
            _, cnt = np.unique(lab, return_counts=True)
            agree = int((cnt * (cnt - 1) // 2).sum())
            total = idx.size * (idx.size - 1) // 2
            per_cluster.append(agree / total)
        if per_cluster:
            scores.append(min(per_cluster))
    mean_ps = float(np.mean(scores)) if scores else 0.0
    return {"k": int(k), "mean_PS": mean_ps,
            "sd_PS": float(np.std(scores)) if scores else 0.0,
            "n_reps": len(scores), "passes_spec": mean_ps >= 0.8}


# =============================================================================
# Layer preparation
# =============================================================================
def _apply_groups(df: pd.DataFrame, groups: Dict[str, List[str]]) -> pd.DataFrame:
    """Collapse named feature groups to a single presence call per group."""
    used = set()
    cols = {}
    for name, members in groups.items():
        present = [m for m in members if m in df.columns]
        if not present:
            continue
        cols[name] = (df[present].sum(axis=1) > 0).astype(int)
        used.update(present)
    for c in df.columns:
        if c not in used:
            cols[str(c)] = df[c].astype(int)
    return pd.DataFrame(cols, index=df.index)


def _load_gated_layers(cfg: Config, ds) -> Tuple[List[str], List[pd.DataFrame],
                                                 pd.Index, List[dict]]:
    """Load, optionally group-collapse, and gate each configured layer.

    Returns ``(names, frames, strain_ids, per_layer_reports)``. Every feature
    removed by the gate is recorded with its prevalence, so nothing disappears
    silently: in AMR surveillance the dropped rare determinant is often the one
    that mattered.
    """
    gate = cfg.trait_cluster.prevalence_gate
    top_var = cfg.trait_cluster.top_variance
    names: List[str] = []
    frames: List[pd.DataFrame] = []
    reports: List[dict] = []

    for role in cfg.trait_cluster.layers:
        spec = cfg.files[role]
        df = ds.binary(role)
        n_in = df.shape[1]
        if cfg.trait_cluster.collapse_feature_groups and spec.groups:
            df = _apply_groups(df, spec.groups)
        col_mean = df.mean(axis=0)
        col_count = df.sum(axis=0)
        protected = set(spec.protected)
        keep_mask = ((col_mean > gate.lo) & (col_mean < gate.hi)
                     & (col_count >= gate.min_count))
        keep_mask |= pd.Series([c in protected for c in df.columns],
                               index=df.columns)
        # A protected feature is exempt from the prevalence gate, but a column
        # present in zero isolates cannot inform a distance, cannot be tested,
        # and would occupy a slot in the FDR denominator. Report it, never use
        # it.
        zero = col_count == 0
        keep_mask &= ~zero
        dropped = [{"feature": str(c), "prevalence": float(col_mean[c]),
                    "count": int(col_count[c]),
                    "reason": ("present in no isolate" if zero[c]
                               else "below the prevalence gate")}
                   for c in df.columns if not keep_mask[c]]
        df_gated = df.loc[:, keep_mask[keep_mask].index]

        truncated = []
        if top_var is not None and df_gated.shape[1] > top_var:
            # For a binary column var = p(1-p); with every prevalence below 0.5
            # this ranking is exactly a prevalence ranking. Say so rather than
            # calling it an information criterion.
            variances = df_gated.var(axis=0)
            chosen = set(variances.sort_values(ascending=False).index[:top_var])
            truncated = [{"feature": str(c), "prevalence": float(df_gated[c].mean())}
                         for c in df_gated.columns if c not in chosen]
            df_gated = df_gated.loc[:, [c for c in df_gated.columns if c in chosen]]

        M = df_gated.values.astype(int)
        reports.append({
            "layer": role, "kind": spec.kind,
            "n_features_in": int(n_in), "n_features_used": int(df_gated.shape[1]),
            "gated_out": dropped, "truncated_by_top_variance": truncated,
            "n_uninformative_rows": int(uninformative_rows(M).sum()) if M.size else 0,
            "frac_uninformative_rows": (float(uninformative_rows(M).mean())
                                        if M.size else 0.0),
            "effective_dimension": effective_dimension(M) if M.size else 0.0,
            "density": float(M.mean()) if M.size else 0.0,
        })
        names.append(role)
        frames.append(df_gated)
    return names, frames, ds.strain_ids, reports


def _layer_distance(cfg: Config, role: str, M: np.ndarray) -> np.ndarray:
    kind = cfg.files[role].kind
    if kind == "one_hot":
        return binary_distance_matrix(M, metric=cfg.distance.one_hot_metric,
                                      undefined_pair=cfg.distance.undefined_pair)
    return binary_distance_matrix(M, metric=cfg.distance.metric,
                                  undefined_pair=cfg.distance.undefined_pair)


def _assign_to_training_clusters(cfg, layer_names, layer_mats, train, held_out,
                                 labels_train, K):
    """Place each held-out isolate in the training cluster it sits closest to.

    The package has no way to assign a new isolate. `spectral_from_similarity`
    is a fit on a precomputed affinity, and `snf_fuse` returns a diffused
    operator on the graph it was given with no out-of-sample extension. The
    self-tuning kernel does extend, because its bandwidth is a mean distance to
    the K nearest neighbours and that mean is taken here over training isolates
    only, while the binary distances are pairwise and so carry nothing from one
    held-out row into another. What cannot be extended is the cross-diffusion,
    so the affinity below is the fused kernel before diffusion rather than the
    matrix the reported partition was cut from. That gap belongs to the fusion
    rather than to this function, and it is why the attribution reports an
    assignment control beside the refitted figure.

    The rule itself, highest mean affinity to a training cluster, is the one
    `prediction_strength` already uses, so the run does not acquire a second
    notion of what it means for an isolate to be near a cluster.
    """
    n_train = len(train)
    Kk = int(np.clip(K, 1, max(n_train - 1, 1)))
    rows = np.concatenate([train, held_out])
    affinities = []
    for role, M in zip(layer_names, layer_mats):
        D = _layer_distance(cfg, role, M[rows])
        D_tt = D[:n_train, :n_train].copy()
        np.fill_diagonal(D_tt, np.inf)
        part = np.partition(D_tt, Kk - 1, axis=1)[:, :Kk]
        with np.errstate(invalid="ignore"):
            scale_train = np.nan_to_num(
                np.nanmean(np.where(np.isfinite(part), part, np.nan), axis=1),
                nan=0.0)
        D_ht = D[n_train:, :n_train]
        scale_held = np.partition(D_ht, Kk - 1, axis=1)[:, :Kk].mean(axis=1)
        eps = (scale_held[:, None] + scale_train[None, :] + D_ht) / 3.0
        eps = np.where(eps < 1e-12, 1e-12, eps)
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            W = np.exp(-(D_ht ** 2) / (cfg.snf.mu * eps))
        affinities.append(np.nan_to_num(W, nan=0.0, posinf=0.0, neginf=0.0))
    affinity = np.mean(affinities, axis=0)
    present = np.unique(labels_train)
    to_cluster = np.column_stack(
        [affinity[:, labels_train == c].mean(axis=1) for c in present])
    out = np.empty(n_train + len(held_out), dtype=int)
    out[train] = labels_train
    out[held_out] = present[to_cluster.argmax(axis=1)]
    return out


# =============================================================================
# k selection
# =============================================================================
def _select_k(cfg: Config, k_select: str, X_combined: np.ndarray,
              W_fused: np.ndarray, rng: np.random.Generator, *,
              ps_n_reps: int, gap_B: int, mdl_complexity: str,
              null_labeler: Callable[[np.ndarray, int], np.ndarray]) -> dict:
    """Choose the number of archetypes, with ``k = 1`` a reachable answer.

    Every criterion here can return 1. That is not a detail: a cluster-number
    criterion whose sweep starts at 2 cannot ever say "there are no clusters",
    and on pure noise such a sweep reports k = 2 with a *negative* MDL gain,
    that is, while its own criterion says the unclustered code is shorter.

    * ``mdl`` - k with the largest description-length gain, but ``k = 1`` when
      no k achieves a positive gain.
    * ``prediction_strength`` - Tibshirani & Walther (2005): the **largest** k
      whose mean prediction strength reaches the threshold, with ``PS(1) = 1``
      by convention so that k = 1 is selected when no k qualifies. (``argmax PS``,
      the obvious alternative, always returns some k >= 2.)
    * ``gap`` - the 1-standard-error rule over a sweep that includes k = 1.
    * ``bic_mixture`` - BIC over Bernoulli-mixture fits including k = 1.
    """
    k_range = list(cfg.trait_cluster.k_range)
    threshold = cfg.trait_cluster.prediction_strength_threshold
    requested = (cfg.trait_cluster.k_select_primary if k_select == "auto"
                 else k_select)
    diag: Dict[str, object] = {"requested_method": requested,
                               "k_select_arg": k_select,
                               "mdl_complexity": mdl_complexity,
                               "k_range_swept": [1] + k_range}

    def run_mdl():
        sweep = []
        for kk in [1] + k_range:
            info = mdl_archetype(X_combined, spectral_from_similarity(W_fused, kk),
                                 complexity=mdl_complexity)
            sweep.append({"k": kk, "gain": info["gain"], "L_null": info["L_null"],
                          "L_archetype": info["L_archetype"]})
        best = max(sweep, key=lambda r: r["gain"])
        chosen = int(best["k"]) if best["gain"] > 0 else 1
        l_null = float(sweep[0]["L_null"]) if sweep else float("nan")
        frac = float(best["gain"] / l_null) if l_null else float("nan")
        return chosen, {
            "mdl_sweep": sweep,
            "mdl_best_gain": float(best["gain"]),
            "mdl_gain_fraction": frac,
            "mdl_gain_bits_per_isolate": float(best["gain"] / X_combined.shape[0]),
            "mdl_prefers_no_clustering": bool(best["gain"] <= 0),
            "mdl_optimum_at_k_range_boundary": bool(
                int(best["k"]) == max(k_range) and best["gain"] > 0),
            "mdl_note": "gain as a fraction of the i.i.d. code length is the "
                        "interpretable scale; a fraction below ~0.01 is a weak "
                        "preference however large the bit count looks",
        }

    def run_ps():
        sweep = [prediction_strength(W_fused, kk, n_reps=ps_n_reps, seed=42)
                 for kk in k_range]
        qualifying = [r for r in sweep if r["mean_PS"] >= threshold]
        chosen = int(max(r["k"] for r in qualifying)) if qualifying else 1
        best = max(sweep, key=lambda r: r["mean_PS"]) if sweep else {"mean_PS": 0.0}
        return chosen, {"prediction_strength_sweep": sweep,
                        "prediction_strength_threshold": threshold,
                        "best_mean_PS": best["mean_PS"],
                        "rule": "largest k with mean_PS >= threshold "
                                "(Tibshirani & Walther 2005); k=1 otherwise"}, \
               best["mean_PS"]

    def run_gap():
        gap = gap_statistic_k(
            X_combined, k_range,
            observed_labeler=lambda kk: spectral_from_similarity(W_fused, kk),
            null_labeler=null_labeler, B=gap_B, seed=42, include_k1=True)
        return int(gap["k_star"] or 1), {"gap_statistic": gap}

    def run_bic():
        sel = _baselines.bernoulli_mixture_select_k(
            X_combined, [1] + k_range, rng=rng, n_init=8)
        return int(sel["k_selected"]), {"bic_mixture": {"table": sel["table"],
                                                        "criterion": "BIC"}}

    if requested == "gap":
        k, d = run_gap()
        diag.update(d, method_used="gap")
    elif requested == "bic_mixture":
        k, d = run_bic()
        diag.update(d, method_used="bic_mixture")
    elif requested == "prediction_strength":
        k, d, best_ps = run_ps()
        diag.update(d)
        if k == 1 and best_ps < threshold:
            k_gap, dgap = run_gap()
            diag.update(dgap, method_used="gap",
                        fallback_from="prediction_strength")
            k = k_gap
        else:
            diag["method_used"] = "prediction_strength"
    else:
        k, d = run_mdl()
        diag.update(d, method_used="mdl")

    diag["selected_k"] = int(k)
    diag["no_structure"] = bool(k <= 1)
    return diag


# =============================================================================
# Public API
# =============================================================================
def run(cfg: Config, *, results_dir=None, seed: int = 42,
        k_select: str = "auto", **options) -> dict:
    """Run the full pipeline on a configured cohort.

    Parameters
    ----------
    cfg : validated :class:`~amr_clonalshare.config.Config`.
    results_dir : where to write ``archetype_profiles.tsv`` (and, from the CLI,
        ``cluster_result.json``).
    seed : master seed; each stochastic stage gets its own spawned generator.
    k_select : ``auto`` | ``mdl`` | ``prediction_strength`` | ``gap`` |
        ``bic_mixture``.
    options : Monte-Carlo budgets (``consensus_B``, ``n_perm``, ``n_boot``,
        ``ps_n_reps``, ``gap_B``, ``recoverability_n``, ``lineage_n_perm``),
        ``mdl_complexity`` (``"rissanen"`` | ``"legacy"``), ``auto_cap_B``, and
        ``run_baselines``.
    """
    rngs = _spawn_rngs(seed)

    ds = load_dataset(cfg)
    layer_names, layer_frames, strain_ids, layer_reports = _load_gated_layers(cfg, ds)
    n_isolates = 0 if strain_ids is None else len(strain_ids)
    if n_isolates < 4:
        raise ValueError(f"need at least 4 aligned isolates; got {n_isolates}")
    for name, df in zip(layer_names, layer_frames):
        if df.shape[1] == 0:
            raise ValueError(
                f"layer {name!r} has no features left after gating; loosen "
                f"trait_cluster.prevalence_gate or add features to 'protected'")

    small = (not bool(options.get("auto_cap_B", True))) or n_isolates <= 300
    consensus_B = int(options.get("consensus_B", 100 if small else 30))
    n_perm = int(options.get("n_perm", 0))
    n_boot = int(options.get("n_boot", 500 if small else 200))
    ps_n_reps = int(options.get("ps_n_reps", 30 if small else 10))
    gap_B = int(options.get("gap_B", 30 if small else 10))
    recoverability_n = int(options.get("recoverability_n", 20))
    lineage_n_perm = int(options.get("lineage_n_perm", 500))
    mdl_complexity = str(options.get("mdl_complexity", "rissanen"))
    run_baselines = bool(options.get("run_baselines", True))

    X_df = pd.concat(layer_frames, axis=1)
    X_combined = X_df.values.astype(int)
    layer_mats = [df.values.astype(int) for df in layer_frames]
    col_splits = []
    off = 0
    for df in layer_frames:
        col_splits.append(np.arange(off, off + df.shape[1]))
        off += df.shape[1]

    K = cfg.snf.K if cfg.snf.K is not None else default_K(n_isolates)

    def kernel_of(role: str, M: np.ndarray) -> np.ndarray:
        return snf_kernel(_layer_distance(cfg, role, M), mu=cfg.snf.mu, K=K)

    def fuse(Ws: List[np.ndarray]) -> np.ndarray:
        return snf_fuse(Ws, K=K, T=cfg.snf.T, alpha=cfg.snf.alpha,
                        tie_policy=cfg.snf.tie_policy, update=cfg.snf.update)

    W_layers = [kernel_of(r, M) for r, M in zip(layer_names, layer_mats)]
    W_fused = fuse(W_layers)

    for rep, W in zip(layer_reports, W_layers):
        rep["knn_ties"] = knn_tie_diagnostics(W, K)

    def refuse_matrix(X_any: np.ndarray) -> np.ndarray:
        """Rebuild the fused affinity from an arbitrary matrix of the same shape.

        Used for the MDL permutation null and for the gap statistic's reference
        samples, so that null data passes through exactly the pipeline the
        observed data did.
        """
        return fuse([snf_kernel(_layer_distance(cfg, r, X_any[:, cols]),
                                mu=cfg.snf.mu, K=K)
                     for r, cols in zip(layer_names, col_splits)])

    def null_labeler(X_null: np.ndarray, kk: int) -> np.ndarray:
        return spectral_from_similarity(refuse_matrix(X_null), kk)

    ksel = _select_k(cfg, k_select, X_combined, W_fused, rngs["baselines"],
                     ps_n_reps=ps_n_reps, gap_B=gap_B,
                     mdl_complexity=mdl_complexity, null_labeler=null_labeler)
    k = ksel["selected_k"]
    labels = spectral_from_similarity(W_fused, k)

    def refit_partition(train, held_out, labels=None):
        """Rebuild the partition from the training isolates of one fold.

        The clustering half is this run's own chain, distance then affinity
        then fusion then spectral, restricted to rows, so the fold-internal
        partition is the same kind of object the run reports and not an
        imitation of it. k is held at the value the whole cohort selected,
        because letting every fold choose its own would mix a change of
        resolution into a measurement about leakage.

        Passing ``labels`` skips the rebuild and places the held-out isolates
        against the partition already fitted on every isolate. That is the
        assignment control: it separates what the placement rule contributes
        from what rebuilding the clustering contributes, and without it the two
        arrive as one number that cannot be taken apart afterwards.
        """
        if labels is None:
            Ws = [snf_kernel(_layer_distance(cfg, r, M[train]), mu=cfg.snf.mu,
                             K=K)
                  for r, M in zip(layer_names, layer_mats)]
            labels_train = spectral_from_similarity(fuse(Ws), k)
        else:
            labels_train = np.asarray(labels)[train]
        return _assign_to_training_clusters(cfg, layer_names, layer_mats,
                                            train, held_out, labels_train, K)

    # The MDL permutation null re-runs distance -> fusion -> spectral for every
    # permutation, which is O(n_perm * T * m * n^3): at n = 1500 a p-value with
    # any resolution costs hours, and a default of 5 permutations would
    # produce a "p = 0.0" that no number of permutations could justify. It is
    # therefore opt-in (`n_perm=`), and the post-clustering test in
    # `post_clustering_inference` is the default instrument for "is this real".
    if n_perm > 0:
        mdl_calib = mdl_null_calibration(X_combined, k, W_fused, rngs["mdl_null"],
                                         refuse_matrix, n_perm=n_perm,
                                         complexity=mdl_complexity)
    else:
        mdl_calib = {
            "n_perm": 0,
            "p_MDL": None,
            "reject_H0": None,
            "gain_obs": float(mdl_archetype(X_combined, labels,
                                            complexity=mdl_complexity)["gain"]),
            "note": "permutation null not run (n_perm=0). Re-run with "
                    "run(..., n_perm=199) for a p-value; cost is O(n_perm * "
                    "T * m * n^3). The default evidence for real structure is "
                    "result['post_clustering_inference'].",
        }
    mdl_contrib = mdl_per_archetype(X_combined, labels)

    M = consensus_matrix(W_fused, k=k, rng=rngs["consensus"], B=consensus_B, rho=0.8)
    conf = assignment_confidence(M, labels)
    stab = label_stability(conf, cfg.trait_cluster.stability.tau_robust,
                           cfg.trait_cluster.stability.tau_fragile)

    profiles = archetype_profiles(
        X_df, labels, rngs["bootstrap"],
        q_fdr=cfg.trait_cluster.defining.q_fdr, n_boot=n_boot,
        min_abs_h=cfg.trait_cluster.defining.min_abs_cohens_h,
        min_support=cfg.trait_cluster.defining.min_support)

    # ---- diagnostics that decide whether any of the above means anything ----
    artifacts = _artifact_diagnostics(layer_names, layer_mats, labels, layer_reports)

    infl = None
    if cfg.influence.enabled:
        # The leave-one-layer-out influence costs m extra fusions and is always
        # computed. The permutation test for "does this layer contribute
        # anything beyond its own marginals" costs 2 * n_perm * m fusions, which
        # is tens of minutes at n = 1500; cap it by cohort size rather than
        # letting a config written for a small cohort silently run for an hour,
        # and report that the cap fired.
        infl_perm = cfg.influence.n_perm
        if n_isolates > 500 and infl_perm > 0:
            infl_perm = int(options.get("influence_n_perm", 0))
        infl = _layer_influence(cfg, layer_names, layer_mats, fuse, k, labels,
                                rngs["influence"], K, n_perm=infl_perm)
        infl["n_perm_requested"] = int(cfg.influence.n_perm)
        if infl_perm != cfg.influence.n_perm:
            infl["n_perm_note"] = (
                f"permutation test skipped at n = {n_isolates} (> 500): it costs "
                f"2 * n_perm * n_layers network fusions. Force it with "
                f"run(..., influence_n_perm={cfg.influence.n_perm}).")

    # ---- valid post-clustering inference -----------------------------------
    inference_block = None
    if cfg.tva.enabled:
        inference_block = _run_inference(cfg, layer_names, layer_frames, X_df, k,
                                         rngs["tva"], W_fused, K,
                                         reported_labels=labels)

    frag = fragility_analysis(X_combined, labels)
    tati = tati_graph(X_combined, labels, feature_names=list(X_df.columns))
    recover = label_recoverability(X_df, labels, rngs["recoverability"],
                                   frac_observed=0.5, n_demo=recoverability_n)

    layer_only = [spectral_from_similarity(W, k) for W in W_layers]
    agree = layer_agreement([labels] + layer_only, ["fused"] + layer_names)
    layer_contrib = _layer_contribution(cfg, layer_names, layer_mats, labels, k)

    meta_block = _metadata_diagnostics(cfg, ds, strain_ids, labels,
                                       rngs["lineage"], lineage_n_perm,
                                       X_df=X_df, refit=refit_partition)
    _md = getattr(ds, "metadata", None)
    _ext = None
    if _md is not None and cfg.dataset.external_columns:
        _mdx = _md.reindex(strain_ids)
        _ext = {c: _mdx[c].tolist() for c in cfg.dataset.external_columns
                if c in _mdx.columns}
    pheno_block = _phenotype_validation(cfg, ds, strain_ids, labels, layer_names,
                                        layer_mats, _ext)

    base_block = None
    if run_baselines:
        externals = None
        md = getattr(ds, "metadata", None)
        if md is not None and cfg.dataset.external_columns:
            md = md.reindex(strain_ids)
            externals = {c: md[c].tolist() for c in cfg.dataset.external_columns
                         if c in md.columns}
        base_block = _run_baselines(cfg, layer_names, layer_mats, X_combined,
                                    labels, k, rngs["baselines"], K,
                                    externals=externals)

    if results_dir is not None:
        out_dir = Path(results_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        profiles.to_csv(out_dir / "archetype_profiles.tsv", sep="\t", index=False)

    ids = list(strain_ids)
    delta_names = ["fused"] + layer_names
    # layer_conflict compares an isolate's per-layer cluster assignments. For an
    # isolate carrying none of a layer's features that assignment is the empty
    # stratum, so the comparison measures which empty stratum the isolate is in
    # rather than any biological discordance: with two layers the raw statistic
    # is exactly 1{positive in exactly one layer}, which is 0 both for isolates
    # positive in neither layer and for isolates positive in both. Those are the
    # two groups a user most needs to tell apart, so the statistic is emitted as
    # null unless the isolate is informative in at least two layers.
    positive = np.array([~uninformative_rows(M) for M in layer_mats])
    n_informative = positive.sum(axis=0)
    assignment = [
        {"isolate_id": ids[i], "cluster": int(labels[i]),
         "confidence": float(conf[i]), "stability": stab[i],
         "fragility": float(frag["F_isolate"][i]),
         "positive_layers": [layer_names[j] for j in range(len(layer_names))
                             if positive[j, i]],
         "n_informative_layers": int(n_informative[i]),
         "layer_conflict": (float(agree["conflict_per_isolate"][i])
                            if n_informative[i] >= 2 else None),
         "discordant_layer": (
             delta_names[int(agree["discordant_layer_per_isolate"][i])]
             if n_informative[i] >= 2 else None)}
        for i in range(len(ids))
    ]

    if results_dir is not None:
        # The per-isolate table is the join a surveillance user needs, and it
        # existed only nested inside the JSON. Every metadata column is carried
        # through, not only the lineage and external columns the analysis used.
        atab = pd.DataFrame(assignment)
        atab["positive_layers"] = atab["positive_layers"].map("+".join)
        md_all = getattr(ds, "metadata", None)
        if md_all is not None:
            extra = md_all.reindex(ids)
            extra = extra[[c for c in extra.columns if c not in atab.columns]]
            atab = pd.concat([atab.set_index("isolate_id"), extra],
                             axis=1).rename_axis("isolate_id").reset_index()
        atab.to_csv(Path(results_dir) / "assignment.tsv", sep="\t", index=False)

    return {
        "schema_version": "1.0",
        "seed": int(seed),
        "config": {
            "distance": {"metric": cfg.distance.metric,
                         "undefined_pair": cfg.distance.undefined_pair,
                         "one_hot_metric": cfg.distance.one_hot_metric},
            "snf": {"K": int(K), "mu": cfg.snf.mu, "T": cfg.snf.T,
                    "alpha": cfg.snf.alpha, "tie_policy": cfg.snf.tie_policy,
                    "update": cfg.snf.update},
            "prevalence_gate": {"lo": cfg.trait_cluster.prevalence_gate.lo,
                                "hi": cfg.trait_cluster.prevalence_gate.hi,
                                "min_count": cfg.trait_cluster.prevalence_gate.min_count},
            "mdl_complexity": mdl_complexity,
            "monte_carlo": {"consensus_B": consensus_B, "n_perm": n_perm,
                            "n_boot": n_boot, "ps_n_reps": ps_n_reps,
                            "gap_B": gap_B},
        },
        "selected_k": int(k),
        "k_selection": ksel,
        "layers": layer_names,
        "layer_reports": layer_reports,
        "n_isolates": len(ids),
        "cluster_sizes": {int(c): int((labels == c).sum()) for c in np.unique(labels)},
        "assignment": assignment,
        "descriptive_archetype_profiles": profiles.to_dict("records"),
        "n_defining_features_descriptive": int(profiles["is_defining_descriptive"].sum())
        if len(profiles) else 0,
        "mdl_calibration": mdl_calib,
        "mdl_contribution": mdl_contrib.to_dict("records"),
        "pac": pac_metric(M),
        "cdf_area": cdf_area(M),
        "artifact_diagnostics": artifacts,
        "layer_influence": infl,
        "layer_contribution": layer_contrib,
        "layer_agreement": {
            "pair_agreement": agree["pair_agreement"],
            "fused_vs_layer": agree["fused_vs_layer"],
            "conflict_mean_layer_pairs_only": agree["conflict_mean_layer_pairs_only"],
            "n_layer_pairs": agree["n_layer_pairs"],
        },
        "baselines": base_block,
        "metadata_diagnostics": meta_block,
        "phenotype_validation": pheno_block,
        "fragility": {"F_isolate": frag["F_isolate"].tolist(),
                      "F_feature": frag["F_feature"].tolist(),
                      "boundary_matrix": frag["boundary_matrix"].tolist()},
        "tati_edges": tati["edges"],
        "label_recoverability": recover,
        "post_clustering_inference": inference_block,
    }


# --------------------------------------------------------------------------- #
# run() helpers
# --------------------------------------------------------------------------- #
def _layer_influence(cfg, layer_names, layer_mats, fuse, k, labels, rng, K,
                     *, n_perm=None):
    """Adapt :func:`influence.layer_influence` to per-layer distance metrics.

    Layers can carry different ``kind``s and therefore different distances, so
    the kernel is tagged with the layer role rather than inferred.
    """
    tagged = [_TaggedMatrix(M, role) for M, role in zip(layer_mats, layer_names)]

    def kernel_fn(M_):
        role = getattr(M_, "layer_role", layer_names[0])
        return snf_kernel(_layer_distance(cfg, role, np.asarray(M_)),
                          mu=cfg.snf.mu, K=K)

    return _influence.layer_influence(
        tagged, layer_names, kernel_fn=kernel_fn, fuse_fn=fuse,
        cluster_fn=spectral_from_similarity, k=k, fused_labels=labels,
        rng=rng,
        n_perm=cfg.influence.n_perm if n_perm is None else n_perm,
        collapse_threshold=cfg.influence.collapse_threshold)


class _TaggedMatrix(np.ndarray):
    """ndarray subclass carrying the layer role through slicing/permutation."""

    def __new__(cls, arr, role):
        obj = np.asarray(arr).view(cls)
        obj.layer_role = role
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.layer_role = getattr(obj, "layer_role", None)


def _artifact_diagnostics(layer_names, layer_mats, labels, layer_reports) -> dict:
    """Is the partition an artefact of the empty stratum or of tied affinities?

    ``empty_stratum_ari`` is the ARI between the partition and the binary
    indicator "this isolate carries none of this layer's retained features". A
    high value means the discovered clusters are, at top level, the trait-absent
    stratum produced by the empty-union convention rather than a multi-trait
    structure.
    """
    per_layer = []
    worst = 0.0
    worst_layer = None
    for name, M, rep in zip(layer_names, layer_mats, layer_reports):
        empty = uninformative_rows(M)
        if 0 < empty.sum() < len(empty):
            a = float(adjusted_rand_score(empty.astype(int), labels))
            ami = float(adjusted_mutual_info_score(empty.astype(int), labels))
        else:
            a = float("nan")
            ami = float("nan")
        per_layer.append({
            "layer": name,
            "n_uninformative_rows": int(empty.sum()),
            "frac_uninformative_rows": float(empty.mean()),
            "empty_stratum_ari": a,
            "empty_stratum_ami": ami,
            "effective_dimension": rep["effective_dimension"],
            "n_features_used": rep["n_features_used"],
            "knn_tie_inflation": rep.get("knn_ties", {}).get("tie_inflation"),
        })
        if np.isfinite(a) and a > worst:
            worst, worst_layer = a, name
    per_cluster = {}
    for c in np.unique(labels):
        mask = labels == c
        per_cluster[str(int(c))] = {
            name: float(uninformative_rows(M)[mask].mean())
            for name, M in zip(layer_names, layer_mats)
        }
    # Which combination of layers is an isolate positive in? On the Klebsiella
    # cohort the smaller cluster is 70 % positive in *neither* layer, which
    # "frac_uninformative_per_cluster" reports one layer at a time and therefore
    # does not show. The joint pattern is what tells a reader whether a cluster
    # is a trait archetype or a double empty stratum, and the doubly-positive
    # count is the one a surveillance user is looking for.
    pos = np.array([~uninformative_rows(M) for M in layer_mats])
    pattern_per_cluster = {}
    for c in np.unique(labels):
        mask = labels == c
        counts: Dict[str, int] = {}
        for i in np.flatnonzero(mask):
            key = "+".join([n for j, n in enumerate(layer_names) if pos[j, i]]) or "none"
            counts[key] = counts.get(key, 0) + 1
        pattern_per_cluster[str(int(c))] = dict(
            sorted(counts.items(), key=lambda kv: -kv[1]))
    return {
        "per_layer": per_layer,
        "frac_uninformative_per_cluster": per_cluster,
        "layer_positivity_pattern_per_cluster": pattern_per_cluster,
        "frac_positive_in_no_layer_per_cluster": {
            str(int(c)): float((~pos[:, labels == c].any(axis=0)).mean())
            for c in np.unique(labels)},
        "frac_positive_in_all_layers_per_cluster": {
            str(int(c)): float(pos[:, labels == c].all(axis=0).mean())
            for c in np.unique(labels)},
        "max_empty_stratum_ari": float(worst),
        "empty_stratum_layer": worst_layer,
        "empty_stratum_warning": bool(worst >= 0.5),
        "note": "empty_stratum_ari compares the partition to the indicator "
                "'carries none of this layer's features'; a high value means "
                "the clusters are the trait-absent stratum, an artefact of the "
                "empty-union distance convention rather than trait structure",
    }


def _layer_contribution(cfg, layer_names, layer_mats, fused_labels, k) -> List[dict]:
    rows = []
    for name, M in zip(layer_names, layer_mats):
        D = _layer_distance(cfg, name, M)
        solo = spectral_from_similarity(snf_kernel(D, mu=cfg.snf.mu), k)
        try:
            sil = float(silhouette_score(np.nan_to_num(D, nan=1.0), fused_labels,
                                         metric="precomputed"))
        except Exception:
            sil = float("nan")
        rows.append({"layer": name,
                     "ari_solo_vs_fused": float(adjusted_rand_score(fused_labels, solo)),
                     "silhouette_of_fused_in_layer": sil})
    return rows


_REPORTED_LABELS: Dict[str, object] = {}


def _run_inference(cfg, layer_names, layer_frames, X_df, k, rng, W_fused,
                   K, reported_labels=None) -> dict:
    """Valid post-clustering inference, and an honest statement of its scope.

    Two questions, deliberately separated because a single p-value cannot answer
    both:

    * **Is there structure?** Multi-split feature splitting
      (:func:`inference.feature_split_report`): cluster on a random half of the
      feature *blocks*, test the held-out half, merge the split p-values with a
      rule valid for exchangeable, arbitrarily dependent p-values. Splitting by
      block rather than by feature is required - a per-feature split has
      family-wise error 1.000 under a no-cluster null with correlated features.

    * **Is that structure discrete?** It does not follow. The split test rejects
      on a single continuous gradient as readily as on discrete groups, because
      both make the held-out features depend on the discovery labels.
      :func:`inference.continuum_null_test` answers the second question with a
      parametric bootstrap from a fitted one-factor latent-trait model, so the
      null preserves the feature-feature dependence.

    Count thinning is reported as a comparison where a count model is
    defensible, together with the realised corr(X1, X2) that checks it.
    """
    if k <= 1:
        return {"status": "skipped", "reason": "k = 1: no partition to test",
                "structure_detected": False, "discreteness": None}
    _REPORTED_LABELS["labels"] = reported_labels

    from . import inference as _inf

    # The discovery block is clustered the way the reported partition was:
    # per-layer distances (one_hot layers get the matching coefficient), fused.
    col_layer = []
    for name, df in zip(layer_names, layer_frames):
        col_layer.extend([name] * df.shape[1])
    col_layer = np.asarray(col_layer, dtype=object)

    def cluster_block(block, kk, seed, cols=None):
        block = np.asarray(block)
        roles = col_layer[cols] if cols is not None else None
        if roles is None or len(set(roles.tolist())) < 2:
            role = roles[0] if roles is not None and len(roles) else layer_names[0]
            W = snf_kernel(_layer_distance(cfg, role, block), mu=cfg.snf.mu, K=K)
        else:
            Ws = []
            for role in layer_names:
                sel = np.flatnonzero(roles == role)
                if sel.size:
                    Ws.append(snf_kernel(_layer_distance(cfg, role, block[:, sel]),
                                         mu=cfg.snf.mu, K=K))
            W = snf_fuse(Ws, K=K, T=cfg.snf.T, alpha=cfg.snf.alpha,
                         tie_policy=cfg.snf.tie_policy)
        return spectral_from_similarity(W, kk, random_state=seed)

    split_groups = _feature_split_groups(cfg, layer_names, layer_frames, X_df)
    # A one_hot block is one hypothesis, not one per level: its indicators are
    # mutually exclusive, so once some levels mark a cluster the rest mark its
    # complement by construction. Excluded from the per-feature FDR set and
    # reported separately.
    onehot_cols = [c for name, df in zip(layer_names, layer_frames)
                   if cfg.files[name].kind == "one_hot" for c in df.columns]
    testable = [c for c in X_df.columns if c not in set(onehot_cols)]

    primary = _inf.feature_split_report(
        X_df, k, cluster_fn=cluster_block, rng=rng,
        n_splits=cfg.tva.n_splits, frac_discovery=0.5,
        merge=cfg.tva.merge, q_fdr=cfg.tva.q_fdr, binary=True,
        groups=split_groups, fdr_features=testable,
        reference_labels=_REPORTED_LABELS.get("labels"))

    primary["discreteness"] = _inf.continuum_null_test(
        X_df.values, k, rng=rng, n_boot=int(cfg.tva.n_boot_continuum),
        q_max=int(cfg.tva.continuum_q_max))
    primary["count_thinning_comparison"] = _count_thinning_comparison(
        cfg, layer_names, layer_frames, k, rng)
    primary["excluded_from_fdr"] = {
        "one_hot_columns": len(onehot_cols),
        "reason": "mutually exclusive indicators are one hypothesis, not p",
    }
    return primary


def _feature_split_groups(cfg, layer_names, layer_frames, X_df) -> List[str]:
    """One split-unit label per column of ``X_df``.

    A feature belongs to its declared group if the config names one; otherwise
    to a singleton unit within its layer, except for ``one_hot`` layers whose
    columns are mutually exclusive by construction and therefore form one unit.

    Declared groups only cover the co-inheritance the analyst thought of. On the
    *Klebsiella* cohort the config declares the five virulence loci and nothing
    for AMR, yet the acquired-resistance classes are strongly co-selected
    (aminoglycoside-sulphonamide r = 0.81, sulphonamide-trimethoprim 0.80) and
    the KpVP-1 plasmid loci correlate across declared groups (iro-rmp r = 0.93).
    An all-singleton design over correlated features is the arm the calibration
    study puts at family-wise error 1.000, so the declared groups are unioned
    with the connected components of the thresholded correlation graph
    (:func:`inference.correlation_blocks`). The blocks are a function of the
    feature matrix alone, never of the labels.
    """
    labels: List[str] = []
    for name, df in zip(layer_names, layer_frames):
        spec = cfg.files[name]
        member_of = {}
        for grp, members in (spec.groups or {}).items():
            for m in members:
                member_of[m] = f"{name}:{grp}"
        for col in df.columns:
            if spec.kind == "one_hot":
                labels.append(f"{name}:onehot")
            else:
                labels.append(member_of.get(col, f"{name}:{col}"))
    if len(labels) != X_df.shape[1]:      # defensive: fall back to per-layer
        labels = []
        for name, df in zip(layer_names, layer_frames):
            labels.extend([name] * df.shape[1])

    from . import inference as _inf
    blocks = _inf.correlation_blocks(X_df.values,
                                     threshold=cfg.tva.block_threshold,
                                     groups=labels)
    return [f"block{b}" for b in blocks]


def _count_thinning_comparison(cfg, layer_names, layer_frames, k, rng) -> dict:
    """Negative-binomial count thinning on per-layer aggregates, with its check.

    Reported alongside the feature-split result, never instead of it. Columns
    that cannot plausibly be negative binomial are screened out with a stated
    reason, and the realised ``corr(X1, X2)`` is reported for those that remain:
    a systematically non-zero value means the independence the method rests on
    did not hold (Neufeld et al. 2024, JMLR 25(57), Prop. 11).
    """
    from . import inference as _inf
    from . import tva as _tva

    per_layer = []
    excluded = []
    for name, df in zip(layer_names, layer_frames):
        spec = cfg.files[name]
        if spec.kind == "one_hot":
            excluded.append({"layer": name,
                             "reason": "one_hot layer aggregates to a constant"})
            continue
        groups = spec.groups or _tva.auto_prefix_groups(list(df.columns),
                                                        min_group=2, prefix_len=3)
        if groups:
            g = _tva.binary_layer_to_counts(df, layer_name=name, groups=groups)
            g.columns = [f"{name}_{c}" for c in g.columns]
            per_layer.append(g)
        per_layer.append(_tva.binary_layer_to_counts(df, layer_name=name))
    if not per_layer:
        return {"status": "skipped", "reason": "no count-bearing layers",
                "excluded_layers": excluded}
    counts = pd.concat(per_layer, axis=1)
    screening = _tva.screen_thinnable(counts, min_distinct=cfg.tva.min_distinct)
    if screening["n_kept"] < cfg.tva.min_columns:
        return {"status": "skipped",
                "reason": f"only {screening['n_kept']} of {counts.shape[1]} count "
                          f"columns are NB-thinnable (need >= {cfg.tva.min_columns})",
                "screening": screening, "excluded_layers": excluded}
    kept = counts.loc[:, screening["keep"]]
    Xi = kept.values.astype(np.int64)
    r_arr = _tva.estimate_dispersion_per_feature(
        Xi, method="mle" if cfg.tva.dispersion in ("mle", "pooled_mle") else "mom")
    X1, X2 = _inf.nb_thin(Xi, r=r_arr, eps=cfg.tva.eps, rng=rng)
    dep = _inf.thinning_dependence(Xi, X1, X2)
    block = _tva.tva_report(kept, k=max(k, 2), rng=rng, eps=cfg.tva.eps,
                            n_splits=cfg.tva.n_splits, q_fdr=cfg.tva.q_fdr,
                            dispersion=cfg.tva.dispersion, merge=cfg.tva.merge,
                            min_distinct=cfg.tva.min_distinct,
                            min_columns=cfg.tva.min_columns)
    block["excluded_layers"] = excluded
    block["realised_dependence"] = {
        c: (None if not np.isfinite(v) else float(v))
        for c, v in zip(screening["keep"], dep)}
    finite = dep[np.isfinite(dep)]
    block["max_abs_realised_dependence"] = (float(np.abs(finite).max())
                                            if finite.size else None)
    dependence_diagnostic_passed = bool(
        finite.size and np.abs(finite).max() < 0.05
    )
    block["dependence_diagnostic_passed"] = dependence_diagnostic_passed
    block["guarantee_holds"] = bool(
        block.get("nominal_inference_valid") and dependence_diagnostic_passed
    )
    block["note"] = ("corr(X1, X2) should be ~0; a systematically non-zero value "
                     "means the assumed negative-binomial family or its estimated "
                     "dispersion is wrong; a near-zero value alone does not validate "
                     "parameters estimated from the same matrix "
                     "(Neufeld et al. 2024 JMLR 25(57) Prop. 11)")
    return block


def _run_baselines(cfg, layer_names, layer_mats, X_combined, labels, k, rng, K,
                   externals=None) -> dict:
    """Every candidate partition, scored on the same external criteria.

    Agreement with the fused partition answers "is this the same answer", which
    cannot rank candidates: it measures agreement with the tool's own output.
    What ranks them is agreement with something external - a published score, a
    domain rule, a lineage label - plus a model-fit criterion. Both are reported
    here for the fused partition, the naive concatenation, every single layer, a
    Bernoulli-mixture latent class model, and any domain rules supplied.

    If a single layer or a one-gene rule scores better than the fusion on the
    external criteria, that is the result, and it belongs in the paper.
    """
    D_cat = binary_distance_matrix(np.hstack(layer_mats), metric=cfg.distance.metric,
                                   undefined_pair=cfg.distance.undefined_pair)
    candidates = {
        "fused_snf": np.asarray(labels),
        "concatenation": spectral_from_similarity(
            snf_kernel(D_cat, mu=cfg.snf.mu, K=K), k),
    }
    for name, M in zip(layer_names, layer_mats):
        candidates[f"solo_{name}"] = spectral_from_similarity(
            snf_kernel(_layer_distance(cfg, name, M), mu=cfg.snf.mu, K=K), k)
    bmm_sel = _baselines.bernoulli_mixture_select_k(
        X_combined, [1] + list(cfg.trait_cluster.k_range), rng=rng, n_init=8)
    bmm_k = _baselines.bernoulli_mixture(X_combined, k, rng=rng, n_init=8)
    candidates["bernoulli_mixture"] = np.asarray(bmm_k["labels"])

    table = []
    for name, lab in candidates.items():
        row = {
            "partition": name,
            "n_clusters": int(np.unique(lab).size),
            "cluster_sizes": {int(c): int((lab == c).sum()) for c in np.unique(lab)},
            "ari_vs_fused": float(adjusted_rand_score(labels, lab)),
            "mdl_gain": float(mdl_archetype(X_combined, lab)["gain"]),
        }
        if externals:
            row["external_agreement"] = _baselines.external_agreement(lab, externals)
        table.append(row)
    return {
        "candidates": table,
        "bernoulli_mixture_k_by_bic": bmm_sel["k_selected"],
        "bernoulli_mixture_bic_table": bmm_sel["table"],
        "note": "ari_vs_fused only says whether a candidate agrees with the "
                "tool's own answer. Rank candidates by external_agreement and "
                "mdl_gain; if a single layer or a domain rule wins, report that",
    }


def _phenotype_validation(cfg, ds, strain_ids, labels, layer_names, layer_mats,
                          externals) -> Optional[dict]:
    """Score the partition against measured susceptibility, if any was supplied.

    The competing rules are built generically: one "carries at least one feature
    in this layer" indicator per layer, plus every declared external column
    binarised at its median. The layer indicators are the ones that matter,
    because a fused partition that cannot beat "carries anything in the AMR
    layer" on a laboratory-measured phenotype has not earned the word
    "archetype".
    """
    path = getattr(cfg.dataset, "phenotype", None)
    if not path:
        return None
    from . import phenotype as _ph
    fp = Path(path)
    if not fp.is_absolute():
        fp = cfg.data_root / path
    if not fp.exists():
        return {"status": "skipped", "reason": f"phenotype file not found: {fp}"}
    long_df = pd.read_csv(fp)
    wide = _ph.to_non_susceptible(
        long_df,
        id_column=cfg.dataset.phenotype_id_column,
        antibiotic_column=cfg.dataset.phenotype_antibiotic_column,
        call_column=cfg.dataset.phenotype_call_column,
        intermediate=cfg.dataset.phenotype_intermediate)

    rules = {f"any_feature_in_{nm}": (np.asarray(M).sum(axis=1) > 0).astype(int)
             for nm, M in zip(layer_names, layer_mats)}
    for name, vals in (externals or {}).items():
        v = pd.to_numeric(pd.Series(vals), errors="coerce").to_numpy(dtype=float)
        if np.isfinite(v).sum() < 4 or np.nanmin(v) == np.nanmax(v):
            continue
        rules[f"external_{name}_above_median"] = (
            v > np.nanmedian(v)).astype(int)
    # Stratify by whatever categorical metadata column is most collection-like,
    # so the leave-one-collection-out check has something to hold out.
    strata = None
    md = getattr(ds, "metadata", None)
    col = getattr(cfg.dataset, "phenotype_stratum_column", None) or "Country"
    if md is not None and col in md.columns:
        strata = md.reindex(strain_ids)[col].tolist()
    return _ph.phenotype_concordance(
        labels, list(strain_ids), wide, competing_rules=rules,
        intermediate=cfg.dataset.phenotype_intermediate, strata=strata)


def _metadata_diagnostics(cfg, ds, strain_ids, labels, rng, n_perm,
                          X_df=None, refit=None) -> Optional[dict]:
    md = getattr(ds, "metadata", None)
    if md is None or cfg.dataset.lineage_column is None:
        return None
    md = md.reindex(strain_ids)
    out: Dict[str, object] = {}
    col = cfg.dataset.lineage_column
    if col in md.columns:
        out["lineage_column"] = col
        raw = md[col].tolist()
        out["lineage_concordance"] = _lineage.lineage_concordance(
            labels, raw, n_perm=n_perm, rng=rng)
        # The concordance test answers "is there any lineage signal" and on a
        # real cohort the answer is always yes, more emphatically the larger
        # the cohort. What decides whether the partition may be read as trait
        # structure is how much, so the z is reported beside an out-of-sample
        # variance attribution that does not grow with n, and it is the
        # attribution that carries the threshold.
        att = getattr(cfg, "attribution", None)
        if (att is not None and getattr(att, "enabled", True)
                and X_df is not None and len(X_df.columns)):
            Xa = X_df.to_numpy(dtype=float)
            # The partition the gate judges was fitted on every isolate, so a
            # fold that holds an isolate out of the group means still scores it
            # against clusters it helped draw. Rebuilding the partition inside
            # the fold is what makes r2_partition an out-of-sample quantity,
            # and the bootstrap is set aside there because it duplicates
            # isolates and a refit cannot be done on a duplicated cohort.
            use_refit = refit if getattr(att, "refit_in_fold", True) else None
            out["lineage_attribution"] = _attribution.attribute_partition(
                Xa, labels, raw, folds=att.folds, repeats=att.repeats,
                n_boot=0 if use_refit is not None else att.n_boot,
                n_perm=att.n_perm, refit=use_refit, seed=rng).as_dict()
            out["lineage_attribution"]["lambda_gate"] = att.lambda_gate
            out["clonal_share"] = {
                str(feat): _attribution.clonal_share(
                    X_df[feat].to_numpy(dtype=float), raw, folds=att.folds,
                    repeats=att.repeats, n_boot=att.n_boot,
                    n_perm=att.n_perm, seed=rng).as_dict()
                for feat in X_df.columns
            }
            out["clonal_share_by_layer"] = {
                "all_clustering_features": _attribution.layer_clonal_share(
                    Xa, raw, folds=att.folds, repeats=att.repeats,
                    n_boot=att.n_boot, n_perm=att.n_perm, seed=rng).as_dict()
            }
        # The share says how much. Whether there is anything there at all is a
        # separate question, and surveillance asks it again every year on the
        # same panel. A p-value recomputed at each of an unplanned number of
        # looks controls nothing; an e-value may be read as often as wanted.
        ev = getattr(cfg, "evidence", None)
        if (ev is not None and getattr(ev, "enabled", True)
                and X_df is not None and len(X_df.columns)):
            per_feature = {
                str(feat): _evalues.e_process(
                    X_df[feat].to_numpy(dtype=float), raw, folds=ev.folds,
                    repeats=ev.repeats, seed=rng).as_dict()
                for feat in X_df.columns
            }
            names = list(per_feature)
            decision = _evalues.e_bh(
                [per_feature[k]["e_value"] for k in names], alpha=ev.alpha)
            out["lineage_evidence"] = {
                "per_feature": per_feature,
                "e_bh": dict(decision, rejected_features=[
                    names[i] for i in decision["rejected"]]),
                "note": ("e-value in the betting sense of Vovk and Wang, not "
                         "the BLAST expectation value and not the E-value of "
                         "VanderWeele and Ding"),
            }
        cen = getattr(cfg, "censored", None)
        if cen is not None and getattr(cen, "enabled", True):
            found = _censored_diagnostics(cfg, ds, strain_ids, raw, cen)
            if found:
                out["censored_share"] = found
        # Prevalence per isolate and per lineage are different estimands and a
        # surveillance report quotes only the first. They separate exactly when
        # sampling across lineages is uneven, which is the normal state of a
        # submission-driven collection, so both are emitted per feature.
        surv = getattr(cfg, "surveillance", None)
        n_boot_s = getattr(surv, "n_boot", n_perm) if surv else n_perm
        n_perm_s = getattr(surv, "n_perm", n_perm) if surv else n_perm
        surv_on = getattr(surv, "enabled", True) if surv else True
        if surv_on and X_df is not None and len(X_df.columns):
            out["lineage_resolved_prevalence"] = {
                str(feat): _clonality.lineage_resolved_prevalence(
                    X_df[feat].to_numpy(), raw, n_boot=n_boot_s, rng=rng)
                for feat in X_df.columns
            }
            # How many lineages effectively carry each trait, and whether that
            # number departs from what the cohort's own lineage abundances
            # would produce by chance. Direction is reported because the
            # departure statistic fires on dispersion as readily as on
            # clonality.
            out["trait_concentration"] = {
                str(feat): _clonality.trait_concentration(
                    X_df[feat].to_numpy(), raw, n_perm=n_perm_s, rng=rng)
                for feat in X_df.columns
            }
        cc = cfg.dataset.contrast_column
        if (surv_on and X_df is not None and len(X_df.columns) and cc
                and cc in md.columns
                and len(cfg.dataset.contrast_levels) == 2):
            # A prevalence difference between two collections splits into a
            # change in the lineage mix and a change in the within-lineage
            # rate. The two imply opposite interventions, and the reported
            # prevalence cannot separate them.
            lv_a, lv_b = (str(v) for v in cfg.dataset.contrast_levels)
            side = md[cc].astype(str)
            in_a = (side == lv_a).to_numpy()
            in_b = (side == lv_b).to_numpy()
            lin_arr = np.asarray(raw, dtype=object)
            panel = _clonality.decompose_panel(
                X_df.loc[in_a], lin_arr[in_a],
                X_df.loc[in_b], lin_arr[in_b],
                n_boot=n_boot_s, rng=rng,
                q=getattr(surv, "q_fdr", 0.05) if surv else 0.05,
                min_shared_support=(getattr(surv, "min_shared_support", 0.8)
                                    if surv else 0.8))
            out["prevalence_decomposition"] = {
                "contrast_column": cc,
                "levels": [lv_a, lv_b],
                "n_a": int(in_a.sum()),
                "n_b": int(in_b.sum()),
                "family": panel["family"],
                "per_feature": panel["per_agent"],
            }
        # MLST variant notation (ST17-1LV) makes a variant of a common clone a
        # different lineage string, which shrinks the number of informative
        # pairs and biases the diagnostic towards "not confounded". Report the
        # collapsed version too whenever collapsing changes anything.
        collapsed = _lineage.collapse_st_variants(raw)
        if not np.array_equal(collapsed, np.asarray(raw, dtype=object)):
            out["lineage_concordance_variants_collapsed"] = \
                _lineage.lineage_concordance(labels, collapsed, n_perm=n_perm,
                                             rng=rng, collapse_variants=True)
    ext_cols = [c for c in cfg.dataset.external_columns if c in md.columns]
    if ext_cols:
        out["external_agreement"] = _baselines.external_agreement(
            labels, {c: md[c].tolist() for c in ext_cols})
    # Every metadata column, not just the lineage and external ones. Country,
    # year and source were read in and never emitted, so the three questions a
    # surveillance user always asks - is this cluster growing, where is it, is
    # it in animals - were unanswerable from the output. `top` bounds the size.
    comp_cols = [c for c in ([col] + list(cfg.dataset.external_columns)
                             + list(md.columns))
                 if c and c in md.columns]
    comp_cols = list(dict.fromkeys(comp_cols))
    if comp_cols:
        out["cluster_composition"] = _lineage.cluster_composition(labels, md, comp_cols)
    return out or None


# =============================================================================
# Synthetic validation harness
# =============================================================================
def validate(cfg=None, *, seed: int = 42, overlap: float = 0.05,
             n: int = 150, k_true: int = 3, **options) -> dict:
    """Planted-truth recovery check.

    Generates two binary layers with ``k_true`` planted archetypes at separation
    ``1 - overlap``, runs the distance -> fusion -> spectral chain, and reports
    the adjusted Rand index against the planted partition together with the k
    recovered by the MDL sweep. ``overlap`` is exposed so the harness can be run
    as a difficulty sweep rather than at one easy setting; see
    ``benchmarks/recovery_sweep.py``.
    """
    ari_min = cfg.validation.cluster_ari_min if cfg is not None else 0.80
    k_range = list(cfg.trait_cluster.k_range) if cfg is not None else [2, 3, 4, 5]
    mu = cfg.snf.mu if cfg is not None else 0.5
    T = cfg.snf.T if cfg is not None else 20

    amr, vir, truth = synth_cluster_archetypes(
        n=n, p_amr=21, p_vir=60, k_true=k_true, overlap=overlap, seed=seed)
    W = snf_fuse([snf_kernel(binary_distance_matrix(amr.values), mu=mu),
                  snf_kernel(binary_distance_matrix(vir.values), mu=mu)], T=T)
    pred = spectral_from_similarity(W, k_true)
    X = np.hstack([amr.values, vir.values]).astype(int)
    sweep = [(int(k), float(mdl_archetype(X, spectral_from_similarity(W, k))["gain"]))
             for k in k_range]
    return {
        "ari": float(adjusted_rand_score(truth, pred)),
        "nmi": float(normalized_mutual_info_score(truth, pred)),
        "selected_k": int(max(sweep, key=lambda x: x[1])[0]),
        "k_true": int(k_true), "n": int(n), "overlap": float(overlap),
        "threshold": float(ari_min),
        "passed": bool(adjusted_rand_score(truth, pred) >= ari_min),
        "mdl_gain_sweep": sweep,
    }

def _censored_diagnostics(cfg, ds, strain_ids, lineage, cen) -> Optional[dict]:
    """The clonal share of each antimicrobial read from the dilution panel.

    A dichotomised call and a recorded MIC are the same likelihood at two
    interval widths, so this reports the same quantity as the binary share on
    a scale that has not thrown away where in the panel each isolate fell. The
    panel geometry is reported first and decides which modes the data support:
    a panel with fewer than three tested wells cannot carry an interval, and a
    panel with a fifth of its readings piled on an end well cannot be read as
    a set of point values.
    """
    mic = getattr(ds, "mic", None)
    if mic is None or not len(mic):
        return None
    idc = cfg.dataset.mic_id_column
    abc = cfg.dataset.mic_antibiotic_column
    vc = cfg.dataset.mic_value_column
    opc = cfg.dataset.mic_operator_column
    ids = pd.Index(strain_ids).astype(str)
    lineage = pd.Series(list(lineage), index=ids)

    out: Dict[str, object] = {"join": getattr(ds, "mic_join", None),
                              "per_agent": {}}
    for agent, block in mic.groupby(mic[abc].astype(str), sort=True):
        block = block.drop_duplicates(subset=[idc], keep="first")
        block = block.set_index(block[idc].astype(str))
        values = pd.to_numeric(block[vc], errors="coerce").reindex(ids)
        present = values.notna().to_numpy()
        if present.sum() < 8:
            continue
        v = values.to_numpy(dtype=float)[present]
        lin = lineage.to_numpy(dtype=object)[present]
        ops = (block[opc].reindex(ids).to_numpy(dtype=object)[present]
               if opc else None)
        geometry = _censored.panel_geometry(v)
        lo, hi = _censored.intervals_from_mic(
            v, operators=ops,
            treat_end_wells_as_censored=cen.end_wells_censored)
        share = _censored.censored_clonal_share(lo, hi, lin,
                                                n_boot=cen.n_boot, seed=0)
        entry: Dict[str, object] = {
            "n": int(present.sum()),
            "panel": geometry.as_dict(),
            "share": share.as_dict(),
        }
        if cen.sensitivity and "point" not in geometry.admissible_modes:
            entry["end_well_sensitivity"] = _censored.sensitivity_endpoints(
                v, lin, operators=ops, n_boot=max(cen.n_boot // 4, 25), seed=0)
        out["per_agent"][str(agent)] = entry
    return out if out["per_agent"] else None
