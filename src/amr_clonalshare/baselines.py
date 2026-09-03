"""baselines.py — the comparators every multi-layer clustering claim needs.

A multi-view method earns its complexity only if it beats the obvious
alternatives on the same data. The multi-omic clustering literature treats that
comparison as mandatory (Rappoport & Shamir 2018, Nucleic Acids Res
46:10546-10562, doi:10.1093/nar/gky889; Tini et al. 2019, Brief Bioinform
20:1269-1279, doi:10.1093/bib/bbx167), and on the *Klebsiella* example shipped
with this package the alternatives are not obviously worse - which is a result,
not an embarrassment, and is reported as one.

Baselines provided
------------------
``concatenate_cluster``
    Stack every layer's columns and run the same distance + spectral pipeline.
    The thing similarity-network fusion is supposed to improve on.
``single_layer_cluster``
    Each layer on its own. If one of these matches the fused partition, the
    fusion did not fuse.
``bernoulli_mixture``
    A finite mixture of independent Bernoullis, i.e. latent class analysis
    (Lazarsfeld & Henry 1968; Vermunt & Magidson 2002), fitted by EM with
    random restarts. This is the textbook generative model for the object this
    package calls an "archetype": it gives k by BIC, per-isolate posterior
    responsibilities instead of a consensus heuristic, per-class Bernoulli
    parameters instead of a Delta-p table, and it handles all-zero rows through
    the likelihood rather than through a Jaccard convention. It costs
    milliseconds. Any pipeline that does not beat it should say so.
``external_agreement``
    ARI/AMI of a partition against an external labelling (a published score, a
    rule, a lineage) - the honest way to ask "did we learn anything that was
    not already in one column of the input".
"""
from __future__ import annotations

from typing import Callable, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

__all__ = [
    "bernoulli_mixture",
    "bernoulli_mixture_select_k",
    "concatenate_cluster",
    "single_layer_cluster",
    "external_agreement",
]


def _bmm_em(X: np.ndarray, k: int, rng: np.random.Generator, *,
            max_iter: int = 300, tol: float = 1e-6,
            eps: float = 1e-6):
    n, p = X.shape
    pi = np.full(k, 1.0 / k)
    theta = rng.uniform(0.15, 0.85, size=(k, p))
    ll_old = -np.inf
    resp = np.full((n, k), 1.0 / k)
    for _ in range(max_iter):
        lt = np.log(np.clip(theta, eps, 1 - eps))
        lnt = np.log(np.clip(1 - theta, eps, 1 - eps))
        log_lik = X @ lt.T + (1 - X) @ lnt.T + np.log(np.clip(pi, eps, None))
        mx = log_lik.max(axis=1, keepdims=True)
        w = np.exp(log_lik - mx)
        s = w.sum(axis=1, keepdims=True)
        resp = w / s
        ll = float((np.log(s) + mx).sum())
        nk = resp.sum(axis=0) + eps
        pi = nk / n
        theta = (resp.T @ X + eps) / nk[:, None]
        if abs(ll - ll_old) < tol * max(1.0, abs(ll_old)):
            ll_old = ll
            break
        ll_old = ll
    return ll_old, pi, theta, resp


def bernoulli_mixture(X: np.ndarray, k: int, *,
                      rng: Optional[np.random.Generator] = None,
                      n_init: int = 10, max_iter: int = 300) -> Dict[str, object]:
    """Fit a mixture of ``k`` independent Bernoullis (latent class model).

    Returns the best of ``n_init`` EM restarts: log-likelihood, BIC, mixing
    proportions, per-class Bernoulli parameters, posterior responsibilities and
    hard labels.
    """
    if rng is None:
        rng = np.random.default_rng()
    Xf = np.asarray(X, dtype=float)
    n, p = Xf.shape
    best = None
    for _ in range(int(n_init)):
        out = _bmm_em(Xf, k, rng, max_iter=max_iter)
        if best is None or out[0] > best[0]:
            best = out
    ll, pi, theta, resp = best
    n_par = k * p + (k - 1)
    return {
        "k": int(k),
        "loglik": float(ll),
        "n_par": int(n_par),
        "bic": float(-2 * ll + n_par * np.log(n)),
        "aic": float(-2 * ll + 2 * n_par),
        "weights": pi.tolist(),
        "theta": theta,
        "responsibilities": resp,
        "labels": resp.argmax(axis=1).astype(int),
        "max_posterior": resp.max(axis=1),
    }


def bernoulli_mixture_select_k(X: np.ndarray, k_range: Sequence[int], *,
                               rng: Optional[np.random.Generator] = None,
                               n_init: int = 10) -> Dict[str, object]:
    """Fit the Bernoulli mixture over ``k_range`` and select k by BIC."""
    if rng is None:
        rng = np.random.default_rng()
    fits = [bernoulli_mixture(X, int(k), rng=rng, n_init=n_init) for k in k_range]
    best = min(fits, key=lambda f: f["bic"])
    return {
        "k_selected": int(best["k"]),
        "criterion": "BIC",
        "table": [{"k": f["k"], "loglik": f["loglik"], "n_par": f["n_par"],
                   "bic": f["bic"], "aic": f["aic"]} for f in fits],
        "best": best,
    }


def concatenate_cluster(layer_mats: List[np.ndarray], k: int, *,
                        distance_fn: Callable[[np.ndarray], np.ndarray],
                        kernel_fn: Callable[[np.ndarray], np.ndarray],
                        cluster_fn: Callable[[np.ndarray, int], np.ndarray]
                        ) -> np.ndarray:
    """Cluster the column-wise concatenation of all layers."""
    X = np.hstack([np.asarray(M) for M in layer_mats])
    return cluster_fn(kernel_fn(distance_fn(X)), k)


def single_layer_cluster(layer_mats: List[np.ndarray], k: int, *,
                         distance_fn: Callable[[np.ndarray], np.ndarray],
                         kernel_fn: Callable[[np.ndarray], np.ndarray],
                         cluster_fn: Callable[[np.ndarray, int], np.ndarray]
                         ) -> List[np.ndarray]:
    """Cluster each layer independently."""
    return [cluster_fn(kernel_fn(distance_fn(np.asarray(M))), k)
            for M in layer_mats]


def external_agreement(labels: Sequence[int],
                       externals: Dict[str, Sequence]) -> Dict[str, Dict[str, float]]:
    """ARI and AMI of a partition against each external labelling.

    Use this to answer the only question that matters for a new clustering
    method on a well-studied organism: does it recover anything that a single
    published score, or a two-gene rule, does not already give for free?
    """
    lab = np.asarray(list(labels))
    out: Dict[str, Dict[str, float]] = {}
    for name, ext in externals.items():
        ser = pd.Series(list(ext), dtype="object")
        ok = np.array(ser.notna().to_numpy(), dtype=bool, copy=True)
        e = np.asarray(ser.fillna("__missing__").astype(str).to_numpy(), dtype=object)
        ok &= ~np.isin(e, ("nan", "None", "", "__missing__"))
        if ok.sum() < 2 or np.unique(e[ok]).size < 2:
            continue
        out[name] = {
            "ari": float(adjusted_rand_score(e[ok], lab[ok])),
            "ami": float(adjusted_mutual_info_score(e[ok], lab[ok])),
            "n": int(ok.sum()),
        }
    return out
