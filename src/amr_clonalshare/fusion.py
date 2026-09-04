"""fusion.py — similarity-network fusion, faithful and order-invariant.

Implements the affinity kernel and cross-diffusion of Wang et al. (2014),
"Similarity network fusion for aggregating data types on a genomic scale",
Nature Methods 11:333-337, doi:10.1038/nmeth.2810, with two corrections and
one behaviour that the reference implementations leave undefined.

1. **Local scaling is local.** Wang et al. Eq. (1) sets

       W(i, j) = exp( -rho(x_i, x_j)^2 / (mu * eps_ij) ),
       eps_ij  = ( mean_{k in N_i} rho(x_i, x_k)
                 + mean_{k in N_j} rho(x_j, x_k)
                 + rho(x_i, x_j) ) / 3

   where ``N_i`` is the set of ``x_i``'s ``K`` nearest neighbours. The previous
   implementation used the mean distance from ``i`` to **all** ``n`` points,
   which removes the local scaling the kernel exists to provide and widens the
   bandwidth by a large, data-dependent factor. Both reference implementations
   (SNFtool ``affinityMatrix.R``; snfpy ``snf.compute.make_affinity``) use the
   K-nearest-neighbour mean, and this module now does too.

2. **In-loop symmetrisation.** SNFtool re-symmetrises each view's transition
   matrix inside every diffusion iteration; the previous implementation
   symmetrised only the final average. Asymmetry accumulates over ``T``
   iterations, so this is not cosmetic.

3. **Ties are resolved by value, not by row order.** The K-nearest-neighbour
   sparsification ``argsort(-W[i])[:K]`` is ill-defined when more than ``K``
   candidates share the K-th largest affinity, and binary trait data produces
   exactly that: on the *Klebsiella* example the median row of the virulence
   layer has ~1000 candidates tied at the cut. Which of them ``argsort`` keeps
   then depends on the input row order, so the clustering was not equivariant
   under relabelling of samples: permuting the rows of the input CSV and
   un-permuting the labels changed the partition (ARI 0.59-0.98 against the
   unpermuted run). A clustering must be a function of the data, not of the
   file's row order.

   ``sparsify`` therefore keeps **every** neighbour whose affinity ties the
   K-th largest ("inclusive" policy). The result is order-invariant by
   construction. Where ties are pervasive the retained degree exceeds ``K``,
   which is real information about the data rather than something to hide:
   :func:`knn_tie_diagnostics` reports it, and a mean degree far above ``K``
   means the affinity has too few distinct values for a K-nearest-neighbour
   graph to be meaningful.

4. **The diffusion update is the reference update. There are two of them, and
   they are not the same.** An earlier draft of this docstring recorded, as a
   fourth deviation, that ``snf_fuse`` re-applies the full row normalisation at
   every iteration "whereas SNFtool and snfpy add a constant to the diagonal and
   normalise once", and reported that under the reference update the
   *Klebsiella* fusion collapsed onto the virulence layer at ARI 1.0000. Both
   halves of that were wrong. What is true, read off the published sources:

   =====================  ===========================================  ==========
   implementation         in-loop update                               era
   =====================  ===========================================  ==========
   Wang's MATLAB          ``Wall{j} = B0_normalized(Wall{j}, alpha)``   2014
   SNFtool <= 2.2         ``Wall[[j]] = nextW[[j]] + diag(n)``          2014
   snfpy 0.2.2            ``_B0_normalized(aff0, alpha) = W + alpha I`` port of 2014
   SNFtool >= 2.2.1       ``Wall[[j]] <- normalize(nextW[[j]])``        2017-11-24
   SNFtool 2.3.1 (CRAN)   same                                         2021-06-11
   =====================  ===========================================  ==========

   SNFtool 2.2.1 replaced ``+ diag(n)`` with a per-iteration renormalisation and
   at the same time replaced ``normalize(X) = X / rowSums(X)`` with
   ``X / (2 * (rowSums(X) - diag(X)))`` and ``diag(X) <- 0.5`` - which is
   exactly ``_renormalise`` here at ``alpha = 0.5``, and exactly equation (2) of
   the paper. So this package implements the *maintained* reference update, and
   what the previous docstring called "the reference" is the 2014 update that
   SNFtool's own authors replaced nine years ago and that snfpy still ports.

   Both are now available (``update="renormalise"`` / ``"add_identity"``) and
   both are benchmarked in ``benchmarks/snf_update_benchmark.py``, which
   reproduces ``snfpy.snf`` from ``_cross_diffuse(update="add_identity")`` to a
   maximum relative difference of 2e-16. The empirical answer is that the choice
   does not matter here: over 20 replicates at each of four overlap levels the
   planted-truth ARI is identical for the two updates at overlaps 0.05/0.15/0.25
   (paired margin exactly 0.0000, 100 % ties) and indistinguishable at 0.35
   (-0.0010 +/- 0.0136). On the *Klebsiella* cohort the two give the same
   cluster sizes (885/615), the same regime (``complementary``, no collapse) and
   the same effective number of layers (1.847 vs 1.853) - and, contrary to the
   retracted claim, neither collapses onto the virulence layer (ARI -0.0024 and
   -0.0010 against the virulence-alone partition).

   Two smaller infidelities found while settling this are also fixed: ``P0`` is
   now symmetrised after normalisation, and the in-loop symmetrisation now
   happens *after* the update rather than before, both as the references do.
   Row-normalising a symmetric matrix with unequal row sums destroys symmetry,
   so the order is not cosmetic over 20 iterations.

A single-layer call is legal and returns that layer's normalised transition
matrix (the previous implementation raised ``ZeroDivisionError``).
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional

import numpy as np

__all__ = [
    "snf_kernel",
    "snf_fuse",
    "sparsify",
    "knn_tie_diagnostics",
    "default_K",
]


def default_K(n: int) -> int:
    """Wang et al. neighbourhood size heuristic used by this package."""
    return int(max(2, min(math.ceil(n / 5), 30)))


def snf_kernel(D: np.ndarray, *, mu: float = 0.5,
               K: Optional[int] = None) -> np.ndarray:
    """Self-tuning affinity from a distance matrix (Wang et al. 2014, Eq. 1).

    ``eps_ij`` uses the mean distance from ``i`` to its ``K`` nearest
    neighbours (excluding itself), likewise for ``j``, plus ``d_ij``.
    """
    D = np.asarray(D, dtype=float)
    n = D.shape[0]
    if n == 0:
        return np.zeros((0, 0))
    if K is None:
        K = default_K(n)
    K = int(np.clip(K, 1, max(n - 1, 1)))

    # Mean distance to the K nearest neighbours, self-distance excluded.
    Dm = D.copy()
    np.fill_diagonal(Dm, np.inf)
    if n > 1:
        part = np.partition(Dm, K - 1, axis=1)[:, :K]
        knn_mean = np.where(np.isfinite(part), part, np.nan)
        with np.errstate(invalid="ignore"):
            mean_d = np.nanmean(knn_mean, axis=1)
        mean_d = np.nan_to_num(mean_d, nan=0.0)
    else:
        mean_d = np.zeros(1)

    eps = (mean_d[:, None] + mean_d[None, :] + D) / 3.0
    eps = np.where(eps < 1e-12, 1e-12, eps)
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        W = np.exp(-(D ** 2) / (mu * eps))
    W = np.nan_to_num(W, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(W, 0.0)
    return (W + W.T) / 2.0


def sparsify(W: np.ndarray, K: int, *, tie_policy: str = "inclusive"
             ) -> np.ndarray:
    """Row-normalised K-nearest-neighbour sparsification of an affinity matrix.

    ``tie_policy="inclusive"`` (default) keeps every neighbour whose affinity
    is greater than or equal to the K-th largest in that row, making the result
    invariant to row order. ``tie_policy="strict"`` reproduces the old
    order-dependent ``argsort`` behaviour and exists only so the regression
    test can demonstrate the difference.
    """
    W = np.asarray(W, dtype=float)
    n = W.shape[0]
    if n == 0:
        return W.copy()
    K = int(np.clip(K, 1, max(n - 1, 1)))
    S = np.zeros_like(W)
    if tie_policy == "strict":
        for i in range(n):
            topk = np.argsort(-W[i])[:K]
            S[i, topk] = W[i, topk]
    elif tie_policy == "inclusive":
        # K-th largest value per row (diagonal is 0 and excluded by ordering)
        thresh = -np.partition(-W, K - 1, axis=1)[:, K - 1]
        S = np.where(W >= thresh[:, None], W, 0.0)
        np.fill_diagonal(S, 0.0)
    else:
        raise ValueError(f"unknown tie_policy {tie_policy!r}")
    rs = S.sum(axis=1, keepdims=True)
    return np.divide(S, np.where(rs > 0, rs, 1.0))


def knn_tie_diagnostics(W: np.ndarray, K: int) -> Dict[str, float]:
    """How badly is the K-nearest-neighbour graph degenerate through ties?

    ``mean_degree`` is the average number of neighbours retained by the
    inclusive policy. ``tie_inflation = mean_degree / K``; a value near 1 means
    the graph is well defined, and a large value means the affinity takes too
    few distinct values for a KNN graph to carry information.
    """
    W = np.asarray(W, dtype=float)
    n = W.shape[0]
    if n < 2:
        return {"K": float(K), "mean_degree": 0.0, "tie_inflation": 1.0,
                "frac_rows_tied": 0.0, "n_distinct_affinities": 0.0}
    K = int(np.clip(K, 1, n - 1))
    thresh = -np.partition(-W, K - 1, axis=1)[:, K - 1]
    deg = (W >= thresh[:, None]).sum(axis=1) - (np.diag(W) >= thresh).astype(int)
    return {
        "K": float(K),
        "mean_degree": float(deg.mean()),
        "tie_inflation": float(deg.mean() / K),
        "frac_rows_tied": float((deg > K).mean()),
        "n_distinct_affinities": float(np.unique(np.round(W, 12)).size),
    }


def _renormalise(W: np.ndarray, alpha: float, idx: np.ndarray) -> np.ndarray:
    """SNFtool >= 2.2.1 ``normalize()``: off-diagonal mass ``1 - alpha``, diagonal ``alpha``."""
    P = np.array(W, dtype=float, copy=True)
    P[idx, idx] = 0.0
    row = P.sum(axis=1, keepdims=True)
    np.divide(P, np.where(row > 0, row, 1.0), out=P)
    P *= (1.0 - alpha)
    P[idx, idx] += alpha
    return P


def _rowstochastic(W: np.ndarray) -> np.ndarray:
    """SNFtool <= 2.2 / snfpy ``normalize()``: ``X / rowSums(X)``, diagonal kept."""
    P = np.array(W, dtype=float, copy=True)
    row = P.sum(axis=1, keepdims=True)
    np.divide(P, np.where(row != 0, row, 1.0), out=P)
    return P


def _cross_diffuse(P: List[np.ndarray], S: List[np.ndarray], *, T: int,
                   alpha: float, update: str, idx: np.ndarray) -> List[np.ndarray]:
    """``T`` iterations of the cross-diffusion loop, shared by both update rules.

    ``P`` are the normalised per-view transition matrices, ``S`` their sparsified
    counterparts. The recursion itself,
    ``P_v <- S_v (sum_{u != v} P_u / (m - 1)) S_v^T``, is common to every
    published implementation; ``update`` selects what is done to the result.
    Both references then symmetrise, and both do so *after* the update rather
    than before (SNFtool ``(Wall[[j]] + t(Wall[[j]]))/2`` after ``normalize``;
    snfpy ``check_symmetric`` inside ``_B0_normalized``). The order matters:
    row-normalising a symmetric matrix with unequal row sums destroys symmetry,
    so normalising second would leave an asymmetric transition matrix to
    accumulate over ``T`` iterations.
    """
    m = len(P)
    # `others` is the mean of the other views' transition matrices. Computing it
    # as (total - P[v]) / (m - 1) from one running total is O(m) per iteration
    # instead of O(m^2), which matters because this loop dominates the runtime.
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        for _ in range(T):
            total = P[0].copy()
            for u in range(1, m):
                total += P[u]
            P_new = []
            for v in range(m):
                others = (total - P[v]) / (m - 1)
                Pv = S[v] @ others @ S[v].T
                if update == "renormalise":
                    Pv = _renormalise(Pv, alpha, idx)
                else:
                    Pv[idx, idx] += alpha
                Pv += Pv.T                       # in-loop symmetrisation
                Pv *= 0.5
                P_new.append(Pv)
            P = P_new
    return P


def snf_fuse(W_layers: List[np.ndarray], *, K: Optional[int] = None,
             T: int = 20, alpha: float = 0.5,
             tie_policy: str = "inclusive",
             update: str = "renormalise") -> np.ndarray:
    """Cross-diffusion fusion of per-layer affinity matrices (Wang et al. 2014).

    Parameters
    ----------
    W_layers : list of (n, n) affinity matrices. A list of length 1 is legal
        and returns that layer's normalised transition matrix.
    K : neighbourhood size for the sparsification; defaults to
        ``min(ceil(n/5), 30)``.
    T : diffusion iterations.
    alpha : self-weight. Under ``update="renormalise"`` it is the diagonal of
        the row-normalised transition matrix (0.5 in Wang et al. and in
        SNFtool); under ``update="add_identity"`` it is the constant added to
        the diagonal at every iteration (1.0 in snfpy and in SNFtool <= 2.2).
    tie_policy : passed to :func:`sparsify`.
    update : ``"renormalise"`` (default) re-applies the full row normalisation
        after every diffusion step, as SNFtool >= 2.2.1 does.
        ``"add_identity"`` adds ``alpha * I`` and lets the row mass float, as
        Wang's MATLAB code, SNFtool <= 2.2 and snfpy do. See the module
        docstring for which reference uses which and what the difference costs.
    """
    if not W_layers:
        raise ValueError("snf_fuse needs at least one affinity matrix")
    if update not in ("renormalise", "add_identity"):
        raise ValueError(f"unknown update {update!r}")
    n = W_layers[0].shape[0]
    for W in W_layers:
        if W.shape != (n, n):
            raise ValueError("all affinity matrices must be (n, n) and equal-sized")
    if K is None:
        K = default_K(n)
    m = len(W_layers)
    idx = np.arange(n)

    if update == "renormalise":
        P = [_renormalise(W, alpha, idx) for W in W_layers]
    else:
        P = [_rowstochastic(W) for W in W_layers]
    P = [(Pv + Pv.T) / 2.0 for Pv in P]          # both references symmetrise P0
    if m == 1:
        return P[0]

    S = [sparsify(W, K, tie_policy=tie_policy) for W in W_layers]
    P = _cross_diffuse(P, S, T=T, alpha=alpha, update=update, idx=idx)

    P_fused = P[0].copy()
    for v in range(1, m):
        P_fused += P[v]
    P_fused /= m
    if update == "renormalise":
        P_fused = _renormalise(P_fused, alpha, idx)
    else:
        P_fused = _rowstochastic(P_fused)
        P_fused[idx, idx] += 1.0
    return (P_fused + P_fused.T) / 2.0
