"""influence.py — layer-influence diagnostics for similarity-network fusion.

Motivation
----------
Similarity Network Fusion (Wang et al. 2014, Nat Methods 11:333-337,
doi:10.1038/nmeth.2810) is routinely presented as a way to let *every* input
layer contribute to a joint partition. Nothing in the algorithm guarantees
this. Cross-diffusion is a fixed-point iteration in which a layer whose
affinity matrix is sharper (larger within/between contrast) can drive the
fused matrix to its own structure while the remaining layers act as near-noise.
The user is never told when this happens: the fused partition looks like a
consensus but is a single layer wearing a consensus costume.

This module makes that failure mode observable. It computes, for a fusion of
``m`` layers:

``ari_solo[l]``
    ARI between the partition obtained from layer ``l`` alone and the fused
    partition. High for a dominating layer. This is descriptive only: a layer
    can agree with the fused partition without having *caused* it.

``delta_loo[l]`` (leave-one-layer-out influence)
    ``1 - ARI(fuse(all layers), fuse(all layers except l))``. This is causal in
    the interventional sense: it is the change in the output produced by
    removing the layer. ``delta_loo == 0`` means the layer is inert - removing
    it changes nothing.

``weight[l]``
    ``delta_loo[l]`` normalised to sum to 1 (uniform if all deltas are 0).

``n_eff``
    The Hill number of order 1 of ``weight``: ``exp(-sum w log w)``, i.e. the
    *effective number of contributing layers*. ``n_eff -> 1`` means the fusion
    has collapsed onto one layer; ``n_eff -> m`` means all layers contribute
    equally. The same participation-ratio idea is standard for quantifying
    effective diversity (Hill 1973, Ecology 54:427-432, doi:10.2307/1934352)
    and effective dimensionality; the inverse Simpson variant
    ``1 / sum w^2`` is reported alongside as ``n_eff_simpson``.

``p_value[l]`` (optional, ``n_perm > 0``)
    A calibrated test of "layer ``l`` contributes nothing beyond its marginal
    structure". The null permutes the *rows* of layer ``l`` only, which
    destroys that layer's alignment with the other layers while preserving its
    own feature marginals and its own internal correlation structure exactly.
    Re-fusing under the null gives a null distribution of ``delta_loo[l]``;
    the p-value is the Phipson-Smyth corrected upper tail. Row permutation is
    the standard null for "is this data view informative about the joint
    structure" in multi-view integration and is the same device used for the
    MDL calibration elsewhere in this package.

``regime``
    ``n_eff`` on its own is not a measure of contribution: it approaches ``m``
    both when every layer says the same thing (so removing one changes nothing)
    and when no layer says anything (so removing one flips an unstable
    partition). The solo agreements resolve the ambiguity, and ``regime``
    reports the resolution as one of ``complementary``, ``collapsed``
    (one layer drives it), ``redundant`` (all agree; the fusion is unnecessary)
    or ``unstructured`` (large influences, no layer reproduces the partition -
    what a fusion of noise looks like). ``collapse`` is ``True``, and the CLI
    exits non-zero, for the last two of those and for ``uninformative``.

Nothing here is organism- or panel-specific; the module takes affinity
matrices and a fusion callable.
"""
from __future__ import annotations

import math
from typing import Any, Callable, Dict, List, Optional, Sequence

import numpy as np
from sklearn.metrics import adjusted_rand_score

__all__ = [
    "layer_influence",
    "effective_n_layers",
    "hill_number",
]


def hill_number(weights: "Sequence[float] | np.ndarray", order: int = 1) -> float:
    """Hill number of a weight vector (Hill 1973). order=1 -> exp(Shannon)."""
    w = np.asarray(weights, dtype=float)
    w = w[w > 0]
    if w.size == 0:
        return 0.0
    w = w / w.sum()
    if order == 1:
        return float(math.exp(-np.sum(w * np.log(w))))
    if order == 2:
        return float(1.0 / np.sum(w ** 2))
    return float((np.sum(w ** order)) ** (1.0 / (1.0 - order)))


def effective_n_layers(deltas: "Sequence[float] | np.ndarray") -> Dict[str, Any]:
    """Effective number of contributing layers from leave-one-out influences."""
    d = np.asarray(deltas, dtype=float)
    d = np.clip(d, 0.0, None)
    if d.sum() <= 0:
        w = np.full(len(d), 1.0 / max(len(d), 1))
    else:
        w = d / d.sum()
    return {
        "weights": w.tolist(),
        "n_eff": hill_number(w, order=1),
        "n_eff_simpson": hill_number(w, order=2),
    }


def layer_influence(
    layer_mats: List[np.ndarray],
    layer_names: List[str],
    *,
    kernel_fn: Callable[[np.ndarray], np.ndarray],
    fuse_fn: Callable[[List[np.ndarray]], np.ndarray],
    cluster_fn: Callable[[np.ndarray, int], np.ndarray],
    k: int,
    fused_labels: Optional[np.ndarray] = None,
    rng: Optional[np.random.Generator] = None,
    n_perm: int = 0,
    collapse_threshold: float = 1.5,
    unstructured_solo_max: float = 0.35,
) -> Dict[str, Any]:
    """Quantify how much each layer actually drives the fused partition.

    Parameters
    ----------
    layer_mats : list of (n, p_l) binary/numeric matrices, one per layer.
    layer_names : names, same order.
    kernel_fn : matrix -> affinity matrix for one layer.
    fuse_fn : list of affinity matrices -> fused affinity matrix. Must accept a
        list of length 1 (see ``core.snf_fuse``, which returns the single
        normalised matrix in that case).
    cluster_fn : (affinity, k) -> labels.
    k : number of clusters.
    fused_labels : labels from the full fusion; recomputed if None.
    rng : generator for the permutation null.
    n_perm : permutations per layer for the contribution test (0 = skip).
    collapse_threshold : n_eff below this flags a collapsed fusion.
    unstructured_solo_max : maximum solo-vs-fused ARI still compatible with the
        ``unstructured`` verdict (see the inline note on its calibration).

    Returns
    -------
    dict with ``per_layer`` records, ``n_eff``, ``n_eff_simpson``, ``collapse``,
    ``n_layers`` and ``n_perm``.
    """
    m = len(layer_mats)
    if m == 0:
        raise ValueError("layer_influence needs at least one layer")
    if len(layer_names) != m:
        raise ValueError("layer_names must match layer_mats")
    if rng is None:
        rng = np.random.default_rng()

    W = [kernel_fn(X) for X in layer_mats]
    if fused_labels is None:
        fused_labels = cluster_fn(fuse_fn(W), k)

    solo_labels = [cluster_fn(Wl, k) for Wl in W]
    ari_solo = [float(adjusted_rand_score(fused_labels, sl)) for sl in solo_labels]

    delta_loo: List[float] = []
    for l in range(m):
        if m == 1:
            delta_loo.append(1.0)
            continue
        W_minus = [W[u] for u in range(m) if u != l]
        labels_minus = cluster_fn(fuse_fn(W_minus), k)
        # ARI can go negative, so 1 - ARI can exceed 1; clip so the weights are
        # a genuine partition of influence rather than an unbounded score.
        delta_loo.append(float(np.clip(
            1.0 - adjusted_rand_score(fused_labels, labels_minus), 0.0, 1.0)))

    eff: Dict[str, Any] = effective_n_layers(delta_loo)

    pvals: List[Optional[float]] = [None] * m
    null_means: List[Optional[float]] = [None] * m
    if n_perm and m >= 2:
        n = layer_mats[0].shape[0]
        for l in range(m):
            null_deltas = []
            for _ in range(n_perm):
                perm = rng.permutation(n)
                W_perm = list(W)
                W_perm[l] = kernel_fn(layer_mats[l][perm])
                labels_perm_full = cluster_fn(fuse_fn(W_perm), k)
                W_perm_minus = [W_perm[u] for u in range(m) if u != l]
                labels_perm_minus = cluster_fn(fuse_fn(W_perm_minus), k)
                null_deltas.append(
                    1.0 - adjusted_rand_score(labels_perm_full, labels_perm_minus)
                )
            deltas = np.asarray(null_deltas, dtype=float)
            # Phipson & Smyth (2010) corrected permutation p-value:
            # never zero, (b + 1) / (n_perm + 1).
            b = int((deltas >= delta_loo[l]).sum())
            pvals[l] = float((b + 1) / (n_perm + 1))
            null_means[l] = float(deltas.mean())

    per_layer = [
        {
            "layer": layer_names[l],
            "ari_solo_vs_fused": ari_solo[l],
            "delta_loo": delta_loo[l],
            "weight": eff["weights"][l],
            "p_contribution": pvals[l],
            "null_mean_delta_loo": null_means[l],
            "inert": bool(delta_loo[l] <= 1e-12),
        }
        for l in range(m)
    ]

    # delta_loo alone cannot tell "every layer says the same thing, so removing
    # any one changes nothing" from "no layer says anything". Both give
    # delta_loo ~ 0 and a flat n_eff. The solo agreements separate them: in the
    # first case every layer reproduces the fused partition on its own, in the
    # second none does.
    all_inert = bool(m >= 2 and max(delta_loo) <= 0.05)
    max_solo = max(ari_solo) if ari_solo else 0.0
    if all_inert and max_solo >= 0.5:
        regime = "redundant"          # layers agree; the fusion is not needed
    elif all_inert:
        regime = "uninformative"      # no layer carries the fused partition
    elif m >= 2 and eff["n_eff"] < collapse_threshold:
        regime = "collapsed"          # one layer drives the result
    elif m >= 2 and max(delta_loo) > 0.5 and max_solo < unstructured_solo_max:
        # Every layer moves the answer a lot and none reproduces it on its own.
        # That is what pure noise looks like: the partition is unstable, so
        # removing any layer flips it, and n_eff approaches m for the worst
        # possible reason. Without this branch a fusion of three noise layers is
        # labelled "complementary" with n_eff 2.81 of 3.
        #
        # The 0.35 cut-off is a heuristic, calibrated on the constructed cases in
        # tests/test_diagnostics.py: three noise layers give max solo agreement
        # 0.28, whereas the Klebsiella case study gives 0.63 and every planted
        # case gives 1.00. It is exposed as `unstructured_solo_max` because a
        # panel with genuinely complementary but individually weak layers could
        # sit below it, and such a run should be inspected rather than trusted
        # either way.
        regime = "unstructured"
    else:
        regime = "complementary"
    return {
        "per_layer": per_layer,
        "n_layers": m,
        "n_eff": eff["n_eff"],
        "n_eff_simpson": eff["n_eff_simpson"],
        "max_ari_solo_vs_fused": float(max_solo),
        "regime": regime,
        "regime_note": {
            "redundant": "every layer reproduces the fused partition alone; the "
                         "fusion adds nothing over the best single layer",
            "uninformative": "no layer reproduces the fused partition and "
                             "removing any layer changes nothing; the partition "
                             "is not supported by any input view",
            "collapsed": "one layer drives the fused partition; this is a "
                         "single-layer clustering presented as an integrated one",
            "unstructured": "removing any layer changes the partition and no "
                            "layer reproduces it on its own; the fused partition "
                            "is unstable rather than jointly supported, which is "
                            "what a fusion of uninformative layers looks like",
            "complementary": "layers contribute differently and jointly",
        }[regime],
        "collapse": bool(regime in ("collapsed", "uninformative", "unstructured")),
        "collapse_threshold": float(collapse_threshold),
        "n_perm": int(n_perm),
        "dominant_layer": layer_names[int(np.argmax(delta_loo))] if m else None,
    }
