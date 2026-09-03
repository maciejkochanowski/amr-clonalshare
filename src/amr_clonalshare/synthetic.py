"""synthetic.py — planted-truth generators for validation and calibration.

Two generators, deliberately different in what they are for:

``synth_cluster_archetypes``
    Clean planted archetypes with a tunable ``overlap``. Used by
    :func:`amr_clonalshare.core.validate` and by the recovery sweep. Setting
    ``k_true = 1`` gives a **null** cohort with no clusters at all, which is
    the case a k-selector and a post-clustering test must both get right and
    which a harness of planted clusters alone never generates.

``synth_lineage_cohort``
    A cohort whose trait content is generated *on a clonal frame*: isolates
    belong to lineages, traits are acquired at the lineage level with only a
    little within-lineage variation, and lineages are sampled unevenly. Every
    apparent "archetype" is then a lineage effect. This is the adversarial case
    for the whole approach - a pipeline that reports confident archetypes here
    is reporting sampling depth - and it is what
    :mod:`amr_clonalshare.lineage` exists to detect.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np
import pandas as pd

__all__ = ["synth_cluster_archetypes", "synth_lineage_cohort"]


def synth_cluster_archetypes(
    n: int = 150,
    p_amr: int = 20,
    p_vir: int = 40,
    k_true: int = 3,
    overlap: float = 0.05,
    seed: int = 42,
) -> Tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    """Two binary layers with ``k_true`` planted archetypes.

    Each archetype owns a disjoint block of features in each layer, present with
    probability ``1 - overlap``; every feature also fires at background rate
    ``overlap``. ``overlap -> 0.5`` erases the structure, so sweeping it gives a
    difficulty curve rather than a single easy setting.

    ``k_true = 1`` produces a null cohort: every feature is independent
    Bernoulli(``overlap``) and the returned labels are all zero.
    """
    rng = np.random.default_rng(seed)
    if k_true < 1:
        raise ValueError("k_true must be >= 1")
    idx = [f"SC_{i:04d}" for i in range(n)]
    amr_cols = [f"amr_{j:02d}" for j in range(p_amr)]
    vir_cols = [f"vir_{j:02d}" for j in range(p_vir)]

    if k_true == 1:
        amr = rng.binomial(1, overlap, (n, p_amr))
        vir = rng.binomial(1, overlap, (n, p_vir))
        labels = np.zeros(n, dtype=int)
        return (pd.DataFrame(amr, index=idx, columns=amr_cols),
                pd.DataFrame(vir, index=idx, columns=vir_cols), labels)

    sizes = np.full(k_true, n // k_true)
    sizes[-1] += n - sizes.sum()
    labels = np.concatenate([np.full(s, c) for c, s in enumerate(sizes)])
    amr = np.zeros((n, p_amr), dtype=int)
    vir = np.zeros((n, p_vir), dtype=int)
    amr_per, vir_per = p_amr // k_true, p_vir // k_true
    for c in range(k_true):
        mask = labels == c
        amr[mask, c * amr_per:(c + 1) * amr_per] = rng.binomial(
            1, 1 - overlap, (int(mask.sum()), amr_per))
        vir[mask, c * vir_per:(c + 1) * vir_per] = rng.binomial(
            1, 1 - overlap, (int(mask.sum()), vir_per))
    amr = np.clip(amr + rng.binomial(1, overlap, amr.shape), 0, 1)
    vir = np.clip(vir + rng.binomial(1, overlap, vir.shape), 0, 1)
    return (pd.DataFrame(amr, index=idx, columns=amr_cols),
            pd.DataFrame(vir, index=idx, columns=vir_cols), labels)


def synth_lineage_cohort(
    n: int = 400,
    n_lineages: int = 40,
    p_amr: int = 15,
    p_vir: int = 20,
    within_lineage_noise: float = 0.03,
    sampling_skew: float = 2.0,
    seed: int = 0,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """A cohort in which all trait structure is clonal.

    Each lineage draws one random trait profile; isolates of that lineage carry
    it with a small per-feature flip probability. Lineage sizes follow a
    power law with exponent ``sampling_skew``, mimicking a surveillance
    collection in which a few outbreak clones are sampled heavily.

    There are no archetypes here, only clones. Returns ``(amr, vir, metadata)``
    where ``metadata`` carries the ``lineage`` column.
    """
    rng = np.random.default_rng(seed)
    w = 1.0 / np.power(np.arange(1, n_lineages + 1), sampling_skew)
    w /= w.sum()
    assign = rng.choice(n_lineages, size=n, p=w)
    prof_amr = rng.binomial(1, 0.3, (n_lineages, p_amr))
    prof_vir = rng.binomial(1, 0.2, (n_lineages, p_vir))
    amr = prof_amr[assign]
    vir = prof_vir[assign]
    amr = np.where(rng.random(amr.shape) < within_lineage_noise, 1 - amr, amr)
    vir = np.where(rng.random(vir.shape) < within_lineage_noise, 1 - vir, vir)
    idx = [f"LC_{i:04d}" for i in range(n)]
    meta = pd.DataFrame({"lineage": [f"L{a:03d}" for a in assign]}, index=idx)
    return (pd.DataFrame(amr, index=idx, columns=[f"amr_{j:02d}" for j in range(p_amr)]),
            pd.DataFrame(vir, index=idx, columns=[f"vir_{j:02d}" for j in range(p_vir)]),
            meta)
