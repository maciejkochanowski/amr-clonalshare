"""lineage.py — population-structure diagnostics for discovered partitions.

Bacterial isolates are not independent draws. Accessory-genome content is
strongly structured by clonal descent, so any partition of isolates by trait
content will partly recapitulate the phylogeny rather than recurrent trait
combinations. Separating a lineage effect from a locus effect is the defining
problem of bacterial genome-wide association analysis (Earle et al. 2016, Nat
Microbiol 1:16041, doi:10.1038/nmicrobiol.2016.41; Collins & Didelot 2018, PLoS
Comput Biol 14:e1005958, doi:10.1371/journal.pcbi.1005958; Lees et al. 2016,
Nat Commun 7:12797, doi:10.1038/ncomms12797), and an unsupervised
archetype analysis inherits it in full.

Two things must be reported, and neither is optional:

``lineage_concordance``
    A **pair-based** test. For every pair of isolates from the same lineage,
    how often are they placed in the same cluster, against a permutation null
    that reshuffles the cluster labels? Mutual-information scores such as AMI
    are the wrong instrument here: with hundreds of singleton lineages the
    chance-correction dominates and AMI is small even when same-lineage pairs
    are co-clustered far more often than chance. On the *Klebsiella* example
    AMI(cluster, ST) = 0.039 - which invites the conclusion "not a clonal
    artefact" - while the same-ST pair concordance is 0.701 against a null of
    0.524 (z = +14.0). The pair statistic is the one that answers the question.
    (Those four numbers are ``metadata_diagnostics.lineage_concordance`` of the
    committed ``examples/klebsiella/expected/`` run.)

    The statistic's *coverage* is reported alongside it and should be read
    first. Lineage labels are taken verbatim, so ``ST17-1LV`` is a different
    lineage from ``ST17`` and isolates in singleton lineages contribute no
    pairs at all: on this heavily de-replicated cohort the whole statistic
    rests on 1600 same-ST pairs, 0.14 % of the 1 124 250 pairs in the cohort.
    ``collapse_variants=True`` folds single- and double-locus variants into
    their parent ST, which here raises the coverage to 1978 pairs and the
    concordance to 0.718 against the same null - i.e. the verbatim labelling is
    conservative rather than wrong, but a reader has to be able to check that.

``dereplication_sensitivity``
    The direct sensitivity analysis: re-run the whole analysis on one isolate
    per lineage and report how much the answer moves. If the number of
    archetypes or the partition changes materially, the original result is a
    function of how densely each lineage happened to be sampled - which, for a
    surveillance collection, is an artefact of submission behaviour.

Neither diagnostic requires a phylogeny; a sequence-type or clonal-group label
is enough. For a tree-aware test of whether a *trait combination* recurs
independently, see :mod:`amr_clonalshare.archephy`.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

import re
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

from .stats import permutation_pvalue

__all__ = ["lineage_concordance", "cluster_composition", "dereplicate_index",
           "collapse_st_variants"]


def _same_group_pairs(codes: np.ndarray) -> int:
    _, counts = np.unique(codes, return_counts=True)
    return int((counts * (counts - 1) // 2).sum())


def _concordance(codes: np.ndarray, labels: np.ndarray) -> float:
    """Fraction of same-lineage pairs that are also same-cluster."""
    total = 0
    same = 0
    for g in np.unique(codes):
        idx = np.where(codes == g)[0]
        if idx.size < 2:
            continue
        lab = labels[idx]
        _, cnt = np.unique(lab, return_counts=True)
        n = idx.size
        total += n * (n - 1) // 2
        same += int((cnt * (cnt - 1) // 2).sum())
    return float(same / total) if total else float("nan")


_VARIANT_SUFFIX = re.compile(r"[-_ ]?\d*(SLV|DLV|LV)$", re.IGNORECASE)


def collapse_st_variants(lineage: "Sequence | np.ndarray") -> np.ndarray:
    """Fold single- and double-locus variants onto their parent sequence type.

    ``ST17-1LV`` and ``ST17-2LV`` become ``ST17``. MLST variant notation is what
    makes a verbatim lineage label under-count: a variant of a common clone is a
    different string but the same lineage for the purpose of asking whether the
    partition is clonal, so treating them separately shrinks the number of
    informative pairs and biases the diagnostic towards "not confounded".
    """
    return np.array([_VARIANT_SUFFIX.sub("", str(x)) if x is not None else x
                     for x in lineage], dtype=object)


def lineage_concordance(labels: "Sequence[int] | np.ndarray", lineage: "Sequence | np.ndarray",
                        *, n_perm: int = 500,
                        rng: Optional[np.random.Generator] = None,
                        collapse_variants: bool = False
                        ) -> Dict[str, Any]:
    """Test whether same-lineage isolates are co-clustered more than chance.

    Parameters
    ----------
    labels : cluster assignment per isolate.
    lineage : lineage label per isolate (ST, clonal group, ...). Missing values
        are dropped.
    n_perm : permutations of the cluster labels for the null.
    collapse_variants : fold ``ST17-1LV`` onto ``ST17`` first
        (:func:`collapse_st_variants`).

    Returns
    -------
    dict with the observed same-lineage pair concordance, the null mean and sd,
    a z-score, a Phipson-Smyth permutation p-value, the number of informative
    pairs *and the fraction of all pairs they represent*, and AMI/ARI against
    lineage for reference (with a note that they are not the primary
    statistic). ``frac_pairs_informative`` is the coverage of the test: on a
    heavily de-replicated cohort it is well under 1 %, and a diagnostic that
    inspects 0.14 % of the pairs should say so next to its own z-score.
    """
    if rng is None:
        rng = np.random.default_rng()
    lab = np.asarray(list(labels))
    lin = pd.Series(list(collapse_st_variants(lineage) if collapse_variants
                         else lineage), dtype="object")
    ok = np.array(lin.notna().to_numpy(), dtype=bool, copy=True)
    lin_s = np.asarray(lin.fillna("__missing__").astype(str).to_numpy(), dtype=object)
    ok &= ~np.isin(lin_s, ("nan", "None", "", "__missing__"))
    lab, lin = lab[ok], lin_s[ok]
    if lab.size == 0:
        return {"status": "skipped", "reason": "no lineage labels"}
    codes = pd.factorize(lin)[0]
    n_pairs = _same_group_pairs(codes)
    if n_pairs == 0:
        return {"status": "skipped",
                "reason": "every lineage is a singleton; no same-lineage pairs",
                "n_isolates": int(lab.size),
                "n_lineages": int(np.unique(codes).size)}

    obs = _concordance(codes, lab)
    null = np.array([_concordance(codes, rng.permutation(lab))
                     for _ in range(int(n_perm))], dtype=float)
    sd = float(null.std(ddof=1)) if null.size > 1 else 0.0
    return {
        "status": "ok",
        "n_isolates": int(lab.size),
        "n_lineages": int(np.unique(codes).size),
        "n_same_lineage_pairs": int(n_pairs),
        "n_pairs_total": int(lab.size * (lab.size - 1) // 2),
        "frac_pairs_informative": float(
            n_pairs / (lab.size * (lab.size - 1) / 2)) if lab.size > 1 else 0.0,
        "collapse_variants": bool(collapse_variants),
        "concordance_observed": obs,
        "concordance_null_mean": float(null.mean()),
        "concordance_null_sd": sd,
        "z": float((obs - null.mean()) / sd) if sd > 0 else float("nan"),
        "p_value": permutation_pvalue(null, obs, tail="greater"),
        "n_perm": int(n_perm),
        "ami_vs_lineage": float(adjusted_mutual_info_score(lin, lab)),
        "ari_vs_lineage": float(adjusted_rand_score(lin, lab)),
        "note": "the pair concordance, not AMI/ARI, is the primary statistic: "
                "with many singleton lineages AMI is chance-corrected towards "
                "zero even under strong clonal confounding",
    }


def cluster_composition(labels: Sequence[int], metadata: pd.DataFrame,
                        columns: Sequence[str], *, top: int = 5
                        ) -> Dict[str, Any]:
    """Per-cluster composition over categorical metadata columns.

    Returns, for each cluster and each requested column, the top categories with
    counts and the maximum single-category share. A cluster whose top lineage
    share is high is a clone, not an archetype.
    """
    lab = np.asarray(list(labels))
    out: Dict[str, Any] = {}
    for c in np.unique(lab):
        mask = lab == c
        entry: Dict[str, Any] = {"n": int(mask.sum())}
        for col in columns:
            if col not in metadata.columns:
                continue
            vc = (metadata.loc[mask, col].fillna("__missing__")
                  .astype(str).value_counts())
            entry[col] = {
                "top": {str(k): int(v) for k, v in vc.head(top).items()},
                "n_distinct": int(vc.size),
                "max_share": float(vc.iloc[0] / mask.sum()) if vc.size else float("nan"),
            }
        out[str(int(c))] = entry
    return out


def dereplicate_index(lineage: Sequence, *, rng: Optional[np.random.Generator] = None,
                      strategy: str = "first") -> np.ndarray:
    """Row indices retaining one isolate per lineage.

    ``strategy="first"`` is deterministic and reproducible; ``"random"`` draws
    one member per lineage with ``rng`` and is appropriate when averaging the
    sensitivity analysis over several de-replications.
    """
    lin = (pd.Series(list(lineage), dtype="object")
           .fillna("__missing__").astype(str).to_numpy())
    keep: List[int] = []
    for g in pd.unique(lin):
        idx = np.where(lin == g)[0]
        if strategy == "random":
            if rng is None:
                rng = np.random.default_rng()
            keep.append(int(rng.choice(idx)))
        else:
            keep.append(int(idx[0]))
    return np.sort(np.asarray(keep, dtype=int))
