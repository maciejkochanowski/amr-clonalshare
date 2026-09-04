"""phenotype.py — score a discovered partition against measured phenotype.

The validation gap this closes
------------------------------
A partition scored against a summary of its own input is not an external
validation. Kleborate's ``virulence_score``
and ``resistance_score`` are reconstructible from the shipped layers at 99.7 %
and 97.8 %, so "the virulence layer alone recovers the published virulence score
at ARI 0.875" is close to a tautology. A partition of genotypes scored against a
summary of the same genotypes tells you the arithmetic is consistent; it tells
you nothing about whether the partition means anything in a laboratory.

Measured antimicrobial susceptibility is not a function of the acquired-gene
panel. It is produced by a different instrument, on a different day, and it is
sensitive to mechanisms the panel cannot see -- *mgrB* inactivation for
colistin, *ramR*/*tet(A)* efflux for tigecycline, *gyrA*/*parC* for
fluoroquinolones, porin loss for carbapenems. That makes it the anchor.

What is actually being tested
-----------------------------
Not "is the partition associated with resistance" -- with 918 versus 582
isolates split almost exactly on "carries any acquired determinant", it must be.
The question worth asking is **whether the multi-layer partition beats the
one-bit rule it is nearly identical to**. So every statistic here is computed
for the partition *and* for each competing rule supplied by the caller, on
exactly the same isolates, and the comparison is reported as a difference.

Two conventions, both stated rather than defaulted silently:

* ``intermediate`` decides whether CLSI/EUCAST "Intermediate" counts as
  non-susceptible. ``"non_susceptible"`` (the surveillance convention, and the
  default) groups I with R; ``"drop"`` excludes those isolates; ``"susceptible"``
  groups I with S. The choice moves prevalence, so it is echoed into the result.
* Balanced accuracy, not accuracy. Resistance prevalence here ranges from 0.011
  (tigecycline) to 1.000 (ampicillin, which is intrinsic in *K. pneumoniae* and
  therefore uninformative). Plain accuracy on a 1 %-prevalence antibiotic is
  0.99 for the rule "call everything susceptible".

An antibiotic on which every isolate has the same phenotype carries no
information about any partition and is reported with ``status="uninformative"``
rather than being silently dropped or scored.
"""
from __future__ import annotations

from typing import Any, Dict, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from .stats import benjamini_hochberg, fisher_exact_p

__all__ = ["to_non_susceptible", "phenotype_concordance"]

_SIR = {"resistant": 1, "intermediate": None, "susceptible": 0,
        "nonsusceptible": 1, "non-susceptible": 1}


def to_non_susceptible(long_df: pd.DataFrame, *, id_column: str = "Strain_ID",
                       antibiotic_column: str = "antibiotic",
                       call_column: str = "resistant_phenotype",
                       intermediate: str = "non_susceptible") -> pd.DataFrame:
    """Long AST records -> wide isolate x antibiotic matrix of 0/1/NaN.

    ``NaN`` means "not tested", which is the overwhelming majority of cells and
    must stay distinguishable from "tested and susceptible".
    """
    if intermediate not in ("non_susceptible", "drop", "susceptible"):
        raise ValueError(f"unknown intermediate policy {intermediate!r}")
    df = long_df[[id_column, antibiotic_column, call_column]].copy()
    df[call_column] = df[call_column].astype(str).str.strip().str.lower()
    mapped = df[call_column].map(_SIR)
    inter = df[call_column].eq("intermediate")
    if intermediate == "non_susceptible":
        mapped = mapped.where(~inter, 1.0)
    elif intermediate == "susceptible":
        mapped = mapped.where(~inter, 0.0)
    df["_ns"] = mapped
    df = df.dropna(subset=["_ns"])
    # An isolate can carry several records for one drug (repeat testing or
    # several typing methods). Any non-susceptible result wins, which is the
    # conservative reading for surveillance.
    wide = df.pivot_table(index=id_column, columns=antibiotic_column,
                          values="_ns", aggfunc="max")
    return wide


def _balanced_accuracy(rule: np.ndarray, truth: np.ndarray) -> float:
    """Best of the two label orientations; a partition has no intrinsic polarity.

    Taking the maximum over orientations is the right thing -- cluster 0 is not
    "the resistant one" by construction -- but it is a selection, and it lifts
    the value achievable under a pure null from 0.50 to about 0.55. That floor is
    measured, not assumed: :func:`orientation_null_floor` permutes the labels
    within tested isolates and reports what this statistic returns on noise, and
    every result carries it.
    """
    best = 0.0
    for r in (rule, 1 - rule):
        tp = float(((r == 1) & (truth == 1)).sum())
        fn = float(((r == 0) & (truth == 1)).sum())
        tn = float(((r == 0) & (truth == 0)).sum())
        fp = float(((r == 1) & (truth == 0)).sum())
        if (tp + fn) == 0 or (tn + fp) == 0:
            continue
        best = max(best, 0.5 * (tp / (tp + fn) + tn / (tn + fp)))
    return best


def orientation_null_floor(rule: np.ndarray, truth: np.ndarray, *,
                           n_perm: int = 400, seed: int = 0) -> float:
    """What :func:`_balanced_accuracy` returns when the rule carries no signal.

    Permutes the rule within the tested isolates, preserving n, prevalence and
    the rule's own marginal. The result is well above 0.5 because of the
    orientation maximum, so a raw balanced accuracy of 0.55 is *chance*.
    """
    rng = np.random.default_rng(seed)
    r = np.asarray(rule)
    return float(np.mean([_balanced_accuracy(rng.permutation(r), truth)
                          for _ in range(int(n_perm))]))


def _mean_ba(rules: Dict[str, np.ndarray], pheno: np.ndarray,
             cols: Sequence[int], rows: np.ndarray) -> Dict[str, float]:
    """Mean balanced accuracy per rule over ``cols``, on the isolates ``rows``."""
    out = {}
    for name, r in rules.items():
        vals = []
        for j in cols:
            col = pheno[rows, j]
            tested = np.isfinite(col)
            truth = col[tested].astype(int)
            if truth.size < 4 or truth.min() == truth.max():
                continue
            vals.append(_balanced_accuracy(r[rows][tested], truth))
        out[name] = float(np.mean(vals)) if vals else float("nan")
    return out


def _bootstrap_head_to_head(rules: Dict[str, np.ndarray], pheno: np.ndarray,
                            cols: Sequence[int], n_boot: int, seed: int
                            ) -> Dict[str, Dict[str, Any]]:
    """Percentile CI for mean(BA[competitor]) - mean(BA[partition]).

    The resampling unit is the **isolate**, not the antibiotic: antibiotics are
    measured on the same isolates and co-resistance is the rule, so treating
    them as independent would understate the uncertainty badly. Without this,
    "the partition loses by 0.012" is a point estimate presented as a fact, and
    it is what the CLI gate fires on.
    """
    n = pheno.shape[0]
    rng = np.random.default_rng(seed)
    diffs: Dict[str, list] = {k: [] for k in rules if k != "partition"}
    for _ in range(int(n_boot)):
        rows = rng.integers(0, n, size=(n,))
        m = _mean_ba(rules, pheno, cols, rows)
        for k in diffs:
            if np.isfinite(m.get(k, np.nan)) and np.isfinite(m.get("partition", np.nan)):
                diffs[k].append(m[k] - m["partition"])
    out: Dict[str, Dict[str, Any]] = {}
    for k, v in diffs.items():
        a = np.asarray(v, dtype=float)
        if a.size < 10:
            out[k] = {"status": "insufficient_bootstrap"}
            continue
        lo, hi = np.percentile(a, [2.5, 97.5])
        out[k] = {
            "mean_difference_minus_partition": float(a.mean()),
            "ci95_low": float(lo), "ci95_high": float(hi),
            "prob_competitor_better": float((a > 0).mean()),
            "significantly_better_than_partition": bool(lo > 0.0),
            "n_boot": int(a.size),
        }
    return out


def phenotype_concordance(labels: "Sequence[int] | np.ndarray", strain_ids: Sequence[str],
                          phenotype: pd.DataFrame, *,
                          competing_rules: Optional[Mapping[str, "Sequence[int] | np.ndarray"]] = None,
                          min_tested: int = 20,
                          q_fdr: float = 0.05,
                          intermediate: str = "non_susceptible",
                          n_boot: int = 300,
                          n_perm_floor: int = 200,
                          seed: int = 0,
                          strata: Optional[Sequence] = None
                          ) -> Dict[str, Any]:
    """Does the partition predict measured non-susceptibility -- better than what?

    Parameters
    ----------
    labels : cluster assignment per isolate.
    strain_ids : isolate identifiers, aligned with ``labels``.
    phenotype : wide isolate x antibiotic matrix of 0/1/NaN
        (:func:`to_non_susceptible`). Indexed by isolate id; isolates absent
        from the index are simply untested.
    competing_rules : named binary rules over the same isolates, each of which
        the partition should be compared against. The rule that matters is
        "carries at least one acquired determinant": if the partition does not
        beat it, the fusion has not bought anything a single column could not.
    min_tested : antibiotics with fewer tested isolates are reported but not
        included in the FDR family.

    Returns
    -------
    dict with a per-antibiotic table, the FDR-significant set, and a head-to-head
    summary of mean balanced accuracy for the partition and each competing rule.
    """
    lab = np.asarray(list(labels))
    ids = pd.Index([str(s) for s in strain_ids])
    if lab.size != ids.size:
        raise ValueError("labels and strain_ids must be the same length")
    if not isinstance(phenotype, pd.DataFrame) or phenotype.empty:
        return {"status": "skipped", "reason": "no phenotype matrix supplied"}

    pheno = phenotype.reindex(ids)
    rules = {"partition": lab}
    for name, given in (competing_rules or {}).items():
        arr = np.asarray(list(given))
        if arr.size != lab.size:
            raise ValueError(f"competing rule {name!r} has the wrong length")
        rules[name] = arr

    rows = []
    for ab in pheno.columns:
        col = pheno[ab].to_numpy(dtype=float)
        tested = np.isfinite(col)
        n = int(tested.sum())
        truth = col[tested].astype(int)
        if n < 4 or truth.min() == truth.max():
            rows.append({"antibiotic": str(ab), "n_tested": n,
                         "n_non_susceptible": int(truth.sum()) if n else 0,
                         "status": "uninformative",
                         "reason": ("every tested isolate has the same phenotype"
                                    if n >= 4 else "too few tested isolates")})
            continue
        sub = lab[tested]
        # 2 x 2 of the largest cluster against the rest
        uq, cnt = np.unique(sub, return_counts=True)
        top = (sub == uq[np.argmax(cnt)]).astype(int)
        a = int(((top == 1) & (truth == 1)).sum())
        b = int(((top == 1) & (truth == 0)).sum())
        c = int(((top == 0) & (truth == 1)).sum())
        d = int(((top == 0) & (truth == 0)).sum())
        rec = {"antibiotic": str(ab), "n_tested": n,
               "n_non_susceptible": int(truth.sum()),
               "frac_non_susceptible": float(truth.mean()),
               "status": "ok",
               "p_fisher": float(fisher_exact_p(a, b, c, d)),
               "in_fdr_family": bool(n >= min_tested)}
        for name, r in rules.items():
            rsub = r[tested]
            uq2, cnt2 = np.unique(rsub, return_counts=True)
            binr = (rsub == uq2[np.argmax(cnt2)]).astype(int)
            rec[f"balanced_accuracy_{name}"] = _balanced_accuracy(binr, truth)
        rows.append(rec)

    table = pd.DataFrame(rows)
    fam = table[(table.get("status") == "ok") & table.get("in_fdr_family", False)]
    q = np.full(len(table), np.nan)
    if len(fam):
        qv, rej = benjamini_hochberg(fam["p_fisher"].to_numpy(), q=q_fdr)
        q[fam.index.to_numpy()] = qv
        table["q_bh"] = q
        table["fdr_significant"] = False
        table.loc[fam.index, "fdr_significant"] = rej
    else:
        table["q_bh"] = np.nan
        table["fdr_significant"] = False

    # Two families, reported separately because they give different answers
    # and averaging over one while naming the other is an easy error.
    # `fdr_family` is every antibiotic with enough tested isolates;
    # `fdr_significant` is the subset where the partition itself is significant,
    # which is a *selected* set and therefore flatters the partition.
    sig = table[table["fdr_significant"].astype(bool)] if len(table) else table
    head_to_head, head_to_head_significant = {}, {}
    if len(fam):
        for name in rules:
            col = f"balanced_accuracy_{name}"
            head_to_head[name] = float(table.loc[fam.index, col].mean())
            if len(sig):
                head_to_head_significant[name] = float(sig[col].mean())

    pheno_arr = pheno.to_numpy(dtype=float)
    fam_cols = [list(pheno.columns).index(r["antibiotic"])
                for _, r in table.loc[fam.index].iterrows()]
    boot = (_bootstrap_head_to_head(rules, pheno_arr, fam_cols, n_boot, seed)
            if len(fam) and n_boot else {})

    floors = {}
    if len(fam):
        for name, r in rules.items():
            vals = []
            for j in fam_cols:
                col = pheno_arr[:, j]
                tested = np.isfinite(col)
                truth = col[tested].astype(int)
                if truth.size >= 4 and truth.min() != truth.max():
                    vals.append(orientation_null_floor(r[tested], truth,
                                                       n_perm=n_perm_floor,
                                                       seed=seed))
            floors[name] = float(np.mean(vals)) if vals else float("nan")

    # "The partition loses" is only a finding if the loss survives resampling.
    beats = {n: (not boot.get(n, {}).get("significantly_better_than_partition",
                                         False))
             for n in head_to_head if n != "partition"}

    # Leave-one-collection-out. Public AST is deposited in blocks: on this
    # cohort a single 67-isolate Caribbean collection is 32 % of the phenotyped
    # subset and 0 % of the rest, and dropping it takes the number of
    # FDR-significant antibiotics from 12 to 1. A comparison that one collection
    # can reverse is a statement about that collection, so if strata are
    # supplied every head-to-head is recomputed with each stratum held out and
    # the worst case is reported next to the headline.
    loco: Dict[str, Any] = {}
    if strata is not None:
        st = pd.Series([str(x) for x in strata], index=ids)
        for g, cnt in st.value_counts().items():
            if cnt < 5:
                continue
            held = np.flatnonzero((st != g).to_numpy())
            if held.size < 20:
                continue
            m = _mean_ba(rules, pheno_arr, fam_cols, held)
            loco[str(g)] = {
                "n_held_out": int(cnt),
                "mean_balanced_accuracy": {k: (float(v) if np.isfinite(v) else None)
                                           for k, v in m.items()},
                "partition_still_loses_to": [
                    k for k, v in m.items()
                    if k != "partition" and np.isfinite(v)
                    and np.isfinite(m["partition"]) and v > m["partition"]],
            }
        flips = [g for g, r in loco.items()
                 if set(r["partition_still_loses_to"]) != set(
                     k for k in head_to_head
                     if k != "partition" and head_to_head[k] > head_to_head["partition"])]
        loco = {"per_stratum": loco,
                "strata_whose_removal_changes_the_verdict": flips,
                "verdict_is_stable": not flips,
                "note": ("public AST is deposited in collections, so a "
                         "head-to-head that one collection can reverse is a "
                         "statement about that collection rather than about "
                         "the cohort")}

    return {
        "status": "ok",
        "intermediate_policy": intermediate,
        "n_isolates_with_any_phenotype": int(np.isfinite(
            pheno.to_numpy(dtype=float)).any(axis=1).sum()),
        "n_isolates_total": int(ids.size),
        "n_antibiotics": int(pheno.shape[1]),
        "min_tested_for_fdr": int(min_tested),
        "per_antibiotic": table.to_dict("records"),
        "n_in_fdr_family": int(len(fam)),
        "n_fdr_significant": int(table["fdr_significant"].sum()),
        "mean_balanced_accuracy": head_to_head,
        "mean_balanced_accuracy_over_fdr_family": head_to_head,
        "mean_balanced_accuracy_over_significant_only": head_to_head_significant,
        "orientation_null_floor": floors,
        "head_to_head_bootstrap": boot,
        "leave_one_stratum_out": loco or None,
        "partition_beats": beats,
        "how_to_read_this": (
            "mean_balanced_accuracy is over the FDR *family* (every antibiotic "
            "with >= min_tested tested isolates), not over the significant "
            "subset; the significant subset is selected on the partition's own "
            "p-value and flatters it, and both are reported. "
            "orientation_null_floor is what the statistic returns on permuted "
            "labels -- around 0.55, not 0.50, because balanced accuracy is "
            "maximised over both label orientations -- so a raw value near 0.55 "
            "is chance. partition_beats is False only when the competitor's "
            "advantage survives an isolate-level bootstrap; a point difference "
            "of a couple of percent does not."),
        "what_this_tests": (
            "whether the discovered partition predicts laboratory-measured "
            "non-susceptibility, and whether it does so better than the "
            "competing one-column rules. A partition that does not beat "
            "'carries any acquired determinant' has added nothing that a single "
            "indicator column did not already provide."),
    }
