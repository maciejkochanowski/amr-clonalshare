"""Input quality control: what the loader found, and what the estimators will do with it.

The record answers the question a user asks before the run: is this cohort
large enough. There is no single number. The clonal-share estimator learns
one rate per lineage inside each fold, so what binds is the number of
lineages with at least two isolates (they alone carry within-lineage
information) and the count of the rarer outcome per antimicrobial. Both are
reported per group and per antimicrobial, against the thresholds the package
already ships, so the same rule that gates the estimate is visible before the
estimate is made.
"""
from __future__ import annotations

import html
from fractions import Fraction

from typing import Dict, Optional

import numpy as np
import pandas as pd

from .attribution import SUPPORT_THRESHOLD, _codes

__all__ = ["MIN_MINOR_COUNT", "MIN_GROUP_SIZE", "group_adequacy",
           "trait_adequacy", "input_qc", "render_markdown", "support_pairing_plan"]

#: Fewest isolates of the rarer outcome for an antimicrobial to be reported
#: as adequate. This is a reporting convention fixed before the real-data
#: campaigns, where the same rule declares a cell void. It
#: is a reporting threshold, not a gate: the
#: estimator still runs below it and its own bootstrap interval says how
#: little the data supports.
MIN_MINOR_COUNT = 20
#: Isolates a lineage needs before it carries within-lineage information. A
#: singleton lineage has one isolate and one outcome; it cannot say how much
#: the trait varies inside the lineage. ``support`` is the share of isolates
#: in lineages of at least this size, and the estimator refuses a cell whose
#: support falls below :data:`amr_clonalshare.attribution.SUPPORT_THRESHOLD`.
MIN_GROUP_SIZE = 2
_LIST_CAP = 10

FEASIBILITY_SCOPE = (
    "Each row uses only isolates with both a readable result for that agent and "
    "a recorded lineage. The last column is an algebraic support scenario: one "
    "genuinely new tested isolate in each of that many distinct singleton lineages. "
    "It does not guarantee interval precision, coverage, or overall estimability; "
    "it cannot repair a constant trait or a missing lineage contrast. Duplicating "
    "existing rows adds no evidence. The rarer-outcome count is a warning threshold, "
    "not an additional estimator gate."
)
FEASIBILITY_HEADERS = ("Agent", "Retained", "Support", "Input failures",
                       "Rarer outcome", "New isolates for support")


def support_pairing_plan(n: int, singletons: int,
                         threshold: float = SUPPORT_THRESHOLD) -> dict:
    """Minimal singleton pairings to meet support, not a precision calculation.

    For k different current singleton groups each gaining one new isolate,
    support becomes (n - singletons + 2*k)/(n + k). Decimal-rational arithmetic
    avoids a spurious extra isolate at exact threshold boundaries.
    """
    if n < 0 or singletons < 0 or singletons > n or not 0 < threshold < 1:
        raise ValueError("require 0 <= singletons <= n and 0 < threshold < 1")
    if not n:
        return {"additional_isolates": None, "projected_support": None}
    target = Fraction(str(threshold))
    needed = (singletons - (1 - target) * n) / (2 - target)
    k = max(0, -(-needed.numerator // needed.denominator))
    return {"additional_isolates": k,
            "projected_support": float(Fraction(n - singletons + 2*k, n + k))}


def _agent_feasibility(panel: pd.DataFrame, lineage: pd.Series) -> Dict[str, dict]:
    """Input conditions on the same finite-then-typed subset as clonal_share."""
    labels = lineage.to_numpy(dtype=object)
    typed = np.array([not _is_missing(v) for v in labels], dtype=bool)
    out = {}
    for name in panel:
        y = panel[name].to_numpy(dtype=float)
        finite = np.isfinite(y)
        keep = finite & typed
        values = y[keep]
        group = group_adequacy(pd.Series(labels[keep]))
        n = int(keep.sum())
        minor = int(min(values.sum(), n - values.sum())) if n else 0
        failures = []
        if n < 2:
            failures.append("fewer than 2 tested and typed isolates")
        if group["n_groups"] < 2:
            failures.append("fewer than 2 lineages")
        if n and np.ptp(values) == 0:
            failures.append("constant trait")
        if n >= 2 and group["support"] < SUPPORT_THRESHOLD:
            failures.append(f"support below {SUPPORT_THRESHOLD:.2f}")
        plan = support_pairing_plan(n, group["n_singletons"])
        out[str(name)] = {
            # both counts describe the cohort as supplied, as the estimators
            # count them; an unread and untyped isolate is in both
            "n_retained": n, "n_dropped_non_finite": int((~finite).sum()),
            "n_dropped_untyped": int((~typed).sum()),
            "prevalence": float(values.mean()) if n else None,
            "minor_count": minor, "n_groups": group["n_groups"],
            "n_singletons": group["n_singletons"],
            "support": group["support"] if n >= 2 else None,
            "support_threshold": SUPPORT_THRESHOLD,
            "input_feasible": not failures, "failure_reasons": failures,
            "warnings": ([f"fewer than {MIN_MINOR_COUNT} isolates of the rarer outcome"]
                         if minor < MIN_MINOR_COUNT else []),
            "support_pairings_needed": plan["additional_isolates"],
            "projected_support": plan["projected_support"],
        }
    return out


def feasibility_rows(records: dict) -> list:
    """Shared display rows for the input check and both run report formats."""
    return [(name, str(r["n_retained"]),
             _pct(r["support"]) if r["support"] is not None else "not defined",
             "; ".join(r["failure_reasons"]) or "none at input level",
             str(r["minor_count"]) + (f" (below {MIN_MINOR_COUNT})" if r["warnings"] else ""),
             str(r["support_pairings_needed"]) if r["support_pairings_needed"] is not None else "not applicable")
            for name, r in records.items()]


def group_adequacy(lineage: pd.Series) -> dict:
    """Per-group sizes and what the clonal-share estimator will do with them.

    Untyped isolates (empty, ``NaN``, ``"NA"``, ``"none"``) are set aside
    before the sizes are read, exactly as the estimator sets them aside before
    it fits, so ``support`` is the estimator's own quantity: the share of the
    typed isolates in lineages of at least :data:`MIN_GROUP_SIZE` members.
    The count of untyped isolates is reported beside it, because a label that
    a third of the cohort lacks is a finding of its own.
    """
    values = np.asarray(lineage, dtype=object)
    n_all = int(values.size)
    missing = np.array([_is_missing(v) for v in values], dtype=bool)
    n_untyped = int(missing.sum())
    typed = values[~missing]
    n = int(typed.size)
    codes = _codes(typed) if n else np.zeros(0, dtype=int)
    counts = np.bincount(codes) if n else np.zeros(0, dtype=int)
    names: Dict[int, str] = {}
    for c, v in zip(codes, typed):
        names.setdefault(int(c), str(v))
    sizes = pd.Series({names[i]: int(counts[i]) for i in range(counts.size)},
                      dtype=int)
    informative = sizes[sizes >= MIN_GROUP_SIZE]
    support = float(informative.sum() / n) if n else 0.0
    n_groups = int(sizes.size)
    n0 = ((n - float((sizes ** 2).sum()) / n) / (n_groups - 1)
          if n_groups > 1 and n else float("nan"))
    return {
        "n": n_all,
        "n_typed": n,
        "n_untyped": n_untyped,
        "n_groups": n_groups,
        "n_singletons": int((sizes == 1).sum()),
        "n_informative_groups": int(informative.size),
        "support": support,
        "support_threshold": SUPPORT_THRESHOLD,
        "estimable": bool(n_groups >= 2 and support >= SUPPORT_THRESHOLD),
        "effective_group_size": float(n0),
        "group_sizes": {str(k): int(v) for k, v in sizes.items()},
        "smallest_groups": {str(k): int(v)
                            for k, v in sizes.sort_values().head(_LIST_CAP).items()},
    }


def _is_missing(v) -> bool:
    return bool(pd.isna(v) or str(v).strip().lower() in ("", "nan", "na", "none", "n/a", "<na>"))


def trait_adequacy(X: pd.DataFrame) -> Dict[str, dict]:
    """Per-trait counts of the rarer outcome against :data:`MIN_MINOR_COUNT`."""
    out = {}
    for col in X.columns:
        y = X[col].to_numpy(dtype=float)
        # an untested cell is not a susceptible one: the counts are taken over
        # the isolates that carry a result for this agent
        y = y[np.isfinite(y)]
        n = int(y.size)
        ones = int(y.sum())
        minor = min(ones, n - ones)
        out[str(col)] = {"n": n, "n_missing": len(X) - n, "n_positive": ones,
                         "prevalence": float(ones / n) if n else None,
                         "minor_count": int(minor), "constant": bool(n > 0 and minor == 0),
                         "empty": n == 0,
                         "adequate": bool(minor >= MIN_MINOR_COUNT)}
    return out


def input_qc(panel: pd.DataFrame, *,
             metadata: Optional[pd.DataFrame] = None,
             lineage_column: Optional[str] = None,
             mic_join: Optional[dict] = None) -> dict:
    """The full input record written with every run and by ``--check-input``."""
    ids = panel.index
    qc = {"n_isolates": int(len(ids)),
          "n_antimicrobials": int(panel.shape[1]),
          "n_empty_agents": int(panel.isna().all(axis=0).sum()),
          "traits": trait_adequacy(panel),
          "min_minor_count": MIN_MINOR_COUNT, "min_group_size": MIN_GROUP_SIZE}
    if mic_join:
        qc["mic_join"] = mic_join
    if metadata is not None:
        sid = ids
        meta = metadata
        joined = sid.isin(meta.index)
        qc["metadata_join"] = {"n_joined": int(joined.sum()),
                               "share_joined": float(joined.mean()) if len(sid) else 0.0,
                               "unjoined_examples": [str(i) for i in sid[~joined][:_LIST_CAP]]}
        if lineage_column and lineage_column in meta.columns:
            lin = meta.reindex(sid)[lineage_column]
            qc["lineage"] = {"column": lineage_column, **group_adequacy(lin)}
            qc["agent_feasibility"] = _agent_feasibility(panel, lin)
    return qc


def _pct(x: float) -> str:
    return f"{100 * x:.1f} %"


def render_markdown(qc: dict) -> str:
    """The same record in plain language, for a reader who is not a statistician."""
    lines = ["# Input check", ""]
    lines.append(f"Isolates in the raw call/MIC cohort: {qc['n_isolates']}; "
                 f"antimicrobials recorded in the call table: {qc['n_antimicrobials']}.")
    reading = qc.get("phenotype_reading") or {}
    if reading:
        lines.append(
            f"Susceptibility rows read: {reading['rows_read']}; "
            f"{reading['n_unrecognized_calls']} unrecognized call(s), "
            f"{reading['n_missing_calls']} missing call(s), and "
            f"{reading['n_intermediate_dropped']} intermediate call(s) set aside. "
            f"Isolates without any readable call: {reading['isolates_without_readable_calls']}.")
    if "mic_join" in qc:
        j = qc["mic_join"]
        lines.append(f"Recorded dilutions: {j['rows_joined']} rows for "
                     f"{j['strains_with_mic']} of {j['strains_aligned']} isolates, "
                     f"{len(j['antimicrobials'])} antimicrobials.")
    lines.append("")
    traits = qc["traits"]
    if traits:
        lines.append("## Antimicrobials")
        lines.append("")
        empty = [c for c, t in traits.items() if t.get("empty")]
        weak = [c for c, t in traits.items() if not t["adequate"] and not t["constant"] and not t.get("empty")]
        const = [c for c, t in traits.items() if t["constant"]]
        ok = len(traits) - len(weak) - len(const) - len(empty)
        lines.append(
            f"Antimicrobials with at least {qc['min_minor_count']} isolates of the "
            f"rarer outcome: {ok} of {len(traits)}. This count is a warning threshold; "
            "each method reports whether its own input requirements are met.")
        if empty:
            lines.append(f"Without any readable call: {', '.join(empty[:_LIST_CAP])}.")
        if weak:
            lines.append(f"Below the threshold: {', '.join(weak[:_LIST_CAP])}"
                         + (" and others" if len(weak) > _LIST_CAP else "") + ".")
        if const:
            lines.append(f"With a single value (nothing to explain): "
                         f"{', '.join(const[:_LIST_CAP])}"
                         + (" and others" if len(const) > _LIST_CAP else "") + ".")
        lines.append("")
    if "metadata_join" in qc:
        j = qc["metadata_join"]
        lines.append("## Metadata")
        lines.append("")
        lines.append(f"{j['n_joined']} isolates ({_pct(j['share_joined'])}) have a "
                     "metadata row." + ("" if j["share_joined"] == 1.0 else
                     f" Identifiers without one include {j['unjoined_examples']}."))
        lines.append("")
    if "lineage" in qc:
        g = qc["lineage"]
        lines.append(f"## Lineage groups (`{g['column']}`)")
        lines.append("")
        lines.append(
            f"{g['n_groups']} lineages among {g['n'] - g['n_untyped']} typed isolates "
            f"({g['n_untyped']} untyped). {g['n_singletons']} lineages hold a single "
            "isolate; a single isolate cannot show how much a trait varies inside "
            "its lineage, so those isolates do not contribute to that part of "
            "the estimate.")
        lines.append(
            f"Share of isolates in lineages of at least {qc['min_group_size']} "
            f"(support): {_pct(g['support'])}; the estimator needs "
            f"{_pct(g['support_threshold'])}. "
            + ("The clonal share can be estimated on this cohort."
               if g["estimable"] else
               "The clonal share needs at least two distinct lineages; there "
               "is no between-lineage comparison on this cohort."
               if g["n_groups"] < 2 else
               "The clonal share will be reported as not estimable at this typing "
               "resolution; a coarser lineage definition (for example serovar "
               "instead of SNP cluster) raises support."))
        size = g.get("effective_group_size")
        effective = (f"{size:.1f}" if isinstance(size, (int, float))
                     and np.isfinite(size) else "not defined")
        lines.append(f"Effective lineage size for the between-lineage variance: {effective}.")
        lines.append("")
    feasibility = qc.get("agent_feasibility") or {}
    if feasibility:
        lines += ["", "## Per-agent feasibility", "", FEASIBILITY_SCOPE, "",
                  "| " + " | ".join(FEASIBILITY_HEADERS) + " |",
                  "| " + " | ".join(["---"] * len(FEASIBILITY_HEADERS)) + " |"]
        for row in feasibility_rows(feasibility):
            lines.append("| " + " | ".join(str(v).replace("|", "\\|").replace("\n", " ")
                                            for v in row) + " |")
    return html.escape("\n".join(lines), quote=False)
