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

from typing import Dict, Optional

import numpy as np
import pandas as pd

from .attribution import SUPPORT_THRESHOLD, _codes

__all__ = ["MIN_MINOR_COUNT", "MIN_GROUP_SIZE", "group_adequacy",
           "trait_adequacy", "input_qc", "render_markdown"]

#: Fewest isolates of the rarer outcome for an antimicrobial to be reported
#: as adequate. This is a reporting convention fixed before the real-data
#: campaigns of the article, where the same rule declares a cell void. It
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
        "estimable": bool(n > 0 and support >= SUPPORT_THRESHOLD),
        "effective_group_size": float(n0),
        "group_sizes": {str(k): int(v) for k, v in sizes.items()},
        "smallest_groups": {str(k): int(v)
                            for k, v in sizes.sort_values().head(_LIST_CAP).items()},
    }


def _is_missing(v) -> bool:
    return (v is None or v != v or str(v).strip() == ""
            or str(v).lower() in ("nan", "na", "none"))


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
        out[str(col)] = {"n": n, "prevalence": float(ones / n) if n else 0.0,
                         "minor_count": int(minor), "constant": bool(minor == 0),
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
          "traits": trait_adequacy(panel),
          "min_minor_count": MIN_MINOR_COUNT, "min_group_size": MIN_GROUP_SIZE}
    if mic_join:
        qc["mic_join"] = mic_join
    if metadata is not None:
        sid = ids.astype(str)
        meta = metadata.set_axis(metadata.index.astype(str))
        joined = sid.isin(meta.index)
        qc["metadata_join"] = {"n_joined": int(joined.sum()),
                               "share_joined": float(joined.mean()) if len(sid) else 0.0,
                               "unjoined_examples": [str(i) for i in sid[~joined][:_LIST_CAP]]}
        if lineage_column and lineage_column in meta.columns:
            lin = meta.reindex(sid)[lineage_column]
            qc["lineage"] = {"column": lineage_column, **group_adequacy(lin)}
    return qc


def _pct(x: float) -> str:
    return f"{100 * x:.1f} %"


def render_markdown(qc: dict) -> str:
    """The same record in plain language, for a reader who is not a statistician."""
    lines = ["# Input check", ""]
    lines.append(f"Isolates in the phenotype table: {qc['n_isolates']}; "
                 f"antimicrobials with a readable call: {qc['n_antimicrobials']}.")
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
        weak = [c for c, t in traits.items() if not t["adequate"] and not t["constant"]]
        const = [c for c, t in traits.items() if t["constant"]]
        ok = len(traits) - len(weak) - len(const)
        lines.append(
            f"Antimicrobials with at least {qc['min_minor_count']} isolates of the "
            f"rarer outcome: {ok} of {len(traits)}. Below that number the estimate "
            "is still computed, but its interval will be wide and the report says so.")
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
               "The clonal share will be reported as not estimable at this typing "
               "resolution; a coarser lineage definition (for example serovar "
               "instead of SNP cluster) raises support."))
        lines.append(f"Effective lineage size for the between-lineage variance: "
                     f"{g['effective_group_size']:.1f}.")
        lines.append("")
    return "\n".join(lines)
