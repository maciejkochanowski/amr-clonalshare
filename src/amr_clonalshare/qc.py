"""Input quality control: what the loader found, and what the estimators will do with it.

A binary layer has two states, present and absent, and a cohort has a third:
not measured. Nothing downstream can tell a 0 that was measured from a 0 that
was written to fill a gap, so the distinction has to be made here, before any
matrix is cast to integers. This module records every missing cell, applies
the policy the configuration names, and refuses to impute.

The second half of the record answers the question a user asks before the
run: is this cohort large enough. There is no single number. The clonal-share
estimator learns one rate per lineage inside each fold, so what binds is the
number of lineages with at least two isolates (they alone carry within-lineage
information) and the count of the rarer outcome per trait. Both are reported
per group and per trait, against the thresholds the package already ships, so
the same rule that gates the estimate is visible before the estimate is made.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .attribution import SUPPORT_THRESHOLD, _codes
from .config import ConfigError

__all__ = ["MIN_MINOR_COUNT", "MIN_GROUP_SIZE", "MissingRecord", "layer_missing",
           "apply_missing_policy", "non_binary_values", "group_adequacy",
           "trait_adequacy", "input_qc", "render_markdown"]

#: Fewest isolates of the rarer outcome for a trait to be reported as adequate.
#: This is the project's convention, fixed before every real-data campaign in
#: the validation ledger, and the same rule that declares a cell void there.
#: It is a reporting threshold, not a gate: the estimator still runs below it
#: and its own bootstrap interval says how little the data supports.
MIN_MINOR_COUNT = 20
#: Isolates a lineage needs before it carries within-lineage information. A
#: singleton lineage has one isolate and one outcome; it cannot say how much
#: the trait varies inside the lineage. ``support`` is the share of isolates
#: in lineages of at least this size, and the estimator refuses a cell whose
#: support falls below :data:`amr_clonalshare.attribution.SUPPORT_THRESHOLD`.
MIN_GROUP_SIZE = 2
_LIST_CAP = 10


@dataclass
class MissingRecord:
    """Where the missing cells were and what was done about them."""
    policy: str
    layers: Dict[str, dict] = field(default_factory=dict)
    dropped_rows: List[str] = field(default_factory=list)
    dropped_columns: Dict[str, List[str]] = field(default_factory=dict)

    @property
    def any_missing(self) -> bool:
        return any(v["cells"] > 0 for v in self.layers.values())

    def as_dict(self) -> dict:
        return {"policy": self.policy, "layers": self.layers,
                "dropped_rows": self.dropped_rows,
                "dropped_columns": self.dropped_columns}


def layer_missing(df: pd.DataFrame) -> dict:
    """Count of empty cells in one layer, with the rows and columns that hold them."""
    mask = df.isna()
    cells = int(mask.to_numpy().sum())
    rows = mask.any(axis=1)
    cols = mask.any(axis=0)
    return {"cells": cells,
            "share": float(cells / mask.size) if mask.size else 0.0,
            "rows_with_missing": int(rows.sum()),
            "columns_with_missing": int(cols.sum()),
            "example_rows": [str(i) for i in df.index[rows][:_LIST_CAP]],
            "example_columns": [str(c) for c in df.columns[cols][:_LIST_CAP]]}


def apply_missing_policy(frames: Dict[str, pd.DataFrame], policy: str
                         ) -> Tuple[Dict[str, pd.DataFrame], MissingRecord]:
    """Apply ``dataset.missing_policy`` to aligned layers. Never imputes.

    ``refuse`` raises :class:`ConfigError` naming the layer, the count and
    example positions. ``drop_rows`` removes an isolate from every layer if
    any layer has no value for it, so the layers stay aligned. ``drop_columns``
    removes the affected features from the layer that holds them.
    """
    record = MissingRecord(policy=policy)
    for role, df in frames.items():
        record.layers[role] = layer_missing(df)
    if not record.any_missing:
        return frames, record
    if policy == "refuse":
        parts = []
        for role, m in record.layers.items():
            if m["cells"]:
                parts.append(
                    f"layer {role!r}: {m['cells']} empty cells in "
                    f"{m['rows_with_missing']} rows and "
                    f"{m['columns_with_missing']} columns (rows such as "
                    f"{m['example_rows']}, columns such as {m['example_columns']})")
        raise ConfigError(
            "binary layers hold cells with no value; a 0 written there would "
            "read as absence of the trait, so the run stops. " + "; ".join(parts)
            + ". Set dataset.missing_policy to 'drop_rows' or 'drop_columns' "
            "to remove them, or repair the input.")
    if policy == "drop_rows":
        bad = pd.Index([])
        for df in frames.values():
            bad = bad.union(df.index[df.isna().any(axis=1)])
        record.dropped_rows = [str(i) for i in bad]
        frames = {role: df.drop(index=bad, errors="ignore")
                  for role, df in frames.items()}
        return frames, record
    if policy == "drop_columns":
        out = {}
        for role, df in frames.items():
            cols = df.columns[df.isna().any(axis=0)]
            record.dropped_columns[role] = [str(c) for c in cols]
            out[role] = df.drop(columns=cols)
        return out, record
    raise ConfigError(f"unknown missing_policy {policy!r}")


def non_binary_values(df: pd.DataFrame) -> Dict[str, int]:
    """Values other than 0 and 1 (empty cells excluded), with their counts.

    ``True`` and ``False`` count as other values: pandas reads them from a
    CSV as booleans, ``True == 1`` would let them through, and a layer coded
    that way was not coded as the contract states.
    """
    vals = pd.Series(df.to_numpy().ravel())
    vals = vals[vals.notna()]
    is_bool = vals.map(lambda v: isinstance(v, (bool, np.bool_)))
    bad = vals[is_bool | ~vals.isin([0, 1, 0.0, 1.0, "0", "1"])]
    return {str(k): int(v) for k, v in bad.value_counts().head(_LIST_CAP).items()}


def group_adequacy(lineage: pd.Series) -> dict:
    """Per-group sizes and what the clonal-share estimator will do with them.

    The groups are coded by :func:`amr_clonalshare.attribution._codes`, the
    same routine the estimator uses, so an untyped isolate (empty, ``NaN``,
    ``"NA"``, ``"none"``) joins one ``__missing__`` level here exactly as it
    does there, and ``support`` is the estimator's own quantity: the share of
    isolates in levels of at least :data:`MIN_GROUP_SIZE` members, the missing
    level included. The count of untyped isolates is reported beside it so a
    reader can see how much of the support that level supplies.
    """
    values = np.asarray(lineage, dtype=object)
    n = int(values.size)
    codes = _codes(values) if n else np.zeros(0, dtype=int)
    counts = np.bincount(codes) if n else np.zeros(0, dtype=int)
    names: Dict[int, str] = {}
    for c, v in zip(codes, values):
        names.setdefault(int(c), _MISSING_NAME if _is_missing(v) else str(v))
    sizes = pd.Series({names[i]: int(counts[i]) for i in range(counts.size)})
    n_untyped = int(sum(_is_missing(v) for v in values))
    informative = sizes[sizes >= MIN_GROUP_SIZE]
    support = float(informative.sum() / n) if n else 0.0
    n_groups = int(sizes.size)
    n0 = ((n - float((sizes ** 2).sum()) / n) / (n_groups - 1)
          if n_groups > 1 and n else float("nan"))
    return {
        "n": n,
        "n_untyped": n_untyped,
        "n_groups": n_groups,
        "n_singletons": int((sizes == 1).sum()),
        "n_informative_groups": int(informative.size),
        "support": support,
        "support_threshold": SUPPORT_THRESHOLD,
        "estimable": bool(support >= SUPPORT_THRESHOLD),
        "effective_group_size": float(n0),
        "group_sizes": {str(k): int(v) for k, v in sizes.items()},
        "smallest_groups": {str(k): int(v)
                            for k, v in sizes.sort_values().head(_LIST_CAP).items()},
    }


_MISSING_NAME = "__missing__"


def _is_missing(v) -> bool:
    return (v is None or v != v or str(v).strip() == ""
            or str(v).lower() in ("nan", "na", "none"))


def trait_adequacy(X: pd.DataFrame) -> Dict[str, dict]:
    """Per-trait counts of the rarer outcome against :data:`MIN_MINOR_COUNT`."""
    out = {}
    for col in X.columns:
        y = X[col].to_numpy(dtype=float)
        n = int(y.size)
        ones = int(np.nansum(y))
        minor = min(ones, n - ones)
        out[str(col)] = {"n": n, "prevalence": float(ones / n) if n else 0.0,
                         "minor_count": int(minor), "constant": bool(minor == 0),
                         "adequate": bool(minor >= MIN_MINOR_COUNT)}
    return out


def input_qc(frames: Dict[str, pd.DataFrame], missing: MissingRecord, *,
             metadata: Optional[pd.DataFrame] = None,
             lineage_column: Optional[str] = None) -> dict:
    """The full input record written with every run and by ``--qc-only``."""
    ids = None
    layers = {}
    for role, df in frames.items():
        ids = df.index if ids is None else ids
        layers[role] = {"n_isolates": int(len(df)), "n_features": int(df.shape[1]),
                        "non_binary_values": non_binary_values(df),
                        "traits": trait_adequacy(df)}
    qc = {"missing": missing.as_dict(), "layers": layers,
          "n_isolates_aligned": int(len(ids)) if ids is not None else 0,
          "min_minor_count": MIN_MINOR_COUNT, "min_group_size": MIN_GROUP_SIZE}
    if metadata is not None and ids is not None:
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
    m = qc["missing"]
    lines.append(f"Isolates aligned across layers: {qc['n_isolates_aligned']}.")
    total_missing = sum(v["cells"] for v in m["layers"].values())
    if total_missing == 0:
        lines.append("No empty cells were found in any binary layer.")
    else:
        lines.append(f"Empty cells were found ({total_missing} in total) and the "
                     f"policy `{m['policy']}` was applied. Nothing was filled in: "
                     "an empty cell is a measurement that was not made, not a "
                     "trait that is absent.")
        if m["dropped_rows"]:
            lines.append(f"Isolates removed: {len(m['dropped_rows'])}.")
        for role, cols in m["dropped_columns"].items():
            if cols:
                lines.append(f"Features removed from layer `{role}`: {len(cols)}.")
    lines.append("")
    for role, L in qc["layers"].items():
        lines.append(f"## Layer `{role}`")
        lines.append("")
        lines.append(f"{L['n_isolates']} isolates, {L['n_features']} features.")
        traits = L["traits"]
        weak = [c for c, t in traits.items() if not t["adequate"] and not t["constant"]]
        const = [c for c, t in traits.items() if t["constant"]]
        ok = len(traits) - len(weak) - len(const)
        lines.append(
            f"Traits with at least {qc['min_minor_count']} isolates of the rarer "
            f"outcome: {ok} of {len(traits)}. Below that number the estimate is "
            "still computed, but its interval will be wide and the report says so.")
        if weak:
            lines.append(f"Traits below the threshold: {', '.join(weak[:_LIST_CAP])}"
                         + (" and others" if len(weak) > _LIST_CAP else "") + ".")
        if const:
            lines.append(f"Traits with a single value (nothing to explain): "
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
