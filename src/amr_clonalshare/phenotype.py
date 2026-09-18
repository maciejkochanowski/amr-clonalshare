"""Explicit categorical phenotype semantics with a lossless raw-record inventory."""
from __future__ import annotations

import warnings
import numpy as np
import pandas as pd

__all__ = ["to_non_susceptible"]

_MISSING = {"", "nan", "none", "na", "n/a", "<na>"}
_SIR = {"s": "s", "susceptible": "s", "i": "i", "intermediate": "i",
        "susceptible, increased exposure": "i", "susceptible increased exposure": "i",
        "r": "r", "resistant": "r"}
_LEGACY_SIR = {"resistant": "r", "intermediate": "i", "susceptible": "s",
               "nonsusceptible": "r", "non-susceptible": "r"}


def _token(value):
    return "" if pd.isna(value) else str(value).strip().lower()


def resolve_duplicates(df, keys, values, policy="error", *, what="table", legacy_max=None):
    """Collapse equal records; conflicts are explicit, including missing vs known."""
    if policy not in ("error", "drop_conflicts", "legacy"):
        raise ValueError(f"unknown duplicate policy {policy!r}")
    duplicates = df.duplicated(keys, keep=False)
    conflicts, conflict_rows = [], 0
    for key, group in df.loc[duplicates].groupby(keys, sort=False, dropna=False):
        if len(group[values].drop_duplicates()) > 1:
            conflicts.append(key if isinstance(key, tuple) else (key,))
            conflict_rows += len(group)
    n_identical = len(df) - len(df.drop_duplicates(keys + values))
    report = {"duplicate_policy": policy, "n_duplicate_rows": int(df.duplicated(keys).sum()),
              "n_identical_duplicates_collapsed": n_identical, "n_conflicting_keys": len(conflicts),
              "n_conflicting_rows": conflict_rows,
              "conflicting_key_examples": [list(k) for k in conflicts[:5]],
              "n_conflict_rows_dropped": conflict_rows if policy == "drop_conflicts" else 0}
    if conflicts and policy == "error":
        raise ValueError(f"{what} has conflicting records for {len(conflicts)} key(s): {conflicts[:5]}; correct the source or explicitly select drop_conflicts or legacy")
    if policy == "drop_conflicts" and conflicts:
        conflict_set = set(conflicts)
        df = df.loc[[tuple(row) not in conflict_set for row in df[keys].itertuples(index=False, name=None)]]
    if policy == "legacy":
        warnings.warn(f"{what}: explicit legacy duplicate policy; conflicting call positives win and other tables keep the first record", RuntimeWarning, stacklevel=2)
        if legacy_max is not None:
            return df.groupby(keys, sort=False, as_index=False, dropna=False)[legacy_max].max(), report
    return df.drop_duplicates(keys).copy(), report


def to_non_susceptible(long_df: pd.DataFrame, *, id_column="Strain_ID",
                       antibiotic_column="antibiotic", call_column="resistant_phenotype",
                       intermediate=None, kind="legacy", duplicate_policy="error",
                       positive_definition=None, source=None, ast_standard=None, ast_version=None):
    """Return a 0/1/NaN panel preserving all raw isolate and agent keys.

    Clinical defaults to R-only; WT/NWT and binary are separate vocabularies.
    The legacy name remains an API alias; it never makes NWT a clinical call.
    ``positive_definition`` supplies a label, not a different recoding rule.
    """
    if kind not in ("legacy", "clinical_sir", "wt_nwt", "binary"):
        raise ValueError(f"unknown phenotype kind {kind!r}")
    if intermediate not in (None, "non_susceptible", "drop", "susceptible"):
        raise ValueError(f"unknown intermediate policy {intermediate!r}")
    if kind in ("wt_nwt", "binary") and intermediate is not None:
        raise ValueError("intermediate policy applies only to clinical_sir or legacy")
    if kind == "legacy":
        warnings.warn("legacy phenotype vocabulary: declare clinical_sir, wt_nwt or binary for explicit interpretation", RuntimeWarning, stacklevel=2)
    policy = intermediate or ("non_susceptible" if kind == "legacy" else "susceptible")
    df = long_df[[id_column, antibiotic_column, call_column]].copy()
    for key in (id_column, antibiotic_column):
        if df[key].isna().any() or df[key].astype(str).str.strip().eq("").any():
            raise ValueError(f"missing key in {key}")
    ids = pd.Index(df[id_column].unique(), name=id_column)
    agents = pd.Index(df[antibiotic_column].unique(), name=antibiotic_column)
    token = df[call_column].map(_token)
    missing = token.isin(_MISSING)
    if kind in ("legacy", "clinical_sir"):
        vocab = _LEGACY_SIR if kind == "legacy" else _SIR
        canonical = token.map(vocab)
        mapped = canonical.map({"s": 0., "r": 1., "i": {"non_susceptible": 1., "susceptible": 0., "drop": np.nan}[policy]})
        definition = "R or I" if policy == "non_susceptible" else "R (resistant)"
        dropped = canonical.eq("i") & (policy == "drop")
    else:
        vocab = ({"wt": "0", "wild type": "0", "wild-type": "0", "nwt": "1", "non-wild-type": "1", "non wild type": "1"}
                 if kind == "wt_nwt" else {"0": "0", "0.0": "0", "1": "1", "1.0": "1"})
        canonical = token.map(vocab)
        mapped = canonical.map({"0": 0., "1": 1.})
        definition = "NWT (non-wild-type)" if kind == "wt_nwt" else "1 (positive)"
        dropped = pd.Series(False, index=df.index)
    unread = ~missing & canonical.isna()
    # Compare semantic calls before dichotomization: I and R remain distinct records.
    df['_canonical'] = canonical.where(canonical.notna(), 'unknown:' + token).where(~missing, 'missing')
    df['_ns'] = mapped
    registry = {}
    for agent in agents:
        sel = df[antibiotic_column].eq(agent)
        registry[str(agent)] = {"rows_read": int(sel.sum()), "n_recorded_isolates": int(df.loc[sel, id_column].nunique()),
                                "n_missing_calls": int((sel & missing).sum()), "n_unrecognized_calls": int((sel & unread).sum()),
                                "n_intermediate_dropped": int((sel & dropped).sum())}
    resolved, counts = resolve_duplicates(df, [id_column, antibiotic_column], ['_canonical'], duplicate_policy,
                                          what='phenotype table', legacy_max='_ns')
    wide = resolved.pivot(index=id_column, columns=antibiotic_column, values='_ns').reindex(index=ids, columns=agents).astype(float)
    for agent in agents:
        registry[str(agent)]['n_readable_isolates'] = int(wide[agent].notna().sum())
    wide.attrs['phenotype_reading'] = {"rows_read": len(df), "n_unrecognized_calls": int(unread.sum()),
        "unrecognized_calls": {str(k): int(v) for k, v in token[unread].value_counts().items()},
        "n_missing_calls": int(missing.sum()), "n_intermediate_dropped": int(dropped.sum()),
        "isolates_without_readable_calls": int(wide.isna().all(axis=1).sum()),
        "phenotype_kind": kind, "positive_definition": positive_definition or definition,
        "positive_coding": definition, "intermediate_policy": policy if kind in ('legacy', 'clinical_sir') else None,
        "source": source, "ast_standard": ast_standard, "ast_version": ast_version,
        "agents": registry, **counts}
    return wide
