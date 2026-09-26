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
_UNDECLARED_SIR = {"resistant": "r", "intermediate": "i", "susceptible": "s",
               "nonsusceptible": "r", "non-susceptible": "r"}
_WT_NWT = {"wt": "0", "wild type": "0", "wild-type": "0", "nwt": "1",
           "non-wild-type": "1", "non wild type": "1"}
_BINARY = {"0": "0", "0.0": "0", "1": "1", "1.0": "1"}
_VOCABULARIES = {"clinical_sir": _SIR, "undeclared": _UNDECLARED_SIR,
                 "wt_nwt": _WT_NWT, "binary": _BINARY}


def _kinds_reading(tokens, declared):
    """The other declared kinds whose vocabulary reads every unreadable token."""
    return [name for name, vocab in _VOCABULARIES.items()
            if name != declared and all(t in vocab for t in tokens)]


def _token(value):
    return "" if pd.isna(value) else str(value).strip().lower()


def resolve_duplicates(df, keys, values, policy="error", *, what="table", positive_max=None):
    """Collapse equal records; conflicts are explicit, including missing vs known."""
    if policy not in ("error", "drop_conflicts", "positive_wins"):
        raise ValueError(f"unknown duplicate policy {policy!r}")
    duplicates = df.duplicated(keys, keep=False)
    conflicts, conflict_rows, varying = [], 0, set()
    for key, group in df.loc[duplicates].groupby(keys, sort=False, dropna=False):
        if len(group[values].drop_duplicates()) > 1:
            conflicts.append(key if isinstance(key, tuple) else (key,))
            conflict_rows += len(group)
            varying.update(c for c in values if group[c].nunique(dropna=False) > 1)
    n_identical = len(df) - len(df.drop_duplicates(keys + values))
    report = {"duplicate_policy": policy, "n_duplicate_rows": int(df.duplicated(keys).sum()),
              "n_identical_duplicates_collapsed": n_identical, "n_conflicting_keys": len(conflicts),
              "n_conflicting_rows": conflict_rows,
              "conflicting_key_examples": [list(k) for k in conflicts[:5]],
              "n_conflict_rows_dropped": conflict_rows if policy == "drop_conflicts" else 0}
    if conflicts and policy == "error":
        where = ", ".join(repr(c) for c in sorted(varying))
        raise ValueError(f"{what} has conflicting records for {len(conflicts)} key(s) in column(s) {where}: {conflicts[:5]}; correct the source or explicitly select drop_conflicts or positive_wins")
    if policy == "drop_conflicts" and conflicts:
        conflict_set = set(conflicts)
        df = df.loc[[tuple(row) not in conflict_set for row in df[keys].itertuples(index=False, name=None)]]
    if policy == "positive_wins":
        warnings.warn(f"{what}: explicit positive_wins duplicate policy; conflicting call positives win and other tables keep the first record", RuntimeWarning, stacklevel=2)
        if positive_max is not None:
            return df.groupby(keys, sort=False, as_index=False, dropna=False)[positive_max].max(), report
    return df.drop_duplicates(keys).copy(), report


def to_non_susceptible(long_df: pd.DataFrame, *, id_column="Strain_ID",
                       antibiotic_column="antibiotic", call_column="resistant_phenotype",
                       intermediate=None, kind="undeclared", duplicate_policy="error",
                       positive_definition=None, source=None, ast_standard=None, ast_version=None):
    """Return a 0/1/NaN panel preserving all raw isolate and agent keys.

    Clinical defaults to R-only; WT/NWT and binary are separate vocabularies.
    The undeclared kind reads resistance words only; it never makes NWT a
    clinical call.
    ``positive_definition`` supplies a label, not a different recoding rule.
    """
    if kind not in ("undeclared", "clinical_sir", "wt_nwt", "binary"):
        raise ValueError(f"unknown phenotype kind {kind!r}")
    if intermediate not in (None, "non_susceptible", "drop", "susceptible"):
        raise ValueError(f"unknown intermediate policy {intermediate!r}")
    if kind in ("wt_nwt", "binary") and intermediate is not None:
        raise ValueError("intermediate policy applies only to clinical_sir or undeclared")
    if kind == "undeclared":
        warnings.warn("undeclared phenotype vocabulary: declare clinical_sir, wt_nwt or binary for explicit interpretation", RuntimeWarning, stacklevel=2)
    policy = intermediate or ("non_susceptible" if kind == "undeclared" else "susceptible")
    df = long_df[[id_column, antibiotic_column, call_column]].copy()
    for key in (id_column, antibiotic_column):
        if df[key].isna().any() or df[key].astype(str).str.strip().eq("").any():
            raise ValueError(f"missing key in {key}")
    # The panel is a function of the records, not of the order in which the
    # rows were written: a table sorted differently is the same table, and one
    # analysis of it must give the same numbers at the same seed. Identifiers
    # and agent names are ordered as strings, as the matched-record comparison
    # already orders its records.
    ids = pd.Index(sorted(df[id_column].unique().tolist(), key=str), name=id_column)
    folded: dict[str, list[str]] = {}
    for name in df[antibiotic_column].dropna().unique().tolist():
        folded.setdefault(str(name).strip().lower(), []).append(str(name))
    clashes = [sorted(v) for v in folded.values() if len(v) > 1]
    if clashes:
        warnings.warn(f"agent names that differ only in case or spaces are read as different agents: "
                      f"{clashes[:3]}; write each agent one way", RuntimeWarning, stacklevel=2)
    agents = pd.Index(sorted(df[antibiotic_column].unique().tolist(), key=str),
                      name=antibiotic_column)
    token = df[call_column].map(_token)
    missing = token.isin(_MISSING)
    if kind in ("undeclared", "clinical_sir"):
        vocab = _UNDECLARED_SIR if kind == "undeclared" else _SIR
        canonical = token.map(vocab)
        mapped = canonical.map({"s": 0., "r": 1., "i": {"non_susceptible": 1., "susceptible": 0., "drop": np.nan}[policy]})
        definition = "R or I" if policy == "non_susceptible" else "R (resistant)"
        dropped = canonical.eq("i") & (policy == "drop")
    else:
        vocab = _WT_NWT if kind == "wt_nwt" else _BINARY
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
                                          what='phenotype table', positive_max='_ns')
    wide = resolved.pivot(index=id_column, columns=antibiotic_column, values='_ns').reindex(index=ids, columns=agents).astype(float)
    for agent in agents:
        registry[str(agent)]['n_readable_isolates'] = int(wide[agent].notna().sum())
    reading = {"rows_read": len(df), "n_unrecognized_calls": int(unread.sum()),
        "unrecognized_calls": {str(k): int(v) for k, v in token[unread].value_counts().items()},
        "n_missing_calls": int(missing.sum()), "n_intermediate_dropped": int(dropped.sum()),
        "isolates_without_readable_calls": int(wide.isna().all(axis=1).sum()),
        "phenotype_kind": kind, "positive_definition": positive_definition or definition,
        "positive_coding": definition, "intermediate_policy": policy if kind in ('undeclared', 'clinical_sir') else None,
        "source": source, "ast_standard": ast_standard, "ast_version": ast_version,
        "agents": registry, **counts}
    if reading["n_unrecognized_calls"] and not wide.notna().to_numpy().any():
        # Nothing was read at all: the recorded vocabulary does not match the
        # declared kind. That is an input error, not an empty estimate.
        seen = list(reading["unrecognized_calls"].items())
        shown = ", ".join(f"{value!r} ({n})" for value, n in seen[:6])
        if len(seen) > 6:
            shown += f" and {len(seen) - 6} more"
        fits = _kinds_reading(reading["unrecognized_calls"], kind)
        raise ValueError(
            f"no call in the phenotype table could be read under phenotype_kind "
            f"{kind!r}: {reading['n_unrecognized_calls']} unreadable call(s), {shown}; "
            + (f"these values are read by phenotype_kind: {' or '.join(fits)}"
               if fits else
               "declare the phenotype_kind that matches the recorded vocabulary"))
    wide.attrs['phenotype_reading'] = reading
    return wide
