"""A first configuration read off the columns of the supplied tables.

The mapping between a laboratory's column names and the settings a run needs
is the step that costs a new user the most, and it is mechanical: the names
carry the meaning. This module guesses it, marks every guess, and leaves the
correction to the reader. It writes nothing and runs no analysis.
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .phenotype import _VOCABULARIES, _MISSING, _token

__all__ = ["draft_config"]

# The first pattern that a column name contains, in this order, wins.
_ROLES: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("strain_id_column", ("isolate", "strain", "sample", "accession", "id")),
    ("lineage_column", ("lineage", "sequence_type", "mlst", "clonal", "clone",
                        "cluster", "serovar", "serotype", "st", "cc")),
    ("phenotype_antibiotic_column", ("antibiotic", "antimicrobial", "agent",
                                     "drug", "compound")),
    ("phenotype_call_column", ("resistant_phenotype", "interpretation", "sir",
                               "phenotype", "result", "call")),
)


def _is_workbook(path: Path) -> bool:
    return str(path).lower().endswith((".xlsx", ".xlsm"))


def _sheet(path: Path, limit: int = 500) -> List[dict]:
    try:
        import pandas as pd
    except ImportError as exc:  # pragma: no cover - pandas is a hard dependency
        raise ValueError(f"reading {path} needs pandas") from exc
    try:
        frame = pd.read_excel(path, dtype=str, keep_default_na=False, nrows=limit)
    except ImportError as exc:
        raise ValueError(
            f"{path} is a workbook; reading one needs openpyxl, which the core "
            "install leaves out: pip install 'amr-clonalshare[excel]'") from exc
    return [{str(k): ("" if v is None else str(v)) for k, v in row.items()}
            for row in frame.to_dict("records")]


def _columns(path: Path) -> List[str]:
    if _is_workbook(path):
        rows = _sheet(path, limit=1)
        return list(rows[0]) if rows else []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return next(csv.reader(handle), [])


def _rows(path: Path, limit: int = 500) -> List[dict]:
    if _is_workbook(path):
        return _sheet(path, limit)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        return [row for _, row in zip(range(limit), reader)]


def _match(columns: Sequence[str], patterns: Sequence[str]) -> Optional[str]:
    lowered = [(c, c.strip().lower()) for c in columns]
    for pattern in patterns:
        for original, name in lowered:
            if name == pattern:
                return original
        for original, name in lowered:
            if pattern in name:
                return original
    return None


def _kind(values: Sequence[str]) -> Optional[str]:
    """The declared vocabulary that reads every value it is shown."""
    seen = {_token(v) for v in values} - _MISSING
    if not seen:
        return None
    fits = [name for name, vocab in _VOCABULARIES.items()
            if name != "undeclared" and seen <= set(vocab)]
    return fits[0] if len(fits) == 1 else None


def draft_config(paths: Sequence[Path]) -> str:
    """A YAML configuration for these tables, with every guess marked."""
    if not paths:
        raise ValueError("name at least one CSV table")
    tables: Dict[Path, List[str]] = {}
    for path in paths:
        if not path.exists():
            raise ValueError(f"{path} does not exist")
        columns = _columns(path)
        if not columns:
            raise ValueError(f"{path} has no header row")
        tables[path] = columns

    def find(role: str) -> Tuple[Optional[Path], Optional[str]]:
        patterns = dict(_ROLES)[role]
        for path, columns in tables.items():
            hit = _match(columns, patterns)
            if hit is not None:
                return path, hit
        return None, None

    calls_path, agent = find("phenotype_antibiotic_column")
    _, call = find("phenotype_call_column")
    meta_path, lineage = find("lineage_column")
    _, ident = find("strain_id_column")
    if calls_path is None:
        calls_path = next(iter(tables))
    if meta_path is None:
        meta_path = calls_path

    kind = None
    agent_columns: List[str] = []
    if call is not None:
        kind = _kind([row.get(call, "") for row in _rows(calls_path)])
    else:
        # No column names a call, so the sheet may be wide: one column per
        # agent. A column counts as one when its own values fit a vocabulary.
        rows = _rows(calls_path)
        taken = {ident, lineage} - {None}
        fits: Dict[str, Optional[str]] = {
            column: _kind([row.get(column, "") for row in rows])
            for column in tables[calls_path] if column not in taken}
        vocabularies = {k for k in fits.values() if k}
        if len(vocabularies) == 1:
            kind = vocabularies.pop()
            agent_columns = [c for c, k in fits.items() if k == kind]

    unknown = "REPLACE_ME"
    # data_dir is absolute, so the draft runs wherever it is saved; the tables
    # are named relative to it. Every value is quoted, so a column called
    # "Isolate #", "yes" or "2019" reads back as the same string.
    directory = calls_path.parent.resolve()
    relative = (lambda p: Path(os.path.relpath(p.resolve(), directory)).as_posix())
    q = json.dumps
    lines = [
        "# Drafted by amr-clonalshare --init from the column names of:",
        *[f"#   {Path(p).as_posix()}" for p in tables],
        "# Every value below is a guess. Read it before running the analysis;",
        "# a line marked REPLACE_ME could not be guessed at all.",
        "dataset:",
        "  name: my_collection",
        f"  data_dir: {q(directory.as_posix())}",
        f"  metadata: {q(relative(meta_path))}",
        f"  strain_id_column: {q(ident or unknown)}",
        f"  lineage_column: {q(lineage or unknown)}",
        f"  phenotype: {q(relative(calls_path))}",
        f"  phenotype_id_column: {q(ident or unknown)}",
    ]
    if agent_columns:
        lines.append("  # one column per agent; the long shape names the two "
                     "columns instead")
        lines.append("  phenotype_agent_columns: [%s]" % ", ".join(q(c) for c in agent_columns))
    else:
        lines.append(f"  phenotype_antibiotic_column: {q(agent or unknown)}")
        lines.append(f"  phenotype_call_column: {q(call or unknown)}")
    if kind:
        lines.append(f"  phenotype_kind: {kind}")
    else:
        lines.append(f"  phenotype_kind: {unknown}"
                     "   # clinical_sir, wt_nwt or binary; the recorded values fit none of them")
    return "\n".join(lines) + "\n"
