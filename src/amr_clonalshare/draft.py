"""A first configuration read off the columns of the supplied tables.

The mapping between a laboratory's column names and the settings a run needs
is the step that costs a new user the most, and it is mechanical: the names
carry the meaning. This module guesses it, marks every guess, and leaves the
correction to the reader. It writes nothing and runs no analysis.
"""
from __future__ import annotations

import csv
import json
import math
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .phenotype import _VOCABULARIES, _MISSING, _token

__all__ = ["draft_config"]

# The first pattern that a column name contains, in this order, wins.
_ROLES: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("strain_id_column", ("isolate", "strain", "sample", "accession", "id")),
    ("lineage_column", ("lineage", "sequence_type", "st", "mlst", "clonal", "clone",
                        "cluster", "serovar", "serotype", "cc", "phylogroup")),
    ("phenotype_antibiotic_column", ("antibiotic", "antimicrobial", "agent",
                                     "drug", "compound")),
    ("phenotype_call_column", ("resistant_phenotype", "interpretation", "sir",
                               "phenotype", "result", "call")),
)


def _is_workbook(path: Path) -> bool:
    from .io import is_workbook
    return is_workbook(path)


def _sheet(path: Path, limit: int = 500) -> List[dict]:
    """The first rows of the one sheet of a workbook, read as the run reads
    it, so that a workbook the run would refuse is refused here too."""
    from .io import read_workbook
    frame = read_workbook(path, "table", dtype=str, nrows=limit)
    return [{str(k): ("" if v is None else str(v)) for k, v in row.items()}
            for row in frame.to_dict("records")]


def _columns(path: Path) -> List[str]:
    if _is_workbook(path):
        rows = _sheet(path, limit=1)
        return [str(c).strip() for c in rows[0]] if rows else []
    from .io import text_separator
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return [c.strip() for c in next(csv.reader(handle, delimiter=text_separator(path)), [])]


def _rows(path: Path, limit: int = 500) -> List[dict]:
    """The first rows of a table as the run reads them: the separator of the
    file, column names without surrounding spaces, no empty rows."""
    if _is_workbook(path):
        rows = _sheet(path, limit)
    else:
        from .io import text_separator
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle, delimiter=text_separator(path))
            rows = [{(k or "").strip(): (v or "") for k, v in row.items()}
                    for _, row in zip(range(limit), reader)]
    return [row for row in rows if any(str(v).strip() for v in row.values())]


def _match(columns: Sequence[str], patterns: Sequence[str]) -> Optional[str]:
    """The first column named exactly as one of the patterns, in their order;
    failing that, the first column whose name contains one, where a pattern of
    one or two letters must be a whole word of the name (``genome_id`` holds
    ``id``; ``isolate_id`` does not hold ``st``). A column called ``st`` is
    therefore preferred to ``serotype`` in a table that holds both."""
    lowered = [(c, c.strip().lower()) for c in columns]
    for pattern in patterns:
        for original, name in lowered:
            if name == pattern:
                return original
    for pattern in patterns:
        for original, name in lowered:
            words = re.split(r"[^a-z0-9]+", name)
            if (pattern in name) if len(pattern) >= 3 else (pattern in words):
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


_SIGN = re.compile(r"^(<=|>=|<|>|\u2264|\u2265)\s*")
_MIC_VALUE_PATTERNS = ("measurement", "mic", "value", "concentration", "dilution")
_MIC_OPERATOR_PATTERNS = ("operator", "sign", "censor", "modifier")
_PANEL_PATTERNS = ("laboratory", "testing_lab", "lab", "panel", "plate")


def _concentration(value: str) -> Optional[float]:
    """The number of a recorded concentration without its sign, or None."""
    text = _SIGN.sub("", _token(value))
    text = re.sub(r"^(\d+),(\d+)$", r"\1.\2", text)
    text = re.sub(r"^([0-9.eE+-]+)\s*/\s*[0-9.eE+-]+$", r"\1", text)
    try:
        number = float(text)
    except ValueError:
        return None
    return number if math.isfinite(number) and number > 0 else None


def _mic_like(values: Sequence[str]) -> bool:
    """Whether a column holds recorded concentrations: nearly every value a
    positive number on the twofold dilution series, and at least one of them
    censored, below 1 or one of four distinct steps, so that a count or a
    year column is not taken for one."""
    seen = [v for v in values if _token(v) not in _MISSING]
    if len(seen) < 4:
        return False
    numbers = [_concentration(v) for v in seen]
    good = [n for n in numbers if n is not None]
    if len(good) < 0.9 * len(seen):
        return False
    # 0.015, 0.03, 0.06 and 0.12 are 1/64 to 1/8 mg/L as a plate sheet writes
    # them: a value within a quarter of a doubling of a power of two is on the series
    on_series = [abs(math.log2(n) - round(math.log2(n))) <= 0.25 for n in good]
    if sum(on_series) < 0.9 * len(good):
        return False
    # a column of years sits on one step of the series; a panel spans several
    if len({round(math.log2(n)) for n in good}) < 2:
        return False
    signed = any(_SIGN.match(_token(v)) for v in seen)
    return signed or min(good) < 1 or len(set(good)) >= 4


def observed_wells(values: Sequence[str]) -> List[float]:
    """The distinct concentrations of a column, in order."""
    return sorted({n for n in (_concentration(v) for v in values) if n is not None})


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
    rows: Dict[Path, List[dict]] = {path: _rows(path, limit=10 ** 7) for path in tables}

    def find(role: str, patterns: Optional[Sequence[str]] = None,
             among: Optional[Sequence[Path]] = None) -> Tuple[Optional[Path], Optional[str]]:
        patterns = patterns or dict(_ROLES)[role]
        for path in (among if among is not None else tables):
            hit = _match(tables[path], patterns)
            if hit is not None:
                return path, hit
        return None, None

    meta_path, lineage = find("lineage_column")
    _, ident = find("strain_id_column")
    taken = {ident, lineage} - {None}

    # The MIC table: a column of concentrations beside an agent column (the
    # long shape), or several concentration columns (one per agent).
    mic_path: Optional[Path] = None
    mic_value = mic_agent = mic_operator = mic_panel = mic_unit = None
    mic_agents: List[str] = []
    for path, columns in tables.items():
        fits = [c for c in columns if c not in taken and _mic_like([r.get(c, "") for r in rows[path]])]
        agent = _match([c for c in columns if c not in fits], dict(_ROLES)["phenotype_antibiotic_column"])
        if agent and len(fits) == 1:
            mic_path, mic_agent, mic_value = path, agent, fits[0]
            break
        if fits and not agent:
            mic_path, mic_agents = path, fits
            break
    if mic_path is not None:
        mic_operator = _match([c for c in tables[mic_path] if c not in (mic_value, mic_agent)], _MIC_OPERATOR_PATTERNS)
        mic_panel = _match([c for c in tables[mic_path] if c not in (mic_value, mic_agent, mic_operator)], _PANEL_PATTERNS)
        mic_unit = _match([c for c in tables[mic_path] if c not in (mic_value, mic_agent, mic_operator, mic_panel)], ("unit",))
        if mic_agents:
            mic_operator = None

    # The call table: any other table with an agent column, or with columns
    # whose values fit one vocabulary.
    others = [p for p in tables if p != mic_path]
    calls_path, agent = find("phenotype_antibiotic_column", among=others)
    _, call = find("phenotype_call_column", among=others)
    kind = None
    agent_columns: List[str] = []
    if calls_path is not None and call is not None:
        kind = _kind([row.get(call, "") for row in rows[calls_path]])
    elif calls_path is not None:
        call = None
    if calls_path is None:
        # No column names an agent, so a sheet may be wide: one column per
        # agent. A column counts as one when its own values fit a vocabulary.
        for path in others:
            kinds = {column: _kind([row.get(column, "") for row in rows[path]])
                     for column in tables[path] if column not in taken}
            vocabularies = {k for k in kinds.values() if k}
            if len(vocabularies) == 1:
                kind = vocabularies.pop()
                agent_columns = [c for c, k in kinds.items() if k == kind]
                calls_path = path
                break
    anchor = calls_path or mic_path or next(iter(tables))
    if calls_path is None and mic_path is None:
        calls_path = anchor
    if meta_path is None:
        meta_path = calls_path or mic_path

    unknown = "REPLACE_ME"
    # data_dir is absolute, so the draft runs wherever it is saved; the tables
    # are named relative to it. Every value is quoted, so a column called
    # "Isolate #", "yes" or "2019" reads back as the same string.
    directory = anchor.parent.resolve()
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
    ]
    if lineage and meta_path in (calls_path, mic_path):
        lines.insert(-3, "  # One sheet carries the lineage beside the readings, so it is the metadata too.")
    if calls_path is not None:
        lines += [
            "",
            "  # Susceptibility calls: S/I/R, wild-type/non-wild-type or 0/1.",
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
    if mic_path is not None:
        lines += [
            "",
            "  # Recorded MICs: the lineage share of the recorded MIC ordering and the",
            "  # bounds on the latent MIC ordering are read from this table.",
            f"  mic: {q(relative(mic_path))}",
            f"  mic_id_column: {q(ident or unknown)}",
        ]
        if mic_agents:
            lines.append("  # one column per agent, each cell a concentration with its sign (<=0.5, 8, >16)")
            lines.append("  mic_agent_columns: [%s]" % ", ".join(q(c) for c in mic_agents))
            wells = {c: observed_wells([r.get(c, "") for r in rows[mic_path]]) for c in mic_agents}
        else:
            lines.append(f"  mic_antibiotic_column: {q(mic_agent)}")
            lines.append(f"  mic_value_column: {q(mic_value)}")
            if mic_operator:
                lines.append(f"  mic_operator_column: {q(mic_operator)}")
            readings: Dict[str, List[str]] = {}
            for row in rows[mic_path]:
                readings.setdefault(str(row.get(mic_agent, "")).strip(), []).append(row.get(mic_value, ""))
            wells = {a: observed_wells(v) for a, v in readings.items() if a}
        if mic_unit:
            lines.append(f"  mic_unit_column: {q(mic_unit)}   # rows in a unit other than mg/L are set aside")
        if mic_panel:
            lines.append(f"  mic_panel_column: {q(mic_panel)}   # laboratory or panel; readings are ordered within its levels")
        lines += [
            "  # The tested concentrations of every agent. The values below are the",
            "  # concentrations OBSERVED in the table, not the tested",
            "  # range: replace them with the dilution range of the panel, lowest to",
            "  # highest well, so that a reading on an end well is read as the interval",
            "  # it is. Delete the block to take the range from the readings themselves.",
            "  mic_wells:",
        ]
        for agent_name, values in wells.items():
            shown = ", ".join(str(int(v)) if float(v).is_integer() else f"{v:g}" for v in values)
            lines.append(f"    {q(agent_name)}: [{shown}]")
    lines += [
        "",
        "  # Optional settings, each a column of the metadata table; remove the",
        "  # leading '# ' to use one.",
        "  # phenotype_covariate_column: \"country\"   # strata of the calls: lineage labels are compared within its levels",
        "  # mic_covariate_columns: [\"country\"]      # strata of the MIC readings, crossed with the panel",
        "  # unit_column: \"farm\"                      # sampling unit: labels are exchanged within units",
        "  # stratify_by: \"period\"                    # repeat every analysis within each level",
        "  # contrast_column: \"period\"                # two collections whose prevalence difference is decomposed",
        "  # contrast_levels: [\"earlier\", \"later\"]",
        "",
        "# censored:",
        "#   end_wells_censored: true   # a reading on an end well is an interval (the default)",
    ]
    return "\n".join(lines) + "\n"
