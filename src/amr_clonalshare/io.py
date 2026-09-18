"""
io.py — config-driven data loading for amr-clonalshare.

The run takes a susceptibility table and a lineage column as input and does
not read sequence data: a phenotype per isolate, as a susceptibility call or
a recorded dilution, and a metadata table that holds the lineage label. The raw
cohort is the union of call-table and dilution-table identifiers; metadata is joined on the
identifier as written, and the join is measured and reported rather than
repaired.
"""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd

from .config import Config, ConfigError
from .phenotype import to_non_susceptible, resolve_duplicates, _token, _MISSING
from .qc import input_qc


# Historical CSV missing markers, normalized only in configured metadata
# analysis columns. IDs and unrelated metadata retain their exact spelling.
_METADATA_MISSING = frozenset({"", "#n/a", "#n/a n/a", "#na", "-1.#ind",
    "-1.#qnan", "-nan", "1.#ind", "1.#qnan", "<na>", "n/a", "na", "nan", "none", "null"})


#: Month abbreviations a spreadsheet writes when it reads a value as a date.
_MONTHS = ("jan", "feb", "mar", "apr", "may", "jun",
           "jul", "aug", "sep", "oct", "nov", "dec")
#: ``1-Feb`` and ``Feb-1`` only. The all-numeric form a spreadsheet also
#: produces, ``01/02/2021``, is deliberately not matched: it is the shape of
#: real strain identifiers, and the shipped cohort contains ``89-88-422``. A
#: check that fires on a strain name is a check the next analyst switches off,
#: so this one keeps to the form that has no other reading.
_DATE_LIKE = re.compile(
    r"^\s*(?:\d{1,2}[-/](?:" + "|".join(_MONTHS) + r")"
    r"|(?:" + "|".join(_MONTHS) + r")[-/]\d{1,2})\s*$",
    re.IGNORECASE)


def date_coerced_values(series: "pd.Series") -> Dict[str, int]:
    """Values in a categorical column that a spreadsheet has turned into dates.

    A serotype written ``1/2`` becomes ``1-Feb``; a gene called ``SEPT9``
    becomes ``9-Sep``. The error is old and well documented - Ziemann, Eren and
    El-Osta (2016), *Genome Biology* 17:177, found it in a fifth of published
    genomics supplements - and it is still reaching public surveillance
    databases. The shipped *Streptococcus suis* metadata carries 39 isolates
    whose serotype 1/2 arrives from BV-BRC as the string ``1-Feb``, and the
    published supplement they both derive from stores the same cells as the
    Excel serial 44228.

    Reported rather than raised. Metadata belongs to whoever supplied it and is
    often not the analyst's to repair, but a silently mistyped serotype is a
    stratum that will not group with itself, so the run says so.
    """
    values = series.dropna().astype(str)
    hits: Dict[str, int] = {}
    for value, count in values.value_counts().items():
        if _DATE_LIKE.match(value):
            hits[value] = int(count)
    return hits


def metadata_quality(metadata: "pd.DataFrame") -> Dict[str, Any]:
    """Data-quality notes on the metadata, emitted with the run."""
    if metadata is None or not len(metadata.columns):
        return {}
    coerced = {}
    for column in metadata.columns:
        found = date_coerced_values(metadata[column])
        if found:
            coerced[str(column)] = found
    notes: Dict[str, Any] = {}
    if coerced:
        notes["date_coerced_values"] = coerced
        notes["date_coerced_note"] = (
            "these values look like dates a spreadsheet produced from "
            "something else, most often a value written as a fraction such as "
            "a serotype 1/2. Check them against the source before reading any "
            "stratum built from this column")
    return notes


@dataclass
class Dataset:
    """The loaded cohort: its isolates, panel, metadata and dilutions."""
    cfg: Config
    strain_ids: Optional[pd.Index] = None
    #: Isolate-by-antimicrobial 0/1 matrix read from the phenotype table, or
    #: an empty frame on the isolates of the dilution table when no phenotype
    #: table is named.
    panel: Optional[pd.DataFrame] = None
    metadata: Optional[pd.DataFrame] = None
    #: Long-format minimum inhibitory concentrations, one row per (isolate,
    #: antimicrobial), restricted to the aligned strains. ``None`` unless
    #: ``dataset.mic`` names a table.
    mic: Optional[pd.DataFrame] = None
    #: Share of the raw cohort carrying at least one raw MIC row, and the
    #: identifiers that did not join. Both are emitted with the run: a join
    #: that quietly matched half the cohort is the failure mode this reports.
    mic_join: Optional[Dict[str, Any]] = None
    #: The input-check record from :func:`amr_clonalshare.qc.input_qc`:
    #: per-antimicrobial adequacy, the metadata join, and per-lineage group
    #: sizes with the estimator's verdict.
    input_qc: Optional[Dict[str, Any]] = None

    @property
    def n(self) -> int:
        return 0 if self.strain_ids is None else len(self.strain_ids)


def _read_table(path, what: str, **kwargs) -> pd.DataFrame:
    """``pd.read_csv`` whose every failure is a configuration error.

    A missing file, a table that is not CSV, an empty file and a file in an
    encoding other than UTF-8 are all properties of the input the user
    supplied, so each is reported as such and the run exits with the
    configuration code rather than a traceback.
    """
    try:
        return pd.read_csv(path, keep_default_na=False, **kwargs)
    except UnicodeDecodeError as exc:
        raise ConfigError(
            f"{what} {path} is not UTF-8 text (byte {exc.start}: {exc.reason}); "
            f"save it as UTF-8 CSV") from exc
    except (FileNotFoundError, pd.errors.ParserError,
            pd.errors.EmptyDataError, OSError) as exc:
        raise ConfigError(f"{what} {path}: {exc}") from exc


def _require_keys(table: pd.DataFrame, columns, what: str) -> None:
    """Identifiers and agent names cannot be missing from a joinable record."""
    for column in columns:
        missing = table[column].isna() | table[column].astype(str).str.strip().eq("")
        if missing.any():
            raise ConfigError(f"{what} has {int(missing.sum())} missing value(s) in key column {column!r}")


def read_panel(cfg: Config) -> pd.DataFrame:
    """The susceptibility panel as an isolate-by-antimicrobial 0/1 matrix."""
    fp = cfg.phenotype_path
    if fp is None:
        raise ConfigError("dataset.phenotype is not configured")
    df = _read_table(fp, "phenotype table",
                     dtype=str)
    for column in (cfg.dataset.phenotype_id_column,
                   cfg.dataset.phenotype_antibiotic_column,
                   cfg.dataset.phenotype_call_column):
        if column not in df.columns:
            raise ConfigError(
                f"phenotype table {fp} has no column {column!r}; it has "
                f"{list(df.columns)[:10]}")
    _require_keys(df, (cfg.dataset.phenotype_id_column,
                       cfg.dataset.phenotype_antibiotic_column), f"phenotype table {fp}")
    try:
        panel = to_non_susceptible(
            df, id_column=cfg.dataset.phenotype_id_column,
            antibiotic_column=cfg.dataset.phenotype_antibiotic_column,
            call_column=cfg.dataset.phenotype_call_column,
            intermediate=cfg.dataset.phenotype_intermediate,
            kind=cfg.dataset.phenotype_kind, duplicate_policy=cfg.dataset.duplicate_policy,
            positive_definition=cfg.dataset.phenotype_positive_definition,
            source=cfg.dataset.phenotype_source, ast_standard=cfg.dataset.ast_standard,
            ast_version=cfg.dataset.ast_version)
    except ValueError as exc:
        raise ConfigError(str(exc)) from exc
    panel.index.name = cfg.dataset.strain_id_column
    reading = panel.attrs.get("phenotype_reading", {})
    if reading.get("n_unrecognized_calls"):
        warnings.warn(
            f"phenotype table {fp} has {reading['n_unrecognized_calls']} unrecognized "
            f"call(s), set aside: {reading['unrecognized_calls']}. "
            "See phenotype_reading in the input check.", RuntimeWarning, stacklevel=2)
    return panel


def _dilution_index(cfg: Config) -> pd.Index:
    """The isolates of the dilution table, for a run with no call table."""
    fp = cfg.mic_path
    idc = cfg.dataset.mic_id_column
    df = _read_table(fp, "MIC table", dtype=str)
    if idc not in df.columns:
        raise ConfigError(f"MIC table {fp} has no column {idc!r}")
    _require_keys(df, (idc,), f"MIC table {fp}")
    ids = pd.Index(df[idc].astype(str).unique(),
                   name=cfg.dataset.strain_id_column)
    return ids


def load_dataset(cfg: Config) -> Dataset:
    """Load the union of raw call/MIC IDs; each method selects its own subset."""
    cfg.dataset.validate()
    strain_col = cfg.dataset.strain_id_column
    panel = read_panel(cfg) if cfg.dataset.phenotype else pd.DataFrame(index=pd.Index([], name=strain_col))
    call_ids = panel.index
    mic_ids = _dilution_index(cfg) if cfg.dataset.mic else pd.Index([])
    core = call_ids.union(mic_ids, sort=False).rename(strain_col)
    panel = panel.reindex(core)
    metadata = None
    metadata_join = None
    if cfg.metadata_path is not None:
        mp = cfg.metadata_path
        metadata = _read_table(mp, "metadata", dtype=str)
        for key, column in (("strain_id_column", strain_col), ("lineage_column", cfg.dataset.lineage_column),
                            ("contrast_column", cfg.dataset.contrast_column), ("batch_column", cfg.dataset.batch_column),
                            ("stratify_by", cfg.dataset.stratify_by)):
            if column and column not in metadata.columns:
                raise ConfigError(f"metadata {mp} has no column {column!r} (dataset.{key})")
        _require_keys(metadata, (strain_col,), f"metadata {mp}")
        analysis_columns = ({cfg.dataset.lineage_column, cfg.dataset.contrast_column,
                             cfg.dataset.batch_column, cfg.dataset.stratify_by,
                             cfg.dataset.mic_panel_column, *cfg.dataset.mic_covariate_columns}
                            & set(metadata.columns)) - {None, strain_col}
        for column in analysis_columns:
            missing = metadata[column].map(_token).isin(_METADATA_MISSING)
            metadata[column] = metadata[column].mask(missing, np.nan)
        raw_metadata_ids = pd.Index(metadata[strain_col].unique())
        n_rows = len(metadata)
        try:
            metadata, duplicates = resolve_duplicates(metadata, [strain_col],
                [c for c in metadata if c != strain_col], cfg.dataset.duplicate_policy, what='metadata table')
        except ValueError as exc:
            raise ConfigError(str(exc)) from exc
        metadata = metadata.set_index(strain_col)
        metadata_join = {"rows_read": n_rows, "n_duplicate_ids": duplicates['n_duplicate_rows'],
            "strains_joined": int(core.isin(metadata.index).sum()), "strains_aligned": len(core),
            "n_metadata_only": int((~raw_metadata_ids.isin(core)).sum()),
            "metadata_only_examples": raw_metadata_ids[~raw_metadata_ids.isin(core)].tolist()[:5], **duplicates}
        if cfg.dataset.contrast_column and len(metadata):
            levels = metadata[cfg.dataset.contrast_column]
            for level in cfg.dataset.contrast_levels or ():
                if not levels.eq(str(level)).any():
                    raise ConfigError(f"dataset.contrast_levels names {level!r}, which no isolate carries in column {cfg.dataset.contrast_column!r}")
        for column, count in (metadata_quality(metadata).get("date_coerced_values") or {}).items():
            warnings.warn(f"metadata column {column!r} holds values that look like dates a spreadsheet produced from something else: {count}", RuntimeWarning, stacklevel=2)
    mic = mic_join = None
    if cfg.mic_path is not None:
        mic, mic_join = _load_mic(cfg, cfg.mic_path, core)
        named = [("mic_panel_column", cfg.dataset.mic_panel_column)]
        named += [("mic_covariate_columns", c) for c in cfg.dataset.mic_covariate_columns]
        for key, column in named:
            if column and column not in mic.columns and (metadata is None or column not in metadata.columns):
                raise ConfigError(f"dataset.{key} = {column!r} is a column of neither the MIC table nor the metadata")
        unknown = sorted(set(cfg.censored.calibration_agents) - set(mic_join['antimicrobials']))
        if unknown:
            raise ConfigError(f"censored.calibration_agents names {unknown}, which the MIC table does not hold; "
                              f"its agents are {mic_join['antimicrobials']}")
    qc = input_qc(panel, metadata=metadata, lineage_column=cfg.dataset.lineage_column, mic_join=mic_join)
    reading = panel.attrs.get('phenotype_reading')
    if reading is not None:
        reading['call_table_isolates_without_readable_calls'] = reading['isolates_without_readable_calls']
        reading['isolates_without_readable_calls'] = int(panel.isna().all(axis=1).sum())
        for r in reading['agents'].values():
            r['n_absent_records'] = len(core) - r['n_recorded_isolates']
        qc['phenotype_reading'] = reading
    qc['raw_universe'] = {'definition': 'union of raw call and MIC record identifiers',
        'n_call_isolates': len(call_ids), 'n_mic_isolates': len(mic_ids),
        'n_mic_only_isolates': int((~mic_ids.isin(call_ids)).sum()),
        'call_agents': list(panel.columns), 'mic_agents': mic_join['antimicrobials'] if mic_join else []}
    if metadata_join is not None:
        qc['metadata_join'].update(metadata_join)
    return Dataset(cfg=cfg, strain_ids=core, panel=panel, metadata=metadata,
                   mic=mic, mic_join=mic_join, input_qc=qc)


def _load_mic(cfg: Config, path, strain_ids: "pd.Index"):
    """Read dilutions without rewriting IDs, overwriting conflicts, or snapping wells."""
    ds = cfg.dataset
    idc, abc, vc = ds.mic_id_column, ds.mic_antibiotic_column, ds.mic_value_column
    df = _read_table(path, 'MIC table', dtype=str)
    for column in (idc, abc, vc):
        if column not in df:
            raise ConfigError(f"MIC table {path} has no column {column!r}")
    _require_keys(df, (idc, abc), f'MIC table {path}')
    opc = ds.mic_operator_column
    if opc and opc not in df:
        raise ConfigError(f"dataset.mic_operator_column = {opc!r} is not a column of {path}")
    raw = df[vc].map(_token)
    missing = raw.isin(_MISSING)
    sign = raw.str.extract(r'^(<=|>=|<|>)')[0].fillna('')
    numeric = pd.to_numeric(raw.str.replace(r'^(<=|>=|<|>)\s*', '', regex=True), errors='coerce')
    from .censored import _normalise_operators
    if opc:
        try:
            recorded = pd.Series(_normalise_operators(df[opc]), index=df.index)
        except ValueError as exc:
            raise ConfigError(f'MIC table {path}: {exc}') from exc
        disagree = recorded.ne('') & sign.ne('') & recorded.ne(sign)
        if disagree.any():
            raise ConfigError('MIC table has conflicting recorded and inline censoring operators')
        df[opc] = recorded.where(recorded.ne(''), sign)
    elif sign.ne('').any():
        opc = '_operator_from_value'
        df[opc] = sign
    invalid = ~missing & (~np.isfinite(numeric) | numeric.le(0))
    numeric = numeric.where(~invalid)
    for agent, wells in ds.mic_wells.items():
        sel = df[abc].eq(agent) & numeric.notna()
        values = numeric.loc[sel].to_numpy(dtype=float)
        if len(values) and not np.isclose(values[:, None], np.asarray(wells)[None, :], rtol=1e-10, atol=0).any(axis=1).all():
            raise ConfigError(f'MIC values for {agent!r} are outside configured wells; correct the data or panel specification (no snapping)')
    agents = list(df[abc].unique())
    registry = {}
    for agent in agents:
        sel = df[abc].eq(agent)
        registry[str(agent)] = {'rows_read': int(sel.sum()), 'n_recorded_isolates': int(df.loc[sel,idc].nunique()),
            'n_missing_values': int((sel & missing).sum()), 'n_unparseable_values': int((sel & invalid).sum())}
    rows_read = len(df)
    raw_ids = pd.Index(df[idc].unique())
    df[vc] = numeric
    # The panel and covariate labels of a reading are part of the reading: a
    # label that is blank is refused, and two rows that differ only in a label
    # are a conflict, not a duplicate.
    labels = [c for c in (ds.mic_panel_column, *ds.mic_covariate_columns) if c and c in df.columns]
    for column in labels:
        blank = df[column].map(_token).isin(_METADATA_MISSING) & numeric.notna()
        if blank.any():
            raise ConfigError(f"MIC table {path}: {int(blank.sum())} reading(s) without a value in column "
                              f"{column!r}; every reading needs one")
    # Preserve unreadable tokens in duplicate comparison; distinct errors are not equal observations.
    df['_duplicate_value'] = numeric.astype(str).where(numeric.notna(), 'unknown:' + raw).where(~missing, 'missing')
    try:
        df, duplicates = resolve_duplicates(df, [idc, abc], ['_duplicate_value'] + ([opc] if opc else []) + labels,
                                            ds.duplicate_policy, what='MIC table')
    except ValueError as exc:
        raise ConfigError(str(exc)) from exc
    df = df.drop(columns=['_duplicate_value'])
    wanted = pd.Index(strain_ids)
    joined = df.loc[df[idc].isin(wanted)].copy()
    matched = raw_ids[raw_ids.isin(wanted)]
    for agent in agents:
        r = registry[str(agent)]
        r['n_absent_records'] = len(wanted) - r['n_recorded_isolates']
        r['n_readable_isolates'] = int(joined.loc[joined[abc].eq(agent) & joined[vc].notna(), idc].nunique())
    try:
        shown = str(Path(path).resolve().relative_to(cfg.data_root.resolve()))
    except ValueError:
        shown = str(path)
    unparseable = {a: r['n_unparseable_values'] for a, r in registry.items() if r['n_unparseable_values']}
    report = {'path': shown, 'rows_read': rows_read, 'rows_joined': len(joined),
        'strains_with_mic': len(matched), 'strains_aligned': len(wanted),
        'join_rate': len(matched) / len(wanted) if len(wanted) else 0.,
        'unmatched_examples': wanted[~wanted.isin(matched)].tolist()[:5], 'antimicrobials': agents,
        'n_unparseable_values': int(invalid.sum()), 'n_missing_values': int(missing.sum()),
        'unparseable_values_by_agent': unparseable, 'operator_read_from_value': opc == '_operator_from_value',
        'mic_wells': ds.mic_wells, 'mic_units': ds.mic_units, 'agents': registry, **duplicates}
    if invalid.any():
        warnings.warn(f'MIC table {path} holds {int(invalid.sum())} value(s) that do not parse as a positive finite concentration and are set aside: {unparseable}. See mic_join.', RuntimeWarning, stacklevel=2)
    if len(wanted) and report['join_rate'] < 0.5:
        warnings.warn(f"MIC table joined {report['join_rate']:.1%} of the raw cohort; the censored share describes its readable subset. See mic_join.", RuntimeWarning, stacklevel=2)
    joined.attrs['mic_join'] = report
    return joined, report
