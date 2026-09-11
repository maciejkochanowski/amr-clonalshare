"""
io.py — config-driven data loading for amr-clonalshare.

The run takes a susceptibility table and a lineage column as input and does
not read sequence data: a phenotype per isolate, as a susceptibility call or
a recorded dilution, and a metadata table that holds the lineage label. The isolates of the run are the isolates the phenotype
table holds; the metadata and the dilution table are joined to them on the
identifier as written, and the join is measured and reported rather than
repaired.
"""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd

from .config import Config, ConfigError
from .phenotype import to_non_susceptible
from .qc import input_qc


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
    #: Share of the aligned strains carrying at least one MIC row, and the
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
        return pd.read_csv(path, **kwargs)
    except UnicodeDecodeError as exc:
        raise ConfigError(
            f"{what} {path} is not UTF-8 text (byte {exc.start}: {exc.reason}); "
            f"save it as UTF-8 CSV") from exc
    except (FileNotFoundError, pd.errors.ParserError,
            pd.errors.EmptyDataError, OSError) as exc:
        raise ConfigError(f"{what} {path}: {exc}") from exc


def read_panel(cfg: Config) -> pd.DataFrame:
    """The susceptibility panel as an isolate-by-antimicrobial 0/1 matrix."""
    fp = cfg.phenotype_path
    if fp is None:
        raise ConfigError("dataset.phenotype is not configured")
    df = _read_table(fp, "phenotype table",
                     dtype={cfg.dataset.phenotype_id_column: str})
    for column in (cfg.dataset.phenotype_id_column,
                   cfg.dataset.phenotype_antibiotic_column,
                   cfg.dataset.phenotype_call_column):
        if column not in df.columns:
            raise ConfigError(
                f"phenotype table {fp} has no column {column!r}; it has "
                f"{list(df.columns)[:10]}")
    panel = to_non_susceptible(
        df, id_column=cfg.dataset.phenotype_id_column,
        antibiotic_column=cfg.dataset.phenotype_antibiotic_column,
        call_column=cfg.dataset.phenotype_call_column,
        intermediate=cfg.dataset.phenotype_intermediate)
    panel.index = panel.index.astype(str)
    panel.index.name = cfg.dataset.strain_id_column
    if not len(panel.index):
        raise ConfigError(f"phenotype table {fp} names no isolates with a "
                          f"readable call")
    return panel


def _dilution_index(cfg: Config) -> pd.Index:
    """The isolates of the dilution table, for a run with no call table."""
    fp = cfg.mic_path
    idc = cfg.dataset.mic_id_column
    df = _read_table(fp, "MIC table", dtype={idc: str})
    if idc not in df.columns:
        raise ConfigError(f"MIC table {fp} has no column {idc!r}")
    ids = pd.Index(df[idc].astype(str).dropna().unique(),
                   name=cfg.dataset.strain_id_column)
    if not len(ids):
        raise ConfigError(f"MIC table {fp} names no isolates")
    return ids


def load_dataset(cfg: Config) -> Dataset:
    """Read the phenotype table, the metadata and the dilutions and join them."""
    strain_col = cfg.dataset.strain_id_column
    if cfg.dataset.phenotype:
        panel = read_panel(cfg)
        core = panel.index
    else:
        core = _dilution_index(cfg)
        panel = pd.DataFrame(index=core)

    metadata = None
    metadata_join: Optional[Dict[str, Any]] = None
    mp = cfg.metadata_path
    if mp is not None and mp.exists():
        metadata = _read_table(mp, "metadata", dtype=str)
        if strain_col not in metadata.columns:
            raise ConfigError(
                f"metadata {mp} has no column {strain_col!r} "
                f"(dataset.strain_id_column); it has {list(metadata.columns)[:10]}")
        for key, column in (("lineage_column", cfg.dataset.lineage_column),
                            ("contrast_column", cfg.dataset.contrast_column),
                            ("batch_column", cfg.dataset.batch_column)):
            if column and column not in metadata.columns:
                raise ConfigError(
                    f"metadata {mp} has no column {column!r} "
                    f"(dataset.{key}); it has {list(metadata.columns)[:10]}")
        n_rows = int(len(metadata))
        metadata = metadata.drop_duplicates(subset=[strain_col]).set_index(strain_col)
        n_duplicate_ids = n_rows - int(len(metadata))
        if n_duplicate_ids:
            warnings.warn(
                f"metadata {mp} repeats {n_duplicate_ids} isolate identifier(s); "
                f"the first row of each is kept. See metadata_join in the run "
                f"output.", RuntimeWarning, stacklevel=2)
        joined_ids = int(core.astype(str).isin(metadata.index.astype(str)).sum())
        metadata_join = {"rows_read": n_rows, "n_duplicate_ids": n_duplicate_ids,
                         "strains_joined": joined_ids,
                         "strains_aligned": int(len(core))}
        if not joined_ids:
            raise ConfigError(
                f"metadata {mp} joined none of the {len(core)} isolates of the "
                f"phenotype table on column {strain_col!r}. Phenotype "
                f"identifiers look like {list(core.astype(str)[:3])}, metadata "
                f"identifiers like {list(metadata.index.astype(str)[:3])}. Align "
                f"them in the source table rather than here.")
        if cfg.dataset.contrast_column:
            levels = metadata[cfg.dataset.contrast_column].astype(str)
            for level in cfg.dataset.contrast_levels or ():
                if not int((levels == str(level)).sum()):
                    raise ConfigError(
                        f"dataset.contrast_levels names {level!r}, which no "
                        f"isolate carries in column "
                        f"{cfg.dataset.contrast_column!r}; the levels present "
                        f"are {sorted(levels.dropna().unique())[:10]}")
        quality = metadata_quality(metadata)
        for value, count in (quality.get("date_coerced_values") or {}).items():
            warnings.warn(
                f"metadata column {value!r} holds values that look like dates "
                f"a spreadsheet produced from something else: {count}; "
                f"amr_clonalshare.io.metadata_quality lists them.", RuntimeWarning,
                stacklevel=2)

    mic = mic_join = None
    qp = cfg.mic_path
    if qp is not None and qp.exists():
        mic, mic_join = _load_mic(cfg, qp, core)

    if len(core) < 4:
        raise ConfigError(
            f"the tables name {len(core)} isolate(s); a share needs at least "
            f"four, two lineages of two, before anything can be held out")

    qc = input_qc(panel, metadata=metadata,
                  lineage_column=cfg.dataset.lineage_column, mic_join=mic_join)
    if metadata_join is not None:
        qc["metadata_join"].update(
            {k: v for k, v in metadata_join.items() if k not in qc["metadata_join"]})
    return Dataset(cfg=cfg, strain_ids=core, panel=panel, metadata=metadata,
                   mic=mic, mic_join=mic_join, input_qc=qc)


def _load_mic(cfg: Config, path, strain_ids: "pd.Index"):
    """Read the long-format MIC table and restrict it to the aligned strains.

    The identifier column is joined as given. Nothing is stripped, prefixed or
    case-folded: a loader that repaired identifiers would be a loader that
    could join the wrong isolates without saying so, and the cost of that
    error is a variance component attributed to the wrong lineage. What the
    loader does instead is measure the join and report it, and refuse only
    when nothing matched at all, which is always a configuration error rather
    than a property of the cohort.
    """
    idc = cfg.dataset.mic_id_column
    abc = cfg.dataset.mic_antibiotic_column
    vc = cfg.dataset.mic_value_column
    df = _read_table(path, "MIC table", dtype={idc: str})
    for column in (idc, abc, vc):
        if column not in df.columns:
            raise ConfigError(
                f"MIC table {path} has no column {column!r}; it has "
                f"{sorted(df.columns)}")
    opc = cfg.dataset.mic_operator_column
    if opc and opc not in df.columns:
        raise ConfigError(
            f"dataset.mic_operator_column = {opc!r} is not a column of {path}")
    # A laboratory export writes a censored reading as "<=0.5" or ">32" in
    # the value column. The sign is read off and kept as the operator where
    # no operator column is configured; a value that still does not parse is
    # counted per agent and reported, never silently dropped.
    raw = df[vc].astype(str).str.strip()
    sign = raw.str.extract(r"^(<=|>=|<|>)")[0]
    stripped = raw.str.replace(r"^(<=|>=|<|>)\s*", "", regex=True)
    numeric = pd.to_numeric(stripped, errors="coerce")
    if not opc and sign.notna().any():
        df["_operator_from_value"] = sign.where(sign.notna(), "=")
        opc = "_operator_from_value"
    df[vc] = numeric
    present = df[vc].isna() & raw.ne("") & raw.str.lower().ne("nan")
    unparseable = (df.loc[present, abc].astype(str).value_counts().to_dict()
                   if present.any() else {})
    if unparseable:
        warnings.warn(
            f"MIC table {path} holds {int(present.sum())} value(s) that do not "
            f"parse as a concentration and are set aside: "
            f"{dict(list(unparseable.items())[:5])}. See mic_join in the run "
            f"output.", RuntimeWarning, stacklevel=2)
    wanted = pd.Index(strain_ids).astype(str)
    keep = df[idc].astype(str).isin(set(wanted))
    joined = df.loc[keep].copy()
    matched = sorted(set(joined[idc].astype(str)))
    unmatched = [s for s in wanted if s not in set(matched)]
    # The path is recorded relative to the configuration, so that a record
    # written on one machine names the same file on another.
    try:
        shown = str(Path(path).resolve().relative_to(cfg.data_root.resolve()))
    except ValueError:
        shown = str(path)
    report: Dict[str, Any] = {
        "path": shown,
        "rows_read": int(len(df)),
        "rows_joined": int(len(joined)),
        "strains_with_mic": len(matched),
        "strains_aligned": int(len(wanted)),
        "join_rate": (len(matched) / len(wanted)) if len(wanted) else 0.0,
        "unmatched_examples": unmatched[:5],
        "antimicrobials": sorted(str(a) for a in joined[abc].dropna().unique()),
        "n_unparseable_values": int(present.sum()),
        "unparseable_values_by_agent": {str(k): int(v) for k, v in unparseable.items()},
        "operator_read_from_value": bool(opc == "_operator_from_value"),
    }
    if not matched:
        raise ConfigError(
            f"MIC table {path} joined none of the {len(wanted)} aligned "
            f"strains on column {idc!r}. Phenotype identifiers look like "
            f"{list(wanted[:3])}, MIC identifiers like "
            f"{list(df[idc].astype(str)[:3])}. Align them in the source table "
            f"rather than here.")
    if report["join_rate"] < 0.5:
        warnings.warn(
            f"MIC table joined {report['join_rate']:.1%} of the aligned "
            f"strains; the censored share will describe that subset. See "
            f"mic_join in the run output.", RuntimeWarning, stacklevel=2)
    return joined, report
