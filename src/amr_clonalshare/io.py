"""
io.py — config-driven data loading for amr-clonalshare.

The trait-clustering pipeline is built from one or more ``wide_binary`` layers
(strains × features, 0/1 presence/absence). This loader reads exactly those
layers, aligns them on a shared strain index according to
``dataset.strain_alignment_policy``, and ignores any other file kinds that may
appear in a shared project config (phenotype/series/long/newick) — they belong
to sibling tools, not to the clustering.

Alignment policy (config.dataset.strain_alignment_policy):
  * intersect_core — align on strains present in ALL loaded binary layers (default)
  * union          — union of all layer indices
  * strict_n       — like intersect_core but raise unless exactly expected_n align
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
import re
from typing import Dict, List, Optional

import pandas as pd

from .config import Config, ConfigError
from .qc import apply_missing_policy, input_qc, non_binary_values


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


def metadata_quality(metadata: "pd.DataFrame") -> Dict[str, object]:
    """Data-quality notes on the metadata, emitted with the run."""
    if metadata is None or not len(metadata.columns):
        return {}
    coerced = {}
    for column in metadata.columns:
        found = date_coerced_values(metadata[column])
        if found:
            coerced[str(column)] = found
    notes: Dict[str, object] = {}
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
    """Loaded, aligned binary layers for one cohort. Role keys come from the config."""
    cfg: Config
    frames: Dict[str, object] = field(default_factory=dict)  # role -> DataFrame
    strain_ids: Optional[pd.Index] = None
    metadata: Optional[pd.DataFrame] = None
    #: Long-format minimum inhibitory concentrations, one row per (isolate,
    #: antimicrobial), restricted to the aligned strains. ``None`` unless
    #: ``dataset.mic`` names a table.
    mic: Optional[pd.DataFrame] = None
    #: Share of the aligned strains carrying at least one MIC row, and the
    #: identifiers that did not join. Both are emitted with the run: a join
    #: that quietly matched half the cohort is the failure mode this reports.
    mic_join: Optional[Dict[str, object]] = None
    #: The input-check record from :func:`amr_clonalshare.qc.input_qc`:
    #: missing cells and the policy applied, per-trait adequacy, the metadata
    #: join, and per-lineage group sizes with the estimator's verdict.
    input_qc: Optional[Dict[str, object]] = None

    def get(self, role: str):
        if role not in self.frames:
            raise ConfigError(f"role {role!r} not loaded (available: {sorted(self.frames)})")
        return self.frames[role]

    def has(self, role: str) -> bool:
        return role in self.frames

    @property
    def wide_binary_roles(self) -> List[str]:
        return [r for r, spec in self.cfg.files.items()
                if spec.kind == "wide_binary" and r in self.frames]

    def binary(self, role: str) -> pd.DataFrame:
        df = self.get(role)
        if not isinstance(df, pd.DataFrame):
            raise ConfigError(f"role {role!r} is not a wide matrix")
        return df

    @property
    def n(self) -> int:
        return 0 if self.strain_ids is None else len(self.strain_ids)


def _read_wide(path, strain_col: str, *, strip_columns: bool) -> pd.DataFrame:
    try:
        df = pd.read_csv(path)
    except (FileNotFoundError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        raise ConfigError(f"layer file {path}: {exc}") from exc
    if strain_col not in df.columns:
        raise ConfigError(
            f"layer file {path} has no column {strain_col!r} "
            f"(dataset.strain_id_column); it has {list(df.columns)[:10]}")
    df = df.set_index(strain_col)
    if strip_columns:
        df.columns = df.columns.str.strip()
    if df.shape[1] == 0:
        raise ConfigError(f"layer file {path} holds no feature columns")
    return df


def load_dataset(cfg: Config) -> Dataset:
    """Load every ``wide_binary`` layer declared in ``cfg`` and align them.

    Non-binary file kinds are skipped (they belong to other tools).
    """
    strain_col = cfg.dataset.strain_id_column
    frames: Dict[str, object] = {}
    indices: List[pd.Index] = []

    for role, spec in cfg.files.items():
        if spec.kind not in ("wide_binary", "one_hot"):
            continue  # not consumed by the trait-clustering pipeline
        df = _read_wide(cfg.file_path(role), strain_col, strip_columns=spec.strip_columns)
        if df.index.has_duplicates:
            dup = df.index[df.index.duplicated()].unique()[:5].tolist()
            raise ConfigError(
                f"layer {role!r} has duplicate strain IDs (e.g. {dup}); "
                f"de-duplicate the input before clustering")
        frames[role] = df
        indices.append(df.index)

    if not indices:
        raise ConfigError("no wide_binary layers configured for trait clustering")

    policy = cfg.dataset.strain_alignment_policy
    core = indices[0]
    if policy == "union":
        for idx in indices[1:]:
            core = core.union(idx)
    else:  # intersect_core or strict_n
        for idx in indices[1:]:
            core = core.intersection(idx)

    if policy == "strict_n":
        if cfg.dataset.expected_n is not None and len(core) != cfg.dataset.expected_n:
            raise ConfigError(
                f"strict_n policy: aligned {len(core)} strains, expected "
                f"{cfg.dataset.expected_n}"
            )
    elif cfg.dataset.expected_n is not None and len(core) != cfg.dataset.expected_n:
        warnings.warn(
            f"aligned {len(core)} core strains, config expected_n="
            f"{cfg.dataset.expected_n} (policy={policy})",
            stacklevel=2,
        )

    # Reindexing on the union creates rows with no value; an empty CSV field
    # arrives the same way. Neither is a 0, so the policy is applied before
    # the cast and recorded with the run. See :mod:`amr_clonalshare.qc`.
    for role in list(frames):
        frames[role] = frames[role].reindex(core)
    frames, missing = apply_missing_policy(frames, cfg.dataset.missing_policy)
    if missing.dropped_rows:
        core = core[~core.astype(str).isin(missing.dropped_rows)]
    if len(core) == 0:
        raise ConfigError(
            f"no isolates remain after alignment ({policy}) and the missing "
            f"policy ({cfg.dataset.missing_policy}); nothing can be estimated")
    for role in list(frames):
        if frames[role].shape[1] == 0:
            raise ConfigError(
                f"layer {role!r} holds no feature columns after the missing "
                f"policy ({cfg.dataset.missing_policy}) removed them")
        bad = non_binary_values(frames[role])
        if bad:
            raise ConfigError(
                f"layer {role!r} holds values other than 0 and 1: {bad}; a "
                f"binary layer cannot carry them and the run stops")
        frames[role] = frames[role].loc[core].astype(int)

    # one_hot layers must actually be one-hot; a silent violation would make the
    # matching-coefficient distance and the count screening meaningless.
    for role, spec in cfg.files.items():
        if spec.kind == "one_hot" and role in frames:
            rs = frames[role].sum(axis=1)
            if not ((rs >= 0) & (rs <= 1)).all():
                bad = int((rs > 1).sum())
                raise ConfigError(
                    f"layer {role!r} is declared kind='one_hot' but {bad} rows "
                    f"have more than one positive column; split it into one "
                    f"layer per categorical variable")

    metadata = None
    mp = cfg.metadata_path
    if mp is not None and mp.exists():
        metadata = pd.read_csv(mp, dtype=str)
        if strain_col not in metadata.columns:
            raise ConfigError(
                f"metadata {mp} has no column {strain_col!r} "
                f"(dataset.strain_id_column); it has {list(metadata.columns)[:10]}")
        lineage_col = cfg.dataset.lineage_column
        if lineage_col and lineage_col not in metadata.columns:
            raise ConfigError(
                f"metadata {mp} has no column {lineage_col!r} "
                f"(dataset.lineage_column); it has {list(metadata.columns)[:10]}")
        metadata = metadata.drop_duplicates(subset=[strain_col]).set_index(strain_col)
        quality = metadata_quality(metadata)
        for value, count in (quality.get("date_coerced_values") or {}).items():
            warnings.warn(
                f"metadata column {value!r} holds values that look like dates "
                f"a spreadsheet produced from something else: {count}. See "
                f"metadata_quality in the run output.", RuntimeWarning,
                stacklevel=2)

    mic = mic_join = None
    qp = cfg.mic_path
    if qp is not None and qp.exists():
        mic, mic_join = _load_mic(cfg, qp, core)

    ds = Dataset(cfg=cfg, frames=frames, strain_ids=core, metadata=metadata,
                 mic=mic, mic_join=mic_join,
                 input_qc=input_qc(frames, missing, metadata=metadata,
                                   lineage_column=cfg.dataset.lineage_column))

    # Fail-loud pandera validation (G.0 Warstwa 1): every loaded wide_binary
    # layer must be strictly 0/1 with a unique strain index. Disable via
    # cfg.validation.schema_check=False when running on partially-cleaned data.
    if getattr(cfg.validation, "schema_check", True):
        try:
            from .schemas import validate_dataset as _validate_dataset
        except ImportError:  # pragma: no cover - pandera not installed
            _validate_dataset = None
        if _validate_dataset is not None:
            _validate_dataset(ds)

    return ds

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
    df = pd.read_csv(path, dtype={idc: str})
    for column in (idc, abc, vc):
        if column not in df.columns:
            raise ConfigError(
                f"MIC table {path} has no column {column!r}; it has "
                f"{sorted(df.columns)}")
    opc = cfg.dataset.mic_operator_column
    if opc and opc not in df.columns:
        raise ConfigError(
            f"dataset.mic_operator_column = {opc!r} is not a column of {path}")
    df[vc] = pd.to_numeric(df[vc], errors="coerce")
    wanted = pd.Index(strain_ids).astype(str)
    keep = df[idc].astype(str).isin(set(wanted))
    joined = df.loc[keep].copy()
    matched = sorted(set(joined[idc].astype(str)))
    unmatched = [s for s in wanted if s not in set(matched)]
    report = {
        "path": str(path),
        "rows_read": int(len(df)),
        "rows_joined": int(len(joined)),
        "strains_with_mic": len(matched),
        "strains_aligned": int(len(wanted)),
        "join_rate": (len(matched) / len(wanted)) if len(wanted) else 0.0,
        "unmatched_examples": unmatched[:5],
        "antimicrobials": sorted(str(a) for a in joined[abc].dropna().unique()),
    }
    if not matched:
        raise ConfigError(
            f"MIC table {path} joined none of the {len(wanted)} aligned "
            f"strains on column {idc!r}. Layer identifiers look like "
            f"{list(wanted[:3])}, MIC identifiers like "
            f"{list(df[idc].astype(str)[:3])}. Align them in the source table "
            f"rather than here.")
    if report["join_rate"] < 0.5:
        warnings.warn(
            f"MIC table joined {report['join_rate']:.1%} of the aligned "
            f"strains; the censored share will describe that subset. See "
            f"mic_join in the run output.", RuntimeWarning, stacklevel=2)
    return joined, report
