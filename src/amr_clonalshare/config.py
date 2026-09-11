"""config.py — declarative configuration for amr-clonalshare.

Everything that changes the science is a config key, and every default is
stated here rather than buried in a function signature. Sections:

``dataset``       cohort identity, the phenotype table, the metadata table with
                  its lineage column, and the optional dilution table, contrast
                  and intake columns
``attribution``   the budgets of the lineage-membership share
``surveillance``  the budgets and gates of the prevalence readings and the
                  decomposition
``censored``      the reading of a dilution panel
``evidence``      the e-value and its rejection level

``data_dir`` is resolved relative to the config file, so a config can sit next
to its data. Validation is eager: :func:`load_config` raises :class:`ConfigError`
with an actionable message instead of failing deep inside the analysis.
"""
from __future__ import annotations

import difflib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "PyYAML is required to load configs. Install with `pip install pyyaml` "
        "or `pip install -e .` from the package root."
    ) from exc


class ConfigError(ValueError):
    """Raised when a config is structurally invalid or internally inconsistent."""


@dataclass(frozen=True)
class DatasetConfig:
    name: str
    strain_id_column: str = "Strain_ID"
    data_dir: str = "."
    metadata: Optional[str] = None          # path to a metadata CSV
    lineage_column: Optional[str] = None    # ST / clonal group column therein
    # Two collections to compare when decomposing a prevalence difference into
    # a lineage-mix component and a within-lineage rate component. A period, a
    # country, a host. Left unset, the decomposition is simply not run.
    contrast_column: Optional[str] = None
    contrast_levels: Sequence[str] = ()
    #: Metadata column that names the intake an isolate arrived in, a year of
    #: collection or a batch number. Set, the run also reads each trait as a
    #: programme of looks: the batches are taken in ascending order of this
    #: column's values, so the column must sort into arrival order, and the
    #: running product of their e-values is a test supermartingale that may be
    #: read after any intake. Left unset, only the single-look e-value is
    #: computed.
    batch_column: Optional[str] = None
    #: Long-format CSV of susceptibility calls, one row per (isolate,
    #: antibiotic). The isolates of the run are the isolates this table
    #: holds. See :mod:`amr_clonalshare.phenotype` for the vocabulary read.
    phenotype: Optional[str] = None
    phenotype_id_column: str = "Strain_ID"
    phenotype_antibiotic_column: str = "antibiotic"
    phenotype_call_column: str = "resistant_phenotype"
    phenotype_intermediate: str = "non_susceptible"
    #: Optional long-format CSV of measured minimum inhibitory
    #: concentrations, one row per (isolate, antimicrobial). A dichotomised
    #: call keeps only which side of a cut-off an isolate fell; the recorded
    #: dilution keeps where in the panel it fell, and the difference is
    #: measurable rather than rhetorical. See
    #: :mod:`amr_clonalshare.censored`. The identifier column must already
    #: match the identifiers of the phenotype table: this loader joins, it
    #: does not repair identifiers, because a rule that silently rewrote them
    #: would be a rule that silently joined the wrong isolates.
    mic: Optional[str] = None
    mic_id_column: str = "Strain_ID"
    mic_antibiotic_column: str = "antibiotic"
    mic_value_column: str = "measurement"
    #: Column holding a recorded censoring operator (``<``, ``<=``, ``>``,
    #: ``>=``). Where it is present it decides; where it is absent the panel
    #: geometry decides, and that is an assumption the run reports.
    mic_operator_column: Optional[str] = None

    def validate(self) -> None:
        if not (self.phenotype or self.mic):
            raise ConfigError(
                "dataset.phenotype or dataset.mic must name the table the "
                "share is read from")
        if not self.metadata:
            raise ConfigError("dataset.metadata must name the table that "
                              "holds the lineage column")
        if not self.lineage_column:
            raise ConfigError("dataset.lineage_column must name the lineage "
                              "label: the share is the share that label carries")
        if self.contrast_column and len(self.contrast_levels) != 2:
            raise ConfigError(
                "dataset.contrast_levels must name exactly two values of "
                f"{self.contrast_column!r}; got "
                f"{list(self.contrast_levels)!r}")
        if self.phenotype_intermediate not in ("non_susceptible", "drop",
                                               "susceptible"):
            raise ConfigError(
                "dataset.phenotype_intermediate must be one of "
                "'non_susceptible', 'drop', 'susceptible'; got "
                f"{self.phenotype_intermediate!r}")


@dataclass(frozen=True)
class AttributionConfig:
    """Budgets of the lineage-membership share.

    Every figure is a cross-validated statistic averaged over ``repeats``
    fold draws, so the cost is ``repeats + n_perm + n_boot`` model fits per
    antimicrobial rather than one, and 2000 of each would put a single run
    into hours for no gain in the third decimal place. The surveillance
    readings keep their own budgets because they are different resamplings
    and do not want the same count.
    """
    enabled: bool = True
    folds: int = 5
    repeats: int = 20
    n_boot: int = 400
    n_perm: int = 200

    def validate(self) -> None:
        if self.folds < 2:
            raise ConfigError("attribution.folds must be >= 2")
        if self.repeats < 1:
            raise ConfigError("attribution.repeats must be >= 1")
        if self.n_boot < 0 or self.n_perm < 1:
            raise ConfigError("attribution.n_boot must be >= 0 and n_perm >= 1; "
                              "the permuted-label run is what the share is "
                              "corrected against")


@dataclass(frozen=True)
class SurveillanceConfig:
    """Budgets and thresholds for the lineage-resolved prevalence outputs.

    They are separate from the attribution budgets. A percentile interval
    and a permutation tail are not the same resampling and do not want the
    same count: 500 replicates put the 2.5th percentile on the 13th draw, and 500
    permutations put a p-value floor at 1/501, which a panel of thirteen agents
    then reports thirteen times.
    """
    enabled: bool = True
    n_boot: int = 2000
    n_perm: int = 2000
    q_fdr: float = 0.05
    min_shared_support: float = 0.8
    label_alpha: float = 0.05

    def validate(self) -> None:
        if self.n_boot < 1 or self.n_perm < 1:
            raise ConfigError("surveillance.n_boot and n_perm must be >= 1")
        if not 0 < self.q_fdr < 1:
            raise ConfigError("surveillance.q_fdr must be in (0, 1)")
        if not 0 <= self.min_shared_support <= 1:
            raise ConfigError(
                "surveillance.min_shared_support must be in [0, 1]")
        if not 0 < self.label_alpha < 1:
            raise ConfigError("surveillance.label_alpha must be in (0, 1)")


@dataclass(frozen=True)
class CensoredConfig:
    """The interval-censored reading of a dilution panel.

    Has no effect unless ``dataset.mic`` names a table, because there is
    nothing to read without one. ``n_boot`` is the only budget worth setting:
    the estimator is one expectation-maximisation fit per resample rather than
    a cross-validated average, so a run costs one fit plus ``n_boot``, which
    is why the default is affordable where the attribution defaults are not.

    ``end_wells_censored`` states what a reading on the lowest or highest
    tested well means. Treating it as censored is the coarsened-at-random
    reading of Heitjan and Rubin; treating it as exact is the alternative.
    With ``sensitivity`` on, both are computed and both are reported, so the
    width of the bracket says how much of the answer rests on the choice.
    """
    enabled: bool = True
    n_boot: int = 200
    end_wells_censored: bool = True
    sensitivity: bool = True

    def validate(self) -> None:
        if self.n_boot < 0:
            raise ConfigError("censored.n_boot must be >= 0")


@dataclass(frozen=True)
class EvidenceConfig:
    """Anytime-valid evidence that a trait depends on lineage.

    Surveillance re-reads the same panel every time a year of isolates
    arrives. A p-value recomputed at each look has no error guarantee, because
    the number of looks is not fixed in advance; an e-value does, and may be
    inspected as often as wanted. ``alpha`` sets the rejection threshold at
    ``1 / alpha`` by Ville's inequality, not by a tail probability.

    ``folds`` and ``repeats`` mirror the attribution block because the
    evidence reuses that split structure. They are separate settings so that a
    programme can spend on the quantity it re-inspects without paying the same
    on the one it does not.
    """
    enabled: bool = True
    folds: int = 5
    repeats: int = 20
    alpha: float = 0.05

    def validate(self) -> None:
        if self.folds < 2:
            raise ConfigError("evidence.folds must be >= 2")
        if self.repeats < 1:
            raise ConfigError("evidence.repeats must be >= 1")
        if not 0 < self.alpha < 1:
            raise ConfigError("evidence.alpha must be in (0, 1)")


@dataclass(frozen=True)
class Config:
    """Fully-resolved, validated configuration for one cohort."""
    dataset: DatasetConfig
    attribution: AttributionConfig = field(
        default_factory=AttributionConfig)
    surveillance: SurveillanceConfig = field(
        default_factory=SurveillanceConfig)
    censored: CensoredConfig = field(default_factory=CensoredConfig)
    evidence: EvidenceConfig = field(default_factory=EvidenceConfig)
    config_path: Optional[Path] = None

    @property
    def data_root(self) -> Path:
        base = self.config_path.parent if self.config_path else Path.cwd()
        return (base / self.dataset.data_dir).resolve()

    @property
    def metadata_path(self) -> Optional[Path]:
        if not self.dataset.metadata:
            return None
        return self.data_root / self.dataset.metadata

    @property
    def phenotype_path(self) -> Optional[Path]:
        if not self.dataset.phenotype:
            return None
        return self.data_root / self.dataset.phenotype

    @property
    def mic_path(self) -> Optional[Path]:
        if not self.dataset.mic:
            return None
        return self.data_root / self.dataset.mic

    def validate(self, check_files_exist: bool = True) -> "Config":
        self.dataset.validate()
        self.attribution.validate()
        self.surveillance.validate()
        self.censored.validate()
        self.evidence.validate()
        if check_files_exist:
            missing = [f"{name} -> {p}" for name, p in (
                ("metadata", self.metadata_path),
                ("phenotype", self.phenotype_path),
                ("mic", self.mic_path)) if p is not None and not p.exists()]
            if missing:
                raise ConfigError(
                    "configured input files not found:\n  " + "\n  ".join(missing))
        return self


def _sub(raw: Dict[str, Any], key: str) -> Dict[str, Any]:
    val = (raw or {}).get(key) or {}
    if not isinstance(val, dict):
        raise ConfigError(f"config section {key!r} must be a mapping")
    return val


# --------------------------------------------------------------------------- #
#  Strict key checking
# --------------------------------------------------------------------------- #
#: Every key the loader reads, by section. A key that is not here is refused.
#:
#: A misspelt key is the dangerous case, not an invented one. ``lineage_colum``
#: silently disables the lineage diagnostic, and the run then reports no
#: population-structure gate for the same reason it would report none on a
#: cohort with no lineage labels at all: the two are indistinguishable in the
#: output. Refusing the key is the only way a configuration can be a contract
#: rather than a suggestion.
_KNOWN_KEYS: Dict[str, frozenset] = {
    "": frozenset({"dataset", "attribution", "surveillance", "censored",
                   "evidence"}),
    "dataset": frozenset({
        "name", "strain_id_column", "data_dir", "metadata", "lineage_column",
        "contrast_column", "contrast_levels", "batch_column", "phenotype",
        "phenotype_id_column", "phenotype_antibiotic_column",
        "phenotype_call_column", "phenotype_intermediate", "mic",
        "mic_id_column", "mic_antibiotic_column", "mic_value_column",
        "mic_operator_column"}),
    "attribution": frozenset({"enabled", "folds", "repeats", "n_boot",
                              "n_perm"}),
    "surveillance": frozenset({"enabled", "n_boot", "n_perm", "q_fdr",
                               "min_shared_support", "label_alpha"}),
    "censored": frozenset({"enabled", "n_boot", "end_wells_censored",
                           "sensitivity"}),
    "evidence": frozenset({"enabled", "folds", "repeats", "alpha"}),
}


def _reject_unknown_keys(mapping: Any, section: str = "") -> None:
    """Refuse a key the loader does not read, naming the nearest known one."""
    if not isinstance(mapping, dict):
        return
    allowed = _KNOWN_KEYS.get(section)
    if allowed is None:
        return
    where = f"section {section!r}" if section else "the config root"
    for key in mapping:
        if key in allowed:
            continue
        close = difflib.get_close_matches(str(key), sorted(allowed), n=1,
                                          cutoff=0.7)
        suggestion = f"; did you mean {close[0]!r}?" if close else "."
        raise ConfigError(
            f"unknown key {key!r} in {where}{suggestion} A key that is not "
            f"read is a silently disabled analysis, so it is refused rather "
            f"than ignored. Known keys: {sorted(allowed)}")
    for key in mapping:
        child = f"{section}.{key}" if section else str(key)
        if child in _KNOWN_KEYS:
            _reject_unknown_keys(mapping[key], child)


def _levels(value) -> tuple:
    """``contrast_levels`` as a list of two labels; a bare string is refused
    rather than read letter by letter."""
    if value is None:
        return ()
    if isinstance(value, (str, bytes)) or not isinstance(value, (list, tuple)):
        raise ConfigError(
            f"dataset.contrast_levels must be a list of two labels, got {value!r}")
    return tuple(str(v) for v in value)


def _int(section: str, raw: Dict[str, Any], key: str, default: int) -> int:
    """An integer option, read strictly: a float with a fraction or a string
    that is not a whole number is a configuration error, not a rounding."""
    value = raw.get(key, default)
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ConfigError(f"{section}.{key} must be a whole number, got {value!r}")
    try:
        number = float(value)
    except ValueError as exc:
        raise ConfigError(f"{section}.{key} must be a whole number, got {value!r}") from exc
    if number != int(number):
        raise ConfigError(f"{section}.{key} must be a whole number, got {value!r}")
    return int(number)


def _float(section: str, raw: Dict[str, Any], key: str, default: float) -> float:
    value = raw.get(key, default)
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ConfigError(f"{section}.{key} must be a number, got {value!r}")
    try:
        return float(value)
    except ValueError as exc:
        raise ConfigError(f"{section}.{key} must be a number, got {value!r}") from exc


def _bool(section: str, raw: Dict[str, Any], key: str, default: bool) -> bool:
    """A yes/no option. YAML already reads ``true``/``false``; a quoted
    ``"false"`` is a common slip and is read as the word it is rather than as
    a non-empty string."""
    value = raw.get(key, default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.strip().lower() in ("true", "yes", "on"):
        return True
    if isinstance(value, str) and value.strip().lower() in ("false", "no", "off"):
        return False
    raise ConfigError(f"{section}.{key} must be true or false, got {value!r}")


def from_dict(raw: Dict[str, Any], config_path: Optional[Path] = None) -> Config:
    """Build a validated :class:`Config` from a parsed YAML/JSON mapping."""
    if "dataset" not in raw:
        raise ConfigError("config must define at least: dataset")
    _reject_unknown_keys(raw, "")

    ds = raw["dataset"]
    if not isinstance(ds, dict):
        raise ConfigError(f"config section 'dataset' must be a mapping, got "
                          f"{type(ds).__name__}; it names the cohort and the "
                          f"tables the run reads")
    if "name" not in ds:
        raise ConfigError("dataset.name must name the cohort; it is the "
                          "identifier the record and both reports carry")
    dataset = DatasetConfig(
        name=ds["name"],
        strain_id_column=ds.get("strain_id_column", "Strain_ID"),
        data_dir=ds.get("data_dir", "."),
        metadata=ds.get("metadata"),
        lineage_column=ds.get("lineage_column"),
        contrast_column=ds.get("contrast_column"),
        contrast_levels=_levels(ds.get("contrast_levels")),
        batch_column=ds.get("batch_column"),
        phenotype=ds.get("phenotype"),
        phenotype_id_column=ds.get("phenotype_id_column", "Strain_ID"),
        phenotype_antibiotic_column=ds.get("phenotype_antibiotic_column", "antibiotic"),
        phenotype_call_column=ds.get("phenotype_call_column", "resistant_phenotype"),
        phenotype_intermediate=ds.get("phenotype_intermediate", "non_susceptible"),
        mic=ds.get("mic"),
        mic_id_column=ds.get("mic_id_column", "Strain_ID"),
        mic_antibiotic_column=ds.get("mic_antibiotic_column", "antibiotic"),
        mic_value_column=ds.get("mic_value_column", "measurement"),
        mic_operator_column=ds.get("mic_operator_column"),
    )

    attr_raw = _sub(raw, "attribution")
    attribution = AttributionConfig(
        enabled=_bool("attribution", attr_raw, "enabled", True),
        folds=_int("attribution", attr_raw, "folds", 5),
        repeats=_int("attribution", attr_raw, "repeats", 20),
        n_boot=_int("attribution", attr_raw, "n_boot", 400),
        n_perm=_int("attribution", attr_raw, "n_perm", 200),
    )

    surv_raw = _sub(raw, "surveillance")
    surveillance = SurveillanceConfig(
        enabled=_bool("surveillance", surv_raw, "enabled", True),
        n_boot=_int("surveillance", surv_raw, "n_boot", 2000),
        n_perm=_int("surveillance", surv_raw, "n_perm", 2000),
        q_fdr=_float("surveillance", surv_raw, "q_fdr", 0.05),
        min_shared_support=_float("surveillance", surv_raw, "min_shared_support", 0.8),
        label_alpha=_float("surveillance", surv_raw, "label_alpha", 0.05),
    )

    cen_raw = _sub(raw, "censored")
    censored = CensoredConfig(
        enabled=_bool("censored", cen_raw, "enabled", True),
        n_boot=_int("censored", cen_raw, "n_boot", 200),
        end_wells_censored=_bool("censored", cen_raw, "end_wells_censored", True),
        sensitivity=_bool("censored", cen_raw, "sensitivity", True),
    )

    ev_raw = _sub(raw, "evidence")
    evidence = EvidenceConfig(
        enabled=_bool("evidence", ev_raw, "enabled", True),
        folds=_int("evidence", ev_raw, "folds", 5),
        repeats=_int("evidence", ev_raw, "repeats", 20),
        alpha=_float("evidence", ev_raw, "alpha", 0.05),
    )

    return Config(
        dataset=dataset,
        attribution=attribution,
        surveillance=surveillance,
        censored=censored,
        evidence=evidence,
        config_path=Path(config_path).resolve() if config_path else None,
    )


def load_config(path: "str | Path", check_files_exist: bool = True) -> Config:
    """Load and validate a YAML config file. Returns a frozen, validated Config."""
    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise ConfigError(f"config file not found: {p}")
    with p.open("r", encoding="utf-8") as fh:
        try:
            raw = yaml.safe_load(fh)
        except yaml.YAMLError as exc:
            raise ConfigError(f"config file {p} is not valid YAML: {exc}") from exc
    if not isinstance(raw, dict):
        raise ConfigError(f"config root must be a mapping, got {type(raw).__name__}")
    return from_dict(raw, config_path=p).validate(check_files_exist=check_files_exist)
