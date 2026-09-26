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
``population_probit`` optional population-model ICC profile, disabled by default

``data_dir`` is resolved relative to the config file, so a config can sit next
to its data. Validation is eager: :func:`load_config` raises :class:`ConfigError`
with an actionable message instead of failing deep inside the analysis.
"""
from __future__ import annotations

import difflib
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple

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
    #: antibiotic). The raw cohort is the union of this table's IDs and MIC
    #: IDs. See :mod:`amr_clonalshare.phenotype` for the vocabulary read.
    phenotype: Optional[str] = None
    phenotype_id_column: str = "Strain_ID"
    phenotype_antibiotic_column: str = "antibiotic"
    phenotype_call_column: str = "resistant_phenotype"
    #: Names the antimicrobial columns of a wide call table, one column per
    #: agent and one row per isolate, as a laboratory spreadsheet is usually
    #: kept. The table is melted to the long shape before it is read, so the
    #: vocabulary, the duplicate policy and every count are unchanged. Leave it
    #: unset for the long shape.
    phenotype_agent_columns: Tuple[str, ...] = ()
    phenotype_intermediate: Optional[str] = None
    phenotype_kind: str = "undeclared"
    phenotype_positive_definition: Optional[str] = None
    phenotype_source: Optional[str] = None
    ast_standard: Optional[str] = None
    ast_version: Optional[str] = None
    duplicate_policy: str = "error"
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
    #: Column naming the panel a reading was made on, in the MIC table or,
    #: per isolate, in the metadata: the laboratory, the panel product, or
    #: any label under which every isolate was tested on the same wells. The
    #: tested wells, and with them the end wells, are then inferred within
    #: each label, so that the lowest well of one laboratory is not read as an
    #: interior dilution because another laboratory tested lower.
    mic_panel_column: Optional[str] = None
    #: Columns naming categorical covariates of the MIC model, in the MIC
    #: table or, per isolate, in the metadata: the testing laboratory, the
    #: country, the year. Their levels enter the calibrated interval as
    #: additive fixed effects, estimated with the other parameters and
    #: re-estimated in every simulated dataset, so that a shift between
    #: levels is removed from the variance components instead of being
    #: absorbed into whichever of them it is aligned with. The moment
    #: estimator does not use them. ``mic_covariate_column`` names one.
    mic_covariate_columns: Tuple[str, ...] = ()
    #: Metadata column whose levels define strata: every analysis is then
    #: repeated on the isolates of each level, with its own denominators,
    #: gates and intervals, and the run adds one table across the strata.
    #: The whole-collection analysis is unchanged.
    stratify_by: Optional[str] = None
    mic_wells: Dict[str, Sequence[float]] = field(default_factory=dict)
    mic_units: Dict[str, str] = field(default_factory=dict)

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
        if self.contrast_column and len(set(self.contrast_levels)) != 2:
            raise ConfigError("dataset.contrast_levels must name two distinct values")
        if self.contrast_levels and not self.contrast_column:
            raise ConfigError("dataset.contrast_levels needs dataset.contrast_column to name the column they are values of")
        if not isinstance(self.name, str) or not self.name.strip():
            raise ConfigError("dataset.name must be a non-empty string")
        for key in ("mic_panel_column", "stratify_by"):
            value = getattr(self, key)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ConfigError(f"dataset.{key} must name one column")
        if self.phenotype_intermediate not in (None, "non_susceptible", "drop",
                                               "susceptible"):
            raise ConfigError(
                "dataset.phenotype_intermediate must be one of "
                "'non_susceptible', 'drop', 'susceptible'; got "
                f"{self.phenotype_intermediate!r}")
        if self.phenotype_agent_columns:
            if not self.phenotype:
                raise ConfigError("dataset.phenotype_agent_columns names columns of "
                                  "dataset.phenotype, which is not set")
            if len(set(self.phenotype_agent_columns)) != len(self.phenotype_agent_columns):
                raise ConfigError("dataset.phenotype_agent_columns repeats a column")
            if self.phenotype_id_column in self.phenotype_agent_columns:
                raise ConfigError("dataset.phenotype_agent_columns must not name "
                                  "dataset.phenotype_id_column")
        if self.phenotype_kind not in ("undeclared", "clinical_sir", "wt_nwt", "binary"):
            raise ConfigError("dataset.phenotype_kind must be undeclared, clinical_sir, wt_nwt or binary")
        if self.phenotype_kind in ("wt_nwt", "binary") and self.phenotype_intermediate is not None:
            raise ConfigError("dataset.phenotype_intermediate applies only to clinical_sir or undeclared")
        if self.duplicate_policy not in ("error", "drop_conflicts", "positive_wins"):
            raise ConfigError("dataset.duplicate_policy must be error, drop_conflicts or positive_wins")
        for name in ("phenotype_positive_definition", "phenotype_source", "ast_standard", "ast_version"):
            value = getattr(self, name)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ConfigError(f"dataset.{name} must be a nonempty string or null")
        if not isinstance(self.mic_wells, dict):
            raise ConfigError("dataset.mic_wells must map agent names to tested concentrations")
        for agent, wells in self.mic_wells.items():
            valid = isinstance(agent, str) and bool(agent.strip()) and isinstance(wells, (list, tuple)) and len(wells) > 0
            if valid:
                valid = all(not isinstance(v, bool) and isinstance(v, (int, float)) and math.isfinite(v) and v > 0 for v in wells)
            if valid:
                valid = all(a < b for a, b in zip(wells, wells[1:]))
            if not valid:
                raise ConfigError("dataset.mic_wells requires agent names and strictly increasing positive finite concentrations")
        if not isinstance(self.mic_units, dict) or any(not isinstance(k, str) or not k.strip() or not isinstance(v, str) or not v.strip() for k, v in self.mic_units.items()):
            raise ConfigError("dataset.mic_units must map agent names to nonempty unit strings")


@dataclass(frozen=True)
class AttributionConfig:
    """Budgets of the lineage-membership share.

    Every figure is a cross-validated statistic averaged over ``repeats``
    fold draws, so the cost is ``repeats + n_perm + n_boot`` model fits per
    antimicrobial rather than one, and 2000 of each would put a single run
    into hours for no gain in the third decimal place. The surveillance
    readings keep their own budgets because they are different resamplings
    and do not want the same count. ``n_perm`` is larger because the
    permutation p-values select the agents with a lineage effect
    (``lineage_selection``): the smallest attainable p-value, 1/(n_perm + 1),
    has to fall below the step-up's first threshold for a single agent to be
    selectable, and 999 permutations clear it for a panel of up to 15
    agents; a larger panel needs more.
    """
    enabled: bool = True
    folds: int = 5
    repeats: int = 20
    n_boot: int = 400
    n_perm: int = 999

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
    """Budget and thresholds for the decomposition of a prevalence difference.

    The budget is separate from the attribution budgets: 500 replicates put
    the 2.5th percentile on the 13th draw, so the percentile limits of the two
    components want more draws than a permutation tail would.
    """
    enabled: bool = True
    n_boot: int = 2000
    q_fdr: float = 0.05
    min_shared_support: float = 0.9
    label_alpha: float = 0.05

    def validate(self) -> None:
        if self.n_boot < 0:
            raise ConfigError("surveillance.n_boot must be >= 0")
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
    nothing to read without one. The moment estimate and its approximate F
    interval cost one fit per agent and have no budget.

    ``end_wells_censored`` states what a reading on the lowest or highest
    tested well means. Treating it as censored is the coarsened-at-random
    reading of Heitjan and Rubin; treating it as exact is the alternative.

    ``calibrated_interval`` adds the validated population interval, the
    null-wise parametric-bootstrap likelihood-ratio inversion with
    ``calibration_n_boot`` draws per tested value. (A run's readings are
    always dilution intervals; the exact generalized F pivot for exact
    readings is reached through ``amr-clonalshare-mic --method exact`` and
    the Python API.) It is the slow step;
    ``workers`` sets the processes used (0 means every available CPU) and
    does not change the result. A screening run sets ``screening_n_boot``
    (fewer draws, coarser p-values) for every agent except those listed in
    ``calibration_agents``, which keep ``calibration_n_boot``; the draws used
    are recorded with every interval.
    """
    enabled: bool = True
    end_wells_censored: bool = True
    calibrated_interval: bool = True
    calibration_n_boot: int = 199
    workers: int = 0
    screening_n_boot: Optional[int] = None
    calibration_agents: Tuple[str, ...] = ()

    def draws_for(self, agent: str) -> int:
        """The bootstrap draws the calibrated interval of ``agent`` uses."""
        if self.screening_n_boot is None or str(agent) in self.calibration_agents:
            return self.calibration_n_boot
        return self.screening_n_boot

    def validate(self) -> None:
        if self.calibration_n_boot < 19:
            raise ConfigError("censored.calibration_n_boot must be >= 19 (199 was validated)")
        if self.screening_n_boot is not None and self.screening_n_boot < 19:
            raise ConfigError("censored.screening_n_boot must be >= 19")
        if self.calibration_agents and self.screening_n_boot is None:
            raise ConfigError("censored.calibration_agents needs censored.screening_n_boot")
        if self.workers < 0:
            raise ConfigError("censored.workers must be >= 0 (0 uses every available CPU)")


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
class PopulationProbitConfig:
    """Separate population methods, including an opt-in general computation route."""
    enabled: bool = False
    interval_method: str = "fixed_cutoff"
    memory_budget_mb: int = 512
    table_cache_mb: int = 128
    cache_dir: Optional[str] = None

    def validate(self) -> None:
        if not isinstance(self.enabled, bool):
            raise ConfigError("population_probit.enabled must be a boolean")
        if not isinstance(self.interval_method, str) or self.interval_method not in ("fixed_cutoff", "general"):
            raise ConfigError("population_probit.interval_method must be fixed_cutoff or general")
        for name in ("memory_budget_mb", "table_cache_mb"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ConfigError(f"population_probit.{name} must be a positive integer")
        if self.table_cache_mb >= self.memory_budget_mb:
            raise ConfigError("population_probit.table_cache_mb must be smaller than memory_budget_mb")
        if self.cache_dir is not None and (not isinstance(self.cache_dir, str) or not self.cache_dir.strip()):
            raise ConfigError("population_probit.cache_dir must be a nonempty path string or null")


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
    population_probit: PopulationProbitConfig = field(default_factory=PopulationProbitConfig)

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
        self.population_probit.validate()
        if check_files_exist:
            missing = [f"{name} -> {p}" for name, p in (
                ("metadata", self.metadata_path),
                ("phenotype", self.phenotype_path),
                ("mic", self.mic_path)) if p is not None and not p.exists()]
            if missing:
                raise ConfigError(
                    "configured input files not found:\n  " + "\n  ".join(missing))
        return self


def _covariate_columns(ds) -> Tuple[str, ...]:
    """``mic_covariate_columns`` (a list) and ``mic_covariate_column`` (one name), in that order."""
    columns = ds.get("mic_covariate_columns") or ()
    if isinstance(columns, str) or not isinstance(columns, (list, tuple)):
        raise ConfigError("dataset.mic_covariate_columns must be a list of column names")
    single = ds.get("mic_covariate_column")
    names = [str(c) for c in columns] + ([str(single)] if single else [])
    if any(not c.strip() for c in names) or len(set(names)) != len(names):
        raise ConfigError("dataset.mic_covariate_columns must name distinct, nonempty columns")
    return tuple(names)


def _sub(raw: Dict[str, Any], key: str) -> Dict[str, Any]:
    val = raw.get(key, {})
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
                   "evidence", "population_probit"}),
    "dataset": frozenset({
        "name", "strain_id_column", "data_dir", "metadata", "lineage_column",
        "contrast_column", "contrast_levels", "batch_column", "phenotype",
        "phenotype_id_column", "phenotype_antibiotic_column",
        "phenotype_call_column", "phenotype_agent_columns",
        "phenotype_intermediate", "mic",
        "mic_id_column", "mic_antibiotic_column", "mic_value_column",
        "mic_operator_column", "mic_panel_column", "mic_covariate_column", "mic_covariate_columns", "mic_wells", "mic_units", "phenotype_kind", "stratify_by",
        "phenotype_positive_definition", "phenotype_source", "ast_standard", "ast_version", "duplicate_policy"}),
    "attribution": frozenset({"enabled", "folds", "repeats", "n_boot",
                              "n_perm"}),
    "surveillance": frozenset({"enabled", "n_boot", "q_fdr",
                               "min_shared_support", "label_alpha"}),
    "censored": frozenset({"enabled", "end_wells_censored", "calibrated_interval",
                           "calibration_n_boot", "workers", "screening_n_boot",
                           "calibration_agents"}),
    "evidence": frozenset({"enabled", "folds", "repeats", "alpha"}),
    "population_probit": frozenset({"enabled", "interval_method", "memory_budget_mb", "table_cache_mb", "cache_dir"}),
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


def _names(key, value) -> tuple:
    """A list of labels; a bare string is refused rather than read letter by
    letter."""
    if value is None:
        return ()
    if isinstance(value, (str, bytes)) or not isinstance(value, (list, tuple)):
        raise ConfigError(f"{key} must be a list of labels, got {value!r}")
    return tuple(str(v) for v in value)


def _levels(value) -> tuple:
    """``contrast_levels`` as a list of two labels."""
    return _names("dataset.contrast_levels", value)


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
    if not math.isfinite(number) or number != int(number):
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


def _beside(path, config_path: Optional[Path]):
    """A relative path read against the configuration file, as data_dir is."""
    if not isinstance(path, str) or not path.strip() or config_path is None or Path(path).is_absolute():
        return path
    return str((Path(config_path).parent / path).resolve())


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
        phenotype_agent_columns=_names("dataset.phenotype_agent_columns",
                                       ds.get("phenotype_agent_columns")),
        phenotype_id_column=ds.get("phenotype_id_column", "Strain_ID"),
        phenotype_antibiotic_column=ds.get("phenotype_antibiotic_column", "antibiotic"),
        phenotype_call_column=ds.get("phenotype_call_column", "resistant_phenotype"),
        phenotype_intermediate=ds.get("phenotype_intermediate"),
        phenotype_kind=ds.get("phenotype_kind", "undeclared"),
        phenotype_positive_definition=ds.get("phenotype_positive_definition"),
        phenotype_source=ds.get("phenotype_source"),
        ast_standard=ds.get("ast_standard"),
        ast_version=ds.get("ast_version"),
        duplicate_policy=ds.get("duplicate_policy", "error"),
        mic=ds.get("mic"),
        mic_id_column=ds.get("mic_id_column", "Strain_ID"),
        mic_antibiotic_column=ds.get("mic_antibiotic_column", "antibiotic"),
        mic_value_column=ds.get("mic_value_column", "measurement"),
        mic_operator_column=ds.get("mic_operator_column"),
        mic_panel_column=ds.get("mic_panel_column"),
        mic_covariate_columns=_covariate_columns(ds),
        stratify_by=ds.get("stratify_by"),
        mic_wells=ds.get("mic_wells", {}),
        mic_units=ds.get("mic_units", {}),
    )

    attr_raw = _sub(raw, "attribution")
    attribution = AttributionConfig(
        enabled=_bool("attribution", attr_raw, "enabled", True),
        folds=_int("attribution", attr_raw, "folds", 5),
        repeats=_int("attribution", attr_raw, "repeats", 20),
        n_boot=_int("attribution", attr_raw, "n_boot", 400),
        n_perm=_int("attribution", attr_raw, "n_perm", 999),
    )

    surv_raw = _sub(raw, "surveillance")
    surveillance = SurveillanceConfig(
        enabled=_bool("surveillance", surv_raw, "enabled", True),
        n_boot=_int("surveillance", surv_raw, "n_boot", 2000),
        q_fdr=_float("surveillance", surv_raw, "q_fdr", 0.05),
        min_shared_support=_float("surveillance", surv_raw, "min_shared_support", 0.9),
        label_alpha=_float("surveillance", surv_raw, "label_alpha", 0.05),
    )

    cen_raw = _sub(raw, "censored")
    censored = CensoredConfig(
        enabled=_bool("censored", cen_raw, "enabled", True),
        end_wells_censored=_bool("censored", cen_raw, "end_wells_censored", True),
        calibrated_interval=_bool("censored", cen_raw, "calibrated_interval", True),
        calibration_n_boot=_int("censored", cen_raw, "calibration_n_boot", 199),
        workers=_int("censored", cen_raw, "workers", 0),
        screening_n_boot=(None if cen_raw.get("screening_n_boot") is None
                          else _int("censored", cen_raw, "screening_n_boot", 99)),
        calibration_agents=_names("censored.calibration_agents", cen_raw.get("calibration_agents")),
    )

    ev_raw = _sub(raw, "evidence")
    evidence = EvidenceConfig(
        enabled=_bool("evidence", ev_raw, "enabled", True),
        folds=_int("evidence", ev_raw, "folds", 5),
        repeats=_int("evidence", ev_raw, "repeats", 20),
        alpha=_float("evidence", ev_raw, "alpha", 0.05),
    )

    profile_raw = _sub(raw, "population_probit")
    population_probit = PopulationProbitConfig(
        enabled=_bool("population_probit", profile_raw, "enabled", False),
        interval_method=profile_raw.get("interval_method", "fixed_cutoff"),
        memory_budget_mb=_int("population_probit", profile_raw, "memory_budget_mb", 512),
        table_cache_mb=_int("population_probit", profile_raw, "table_cache_mb", 128),
        cache_dir=_beside(profile_raw.get("cache_dir"), config_path))

    return Config(
        dataset=dataset,
        attribution=attribution,
        surveillance=surveillance,
        censored=censored,
        evidence=evidence,
        population_probit=population_probit,
        config_path=Path(config_path).resolve() if config_path else None,
    ).validate(check_files_exist=False)


class _UniqueKeyLoader(yaml.SafeLoader):
    """Safe YAML with duplicate keys rejected instead of silently overwritten."""


def _unique_mapping(loader, node, deep=False):
    mapping = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        try:
            duplicate = key in mapping
        except TypeError as exc:
            raise ConfigError("config mapping keys must be scalar values") from exc
        if duplicate:
            raise ConfigError(f"duplicate config key {key!r} at line {key_node.start_mark.line + 1}")
        mapping[key] = loader.construct_object(value_node, deep=deep)
    return mapping


_UniqueKeyLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _unique_mapping)


def load_config(path: "str | Path", check_files_exist: bool = True) -> Config:
    """Load and validate a YAML config file. Returns a frozen, validated Config."""
    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise ConfigError(f"config file not found: {p}")
    try:
        with p.open("r", encoding="utf-8") as fh:
            raw = yaml.load(fh, Loader=_UniqueKeyLoader)
    except yaml.YAMLError as exc:
        raise ConfigError(f"config file {p} is not valid YAML: {exc}") from exc
    except UnicodeDecodeError as exc:
        raise ConfigError(f"config file {p} is not UTF-8 text (byte {exc.start}: {exc.reason})") from exc
    except OSError as exc:
        raise ConfigError(f"config file {p} cannot be read: {exc.strerror or exc}") from exc
    if not isinstance(raw, dict):
        raise ConfigError(f"config root must be a mapping, got {type(raw).__name__}")
    return from_dict(raw, config_path=p).validate(check_files_exist=check_files_exist)
