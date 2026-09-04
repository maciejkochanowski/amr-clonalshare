"""config.py — declarative configuration for amr-clonalshare.

Everything that changes the science is a config key, and every default is
stated here rather than buried in a function signature. Sections:

``dataset``       cohort identity, strain-ID column, alignment policy, optional
                  metadata file with a lineage column
``files``         one entry per input layer; ``kind`` selects the parser and the
                  distance treatment (``wide_binary`` vs ``one_hot``)
``distance``      binary coefficient and, explicitly, what an empty union means
``snf``           fusion hyperparameters (K, mu, T, alpha, tie policy)
``trait_cluster`` gating, feature grouping, k range and k selector, thresholds
``tva``           multi-split thinning inference
``influence``     layer-influence / fusion-collapse diagnostics
``validation``    acceptance thresholds for the synthetic harness

``data_dir`` is resolved relative to the config file, so a config can sit next
to its data. Validation is eager: :func:`load_config` raises :class:`ConfigError`
with an actionable message instead of failing deep inside the analysis.
"""
from __future__ import annotations

import difflib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "PyYAML is required to load configs. Install with `pip install pyyaml` "
        "or `pip install -e .` from the package root."
    ) from exc


class ConfigError(ValueError):
    """Raised when a config is structurally invalid or internally inconsistent."""


# --------------------------------------------------------------------------- #
#  Sections
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class FileSpec:
    """One input layer.

    ``kind``:
      * ``wide_binary`` - independent 0/1 presence calls; Jaccard-family distance.
      * ``one_hot`` - a categorical variable expanded to indicator columns
        (capsule K locus, O antigen, MLST). Exactly one column is 1 per row.
        These are validated at load time, excluded from count aggregation for
        thinning (a one-hot layer aggregates to a constant, not a count) and
        given a matching-coefficient distance rather than Jaccard, because
        Jaccard on a one-hot block takes only a handful of distinct values.
    """
    path: str
    kind: str = "wide_binary"
    value_column: Optional[str] = None
    encoding: Optional[Dict[str, str]] = None
    strip_columns: bool = False
    groups: Optional[Dict[str, List[str]]] = None   # locus-level grouping
    protected: Sequence[str] = ()                   # never dropped by the gate

    _KINDS = ("wide_binary", "one_hot", "wide_phenotype", "series", "long", "newick")
    _MATRIX_KINDS = ("wide_binary", "one_hot")

    def validate(self, role: str) -> None:
        if self.kind not in self._KINDS:
            raise ConfigError(
                f"files.{role}.kind = {self.kind!r} is not one of {self._KINDS}")
        if self.kind in ("series", "long") and not self.value_column:
            raise ConfigError(
                f"files.{role}.value_column is required for kind={self.kind!r}")


@dataclass(frozen=True)
class DatasetConfig:
    name: str
    strain_id_column: str = "Strain_ID"
    data_dir: str = "."
    expected_n: Optional[int] = None
    strain_alignment_policy: str = "intersect_core"
    metadata: Optional[str] = None          # path to a metadata CSV
    lineage_column: Optional[str] = None    # ST / clonal group column therein
    external_columns: Sequence[str] = ()    # columns to score the partition against
    # Two collections to compare when decomposing a prevalence difference into
    # a lineage-mix component and a within-lineage rate component. A period, a
    # country, a host. Left unset, the decomposition is simply not run.
    contrast_column: Optional[str] = None
    contrast_levels: Sequence[str] = ()
    #: Optional long-format CSV of measured antimicrobial susceptibility, one
    #: row per (isolate, antibiotic). This is the only truly external criterion
    #: available for a genotype panel: it is produced by a different instrument
    #: and is sensitive to mechanisms the acquired-gene panel cannot see. Not
    #: used for clustering; see :mod:`amr_clonalshare.phenotype`.
    phenotype: Optional[str] = None
    phenotype_id_column: str = "Strain_ID"
    phenotype_antibiotic_column: str = "antibiotic"
    phenotype_call_column: str = "resistant_phenotype"
    phenotype_intermediate: str = "non_susceptible"
    #: Metadata column identifying the depositing collection, for the
    #: leave-one-collection-out check on the phenotype comparison.
    phenotype_stratum_column: str = "Country"
    #: Optional long-format CSV of measured minimum inhibitory
    #: concentrations, one row per (isolate, antimicrobial). A dichotomised
    #: call keeps only which side of a cut-off an isolate fell; the recorded
    #: dilution keeps where in the panel it fell, and the difference is
    #: measurable rather than rhetorical. See
    #: :mod:`amr_clonalshare.censored`. The identifier column must already
    #: match the strain identifiers of the layers: this loader joins, it does
    #: not repair identifiers, because a rule that silently rewrote them would
    #: be a rule that silently joined the wrong isolates.
    mic: Optional[str] = None
    mic_id_column: str = "Strain_ID"
    mic_antibiotic_column: str = "antibiotic"
    mic_value_column: str = "measurement"
    #: Column holding a recorded censoring operator (``<``, ``<=``, ``>``,
    #: ``>=``). Where it is present it decides; where it is absent the panel
    #: geometry decides, and that is an assumption the run reports.
    mic_operator_column: Optional[str] = None
    #: What to do with a cell of a binary layer that holds no value: an empty
    #: field in the CSV, or an isolate absent from one layer under the
    #: ``union`` alignment policy. A binary layer has no state for "not
    #: measured", and writing 0 there would turn absence of a measurement into
    #: absence of the trait. ``refuse`` (default) stops the run with the
    #: positions listed; ``drop_rows`` removes every isolate with a missing
    #: cell; ``drop_columns`` removes every feature with a missing cell. Each
    #: choice is reported in the input QC record. Missing values are never
    #: imputed.
    missing_policy: str = "refuse"

    _POLICIES = ("intersect_core", "union", "strict_n")
    _MISSING_POLICIES = ("refuse", "drop_rows", "drop_columns")

    def validate(self) -> None:
        if self.strain_alignment_policy not in self._POLICIES:
            raise ConfigError(
                f"dataset.strain_alignment_policy = {self.strain_alignment_policy!r} "
                f"not one of {self._POLICIES}")
        if self.missing_policy not in self._MISSING_POLICIES:
            raise ConfigError(
                f"dataset.missing_policy = {self.missing_policy!r} not one of "
                f"{self._MISSING_POLICIES}")
        if self.strain_alignment_policy == "strict_n" and self.expected_n is None:
            raise ConfigError(
                "dataset.expected_n must be set when strain_alignment_policy='strict_n'")
        if self.lineage_column and not self.metadata:
            raise ConfigError(
                "dataset.lineage_column requires dataset.metadata to be set")
        if self.contrast_column:
            if not self.metadata:
                raise ConfigError(
                    "dataset.contrast_column requires dataset.metadata to be set")
            if not self.lineage_column:
                raise ConfigError(
                    "dataset.contrast_column requires dataset.lineage_column: "
                    "the decomposition splits a difference by lineage")
            if len(self.contrast_levels) != 2:
                raise ConfigError(
                    "dataset.contrast_levels must name exactly two values of "
                    f"{self.contrast_column!r}; got "
                    f"{list(self.contrast_levels)!r}")
        if self.mic and not self.lineage_column:
            raise ConfigError(
                "dataset.mic requires dataset.lineage_column: the censored "
                "share divides the latent concentration between lineage and "
                "everything else, and without a lineage there is nothing to "
                "divide it by")
        if self.phenotype_intermediate not in ("non_susceptible", "drop",
                                               "susceptible"):
            raise ConfigError(
                "dataset.phenotype_intermediate must be one of "
                "'non_susceptible', 'drop', 'susceptible'; got "
                f"{self.phenotype_intermediate!r}")


@dataclass(frozen=True)
class DistanceConfig:
    """Binary distance coefficient and the empty-union convention.

    ``undefined_pair`` decides what two isolates carrying none of a layer's
    features mean. ``"identical"`` is the classical Jaccard convention and
    merges them into one block; on sparse accessory-genome layers that block can
    be the majority of the cohort and can become the largest reported
    "archetype". ``"distinct"`` refuses to merge them. There is no neutral
    choice, so there is no silent default: the value used is echoed into the
    result JSON and both must be reported.
    """
    metric: str = "jaccard"
    undefined_pair: str = "identical"
    one_hot_metric: str = "simple_matching"

    _METRICS = ("jaccard", "dice", "simple_matching", "hamming")
    _UNDEFINED = ("identical", "distinct", "nan")

    def validate(self) -> None:
        for name, val, allowed in (("metric", self.metric, self._METRICS),
                                   ("one_hot_metric", self.one_hot_metric, self._METRICS),
                                   ("undefined_pair", self.undefined_pair, self._UNDEFINED)):
            if val not in allowed:
                raise ConfigError(f"distance.{name} = {val!r} not one of {allowed}")


@dataclass(frozen=True)
class SNFConfig:
    """Similarity-network-fusion hyperparameters (Wang et al. 2014).

    Defaults follow the paper and the reference implementations: ``K`` in
    [10, 30], ``mu`` in [0.3, 0.8], ``T = 20``, self-loop ``alpha = 0.5``.
    ``K = None`` uses ``min(ceil(n/5), 30)``.

    ``tie_policy="inclusive"`` keeps every neighbour tied with the K-th largest
    affinity, which makes the fusion invariant to input row order. The
    ``"strict"`` policy reproduces the order-dependent ``argsort`` behaviour and
    exists only for the regression test that demonstrates the difference.

    ``update="renormalise"`` is the cross-diffusion update of SNFtool >= 2.2.1
    (the current CRAN release); ``update="add_identity"`` is the 2014 update of
    Wang's MATLAB code, SNFtool <= 2.2 and snfpy. ``benchmarks/
    snf_update_benchmark.py`` shows the two agree on every control and on the
    case study; the option exists so that a user reproducing snfpy can.
    """
    K: Optional[int] = None
    mu: float = 0.5
    T: int = 20
    alpha: float = 0.5
    tie_policy: str = "inclusive"
    update: str = "renormalise"

    def validate(self) -> None:
        if not (0 < self.mu):
            raise ConfigError(f"snf.mu must be positive; got {self.mu}")
        if self.T < 1:
            raise ConfigError(f"snf.T must be >= 1; got {self.T}")
        # Under "renormalise", alpha is the diagonal of a row-stochastic matrix
        # and must be a proper fraction; under "add_identity" it is the constant
        # added to the diagonal each iteration, whose reference value is 1.0.
        if self.update == "renormalise" and not (0.0 <= self.alpha < 1.0):
            raise ConfigError(f"snf.alpha must lie in [0, 1); got {self.alpha}")
        if self.alpha < 0.0:
            raise ConfigError(f"snf.alpha must be non-negative; got {self.alpha}")
        if self.tie_policy not in ("inclusive", "strict"):
            raise ConfigError(
                f"snf.tie_policy must be 'inclusive' or 'strict'; got {self.tie_policy!r}")
        if self.update not in ("renormalise", "add_identity"):
            raise ConfigError(
                "snf.update must be 'renormalise' or 'add_identity'; "
                f"got {self.update!r}")
        if self.K is not None and self.K < 1:
            raise ConfigError(f"snf.K must be >= 1 or null; got {self.K}")


@dataclass(frozen=True)
class PrevalenceGate:
    """Which features enter the distance computation.

    ``lo``/``hi`` are prevalence fractions; ``min_count`` is an absolute
    alternative (a feature is kept if it is present in at least this many
    isolates). The default ``lo = 0.0`` keeps rare determinants: in AMR
    surveillance the emergent low-prevalence resistance is the signal, and a
    2% floor deletes exactly the colistin/tigecycline/carbapenemase columns the
    analysis exists to find. Gated features are always reported, and any
    feature named in a layer's ``protected`` list is never gated out.
    """
    lo: float = 0.0
    hi: float = 1.0
    min_count: int = 2

    def validate(self, where: str) -> None:
        if not (0.0 <= self.lo < self.hi <= 1.0):
            raise ConfigError(
                f"{where}: require 0 <= lo < hi <= 1, got lo={self.lo}, hi={self.hi}")
        if self.min_count < 0:
            raise ConfigError(f"{where}: min_count must be >= 0")


@dataclass(frozen=True)
class DefiningConfig:
    """Thresholds for calling a feature 'defining' in the descriptive profile."""
    min_abs_cohens_h: float = 0.5
    q_fdr: float = 0.05
    min_support: int = 3

    def validate(self) -> None:
        if self.min_abs_cohens_h < 0:
            raise ConfigError("defining.min_abs_cohens_h must be >= 0")
        if not (0 < self.q_fdr < 1):
            raise ConfigError("defining.q_fdr must lie in (0, 1)")


@dataclass(frozen=True)
class StabilityConfig:
    """Consensus-confidence cutoffs for {robust, intermediate, fragile}."""
    tau_robust: float = 0.4
    tau_fragile: float = 0.1

    def validate(self) -> None:
        if not (self.tau_fragile < self.tau_robust):
            raise ConfigError("stability.tau_fragile must be < stability.tau_robust")


@dataclass(frozen=True)
class TraitClusterConfig:
    prevalence_gate: PrevalenceGate = field(default_factory=PrevalenceGate)
    layers: List[str] = field(default_factory=lambda: ["amr", "vir"])
    k_range: List[int] = field(default_factory=lambda: [2, 3, 4, 5, 6])
    k_select_primary: str = "mdl"
    prediction_strength_threshold: float = 0.8
    top_variance: Optional[int] = None
    collapse_feature_groups: bool = False
    defining: DefiningConfig = field(default_factory=DefiningConfig)
    stability: StabilityConfig = field(default_factory=StabilityConfig)

    _K_SELECT = ("mdl", "prediction_strength", "gap", "bic_mixture")

    def validate(self) -> None:
        self.prevalence_gate.validate("trait_cluster.prevalence_gate")
        self.defining.validate()
        self.stability.validate()
        if self.k_select_primary not in self._K_SELECT:
            raise ConfigError(
                f"trait_cluster.k_select_primary = {self.k_select_primary!r} "
                f"not one of {self._K_SELECT}")
        if len(self.k_range) < 2 or any(k < 2 for k in self.k_range):
            raise ConfigError("trait_cluster.k_range must list >=2 values, each >=2")
        if self.top_variance is not None and self.top_variance < 1:
            raise ConfigError("trait_cluster.top_variance must be >= 1 or null")


@dataclass(frozen=True)
class TVAConfig:
    """Multi-split data-thinning post-clustering inference."""
    enabled: bool = True
    n_splits: int = 9
    eps: float = 0.5
    dispersion: str = "mom"                 # mom | mle | pooled_mom | pooled_mle
    merge: str = "exchangeable_ruger"       # exchangeable_ruger | ruger | twice_mean
    q_fdr: float = 0.05
    min_distinct: int = 3
    min_columns: int = 3
    n_boot_continuum: int = 99
    # BIC starts at q=1..3 and extends to continuum_q_max only when the current
    # optimum is at the boundary. Four is the largest product Gauss-Hermite grid
    # that the dense implementation supports safely; q=4 itself remains a
    # boundary, not a pass.
    continuum_q_max: int = 4
    #: Features correlating at or above this are forced into the same split
    #: unit. The feature split is only valid if the held-out half is
    #: independent of the discovery half under the null; per-feature splitting
    #: on a correlated panel rejected in 3999 of 4000 replicates (calibration
    #: study, experiment B). Experiment G sweeps this threshold: type-I 0.0163
    #: at 0.5, 0.0195 at 0.7 and 0.7897 at 0.85, so 0.5 and 0.7 both hold the
    #: nominal 0.05 and the level alone does not separate them
    #: (``benchmarks/results_n4000_2026-08-14/``). The default is chosen on the
    #: shipped cohorts instead: over 0.30 to 0.99 the split design is adequate
    #: for 0.30-0.75 on the Klebsiella case study and 0.49-0.72 on the clonal
    #: control, and inside that intersection 0.7 yields the finer units and
    #: more margin to the chaining bound than 0.5
    #: (``benchmarks/results_adequacy_2026-08-14/``). Power is not compared.
    #: Every run additionally reports and enforces ``split_design_adequate`` on
    #: the design actually produced.
    #: Set to 0 to use only the groups declared in `files`.
    block_threshold: float = 0.7

    _DISPERSION = ("mom", "mle", "pooled_mom", "pooled_mle")
    _MERGE = ("exchangeable_ruger", "ruger", "twice_mean")

    def validate(self) -> None:
        if self.n_splits < 1:
            raise ConfigError("tva.n_splits must be >= 1")
        if self.n_boot_continuum < 20:
            raise ConfigError(
                "tva.n_boot_continuum must be >= 20: with the corrected "
                "bootstrap p-value 1/(B+1), fewer replicates cannot attain "
                "the public alpha = 0.05")
        if not (1 <= self.continuum_q_max <= 4):
            raise ConfigError(
                "tva.continuum_q_max must lie in [1, 4]; q=4 is the current "
                "computational safety boundary")
        if not (0.0 <= self.block_threshold <= 1.0):
            raise ConfigError("tva.block_threshold must lie in [0, 1]; "
                              f"got {self.block_threshold}")
        if not (0 < self.eps < 1):
            raise ConfigError("tva.eps must lie in (0, 1)")
        if self.dispersion not in self._DISPERSION:
            raise ConfigError(f"tva.dispersion not one of {self._DISPERSION}")
        if self.merge not in self._MERGE:
            raise ConfigError(f"tva.merge not one of {self._MERGE}")


@dataclass(frozen=True)
class InfluenceConfig:
    """Layer-influence / fusion-collapse diagnostics."""
    enabled: bool = True
    n_perm: int = 0
    collapse_threshold: float = 1.5

    def validate(self) -> None:
        if self.n_perm < 0:
            raise ConfigError("influence.n_perm must be >= 0")


@dataclass(frozen=True)
class AttributionConfig:
    """Budgets and the gate threshold for the lineage attribution.

    Kept apart from ``surveillance`` for the reason that block gives for
    keeping itself apart from ``influence``: these are different resamplings
    and do not want the same count. Every attribution figure is a
    cross-validated statistic averaged over ``repeats`` fold draws, so the cost
    is ``repeats + n_perm + n_boot`` model fits per quantity rather than one,
    and 2000 of each would put a single run into hours for no gain in the third
    decimal place.

    ``lambda_gate`` is the share of a partition's explanatory power that may be
    attributable to the lineage label before the partition stops being readable
    as trait structure. It is a magnitude, and it replaces the significance
    test that preceded it as the gate: the concordance z is still computed and
    still reported, but it grows with the number of pairs at fixed structure
    and so cannot carry a threshold.

    ``refit_in_fold`` rebuilds the partition from the training isolates of
    every fold instead of scoring the one partition the run fitted on all of
    them. It is on by default because without it ``r2_partition`` is out of
    sample over isolates but not over the clustering, so a held-out isolate
    helped draw the cluster it is then scored against, and the share that comes
    back is the wrong side of the question the gate asks. The cost is real and
    is stated here so that switching it off is an informed choice: it is
    ``folds * repeats`` extra clusterings, measured at 94 seconds on the
    shipped *S. suis* cohort and 243 seconds on *Klebsiella*, against 0.6 and
    1.7 seconds for the fixed-partition estimate. It also sets ``n_boot`` aside
    for the run, because the lineage bootstrap duplicates isolates and a
    partition cannot honestly be rebuilt on a cohort holding the same isolate
    twice.
    """
    enabled: bool = True
    folds: int = 5
    repeats: int = 20
    n_boot: int = 400
    n_perm: int = 200
    lambda_gate: float = 0.5
    refit_in_fold: bool = True

    def validate(self) -> None:
        if self.folds < 2:
            raise ConfigError("attribution.folds must be >= 2")
        if self.repeats < 1:
            raise ConfigError("attribution.repeats must be >= 1")
        if self.n_boot < 0 or self.n_perm < 0:
            raise ConfigError("attribution.n_boot and n_perm must be >= 0")
        if not 0 < self.lambda_gate <= 1:
            raise ConfigError("attribution.lambda_gate must be in (0, 1]")


@dataclass(frozen=True)
class SurveillanceConfig:
    """Budgets and thresholds for the lineage-resolved prevalence outputs.

    They are separate from ``influence.n_perm``, which sets the permutation
    budget of a different diagnostic. A percentile interval and a
    permutation tail are not the same resampling and do not want the same
    count: 500 replicates put the 2.5th percentile on the 13th draw, and 500
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
        if self.n_boot < 0 or self.n_perm < 0:
            raise ConfigError("surveillance.n_boot and n_perm must be >= 0")
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
class ValidationConfig:
    cluster_ari_min: float = 0.80
    schema_check: bool = True


@dataclass(frozen=True)
class Config:
    """Fully-resolved, validated configuration for one cohort."""
    dataset: DatasetConfig
    files: Dict[str, FileSpec]
    trait_cluster: TraitClusterConfig
    validation: ValidationConfig
    distance: DistanceConfig = field(default_factory=DistanceConfig)
    snf: SNFConfig = field(default_factory=SNFConfig)
    tva: TVAConfig = field(default_factory=TVAConfig)
    influence: InfluenceConfig = field(default_factory=InfluenceConfig)
    surveillance: SurveillanceConfig = field(
        default_factory=SurveillanceConfig)
    attribution: AttributionConfig = field(
        default_factory=AttributionConfig)
    censored: CensoredConfig = field(default_factory=CensoredConfig)
    evidence: EvidenceConfig = field(default_factory=EvidenceConfig)
    config_path: Optional[Path] = None

    @property
    def data_root(self) -> Path:
        base = self.config_path.parent if self.config_path else Path.cwd()
        return (base / self.dataset.data_dir).resolve()

    def file_path(self, role: str) -> Path:
        if role not in self.files:
            raise ConfigError(f"no file configured for role {role!r}")
        return self.data_root / self.files[role].path

    @property
    def metadata_path(self) -> Optional[Path]:
        if not self.dataset.metadata:
            return None
        return self.data_root / self.dataset.metadata

    @property
    def mic_path(self) -> Optional[Path]:
        if not self.dataset.mic:
            return None
        return self.data_root / self.dataset.mic

    def validate(self, check_files_exist: bool = True) -> "Config":
        self.dataset.validate()
        self.trait_cluster.validate()
        self.distance.validate()
        self.snf.validate()
        self.tva.validate()
        self.influence.validate()
        self.surveillance.validate()
        self.attribution.validate()
        self.censored.validate()
        self.evidence.validate()
        for role, spec in self.files.items():
            spec.validate(role)
        for role in self.trait_cluster.layers:
            if role not in self.files:
                raise ConfigError(
                    f"trait_cluster.layers references {role!r}, which has no files.* entry")
            if self.files[role].kind not in FileSpec._MATRIX_KINDS:
                raise ConfigError(
                    f"trait_cluster.layers references {role!r} of kind "
                    f"{self.files[role].kind!r}; clustering layers must be "
                    f"{FileSpec._MATRIX_KINDS}")
        if check_files_exist:
            missing = [f"{role} -> {self.file_path(role)}"
                       for role in self.files if not self.file_path(role).exists()]
            mp = self.metadata_path
            if mp is not None and not mp.exists():
                missing.append(f"metadata -> {mp}")
            qp = self.mic_path
            if qp is not None and not qp.exists():
                missing.append(f"mic -> {qp}")
            if missing:
                raise ConfigError(
                    "configured data files not found (data_root="
                    f"{self.data_root}):\n  " + "\n  ".join(missing))
        return self


# --------------------------------------------------------------------------- #
#  Loading
# --------------------------------------------------------------------------- #
def _files_from_dict(raw: Dict[str, Any]) -> Dict[str, FileSpec]:
    out: Dict[str, FileSpec] = {}
    for role, spec in (raw or {}).items():
        if not isinstance(spec, dict):
            raise ConfigError(
                f"files.{role} must be a mapping, got {type(spec).__name__}")
        groups = spec.get("groups")
        if groups is not None and not isinstance(groups, dict):
            raise ConfigError(f"files.{role}.groups must be a mapping name -> [columns]")
        out[role] = FileSpec(
            path=spec["path"],
            kind=spec.get("kind", "wide_binary"),
            value_column=spec.get("value_column"),
            encoding=spec.get("encoding"),
            strip_columns=bool(spec.get("strip_columns", False)),
            groups={k: list(v) for k, v in groups.items()} if groups else None,
            protected=tuple(spec.get("protected", ()) or ()),
        )
    return out


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
    "": frozenset({"dataset", "files", "trait_cluster", "distance", "snf",
                   "tva", "influence", "surveillance", "attribution",
                   "censored", "evidence", "validation"}),
    "dataset": frozenset({
        "name", "strain_id_column", "data_dir", "expected_n",
        "strain_alignment_policy", "metadata", "lineage_column",
        "external_columns", "contrast_column", "contrast_levels", "phenotype",
        "phenotype_id_column", "phenotype_antibiotic_column",
        "phenotype_call_column", "phenotype_intermediate",
        "phenotype_stratum_column", "mic", "mic_id_column",
        "mic_antibiotic_column", "mic_value_column", "mic_operator_column",
        "missing_policy"}),
    "trait_cluster": frozenset({
        "prevalence_gate", "defining", "stability", "layers", "k_range",
        "k_select_primary", "prediction_strength_threshold", "top_variance",
        "vir_top_variance", "collapse_feature_groups"}),
    "trait_cluster.prevalence_gate": frozenset({"lo", "hi", "min_count"}),
    "trait_cluster.defining": frozenset({"min_abs_cohens_h", "q_fdr",
                                         "min_support"}),
    "trait_cluster.stability": frozenset({"tau_robust", "tau_fragile"}),
    "distance": frozenset({"metric", "undefined_pair", "one_hot_metric"}),
    "snf": frozenset({"K", "mu", "T", "alpha", "tie_policy", "update"}),
    "tva": frozenset({
        "enabled", "n_splits", "eps", "dispersion", "merge", "q_fdr",
        "min_distinct", "min_columns", "n_boot_continuum", "continuum_q_max",
        "block_threshold"}),
    "influence": frozenset({"enabled", "n_perm", "collapse_threshold"}),
    "attribution": frozenset({"enabled", "folds", "repeats", "n_boot",
                              "n_perm", "lambda_gate", "refit_in_fold"}),
    "surveillance": frozenset({"enabled", "n_boot", "n_perm", "q_fdr",
                               "min_shared_support", "label_alpha"}),
    "censored": frozenset({"enabled", "n_boot", "end_wells_censored",
                           "sensitivity"}),
    "evidence": frozenset({"enabled", "folds", "repeats", "alpha"}),
    "validation": frozenset({"cluster_ari_min", "schema_check"}),
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
        suggestion = f"; did you mean {close[0]!r}?" if close else ""
        raise ConfigError(
            f"unknown key {key!r} in {where}{suggestion} A key that is not "
            f"read is a silently disabled analysis, so it is refused rather "
            f"than ignored. Known keys: {sorted(allowed)}")
    for key in mapping:
        child = f"{section}.{key}" if section else str(key)
        if child in _KNOWN_KEYS:
            _reject_unknown_keys(mapping[key], child)


def from_dict(raw: Dict[str, Any], config_path: Optional[Path] = None) -> Config:
    """Build a validated :class:`Config` from a parsed YAML/JSON mapping."""
    if "dataset" not in raw or "files" not in raw:
        raise ConfigError("config must define at least: dataset, files")
    _reject_unknown_keys(raw, "")

    ds = raw["dataset"]
    dataset = DatasetConfig(
        name=ds["name"],
        strain_id_column=ds.get("strain_id_column", "Strain_ID"),
        data_dir=ds.get("data_dir", "."),
        expected_n=ds.get("expected_n"),
        strain_alignment_policy=ds.get("strain_alignment_policy", "intersect_core"),
        metadata=ds.get("metadata"),
        lineage_column=ds.get("lineage_column"),
        external_columns=tuple(ds.get("external_columns", ()) or ()),
        contrast_column=ds.get("contrast_column"),
        contrast_levels=tuple(str(v) for v in (ds.get("contrast_levels") or ())),
        phenotype=ds.get("phenotype"),
        phenotype_id_column=ds.get("phenotype_id_column", "Strain_ID"),
        phenotype_antibiotic_column=ds.get("phenotype_antibiotic_column", "antibiotic"),
        phenotype_call_column=ds.get("phenotype_call_column", "resistant_phenotype"),
        phenotype_intermediate=ds.get("phenotype_intermediate", "non_susceptible"),
        phenotype_stratum_column=ds.get("phenotype_stratum_column", "Country"),
        mic=ds.get("mic"),
        mic_id_column=ds.get("mic_id_column", "Strain_ID"),
        mic_antibiotic_column=ds.get("mic_antibiotic_column", "antibiotic"),
        mic_value_column=ds.get("mic_value_column", "measurement"),
        mic_operator_column=ds.get("mic_operator_column"),
        missing_policy=ds.get("missing_policy", "refuse"),
    )

    tc = _sub(raw, "trait_cluster")
    g = _sub(tc, "prevalence_gate")
    gate_default = PrevalenceGate()
    gate = PrevalenceGate(
        lo=float(g.get("lo", gate_default.lo)),
        hi=float(g.get("hi", gate_default.hi)),
        min_count=int(g.get("min_count", gate_default.min_count)),
    )
    d = _sub(tc, "defining")
    dd = DefiningConfig(
        min_abs_cohens_h=float(d.get("min_abs_cohens_h", 0.5)),
        q_fdr=float(d.get("q_fdr", 0.05)),
        min_support=int(d.get("min_support", 3)),
    )
    s = _sub(tc, "stability")
    ss = StabilityConfig(
        tau_robust=float(s.get("tau_robust", 0.4)),
        tau_fragile=float(s.get("tau_fragile", 0.1)),
    )
    # `vir_top_variance` is an accepted alias for the same setting.
    top_var = tc.get("top_variance", tc.get("vir_top_variance"))
    trait_cluster = TraitClusterConfig(
        prevalence_gate=gate,
        layers=list(tc.get("layers", ["amr", "vir"])),
        k_range=list(tc.get("k_range", [2, 3, 4, 5, 6])),
        k_select_primary=tc.get("k_select_primary", "mdl"),
        prediction_strength_threshold=float(
            tc.get("prediction_strength_threshold", 0.8)),
        top_variance=None if top_var is None else int(top_var),
        collapse_feature_groups=bool(tc.get("collapse_feature_groups", False)),
        defining=dd,
        stability=ss,
    )

    dist_raw = _sub(raw, "distance")
    distance = DistanceConfig(
        metric=dist_raw.get("metric", "jaccard"),
        undefined_pair=dist_raw.get("undefined_pair", "identical"),
        one_hot_metric=dist_raw.get("one_hot_metric", "simple_matching"),
    )

    snf_raw = _sub(raw, "snf")
    snf = SNFConfig(
        K=snf_raw.get("K"),
        mu=float(snf_raw.get("mu", 0.5)),
        T=int(snf_raw.get("T", 20)),
        alpha=float(snf_raw.get("alpha", 0.5)),
        tie_policy=snf_raw.get("tie_policy", "inclusive"),
        update=snf_raw.get("update", "renormalise"),
    )

    tva_raw = _sub(raw, "tva")
    tva = TVAConfig(
        enabled=bool(tva_raw.get("enabled", True)),
        n_splits=int(tva_raw.get("n_splits", 9)),
        eps=float(tva_raw.get("eps", 0.5)),
        dispersion=tva_raw.get("dispersion", "mom"),
        merge=tva_raw.get("merge", "exchangeable_ruger"),
        q_fdr=float(tva_raw.get("q_fdr", 0.05)),
        min_distinct=int(tva_raw.get("min_distinct", 3)),
        min_columns=int(tva_raw.get("min_columns", 3)),
        n_boot_continuum=int(tva_raw.get("n_boot_continuum", 99)),
        continuum_q_max=int(tva_raw.get("continuum_q_max", 4)),
        block_threshold=float(tva_raw.get("block_threshold", 0.7)),
    )

    inf_raw = _sub(raw, "influence")
    influence = InfluenceConfig(
        enabled=bool(inf_raw.get("enabled", True)),
        n_perm=int(inf_raw.get("n_perm", 0)),
        collapse_threshold=float(inf_raw.get("collapse_threshold", 1.5)),
    )

    surv_raw = _sub(raw, "surveillance")
    surveillance = SurveillanceConfig(
        enabled=bool(surv_raw.get("enabled", True)),
        n_boot=int(surv_raw.get("n_boot", 2000)),
        n_perm=int(surv_raw.get("n_perm", 2000)),
        q_fdr=float(surv_raw.get("q_fdr", 0.05)),
        min_shared_support=float(surv_raw.get("min_shared_support", 0.8)),
        label_alpha=float(surv_raw.get("label_alpha", 0.05)),
    )

    attr_raw = _sub(raw, "attribution")
    attribution = AttributionConfig(
        enabled=bool(attr_raw.get("enabled", True)),
        folds=int(attr_raw.get("folds", 5)),
        repeats=int(attr_raw.get("repeats", 20)),
        n_boot=int(attr_raw.get("n_boot", 400)),
        n_perm=int(attr_raw.get("n_perm", 200)),
        lambda_gate=float(attr_raw.get("lambda_gate", 0.5)),
        refit_in_fold=bool(attr_raw.get("refit_in_fold", True)),
    )

    cen_raw = _sub(raw, "censored")
    censored = CensoredConfig(
        enabled=bool(cen_raw.get("enabled", True)),
        n_boot=int(cen_raw.get("n_boot", 200)),
        end_wells_censored=bool(cen_raw.get("end_wells_censored", True)),
        sensitivity=bool(cen_raw.get("sensitivity", True)),
    )

    ev_raw = _sub(raw, "evidence")
    evidence = EvidenceConfig(
        enabled=bool(ev_raw.get("enabled", True)),
        folds=int(ev_raw.get("folds", 5)),
        repeats=int(ev_raw.get("repeats", 20)),
        alpha=float(ev_raw.get("alpha", 0.05)),
    )

    vv = _sub(raw, "validation")
    validation = ValidationConfig(
        cluster_ari_min=float(vv.get("cluster_ari_min", 0.80)),
        schema_check=bool(vv.get("schema_check", True)),
    )

    return Config(
        dataset=dataset,
        files=_files_from_dict(raw["files"]),
        trait_cluster=trait_cluster,
        validation=validation,
        distance=distance,
        snf=snf,
        tva=tva,
        influence=influence,
        surveillance=surveillance,
        attribution=attribution,
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
