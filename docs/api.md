# API reference

The public surface is the eight modules listed in `amr_clonalshare.__all__`
together with the two that run a configuration, ten blocks in all, each loaded
on first use. Everything below is read from the docstrings in the
source, so the reference and the code cannot disagree. Names that begin with
an underscore are internal and are not documented here.

## Estimators

### attribution

::: amr_clonalshare.attribution
    options:
      members:
        - clonal_share
        - layer_clonal_share
        - ShareResult
        - SUPPORT_THRESHOLD
        - FEW_LINEAGES

### realised

::: amr_clonalshare.realised
    options:
      members:
        - realised_share
        - realised_interval
        - superpopulation_interval
        - RealisedShare
        - KURTOSIS_LIMIT
        - MIN_GROUPS

### latent

::: amr_clonalshare.latent
    options:
      members:
        - latent_share
        - observed_share
        - latent_interval
        - bivariate_normal_cdf_equal

### censored

::: amr_clonalshare.censored
    options:
      members:
        - censored_clonal_share
        - profile_interval
        - sensitivity_endpoints
        - intervals_from_mic
        - intervals_from_binary
        - panel_geometry
        - scale_is_identified
        - marginal_loglik
        - PanelGeometry
        - CensoredShare

### clonality

::: amr_clonalshare.clonality
    options:
      members:
        - decompose_prevalence_difference
        - decompose_panel
        - lineage_resolved_prevalence
        - trait_concentration

### evalues

::: amr_clonalshare.evalues
    options:
      members:
        - e_process
        - sequential_e_process
        - e_bh
        - e_bh_log
        - combine_independent
        - combine_within_cohort
        - EResult
        - SequentialEResult

### phenotype

::: amr_clonalshare.phenotype
    options:
      members:
        - to_non_susceptible

### stats

::: amr_clonalshare.stats
    options:
      members:
        - permutation_pvalue
        - benjamini_hochberg
        - fisher_exact_p
        - effective_dimension

## Running a configuration

::: amr_clonalshare.core
    options:
      members:
        - run

::: amr_clonalshare.config
    options:
      members:
        - load_config
        - Config
        - DatasetConfig
        - AttributionConfig
        - SurveillanceConfig
        - CensoredConfig
        - EvidenceConfig
        - ConfigError
