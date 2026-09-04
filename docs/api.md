# API reference

The public surface is the eleven modules listed in `amr_clonalshare.__all__`,
each loaded on first use. Everything below is read from the docstrings in the
source, so the reference and the code cannot disagree. Names that begin with
an underscore are internal and are not documented here.

## Estimators

### attribution

::: amr_clonalshare.attribution
    options:
      members:
        - clonal_share
        - layer_clonal_share
        - attribute_partition
        - ShareResult
        - AttributionResult
        - SUPPORT_THRESHOLD

### realised

::: amr_clonalshare.realised
    options:
      members:
        - realised_share
        - realised_interval
        - superpopulation_interval
        - RealisedShare

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
        - e_bh

## Partition diagnostics

### inference

::: amr_clonalshare.inference
    options:
      members:
        - feature_split_test
        - merge_pvalues
        - merge_order_sensitivity
        - correlation_blocks
        - continuum_null_test
        - select_latent_dimension
        - poisson_thin
        - nb_thin
        - binomial_thin
        - thinning_dependence

### influence

::: amr_clonalshare.influence

### lineage

::: amr_clonalshare.lineage

### tva

::: amr_clonalshare.tva

### archephy

::: amr_clonalshare.archephy
    options:
      members:
        - load_tree
        - fitch_change_edges
        - fitch_asr
        - joint_homoplasy_pvalue
        - archephy_cs_test

### baselines

::: amr_clonalshare.baselines

## Running a configuration

::: amr_clonalshare.core
    options:
      members:
        - run
        - validate

::: amr_clonalshare.config
    options:
      members:
        - load_config
        - Config
        - ConfigError
