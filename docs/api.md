# Python API

## Run a complete workflow

```python
from amr_clonalshare import load_config, run

cfg = load_config("examples/workflows/calls.yaml")
record = run(cfg, results_dir="out/calls_api", seed=20260913)
```

Substitute `mic.yaml`, `contrasts.yaml` or `population.yaml` to run the other recipes. The command line and Python use the same configuration and scientific path. `run` also accepts `overwrite=False` and `progress=None`. Use a new directory by default; a deliberate overwrite replaces the earlier results files. The progress callback is optional and should not modify analytical state.

## Read saved results

```python
from amr_clonalshare.outputs import read_result

record = read_result("out/calls_api/clonal_share_result.json")
```

The canonical reader accepts records of schema 1.0 and refuses any other. A record stores collection bounds and method results; inspect statuses and retained cohorts before using a number. See [result interpretation](manual/05-results.md).

## Exact bounds in a recorded frame

```python
from amr_clonalshare.missingness import finite_collection_bounds

bounds = finite_collection_bounds([1, 0, None], total_count=5)
```

This describes the range allowed by missing binary outcomes in the stated finite frame. It is not a confidence interval, does not impute outcomes and does not establish a population parameter. Use only a defensible recorded-frame denominator.

## Numerical functions

`attribution.clonal_share`, `realised.realised_share`, `censored.censored_clonal_share`, `clonality.decompose_prevalence_difference` and the evidence functions retain their distinct targets. Optional `population_probit_icc` and `general_probit_icc` target a Gaussian latent-liability population ICC. Requesting the general method does not change the default or validate other models. Research prototypes remain separate from the normal analysis path unless their documented promotion gate is satisfied.

::: amr_clonalshare.core.run

::: amr_clonalshare.config.load_config

::: amr_clonalshare.missingness.finite_collection_bounds

::: amr_clonalshare.censored.censored_clonal_share
