# Python API

## Run a complete workflow

```python
from amr_clonalshare import load_config, run

cfg = load_config("examples/workflows/calls.yaml")
record = run(cfg, results_dir="out/calls_api", seed=20260913)
```

Substitute `mic.yaml` or `contrasts.yaml` to run the other recipes. The command line and Python use the same configuration and scientific path. `run` also accepts `overwrite=False`, `progress=None` and `check_files_exist=True`. Use a new directory by default; a deliberate overwrite replaces the earlier results files. The progress callback is optional and should not modify analytical state.

## Read saved results

```python
from amr_clonalshare.outputs import read_result

record = read_result("out/calls_api/clonal_share_result.json")
```

The canonical reader accepts records of schema 1.0 and refuses any other. A record stores collection bounds and method results; inspect statuses and retained collections before using a number. See [result interpretation](manual/05-results.md).

## Exact bounds in a recorded frame

```python
from amr_clonalshare.missingness import finite_collection_bounds

bounds = finite_collection_bounds([1, 0, None], total_count=5)
```

This describes the range allowed by missing binary outcomes in the stated finite frame. It is not a confidence interval, does not impute outcomes and does not establish a population parameter. Use only a defensible recorded-frame denominator.

## Numerical functions

`attribution.clonal_share` and `attribution.layer_clonal_share` return a `ShareResult`: the lineage share of one trait among the represented lineages (`kappa_adj`), its 95% interval for those lineages (`observed_low`, `observed_high`) with the standard error it is studentized by (`observed_se`), the permuted control (`null_mean`) and the permutation p-value. `mic_order.mic_order_share` returns the same fields for the within-stratum MIC ordering as a record, together with the bounds of `latent_order.order_bounds` on the lineage share of the latent ordering; `latent_order.call_bounds` gives those bounds for a binary call. `clonality.decompose_prevalence_difference` and `clonality.decompose_panel` split a prevalence difference between two collections, `comparison.compare_lineage_definitions` compares two lineage definitions on matched records, and the functions of `evalues` give fixed-look and sequential e-values and the panel rules that read them. Each has its own target, stated with it; none is a substitute for another.

::: amr_clonalshare.core.run

::: amr_clonalshare.config.load_config

::: amr_clonalshare.missingness.finite_collection_bounds

::: amr_clonalshare.missingness.difference_bounds

::: amr_clonalshare.attribution.clonal_share

::: amr_clonalshare.attribution.layer_clonal_share

::: amr_clonalshare.attribution.ShareResult

::: amr_clonalshare.mic_order.mic_order_share

::: amr_clonalshare.mic_order.order_scores

::: amr_clonalshare.latent_order.order_bounds

::: amr_clonalshare.latent_order.call_bounds

::: amr_clonalshare.censored.intervals_from_mic

::: amr_clonalshare.censored.intervals_by_panel

::: amr_clonalshare.censored.panel_geometry

::: amr_clonalshare.clonality.decompose_prevalence_difference

::: amr_clonalshare.clonality.decompose_panel

::: amr_clonalshare.comparison.compare_lineage_definitions

::: amr_clonalshare.evalues.e_process

::: amr_clonalshare.evalues.sequential_e_process

::: amr_clonalshare.evalues.e_bh

::: amr_clonalshare.evalues.anytime_bonferroni

::: amr_clonalshare.qc.input_qc

::: amr_clonalshare.outputs.read_result
