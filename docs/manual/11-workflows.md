# Four workflows using one small fictional collection

These CSV files contain 96 fictional isolates, 12 arbitrary lineage labels and one invented agent. They demonstrate formats and software behavior; they are not biological observations, clinical recommendations or performance evidence. No Python programming is needed for the command-line routes. Run from the repository root after installing the package.

## 1 Calls and lineage labels

```bash
amr-clonalshare --config examples/workflows/calls.yaml --check-input
amr-clonalshare --config examples/workflows/calls.yaml --results-dir out/calls
```

Read `data/calls.csv` with `data/metadata.csv`. The declared meaning is R only under `clinical_sir`. Inspect input QC, the observed-scale attribution, support, uncertainty and any method refusal. Replace the fictional source declaration with your actual source; add AST standard/version only when known.

## 2 Recorded MICs

```bash
amr-clonalshare --config examples/workflows/mic.yaml --check-input
amr-clonalshare --config examples/workflows/mic.yaml --results-dir out/mic
```

This route reads `data/mic.csv` without a calls table. The configuration records agent-specific wells and units. End-well interpretation and the sensitivity result must accompany the dilution-scale output. MIC-only identifiers belong to the recorded frame. The method target differs from the binary attribution share.

## 3 Compare two recorded collections

```bash
amr-clonalshare --config examples/workflows/contrasts.yaml --check-input
amr-clonalshare --config examples/workflows/contrasts.yaml --results-dir out/contrasts
```

Metadata declares `earlier` and `later`; the difference is later minus earlier. Read composition and within-lineage terms with shared support and scope. With missing labels, a supported result describes the typed subset and cannot establish collection generalization. The fixture has complete labels; it does not validate performance under missingness.

## 4 Optional Gaussian population model

```bash
amr-clonalshare --config examples/workflows/population.yaml --check-input
amr-clonalshare --config examples/workflows/population.yaml --results-dir out/population
```

This configuration explicitly selects `interval_method: general`. The population module is otherwise disabled and its enabled default remains `fixed_cutoff`. Read the distinct latent-liability ICC target, computation status and model assumptions. An incomplete computation is not a completed interval. General calculations can take longer than the other demonstrations; the manual explains resources and restart behavior.

## The same four routes in Python

```python
from amr_clonalshare import load_config, run

for recipe in ("calls", "mic", "contrasts", "population"):
    config = load_config(f"examples/workflows/{recipe}.yaml")
    result = run(config, results_dir=f"out/{recipe}_api", seed=20260913)
```

## Preserve and interpret the bundle

Inspect the completion manifest, canonical JSON and method statuses before using the CSV or offline reports. Use new destination names for reruns. Explicit `--overwrite` replaces the earlier results files. Demonstration resampling budgets are intentionally small; they are unsuitable for precision claims. Numerical validation campaigns use their separately recorded protocols and seeds. File-format examples do not replace provenance, sampling assessment or method-specific validation.
