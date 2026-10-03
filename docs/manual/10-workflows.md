# Three workflows using one small fictional collection

These CSV files contain 96 fictional isolates, 12 arbitrary lineage labels and one invented agent. They demonstrate formats and software behavior; they are not biological observations, clinical recommendations or performance evidence. No Python programming is needed for the command-line routes. Run from the repository root after installing the package.

## 1 Calls and lineage labels

```bash
amr-clonalshare --config examples/workflows/calls.yaml --check-input
amr-clonalshare --config examples/workflows/calls.yaml --results-dir out/calls
```

Read `data/calls.csv` with `data/metadata.csv`. The declared meaning is R only under `clinical_sir`. Inspect input QC, the lineage share of the call with its interval for the represented lineages, its support, the bounds the call places on the share of the latent ordering and any method refusal. Replace the fictional source declaration with your actual source; add AST standard/version only when known.

## 2 Recorded MICs

```bash
amr-clonalshare --config examples/workflows/mic.yaml --check-input
amr-clonalshare --config examples/workflows/mic.yaml --results-dir out/mic
```

This route reads `data/mic.csv` without a calls table. The configuration records agent-specific wells and units. The run reports, per agent, the lineage share of the recorded MIC ordering with its interval for the represented lineages and its permutation p-value, and the bounds the readings place on the lineage share of the latent MIC ordering with the one-sided lower confidence limit of the lower bound they establish. The panel geometry and the share of end-well readings must accompany the result. MIC-only identifiers belong to the recorded frame. The share of the MIC ordering is a different quantity from the share of a call.

## 3 Compare two recorded collections

```bash
amr-clonalshare --config examples/workflows/contrasts.yaml --check-input
amr-clonalshare --config examples/workflows/contrasts.yaml --results-dir out/contrasts
```

Metadata declares `earlier` and `later` in that order; the difference is earlier minus later. Read composition and within-lineage terms with shared support and scope. With missing labels, a supported result describes the typed subset and cannot establish collection generalization. The fixture has complete labels; it does not validate performance under missingness.

## The same three routes in Python

```python
from amr_clonalshare import load_config, run

for recipe in ("calls", "mic", "contrasts"):
    config = load_config(f"examples/workflows/{recipe}.yaml")
    result = run(config, results_dir=f"out/{recipe}_api", seed=20260913)
```

## Preserve and interpret the bundle

Inspect the completion manifest, canonical JSON and method statuses before using the CSV or offline reports. Use new destination names for reruns. Explicit `--overwrite` replaces the earlier results files. Demonstration resampling budgets are intentionally small; they are unsuitable for precision claims. The validation campaign uses its separately recorded protocol and seeds. File-format examples do not replace provenance, sampling assessment or method-specific validation.
