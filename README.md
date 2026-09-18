# amr-clonalshare

## Version 1.0.0

This release reports a calibrated 95% interval for the lineage share of MIC variation from dilution intervals: the panel of each testing laboratory is read on its own wells, declared covariates (a laboratory, a country, a year) enter the model as fixed effects, and the interval inverts a likelihood-ratio test calibrated by a null-wise parametric bootstrap (Wald's generalised F pivot for exact readings). In prespecified simulations the interval kept nominal coverage in all 72 main Gaussian designs and in 37 of 39 designs with heavy censoring, two-component residuals and laboratory or year effects; the simulation drivers and results are in `benchmarks/mic_inference/` and `benchmarks/results_mic_release/`. The release also adds a matched-record comparison of two lineage definitions, a run repeated within the strata of a metadata column (`stratify_by`), a screening budget for the calibrated interval, and reports that show the panel per laboratory, the fixed-effect coefficients, the status of every interval and a lineage-by-dilution heat map. The earlier moment estimate and its approximate F interval stay in the output as a labelled comparison.

Describe how recorded antimicrobial phenotypes vary across recorded lineages.
The package reads CSV tables of isolate identifiers, lineage labels and susceptibility calls or minimum inhibitory concentrations (MICs). It does not read genomes. Laboratory analysts can use the command line with a YAML configuration; bioinformaticians can use the same configuration through Python.

**Software version: 1.0.0.** Cite the version-specific Zenodo record; the [concept DOI](https://doi.org/10.5281/zenodo.22306353) groups all versions.

## Install and run

From PyPI or from this source directory:

```bash
python -m pip install amr-clonalshare   # or: python -m pip install .
amr-clonalshare --config examples/workflows/calls.yaml --check-input
amr-clonalshare --config examples/workflows/calls.yaml --results-dir out/calls
```

Python 3.11 or later and NumPy, pandas, SciPy and PyYAML are required. A container image with the pinned stack of the shipped records is published with every release as `ghcr.io/maciejkochanowski/amr-clonalshare:1.0.0` (see the [reproduction guide](docs/manual/08-container.md)).

Start with the [four executable recipes](examples/workflows/README.md): calls, MIC, two collections, and the optional population model. Their small fictional data teach the file format; they are not biological evidence. The [input manual](docs/manual/02-input.md) explains how to substitute your own tables.

## Choose a question

| User question | Route | Interpretation |
|---|---|---|
| How strongly do the recorded labels describe the recorded calls? | Calls recipe | Observed-scale lineage-membership share, support and uncertainty |
| What do recorded MIC intervals show at this label resolution? | MIC recipe | Dilution-scale variance component with a calibrated 95% interval |
| How do two recorded collections differ? | Contrast recipe | Descriptive lineage-composition and within-lineage rate terms |
| What does a Gaussian lineage population model imply? | Population recipe | Separate latent-liability ICC and explicit model assumptions |

Clinical S/I/R, WT/NWT and binary data have separate `phenotype_kind` settings. New configurations should state the positive outcome, source and any applicable AST standard/version. Clinical S/I/R defaults to R only; I is not silently treated as resistant. Legacy configurations remain readable for explicit historical replay.

## Read the result bundle

A completed run saves `clonal_share_result.json`, a CSV summary, Markdown and offline HTML reports, input-QC records and a completion manifest. JSON is the canonical analytical record; the other formats render its values. Inspect the manifest before treating a directory as a complete run. Existing output is protected; explicit `--overwrite` preserves a backup. Use a new output directory for a comparison run.

Read the phenotype definition, each method's retained cohort, method status and reasons before comparing estimates. Small or incomplete data can support a description even when a particular interval is withheld. A missing method is not a zero estimate.

Missingness checks describe observed associations. A nonsignificant check or equal typing fractions cannot establish representativeness. When labels are missing, a supported decomposition describes the typed subset; collection generalization remains unsupported. Exact finite-collection bounds show what missing binary outcomes could change within the recorded frame; these are not confidence intervals and do not extend to unrecorded isolates.

## Python using the same configuration

```python
from amr_clonalshare import load_config, run

config = load_config("examples/workflows/calls.yaml")
result = run(config, results_dir="out/calls_api", seed=20260913)
```

See the [API guide](docs/api.md) for result reading and the [results manual](docs/manual/05-results.md) for schema 2.0. The reader also accepts historical schema 1.0 records without manufacturing current fields.

## Evidence and limits

The empirical examples include a 677-isolate *Streptococcus suis* collection and a 7,049-isolate poultry-meat *Salmonella* input cell; their provenance, exclusions and negative findings are recorded beside the data in `examples/`. The MIC interval calibration and the *S. suis* and *Salmonella* reanalyses of this release are in `benchmarks/results_mic_release`; the other benchmark drivers in `benchmarks/` are not part of the reported results and write their output where they are pointed. [REPRODUCIBILITY.md](REPRODUCIBILITY.md) explains the distinction.

Lineage association does not establish transmission, a resistance mechanism, intervention benefit or clinical utility. Changing the lineage definition changes the question. Population-model intervals depend on their stated assumptions, and passing a diagnostic does not certify those assumptions. The optional general Gaussian-probit route is selected explicitly with `interval_method: general`; it does not change the default. Read its [validation and resource guidance](docs/GENERAL_POPULATION_ICC.md).

## Documentation and licence

[User manual](docs/index.md) · [Methods](docs/methodology.md) · [API](docs/api.md) · [Changes](CHANGELOG.md) · [Code metadata](docs/CODE_METADATA.md)

Code is distributed under the MIT licence; see `LICENSE` and `Licence.txt`. Example-data provenance and licences are stated alongside each collection.

## Compare two lineage definitions on matched records

```bash
amr-clonalshare-compare --input examples/matched_lineages/input.csv \
  --id-column isolate_id --outcome positive --lineage-a lineage_a \
  --lineage-b lineage_b --output out/lineage_comparison
```

The four arms distinguish record restriction from label changes. The two common-frame arms have the same IDs and outcomes; unsupported arms remain diagnostics. Arm-wise bootstrap intervals are not intervals for paired differences. The example is fictional and teaches the format.

## Interval for one MIC table

```bash
amr-clonalshare-mic readings.csv --method bootstrap --seed 1 \
  --panel-edges "[-3,-2,-1,0,1,2,3]" --workers 8 --output rho.json
```

The input has columns `lo`, `hi` (log2 concentration bounds; equal for an exact reading) and `lineage`. When the records were read on more than one panel, name the column that says which with `--panel-column` and pass `--panel-edges` as a JSON object, one array of cut points per panel name; `--covariate-column` names a column whose levels (a laboratory, a country) enter the model as fixed effects. Use `--method exact` when every reading is exact. The main pipeline computes the same interval for every agent unless `censored.calibrated_interval: false`; `censored.workers` sets the number of processes and does not change the result. One interval for about 700 isolates takes several minutes on 16 cores.
