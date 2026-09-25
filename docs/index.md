# amr-clonalshare

Describe how recorded antimicrobial phenotypes vary across recorded lineages.
The package reads CSV tables (or `.xlsx` workbooks with the `excel` extra) of isolate identifiers, lineage labels and susceptibility calls or minimum inhibitory concentrations (MICs). It does not read genomes. Laboratory analysts can use the command line with a YAML configuration; bioinformaticians can use the same configuration through Python.

For MICs recorded as dilution intervals, amr-clonalshare 1.0.0 reports a calibrated 95% interval for the lineage share of MIC variation: the panel of each testing laboratory is read on its own wells, declared covariates (a laboratory, a country, a year) enter the model as fixed effects, and the interval inverts a likelihood-ratio test calibrated by a null-wise parametric bootstrap (Wald's generalised F pivot for exact readings). In simulations the interval met the prespecified coverage criterion in all 72 main Gaussian designs (single panel, no covariates) and in 45 of 51 additional designs (37 of 39 with heavy censoring, two-component residuals or laboratory and year entered as fixed effects); the simulation drivers and results are in `benchmarks/mic_inference/` and `benchmarks/results_mic_release/`. A moment estimate with an approximate F interval is reported beside it as a labelled comparison. The package also compares two lineage definitions on matched records, repeats a run within the strata of a metadata column (`stratify_by`), applies a screening budget to the calibrated interval, and writes reports that show the panel per laboratory, the fixed-effect coefficients, the status of every interval and a lineage-by-dilution heat map.

**Software version: 1.0.0.** Cite the version-specific Zenodo record; the [concept DOI](https://doi.org/10.5281/zenodo.22306353) groups all versions.

## Install and run

From PyPI or from this source directory:

```bash
python -m pip install amr-clonalshare   # or: python -m pip install .
amr-clonalshare --config examples/workflows/calls.yaml --check-input
amr-clonalshare --config examples/workflows/calls.yaml --results-dir out/calls
```

Python 3.11 or later and NumPy, pandas, SciPy and PyYAML are required.

Start with the [four executable recipes](manual/11-workflows.md): calls, MIC, two collections, and the optional population model. Their small fictional data teach the file format; they are not biological evidence. The [input manual](manual/02-input.md) explains how to substitute your own tables.

## Choose a question

| User question | Route | Interpretation |
|---|---|---|
| How strongly do the recorded labels describe the recorded calls? | Calls recipe | Observed-scale lineage-membership share, support and uncertainty |
| What do recorded MIC intervals show at this label resolution? | MIC recipe | Dilution-scale variance component with a calibrated 95% interval |
| How do two recorded collections differ? | Contrast recipe | Descriptive lineage-composition and within-lineage rate terms |
| What does a Gaussian lineage population model imply? | Population recipe | Separate latent-liability ICC and explicit model assumptions |

Clinical S/I/R, WT/NWT and binary data have separate `phenotype_kind` settings. New configurations should state the positive outcome, source and any applicable AST standard/version. With `phenotype_kind: clinical_sir` only R is positive unless the intermediate policy says otherwise. A configuration that declares no kind is read as `undeclared`: resistance words only, with I counted with R, and the run warns; declare the kind.

## Read the result bundle

A completed run saves `clonal_share_result.json`, a CSV summary, Markdown and offline HTML reports, input-QC records and a completion manifest. JSON is the canonical analytical record; the other formats render its values. Inspect the manifest before treating a directory as a complete run. Existing output is protected; explicit `--overwrite` replaces the results files of an earlier run. Use a new output directory for a comparison run.

Read the phenotype definition, each method's retained cohort, method status and reasons before comparing estimates. Small or incomplete data can support a description even when a particular interval is withheld. A missing method is not a zero estimate.

Missingness checks describe observed associations. A nonsignificant check or equal typing fractions cannot establish representativeness. When labels are missing, a supported decomposition describes the typed subset; collection generalization remains unsupported. Exact finite-collection bounds show what missing binary outcomes could change within the recorded frame; these are not confidence intervals and do not extend to unrecorded isolates.

## Python using the same configuration

```python
from amr_clonalshare import load_config, run

config = load_config("examples/workflows/calls.yaml")
result = run(config, results_dir="out/calls_api", seed=20260913)
```

See the [API guide](api.md) for result reading and the [results manual](manual/05-results.md) for the record schema (1.0).

## Evidence and limits

The empirical examples include a 677-isolate *Streptococcus suis* collection and a 7,049-isolate poultry-meat *Salmonella* input cell; their provenance, exclusions and negative findings are recorded beside the data in `examples/`. The MIC interval calibration and the *S. suis* and *Salmonella* reanalyses are in `benchmarks/results_mic_release`. Every result under `benchmarks/results_*` and every validation file the package reads was produced by one validation campaign, run from one commit on one pinned software stack; the repository-root `REPRODUCIBILITY.md` lists its commands and receipts.

Lineage association does not establish transmission, a resistance mechanism, intervention benefit or clinical utility. Changing the lineage definition changes the question. The optional population model was calibrated and validated in a prespecified simulation study (`benchmarks/population_model/PROTOCOL.md`): its fixed-cut-off interval covered the generating value in 95.9% to 100% of datasets across 144 Gaussian designs, and far less often when lineage effects followed a two-point law or lineage sizes grew with the effects. Population-model intervals depend on their stated assumptions, and passing a diagnostic does not certify those assumptions. The optional general Gaussian-probit route is selected explicitly with `interval_method: general`; it does not change the default. Read its [validation and resource guidance](GENERAL_POPULATION_ICC.md).

## Documentation and licence

[User manual](index.md) · [Methods](methodology.md) · [API](api.md) · [Code metadata](CODE_METADATA.md)

Code is distributed under the MIT licence; see `LICENSE` and `Licence.txt`. Example-data provenance and licences are stated alongside each collection.
