# amr-clonalshare

Describe how recorded antimicrobial phenotypes vary across recorded lineages.
The package reads CSV tables (or `.xlsx` workbooks with the `excel` extra) of isolate identifiers, lineage labels and susceptibility calls or minimum inhibitory concentrations (MICs). It does not read genomes. Laboratory analysts can use the command line with a YAML configuration; bioinformaticians can use the same configuration through Python.

For every antimicrobial, amr-clonalshare 1.0.0 reports the lineage share of the calls: the out-of-sample skill of lineage means in predicting held-out isolates, with isolates assigned to folds within lineages and the cost of estimating each lineage mean removed lineage by lineage, set against the same score of permuted lineage labels. Its 95% interval is for the lineage share of the represented lineages: the lineages and their sizes are held fixed, the isolates of every lineage are drawn again from the smoothed empirical distribution of its readings, and every draw is studentized (a bootstrap-t interval); its ends are never closer to the estimate than the studentized permutation distribution allows, and its upper end also covers a share accounted for by one lineage whose readings sit close to the rest. Its p-value is a permutation test of the between-lineage sum of squares, with lineage labels exchanged within strata and, when a sampling unit such as a farm is declared, within units.

For MICs recorded as dilution intervals, it reports the lineage share of the MIC ordering: within each testing laboratory, crossed with any declared covariate, every reading is scored by its mid-distribution position under the nonparametric maximum-likelihood distribution of that stratum's readings, and the share of these scores accounted for by lineage is estimated, given an interval and tested as the share of a call is, with the scores recomputed in every bootstrap draw. It is the analogue, for the observed lineages, of the rank intraclass correlation (Tu, Li, Zeng and Shepherd, 2023) extended to censored readings; with one cut point it is the share of the call within strata. No distribution is assumed for the MICs.

The readings also bound the lineage share of the ordering of the latent MICs themselves, and a call bounds it as a reading with one cut point. The shares the readings allow form an interval: its lower end, the lower bound the readings establish whatever the order inside the wells, is a convex minimum over the base polytope of a submodular function and is computed with a certified error; its upper end is the largest share an order of the lineages attains, found exactly by dynamic programming for up to twenty lineages and, beyond, reported as a share the readings attain with a bound certified by branch and bound. Under the stated assumptions the bounds are a theorem, not a confidence statement. A one-sided 95% lower confidence limit of the lower bound, and so of the share, is estimated on random halves of the isolates, one half choosing the direction the other half tests. A call cut from a panel never establishes a higher lower bound than the panel.

The package also computes held-out e-values per antimicrobial with running evidence over intakes, a Kitagawa decomposition of a prevalence difference between two collections, exact recorded-frame bounds for missing binary outcomes, and a comparison of two lineage definitions on matched records with paired intervals for its terms; it repeats a run within the strata of a metadata column (`stratify_by`) and writes reports that show the panel per laboratory, the status of every result and a lineage-by-dilution heat map.

The claims that can fail on data are tested in a confirmatory campaign whose studies, designs, replicate counts, seeds and decision rules were committed before it runs (`benchmarks/CONFIRMATORY_PROTOCOL.md`): the coverage of the interval for the represented lineages, the level of the permutation test within strata and within sampling units, the bounds in every simulated dataset, the lower confidence limit and the paired intervals of the comparison of lineage definitions. `benchmarks/confirmatory_verdicts.py` applies the rules, and the verdicts of record are written to `benchmarks/results_confirmatory/`; a claim holds only for the designs whose rule it meets. In the campaign of record the interval met its coverage rule in all 224 designs and the paired intervals in all 25 terms, the bounds held in all 587,000 simulated datasets, the lower limit met its rule in all 130 designs, and the test met its level rule in 30 of 31 designs; the one exception, and every design, are listed there.

**Software version: 1.0.0.** Cite the version-specific Zenodo record; the [concept DOI](https://doi.org/10.5281/zenodo.22306353) groups all versions.

## Install and run

Install the software from PyPI:

```bash
python -m pip install amr-clonalshare
```

The executable recipes are part of the repository and of the source distribution; the wheel installs the software, not the example files. Their configurations name their data by paths relative to the recipe, so run them from a checkout:

```bash
git clone https://github.com/maciejkochanowski/amr-clonalshare
cd amr-clonalshare
amr-clonalshare check --config examples/workflows/calls.yaml
amr-clonalshare run --config examples/workflows/calls.yaml --results-dir out/calls
amr-clonalshare example ecoli --results-dir out/ecoli
```

The subcommands are `run`, `init`, `check`, `example`, `compare`, `gui`, `doctor`, `version` and `completion`; the same options work without a subcommand (`amr-clonalshare --config ...`), and `amr-clonalshare-compare` and `amr-clonalshare-gui` remain.

Installing from the checkout rather than from PyPI is `python -m pip install .`.

Python 3.11 or later and NumPy, pandas, SciPy, PyYAML and threadpoolctl are required. A container image with the pinned stack of `requirements-lock.txt` is published with every release as `ghcr.io/maciejkochanowski/amr-clonalshare:1.0.0` (see the [reproduction guide](docs/manual/08-container.md)).

Start with the [three executable recipes](examples/workflows/README.md): calls, MIC and two collections. Their small fictional data teach the file format; they are not biological evidence. The [input manual](docs/manual/02-input.md) explains how to substitute your own tables, and the [tutorial](docs/manual/11-tutorial.md) walks through the published *S. suis* analysis from the raw tables to the report.

For your own data the shortest path is one sheet carrying the isolate, its lineage and one column per antimicrobial, with calls or with recorded MICs. `amr-clonalshare init data/sheet.csv` reads it, prints a configuration with the mapping it guessed, the vocabulary it read from the values and the concentrations observed per agent, and stops; nothing it could not guess is invented.

## Without a command line

`amr-clonalshare-gui` opens the same run as a form in the browser: the tables are chosen or pasted from a spreadsheet, one sheet with the lineage beside one column per agent being enough, the columns are named from lists with the guesses of `--init` already made, the recorded concentrations of every agent are shown with a button that declares them as the tested range or a preset of an EU harmonised panel, and the run writes the same configuration, record and reports into a folder of its own, from which the command line repeats it. The results table of the run, in the form and on the console, ends with one conclusion per agent: whether the collection shows lineage structure in that measurement. The server answers this computer only. For a computer without Python, every release carries `AMR-ClonalShare-1.0.0-windows.zip`: unpack it anywhere and double-click `AMR-ClonalShare.bat`; nothing is installed. See the [manual page](docs/manual/13-local-form.md).

## Choose a question

| User question | Route | Interpretation |
|---|---|---|
| How strongly do the recorded labels describe the recorded calls? | Calls recipe | Lineage share of the call among the represented lineages, its interval, permutation test and support |
| What do recorded MIC intervals show at this label resolution? | MIC recipe | Lineage share of the MIC ordering, read within laboratories |
| What lower bound do the readings establish on the latent MIC ordering? | Calls or MIC recipe | Bounds on the lineage share of the latent ordering and a one-sided 95% lower confidence limit |
| How do two recorded collections differ? | Contrast recipe | Descriptive lineage-composition and within-lineage rate terms |
| Does another lineage definition read the calls differently? | `amr-clonalshare-compare` | Six analyses, five terms of their difference and paired intervals |

Clinical S/I/R, WT/NWT and binary data have separate `phenotype_kind` settings. New configurations should state the positive outcome, source and any applicable AST standard/version. With `phenotype_kind: clinical_sir` only R is positive unless the intermediate policy says otherwise. A configuration that declares no kind is read as `undeclared`: resistance words only, with I counted with R, and the run warns; declare the kind.

## Read the result bundle

A completed run saves `clonal_share_result.json`, a CSV summary with the p- and q-value of every row and a conclusion per agent, Markdown and offline HTML reports, input-QC records and a completion manifest; the console shows the results table, or the machine-readable summary with `--json`. JSON is the canonical analytical record; the other formats render its values. Inspect the manifest before treating a directory as a complete run. Existing output is protected; explicit `--overwrite` replaces the results files of an earlier run. Use a new output directory for a comparison run.

Read the phenotype definition, each method's retained collection, method status and reasons before comparing estimates. Small or incomplete data can support a description even when a particular interval is not computed. A missing method is not a zero estimate.

Missingness checks describe observed associations. A nonsignificant check or equal typing fractions cannot establish representativeness. When labels are missing, a supported decomposition describes the typed subset; collection generalization remains unsupported. Exact finite-collection bounds show what missing binary outcomes could change within the recorded frame; these are not confidence intervals and do not extend to unrecorded isolates.

## Python using the same configuration

```python
from amr_clonalshare import load_config, run

config = load_config("examples/workflows/calls.yaml")
result = run(config, results_dir="out/calls_api", seed=20260913)
```

See the [API guide](docs/api.md) for result reading and the [results manual](docs/manual/05-results.md) for the record schema (1.0).

## Evidence and limits

The empirical examples include a 481-isolate *Escherichia coli* MIC collection from diseased pigs, a 677-isolate *Streptococcus suis* collection and a 7,049-isolate poultry-meat *Salmonella* input cell; their provenance, exclusions and negative findings are recorded beside the data in `examples/`. The validation of record is the confirmatory campaign described above; [REPRODUCIBILITY.md](REPRODUCIBILITY.md) lists its stack, its commands and the receipts every run writes, and [benchmarks/README.md](benchmarks/README.md) the scripts behind it.

Lineage association does not establish transmission, a resistance mechanism, intervention benefit or clinical utility. Changing the lineage definition changes the question. The interval of a lineage share and the permutation test answer for the represented lineages, with their isolates as the random part; they do not describe lineages, farms, countries or years that were not observed.

## Documentation and licence

[User manual](docs/index.md) · [Methods](docs/methodology.md) · [API](docs/api.md) · [Changes](CHANGELOG.md) · [Code metadata](docs/CODE_METADATA.md)

Code is distributed under the MIT licence; see `LICENSE` and `Licence.txt`. The reports and the local form embed Latin Modern and TeX Gyre Heros, under the GUST Font License (`src/amr_clonalshare/fonts/`). Example-data provenance and licences are stated alongside each collection.

## Compare two lineage definitions on matched records

```bash
amr-clonalshare-compare --input examples/matched_lineages/input.csv \
  --id-column isolate_id --outcome positive --lineage-a lineage_a \
  --lineage-b lineage_b --output out/lineage_comparison
```

The six analyses separate record selection, the isolates each definition sets aside as singletons and relabelling on the isolates both definitions score; the two analyses on those common scorable isolates have the same identifiers and outcomes, and the five terms between the six analyses sum to the total difference. Every term has a paired 95% bootstrap-t interval (`paired_intervals`) for these records and labels: both are held fixed and the outcomes are drawn again. When an analysis is not estimable the comparison stays a diagnostic and no paired interval is computed. The example is fictional and teaches the format. `examples/ecoli_swine/input.csv` holds a real one: the ciprofloxacin call of 481 *E. coli* isolates with their sequence type and phylogroup (`--outcome ciprofloxacin_call --lineage-a sequence_type --lineage-b phylogroup`).
