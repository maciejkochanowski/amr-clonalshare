# 4 Run an analysis

```bash
amr-clonalshare --config examples/workflows/calls.yaml --results-dir out/calls
```

Choose a separate output directory for each phenotype definition, lineage resolution or collection. A normal run stages its files and writes a completion manifest. Treat a directory without a successful completion manifest as incomplete. An existing destination is protected. When replacement is intentional, `--overwrite` keeps a backup of the preceding directory.

```bash
amr-clonalshare --config examples/workflows/calls.yaml --results-dir out/calls --overwrite
```

The switches are `--seed` (default 42), `--threads` (see the installation page), `--check-input` (run the input check only), `--quiet` (no progress lines and no summary on standard output; the files are written as usual), `--overwrite` and `--no-check-files` (do not verify that the configured input files exist before loading, for a configuration that is only being validated). Exit status zero includes completed analyses that withhold unsupported methods; inspect result status. Configuration/input errors return status two.

The Python entry point accepts the same loaded configuration and supports `results_dir`, `seed`, `overwrite=False` and an optional `progress` callback. Preserve the configuration, software version, seed, input files and manifest with the output. Historical evidence has its own recorded seeds and numerical budgets; the small recipes use deliberately modest demonstration budgets.

## Analysis budgets

Every stochastic step has a budget in the configuration. The defaults suit a routine run; the values used for the reported reanalyses are given beside them.

| Section | Key (default) | Meaning |
|---|---|---|
| `attribution`, `evidence`, `surveillance`, `censored` | `enabled` (true) | Switch the whole family off; a family that is off writes no rows for its analyses |
| `attribution` | `folds` (5), `repeats` (20) | Cross-validation folds and repeated fold draws for the prediction score |
| `attribution` | `n_perm` (200), `n_boot` (400) | Label permutations for the null correction and the p-value; lineage-bootstrap draws for the limits |
| `evidence` | `folds` (5), `repeats` (20), `alpha` (0.05) | Split-likelihood-ratio e-values and the e-BH level |
| `surveillance` | `n_boot` (2000), `q_fdr` (0.05), `min_shared_support` (0.8), `label_alpha` (0.05) | The decomposition of a prevalence difference between two collections: its bootstrap budget, the false-discovery level and the label-availability tests |
| `censored` | `end_wells_censored` (true) | Whether a reading on an end well of the panel enters as a censored interval (the coarsened-at-random reading) or as an exact value |
| `censored` | `calibrated_interval` (true), `calibration_n_boot` (199), `workers` (0) | The calibrated MIC interval, simulated datasets per tested value of rho, and processes (0: every CPU the process may use) |
| `censored` | `screening_n_boot` (unset), `calibration_agents` ([]) | A smaller simulation budget for every agent except those listed, which keep `calibration_n_boot` |

`calibration_n_boot: 199` resolves the p-value in steps of 0.005 and was the budget of the simulation study; it is the setting for routine reporting. For an interval that will be quoted, `999` narrows the Monte Carlo error of the limits at about five times the run time; the limits themselves are reported to a tolerance of 0.002 in rho and inherit the Monte Carlo error of the test at each tested value. Results do not depend on `workers`.

`screening_n_boot` runs a screening pass: every agent is calibrated with that smaller budget (at least 19 draws) and only the agents named in `calibration_agents` receive the full `calibration_n_boot`. The record stores the budget each agent was given (`n_boot` in its calibrated block) and the report lists it. Use it on a panel of many agents to find the ones worth the full budget, then rerun those with `calibration_n_boot` alone; a screening interval is not the interval to quote.

## What the lineage label is

The lineage column is an operational grouping of isolates: a sequence type, a core-genome or SNP cluster, a serovar or any other label the user supplies. Labels of different schemes have different depths, and the package does not assume a common phylogenetic depth or a genetic mechanism behind a label. The prediction score describes how well membership of the observed lineages predicts the phenotype of further isolates from those lineages; it does not describe prediction for lineages, farms, countries or years that were not observed, and any host, farm, laboratory or period effect that follows lineage is part of it. Read the lineage-structure table of the report (lineages with repeated isolates, effective number of lineages, support) beside every score.

The population module remains disabled unless requested. `population_probit: {enabled: true}` keeps `fixed_v3`; explicitly specify `interval_method: general` to request the general Gaussian-probit computation. That computation can be expensive and can return incomplete status. See the general-method manual before a large run.
