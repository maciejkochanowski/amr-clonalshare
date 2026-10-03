# 4 Run an analysis

```bash
amr-clonalshare run --config examples/workflows/calls.yaml --results-dir out/calls
```

The command has subcommands: `run` (an analysis), `init` (a draft configuration from the tables), `check` (the input check only), `example` (a shipped collection with the seed of its published record), `compare` (two lineage definitions, the same as `amr-clonalshare-compare`), `gui` (the local form, the same as `amr-clonalshare-gui`), `doctor` (the versions, threads and examples this computer offers), `version` and `completion` (a completion script for bash or PowerShell, to source or dot-source). The options below also work without a subcommand, as in `amr-clonalshare --config ...`, and the old forms stay valid.

```bash
amr-clonalshare example ecoli --results-dir out/ecoli
amr-clonalshare example            # lists the shipped collections
amr-clonalshare doctor
amr-clonalshare completion bash > ~/.amr-clonalshare-completion && source ~/.amr-clonalshare-completion
```

Choose a separate output directory for each phenotype definition, lineage resolution or collection. A normal run stages its files and writes a completion manifest. Treat a directory without a successful completion manifest as incomplete. An existing destination is protected. When replacement is intentional, `--overwrite` replaces the results files once the new run has completed and keeps any other file in the directory; it refuses a directory that holds subdirectories or the working directory, since that is a project folder rather than a results folder.

```bash
amr-clonalshare --config examples/workflows/calls.yaml --results-dir out/calls --overwrite
```

When the run ends, the console shows the results table the local form shows, one line per agent with the share, its interval, the p- and q-values, the bounds and the conclusion, followed by the conclusion for every agent in full and the paths of the report, the results table and the record; `--json` prints the machine-readable summary of the run instead; without `--results-dir` nothing is written and the console says so. The switches are `--seed` (default 42), `--threads` (see the installation page), `--check-input` (run the input check only), `--quiet` (no progress lines and nothing on standard output except the input check of `--check-input`; the files are written as usual), `--overwrite`, `--no-check-files` (do not verify that the configured input files exist before loading, for a configuration that is only being validated), `--permutations N`, `--bootstrap N` (also spelled `--bootstraps`) and `--alpha A` (the budgets and the level of the configuration replaced from the command line, checked as the file's values are and written into the record as run) and `--dry-run` (print the configuration as the run would read it, every default filled in, and stop). The progress lines name every agent as it is scored and, from the second agent on, the time taken so far and about how long the rest will take; the same lines, with the warnings of the run, are kept in `run.log` in the results folder, and a run that fails leaves `run.log` with its error there as the only file, so that one file says what happened; a failed run over a finished one with `--overwrite` leaves the finished run as it was and writes `run_failed.log` beside it. Exit status zero includes completed analyses that withhold unsupported methods; inspect result status. Configuration, input and run errors return status two.

The Python entry point accepts the same loaded configuration and supports `results_dir`, `seed`, `overwrite=False` and an optional `progress` callback. Preserve the configuration, software version, seed, input files and manifest with the output. The validation campaign has its own recorded seeds and budgets.

## Analysis budgets

The resampling budgets are set in the configuration. The defaults suit a routine run; the small recipes use deliberately modest budgets that suit a demonstration only.

| Section | Key (default) | Meaning |
|---|---|---|
| `attribution`, `evidence`, `surveillance`, `censored` | `enabled` (true) | Switch the whole family off; a family that is off writes no rows for its analyses |
| `attribution` | `folds` (5), `repeats` (20) | Cross-validation folds and repeated fold assignments for the lineage share, of the calls and of the MIC ordering |
| `attribution` | `n_perm` (999), `n_boot` (999) | Label permutations for the control, the p-value and the variance of the share without lineage differences; draws of the smoothed bootstrap for the interval for the represented lineages (0: no interval) |
| `evidence` | `folds` (5), `repeats` (20), `alpha` (0.05) | Split-likelihood-ratio e-values; `alpha` is the level of e-BH, of the rule for a programme of intakes and of the Benjamini-Yekutieli selection of agents by their permutation p-values |
| `surveillance` | `n_boot` (2000), `q_fdr` (0.05), `min_shared_support` (0.9), `label_alpha` (0.05) | The decomposition of a prevalence difference between two collections: its bootstrap budget, the false-discovery level, the shared-support gate and the level of the label-availability tests |
| `dataset` | `phenotype_agent_columns` (unset) | Names the antimicrobial columns of a wide call table; the table is melted to the long shape before it is read, and the setting replaces `phenotype_antibiotic_column` and `phenotype_call_column` |
| `censored` | `end_wells_censored` (true) | Whether a reading on an end well of the panel is censored (the coarsened-at-random reading) or read as an interior well is, one doubling wide |

The work of the lineage share is linear in its budgets: per antimicrobial, `repeats + n_perm` cross-validated scores and `n_boot + n_perm` averages of the share over the fold assignments, whatever the number of folds (`tests/test_work_budget.py`). The smallest p-value `n_perm` permutations can return is 1/(n_perm + 1); with a wide panel, raise `n_perm` until that floor falls below the first threshold of the Benjamini-Yekutieli step-up, which the default clears for up to 15 agents at 0.05. The bounds on the share of the latent ordering have no budget in the configuration: their lower confidence limit always averages 25 random splits, and their upper end is exact for up to 20 lineages and at least certified beyond. The comparison of two lineage definitions takes its budgets on its own command line (`--folds`, `--repeats`, `--permutations` and `--bootstraps`, with the defaults of the run: 5, 20, 999 and 999).

## What the lineage label is

The lineage column is an operational grouping of isolates: a sequence type, a core-genome or SNP cluster, a serovar or any other label the user supplies. Labels of different schemes have different depths, and the package does not assume a common phylogenetic depth or a genetic mechanism behind a label. The prediction score describes how well membership of the observed lineages predicts the phenotype of further isolates from those lineages; it does not describe prediction for lineages, farms, countries or years that were not observed, and any host, farm, laboratory or period effect that follows lineage is part of it. Read the lineage-structure table of the report (lineages with repeated isolates, effective number of lineages, support) beside every score.
