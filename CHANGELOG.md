# Changelog

## [1.0.0] - 2026-09-16

- Report a calibrated 95% interval for the latent MIC variance ratio: maximum likelihood with Wald's exact generalised F pivot for exact readings, and a null-wise parametric-bootstrap likelihood-ratio inversion for dilution intervals. The pipeline computes it by default (`censored.calibrated_interval`, `calibration_n_boot`, `workers`) and the new `amr-clonalshare-mic` command runs it on one table.
- Keep the moment estimate and its approximate F interval as a labelled comparison; the F interval undercovered in simulation.
- Validate the new intervals in prespecified simulations (Gaussian, off-grid and stress designs); results and receipts are in `benchmarks/results_mic_release`.
- Accept maxima on the bounded nuisance box inside the bootstrap, which resolves full-range results at extreme candidate values; run bootstrap fits in parallel with unchanged results.
- Add the matched-record comparison of two lineage definitions (`amr-clonalshare-compare`).
- Repair tail probabilities, truncated moments, likelihood integration, bootstrap group identities, convergence checks and e-value arithmetic.
- Record a numerical failure of the MIC reading for the affected agent instead of ending the run, render reports when a sequential product of e-values is zero, and give result folders and files the permissions set by the umask.
- Take the cut points of the calibrated MIC interval from the recorded wells when a panel is configured, so that a tested but unobserved dilution is not dropped.
- Enter one categorical covariate of the MIC model as fixed effects (`dataset.mic_covariate_column`, `--covariate-column`): a laboratory or country shift is then removed from both variance components instead of being absorbed into whichever of them it is aligned with; the coefficients are re-estimated under every candidate ρ and in every simulated dataset.
- Read the dilution panel within each laboratory (`dataset.mic_panel_column`): the end wells are inferred per laboratory, the geometry is reported per laboratory, and the calibrated interval simulates each reading on its own panel. The S. suis example names its two testing laboratories; on the pooled lattice the lowest well of one had passed as an interior dilution of the other.
- Accept several covariates (`dataset.mic_covariate_columns`, repeated `--covariate-column`), each with its own reference level and additive offsets; validated in a simulated country-by-year design.
- Repeat the run within the levels of a metadata column (`dataset.stratify_by`) and add the strata to the record, the report and `strata_results.csv`.
- Add a screening budget for the calibrated interval (`censored.screening_n_boot`, `calibration_agents`) and record the budget each agent was given.
- Show in both reports the panel geometry per laboratory, the fixed-effect coefficients with their reference level, the status and reason of every calibrated interval, and a heat map of the readings per lineage and dilution interval.
- Publish the pinned container image to the GitHub Container Registry on every release.
- Remove what the reports never showed and the article does not describe: the lineage bootstrap, the conditional likelihood support and the end-well sensitivity arm of the moment estimate (`censored.n_boot`, `censored.sensitivity`), the effective number of carrying lineages and the lineage-mean prevalence (`surveillance.n_perm`), the bootstrap from the unrestricted fit that the null-wise test replaced, and the finite-reference research route of the population model with its registered surfaces. The moment estimate, its approximate F interval and the realised interval are unchanged.
- Refuse a covariate level that is censored on one side in every reading, whose effect the data cannot identify; treat the panel and covariate labels of a reading as part of it, so a blank label is refused at load and two rows that differ only in a label are a conflict; name every simulated dataset's covariate levels as the observed ones.
- Keep the fitted lineage probabilities of the e-values away from 0 and 1 in the grand-mean fallback, so that a first intake of identical outcomes cannot drive the running product to zero.
- Report the lineage structure behind each prediction score (lineages with repeated isolates, effective number of lineages, largest lineage share) in the record, the report and the matched-record comparison.
- Seed the Python entry points by default, accept `workers: 0` on the command line, honour `--no-check-files` in the run, accept `None` as a missing outcome in the collection bounds, and document the analysis budgets and the exit status 3 of the MIC command.
- Explain non-computed and full-range statuses with recorded reasons and next checks in both reports, without changing statistical methods or result status values.
- Separate clinical S/I/R, WT/NWT and binary phenotype declarations; retain explicit legacy replay and source/AST metadata.
- Collapse identical duplicates, reject conflicts by default and record explicit conflict exclusions.
- Track the union of raw call and MIC identifiers, method-specific retained cohorts, missing calls/labels and empty agents.
- Preserve supported typed-subset decomposition; missingness diagnostics do not certify collection representativeness.
- Add exact recorded-frame bounds for missing binary outcomes, configurable per-agent MIC wells and units, and four CLI/API recipes.
- Use schema 2.0 with historical schema 1.0 reading, canonical exports, protected output destinations, backups, completion manifests and optional progress.
- Reorganize the software documentation around user tasks.
