# Methods and interpretation

amr-clonalshare reports several quantities with different targets. A comparison must preserve phenotype meaning, lineage resolution, retained isolates, units and the target of uncertainty. A confidence interval is not a sampling-representativeness test.

## Collection attribution and conditional components

`attribution.clonal_share` evaluates lineage means on held-out isolates and corrects the score against permuted labels. Its primary observed-scale lineage-membership share differs from the design-corrected realised component ratio. The lineage bootstrap and the model-based population interval have different targets. A score with support below 0.80 is withheld from the report and kept in the result record; few lineages and dominant groups are reported; passing a support diagnostic does not guarantee coverage.

`realised.realised_share` supplies a conditional noncentral-F interval under Gaussian within-lineage errors. The kurtosis criterion may withhold the interval, but it does not establish Gaussianity: in the validation families of `benchmarks/realised_calibration.py`, residuals from a t distribution on four degrees of freedom lowered its coverage to 0.87. The prevalence-based latent restatement of the share is exploratory; do not substitute it for the directly fitted population model.

## MIC measurements

The primary MIC result is a maximum-likelihood estimate of the latent variance ratio with a calibrated 95% interval (`calibrated` in the record, `censored_mic_calibrated` in `results.csv`). When every reading is exact, the interval inverts Wald's generalised F pivot, which is exact for unequal lineage sizes. For dilution intervals, each candidate ratio is tested with a likelihood-ratio statistic whose null distribution is simulated from the fit under that candidate (199 draws, same lineage sizes and panel); the interval is the numerical hull of kept values. In simulations both routes met the prespecified coverage criterion in all 72 main Gaussian designs (single panel, no covariates), and in 15 of 16 other Gaussian designs; contaminated or heavy-tailed residuals lowered coverage. Results and receipts are in `benchmarks/results_mic_release`.

An empirical-Bayes moment iteration with approximate F limits is reported in the `share` block as a labelled comparison. It is not exact EM, ML or REML and undercovered in simulation. Profile fields are conditional likelihood-support diagnostics, not confidence intervals.

The censored likelihood reads recorded dilution intervals on the log2 scale. Explicit operators and configured agent-specific wells define measurement geometry; units are retained without automatic conversion. End-well sensitivity assesses one measurement assumption, not every source of model error. The dilution-scale variance component differs from binary-call attribution. The MIC cohort includes eligible MIC-only isolates.

## Two collections and incomplete observation

The Kitagawa identity decomposes the prevalence of the first level of `contrast_levels` minus that of the second into composition and within-lineage terms. Shared support governs the within-lineage comparison. A supported calculation may describe the typed subset despite missing labels; collection generalization remains unsupported whenever labels are missing. Equal typing fractions and nonsignificant association checks do not establish representativeness. Benjamini–Yekutieli adjustment applies within each reported component family and does not repair selection.

Missing binary outcomes produce exact recorded-frame prevalence bounds from the extreme all-negative and all-positive completions. These are not confidence intervals or imputed estimates. The reported prevalence point concerns observed outcomes; the bounds concern all isolates in the stated recorded frame. Neither establishes a parameter outside that frame.

## Population and sequential methods

The optional grouped Gaussian-probit analysis targets a latent-liability population ICC under independent Gaussian group effects, conditional binomial sampling, noninformative sizes and ignorable selection. It is disabled by default; enabling it retains `fixed_cutoff`, and `general` must be selected explicitly. The [general-method documentation](GENERAL_POPULATION_ICC.md) describes geometry branches, approximate fitted-nuisance inference, resources and its finite-panel validation.

Fixed-look evidence and sequential evidence are distinct. The sequential construction uses past-only prediction and requires the specified conditional null model and actual batch order. It is not an unrestricted permission to reuse fixed-look tests, union rejection sets or search grouping definitions. E-values measure evidence, not the magnitude of a variance component.
