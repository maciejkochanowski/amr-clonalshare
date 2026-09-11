# Benchmarks

The scripts that produced the calibration and application artefacts of the
article, and the results folders a document reads. Every script names the
file it writes; the article's supplement names, by path and digest, the file
each of its numbers was read from.

| script | what it measures | writes |
|---|---|---|
| `attribution_calibration.py` | the clonal share against a known truth; the support gate | `results_attribution_calibration/` |
| `estimator_benchmark.py` | the six estimators on one grid, three truths a cell | the evidence deposit |
| `realised_calibration.py` | the exact interval of the realised share and its kurtosis gate | the evidence deposit |
| `decomposition_calibration.py` | the Kitagawa split against a known truth; the shared-support gate | `results_decomposition_calibration_2026-08-31/` |
| `decomposition_vs_regression.py` | the split against a lineage-adjusted regression | `results_decomposition_vs_regression_2026-08-31/` |
| `censored_calibration.py`, `censored_grid.py` | the dilution reading against a known truth, one design and a grid | the evidence deposit; `results_censored/` |
| `censored_real_cohort.py` | the dilution reading on the shipped *S. suis* panel | `results_censored/` |
| `repeated_looks.py` | false-discovery control when a panel is re-read every year | the evidence deposit |
| `null_uniformity.py` | the permutation p-value under an exact null | `results_null_uniformity_2026-09-04/` |
| `ssuis_mechanism.py`, `ssuis_resolution.py` | the *S. suis* reading against its determinants, and at three lineage definitions | `results_ssuis_mechanism/`, `results_ssuis_resolution/` |
| `comparator_arms.py`, `comparator_report.py` | the population-structure comparators and the restricted-likelihood anchor | the evidence deposit |
| `leverage_check.py` | the leverage of each cut-off on the *S. suis* reading | the evidence deposit |
| `fetch_pathogen_detection.sh`, `pathogen_detection_releases.tsv` | retrieval of the public source at pinned release accessions | `raw/` (not shipped) |
| `atlas_cross_species.py`, `agent_screen.py` | the estimator on fifty organisms of the public source | the evidence deposit |
| `vet_cohorts.py`, `vet_source_taxonomy.py`, `vet_atlas.py`, `resolution_atlas.py`, `verify_vet_claims.py`, `period_split.py` | the veterinary cut of the source, at two lineage definitions, and the poultry cell of the article | the evidence deposit |
