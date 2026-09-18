# Benchmarks

The scripts that produce the calibration and application results, and the
results folders they write. Every script names the file it writes, and every
results folder carries the digests of its files.

| script | what it measures | writes |
|---|---|---|
| `attribution_calibration.py` | the clonal share against a known truth; the support gate | `results_attribution_calibration/` |
| `estimator_benchmark.py` | the six estimators on one grid, three truths a cell | the directory given on the command line (not part of the reported results) |
| `realised_calibration.py` | the exact interval of the realised share and its kurtosis gate | the directory given on the command line (not part of the reported results) |
| `decomposition_calibration.py` | the Kitagawa split against a known truth; the shared-support gate | `results_decomposition_calibration_2026-08-31/` |
| `decomposition_vs_regression.py` | the split against a lineage-adjusted regression | `results_decomposition_vs_regression_2026-08-31/` |
| `repeated_looks.py` | false-discovery control when a panel is re-read every year | the directory given on the command line (not part of the reported results) |
| `null_uniformity.py` | the permutation p-value under an exact null | `results_null_uniformity_2026-09-04/` |
| `ssuis_mechanism.py`, `ssuis_resolution.py` | the *S. suis* reading against its determinants, and at three lineage definitions | `results_ssuis_mechanism/`, `results_ssuis_resolution/` |
| `comparator_arms.py`, `comparator_report.py` | the population-structure comparators and the restricted-likelihood anchor | the directory given on the command line (not part of the reported results) |
| `leverage_check.py` | the leverage of each cut-off on the *S. suis* reading | the directory given on the command line (not part of the reported results) |
| `mic_inference/` (see `PROTOCOL.md`), `empirical/` | the calibration of the MIC variance-ratio interval and the *S. suis* and *Salmonella* reanalyses | `results_mic_release/` |
| `fetch_pathogen_detection.sh`, `pathogen_detection_releases.tsv` | retrieval of the public source at pinned release accessions | `raw/` (not shipped) |
| `atlas_cross_species.py`, `agent_screen.py` | the estimator on fifty organisms of the public source | the directory given on the command line (not part of the reported results) |
| `vet_cohorts.py`, `vet_source_taxonomy.py`, `vet_atlas.py`, `resolution_atlas.py`, `verify_vet_claims.py`, `period_split.py` | the veterinary cut of the source, at two lineage definitions, and the poultry cell | the directory given on the command line (not part of the reported results) |
