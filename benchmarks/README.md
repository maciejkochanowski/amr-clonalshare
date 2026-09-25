# Benchmarks

The scripts that produce the calibration and application results, and the
results folders they write. Every results folder was written by one run of the
validation campaign (`campaign/CAMPAIGN.md`) and carries a receipt with the
sha256 of every file in it.

| script | what it measures | writes |
|---|---|---|
| `attribution_calibration.py` | the clonal share against a known truth; the support gate | `results_attribution_calibration/` |
| `estimator_benchmark.py`, `estimator_grid_summary.py` | the six estimators on one grid, three truths a cell; the validation grid the package reads | `results_estimator_grid/`, `src/amr_clonalshare/validation_grid.json` |
| `realised_calibration.py` | the exact interval of the realised share and its kurtosis gate | `results_realised_calibration/` |
| `censored_gates.py` | the censored-reading gates: share of lineages beyond the panel, prevalence window of a single cut | `results_censored_gates/` |
| `decomposition_calibration.py` | the Kitagawa split against a known truth; the shared-support gate | `results_decomposition_calibration/` |
| `decomposition_vs_regression.py` | the split against a lineage-adjusted regression | `results_decomposition_vs_regression/` |
| `repeated_looks.py` | false-discovery control when a panel is re-read every year | `results_repeated_looks/` |
| `null_uniformity.py` | the permutation p-value under an exact null | `results_null_uniformity/` |
| `ssuis_mechanism.py`, `ssuis_resolution.py` | the *S. suis* reading against its resistance determinants, and at three lineage definitions | `results_ssuis_mechanism/`, `results_ssuis_resolution/` |
| `mic_inference/` (see `PROTOCOL.md`), `empirical/` | the calibration of the MIC variance-ratio interval and the *S. suis* and *Salmonella* reanalyses | `results_mic_release/` |
| `population_model/` (see `PROTOCOL.md`) | calibration and validation of the optional population liability ICC | `results_population_model/`, `src/amr_clonalshare/population_probit_validation.json`, `src/amr_clonalshare/general_probit_protocol.json` |
| `seed_stability.py` | the spread of each agent's share over master seeds | `results_seed_stability/` |
| `profile_run.py` | time and memory of one analysis, by phase | `results_profile/` |
| `fetch_pathogen_detection.sh`, `pathogen_detection_releases.tsv`, `vet_source_taxonomy.py`, `ast_calls.py` | retrieval of the public source at pinned release accessions and the rules that build the *Salmonella* example from it | `raw/` (not shipped), `examples/salmonella_poultry/data/` |
