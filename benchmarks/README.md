# Benchmarks

The scripts that test the package on simulated and recorded data, and the
results folders they write. The validation of record is the confirmatory
campaign of `CONFIRMATORY_PROTOCOL.md`: the protocol was committed before any
run and fixes every study, design, replicate count, root seed and decision
rule; `confirmatory_verdicts.py` applies the rules to the outputs of the runs
and names every design that fails one. Its verdicts, with the receipts of the
steps that wrote them, are written to `results_confirmatory/`. A claim holds
only for the designs whose rule it meets. The repository-root
`REPRODUCIBILITY.md` gives the stack and the commands.

| script | what it measures | writes |
|---|---|---|
| `CONFIRMATORY_PROTOCOL.md`, `confirmatory_verdicts.py` | the protocol of the confirmatory campaign and the rules that decide it: coverage of the interval for the represented lineages and of the paired intervals of the comparison of lineage definitions, level of the permutation test within strata and within sampling units, the bounds on the share of the latent ordering in every dataset, and the lower confidence limit | `results_confirmatory/` |
| `estimator_benchmark.py` | Study 1: the lineage share of calls and traits on a grid of 175 cells, against the in-sample R², the ANOVA variance component and REML | per-cell files and, with `--aggregate`, `estimator_benchmark.json` in the directory given with `--out` |
| `order_calibration.py` | Study 2: coverage and level of the lineage share of the MIC ordering, and the bounds and lower limit against the latent share, on 23 designs | the file given with `--output`, one per chunk of replicates |
| `conditional/` (`designs.py`, `run.py`) | Study 3: designs with the lineages, their sizes and their distributions fixed, so every target is exact: the interval for the represented lineages, the bounds and the lower limit, the test within sampling units, and the paired intervals of the comparison | the file given with `--output`, one per chunk of replicates |
| `competitors/` (`compare_order.py`, `compare_order.R`, `install_r_bench.R`) | Study 4: seven published procedures for a lineage effect on MIC readings, run in R on the datasets of the test of the ordering share | one file per chunk, and a summary |
| `campaign/` | the queues of the confirmatory campaign (`make_confirmatory_queues.py`), the Slurm scripts that run them (`launch_confirmatory.sh`, `queue.sbatch`, `run_queue.py`, `run_summaries.sh`, `env.sh`), the receipt written beside every summary step (`run_logged.py`), and the rerun of the empirical, auxiliary, check and mutation queues against a new commit (`launch_reissue.sh`) with the script that places its outputs (`collect_reissue.py`) | the campaign directory |
| `decomposition_calibration.py`, `decomposition_gate.py` | the Kitagawa split against a known truth, over collection sizes, lineage counts, turnover and prevalence; the shared-support gate read off that grid | `results_decomposition_calibration/` (`decomposition_gate.json`) |
| `decomposition_vs_regression.py` | the split against a lineage-adjusted regression | `results_decomposition_vs_regression/` |
| `repeated_looks.py` | false-discovery control when a panel is re-read every year | `results_repeated_looks/` |
| `null_uniformity.py` | the permutation p-value under an exact null, over its whole distribution | `results_null_uniformity/` |
| `ssuis_mechanism.py`, `ssuis_resolution.py` | the *S. suis* reading against its resistance determinants, and at three lineage definitions | the directory given with `--out` |
| `empirical/` (`ssuis_analysis.py`, `ssuis_panel.py`, `ssuis_sensitivity.py`, `salmonella_analysis.py`) | the *S. suis* and *Salmonella* analyses beyond the configured runs: comparisons of lineage definitions, the ordering share per agent, the share as the panel is coarsened to fewer wells, a decomposition, and the sensitivity of the MIC and call results to the choices the analysis makes; `campaign/empirical_queue.txt` runs them with the configured analyses | the directory given, `results_empirical/` by default |
| `seed_stability.py` | the spread of each agent's share over master seeds | the file given |
| `profile_run.py` | time and memory of one analysis, by phase | the file given |
| `fetch_pathogen_detection.sh`, `pathogen_detection_releases.tsv`, `vet_source_taxonomy.py`, `ast_calls.py` | retrieval of the public source at pinned release accessions and the rules that build the *Salmonella* example from it | `raw/` (not shipped), `examples/salmonella_poultry/data/` |

The four results folders in the table are not rerun by the confirmatory
campaign: the protocol records that the decomposition, the e-values and the
uniformity of the permutation p-value are unchanged in every path the defaults
take. Each folder carries a receipt naming the command and the commit that
wrote it and the sha256 of every file in it.

The truths of the confirmatory designs are written from their definitions in
code that shares nothing with the package, with one exception that the
protocol's wording does not spell out: the population value of ρ_min, the
truth of the rule "lower limit above ρ_min" of Study 3, is the minimum-norm
point of the package (`latent_order._lower`) at the population reading
probabilities of the design (`conditional/designs.py`, `_truth`). That
routine returns its value with a certified duality gap (it stops at 10⁻¹⁰,
`GAP_TOLERANCE`), and `tests/test_latent_order.py` checks it, gap at most
10⁻⁹, against a quadratic programme over the set of lineage means written out
subset by subset, which shares no code with it. ρ itself is computed by quadrature and shares nothing
with the package.

`scripts/reconcile_numbers.py` checks every number the article and its
supplement quote against the records in this directory and in `examples/`
(`scripts/reconcile_numbers.json` names the record behind each one).
