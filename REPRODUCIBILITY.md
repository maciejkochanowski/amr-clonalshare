# Reproduce amr-clonalshare 1.0.0

Install the package with its development extras and run the tests:

```bash
python -m pip install -e '.[dev]'
pytest -q
```

The reported results come from four steps. Each step checks its inputs against stored hashes.

1. S. suis and Salmonella examples: `python benchmarks/empirical/ssuis_analysis.py` and `python benchmarks/empirical/salmonella_analysis.py` (set `AMR_RESULTS_OUT=benchmarks/results_mic_release` to write beside the stored results).
2. MIC interval calibration: `benchmarks/mic_inference/calibrate.py` (exact readings, task table `tasks_confirmation_exact.tsv`) and `calibrate_null.py` (dilution panels, `tasks_null_confirmation.tsv`, root seed 20260916301). The Slurm scripts beside them record the settings used on the PLGrid Helios cluster. `summarize.py` and `summarize_null.py` verify every result file against its receipt before summarising.
3. S. suis MIC intervals: `prepare_empirical.py` writes the inputs and `run_empirical.py` analyses one agent (`helios_empirical.sbatch`).
4. Figures and tables of the MIC calibration are built from the files in `benchmarks/results_mic_release`.

Raw simulation rows, summaries, receipts and job accounting are in `benchmarks/results_mic_release`; the protocol, including the final numerical amendment, is `benchmarks/mic_inference/PROTOCOL.md`. The concept DOI for all versions is https://doi.org/10.5281/zenodo.22306353. Numbers are not promised to be bit-identical across different NumPy or SciPy versions.

## Result environments

Three files carry three different claims and they are deliberately not merged.

`pyproject.toml` states the versions the package supports. It uses ranges,
because a user installing the tool should not be forced onto one stack.

`requirements-lock.txt` states one exact stack, the release environment of
the table below, on which the shipped example records were computed; it is a
pin list read from the installed distribution metadata of that environment.
`uv.lock` is the resolver lock for the uv toolchain, generated from
`pyproject.toml`, and resolves the same ranges to the newest releases at the
time of locking rather than to the pin list; the pin list is what ran. The
stack each artefact was actually produced on is recorded in the artefact
itself and listed below.

`Dockerfile` builds that environment from the pin list on a pinned Python base
image. A difference between a rerun inside the image and a reported number
is therefore a difference in the work rather than in the stack.

## The stacks the results were computed on

Most artefacts record the interpreter and the numerical stack that produced
them in their own provenance block, and where one does not the table below is
the record. Four stacks appear across the repository and the campaigns run
outside it, and the table says which produced what.

| stack | interpreter | numpy | scipy | pandas | produced |
|---|---|---|---|---|---|
| release environment | 3.11.15 | 2.4.4 | 1.17.1 | 3.0.2 | the shipped example records and reports (`examples/*/expected*`), `benchmarks/results_ssuis_mechanism`, `benchmarks/results_ssuis_resolution`, `benchmarks/results_attribution_calibration`, the decomposition tree written by `scripts/run_decomposition_evidence.sh` |
| cluster | 3.11.16 | 2.4.6 | 1.17.1 | 3.0.5 | the 135-cell estimator grid of 2026-09-10 summarised in `src/amr_clonalshare/validation_grid.json`, the repeated-looks campaign, the comparator arms, the cut-off sensitivity envelope, the gate-conditional coverage study, the period split |
| cluster, atlas runs | 3.11.15 | 2.4.6 | 1.17.1 | 3.0.5 | the cross-species and veterinary atlases and the two-resolution atlas of 2026-09-09, the sequential repeated-looks arm |
| workstation | 3.12.13 | 2.5.2 | 1.18.1 | 3.0.5 | `benchmarks/results_null_uniformity_2026-09-04`, the realised-share calibration and the restricted-likelihood anchor of 2026-09-02 |

`requirements-lock.txt` pins the release stack, which is the one the container
image builds, so that the shipped example records reproduce inside it; `uv.lock`
resolves the same ranges for a developer checkout. The decomposition
calibration of 2026-08-31 was run on Python 3.11.5 with numpy 2.4.6, and the
receipt of `benchmarks/results_decomposition_vs_regression_2026-08-31/` records
the same Python 3.11.5 with numpy 2.4.6; neither receipt records a scipy or a
pandas version. Running the same analysis on more
than one stack is not redundancy: a failure that appears on one and not on
the other is an artefact of the machine, and reporting such an artefact as a
defect of the work is the most damaging error available in an audit. The
estimator grid was run on the cluster stack on 2026-09-02 and on the
workstation stack on 2026-09-08 with identical cell-by-cell results for every
estimator both runs hold, so the two stacks produce the same grid.

## Seeds

Every stochastic component is seeded and the seed is written into the
artefact. The estimator seed is 42 in the four atlas and period scripts
(`benchmarks/atlas_cross_species.py`, `vet_atlas.py`, `resolution_atlas.py`,
`period_split.py`), 20260902 in the estimator benchmark and the realised-share
calibration, 20260901 in the censored grid, the censored calibration and the
two *S. suis* arms, and 0 in the censored real-cohort run. The permutation
control carries a seed of its own, so that the real arm and the falsifying arm
cannot share a draw: 20260902 in the veterinary atlas, the two-resolution atlas
and the period split, and 20260901 in the cross-species atlas.

The reference estimator benchmark is the 135-cell
run of 2026-09-10 on the cluster stack, made after the species interval of
`clonal_share` was carried onto the component-ratio scale. Every quantity of
every other estimator, and the point estimate of `clonal_share`, agrees with
the run of 2026-09-09 to the ninth decimal; only the coverage and width of the
species interval moved, by at most 0.03 and 0.10 in any cell.

## Building and running

    docker build -t amr-clonalshare:1.0.0 .
    docker run --rm amr-clonalshare:1.0.0

    # or, without a container
    python -m venv .venv && . .venv/bin/activate
    pip install -r requirements-lock.txt
    pip install --no-deps .

## The S. suis decomposition arms

    bash scripts/run_decomposition_evidence.sh out/decomposition
    python benchmarks/ssuis_mechanism.py --out out/ssuis_mechanism
    python benchmarks/ssuis_resolution.py --out out/ssuis_resolution

The first runs the two input checks and the seven decomposition arms into one
tree with a receipt that records every command and the sha256 of every file
written; the decomposition arms reported in the article were made this way
and reproduce digit for digit. The other two write the per-agent
share against the determinant layer and at three definitions of a lineage.

## Which artefact is the code of this release

Three things carry the version 1.0.0 and one of them is authoritative. The git
tag `v1.0.0` is the source of the release, and the Zenodo concept DOI
https://doi.org/10.5281/zenodo.22306353 always resolves to the version record
deposited from this tag, the one that archives it as it stands, and that
version record is the one to cite for a reproduction;
the wheel published on PyPI is built from the tag and is byte-identical to it
module for module; its file name carries a build number
(`amr_clonalshare-1.0.0-<n>-py3-none-any.whl`), because PyPI never lets a file
name be reused. PyPI carries the wheel of this version and no source
distribution: the source-distribution file name of 1.0.0 is already taken on
the index, so the source of this release is the tag, the GitHub release and
the Zenodo archive. Install the wheel, the tag or the Zenodo archive to reproduce
any reported result. Every campaign records the commit it was run at.
