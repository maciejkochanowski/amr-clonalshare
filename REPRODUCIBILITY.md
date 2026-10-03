# Reproduce amr-clonalshare 1.0.0

## Tests

```bash
python -m pip install -e '.[dev]'
pytest -q
```

## The software stack

Three files carry three different claims.

`pyproject.toml` states the versions the package supports. It uses ranges,
because a user installing the tool should not be forced onto one stack.

`requirements-lock.txt` states one exact stack: Python 3.11.15 with the pinned
versions of NumPy, SciPy, pandas and the other packages listed there. The
shipped example records are computed on this stack, and the validation
campaign runs in a virtual environment built from it.

`Dockerfile` builds that stack on a pinned Python base image, so a difference
between a rerun inside the image and a reported number is a difference in the
work rather than in the stack.

    docker build -t amr-clonalshare:1.0.0 .
    docker run --rm amr-clonalshare:1.0.0 --config examples/ssuis/config.yaml --results-dir /tmp/ssuis

    # or, without a container
    python -m venv .venv && . .venv/bin/activate
    pip install -r requirements-lock.txt
    pip install --no-deps .

Numbers are not promised to be bit-identical on other NumPy or SciPy versions.

## The validation campaign

The validation of record is the confirmatory campaign of
`benchmarks/CONFIRMATORY_PROTOCOL.md`. The protocol fixes, before any run, what
every study measures, on which designs, with how many replicates, from which
root seeds, and the rule that decides it; a result that fails its rule is
reported as a failure, and the claim it would have supported is restricted to
the designs that pass. The runs use the package at the commit that carries the
protocol, and a change to the package after the protocol voids the campaign for
the procedures it touches.

The campaign runs on a Slurm cluster from a campaign directory that holds
`source/` (a checkout of that commit), `env/` (a virtual environment built from
`requirements-lock.txt`) and `rlib/` (the R packages that
`benchmarks/competitors/install_r_bench.R` installs):

    AMR_CAMPAIGN_ROOT=/path/to/campaign benchmarks/campaign/launch_confirmatory.sh \
        ACCOUNT PARTITION "<command that sets up R>"

`benchmarks/campaign/make_confirmatory_queues.py` writes the queues: Studies 1
to 3 (the estimator grid of `estimator_benchmark.py`, the MIC ordering designs
of `order_calibration.py` and the fixed-distribution designs of `conditional/`) in one,
Study 4 (the published procedures of `competitors/`, run in R) in another, and
the merges, the summaries and the verdicts last, run once both have ended well.
Every replicate draws from its own stream of the study's root seed, and the
summary steps run under `benchmarks/campaign/run_logged.py`, which writes a
receipt naming the command, the commit, a digest of every Python file of the
package, the interpreter and the numerical stack, the Slurm job, the start, the
end and the exit status, and the sha256 of every file the step wrote.
`benchmarks/confirmatory_verdicts.py` applies the rules of the protocol to the
outputs and writes `verdicts.json` and `verdicts.csv`, and
`benchmarks/confirmatory_summaries.py` writes the per-design widths, means and
power beside them (`design_summaries.json`) with the SHA-256 of every raw
output file it read (`raw_outputs_sha256.json`). The verdicts of record, with
their receipts, are in `benchmarks/results_confirmatory/`; the raw outputs,
too large for the repository, are archived with the release.
`results_confirmatory/campaign/` lists the jobs and the commit each ran from:
Studies 1 to 4 from the commit that carries the protocol; procedure 5 again,
under the amendment of Section 5 of the protocol, by
`benchmarks/campaign/comparison_rerun.txt`; the undecided cells of Study 1 on
10,000 datasets (`estimator_benchmark.py --replicates 10000`, decided through
`confirmatory_verdicts.py --grid-rerun`); and, outside the protocol, the
follow-up of the one cell whose level rule was not met
(`results_confirmatory/followup/`).

The prevalence-difference decomposition, the e-values and the uniformity of the
permutation p-value are not rerun by that campaign (the protocol, Section 1);
their results folders (`results_decomposition_calibration`,
`results_decomposition_vs_regression`, `results_null_uniformity` and
`results_repeated_looks`) are kept, each with a receipt naming the command and
the commit that produced it.

## The examples

The records under `examples/*/expected*` are the output of the example
configurations at the seed of their published record (42; 20261001 for
*E. coli*), written with `amr-clonalshare --config <config> --results-dir <dir>`,
and of the comparison of `examples/ecoli_swine/input.csv` at the default seed
(`examples/ecoli_swine/README.md`).
`benchmarks/empirical/ssuis_analysis.py`, `ssuis_panel.py`,
`ssuis_sensitivity.py` and `salmonella_analysis.py` recompute the *S. suis* and
*Salmonella* analyses beyond the configured runs; each takes the folder it
writes to as its one argument and records the digests of its inputs and of the
package source beside its results. `benchmarks/campaign/empirical_queue.txt`
runs the eight configured analyses and the comparison, these four scripts, the seed-stability check
and the profile of a run, each under `run_logged.py`, in the campaign directory
of the confirmatory campaign; its results, with their receipts, are placed in
`benchmarks/results_empirical/` and the example records are copied from it.

## The numbers of the article

`scripts/reconcile_numbers.py` reads every number that the article and its
supplement quote from the record it comes from (`scripts/reconcile_numbers.json`
names the record, the field and the format of each) and checks the text
against it:

    python scripts/reconcile_numbers.py --manuscript manuscript.md --supplement supplement.md

It exits non-zero on a number that the records no longer give, so a rebuilt
article is checked in one step after the records are regenerated.

## Which artefact is the code of this release

The git tag `v1.0.0` is the source of the release. The Zenodo version record
deposited from it archives the repository as it stands and is the record to
cite for a reproduction; the concept DOI https://doi.org/10.5281/zenodo.22306353
groups all records. The wheel published on PyPI is built from the tag; its file
name carries a build number (`amr_clonalshare-1.0.0-<n>-py3-none-any.whl`),
because PyPI never lets a file name be reused, and pip installs the highest
build number. The container image `ghcr.io/maciejkochanowski/amr-clonalshare:1.0.0`
is built from the tag; cite it by the digest `docker pull` prints.

The repository is published as one commit. The git history of the work that
preceded it, which holds the commits named by the receipts, the protocol and
the appendices of the article, is attached to the release as `history.bundle`;
`git fetch history.bundle '+refs/heads/*:refs/history/*'` in a clone makes
those commits available, and `git bundle verify history.bundle` checks the file.
