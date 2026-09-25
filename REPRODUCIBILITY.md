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
versions of NumPy, SciPy, pandas and the other packages listed there. Every
reported result, every validation file shipped in the package and every shipped
example record was computed on this stack.

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

Every result under `benchmarks/results_*` and the three validation files the
package reads at run time (`validation_grid.json`,
`population_probit_validation.json`, `general_probit_protocol.json`) were
produced by one campaign, run from one commit of this repository on the stack
above. `benchmarks/campaign/CAMPAIGN.md` lists its commands in order. Each step
wrote a receipt naming the command, the commit, a digest of every Python file of
the package, the library versions, the Slurm job and the sha256 of every file it
wrote; the per-task receipts of the MIC simulations also carry the hash of every
source file and are checked by `summarize.py` and `summarize_null.py` before
anything is summarised. `benchmarks/results_mic_release/accounting.json` records
the CPU-hours of every job of the campaign.

The protocols are `benchmarks/mic_inference/PROTOCOL.md` (MIC interval) and
`benchmarks/population_model/PROTOCOL.md` (population model). Every seed is
fixed in the driver that uses it and written into the result.

## The examples

The shipped records under `examples/*/expected*` are the output of the three
example configurations at the default seed; `benchmarks/empirical/ssuis_analysis.py`
and `benchmarks/empirical/salmonella_analysis.py` write the *S. suis* and
*Salmonella* analyses the article reports (set `AMR_RESULTS_OUT` to write
elsewhere). Each checks its inputs against stored hashes.

## Which artefact is the code of this release

The git tag `v1.0.0` is the source of the release. The Zenodo version record
deposited from it archives the repository as it stands and is the record to
cite for a reproduction; the concept DOI https://doi.org/10.5281/zenodo.22306353
groups all records. The wheel published on PyPI is built from the tag; its file
name carries a build number (`amr_clonalshare-1.0.0-<n>-py3-none-any.whl`),
because PyPI never lets a file name be reused, and pip installs the highest
build number. The container image `ghcr.io/maciejkochanowski/amr-clonalshare:1.0.0`
is built from the tag; cite it by the digest `docker pull` prints.
