# 8. Container and reproducibility

Two files make two separate claims. `pyproject.toml` states the version ranges
the software supports. `requirements-lock.txt` states the exact stack the
shipped example records were computed on, the release stack of
`REPRODUCIBILITY.md`, which also lists the stack every deposited artefact was
produced on. The `Dockerfile` installs that pinned stack, so that a difference
between a rerun and the shipped record is a difference in the work and not in
the stack.

## Build and check

```bash
docker build -t amr-clonalshare:1.0.0 .
docker run --rm amr-clonalshare:1.0.0 --version
docker run --rm amr-clonalshare:1.0.0 \
    --config examples/ssuis/config.yaml --check-input
```

The image runs as an unprivileged user, sets every thread limit to one so
that a run is deterministic, and carries the shipped cohort.

## Run your own data

Mount the directory that holds the configuration and the CSV files it names:

```bash
docker run --rm -v "$PWD":/data amr-clonalshare:1.0.0 \
    --config /data/config.yaml --results-dir /data/out --seed 42
```

Paths inside the configuration are relative to `data_dir`, so a configuration
that works on the host works unchanged inside the container.

## Reproduce the shipped run

```bash
mkdir -p out
docker run --rm -v "$PWD/out":/out amr-clonalshare:1.0.0 \
    --config examples/ssuis/config.yaml --results-dir /out/ssuis --quiet
```

`out/ssuis/clonal_share_result.json` carries the thirteen shares of
`examples/ssuis/expected/clonal_share_result.json`, and on the release stack
the two records agree to the last digit.
The continuous-integration workflow in `.github/workflows/ci.yml` builds this
image and runs the cohort inside it on every push, so the container is
tested, not only shipped.

## Without Docker

`REPRODUCIBILITY.md` in the repository root gives the same pinned environment
as a virtual environment:

```bash
python3.11 -m venv .venv && . .venv/bin/activate
pip install -r requirements-lock.txt
pip install --no-deps .
```

## What is and is not reproduced bit for bit

Every decision, every share, every interval endpoint to the precision reported
in the record reproduces across machines on the pinned stack. Across different
BLAS builds the eigenvalue step behind the effective number of independent
agents can differ in the last digits, and the
continuous-integration contract is therefore a set of properties (row-order
equivariance, control recovery, strict JSON) and not a byte comparison.
