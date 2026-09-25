# 8 Environment and reproducibility

Preserve the software version, input hashes, configuration, seed, dependency versions and completion manifest of every run. Use a new destination for each run; `--overwrite` replaces the earlier results.

## Reproduce a workflow

Install the package in a separate environment and follow the four recipes in `examples/workflows`. Read method-specific retained cohorts and statuses before comparing numbers. Numerical agreement must use an explicitly defined comparison on matched inputs, targets and settings.

## Container recipe

The repository includes a Dockerfile with a pinned dependency recipe. Build it from the checkout with the local tag `amr-clonalshare:1.0.0`, verify the installed version, then run the documented workflows against mounted inputs and a fresh writable output directory. A local tag is a label, not a deposited release identifier.

The same image is built on every published release and pushed to the GitHub Container Registry as `ghcr.io/maciejkochanowski/amr-clonalshare:1.0.0`, with a second tag naming the commit it was built from; cite the image by the digest `docker pull` prints. The image carries the shipped example, which runs without building or mounting anything:

```
docker pull ghcr.io/maciejkochanowski/amr-clonalshare:1.0.0
docker run --rm ghcr.io/maciejkochanowski/amr-clonalshare:1.0.0 \
    --config examples/ssuis/config.yaml --results-dir /tmp/ssuis
```

For your own data, mount the folder that holds your configuration and tables and run as your own user, so that the results written back into it are yours and the container can write there:

```
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/data \
    ghcr.io/maciejkochanowski/amr-clonalshare:1.0.0 \
    --config /data/config.yaml --results-dir /data/out
```

Container use does not by itself guarantee identical numerical-library behaviour across platforms. Compare the environment and the method-specific results before comparing numbers.

## What a result depends on

A record is a function of the data and the seed. It does not depend on the
order of the rows in any input table, on the position of an antimicrobial
among the columns, on how many antimicrobials the file holds beside it, or on
the number of worker processes: sorting a spreadsheet differently, or adding
one more drug to it, leaves every other number where it was. The package
orders the panel by identifier and by agent name before anything reads it, and
each agent draws from a fresh generator on the same prespecified stream. The
relations are tests (`tests/test_metamorphic.py`), so they cannot quietly stop
holding.

Two things are not promised. A different version of NumPy, SciPy or pandas may
move the last digits of a floating-point result, which is why the record names
the versions it was produced with and the report prints a digest that covers
them; in the release's own checks the six artefacts of every fixture were
identical across two such environments, but that is a measurement and not a
guarantee. And the Monte Carlo error of the estimators themselves remains: a
different seed gives a different draw, with a spread the release measures at
0.002 to 0.007 on the shipped example, two orders below the width of the
reported interval (`benchmarks/results_seed_stability`).

## Validation results

The repository-root `REPRODUCIBILITY.md` lists the stacks, seeds and commands behind the validation results in `benchmarks/`. The MIC interval calibration is described in `benchmarks/mic_inference/PROTOCOL.md`.
