# 8 Environment and reproducibility

Preserve the software version, input hashes, configuration, seed, dependency versions and completion manifest of every run. Use a new destination for each run; a deliberate overwrite keeps a backup.

## Reproduce a workflow

Install the package in a separate environment and follow the four recipes in `examples/workflows`. Read method-specific retained cohorts and statuses before comparing numbers. Numerical agreement must use an explicitly defined comparison on matched inputs, targets and settings.

## Container recipe

The repository includes a Dockerfile with a pinned dependency recipe. Build it from the checkout with the local tag `amr-clonalshare:1.0.0`, verify the installed version, then run the documented workflows against mounted inputs and a fresh writable output directory. A local tag is a label, not a deposited release identifier.

The same image is built on every published release and pushed to the GitHub Container Registry as `ghcr.io/maciejkochanowski/amr-clonalshare:1.0.0`, with a second tag naming the commit it was built from. Pull it and run the shipped example without building anything:

```
docker pull ghcr.io/maciejkochanowski/amr-clonalshare:1.0.0
docker run --rm -v "$PWD":/data ghcr.io/maciejkochanowski/amr-clonalshare:1.0.0 \
    --config /data/config.yaml --results-dir /data/out
```

Container use does not by itself guarantee identical numerical-library behaviour across platforms. Compare the environment and the method-specific results before comparing numbers.

## Validation results

The repository-root `REPRODUCIBILITY.md` lists the stacks, seeds and commands behind the validation results in `benchmarks/`. The MIC interval calibration is described in `benchmarks/mic_inference/PROTOCOL.md`.
