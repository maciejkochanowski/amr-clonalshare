# 1 Install

## Install software 1.0.0

From PyPI or from a source checkout:

```bash
python -m pip install amr-clonalshare      # or: python -m pip install .
amr-clonalshare --version
```

Check that the version is `1.0.0`. Python 3.11 or later and NumPy, pandas, SciPy and PyYAML are required. Use a separate environment for each project.

For development tools and tests, install `python -m pip install -e ".[dev]"` from the source checkout. Software tests and statistical validation answer different questions; keep the interpreter and dependency versions with any result.

## Threads

For a command-line run, `--threads 1` requests a single numerical-library thread before the numerical modules are imported:

```bash
amr-clonalshare --config config.yaml --results-dir out --threads 1
```

It also caps the worker processes of the calibrated MIC interval at the same number unless `censored.workers` sets them. Thread settings alone do not guarantee byte-identical output across operating systems or numerical libraries. Use the same recorded environment, source, inputs, seed and settings for a numerical reproduction. See the [reproduction guide](08-container.md).
