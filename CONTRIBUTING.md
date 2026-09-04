# Contributing to amr-clonalshare

Thanks for your interest in improving the tool. This project follows the JOSS
community guidelines for support, reporting, and contribution.

## Getting support / asking questions

Open a **GitHub issue** with the `question` label, or start a discussion. Please
include your Python version, how you installed the package, and a minimal config
that reproduces what you are seeing.

## Reporting a bug

Open a GitHub issue with:

- what you expected vs. what happened;
- a **minimal reproducible example** — ideally a small config plus a few rows of
  synthetic data (please do **not** attach real surveillance data);
- the full traceback and your environment (`python --version`, `pip show amr-clonalshare`).

## Contributing code

1. Fork the repository and create a feature branch.
2. Set up a development environment:
   ```bash
   python -m venv .venv && source .venv/bin/activate
   pip install -e ".[dev]"
   ```
3. Add or update tests under `tests/`. Run them with:
   ```bash
   pytest -m "not slow"     # fast unit + smoke tests
   pytest                   # full suite incl. the synthetic cluster-recovery validation
   ruff check .             # dead code, unused arguments, undefined names
   mypy src/amr_clonalshare # the annotations are checked, not decorative
   ```
   CI runs the same four commands, then the per-module coverage floors in
   `scripts/check_coverage_floors.py`: a module that carries a reported
   number has a floor of its own, so the whole-package figure cannot hide
   an estimator that is barely exercised.
4. Keep statistical changes **provenance-annotated** — methodology lives in
   `docs/methodology.md`; update it alongside the code.
5. Open a pull request describing the change and its scientific rationale.

## Mutation testing

Coverage says a line ran; it does not say a test would notice if the line
were wrong. `mutmut run` (configured in `pyproject.toml`) rewrites one
operator or constant at a time in the estimator modules and reruns the tests
that cover it. A mutant that survives is either equivalent to the original,
which happens with a dtype that was already right or a default that equals
the value passed, or a gap in the tests. The last run on `stats.py` against
its direct tests (`tests/test_stats.py`, `tests/test_properties.py`) is
recorded below; rerun it after changing an estimator, and add a test for any
survivor that is not equivalent.

| date | module | mutants | killed | survived | not covered | score |
|---|---|---|---|---|---|---|
| 2026-09-04 | `stats.py` | 343 | 293 | 42 | 8 | 0.87 of the covered mutants |

Of the 42 survivors inspected on 2026-09-04, those read were equivalent
mutants: a `dtype=float` dropped from an array that was already float, a
`> 0` turned into `>= 0` inside a `where` whose other branch gives the same
value, a default argument removed where the caller passes that value. The
first pass of the same run, before `tests/test_properties.py` gained its
distance, step-up and tail cases, scored 0.70; the difference is what those
tests bought.

## Releasing

A release is a GitHub release, not a tag. Publishing a release archives the
tagged tree on Zenodo through the GitHub integration and runs
`.github/workflows/publish.yml`, which can also be started from the Actions
tab against the release tag: it builds the distribution, checks that the tag, `pyproject.toml` and
`__version__` agree, installs the wheel in a clean environment, and uploads to
PyPI through Trusted Publishing with a PEP 740 provenance attestation on every
file. No API token is stored anywhere.

The trust has to be declared once on PyPI, under the project's *Publishing*
settings: owner and repository as in `pyproject.toml`, workflow file
`publish.yml`, environment `pypi`. The `pypi` environment is created in the
repository settings, with required reviewers if a second pair of eyes is
wanted before an upload. Until both exist the workflow fails at the upload
step and nothing is published.

An installed release can be checked against its attestation:

```bash
pip download --no-deps amr-clonalshare==1.0.0
pypi-attestations verify pypi --repository https://github.com/maciejkochanowski/amr-clonalshare amr_clonalshare-1.0.0-py3-none-any.whl
```

## Code of conduct

Be respectful and constructive. Maintainers may remove comments or contributions
that are abusive or off-topic.
