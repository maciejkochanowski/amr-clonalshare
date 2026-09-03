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
   ```
4. Keep statistical changes **provenance-annotated** — methodology lives in
   `docs/methodology.md`; update it alongside the code.
5. Open a pull request describing the change and its scientific rationale.

## Code of conduct

Be respectful and constructive. Maintainers may remove comments or contributions
that are abusive or off-topic.
