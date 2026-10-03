# Contributing to amr-clonalshare

Thanks for your interest in improving the tool. Support, bug reports and
contributions go through the GitHub repository as described below.

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
   pytest                   # full suite incl. the Monte-Carlo calibration tests
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

## How the package is tested

A test suite of a statistical package can pass while the package reports the
wrong number, because most of what a test usually asserts is that a function
returned something of the right shape. The layers below each answer a question
the others cannot. Run them all with `pytest`; the fast subset is
`pytest -m "not slow"`.

**Units and boundaries.** The ordinary tests, one per element of arithmetic,
plus the input boundaries: a table with a missing column, a call in an
unknown vocabulary, a collection of one.

**Analytic and high-precision oracles** (`tests/test_high_precision_oracles.py`).
Where a quantity has a definition that can be evaluated independently, the
test evaluates it from that definition and compares. The oracle never uses the
package's own formula: the two-sided Fisher test is the hypergeometric sum in
exact rational arithmetic, a Clopper-Pearson limit is checked by the tail
probability it makes equal to alpha/2 at forty digits with `mpmath`, the
closed form of the share averaged over the fold assignments is compared
with every fold assignment of a small collection enumerated in `fractions`, the share of a
smoothed distribution with its definition as a ratio of variances, the exact upper end
of the bounds on the latent ordering with the largest share over every
order of the lineages, and, for two lineages, the largest dual certificate
with the lower bound in closed form. This is the strongest evidence in
the suite that the reported numbers are computing what they say they compute.

**Properties** (`tests/test_properties.py`). Hypothesis states an identity or
an order relation and searches for a counterexample over thousands of inputs,
including the pathological ones a person does not think to write.

**Metamorphic relations** (`tests/test_metamorphic.py`). Nobody can say what
the lineage share of a fictional table ought to be, so no test can compare it
with a known answer. It can compare two runs: change the input in a way whose
effect is known in advance -- relabel the lineages, call the other outcome
positive, change the unit or the origin of the values, sort the rows
differently -- and
check that the answer changed in exactly that way, or not at all. This layer
found that the reported number depended on the order of the rows in the call
table and on the position of an antimicrobial among the columns; both are
fixed, and the relations are now tests.

**Data-centric mutation** (`tests/test_data_mutation.py`). Mutating the source
asks whether a test would notice a changed operator. Mutating the *data* asks
the question that matters for an analysis: a repeated row, a blanked label, an
unreadable call, one flipped outcome, a drug removed from the panel. Each
mutation declares the reaction it must provoke -- a refusal by name, no change
at all, or a change of a named count -- and a mutation that provokes no
reaction where one was declared fails the test.

**Golden artefacts** (`tests/test_golden_artefacts.py`). The whole of what the
user is given -- the record, the two input diagnostics, the result table and
both reports -- compared against a committed copy for four pipeline fixtures,
together with the output of the comparison of two lineage definitions and a
drafted configuration. A renamed key, a reworded heading or a moved number
fails here, and the diff of the golden file is the description of the change. Text is compared exactly; the
record's numbers to the precision a rerun reproduces; the recorded dependency
versions are excluded, because two machines write different strings there
while writing the same numbers.

**The wording of every refusal** (`tests/test_refusal_messages.py`). A refusal
is part of the interface: the person who gets it has to know which table,
which column and which value stopped the run. Every refusal of the measured
modules is called and its whole sentence compared, and a line-tracing check
fails if a refusal exists that no test reaches.

**Conditioning and sensitivity** (`tests/test_sensitivity.py`). Two
implementations agreeing to twelve digits says nothing about whether the
quantity is determined by the collection. These measure the Monte Carlo spread
against the width of the reported interval and the movement caused by
dropping any single isolate. They also fuzz the estimators with NaN, both infinities, denormals and blank labels
and require a named refusal or a result that declares itself not estimable --
never a NaN travelling on as a number.

**The cost** (`tests/test_work_budget.py`). Not a stopwatch, which measures the
machine: the number of scoring passes, which is a closed form in the declared
budgets. It pins the complexity as linear in the repeats, in the permutations
and in the bootstrap replicates, on any machine, with no tolerance to argue
about.

**Mutation of the source** (`mutmut`, below).

### Tolerances

`assert_allclose` and `pytest.approx` check `|a - b| <= atol + rtol * |b|`, so
`atol` decides near zero and `rtol` decides at scale. Write both, and write
beside them the reason for the number: the precision of the reference, the
conditioning of the problem, the error of the scheme. A tolerance loosened
until a test passed is a test that has stopped testing. The oracles in
`tests/test_high_precision_oracles.py` carry their derivations in the file.
One trap is worth naming: convert a float to a high-precision number exactly,
never through `repr`. The shortest decimal that round-trips is not the same
number as the float, and on an interval of width 1e-8 an endpoint moved by one
unit in the last place changes the width by one part in 10^8 -- an oracle
built that way reports an error the implementation does not have.

### Golden artefacts

Regenerate them with `python scripts/update_golden.py` after an intended
change, read the diff, and commit it together with the change that caused it.
Never regenerate to make a test pass.

### Benchmarks

`benchmarks/` holds the campaigns and their receipts. Pin the BLAS thread
count when running one -- `OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1` -- so that the timings are of the code and the last
digits do not move with the order of a threaded reduction.
`benchmarks/profile_run.py` reports where one analysis spends its time and its
memory, by phase and by function; `benchmarks/seed_stability.py` reports how
much of a reported ordering is the Monte Carlo draw.

## Mutation testing

Coverage says a line ran; it does not say a test would notice if the line
were wrong. `mutmut run` (configured in `pyproject.toml`) rewrites one
operator, constant or string at a time and reruns the tests that cover it.

**Where it is applied.** The twelve modules in `source_paths`: `stats.py`,
`clonality.py`, `attribution.py`, `evalues.py`, `missingness.py`, `qc.py`,
`phenotype.py`, `comparison.py`, `draft.py`, `censored.py`, `mic_order.py` and
`latent_order.py`. The rule is two conditions together: the module's
arithmetic carries a number the report prints, *and* its tests can kill a
mutant in seconds.

**Where it is not, and why.** The presentation layer (`report_model.py`,
`report_html.py`, `report.py`, `outputs.py`) is excluded: the golden artefacts
compare its entire output for every fixture, which kills the mutants a
mutation run would find there and does it in seconds.

A mutant that survives is either equivalent to the original, which happens
with a dtype that was already right or a default that equals the value passed,
or a gap in the tests. Equivalent mutants cannot be killed by any test, since
there is nothing to observe; deciding which are equivalent is in general
undecidable, so the survivors are listed with a reason rather than counted and
forgotten. Every survivor of the recorded run is written down in
`tests/mutation_survivors.json`, by class, each class with a verdict and a
reason. A class called equivalent computes what the original computes, exactly
or within the rounding or the stopping tolerance the computation already
admits -- a change in the path of the Newton iteration of `mic_order.py` that
stops at the same maximum, or of the minimum-norm-point algorithm of
`latent_order.py` that stops at a lower end certified to the same tolerance,
is of that kind -- so no test of a reported quantity can kill it, and none
should be written to try. A class called a gap is one a test could kill, and
the record says what that test would have to do. Calling a class a gap rather
than equivalent is deliberate when equivalence is probable but cannot be shown
without reading every member: a class that claims equivalence without an
argument for every member is worth less than an honest list.

The record names the modules its run measured (`run.modules`) and the counts
it found. It holds the run of commit b9a9f50 on the twelve modules above:
9,792 of 11,345 mutants killed (86.3%), 1,543 survivors, 221 of them in
classes argued equivalent. Every class also counts the members that change a
line the value records of `tests/test_golden_values.py` execute
(`n_on_lines_the_value_records_run`, 874 in all): having survived, such a
mutant changes none of the values those records compare, so the gap it names
is an input the fixtures do not contain, not a value a result reports
differently on them.

A run of all mutants at once can kill a mutant through a test that fails for a
reason outside it, so two runs of the same source and tests can differ in a
few mutants. `check_mutation_survivors.py` reports such a mutant as a stale
entry; running it alone (`mutmut run <name>`) decides it.

Survivors in new code are read, not only counted. An earlier run, of commit
656bade, left gaps in newly added code that computed something different:
settings no test passed through, an interval assembled from parts no test
compared, a lower confidence limit no test recomputed. Its tests were written
against those mutants, and reading the comparison of two lineage definitions
for them found that each analysis's `id_sha256` hashed the labels of the
second definition instead of the identifiers. The run of commit 9877288 left
survivors in the certificate of the lower confidence limit whose variance and
degrees of freedom no test recomputed, and in the branch and bound of the upper
end, whose leaves no fixture reached; tests now compare the certificate with
its written form isolate by isolate and with the spread of its estimate over
repeated datasets, run the branch and bound from no incumbent, and pin the
limit where it is positive in the value records. The run of commit 84c69a3 left
two such gaps in code added since, and tests now close both: the certificate is
compared with its written form on readings that overlap within a stratum, and a
within-lineage component that the shared-support gate refuses is carried into the
selection of a panel and shown never to be selected.

`scripts/check_mutation_survivors.py` compares a fresh run against that
record and fails on a survivor it does not name, on an entry whose mutant no
longer survives, and on a class without a verdict or a reason.
`tests/test_mutation_record.py` checks the record itself on every run of the
suite, so it cannot rot between mutation runs.

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

Before the release is published, `scripts/verify_sdist_cleanroom.py` checks
that the source distribution in `dist/` is self-sufficient: it extracts the
sdist, builds a Python 3.11 environment from the shipped `requirements-lock.txt`
so the cleanroom is the release stack, installs the sdist with its `dev` extra
into it and runs the whole test suite from the extracted tree, the slow tests
and the high-precision oracles included, writing
`dist/SDIST_CLEANROOM_RECEIPT.json`. It needs `uv` on the path and network
access to fetch the pinned distributions; run it after `python -m build` and
before publishing the release.

An installed release can be checked against its attestation:

```bash
pip download --no-deps amr-clonalshare==1.0.0
pypi-attestations verify pypi --repository https://github.com/maciejkochanowski/amr-clonalshare amr_clonalshare-1.0.0-*-py3-none-any.whl   # the wheel carries a build number, e.g. -21-
```

## Code of conduct

Be respectful and constructive. Maintainers may remove comments or contributions
that are abusive or off-topic.
