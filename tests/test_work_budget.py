"""What the analysis costs, counted rather than timed.

A wall clock measures the machine as much as the code, so a timing test either
has a threshold so wide that it catches nothing or it fails on a busy runner.
What can be pinned exactly is the *work*: how many times the estimator scores
a held-out split, and how many times the bootstrap evaluates its closed form.
For the clonal share both are closed forms in the declared budgets,

    passes = repeats + n_perm,        fold averages = n_boot + n_perm,

which says in two lines what the resampling plan is: the estimate scores
`repeats` fold draws of the observed labels, each permuted labelling is
scored once for the penalty the estimate is corrected by (the p-value uses the
between-lineage sum of squares of the same labellings and scores nothing),
every bootstrap draw is scored once by the average of the share over the
dealings of the folds, which needs no split at all, each permuted labelling
is averaged so once more for the variance the share has without lineage
differences, and the number of folds does not enter. A change that makes any of these quadratic, that quietly
re-scores the observed draws inside the permutation loop, or that returns to
scoring every bootstrap draw by cross-validation, fails here on any machine
with no tolerance to argue about.

One coarse timing remains, with a very wide margin, because an implementation
can be made slower without being made to do more scoring passes.
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pytest

from amr_clonalshare import attribution

# Captured once: a second patch must wrap the package's function, not the
# counter that the first one installed.
_SKILL = attribution._skill
_FOLD_AVERAGE = attribution._fold_average


def _cohort():
    rng = np.random.default_rng(19)
    code = np.repeat(np.arange(8), 6)
    lineage = np.array([f"L{c}" for c in code], dtype=object)
    return (rng.random(code.size) < 0.4).astype(float), lineage


def _work(monkeypatch, seed=1, **budgets) -> tuple[int, int]:
    counter = {"skill": 0, "average": 0}

    def counted(*args, **kwargs):
        counter["skill"] += 1
        return _SKILL(*args, **kwargs)

    def averaged(*args, **kwargs):
        counter["average"] += 1
        return _FOLD_AVERAGE(*args, **kwargs)

    monkeypatch.setattr(attribution, "_skill", counted)
    monkeypatch.setattr(attribution, "_fold_average", averaged)
    y, lineage = _cohort()
    attribution.clonal_share(y, lineage, seed=seed, **budgets)
    return counter["skill"], counter["average"]


def cost(*, repeats, n_perm, n_boot, **_ignored) -> tuple[int, int]:
    # every bootstrap draw of this cohort varies (48 isolates, prevalence 0.4),
    # so every draw is scored; without permutations there is no control, no
    # corrected share and so no interval to draw
    return repeats + n_perm, n_boot + n_perm if n_perm and n_boot else 0


BUDGETS = [dict(folds=folds, repeats=repeats, n_perm=n_perm, n_boot=n_boot)
           for folds in (2, 3, 5)
           for repeats in (1, 2, 5, 7)
           for n_perm in (0, 9, 19)
           for n_boot in (0, 20)]


@pytest.mark.parametrize("budget", BUDGETS,
                         ids=lambda b: "f{folds}r{repeats}p{n_perm}b{n_boot}".format(**b))
def test_the_work_is_the_declared_cost(monkeypatch, budget):
    assert _work(monkeypatch, **budget) == cost(**budget)


def test_the_number_of_folds_does_not_change_the_work(monkeypatch):
    """One pass scores every fold of a draw; more folds is not more draws."""
    budget = dict(repeats=3, n_perm=9, n_boot=20)
    counts = {folds: _work(monkeypatch, folds=folds, **budget) for folds in (2, 4, 8)}
    assert len(set(counts.values())) == 1, counts


def test_the_work_does_not_depend_on_the_seed(monkeypatch):
    budget = dict(folds=2, repeats=2, n_perm=19, n_boot=20)
    counts = {seed: _work(monkeypatch, seed=seed, **budget) for seed in (0, 1, 7)}
    assert len(set(counts.values())) == 1, f"the plan depends on the draw: {counts}"


@pytest.mark.slow
def test_the_smallest_recipe_still_runs_in_seconds(tmp_path):
    """A margin of more than an order of magnitude over the measured two
    seconds: this catches an analysis that stopped being seconds and became
    hours, not a runner that is busy."""
    from amr_clonalshare import cli
    config = str(Path(__file__).resolve().parent / "golden" / "calls.yaml")
    start = time.perf_counter()
    assert cli.main(["--config", config, "--results-dir", str(tmp_path / "out"),
                     "--quiet"]) == 0
    assert time.perf_counter() - start < 60.0
