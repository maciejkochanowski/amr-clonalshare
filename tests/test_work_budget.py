"""What the analysis costs, counted rather than timed.

A wall clock measures the machine as much as the code, so a timing test either
has a threshold so wide that it catches nothing or it fails on a busy runner.
What can be pinned exactly is the *work*: how many times the estimator scores
a held-out split. For the clonal share that number is a closed form in the
declared budgets,

    passes = (n_perm + 1) * min(repeats, null_repeats)
             + max(repeats - null_repeats, 0)
             + 2 * n_boot,

which says in one line what the resampling plan is: the permutation test
scores `null_repeats` fold draws for the observed labels and for each
permutation, repeats beyond `null_repeats` are scored once more for the
reported mean alone, the bootstrap costs two passes per replicate, and the
number of folds does not enter at all. A change that makes any of these
quadratic, or that quietly re-scores the observed draws inside the
permutation loop, fails here on any machine with no tolerance to argue about.

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


def _cohort():
    rng = np.random.default_rng(19)
    code = np.repeat(np.arange(8), 6)
    lineage = np.array([f"L{c}" for c in code], dtype=object)
    return (rng.random(code.size) < 0.4).astype(float), lineage


def _scoring_passes(monkeypatch, **budgets) -> int:
    counter = {"n": 0}

    def counted(*args, **kwargs):
        counter["n"] += 1
        return _SKILL(*args, **kwargs)

    monkeypatch.setattr(attribution, "_skill", counted)
    y, lineage = _cohort()
    attribution.clonal_share(y, lineage, seed=1, **budgets)
    return counter["n"]


def cost(*, repeats, n_perm, n_boot, null_repeats=5, **_ignored) -> int:
    return ((n_perm + 1) * min(repeats, null_repeats)
            + max(repeats - null_repeats, 0) + 2 * n_boot)


BUDGETS = [dict(folds=folds, repeats=repeats, n_perm=n_perm, n_boot=n_boot,
                null_repeats=null_repeats)
           for folds in (2, 3, 5)
           for repeats in (1, 2, 5, 7)
           for n_perm in (0, 9, 19)
           for n_boot in (0, 20)
           for null_repeats in (1, 5)]


@pytest.mark.parametrize("budget", BUDGETS,
                         ids=lambda b: "f{folds}r{repeats}p{n_perm}b{n_boot}n{null_repeats}"
                         .format(**b))
def test_the_scoring_passes_are_the_declared_cost(monkeypatch, budget):
    assert _scoring_passes(monkeypatch, **budget) == cost(**budget)


def test_the_number_of_folds_does_not_change_the_number_of_passes(monkeypatch):
    """One pass scores every fold of a draw; more folds is not more draws."""
    budget = dict(repeats=3, n_perm=9, n_boot=0, null_repeats=5)
    counts = {folds: _scoring_passes(monkeypatch, folds=folds, **budget)
              for folds in (2, 4, 8)}
    assert len(set(counts.values())) == 1, counts


def test_the_scoring_count_does_not_depend_on_the_seed(monkeypatch):
    budget = dict(folds=2, repeats=2, n_perm=19, n_boot=20)
    counts = {}
    for seed in (0, 1, 7):
        counter = {"n": 0}

        def counted(*args, __counter=counter, **kwargs):
            __counter["n"] += 1
            return _SKILL(*args, **kwargs)

        monkeypatch.setattr(attribution, "_skill", counted)
        y, lineage = _cohort()
        attribution.clonal_share(y, lineage, seed=seed, **budget)
        counts[seed] = counter["n"]
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
