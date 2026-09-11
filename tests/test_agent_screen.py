"""The screening rule, and the accounting identity it exists to protect.

Both atlases used to pass over an agent that failed the tested-isolate
threshold without writing a row. The counts in the shipped files then did not
reconcile: 366 of 907 agents in the veterinary run and 21 in the cross-species
run left no trace, and a reader had no way to tell that apart from selection.
These tests hold the replacement to account.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "benchmarks"))

from agent_screen import (
    MIN_DRUG_N, MIN_LINEAGES, MIN_MINOR, REFUSAL_NO_VARIANCE,
    REFUSAL_TOO_FEW_LINEAGES, REFUSAL_TOO_FEW_TESTED, REFUSALS, reconcile,
    screen_agent)


def _cohort(n_tested, prevalence, n_lineages, n_total=400, seed=0):
    """A cohort with the three quantities the screen reads set on purpose."""
    rng = np.random.default_rng(seed)
    y = np.full(n_total, np.nan)
    k = max(1, int(round(prevalence * n_tested)))
    calls = np.zeros(n_tested)
    calls[:k] = 1.0
    y[:n_tested] = rng.permutation(calls)
    lineage = [f"L{i % max(n_lineages, 1)}" for i in range(n_total)]
    return y, lineage


def test_every_agent_leaves_exactly_one_outcome():
    """The property the accounting rests on: no input reaches a fourth path.

    The screen either analyses or refuses with a named reason, for every
    combination of the three quantities it reads. A path that returned nothing
    is what produced the unreconciled counts, so its absence is asserted
    directly rather than inferred from the totals.
    """
    for n_tested in (0, 1, MIN_DRUG_N - 1, MIN_DRUG_N, 200):
        for prevalence in (0.0, MIN_MINOR / 2, 0.5, 1.0 - MIN_MINOR / 2, 1.0):
            for n_lineages in (1, 2, 40):
                y, lineage = _cohort(n_tested, prevalence, n_lineages)
                screen = screen_agent(y, lineage)
                assert screen.analyse or screen.reason in REFUSALS
                assert screen.analyse != bool(screen.record.get("skipped"))


def test_the_refusals_name_the_condition_that_failed():
    y, lineage = _cohort(MIN_DRUG_N - 1, 0.5, 40)
    assert screen_agent(y, lineage).reason == REFUSAL_TOO_FEW_TESTED

    y, lineage = _cohort(200, 0.0, 40)
    assert screen_agent(y, lineage).reason == REFUSAL_NO_VARIANCE

    y, lineage = _cohort(200, 1.0, 40)
    assert screen_agent(y, lineage).reason == REFUSAL_NO_VARIANCE

    y, lineage = _cohort(200, 0.5, 1)
    screen = screen_agent(y, lineage)
    assert screen.reason == REFUSAL_TOO_FEW_LINEAGES
    assert screen.record["n_lineages_analysed"] < MIN_LINEAGES


def test_an_analysed_agent_carries_the_subset_it_was_screened_on():
    """The analysed subset must be the tested isolates and their own labels.

    Screening on the full vector and estimating on a differently built subset
    is how a permutation control stops being comparable to the cohort it
    controls for.
    """
    y, lineage = _cohort(200, 0.5, 40, n_total=400)
    screen = screen_agent(y, lineage)
    assert screen.analyse
    assert screen.y.size == 200
    assert len(screen.lineage) == 200
    assert np.isfinite(screen.y).all()
    assert list(screen.lineage) == lineage[:200]
    assert screen.record["n"] == 200


def test_the_reconciliation_refuses_a_cohort_that_loses_an_agent():
    """The identity is raised, not reported: a run that cannot account for its
    own units has nothing to report."""
    records = {"a": {"n": 60, "prevalence": 0.4},
               "b": {"n": 10, "skipped": REFUSAL_TOO_FEW_TESTED}}
    reconcile(2, records)

    with pytest.raises(AssertionError, match="every agent must leave one"):
        reconcile(3, records)

    with pytest.raises(AssertionError, match="outside the closed set"):
        reconcile(2, {"a": {"skipped": "because"},
                      "b": {"skipped": REFUSAL_NO_VARIANCE}})


def test_mismatched_lengths_are_refused_rather_than_broadcast():
    with pytest.raises(ValueError, match="must describe the same isolates"):
        screen_agent(np.zeros(10), ["L0"] * 9)
