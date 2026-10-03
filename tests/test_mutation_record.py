"""The record of surviving mutants, checked on every run of the suite.

A mutation run takes an hour and a half, so it is made by hand and its result
is written down. What can be checked in a second is the record itself: that
every class carries a verdict and a reason, that a class calling itself a gap
says what would close it, that no mutant is filed twice, and that the counts
in the file are the counts of the file. Without this the record would be free
to rot between runs.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

RECORD = Path(__file__).resolve().parent / "mutation_survivors.json"
VERDICTS = ("equivalent", "gap")


@pytest.fixture(scope="module")
def record():
    return json.loads(RECORD.read_text(encoding="utf-8"))


def test_every_class_carries_a_verdict_and_a_reason(record):
    for name, entry in record["classes"].items():
        assert entry["verdict"] in VERDICTS, name
        assert len(entry["reason"]) > 40, f"{name}: the reason is not a reason"
        if entry["verdict"] == "gap":
            assert len(entry["what_would_kill_it"]) > 20, name


def test_no_mutant_is_filed_in_two_classes(record):
    seen: dict[str, str] = {}
    for name, entry in record["classes"].items():
        for mutant in entry["mutants"]:
            assert mutant not in seen, f"{mutant}: {seen[mutant]} and {name}"
            seen[mutant] = name


def test_the_counts_are_the_counts_of_the_file(record):
    filed = 0
    for name, entry in record["classes"].items():
        assert entry["n"] == len(entry["mutants"]), name
        filed += entry["n"]
    assert filed == record["run"]["survived"]
    run = record["run"]
    assert run["killed"] + run["survived"] + run["timeout"] == run["mutants"]
    assert run["score"] == pytest.approx(run["killed"] / run["mutants"], abs=5e-5)


def test_every_mutant_names_a_module_the_run_measured(record):
    measured = set(record["run"]["modules"])
    for entry in record["classes"].values():
        for mutant in entry["mutants"]:
            parts = mutant.split(".")
            assert parts[0] == "amr_clonalshare", mutant
            assert parts[1] in measured, mutant


def test_a_claim_of_equivalence_is_argued_rather_than_asserted(record):
    """An equivalent mutant is one no test can kill; saying so needs a reason
    that explains why, not a label."""
    for name, entry in record["classes"].items():
        if entry["verdict"] != "equivalent":
            continue
        reason = entry["reason"].lower()
        assert ("same" in reason or "identical" in reason or "never" in reason), name
