"""Mutating the data instead of the code.

A mutation of the source asks whether a test would notice a changed operator.
Most of those mutants are equivalent, because a statistical program spends
much of its text on bookkeeping that no input can distinguish. A mutation of
the *data* asks the question that matters here: change one record in a way
whose consequence is known -- a repeated row, a blanked label, an unreadable
call, one flipped outcome, a drug removed from the panel -- and check that the
run reacts exactly as declared. A mutant that provokes no reaction where one
was declared is a surviving mutant and fails the test, whether the reaction
that is missing is a refusal, a changed count or a number that should have
moved.
"""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

import pytest

N = 40
LINEAGES = [f"L{i // 4 + 1:02d}" for i in range(N)]
IDS = [f"iso_{i:02d}" for i in range(N)]
# ag_varies has both outcomes, ag_constant has one, ag_unread was never
# interpreted.
PATTERN = {"ag_varies": "RSRSS", "ag_constant": "SSSSS", "ag_unread": "NTNTN"}
CONFIG = """\
dataset:
  name: data_mutation_fixture
  data_dir: data
  metadata: metadata.csv
  strain_id_column: isolate_id
  lineage_column: lineage
  phenotype: calls.csv
  phenotype_id_column: isolate_id
  phenotype_antibiotic_column: agent
  phenotype_call_column: call
  phenotype_kind: clinical_sir
attribution: {folds: 2, repeats: 2, n_boot: 20, n_perm: 19}
surveillance: {n_boot: 20}
evidence: {folds: 2, repeats: 2}
"""


def baseline_tables():
    metadata = [[IDS[i], LINEAGES[i]] for i in range(N)]
    calls = [[IDS[i], agent, PATTERN[agent][i % 5]]
             for agent in sorted(PATTERN) for i in range(N)]
    return calls, metadata


def write(root: Path, calls, metadata) -> Path:
    data = root / "data"
    data.mkdir(parents=True, exist_ok=True)
    with (data / "metadata.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["isolate_id", "lineage"])
        writer.writerows(metadata)
    with (data / "calls.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["isolate_id", "agent", "call"])
        writer.writerows(calls)
    config = root / "config.yaml"
    config.write_text(CONFIG, encoding="utf-8")
    return config


def run(root: Path, calls, metadata) -> dict:
    from amr_clonalshare import cli
    out = root / "out"
    assert cli.main(["--config", str(write(root, calls, metadata)),
                     "--results-dir", str(out), "--quiet"]) == 0
    return json.loads((out / "clonal_share_result.json").read_text(encoding="utf-8"))


def refuse(root: Path, calls, metadata) -> str:
    """Run a table the loader must refuse; return what it told the user."""
    from amr_clonalshare import cli
    import contextlib
    import io

    captured = io.StringIO()
    with contextlib.redirect_stderr(captured):
        status = cli.main(["--config", str(write(root, calls, metadata)),
                           "--results-dir", str(root / "out"), "--quiet"])
    assert status == 2, f"the run ended {status} instead of refusing the table"
    return captured.getvalue()


def reported(record: dict, agent: str, analysis: str) -> dict:
    return next(row for row in record["results"]
                if row["agent"] == agent and row["analysis"] == analysis)


def numbers(record: dict) -> dict:
    diagnostics = record["metadata_diagnostics"]
    out = {}
    for block_name in ("clonal_share", "realised_share"):
        for agent, block in diagnostics.get(block_name, {}).items():
            out[f"{block_name}/{agent}"] = {k: v for k, v in block.items()
                                            if k != "reason"}
    for agent, block in diagnostics["lineage_evidence"]["per_feature"].items():
        out[f"evidence/{agent}"] = dict(block)
    return out


def reading(record: dict) -> dict:
    return record["input_qc"]["phenotype_reading"]


# --------------------------------------------------------------------------- #
# The mutations and the reaction each one must provoke
# --------------------------------------------------------------------------- #

@dataclass
class DataMutant:
    name: str
    mutate: Callable[[list, list], tuple]
    reaction: str                      # "refused", "inert" or "reacts"
    check: Optional[Callable[[dict, dict], None]] = None
    refusal: Optional[str] = None


def _repeat_identical(calls, metadata):
    return calls + [list(calls[0])], metadata


def _repeat_conflicting(calls, metadata):
    row = list(calls[0])
    row[2] = "S" if row[2] == "R" else "R"
    return calls + [row], metadata


def _blank_one_lineage(calls, metadata):
    metadata = [list(row) for row in metadata]
    metadata[3][1] = ""
    return calls, metadata


def _blank_one_call(calls, metadata):
    calls = [list(row) for row in calls]
    index = next(i for i, row in enumerate(calls) if row[1] == "ag_varies")
    calls[index][2] = ""
    return calls, metadata


def _make_one_call_unreadable(calls, metadata):
    calls = [list(row) for row in calls]
    index = next(i for i, row in enumerate(calls) if row[1] == "ag_varies")
    calls[index][2] = "ND"
    return calls, metadata


def _flip_one_call(calls, metadata):
    calls = [list(row) for row in calls]
    index = next(i for i, row in enumerate(calls)
                 if row[1] == "ag_varies" and row[2] == "S")
    calls[index][2] = "R"
    return calls, metadata


def _shuffle(calls, metadata):
    import numpy as np
    calls = [list(row) for row in calls]
    np.random.default_rng(23).shuffle(calls)
    return calls, metadata


def _drop_an_agent(calls, metadata):
    return [row for row in calls if row[1] != "ag_constant"], metadata


def _add_an_isolate_with_nothing(calls, metadata):
    return calls, metadata + [["iso_99", ""]]


def _one_lineage_only(calls, metadata):
    return calls, [[row[0], "L01"] for row in metadata]


def _exchange_two_labels(calls, metadata):
    metadata = [list(row) for row in metadata]
    metadata[0][1], metadata[-1][1] = metadata[-1][1], metadata[0][1]
    return calls, metadata


def _missing_metadata_row(calls, metadata):
    return calls, metadata[:-1]


def _check_untyped(before, after):
    key = "clonal_share/ag_varies"
    assert numbers(after)[key]["n_dropped_untyped"] == \
        numbers(before)[key]["n_dropped_untyped"] + 1


def _check_missing_call(before, after):
    assert reading(after)["n_missing_calls"] == reading(before)["n_missing_calls"] + 1


def _check_unreadable_call(before, after):
    assert reading(after)["n_unrecognized_calls"] == \
        reading(before)["n_unrecognized_calls"] + 1
    assert "nd" in reading(after)["unrecognized_calls"]


def _check_one_more_positive(before, after):
    key = "clonal_share/ag_varies"
    n = numbers(before)[key]["n"]
    assert numbers(after)[key]["prevalence"] == pytest.approx(
        numbers(before)[key]["prevalence"] + 1.0 / n, rel=0, abs=1e-12)


def _check_not_estimable(before, after):
    assert after["metadata_diagnostics"]["clonal_share"]["ag_varies"]["estimable"] is False
    row = reported(after, "ag_varies", "collection_membership")
    assert row["status"] == "unavailable"
    assert row["estimate"] is None
    assert row["reason"] == "Only one retained lineage; no between-lineage comparison"


def _check_join(before, after):
    for key in ("n_joined", "strains_joined"):
        assert after["input_qc"]["metadata_join"][key] == \
            before["input_qc"]["metadata_join"][key] - 1
    assert after["input_qc"]["metadata_join"]["share_joined"] < 1.0


MUTANTS = [
    DataMutant("a repeated row that agrees", _repeat_identical, "inert"),
    DataMutant("a repeated row that disagrees", _repeat_conflicting, "refused",
               refusal="has conflicting records for 1 key(s)"),
    DataMutant("one lineage label blanked", _blank_one_lineage, "reacts",
               check=_check_untyped),
    DataMutant("one call blanked", _blank_one_call, "reacts",
               check=_check_missing_call),
    DataMutant("one call in another vocabulary", _make_one_call_unreadable,
               "reacts", check=_check_unreadable_call),
    DataMutant("one call flipped to the rarer outcome", _flip_one_call, "reacts",
               check=_check_one_more_positive),
    DataMutant("the rows in another order", _shuffle, "inert"),
    DataMutant("an antimicrobial removed from the panel", _drop_an_agent, "inert"),
    DataMutant("an isolate with no label and no call", _add_an_isolate_with_nothing,
               "inert"),
    DataMutant("every isolate in one lineage", _one_lineage_only, "reacts",
               check=_check_not_estimable),
    DataMutant("two labels exchanged", _exchange_two_labels, "reacts"),
    DataMutant("one metadata row missing", _missing_metadata_row, "reacts",
               check=_check_join),
]


@pytest.fixture(scope="module")
def baseline(tmp_path_factory):
    root = tmp_path_factory.mktemp("baseline")
    return run(root, *baseline_tables())


@pytest.mark.parametrize("mutant", MUTANTS, ids=[m.name for m in MUTANTS])
def test_the_data_mutant_provokes_the_declared_reaction(mutant, baseline, tmp_path):
    calls, metadata = mutant.mutate(*baseline_tables())
    if mutant.reaction == "refused":
        assert mutant.refusal in refuse(tmp_path, calls, metadata)
        return
    after = run(tmp_path, calls, metadata)
    before_numbers = numbers(baseline)
    after_numbers = numbers(after)
    if mutant.reaction == "inert":
        # A mutation that removes an antimicrobial removes its rows; what is
        # still reported must be reported unchanged.
        shared = set(before_numbers) & set(after_numbers)
        assert shared, "the mutation left nothing to compare"
        assert {k: after_numbers[k] for k in shared} == \
            {k: before_numbers[k] for k in shared}, (
            "this change of the data may not move a reported number")
        return
    moved = [key for key in before_numbers
             if after_numbers.get(key) != before_numbers[key]]
    assert moved, (
        "the run did not react to this change of the data at all; either the "
        "reaction is missing or the mutation is declared wrongly")
    if mutant.check is not None:
        mutant.check(baseline, after)


def test_every_declared_reaction_is_one_the_harness_can_check():
    assert {m.reaction for m in MUTANTS} <= {"refused", "inert", "reacts"}
    assert all(m.refusal for m in MUTANTS if m.reaction == "refused")
    assert all(m.check is None or m.reaction != "refused" for m in MUTANTS)
