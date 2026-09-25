"""Every refusal, by the sentence it prints.

A refusal is part of the interface: the person who gets it has to know which
table, which column and which value stopped the run. These tests call each
refusal and compare the whole sentence, so a reworded, truncated or emptied
message fails here rather than reaching a user.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from amr_clonalshare import clonality, comparison, evalues, missingness, qc, stats
from amr_clonalshare.attribution import clonal_share, layer_clonal_share
from amr_clonalshare.draft import draft_config
from amr_clonalshare.phenotype import resolve_duplicates, to_non_susceptible
from amr_clonalshare.realised import realised_share


def _long(rows):
    return pd.DataFrame(rows, columns=["Strain_ID", "antibiotic", "resistant_phenotype"])


def _panel(agents=("a",), n=4):
    rows = []
    for i in range(n):
        for agent in agents:
            rows.append((f"i{i}", agent, "Resistant" if i % 2 else "Susceptible"))
    return _long(rows)


CASES = [
    ("stats/q",
     lambda: stats.benjamini_hochberg([0.1], q=0.0),
     "q must satisfy 0 < q < 1"),
    ("stats/nan_policy",
     lambda: stats.benjamini_hochberg([0.1], nan_policy="skip"),
     "nan_policy must be raise or omit"),
    ("stats/dependence",
     lambda: stats.benjamini_hochberg([0.1], dependence="positive"),
     "dependence must be independent or arbitrary"),
    ("stats/shape",
     lambda: stats.benjamini_hochberg([[0.1, 0.2]]),
     "p-values must be a one-dimensional vector"),
    ("stats/range",
     lambda: stats.benjamini_hochberg([1.5]),
     "finite p-values must lie in [0, 1]"),
    ("stats/nonfinite",
     lambda: stats.benjamini_hochberg([0.1, float("nan")]),
     "1 of 2 p-values are not finite (indices [1]); Benjamini-Hochberg would "
     "report no discoveries for the whole family. Fix the tests that produced "
     "them, or pass nan_policy='omit' and report how many were dropped"),
    ("stats/permutation_finite",
     lambda: stats.permutation_pvalue([1.0, float("inf")], 0.5),
     "permutation statistics and observed statistic must be finite"),
    ("stats/tail",
     lambda: stats.permutation_pvalue([1.0], 0.5, tail="two-sided"),
     "tail must be greater or less"),
    ("clonality/length",
     lambda: clonality.decompose_prevalence_difference([0, 1], ["A", "B"], [0], ["A", "B"]),
     "y and lineage must have the same length"),
    ("clonality/agents",
     lambda: clonality.decompose_panel(pd.DataFrame({"a": [0, 1]}), ["L", "L"],
                                       pd.DataFrame({"b": [0, 1]}), ["L", "L"]),
     "both collections must carry the same agent columns"),
    ("attribution/folds",
     lambda: clonal_share([0, 1, 0, 1], ["A", "A", "B", "B"], folds=1),
     "folds must be an integer >= 2; at least two folds are required"),
    ("attribution/repeats",
     lambda: clonal_share([0, 1, 0, 1], ["A", "A", "B", "B"], repeats=0),
     "repeats must be an integer >= 1"),
    ("attribution/rows",
     lambda: layer_clonal_share(np.zeros((2, 5)), ["A", "B", "C"]),
     "X has 2 rows and lineage has 3 entries; they must match"),
    ("realised/alpha",
     lambda: realised_share([0.0, 1.0], ["A", "A"], alpha=0.0),
     "alpha=0.0; the interval needs 0 < alpha < 1"),
    ("realised/length",
     lambda: realised_share([0.0, 1.0], ["A"]),
     "y has 2 entries and lineage has 1; they must match"),
    ("evalues/folds",
     lambda: evalues.e_process([0, 1], ["A", "B"], folds=1),
     "folds must be an integer >= 2"),
    ("evalues/repeats",
     lambda: evalues.e_process([0, 1], ["A", "B"], repeats=0),
     "repeats must be an integer >= 1"),
    ("evalues/length",
     lambda: evalues.e_process([0, 1, 0], ["A", "B"]),
     "lineage has 2 entries for 3 traits"),
    ("evalues/binary",
     lambda: evalues.e_process([0, 2], ["A", "B"]),
     "evidence requires binary outcomes (0 or 1); non-finite values are set aside"),
    ("evalues/batches",
     lambda: evalues.sequential_e_process([[0, 1], [0, 1]], [["A", "B"]]),
     "2 trait batches for 1 lineage batches"),
    ("evalues/batch_length",
     lambda: evalues.sequential_e_process([[0, 1, 0]], [["A", "B"]]),
     "batch 0: 3 traits for 2 labels"),
    ("evalues/batch_binary",
     lambda: evalues.sequential_e_process([[0, 1], [0, 3]], [["A", "B"], ["A", "B"]]),
     "batch 1: evidence requires binary outcomes (0 or 1); non-finite values are set aside"),
    ("evalues/combine_negative",
     lambda: evalues.combine_independent([1.0, -1.0]),
     "e-values must be non-negative and not NaN"),
    ("evalues/combine_zero_times_inf",
     lambda: evalues.combine_independent([0.0, float("inf")]),
     "an infinite e-value multiplied by a zero one is undefined"),
    ("evalues/cohort_empty",
     lambda: evalues.combine_within_cohort([]),
     "e-values must be non-negative, non-empty and not NaN"),
    ("evalues/ebh_negative",
     lambda: evalues.e_bh([1.0, -0.5]),
     "an e-value is non-negative; a negative value was passed"),
    ("evalues/ebh_alpha",
     lambda: evalues.e_bh([1.0], alpha=1.0),
     "alpha must satisfy 0 < alpha < 1"),
    ("missingness/total_count_type",
     lambda: missingness.finite_collection_bounds([1, 0], total_count=2.5),
     "total_count must be a non-negative integer"),
    ("missingness/total_count_small",
     lambda: missingness.finite_collection_bounds([1, 0, 1], total_count=2),
     "total_count cannot be smaller than the number of supplied values"),
    ("missingness/binary",
     lambda: missingness.finite_collection_bounds([0, 1, 2]),
     "values must be binary numeric 0/1 or NaN"),
    ("missingness/numeric",
     lambda: missingness.finite_collection_bounds([0, "resistant"]),
     "values must be binary numeric 0/1 or NaN"),
    ("missingness/bounds_records",
     lambda: missingness.difference_bounds({}, {"lower_bound": 0.0, "upper_bound": 1.0}),
     "difference_bounds requires two non-empty bounds records"),
    ("missingness/bounds_finite",
     lambda: missingness.difference_bounds(
         {"lower_bound": 0.0, "upper_bound": float("inf")},
         {"lower_bound": 0.0, "upper_bound": 1.0}),
     "difference_bounds requires finite bound endpoints"),
    ("qc/support_plan",
     lambda: qc.support_pairing_plan(10, 12),
     "require 0 <= singletons <= n and 0 < threshold < 1"),
    ("phenotype/duplicate_policy",
     lambda: resolve_duplicates(_panel(), ["Strain_ID", "antibiotic"],
                                ["resistant_phenotype"], policy="ignore"),
     "unknown duplicate policy 'ignore'"),
    ("phenotype/conflicts",
     lambda: resolve_duplicates(
         _long([("i0", "a", "Resistant"), ("i0", "a", "Susceptible")]),
         ["Strain_ID", "antibiotic"], ["resistant_phenotype"], what="phenotype table"),
     "phenotype table has conflicting records for 1 key(s) in column(s) "
     "'resistant_phenotype': [('i0', 'a')]; correct the source or explicitly "
     "select drop_conflicts or positive_wins"),
    ("phenotype/kind",
     lambda: to_non_susceptible(_panel(), kind="mic"),
     "unknown phenotype kind 'mic'"),
    ("phenotype/intermediate",
     lambda: to_non_susceptible(_panel(), kind="clinical_sir", intermediate="keep"),
     "unknown intermediate policy 'keep'"),
    ("phenotype/intermediate_scope",
     lambda: to_non_susceptible(_panel(), kind="binary", intermediate="drop"),
     "intermediate policy applies only to clinical_sir or undeclared"),
    ("phenotype/missing_key",
     lambda: to_non_susceptible(_long([("i0", "", "Resistant")]), kind="clinical_sir"),
     "missing key in antibiotic"),
    ("phenotype/unreadable_vocabulary",
     lambda: to_non_susceptible(
         _long([("i0", "a", "S"), ("i1", "a", "R")]), kind="undeclared"),
     "no call in the phenotype table could be read under phenotype_kind "
     "'undeclared': 2 unreadable call(s), 's' (1), 'r' (1); these values are read "
     "by phenotype_kind: clinical_sir"),
    ("comparison/dimensions",
     lambda: comparison.compare_lineage_definitions([0, 1], ["A", "B"], ["A"]),
     "outcome and both lineage vectors must be one-dimensional and equally long"),
    ("comparison/binary",
     lambda: comparison.compare_lineage_definitions([0, 2], ["A", "B"], ["A", "B"]),
     "outcomes must be binary 0/1 or NaN"),
    ("comparison/seed",
     lambda: comparison.compare_lineage_definitions([0, 1], ["A", "B"], ["A", "B"], seed=-1),
     "seed must be a nonnegative integer shared by all four arms"),
    ("comparison/names",
     lambda: comparison.compare_lineage_definitions([0, 1], ["A", "B"], ["A", "B"],
                                                    name_a="X", name_b="X"),
     "lineage names must be nonempty and distinct"),
    ("comparison/id_count",
     lambda: comparison.compare_lineage_definitions([0, 1], ["A", "B"], ["A", "B"],
                                                    ids=["one"]),
     "ids must have one entry per outcome"),
    ("comparison/id_missing",
     lambda: comparison.compare_lineage_definitions([0, 1], ["A", "B"], ["A", "B"],
                                                    ids=["one", " "]),
     "ids must be nonmissing and unique"),
    ("comparison/id_unique",
     lambda: comparison.compare_lineage_definitions([0, 1], ["A", "B"], ["A", "B"],
                                                    ids=["one", "one"]),
     "ids must be unique after string conversion"),
]


@pytest.mark.parametrize("name, call, message", CASES, ids=[c[0] for c in CASES])
def test_the_refusal_prints_the_whole_sentence(name, call, message):
    with pytest.raises(ValueError) as caught:
        call()
    assert str(caught.value) == message


MEASURED = ("stats", "clonality", "attribution", "realised", "evalues",
            "missingness", "qc", "phenotype", "comparison")


def _refusal_lines():
    import ast
    from pathlib import Path

    root = Path(__file__).resolve().parent.parent / "src" / "amr_clonalshare"
    wanted = {}
    for module in MEASURED:
        path = root / f"{module}.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        # Under `mutmut run` the sources on disk carry one copy of every
        # function per mutant. This check is about the repository's source,
        # not about a mutated copy of it, so it steps aside there rather
        # than reporting refusals that cannot be reached by construction.
        if any(isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
               and n.name.endswith("__mutmut_orig") for n in ast.walk(tree)):
            pytest.skip("the sources on disk are a mutation run's copies")
        for node in ast.walk(tree):
            if (isinstance(node, ast.Raise) and isinstance(node.exc, ast.Call)
                    and getattr(node.exc.func, "id", "") == "ValueError"):
                wanted[(str(path), node.lineno)] = f"{module}.py:{node.lineno}"
    return wanted


def test_the_compare_command_names_the_columns_the_table_lacks(tmp_path, capsys):
    table = tmp_path / "input.csv"
    table.write_text("isolate_id,positive\ni0,1\n", encoding="utf-8")
    status = comparison.main(["--input", str(table), "--id-column", "isolate_id",
                              "--outcome", "positive", "--lineage-a", "v1",
                              "--lineage-b", "v2", "--output", str(tmp_path / "out")])
    assert status == 2
    assert capsys.readouterr().err.strip() == (
        "Comparison not completed: CSV columns missing: ['v1', 'v2']")


def test_the_compare_command_names_the_settings_it_refuses(tmp_path, capsys):
    table = tmp_path / "input.csv"
    table.write_text("isolate_id,positive\ni0,1\n", encoding="utf-8")
    status = comparison.main(["--input", str(table), "--id-column", "isolate_id",
                              "--outcome", "positive", "--lineage-a", "v1",
                              "--lineage-b", "v2", "--output", str(tmp_path / "out"),
                              "--folds", "1"])
    assert status == 2
    assert capsys.readouterr().err.strip() == (
        "Comparison not completed: --permutations and --repeats must be at least 1, "
        "--folds at least 2 and --bootstraps not negative")


def test_every_refusal_of_the_measured_modules_is_exercised(tmp_path):
    """A refusal no test reaches cannot have its wording checked."""
    import sys

    wanted = _refusal_lines()
    files = {path for path, _ in wanted}
    seen = set()

    def tracer(frame, event, arg):
        if event == "call":
            return tracer if frame.f_code.co_filename in files else None
        if event == "line":
            seen.add((frame.f_code.co_filename, frame.f_lineno))
        return tracer

    table = tmp_path / "input.csv"
    table.write_text("isolate_id,positive\ni0,1\n", encoding="utf-8")
    # Coverage measures with a tracer of its own. Replacing it and then
    # setting None would switch coverage off for everything that runs after
    # this test, which silently reported whole modules as unexercised; the
    # previous tracer is put back instead.
    previous = sys.gettrace()
    sys.settrace(tracer)
    try:
        for _name, call, _message in CASES:
            with pytest.raises(ValueError):
                call()
        comparison.main(["--input", str(table), "--id-column", "isolate_id",
                         "--outcome", "positive", "--lineage-a", "v1",
                         "--lineage-b", "v2", "--output", str(tmp_path / "out")])
        comparison.main(["--input", str(table), "--id-column", "isolate_id",
                         "--outcome", "positive", "--lineage-a", "v1",
                         "--lineage-b", "v2", "--output", str(tmp_path / "out"),
                         "--permutations", "0"])
    finally:
        sys.settrace(previous)
    missing = sorted(where for key, where in wanted.items() if key not in seen)
    assert not missing, (
        "no test reaches these refusals, so their wording is unchecked: "
        + ", ".join(missing))


def test_the_draft_refuses_a_missing_table_by_name(tmp_path):
    with pytest.raises(ValueError) as caught:
        draft_config([tmp_path / "absent.csv"])
    assert str(caught.value) == f"{tmp_path / 'absent.csv'} does not exist"


def test_the_draft_refuses_an_empty_call_list():
    with pytest.raises(ValueError) as caught:
        draft_config([])
    assert str(caught.value) == "name at least one CSV table"


def test_the_draft_refuses_a_table_without_a_header(tmp_path):
    empty = tmp_path / "empty.csv"
    empty.write_text("", encoding="utf-8")
    with pytest.raises(ValueError) as caught:
        draft_config([empty])
    assert str(caught.value) == f"{empty} has no header row"


def test_the_draft_names_the_missing_workbook_reader(tmp_path, monkeypatch):
    import pandas as pd

    book = tmp_path / "book.xlsx"
    book.write_bytes(b"not a workbook")

    def no_openpyxl(*args, **kwargs):
        raise ImportError("openpyxl")

    monkeypatch.setattr(pd, "read_excel", no_openpyxl)
    with pytest.raises(ValueError) as caught:
        draft_config([book])
    assert str(caught.value) == (
        f"{book} is a workbook; reading one needs openpyxl, which the core "
        "install leaves out: pip install 'amr-clonalshare[excel]'")
