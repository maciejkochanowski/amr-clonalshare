"""Small and partially observed input must yield diagnostics, not exceptions."""
import json

import numpy as np
import pandas as pd
import pytest

from amr_clonalshare.attribution import _debias, clonal_share
from amr_clonalshare.config import from_dict
from amr_clonalshare.evalues import SequentialEResult, sequential_e_process
from amr_clonalshare.io import load_dataset
from amr_clonalshare.jsonio import dumps
from amr_clonalshare.qc import input_qc, render_markdown


@pytest.mark.parametrize("null_mean", [1., 1.01])
def test_debias_refuses_nonpositive_correction_denominator(null_mean):
    assert np.isnan(_debias(0.5, null_mean))


def test_saturated_permutation_null_returns_a_nonestimable_diagnostic():
    # With this fixed seed the only permuted grouping also predicts perfectly.
    result = clonal_share([1, 0, 1, 0], ["a", "b", "a", "b"],
                          n_boot=4, n_perm=1, repeats=1, seed=0)
    assert result.null_mean == 1.
    assert result.kappa == 1.
    assert not result.estimable
    assert np.isnan(result.kappa_adj)
    assert np.isnan(result.ci_low) and np.isnan(result.ci_high)
    assert "null" in result.reason and "denominator" in result.reason


@pytest.mark.parametrize("labels", [["a"], [None, None]])
def test_qc_render_accepts_undefined_effective_size_after_strict_json(labels):
    panel = pd.DataFrame({"Drug": [0.] * len(labels)})
    qc = input_qc(panel, metadata=pd.DataFrame({"lineage": labels}),
                  lineage_column="lineage")
    restored = json.loads(dumps(qc))
    assert restored["lineage"]["effective_group_size"] is None
    text = render_markdown(restored)
    assert "Effective lineage size for the between-lineage variance: not defined" in text
    assert "no between-lineage comparison" in text


def test_missing_call_inventory_uses_raw_call_and_mic_union(tmp_path):
    pd.DataFrame({"Strain_ID": list("abcdefghi"), "ST": ["a"] * 5 + ["b"] * 4})\
        .to_csv(tmp_path / "meta.csv", index=False)
    pd.DataFrame({"Strain_ID": list("ab"), "antibiotic": ["Drug"] * 2,
                  "resistant_phenotype": ["R", "S"]}).to_csv(tmp_path / "calls.csv", index=False)
    pd.DataFrame({"Strain_ID": list("bcdefghi"), "antibiotic": ["Drug"] * 8,
                  "measurement": [1, 2, 4, 8, 1, 2, 4, 8]}).to_csv(tmp_path / "mic.csv", index=False)
    cfg = from_dict({"dataset": {"name": "union", "data_dir": str(tmp_path),
                    "phenotype": "calls.csv", "phenotype_kind": "clinical_sir",
                    "metadata": "meta.csv", "lineage_column": "ST",
                    "mic": "mic.csv"}})
    ds = load_dataset(cfg)
    assert len(ds.strain_ids) == 9
    reading = ds.input_qc["phenotype_reading"]
    assert reading["isolates_without_readable_calls"] == 7
    assert reading["call_table_isolates_without_readable_calls"] == 0
    assert reading["agents"]["Drug"]["n_absent_records"] == 7
    assert "Isolates without any readable call: 7" in render_markdown(ds.input_qc)


@pytest.mark.parametrize("batches,labels,analyzed,scored", [
    ([], [], [], 0),
    ([[np.nan, np.nan]], [["a", "b"]], [0], 0),
    ([[0, 1]], [["a", "b"]], [2], 0),
    ([[0, np.nan], [1, 0]], [["a", "b"], ["a", "b"]], [1, 2], 0),
    ([[0, 1], [1, 0]], [["a", "a"], ["a", "b"]], [2, 2], 0),
    ([[0, 1], [1, 0]], [["a", "b"], ["a", "b"]], [2, 2], 1),
])
def test_sequential_counts_distinguish_analyzed_and_actually_scored_batches(
        batches, labels, analyzed, scored):
    result = sequential_e_process(batches, labels)
    assert result.n_analyzed_per_batch == analyzed
    assert result.n_scored_batches == scored


def test_sequential_new_counts_preserve_archived_sizes_and_overlapping_drops():
    result = sequential_e_process([[0, 1, np.nan], [1, 0, np.nan], [np.nan]],
                                  [["a", "b", None], ["a", "b", "c"], [None]])
    assert result.n_per_batch == [2, 3, 0]
    assert result.n_analyzed_per_batch == [2, 2, 0]
    assert result.n_scored_batches == 1
    assert result.n_dropped_non_finite == 3
    assert result.n_dropped_untyped == 2


def test_sequential_result_positional_construction_has_empty_batch_defaults():
    result = SequentialEResult([0.], [1.], False, False, 1, [2])
    assert result.n_analyzed_per_batch == []
    assert result.n_scored_batches == 0
