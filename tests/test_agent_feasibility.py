"""Actionable input diagnostics use each trait's actual analysed subset."""
import numpy as np
import pandas as pd
import pytest

from amr_clonalshare.cli import _summary
from amr_clonalshare.qc import input_qc, render_markdown
from amr_clonalshare.report import render_report
from amr_clonalshare.report_html import render_html_report


def test_agent_missingness_changes_support_and_retained_outcome_counts():
    panel = pd.DataFrame({"complete": [0, 1] * 5,
                          "sparse": [0, np.nan, 1, np.nan, 0, np.nan, 1, np.nan, 0, np.nan]})
    meta = pd.DataFrame({"lineage": list("aabbccddee")})
    qc = input_qc(panel, metadata=meta, lineage_column="lineage")
    a, b = (qc["agent_feasibility"][name] for name in ("complete", "sparse"))
    assert a["support"] == 1 and a["input_feasible"]
    assert b["support"] == 0 and not b["input_feasible"]
    assert b["n_retained"] == 5 and b["n_dropped_non_finite"] == 5
    assert b["minor_count"] == 2
    assert b["support_pairings_needed"] == 4
    assert b["failure_reasons"] == ["support below 0.80"]
    assert "fewer than 20 isolates of the rarer outcome" in b["warnings"]


def test_typed_tested_intersection_matches_estimator_filter_order():
    panel = pd.DataFrame({"drug": [0, np.nan, 1, 1, 0, np.nan]})
    meta = pd.DataFrame({"lineage": ["a", None, None, "b", "b", "c"]})
    rec = input_qc(panel, metadata=meta, lineage_column="lineage")["agent_feasibility"]["drug"]
    # the second isolate is unread and untyped and is counted in both
    assert (rec["n_dropped_non_finite"], rec["n_dropped_untyped"], rec["n_retained"]) == (2, 2, 3)
    assert rec["support"] == pytest.approx(2/3)
    assert rec["prevalence"] == pytest.approx(1/3)
    from amr_clonalshare.attribution import clonal_share
    fitted = clonal_share(panel["drug"], meta["lineage"], folds=2, repeats=1,
                          n_boot=0, n_perm=0, seed=7)
    assert rec["n_retained"] == fitted.n
    assert rec["n_dropped_non_finite"] == fitted.n_dropped_non_finite
    assert rec["n_dropped_untyped"] == fitted.n_dropped_untyped
    assert rec["support"] == fitted.support
    assert rec["input_feasible"] == fitted.estimable


@pytest.mark.parametrize("n,s", [(10, 1), (11, 2), (20, 3), (100, 10), (100, 11), (50, 50)])
def test_support_plan_is_minimal_and_exact_at_boundary(n, s):
    from amr_clonalshare.qc import support_pairing_plan
    plan = support_pairing_plan(n, s)
    k = plan["additional_isolates"]
    from amr_clonalshare.attribution import SUPPORT_THRESHOLD as t
    assert (n - s + 2*k) / (n + k) >= t
    if k:
        assert (n - s + 2*(k-1)) / (n + k-1) < t
    assert 0 <= k <= s


def test_support_plan_cannot_fix_constant_trait_or_single_lineage():
    panel = pd.DataFrame({"drug": [0] * 6, "missing": [np.nan] * 6})
    meta = pd.DataFrame({"lineage": ["one"] * 6})
    recs = input_qc(panel, metadata=meta, lineage_column="lineage")["agent_feasibility"]
    assert recs["drug"]["support_pairings_needed"] == 0
    assert not recs["drug"]["input_feasible"]
    assert set(recs["drug"]["failure_reasons"]) == {"fewer than 2 lineages", "constant trait"}
    assert recs["missing"]["support_pairings_needed"] is None
    assert recs["missing"]["support"] is None


def test_feasibility_is_visible_in_both_run_reports_and_input_check():
    panel = pd.DataFrame({"drug": [0, 1] * 3})
    meta = pd.DataFrame({"lineage": list("abcdef")})
    qc = input_qc(panel, metadata=meta, lineage_column="lineage")
    record = {"n_isolates": 6, "seed": 1, "metadata_diagnostics": {}}
    outputs = [render_markdown(qc), render_report(record, _summary(record), input_qc=qc),
               render_html_report(record, _summary(record), input_qc=qc)]
    for output in outputs:
        assert "support below 0.80" in output
        assert "distinct singleton lineages" in output
        assert "does not guarantee" in output
        assert "drug" in output
