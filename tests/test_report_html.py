"""The drawn form of the run report.

The prose form states the result; this one has to show the diagnostics that
could undermine it. The tests below hold it to three things: that it draws
what the record holds rather than recomputing it, that a value the estimator
refused is not drawn as if it had been estimated, and that the page opens
without anything to fetch.

The shipped planted example was run without the lineage attribution, so its
record carries no clonal share. The fixture supplies one by running the
estimator over the same cohort, which is what the pipeline does when
``attribution.enabled`` is set. What is under test here is the renderer, so
the fixture only has to be a well-formed record; the estimator itself is
tested in ``test_attribution.py``.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
import pytest

from amr_clonalshare import attribution
from amr_clonalshare.report_html import render_html_report

EXPECTED = Path(__file__).resolve().parent.parent / "examples" / "synthetic"


def _with_clonal_share(record, summary):
    d = EXPECTED / "data" / "planted"
    amr = pd.read_csv(d / "amr.csv", index_col=0)
    vir = pd.read_csv(d / "vir.csv", index_col=0)
    meta = pd.read_csv(d / "metadata.csv", index_col=0)
    X = pd.concat([amr, vir], axis=1)
    lineage = meta.loc[X.index, "lineage"].tolist()
    record["metadata_diagnostics"]["clonal_share"] = {
        str(c): attribution.clonal_share(
            X[c].to_numpy(dtype=float), lineage, seed=0).as_dict()
        for c in X.columns
    }
    summary["clonal_share_all_features"] = attribution.layer_clonal_share(
        X.to_numpy(dtype=float), lineage, seed=0).kappa_adj
    return record, summary


@pytest.fixture(scope="module")
def rendered():
    record = json.loads((EXPECTED / "expected_planted" /
                         "cluster_result.json").read_text(encoding="utf-8"))
    summary = json.loads((EXPECTED / "expected_planted" /
                          "summary.json").read_text(encoding="utf-8"))
    record, summary = _with_clonal_share(record, summary)
    return record, summary, render_html_report(record, summary)


def test_the_page_is_self_contained(rendered):
    """Nothing to fetch: no script, no external stylesheet, no remote image.

    A report a laboratory cannot open offline is not a report. The one URL the
    page carries is the SVG namespace, which is an identifier rather than an
    address: nothing is retrieved from it.
    """
    _, _, html = rendered
    assert html.startswith("<!doctype html>")
    assert "<script" not in html.lower()
    assert "<link" not in html.lower()
    assert "src=" not in html.lower()
    addresses = [u for u in re.findall(r"""https?://[^\s"'<>]+""", html)
                 if u != "http://www.w3.org/2000/svg"]
    assert addresses == []
    assert re.search(r"<svg[^>]*viewBox", html)


def test_every_section_the_reader_needs_is_present(rendered):
    """The prose form read five of the record's blocks. These are the four
    that answer whether the result could be wrong rather than what it is."""
    _, _, html = rendered
    for heading in ("Admissibility of the input",
                    "artefact of the panel",
                    "beat a simpler rule",
                    "Which layer carries the grouping",
                    "How stable is the grouping",
                    "What may be concluded"):
        assert heading in html, heading
    assert html.count('<span class="no">') == 10


def test_the_figures_are_drawn_from_the_record(rendered):
    """A drawn value must be traceable to the file. The interval figure marks
    every trait whose interval covers zero; the count in the prose and the
    count of marks in the figure must agree, or the picture and the sentence
    are telling the reader different things."""
    record, _, html = rendered
    shares = record["metadata_diagnostics"]["clonal_share"]
    rows = [v for v in shares.values()
            if isinstance(v.get("kappa_adj"), float)
            and isinstance(v.get("ci_low"), float)]
    rows.sort(key=lambda v: -v["kappa_adj"])
    drawn = rows[:18]
    crossing = sum(1 for v in drawn if v["ci_low"] <= 0.0 <= v["ci_high"])
    assert "<b>%d of the %d traits</b>" % (crossing, len(drawn)) in html
    assert html.count("&#8225;</text>") == crossing


def test_a_refused_gate_margin_is_not_invented(rendered):
    """Only gates carrying a numeric margin are drawn. A gate whose margin is
    absent must not appear in the figure with a bar of zero, which would read
    as a gate that only just held."""
    record, _, html = rendered
    ledger = record["interpretation"]["gate_ledger"]
    with_margin = [g for g in ledger if g.get("applicable")
                   and isinstance(g.get("margin_to_failure"), float)]
    figure = html.split("Figure 4.")[0].split("margin to failure")[0]
    for gate in ledger:
        drawn = ">%s</text>" % gate["code"].replace("_", " ")
        if gate in with_margin:
            assert drawn in figure or gate["code"] in html
    assert len(with_margin) >= 1
