"""The run report is rendered from the record, in a fixed order, with nothing invented.

A planted cohort is the real case; synthetic records exercise the branches a
run does not reach (a refused agent, no e-values, no lineage column). Both
forms of the report are built from one structure, so what is asserted of the
prose form holds of the page, and the page has to open with nothing to fetch.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pytest

from amr_clonalshare import run
from amr_clonalshare.cli import _summary
from amr_clonalshare.report import render_report
from amr_clonalshare.report_html import render_html_report

ROOT = Path(__file__).resolve().parents[1]
SHIPPED = ROOT / "examples" / "ssuis" / "expected" / "clonal_share_result.json"

FORBIDDEN = (
    r"\bNone\b", r"\bnan\b", r"\bNaN\b", r"\[\]", r"\b0 of the 0\b",
    r"\b0 of 0\b", r"result\[", r"run\(\.\.\.", r"n_perm=", r"\bO\(n",
    r"metadata_diagnostics\.", r"not computed to not computed",
)


def _sections(text):
    return [line[3:] for line in text.splitlines() if line.startswith("## ")]


def _clean(html):
    return re.sub(r"<[^>]+>", " ", html)


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    from conftest import planted_cohort, write_config
    tmp = tmp_path_factory.mktemp("planted")
    _, cfg = write_config(tmp, planted_cohort(tmp))
    record = run(cfg, results_dir=tmp / "out", seed=7)
    raw = (tmp / "out" / "clonal_share_result.json").read_bytes()
    return record, raw


def test_every_number_in_the_report_is_read_from_the_record(planted):
    record, raw = planted
    summary = _summary(record)
    text = render_report(record, summary, record_bytes=raw)
    heads = _sections(text)
    assert heads[0].endswith("Measurement summary")
    assert heads[1].endswith("Admissibility of the input")
    assert heads[2].endswith("trait by trait")
    assert heads[-2].endswith("Provenance and terms")
    assert f"**Isolates:** {record['n_isolates']}" in text
    assert "**Lineage:** lineage" in text
    per = record["metadata_diagnostics"]["clonal_share"]
    for name, d in per.items():
        assert f"| {name} |" in text
        assert f"{d['kappa_adj']:.3f}" in text
    assert hashlib.sha256(raw).hexdigest() in text
    for pat in FORBIDDEN:
        assert not re.search(pat, text), pat


def test_the_page_is_self_contained_and_carries_the_same_numbers(planted):
    """Nothing to fetch: no script, no external stylesheet, no remote image.
    The one URL the page carries is the SVG namespace, an identifier rather
    than an address. And a number in one form is in the other."""
    record, raw = planted
    summary = _summary(record)
    html = render_html_report(record, summary, record_bytes=raw)
    assert html.startswith("<!doctype html>")
    assert "<script" not in html.lower() and "<link" not in html.lower()
    assert "src=" not in html.lower()
    addresses = [u for u in re.findall(r"""https?://[^\s"'<>]+""", html)
                 if u != "http://www.w3.org/2000/svg"]
    assert addresses == []
    assert re.search(r"<svg[^>]*viewBox", html)
    md = render_report(record, summary, record_bytes=raw)
    for d in record["metadata_diagnostics"]["clonal_share"].values():
        token = "%.3f" % d["kappa_adj"]
        assert token in md and token in html
    order = ("Measurement summary", "Admissibility of the input",
             "trait by trait", "survives re-reading", "recorded resolution",
             "carried across lineages", "Provenance and terms")
    positions = [html.index(h) for h in order]
    assert positions == sorted(positions)
    text = _clean(html)
    for pat in FORBIDDEN:
        assert not re.search(pat, text), pat


def test_the_interval_figure_marks_what_the_prose_counts(planted):
    """The figure marks every trait whose interval covers zero; the count in
    the prose and the count of marks must agree."""
    record, raw = planted
    html = render_html_report(record, _summary(record), record_bytes=raw)
    shares = record["metadata_diagnostics"]["clonal_share"]
    rows = sorted(shares.values(), key=lambda v: -v["kappa_adj"])[:18]
    crossing = sum(1 for v in rows if v["ci_low"] <= 0.0 <= v["ci_high"])
    assert "<b>%d of the %d traits</b>" % (crossing, len(rows)) in html
    assert html.count("&#8225;</text>") == crossing


def test_a_refused_agent_is_read_as_refused():
    record = {"n_isolates": 30, "n_traits": 2, "seed": 1,
              "metadata_diagnostics": {
                  "lineage_column": "st",
                  "clonal_share": {
                      "ok": {"kappa_adj": 0.4, "ci_low": 0.2, "ci_high": 0.6,
                             "null_mean": 0.0, "support": 0.95, "estimable": True,
                             "n_groups": 6},
                      "thin": {"kappa_adj": 0.1, "ci_low": -0.1, "ci_high": 0.3,
                               "null_mean": 0.0, "support": 0.5,
                               "estimable": False, "n_groups": 6}}}}
    text = render_report(record, _summary(record))
    assert "1 of 2 antimicrobials estimable" in text
    assert "| thin |" in text and "refused" in text
    assert "†" in text
    assert "No e-values were computed" in text
    assert "No recorded dilutions were supplied" in text
    assert "needs two collections" in text
    for pat in FORBIDDEN:
        assert not re.search(pat, text), pat


def test_the_empty_record_renders_without_asserting_anything():
    record = {"n_isolates": 0, "seed": 0}
    summary = _summary(record)
    for body in (render_report(record, summary),
                 _clean(render_html_report(record, summary))):
        assert "No antimicrobial carried a readable call" in body
        assert "No input check was attached" in body
        for pat in FORBIDDEN:
            assert not re.search(pat, body), pat


def test_formatting_never_invents_a_number():
    from amr_clonalshare.report import _f, _pct
    from amr_clonalshare.report_model import evalue, interval
    assert _f(None) == "not computed" and _f(float("nan")) == "not computed"
    assert _f("withheld") == "withheld" and _f(0.12345) == "0.123"
    assert _pct(None) == "not computed" and _pct(0.5) == "50.0 %"
    assert interval(None, 0.2) == "not computed"
    assert evalue(2.7) == "2.7" and evalue(6.3e18).startswith("6.3 × 10")


@pytest.mark.skipif(not SHIPPED.is_file(), reason="shipped S. suis record not present")
def test_the_shipped_record_renders_cleanly():
    record = json.loads(SHIPPED.read_text(encoding="utf-8"))
    summary = _summary(record)
    md = render_report(record, summary, record_bytes=SHIPPED.read_bytes())
    html = _clean(render_html_report(record, summary,
                                     record_bytes=SHIPPED.read_bytes()))
    for body in (md, html):
        for pat in FORBIDDEN:
            assert not re.search(pat, body), pat
    assert "677" in md and "baps_cluster" in md


def test_the_page_draws_one_figure_per_question_the_record_can_answer():
    raw = SHIPPED.read_bytes()
    record = json.loads(raw)
    summary = _summary(record)
    html = render_html_report(record, summary, input_qc=record.get("input_qc"),
                              record_bytes=raw)
    md = render_report(record, summary, input_qc=record.get("input_qc"),
                       record_bytes=raw)
    labels = re.findall(r"<b>(Figure \d+)\.</b>", html)
    assert labels == ["Figure %d" % i for i in range(1, 7)]
    assert html.count("<figure>") == 6
    # the prose form carries every caption in the same order, drawn or not
    assert [m for m in re.findall(r"\*(Figure \d+) \(drawn on the page\)", md)] == labels
    # Table 6 states what the carriage figure draws
    assert "**Table 6.**" in md and "as many lineages as chance" in md
    text = _clean(html)
    for pat in FORBIDDEN:
        assert not re.search(pat, text), pat
    # a record without dilutions and intakes draws fewer figures, never a blank one
    slim = json.loads(raw)
    slim["metadata_diagnostics"].pop("censored_share")
    slim["metadata_diagnostics"]["lineage_evidence"].pop("sequential")
    html2 = render_html_report(slim, _summary(slim), input_qc=slim.get("input_qc"))
    assert html2.count("<figure>") == 4
    assert "<figure></figure>" not in html2
