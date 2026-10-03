"""Every agent gets one conclusion, derived from the record alone."""
import csv
import json
from pathlib import Path

import pytest

from amr_clonalshare.verdict import agent_verdict, format_interval, format_p, short_label, summary_rows

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"


def _block(**kw):
    base = dict(kappa_adj=0.2, observed_low=0.1, observed_high=0.3, p_value=0.001, p_floor=0.001,
                estimable=True, latent_order_lower=0.0, latent_order_lower_limit=0.0)
    base.update(kw)
    return base


def test_the_four_conclusions():
    assert agent_verdict(_block())["label"] == "lineage structure established"
    assert agent_verdict(_block(observed_low=0.0, p_value=0.2, latent_order_lower_limit=0.02))["label"] == "lineage structure established"
    # a positive lower bound alone, without a positive limit, establishes nothing
    assert agent_verdict(_block(observed_low=0.0, p_value=0.2, latent_order_lower=0.01))["label"] == "not established"
    assert agent_verdict(_block(observed_low=0.0, p_value=0.2))["label"] == "not established"
    # with a selection across agents the q-value decides, not the p-value
    assert agent_verdict(_block(), q=0.2)["label"] == "not established"
    selected = agent_verdict(_block(p_value=0.04), q=0.03)
    assert selected["label"] == "lineage structure established" and "q 0.030" in selected["text"]
    spans = agent_verdict(_block(observed_low=0.0, observed_high=1.0, p_value=0.3))
    assert spans["label"] == "interval spans the whole range" and "spans the whole range" in spans["text"]
    refused = agent_verdict(dict(estimable=False, reason="fewer than two lineages hold repeated readings"))
    assert refused == {"label": "not estimable", "text": "Not estimable: fewer than two lineages hold repeated readings."}
    assert agent_verdict(None)["label"] == "not estimable"


def test_the_smallest_p_value_of_the_permutations_is_written_as_a_bound():
    assert format_p(0.001, 0.001) == "≤ 0.001"
    assert format_p(0.002, 0.001) == "0.002"
    assert format_p(0.05, 0.05) == "≤ 0.050"
    assert format_p(0.3, None) == "0.300"
    assert format_p(None, 0.001) == "—"
    assert format_interval(0.0, 1.0) == "spans the whole range"
    assert format_interval(0.0, 0.314) == "0.000 to 0.314"


@pytest.mark.skipif(not (EXAMPLES / "ecoli_swine" / "expected").is_dir(), reason="examples folder not shipped")
def test_the_shipped_ecoli_record_is_established_for_four_agents():
    record = json.loads((EXAMPLES / "ecoli_swine" / "expected" / "clonal_share_result.json").read_text(encoding="utf-8"))
    rows = {r["agent"]: r for r in summary_rows(record)}
    for agent in ("ciprofloxacin", "cefotaxime", "ceftazidime", "gentamicin"):
        assert rows[agent]["mic_conclusion"]["label"] == "lineage structure established"
    assert rows["ceftazidime"]["mic_p"] == "≤ 0.001"
    assert rows["meropenem"]["mic_interval"] == "spans the whole range"
    assert rows["meropenem"]["mic_conclusion"]["label"] == "interval spans the whole range"
    assert rows["tetracycline"]["mic_conclusion"]["label"] == "not established"
    assert "transmission" not in rows["tetracycline"]["mic_conclusion"]["text"]


def test_results_csv_and_the_console_carry_the_conclusion(tmp_path, capsys):
    from amr_clonalshare.cli import main
    source = EXAMPLES / "workflows"
    if not source.is_dir():
        pytest.skip("examples folder not shipped")
    out = tmp_path / "out"
    assert main(["--config", str(source / "mic.yaml"), "--results-dir", str(out), "--seed", "3"]) == 0
    printed = capsys.readouterr().out
    assert printed.lstrip().startswith("Antimicrobial")
    assert "Conclusion" in printed and "Report: " in printed and "{" not in printed
    assert "demo_agent:" in printed
    with (out / "results.csv").open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    conclusions = {(r["agent"], r["analysis"]): r["conclusion"] for r in rows}
    assert conclusions[("demo_agent", "mic_order")].endswith(").")
    assert all(not r["conclusion"] for r in rows if r["analysis"] not in ("mic_order", "collection_membership"))
    order = next(r for r in rows if r["analysis"] == "mic_order")
    assert order["p_value"] and order["selected"] in ("yes", "no")
    assert main(["--config", str(source / "mic.yaml"), "--results-dir", str(tmp_path / "json"), "--seed", "3", "--json"]) == 0
    printed = capsys.readouterr().out
    assert printed.lstrip().startswith("{")
    report = (out / "report.md").read_text(encoding="utf-8")
    assert "| Conclusion |" in report and "≤ " in report


@pytest.mark.skipif(not (EXAMPLES / "ssuis" / "expected").is_dir(), reason="examples folder not shipped")
def test_the_report_reads_every_share_as_the_form_does_with_its_q_value():
    from amr_clonalshare.report import render_report
    record = json.loads((EXAMPLES / "ssuis" / "expected" / "clonal_share_result.json").read_text(encoding="utf-8"))
    text = render_report(record, {})
    md = record["metadata_diagnostics"]
    rows = {r["agent"]: r for r in summary_rows(record)}

    def table(label):
        block = text.split("**%s.**" % label, 1)[1].split("\n\n")[1]
        lines = [[c.strip() for c in line.strip("|").split("|")] for line in block.splitlines()]
        return lines[0], {cells[0]: cells for cells in lines[2:]}

    head, calls = table("Table 1")
    assert head == ["Trait", "Share", "95 % interval", "p", "q", "Control", "e-value", "Selected", "Conclusion"]
    selection = md["lineage_selection"]
    assert set(calls) == set(md["clonal_share"])
    for agent, cells in calls.items():
        assert cells[3] == rows[agent]["call_p"]
        assert cells[4] == "%.3f" % selection["q_values"][agent]
        assert cells[7] == ("yes" if agent in selection["rejected_features"] else "no")
        assert cells[8] == short_label(rows[agent]["call_conclusion"]["label"])
    head, orders = table("Table 5")
    assert head[4:6] == ["p", "q"]
    order_q = md["censored_share"]["order_selection"]["q_values"]
    for agent, cells in orders.items():
        if cells[2] != "not estimable":
            assert cells[5] == "%.3f" % order_q[agent]
            assert cells[-1] == short_label(rows[agent]["mic_conclusion"]["label"])
    head, structure = table("Table 5b")
    assert head == ["Agent", "Lineages scored", "Singletons set aside", "Effective", "Support"]
    assert set(structure) == set(md["censored_share"]["per_agent"])
