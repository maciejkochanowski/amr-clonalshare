"""The local form drives the same run as the command line.

The server is started on a free loopback port with a working folder of its
own, and the pages are driven with http.client: the tables are uploaded, the
columns are named, a run is started and its report is read back; the record
the form produced is the record the command line produces from the
configuration the form wrote. The form is also checked against the
configuration: every setting of the dataset section is either a field of the
form or named as one the form leaves to the file.
"""
from __future__ import annotations

import dataclasses
import http.client
import json
import re
import threading
import time
import urllib.parse
import uuid

import pandas as pd
import pytest

from amr_clonalshare import gui
from amr_clonalshare.config import DatasetConfig, load_config


@pytest.fixture(scope="module")
def server(tmp_path_factory):
    root = tmp_path_factory.mktemp("gui-root")
    ready = threading.Event()
    box = {}

    def on_ready(address):
        box["address"] = address
        ready.set()

    thread = threading.Thread(target=gui.serve, kwargs=dict(root=root, ready=on_ready), daemon=True)
    thread.start()
    assert ready.wait(10), "the server did not start"
    url = urllib.parse.urlsplit(box["address"])
    token = urllib.parse.parse_qs(url.query)["t"][0]
    return {"host": url.hostname, "port": url.port, "token": token, "root": root}


def _request(server, method, path, body=None, headers=None, query=None):
    query = dict(query or {})
    query.setdefault("t", server["token"])
    target = path + "?" + urllib.parse.urlencode(query)
    conn = http.client.HTTPConnection(server["host"], server["port"], timeout=60)
    conn.request(method, target, body=body, headers=headers or {})
    response = conn.getresponse()
    data = response.read()
    conn.close()
    return response.status, dict(response.getheaders()), data


def _multipart(fields: dict, files: dict):
    boundary = "----amr" + uuid.uuid4().hex
    out = []
    for name, value in fields.items():
        out.append(f"--{boundary}\r\nContent-Disposition: form-data; name=\"{name}\"\r\n\r\n{value}\r\n".encode())
    for name, (filename, payload) in files.items():
        out.append(f"--{boundary}\r\nContent-Disposition: form-data; name=\"{name}\"; filename=\"{filename}\"\r\n"
                   f"Content-Type: application/octet-stream\r\n\r\n".encode() + payload + b"\r\n")
    out.append(f"--{boundary}--\r\n".encode())
    return b"".join(out), {"Content-Type": f"multipart/form-data; boundary={boundary}"}


def _urlencoded(fields: dict):
    items = []
    for key, value in fields.items():
        for v in (value if isinstance(value, list) else [value]):
            items.append((key, v))
    return urllib.parse.urlencode(items).encode(), {"Content-Type": "application/x-www-form-urlencoded"}


def _wait(server, job_url, seconds=240):
    deadline = time.time() + seconds
    while time.time() < deadline:
        status, _, page = _request(server, "GET", job_url)
        assert status == 200
        text = page.decode("utf-8")
        if "Completed" in text or "Not completed" in text:
            return text
        time.sleep(1)
    raise AssertionError("the run did not end in time")


def test_the_pages_need_the_token_of_the_session(server):
    status, _, page = _request(server, "GET", "/", query={"t": "wrong"})
    assert status == 403 and b"Not this address" in page
    status, _, page = _request(server, "GET", "/")
    assert status == 200 and b"Metadata table" in page and b"<script" not in page


def test_every_dataset_setting_is_a_field_of_the_form_or_named_as_file_only():
    settings = {f.name for f in dataclasses.fields(DatasetConfig)}
    covered = {key for key, _ in gui.COLUMN_FIELDS} | set(gui.TEXT_FIELDS) | set(gui.FILE_ONLY_FIELDS)
    assert settings == covered


def test_a_run_from_the_form_is_the_run_of_the_configuration_it_wrote(server, tmp_path):
    from conftest import planted_cohort
    planted_cohort(tmp_path, n_per_lineage=6, n_lineages=6)
    data = tmp_path / "data"
    body, headers = _multipart({"name": "planted"}, {
        "metadata": ("meta.csv", (data / "meta.csv").read_bytes()),
        "calls": ("calls.csv", (data / "calls.csv").read_bytes())})
    status, _, page = _request(server, "POST", "/upload", body, headers)
    assert status == 200, page[:400]
    text = page.decode("utf-8")
    run = re.search(r'action="/run\?run=([^&"]+)', text).group(1)
    # the guesses: the lineage and the identifier were read off the headers
    assert re.search(r'<select id="lineage_column"[^>]*>.*?<option value="lineage" selected>', text, re.S)
    assert re.search(r'<select id="strain_id_column"[^>]*>.*?<option value="iid" selected>', text, re.S)
    fields = {"strain_id_column": "iid", "lineage_column": "lineage", "batch_column": "intake",
              "phenotype_id_column": "iid", "phenotype_antibiotic_column": "antibiotic",
              "phenotype_call_column": "call", "phenotype_kind": "undeclared", "duplicate_policy": "error",
              "n_perm": "40", "n_boot": "40", "folds": "3", "repeats": "2", "surveillance_n_boot": "50",
              "alpha": "0.05", "seed": "7", "action": "run"}
    body, headers = _urlencoded(fields)
    status, response_headers, _ = _request(server, "POST", "/run", body, headers, query={"run": run})
    assert status == 303, status
    job_url = urllib.parse.urlsplit(response_headers["Location"]).path
    page = _wait(server, job_url)
    assert "Completed" in page and "Lineage share of the call" in page, page[-2000:]
    folder = server["root"] / "runs" / run
    record = json.loads((folder / "results" / "clonal_share_result.json").read_text(encoding="utf-8"))
    assert record["n_isolates"] == 36
    # the report is served from the run folder, and the zip holds the run
    status, headers_out, report = _request(server, "GET", f"{job_url}/file/report.html")
    assert status == 200 and headers_out["Content-Type"].startswith("text/html") and b"<svg" in report
    status, headers_out, archive = _request(server, "GET", f"{job_url}/zip")
    assert status == 200 and headers_out["Content-Type"] == "application/zip" and archive[:2] == b"PK"
    # the configuration the form wrote is a configuration of the command line, and gives the same record
    from amr_clonalshare import core
    cfg = load_config(folder / "config.yaml")
    assert cfg.dataset.batch_column == "intake" and cfg.attribution.n_perm == 40 and cfg.evidence.alpha == 0.05
    again = core.run(cfg, results_dir=tmp_path / "again", seed=7)
    for agent, block in record["metadata_diagnostics"]["clonal_share"].items():
        assert again["metadata_diagnostics"]["clonal_share"][agent]["kappa_adj"] == block["kappa_adj"]
    # and the runs page lists it
    status, _, runs = _request(server, "GET", "/runs")
    assert status == 200 and b"planted" in runs


def test_a_refused_configuration_returns_to_the_form_with_the_reason(server, tmp_path):
    from conftest import planted_cohort
    planted_cohort(tmp_path, n_per_lineage=4, n_lineages=4)
    data = tmp_path / "data"
    body, headers = _multipart({"name": "refused"}, {
        "metadata": ("meta.csv", (data / "meta.csv").read_bytes()),
        "calls": ("calls.csv", (data / "calls.csv").read_bytes())})
    status, _, page = _request(server, "POST", "/upload", body, headers)
    run = re.search(r'action="/run\?run=([^&"]+)', page.decode()).group(1)
    body, headers = _urlencoded({"strain_id_column": "iid", "lineage_column": "", "phenotype_id_column": "iid",
                                 "phenotype_antibiotic_column": "antibiotic", "phenotype_call_column": "call"})
    status, _, page = _request(server, "POST", "/run", body, headers, query={"run": run})
    assert status == 400 and b"lineage_column" in page and b'<select id="lineage_column"' in page


def test_the_comparison_form_runs_the_matched_comparison(server, tmp_path):
    rng = __import__("numpy").random.default_rng(3)
    n = 120
    a = [f"A{i % 6}" for i in range(n)]
    b = [f"B{(i // 2) % 8}" for i in range(n)]
    outcome = [int(rng.random() < (0.8 if x in ("A0", "A1") else 0.2)) for x in a]
    table = pd.DataFrame({"isolate": [f"i{i}" for i in range(n)], "outcome": outcome, "st": a, "cluster": b})
    payload = table.to_csv(index=False).encode()
    body, headers = _multipart({"name": "two definitions"}, {"table": ("two.csv", payload)})
    status, _, page = _request(server, "POST", "/compare/upload", body, headers)
    assert status == 200, page[:300]
    text = page.decode()
    run = re.search(r'action="/compare/run\?run=([^&"]+)', text).group(1)
    assert '<option value="isolate" selected>' in text
    body, headers = _urlencoded({"id_column": "isolate", "outcome": "outcome", "lineage_a": "st",
                                 "lineage_b": "cluster", "permutations": "30", "bootstraps": "30",
                                 "folds": "3", "repeats": "2", "seed": "1"})
    status, response_headers, _ = _request(server, "POST", "/compare/run", body, headers, query={"run": run})
    assert status == 303
    job_url = urllib.parse.urlsplit(response_headers["Location"]).path
    page = _wait(server, job_url)
    assert "Completed" in page and "comparison.json" in page, page[-1500:]
    folder = server["root"] / "runs" / run
    result = json.loads((folder / "results" / "comparison.json").read_text(encoding="utf-8"))
    assert result["status"] == "computed" and result["names"] == {"a": "st", "b": "cluster"}
    # a second comparison of the same table keeps the first and gets a folder of its own
    status, response_headers, _ = _request(server, "POST", "/compare/run", body, headers, query={"run": run})
    assert status == 303
    page = _wait(server, urllib.parse.urlsplit(response_headers["Location"]).path)
    assert "Completed" in page, page[-1500:]
    again = json.loads((folder / "results_2" / "comparison.json").read_text(encoding="utf-8"))
    assert again["arms"] == result["arms"]
    assert (folder / "results" / "comparison.json").is_file()
    # a folder name with a separator is refused
    status, _, _ = _request(server, "POST", "/compare/run", body, headers, query={"run": "../" + run})
    assert status == 404


EXAMPLES = gui.examples_folder()


@pytest.fixture(scope="module")
def server_with_examples(tmp_path_factory):
    root = tmp_path_factory.mktemp("gui-examples")
    ready = threading.Event()
    box = {}

    def on_ready(address):
        box["address"] = address
        ready.set()

    thread = threading.Thread(target=gui.serve, kwargs=dict(root=root, ready=on_ready, examples=EXAMPLES), daemon=True)
    thread.start()
    assert ready.wait(10), "the server did not start"
    url = urllib.parse.urlsplit(box["address"])
    token = urllib.parse.parse_qs(url.query)["t"][0]
    return {"host": url.hostname, "port": url.port, "token": token, "root": root}


@pytest.mark.skipif(EXAMPLES is None, reason="the examples folder of the source tree is not beside the package")
def test_the_examples_page_runs_a_shipped_collection_as_the_command_line_does(server_with_examples, tmp_path):
    server = server_with_examples
    status, _, page = _request(server, "GET", "/examples")
    assert status == 200
    text = page.decode("utf-8")
    # every configuration file one folder below examples/ is listed, with its opening comment; so is the comparison table
    assert 'value="workflows/calls.yaml"' in text and 'value="ssuis/config.yaml"' in text
    assert 'value="matched_lineages/input.csv"' in text
    assert "Fictional format demonstration" in text
    # a run of the fictional collection: its tables and configuration copied into a folder of the run
    body, headers = _urlencoded({"example": "workflows/calls.yaml", "seed": "7", "action": "run"})
    status, response_headers, _ = _request(server, "POST", "/examples/run", body, headers)
    assert status == 303, status
    job_url = urllib.parse.urlsplit(response_headers["Location"]).path
    page = _wait(server, job_url)
    assert "Completed" in page and "Lineage share of the call" in page, page[-2000:]
    folder = next(p for p in (server["root"] / "runs").iterdir() if p.name.endswith("workflows_calls"))
    assert (folder / "data" / "calls.csv").is_file() and (folder / "data" / "metadata.csv").is_file()
    record = json.loads((folder / "results" / "clonal_share_result.json").read_text(encoding="utf-8"))
    # the same record the command line writes from the shipped file
    from amr_clonalshare import core
    again = core.run(load_config(EXAMPLES / "workflows" / "calls.yaml"), results_dir=tmp_path / "cli", seed=7)
    assert record["n_isolates"] == again["n_isolates"]
    for agent, block in record["metadata_diagnostics"]["clonal_share"].items():
        assert again["metadata_diagnostics"]["clonal_share"][agent]["kappa_adj"] == block["kappa_adj"]
    # and from the folder's own configuration
    cfg = load_config(folder / "config.yaml")
    assert cfg.data_root == (folder / "data").resolve()
    # an unknown example is refused
    body, headers = _urlencoded({"example": "nowhere/config.yaml", "seed": "7", "action": "run"})
    status, _, _ = _request(server, "POST", "/examples/run", body, headers)
    assert status == 404


@pytest.mark.skipif(EXAMPLES is None, reason="the examples folder of the source tree is not beside the package")
def test_a_configuration_file_of_the_readers_own_runs_in_place(server_with_examples):
    server = server_with_examples
    source = EXAMPLES / "workflows" / "mic.yaml"
    body, headers = _urlencoded({"config_path": str(source), "seed": "3", "action": "check"})
    status, response_headers, _ = _request(server, "POST", "/config/run", body, headers)
    assert status == 303, status
    job_url = urllib.parse.urlsplit(response_headers["Location"]).path
    page = _wait(server, job_url)
    assert "Completed" in page and "input_qc.json" in page, page[-1500:]
    folder = next(p for p in (server["root"] / "runs").iterdir() if p.name.endswith("_workflows"))
    cfg = load_config(folder / "config.yaml")
    # the tables stayed where the file names them
    assert cfg.data_root == (source.parent / "data").resolve()
    assert not any((folder / "data").iterdir())
    body, headers = _urlencoded({"config_path": str(source.parent / "absent.yaml"), "action": "check"})
    status, _, page = _request(server, "POST", "/config/run", body, headers)
    assert status == 400 and b"No configuration file" in page


@pytest.mark.skipif(EXAMPLES is None, reason="the examples folder of the source tree is not beside the package")
def test_the_comparison_example_opens_the_comparison_form_on_its_table(server_with_examples):
    server = server_with_examples
    body, headers = _urlencoded({"example": "matched_lineages/input.csv"})
    status, _, page = _request(server, "POST", "/examples/compare", body, headers)
    assert status == 200, page[:300]
    text = page.decode("utf-8")
    run = re.search(r'action="/compare/run\?run=([^&"]+)', text).group(1)
    assert run.endswith("matched_lineages_comparison")
    assert (server["root"] / "runs" / run / "data" / "input.csv").is_file()
    # the E. coli table opens with its columns already chosen: the call and the two lineage definitions
    body, headers = _urlencoded({"example": "ecoli_swine/input.csv"})
    status, _, page = _request(server, "POST", "/examples/compare", body, headers)
    assert status == 200, page[:300]
    text = page.decode("utf-8")
    for key, column in (("id_column", "isolate_id"), ("outcome", "ciprofloxacin_call"),
                        ("lineage_a", "sequence_type"), ("lineage_b", "phylogroup")):
        chosen = re.search(r'<select id="%s"[^>]*>(.*?)</select>' % key, text, re.S).group(1)
        assert '<option value="%s" selected>' % column in chosen


@pytest.mark.skipif(EXAMPLES is None, reason="the examples folder of the source tree is not beside the package")
def test_the_form_shows_the_first_rows_and_writes_the_tested_wells(server_with_examples):
    server = server_with_examples
    data = EXAMPLES / "workflows" / "data"
    body, headers = _multipart({"name": "wells"}, {
        "metadata": ("metadata.csv", (data / "metadata.csv").read_bytes()),
        "mic": ("mic.csv", (data / "mic.csv").read_bytes())})
    status, _, page = _request(server, "POST", "/upload", body, headers)
    assert status == 200, page[:400]
    text = page.decode("utf-8")
    # the first rows of each table are shown beside the choices, and the budgets are folded away with their defaults
    assert "The first 5 rows of <code>metadata.csv</code>" in text and "example_01_01" in text
    assert "<details" in text and 'name="mic_wells"' in text and 'name="alpha"' in text and 'step="any"' in text
    run = re.search(r'action="/run\?run=([^&"]+)', text).group(1)
    fields = {"strain_id_column": "isolate_id", "lineage_column": "lineage", "mic_id_column": "isolate_id",
              "mic_antibiotic_column": "agent", "mic_value_column": "measurement", "mic_operator_column": "operator",
              "mic_wells": "demo_agent: 0.25, 0.5, 1, 2, 4, 8", "n_perm": "19", "n_boot": "20", "folds": "2",
              "repeats": "2", "seed": "5", "action": "check"}
    body, headers = _urlencoded(fields)
    status, response_headers, page = _request(server, "POST", "/run", body, headers, query={"run": run})
    assert status == 303, page[:600]
    cfg = load_config(server["root"] / "runs" / run / "config.yaml")
    assert cfg.dataset.mic_wells == {"demo_agent": [0.25, 0.5, 1, 2, 4, 8]}
    page = _wait(server, urllib.parse.urlsplit(response_headers["Location"]).path)
    assert "Completed" in page and "Lineage groups" in page  # the input check is shown on the page itself
    # wells that cannot be read send the form back with the reason
    fields["mic_wells"] = "demo_agent 0.25 0.5"
    body, headers = _urlencoded(fields)
    status, _, page = _request(server, "POST", "/run", body, headers, query={"run": run})
    assert status == 400 and b"tested concentrations" in page


@pytest.mark.skipif(EXAMPLES is None, reason="the examples folder of the source tree is not beside the package")
def test_the_fictional_tables_are_served_as_templates(server_with_examples):
    server = server_with_examples
    status, headers, page = _request(server, "GET", "/examples/file", query={"name": "mic.csv"})
    assert status == 200 and page.startswith(b"isolate_id,agent,measurement,operator")
    assert "attachment" in headers["Content-Disposition"]
    status, _, _ = _request(server, "GET", "/examples/file", query={"name": "../config.yaml"})
    assert status == 404


def test_tested_wells_are_read_from_one_line_per_agent():
    assert gui._wells_from_text("cipro: 0.5, 1, 2\n\ngenta: 2 4 8\n") == {"cipro": [0.5, 1, 2], "genta": [2, 4, 8]}
    assert gui._wells_to_text({"cipro": [0.5, 1.0, 2.0]}) == "cipro: 0.5, 1, 2"
    with pytest.raises(gui.ConfigError):
        gui._wells_from_text("cipro: 0.5, one")


def test_a_damaged_or_two_sheet_workbook_is_refused_by_the_form_with_a_message(server, tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    import zipfile
    broken = tmp_path / "broken.xlsx"
    with zipfile.ZipFile(broken, "w") as archive:
        archive.writestr("not_a_workbook.txt", "x")
    calls = b"isolate_id,agent,call\na,tet,R\nb,tet,S\n"
    body, headers = _multipart({"name": "broken"}, {"metadata": ("broken.xlsx", broken.read_bytes()),
                                                    "calls": ("calls.csv", calls)})
    status, _, page = _request(server, "POST", "/upload", body, headers)
    assert status == 400 and b"cannot be read as an .xlsx workbook" in page
    book = openpyxl.Workbook()
    first = book.active
    first.title = "old_export"
    first.append(["isolate_id", "lineage"])
    first.append(["a", "L1"])
    second = book.create_sheet("current_export")
    second.append(["isolate_id", "lineage"])
    second.append(["a", "L2"])
    two = tmp_path / "two.xlsx"
    book.save(two)
    body, headers = _multipart({"name": "two"}, {"metadata": ("two.xlsx", two.read_bytes()),
                                                 "calls": ("calls.csv", calls)})
    status, _, page = _request(server, "POST", "/upload", body, headers)
    assert status == 400 and b"old_export, current_export" in page
    # the server answers the next request as before
    status, _, _ = _request(server, "GET", "/")
    assert status == 200


def _one_sheet(n_lineages=4, n_per=5):
    """One sheet: identifier, lineage and two agents in the wide shape with
    the sign inside the cell, as a laboratory spreadsheet is kept."""
    lines = ["isolate_id\tst\tCIP\tGEN"]
    for lineage in range(n_lineages):
        for k in range(n_per):
            cip = ["<=0.015", "0.03", "0.06", "0.12", "0.25", ">8"][(lineage + k) % 6]
            gen = ["0.5", "1", "2", "4", "8", "16"][(2 * lineage + k) % 6]
            lines.append(f"L{lineage}_{k}\tST{lineage}\t{cip}\t{gen}")
    return "\n".join(lines) + "\n"


def test_one_pasted_sheet_in_the_wide_shape_is_read_with_its_ranges_and_runs(server):
    sheet = _one_sheet()
    body, headers = _multipart({"name": "one sheet", "paste_mic": sheet}, {})
    status, _, page = _request(server, "POST", "/upload", body, headers)
    assert status == 200, page[:600]
    text = page.decode("utf-8")
    run = re.search(r'action="/run\?run=([^&"]+)', text).group(1)
    # the sheet is the metadata too, its agent columns are ticked, and the codes are read as names
    assert re.search(r'<select id="lineage_column"[^>]*>.*?<option value="st" selected>', text, re.S)
    assert 'name="mic_agent_columns" value="CIP" checked' in text
    assert 'name="mic_agent_columns" value="GEN" checked' in text
    assert "<code>CIP</code> read as ciprofloxacin" in text
    # the recorded concentrations, the button and the warning for the agent whose readings crowd the ends
    assert "Use the observed range" in text and "Panel preset" in text
    assert "eu-2020-1729-salmonella-ecoli" in text
    assert re.search(r"<td>CIP</td><td>0\.015, 0\.03, 0\.06, 0\.12, 0\.25, 8</td>", text)
    assert "Many readings lie on the lowest or highest recorded concentration" in text
    fields = {"strain_id_column": "isolate_id", "lineage_column": "st", "mic_id_column": "isolate_id",
              "mic_agent_columns": ["CIP", "GEN"], "mic_panel_preset": "eu-2020-1729-salmonella-ecoli",
              "n_perm": "20", "n_boot": "20", "folds": "2", "repeats": "2", "surveillance_n_boot": "20",
              "alpha": "0.05", "seed": "3", "action": "run"}
    body, headers = _urlencoded(fields)
    status, response_headers, page = _request(server, "POST", "/run", body, headers, query={"run": run})
    assert status == 303, page[:800]
    job_url = urllib.parse.urlsplit(response_headers["Location"]).path
    page = _wait(server, job_url)
    assert "Completed" in page, page[-1500:]
    assert "Conclusion for the MICs" in page and "The conclusion for every agent in full" in page
    folder = server["root"] / "runs" / run
    cfg = load_config(folder / "config.yaml")
    assert cfg.dataset.metadata == cfg.dataset.mic == "pasted_mic.csv"
    assert cfg.dataset.mic_agent_columns == ("CIP", "GEN")
    assert cfg.dataset.mic_wells["gentamicin"][0] == 0.5
    record = json.loads((folder / "results" / "clonal_share_result.json").read_text(encoding="utf-8"))
    agents = record["metadata_diagnostics"]["censored_share"]["per_agent"]
    assert set(agents) == {"ciprofloxacin", "gentamicin"}
    assert agents["ciprofloxacin"]["panel_source"].startswith("preset ")


def test_a_sheet_without_a_lineage_column_still_needs_the_metadata_table(server):
    body, headers = _multipart({"name": "no lineage"}, {
        "mic": ("mic.csv", b"isolate_id,CIP,GEN\nA1,0.5,2\nA2,1,4\nA3,2,8\nA4,0.25,1\n")})
    status, _, page = _request(server, "POST", "/upload", body, headers)
    assert status == 400 and b"carries the lineage column" in page


def test_a_wrong_column_comes_back_to_the_form_before_the_run_starts(server, tmp_path):
    from conftest import planted_cohort
    planted_cohort(tmp_path, n_per_lineage=4, n_lineages=4)
    data = tmp_path / "data"
    body, headers = _multipart({"name": "checked first"}, {
        "metadata": ("meta.csv", (data / "meta.csv").read_bytes()),
        "calls": ("calls.csv", (data / "calls.csv").read_bytes())})
    status, _, page = _request(server, "POST", "/upload", body, headers)
    run = re.search(r'action="/run\?run=([^&"]+)', page.decode()).group(1)
    body, headers = _urlencoded({"strain_id_column": "iid", "lineage_column": "lineage", "phenotype_id_column": "iid",
                                 "phenotype_antibiotic_column": "call", "phenotype_call_column": "antibiotic",
                                 "phenotype_kind": "clinical_sir", "action": "run"})
    status, _, page = _request(server, "POST", "/run", body, headers, query={"run": run})
    assert status == 400 and b"could not be read as named" in page and b'<select id="lineage_column"' in page


def test_the_run_page_offers_the_paragraph_the_folder_and_a_repeat_and_the_runs_page_sets_aside(server, tmp_path):
    from conftest import planted_cohort
    planted_cohort(tmp_path, n_per_lineage=5, n_lineages=5)
    data = tmp_path / "data"
    body, headers = _multipart({"name": "to repeat"}, {
        "metadata": ("meta.csv", (data / "meta.csv").read_bytes()),
        "calls": ("calls.csv", (data / "calls.csv").read_bytes())})
    status, _, page = _request(server, "POST", "/upload", body, headers)
    run = re.search(r'action="/run\?run=([^&"]+)', page.decode()).group(1)
    fields = {"strain_id_column": "iid", "lineage_column": "lineage", "phenotype_id_column": "iid",
              "phenotype_antibiotic_column": "antibiotic", "phenotype_call_column": "call", "phenotype_kind": "undeclared",
              "duplicate_policy": "drop_conflicts", "n_perm": "30", "n_boot": "20", "folds": "3", "repeats": "2",
              "surveillance_n_boot": "50", "alpha": "0.1", "seed": "11", "action": "run"}
    body, headers = _urlencoded(fields)
    status, response_headers, _ = _request(server, "POST", "/run", body, headers, query={"run": run})
    assert status == 303
    job_url = urllib.parse.urlsplit(response_headers["Location"]).path
    page = _wait(server, job_url)
    assert "A paragraph for a laboratory report" in page and "Copy the paragraph" in page
    assert "analysed with AMR-ClonalShare" in page and "Open the folder of this run" in page
    assert f'/?from={run}' in page or f"from={run}" in page
    # the same tables again, with the settings of that run proposed
    status, _, page = _request(server, "GET", "/", query={"from": run})
    assert status == 200 and b"Repeat with new tables" in page
    body, headers = _multipart({"name": "repeated", "from": run}, {
        "metadata": ("meta.csv", (data / "meta.csv").read_bytes()),
        "calls": ("calls.csv", (data / "calls.csv").read_bytes())})
    status, _, page = _request(server, "POST", "/upload", body, headers)
    text = page.decode("utf-8")
    assert status == 200
    assert re.search(r'<select id="duplicate_policy"[^>]*>.*?<option value="drop_conflicts" selected>', text, re.S)
    assert 'name="n_perm" value="30"' in text and 'name="seed" value="11"' in text and 'name="alpha" value="0.1"' in text
    # set the first run aside: its folder moves under _removed and nothing is deleted
    body, headers = _urlencoded({"folder": run})
    status, response_headers, _ = _request(server, "POST", "/runs/remove", body, headers)
    assert status == 303
    assert not (server["root"] / "runs" / run).exists()
    assert (server["root"] / "runs" / "_removed" / run / "config.yaml").is_file()
    status, _, page = _request(server, "GET", "/runs")
    assert status == 200 and b"Set aside" in page


def test_the_published_seed_of_a_shipped_collection_and_the_manual_route():
    assert gui._published_seed("ecoli_swine/config.yaml") == 20261001
    assert gui._published_seed("ssuis/config.yaml") == 42
    assert gui._published_seed("somewhere/else.yaml") == 42


def test_without_a_manual_folder_the_manual_page_says_where_the_manual_is(server):
    status, _, page = _request(server, "GET", "/manual")
    assert status == 404 and b"maciejkochanowski.github.io" in page
