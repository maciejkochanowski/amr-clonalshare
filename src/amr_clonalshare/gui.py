"""A local form for a run, for a laboratory that works without a command line.

``amr-clonalshare-gui`` starts a small web server on the loopback interface
of this computer, opens the browser on it, and lets the reader choose the
tables, name the columns, set the budgets and run the same analysis the
command line runs. Every run writes the same configuration file, the same
record and the same reports into a folder of its own under the working
folder, so that a run made here can be repeated from the command line with
``amr-clonalshare --config <run>/config.yaml``. The "Examples" page runs the
example collections shipped with the package (``--examples`` names their
folder; a source checkout finds its own), and any configuration file already
on this computer, as ``amr-clonalshare --config`` does. Nothing leaves the computer:
the server answers only the loopback address, and every request must carry
the token printed at start, which the page the browser opens carries.

The module uses the standard library only: ``http.server`` for the server,
``email`` for the uploaded forms, ``webbrowser`` for the browser. The pages
carry no script beyond three buttons that copy a text or fill a field; the
page of a running job refreshes itself.
"""
from __future__ import annotations

import argparse
import contextlib
import datetime as _dt
import email
import email.policy
import html
import io
import json
import mimetypes
import os
import re
import secrets
import socketserver
import sys
import threading
import time
import traceback
import urllib.parse
import webbrowser
import zipfile
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import warnings

import yaml

from . import __version__
from .cli import examples_folder
from .config import ConfigError, load_config
from .agents import unify_agents
from .draft import _MIC_OPERATOR_PATTERNS, _MIC_VALUE_PATTERNS, _PANEL_PATTERNS, _ROLES, _columns, _kind, _match, _mic_like, _rows, observed_wells
from .panels import PANEL_PRESETS
from .qc import END_WELL_WARNING
from .verdict import VERDICT_NOTE, short_label, summary_rows

__all__ = ["main", "serve", "Session"]

_PATTERNS = dict(_ROLES)
_MIC_PATTERNS = {
    "mic_id_column": _PATTERNS["strain_id_column"],
    "mic_antibiotic_column": _PATTERNS["phenotype_antibiotic_column"],
    "mic_value_column": _MIC_VALUE_PATTERNS,
    "mic_operator_column": _MIC_OPERATOR_PATTERNS,
    "mic_unit_column": ("unit",),
    "mic_panel_column": _PANEL_PATTERNS,
}
_KINDS = (("undeclared", "read from the words used (resistant, susceptible, ...)"),
          ("clinical_sir", "clinical S / I / R"),
          ("wt_nwt", "wild-type / non-wild-type"),
          ("binary", "0 / 1"))
_INTERMEDIATE = (("", "default of the vocabulary"), ("non_susceptible", "counted as positive"),
                 ("drop", "set aside"), ("susceptible", "counted as negative"))
_DUPLICATES = (("error", "stop the run"), ("drop_conflicts", "set the conflicting key aside"),
               ("positive_wins", "a positive record wins"))
#: Every column setting of the dataset section that the form offers, with
#: the table it is read from. tests/test_gui.py checks this list against
#: DatasetConfig, so a new key cannot be added to the configuration without
#: a decision about the form.
COLUMN_FIELDS: Tuple[Tuple[str, str], ...] = (
    ("strain_id_column", "metadata"),
    ("lineage_column", "metadata"),
    ("batch_column", "metadata"),
    ("contrast_column", "metadata"),
    ("phenotype_covariate_column", "metadata"),
    ("stratify_by", "metadata"),
    ("unit_column", "metadata"),
    ("phenotype_id_column", "calls"),
    ("phenotype_antibiotic_column", "calls"),
    ("phenotype_call_column", "calls"),
    ("phenotype_agent_columns", "calls"),
    ("mic_id_column", "mic"),
    ("mic_antibiotic_column", "mic"),
    ("mic_value_column", "mic"),
    ("mic_operator_column", "mic"),
    ("mic_agent_columns", "mic"),
    ("mic_unit_column", "mic"),
    ("mic_panel_column", "mic"),
    ("mic_covariate_columns", "mic+metadata"),
)
#: Settings of the dataset section entered as text or chosen from a list.
TEXT_FIELDS: Tuple[str, ...] = (
    "name", "contrast_levels", "phenotype_kind", "phenotype_intermediate",
    "duplicate_policy", "phenotype_positive_definition", "phenotype_source",
    "ast_standard", "ast_version", "mic_wells", "mic_panel_preset")
#: Settings the form leaves to the configuration file: the unit per agent is
#: rarely needed, since a unit column covers it, and the manual says how to
#: write it.
FILE_ONLY_FIELDS: Tuple[str, ...] = ("data_dir", "metadata", "phenotype", "mic", "mic_units")


# ---------------------------------------------------------------- the session and its jobs


class Job:
    """One analysis or comparison, run in a thread of its own."""

    def __init__(self, kind: str, name: str, folder: Path) -> None:
        self.id = secrets.token_hex(6)
        self.kind = kind
        self.name = name
        self.folder = folder
        self.results = folder / "results"
        self.status = "queued"
        self.lines: List[str] = []
        self.started = time.time()
        self.finished: Optional[float] = None
        self.error: Optional[str] = None
        self.summary: Dict[str, object] = {}
        self.seed: Optional[int] = None
        self.lock = threading.Lock()

    def log(self, message: str) -> None:
        stamp = _dt.datetime.now().strftime("%H:%M:%S")
        with self.lock:
            self.lines.append(f"{stamp}  {message}")
        try:
            with (self.folder / "gui.log").open("a", encoding="utf-8") as handle:
                handle.write(f"{stamp}  {message}\n")
        except OSError:
            pass

    def finish(self, status: str, error: Optional[str] = None) -> None:
        self.status = status
        self.error = error
        self.finished = time.time()
        record = {"id": self.id, "kind": self.kind, "name": self.name, "status": status,
                  "error": error, "started": self.started, "finished": self.finished,
                  "seed": self.seed, "results": self.results.name, "summary": self.summary}
        try:
            (self.folder / "job.json").write_text(json.dumps(record, indent=1, default=str) + "\n",
                                                  encoding="utf-8")
        except OSError:
            pass

    @property
    def seconds(self) -> float:
        return (self.finished or time.time()) - self.started


class Session:
    """The working folder, the token and the jobs of one server."""
    manual: Optional[Path] = None

    def __init__(self, root: Path, threads: Optional[int] = None, examples: Optional[Path] = None) -> None:
        self.root = root
        self.runs = root / "runs"
        self.runs.mkdir(parents=True, exist_ok=True)
        self.token = secrets.token_urlsafe(24)
        self.threads = threads
        self.examples = examples.resolve() if examples and examples.is_dir() else None
        self.jobs: Dict[str, Job] = {}
        self.jobs_lock = threading.Lock()

    def example_runs(self) -> List[dict]:
        """The example collections: every configuration file one folder below
        the examples folder, and every folder holding an ``input.csv`` for the
        comparison of two lineage definitions, each with the opening comment
        of its file or the opening paragraph of its README as description."""
        found: List[dict] = []
        comparisons: List[dict] = []
        if self.examples is None:
            return found
        # the small fictional collection first, since it runs in seconds; the comparisons last
        folders = sorted((p for p in self.examples.iterdir() if p.is_dir()), key=lambda p: (p.name != "workflows", p.name))
        for folder in folders:
            for path in sorted(folder.glob("*.yaml")):
                found.append({"kind": "analysis", "id": f"{folder.name}/{path.name}", "path": path,
                              "title": _example_title(folder, path), "description": _opening_comment(path),
                              "seed": _published_seed(f"{folder.name}/{path.name}")})
            table = folder / "input.csv"
            if table.is_file():
                comparisons.append({"kind": "comparison", "id": f"{folder.name}/input.csv", "path": table,
                                    "title": _example_title(folder, table), "description": _opening_paragraph(folder / "README.md")})
        return found + comparisons

    def example(self, key: str) -> Optional[dict]:
        for item in self.example_runs():
            if item["id"] == key:
                return item
        return None

    def new_folder(self, kind: str, name: str) -> Path:
        stamp = _dt.datetime.now().strftime("%Y-%m-%d_%H%M%S")
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip()) or kind
        folder = self.runs / f"{stamp}_{safe}"
        n = 1
        while folder.exists():
            n += 1
            folder = self.runs / f"{stamp}_{safe}_{n}"
        (folder / "data").mkdir(parents=True)
        return folder

    def start(self, job: Job, target: Callable[[Job], None]) -> Job:
        with self.jobs_lock:
            self.jobs[job.id] = job

        def run() -> None:
            job.status = "running"
            try:
                target(job)
                job.finish("done")
            except ConfigError as error:
                job.log(f"input error: {error}")
                job.finish("failed", f"input error: {error}")
            except Exception as error:  # a failed run is reported, not raised into the server
                job.log(traceback.format_exc().strip().splitlines()[-1])
                job.finish("failed", f"{type(error).__name__}: {error}")

        threading.Thread(target=run, name=f"amr-clonalshare-{job.id}", daemon=True).start()
        return job

    def live_job(self, folder: Path) -> Optional[Job]:
        """The job still running in ``folder``, if any."""
        with self.jobs_lock:
            jobs = list(self.jobs.values())
        return next((j for j in jobs if j.folder == folder and j.status in ("queued", "running")), None)

    def previous_runs(self) -> List[dict]:
        found = []
        for folder in sorted(self.runs.iterdir(), reverse=True):
            record = folder / "job.json"
            if record.is_file():
                try:
                    item = json.loads(record.read_text(encoding="utf-8"))
                except (OSError, ValueError):
                    continue
                item["folder"] = folder
                found.append(item)
        return found


#: Plain titles for the example collections of the repository; any other folder is named by itself.
_EXAMPLE_TITLES = {
    "workflows": "Fictional teaching collection, 96 isolates",
    "ssuis": "Streptococcus suis, 677 isolates, calls and MICs",
    "ecoli_swine": "Escherichia coli from diseased pigs, 481 isolates, MICs",
    "salmonella_poultry": "Salmonella from poultry meat, 7,049 isolates, calls",
    "matched_lineages": "Two lineage definitions on matched records",
}
_EXAMPLE_FILES = {
    "calls": "susceptibility calls and lineage labels",
    "mic": "recorded MICs without calls",
    "contrasts": "two collections compared",
    "config": "the main analysis",
    "config_contrast": "the same collection at sequence-type resolution, two countries contrasted",
    "config_cluster": "the same isolates at SNP-cluster resolution",
    "config_two_categories": "the same readings merged into two categories per agent",
    "config_without_st410": "the same readings without the most frequent sequence type",
    "config_periods": "calls above a cut per agent by phylogroup, 2014-2017 against 2011-2013",
}


def _published_seed(relative: str) -> int:
    """The seed of the published record of a shipped collection, so that the
    numbers on the page agree with the appendices of the article; 42 for a
    collection without one."""
    from .cli import EXAMPLES
    for config, seed in EXAMPLES.values():
        if config == relative:
            return seed
    return 42


def _example_title(folder: Path, path: Path) -> str:
    head = _EXAMPLE_TITLES.get(folder.name, folder.name)
    tail = _EXAMPLE_FILES.get(path.stem, path.stem) if path.suffix == ".yaml" else "the table of the comparison"
    return f"{head}: {tail}"


def _opening_comment(path: Path) -> str:
    """The first paragraph of the comment a configuration file opens with,
    without the command lines it quotes."""
    lines: List[str] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not raw.startswith("#"):
            break
        line = raw.lstrip("#").strip()
        if not line:
            if lines:
                break
            continue
        if line.startswith(("python ", "amr-clonalshare")):
            continue
        lines.append(line)
    return " ".join(lines)


def _opening_paragraph(path: Path) -> str:
    """The first paragraph of a README after its heading."""
    if not path.is_file():
        return ""
    lines: List[str] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if line.startswith("#"):
            continue
        if not line:
            if lines:
                break
            continue
        lines.append(line)
    return " ".join(lines)


_FILE_KEYS = ("metadata", "phenotype", "mic")


def _safe_name(name: str) -> bool:
    """A folder name of the working folder: no separator, no parent."""
    return bool(re.fullmatch(r"[A-Za-z0-9_.-]+", name or "")) and name not in (".", "..")


def _seed(form: Dict[str, List[str]]) -> int:
    seed = _number(form, "seed", 42)
    if seed < 0:
        raise ConfigError(f"seed must be a non-negative integer, got {seed}")
    return seed


def _run_folder_from_config(session: Session, source: Path, name: str, copy_data: bool, seed: int = 42) -> Path:
    """A run folder holding ``config.yaml`` for the configuration file at
    ``source``. With ``copy_data`` the tables the file names are copied under
    ``data/`` and the folder stands alone, as the form's own runs do; without
    it the configuration keeps pointing at the tables where they are, by the
    absolute path of their folder, so the run reads the reader's files in
    place. Either way the folder's configuration is the one the command line
    repeats the run from."""
    # the file is read as the command line reads it: a duplicated key is refused
    load_config(source, check_files_exist=False)
    text = source.read_text(encoding="utf-8")
    raw = yaml.safe_load(text)
    if not isinstance(raw, dict) or not isinstance(raw.get("dataset"), dict):
        raise ConfigError(f"{source} holds no dataset section")
    dataset = raw["dataset"]
    data_root = (source.parent / str(dataset.get("data_dir") or ".")).resolve()
    folder = session.new_folder("analysis", name)
    header = (f"# Written by amr-clonalshare-gui from {source}; the same run from the command line:\n"
              f"#   amr-clonalshare --config config.yaml --results-dir <folder> --seed {seed}\n")
    if copy_data:
        copied: Dict[str, str] = {}
        for key in _FILE_KEYS:
            value = dataset.get(key)
            if value:
                origin = data_root / str(value)
                if not origin.is_file():
                    raise ConfigError(f"{source} names {key}: {value}, which is not a file under {data_root}")
                if Path(str(value)).name in copied and copied[Path(str(value)).name] != str(origin):
                    raise ConfigError(f"{source} names two tables called {Path(str(value)).name}; "
                                      "rename one so that both can be copied into the run")
                copied[Path(str(value)).name] = str(origin)
                target = folder / "data" / Path(str(value)).name
                target.write_bytes(origin.read_bytes())
                dataset[key] = target.name
        dataset["data_dir"] = "data"
    else:
        dataset["data_dir"] = str(data_root)
    (folder / "config.yaml").write_text(header + yaml.safe_dump(raw, sort_keys=False, allow_unicode=True),
                                        encoding="utf-8")
    load_config(folder / "config.yaml")
    return folder


# ---------------------------------------------------------------- the analysis a form describes


def _thread_limit(threads: Optional[int]):
    """The processor threads of a job: the numerical stack is already
    imported by the form, so the limit is set on its thread pools."""
    if not threads:
        return contextlib.nullcontext()
    from threadpoolctl import threadpool_limits
    return threadpool_limits(limits=threads)


def _analysis_target(session: Session, seed: int, check_only: bool) -> Callable[[Job], None]:
    def target(job: Job) -> None:
        with _thread_limit(session.threads):
            _analysis(job, seed, check_only)
    return target


def _analysis(job: Job, seed: int, check_only: bool) -> None:
    """One analysis of a run folder, or its input check alone."""
    from . import core
    from .cli import _summary
    from .io import load_dataset
    from .outputs import publish_input_check
    cfg = load_config(job.folder / "config.yaml")
    # a second run of the same folder (after an input check, say) replaces
    # the results files of the first
    again = job.results.is_dir() and any(job.results.iterdir())
    if check_only:
        job.log("reading the tables")
        dataset = load_dataset(cfg)
        if dataset.input_qc is None:
            raise ConfigError("the tables gave no input check")
        publish_input_check(dataset.input_qc, job.results, overwrite=again)
        job.summary = {"n_isolates": dataset.input_qc.get("n_isolates"),
                       "n_antimicrobials": dataset.input_qc.get("n_antimicrobials")}
        job.log("input check written")
        return
    log: List[str] = []
    with core.log_warnings(log):
        record = core.run(cfg, results_dir=job.results, seed=seed, progress=job.log,
                          overwrite=again, log=log)
    job.summary = _summary(record)
    job.log("results written")


def _comparison_target(session: Session, arguments: List[str]) -> Callable[[Job], None]:
    def target(job: Job) -> None:
        with _thread_limit(session.threads):
            from .comparison import run_comparison
            job.log("comparing the two lineage definitions")
            run_comparison(arguments + ["--output", str(job.results)])
            job.log("comparison written")
    return target


def _dataset_from_form(form: Dict[str, List[str]], files: Dict[str, str]) -> dict:
    """The dataset section of a configuration, from the fields of the form."""
    def one(key: str) -> Optional[str]:
        values = [v.strip() for v in form.get(key, []) if v.strip()]
        return values[0] if values else None

    def many(key: str) -> List[str]:
        return [v.strip() for v in form.get(key, []) if v.strip()]

    dataset: Dict[str, object] = {"name": one("name") or "collection", "data_dir": "data"}
    dataset["metadata"] = files.get("metadata") or files.get("calls") or files.get("mic")
    if "calls" in files:
        dataset["phenotype"] = files["calls"]
    if "mic" in files:
        dataset["mic"] = files["mic"]
    for key, table in COLUMN_FIELDS:
        if table == "calls" and "calls" not in files:
            continue
        if table.startswith("mic") and "mic" not in files:
            continue
        if key in ("phenotype_agent_columns", "mic_agent_columns", "mic_covariate_columns"):
            values = many(key)
            if values:
                dataset[key] = values
        else:
            value = one(key)
            if value is not None:
                dataset[key] = value
    if dataset.get("phenotype_agent_columns"):
        dataset.pop("phenotype_call_column", None)
        dataset.pop("phenotype_antibiotic_column", None)
    if dataset.get("mic_agent_columns"):
        for key in ("mic_antibiotic_column", "mic_value_column", "mic_operator_column"):
            dataset.pop(key, None)
    for key in ("phenotype_kind", "duplicate_policy"):
        value = one(key)
        if value:
            dataset[key] = value
    for key in ("phenotype_intermediate", "phenotype_positive_definition", "phenotype_source",
                "ast_standard", "ast_version"):
        value = one(key)
        if value:
            dataset[key] = value
    if "mic" in files:
        wells = _wells_from_text(one("mic_wells") or "")
        if wells:
            dataset["mic_wells"] = wells
        preset = one("mic_panel_preset")
        if preset:
            dataset["mic_panel_preset"] = preset
    levels = [v.strip() for v in (one("contrast_levels") or "").split(";") if v.strip()]
    if dataset.get("contrast_column"):
        dataset["contrast_levels"] = levels
    else:
        dataset.pop("contrast_column", None)
    return dataset


def _wells_from_text(text: str) -> Dict[str, List[float]]:
    """The tested concentrations per agent, one agent per line: the agent, a
    colon, and its concentrations separated by commas or spaces, as in
    ``ciprofloxacin: 0.5, 1, 2``."""
    wells: Dict[str, List[float]] = {}
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        agent, sep, rest = line.partition(":")
        if not sep or not agent.strip():
            raise ConfigError(f"tested concentrations: write the agent, a colon and its concentrations; got {line!r}")
        try:
            values = [float(v) for v in re.split(r"[,;\s]+", rest.strip()) if v]
        except ValueError as error:
            raise ConfigError(f"tested concentrations of {agent.strip()}: {rest.strip()!r} is not a list of numbers") from error
        if not values:
            raise ConfigError(f"tested concentrations of {agent.strip()}: none given")
        wells[agent.strip()] = values
    return wells


def _wells_to_text(wells) -> str:
    if isinstance(wells, str):
        return wells
    if not isinstance(wells, dict):
        return ""
    return "\n".join(f"{agent}: {', '.join(_plain(v) for v in values)}" for agent, values in wells.items())


def _plain(value) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    return str(int(number)) if number.is_integer() else f"{number:g}"


def _ranges_block(ranges: Sequence[dict]) -> str:
    """The concentrations every agent was recorded at, a button that writes
    them into the field above, and a warning where many readings sit on the
    lowest or highest of them while no tested range is given."""
    if not ranges:
        return ""
    rows = []
    warned = []
    for r in ranges:
        share = r.get("share_on_ends")
        shown = pct_text(share)
        rows.append(f"<tr><td>{_e(r['agent'])}</td><td>{_e(', '.join(_plain(w) for w in r['wells']))}</td>"
                    f'<td class="num">{r["n"]}</td><td class="num">{shown}</td></tr>')
        if share is not None and share >= END_WELL_WARNING:
            warned.append(f"{r['agent']} ({shown})")
    text = "\n".join(f"{r['agent']}: {', '.join(_plain(w) for w in r['wells'])}" for r in ranges if r["wells"])
    warning = ""
    if warned:
        warning = (f'<div class="error">Many readings lie on the lowest or highest recorded concentration: {_e(", ".join(warned))}. '
                   f"Without the tested range those concentrations are taken as the end wells of the panel and the readings on them as censored. "
                   f"If the panel tested further dilutions, write its range above or choose a panel preset; "
                   f"if the recorded concentrations are the whole panel, use the observed range.</div>")
    return (f'<div class="field"><span class="lab">Recorded concentrations</span><div>'
            f'<table class="ranges"><thead><tr><th>Agent</th><th>Recorded at (mg/L)</th><th class="num">Readings</th>'
            f'<th class="num">On lowest or highest</th></tr></thead><tbody>{"".join(rows)}</tbody></table>'
            f'<button type="button" class="quiet" data-wells="{_e(text)}" '
            f'onclick="document.getElementById(\'mic_wells\').value=this.dataset.wells">Use the observed range</button>'
            f'<div class="note">What the table holds for every agent. "Use the observed range" writes it into the field above; '
            f"the plate sheet of the laboratory may test more dilutions than the readings show.</div>{warning}</div></div>")


def pct_text(share) -> str:
    return "" if share is None else f"{100 * share:.0f} %"


def _preview(path: Path, limit: int = 5) -> str:
    """The first rows of a table, for the reader to see what each column holds."""
    try:
        columns = _columns(path)
        rows = _rows(path, limit=limit)
    except (OSError, ValueError, UnicodeDecodeError):
        return ""
    head = "".join(f"<th>{_e(c)}</th>" for c in columns)
    body = "".join("<tr>" + "".join(f"<td>{_e(_clip(row.get(c, '')))}</td>" for c in columns) + "</tr>"
                   for row in rows[:limit])
    return (f'<div class="preview"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'
            f'<p class="note">The first {min(limit, len(rows))} rows of <code>{_e(path.name)}</code>.</p></div>')


def _clip(value, width: int = 28) -> str:
    text = "" if value is None else str(value)
    return text if len(text) <= width else text[:width - 1] + "…"


def _number(form: Dict[str, List[str]], key: str, default: float, kind=int):
    values = [v.strip() for v in form.get(key, []) if v.strip()]
    if not values:
        return default
    try:
        return kind(values[0])
    except ValueError as error:
        raise ConfigError(f"{key} must be a number; got {values[0]!r}") from error


def _config_from_form(form: Dict[str, List[str]], files: Dict[str, str]) -> dict:
    raw: Dict[str, object] = {"dataset": _dataset_from_form(form, files)}
    raw["attribution"] = {"folds": _number(form, "folds", 5), "repeats": _number(form, "repeats", 20),
                          "n_boot": _number(form, "n_boot", 999), "n_perm": _number(form, "n_perm", 999)}
    raw["surveillance"] = {"n_boot": _number(form, "surveillance_n_boot", 2000)}
    raw["censored"] = {"end_wells_censored": form.get("end_wells_censored", ["on"]) != ["off"]}
    raw["evidence"] = {"alpha": _number(form, "alpha", 0.05, float)}
    return raw


# ---------------------------------------------------------------- guesses for the form


def _guess(columns: Sequence[str], key: str) -> Optional[str]:
    patterns = _PATTERNS.get(key) or _MIC_PATTERNS.get(key)
    return _match(columns, patterns) if patterns else None


def _table_guesses(paths: Dict[str, Path]) -> Dict[str, Any]:
    """Columns per table, one guess per role and the vocabulary of the calls."""
    columns = {table: _columns(path) for table, path in paths.items()}
    guesses: Dict[str, Any] = {}
    meta = columns.get("metadata", [])
    calls = columns.get("calls", [])
    mic = columns.get("mic", [])
    # The identifier is the column every table shares, when there is one;
    # the names alone can prefer a column called "strain" over the join key.
    shared = [c for c in meta if all(c in columns[t] for t in columns)] if len(columns) > 1 else []
    identifier = shared[0] if shared else _guess(meta, "strain_id_column")
    guesses["strain_id_column"] = identifier
    guesses["lineage_column"] = _guess(meta, "lineage_column")
    guesses["batch_column"] = _match(meta, ("year", "intake", "batch", "period"))
    names: Dict[str, str] = {}
    if calls:
        guesses["phenotype_id_column"] = identifier if identifier in calls else _guess(calls, "strain_id_column")
        guesses["phenotype_antibiotic_column"] = _guess(calls, "phenotype_antibiotic_column")
        guesses["phenotype_call_column"] = _guess(calls, "phenotype_call_column")
        kind = None
        agents: List[str] = []
        rows = _rows(paths["calls"])
        if guesses["phenotype_call_column"]:
            kind = _kind([row.get(str(guesses["phenotype_call_column"]), "") for row in rows])
        else:
            taken = {guesses["phenotype_id_column"], guesses["lineage_column"]} - {None}
            kinds = {c: _kind([row.get(c, "") for row in rows]) for c in calls if c not in taken}
            found = {k for k in kinds.values() if k}
            if len(found) == 1:
                kind = found.pop()
                agents = [c for c, k in kinds.items() if k == kind]
        guesses["phenotype_kind"] = kind or "undeclared"
        guesses["phenotype_agent_columns"] = agents
        if agents:
            names.update(unify_agents(agents)[1])
        elif guesses["phenotype_antibiotic_column"]:
            names.update(unify_agents(row.get(str(guesses["phenotype_antibiotic_column"]), "") for row in rows)[1])
    ranges: List[dict] = []
    if mic:
        rows = _rows(paths["mic"], limit=10 ** 7)
        taken = {identifier, guesses["lineage_column"]} - {None}
        fits = [c for c in mic if c not in taken and _mic_like([row.get(c, "") for row in rows])]
        for key in _MIC_PATTERNS:
            guesses[key] = _guess([c for c in mic if c not in fits], key)
        if identifier in mic:
            guesses["mic_id_column"] = identifier
        guesses["mic_agent_columns"] = []
        if guesses["mic_antibiotic_column"] and len(fits) == 1:
            guesses["mic_value_column"] = fits[0]
        elif fits and not guesses["mic_antibiotic_column"]:
            guesses["mic_agent_columns"] = fits
            guesses["mic_value_column"] = None
            guesses["mic_operator_column"] = None
        if guesses["mic_agent_columns"]:
            ranges = [_range(c, [row.get(c, "") for row in rows]) for c in guesses["mic_agent_columns"]]
            names.update(unify_agents(guesses["mic_agent_columns"])[1])
        elif guesses["mic_antibiotic_column"] and guesses["mic_value_column"]:
            by_agent: Dict[str, List[str]] = {}
            for row in rows:
                by_agent.setdefault(str(row.get(str(guesses["mic_antibiotic_column"]), "")).strip(), []).append(
                    row.get(str(guesses["mic_value_column"]), ""))
            ranges = [_range(a, v) for a, v in by_agent.items() if a]
            names.update(unify_agents(by_agent)[1])
    return {"columns": columns, "guesses": guesses, "mic_ranges": ranges, "agent_names": names}


def _carry_settings(config_path: Path, info: Dict[str, Any]) -> None:
    """The settings of an earlier run proposed for new tables: every column
    choice the new tables can carry, the settings that name no column, and
    the budgets and the seed."""
    try:
        raw = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError):
        return
    columns: Dict[str, List[str]] = info["columns"]
    guesses: Dict[str, Any] = info["guesses"]
    dataset = raw.get("dataset") or {}
    for key, table in COLUMN_FIELDS:
        value = dataset.get(key)
        if value is None:
            continue
        available = set(columns.get(table.split("+")[0], [])) | (set(columns.get("metadata", [])) if "+" in table else set())
        if isinstance(value, list):
            kept = [c for c in value if c in available]
            if kept:
                guesses[key] = kept
        elif value in available:
            guesses[key] = value
    for key in ("phenotype_kind", "phenotype_intermediate", "duplicate_policy", "phenotype_positive_definition",
                "phenotype_source", "ast_standard", "ast_version", "mic_wells", "mic_panel_preset"):
        if dataset.get(key) is not None:
            guesses[key] = dataset[key]
    if dataset.get("contrast_levels"):
        guesses["contrast_levels"] = "; ".join(str(v) for v in dataset["contrast_levels"])
    attribution = raw.get("attribution") or {}
    budgets = {k: attribution[k] for k in ("n_perm", "n_boot", "folds", "repeats") if k in attribution}
    if (raw.get("surveillance") or {}).get("n_boot") is not None:
        budgets["surveillance_n_boot"] = raw["surveillance"]["n_boot"]
    if (raw.get("evidence") or {}).get("alpha") is not None:
        budgets["alpha"] = raw["evidence"]["alpha"]
    if (raw.get("censored") or {}).get("end_wells_censored") is False:
        budgets["end_wells_censored"] = "off"
    comment = config_path.read_text(encoding="utf-8", errors="replace")
    seed = re.search(r"--seed (\d+)", comment)
    if seed:
        budgets["seed"] = int(seed.group(1))
    guesses["budgets"] = budgets


def _range(agent: str, values: Sequence[str]) -> dict:
    """The concentrations an agent was recorded at and the share of its
    readings on the lowest or highest of them."""
    from .draft import _concentration
    numbers = [n for n in (_concentration(v) for v in values) if n is not None]
    wells = observed_wells(values)
    ends = sum(1 for n in numbers if wells and (abs(n - wells[0]) < 1e-9 or abs(n - wells[-1]) < 1e-9))
    return {"agent": agent, "wells": wells, "n": len(numbers),
            "share_on_ends": (ends / len(numbers)) if numbers else None}


# ---------------------------------------------------------------- the pages

_CSS = """
:root{--ink:#000;--mid:#3D3D3D;--faint:#6E6E6E;--rule:#BFBFBF;--paper:#FFF;--accent:#0B2545;--flag:#A63D40;--field:#F7F7F5}
*{box-sizing:border-box}
html{background:#F2F2F2}
body{margin:0;background:var(--paper);color:var(--ink);font-family:@TEXT@;font-size:16px;line-height:1.45;font-variant-numeric:lining-nums tabular-nums;-webkit-font-smoothing:antialiased}
.sheet{max-width:46rem;min-height:100vh;margin:0 auto;padding:30px 30px 70px}
header{display:flex;justify-content:space-between;align-items:baseline;border-bottom:.8px solid var(--ink);padding-bottom:8px;margin-bottom:22px}
header h1{font-size:1.05rem;font-weight:600;margin:0;letter-spacing:-.005em}
header h1 span{color:var(--faint);font-weight:400;margin-left:.6em}
nav a{font-size:.8rem;letter-spacing:.09em;text-transform:uppercase;color:var(--mid);text-decoration:none;margin-left:1.4em}
nav a.here{color:var(--accent);font-weight:600}
h2{font-size:1rem;font-weight:600;margin:30px 0 8px;padding-bottom:4px;border-bottom:.8px solid var(--ink)}
h2 .no{color:var(--faint);font-weight:400;margin-right:.55em}
p{margin:0 0 10px}
.lead{color:var(--mid)}
.field{display:grid;grid-template-columns:14em 1fr;gap:4px 16px;padding:7px 0;border-bottom:.5px solid var(--rule);align-items:baseline}
.field:last-child{border-bottom:0}
.field label,.field .lab{font-size:.78rem;letter-spacing:.055em;text-transform:uppercase;color:var(--faint);padding-top:3px}
.field .note{grid-column:2;font-size:.8rem;color:var(--mid);margin-top:1px}
select,input[type=text],input[type=number],textarea{font:inherit;font-size:.92rem;color:var(--ink);background:var(--field);border:.8px solid var(--rule);border-radius:0;padding:5px 7px;width:100%;max-width:28em}
select[multiple]{height:7.5em}
input[type=file]{font-size:.88rem}
input[type=number]{max-width:8em}
textarea{max-width:100%;min-height:3.2em}
.actions{margin-top:26px;display:flex;gap:12px;align-items:center}
button{font:inherit;font-size:.85rem;letter-spacing:.06em;text-transform:uppercase;padding:8px 16px;background:var(--accent);color:#fff;border:.8px solid var(--accent);border-radius:0;cursor:pointer}
button.quiet{background:var(--paper);color:var(--accent)}
button.small{padding:3px 8px;font-size:.72rem}
form.inline{display:inline;margin:0}
a.button{font-size:.85rem;letter-spacing:.06em;text-transform:uppercase;padding:8px 16px;background:var(--paper);color:var(--accent);border:.8px solid var(--accent);text-decoration:none}
.status{margin:0 0 18px}
.status .kicker{font-size:.8rem;letter-spacing:.1em;text-transform:uppercase;color:var(--accent);font-weight:600}
.status .kicker.flag{color:var(--flag)}
.status p{margin:4px 0 0;color:var(--mid)}
table{border-collapse:collapse;width:100%;margin:10px 0 6px;font-size:.86rem}
thead th{border-top:.8px solid var(--ink);border-bottom:.5px solid var(--ink);padding:5px 7px;font-weight:600;text-align:left;vertical-align:bottom}
tbody td{padding:4px 7px;border:0;vertical-align:top}
tbody tr:last-child td{border-bottom:.8px solid var(--ink)}
.num{text-align:right;white-space:nowrap}
code,.mono{font-family:@MONO@;font-size:.92em}
pre{font-family:@MONO@;font-size:.8rem;line-height:1.45;background:var(--field);border:.5px solid var(--rule);padding:10px 12px;overflow:auto;white-space:pre-wrap}
.files a{display:block;padding:5px 0;border-bottom:.5px solid var(--rule);color:var(--accent);text-decoration:none}
.files a:last-child{border-bottom:0}
.files small{color:var(--faint);margin-left:.6em}
.error{border-left:3px solid var(--flag);padding:8px 12px;margin:14px 0;background:var(--field)}
.colophon{margin-top:40px;border-top:.5px solid var(--rule);padding-top:9px;font-size:.76rem;color:var(--faint);line-height:1.5}
a{color:var(--accent)}
details.more{margin:6px 0 2px;border-bottom:.5px solid var(--rule)}
details.more summary{cursor:pointer;font-size:.8rem;letter-spacing:.055em;text-transform:uppercase;color:var(--accent);padding:8px 0}
details.more[open] summary{border-bottom:.5px solid var(--rule)}
.boxes{display:flex;flex-wrap:wrap;gap:4px 14px;max-width:28em}
.boxes .box{font-size:.9rem;color:var(--ink);white-space:nowrap;text-transform:none;letter-spacing:0}
p.note{font-size:.8rem;color:var(--mid)}
.preview{margin:4px 0 10px;overflow:auto}
.preview table{font-size:.78rem;margin:4px 0 2px;width:auto;min-width:50%}
.preview td,.preview th{white-space:nowrap;padding:3px 8px}
.preview .note{font-size:.76rem;color:var(--faint);margin:0}
table.shapes{font-size:.84rem;margin:4px 0 14px}
table.shapes td:first-child{font-weight:600;white-space:nowrap}
form.row{display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin:8px 0 4px}
form.row label.inline{font-size:.78rem;letter-spacing:.055em;text-transform:uppercase;color:var(--faint);margin-left:10px}
form.row input[type=number]{max-width:7em}
.wide{overflow:auto}
.wide table{font-size:.78rem}
.wide td:first-child{white-space:nowrap}
.wide thead th{padding:4px 5px;line-height:1.2;white-space:normal}
.wide tbody td{padding:4px 5px}
.wide td.num{white-space:nowrap}
iframe.report{width:100%;height:60rem;border:.5px solid var(--rule);background:#fff}
"""


def _style() -> str:
    """The style sheet of every page, with its typefaces embedded."""
    from .typeface import MONO, TEXT, font_faces
    return font_faces() + _CSS.replace("@TEXT@", TEXT).replace("@MONO@", MONO)



def _e(value) -> str:
    return html.escape("" if value is None else str(value), quote=True)


class Pages:
    """HTML of the pages, all rendered from the session and its jobs."""

    def __init__(self, session: Session) -> None:
        self.session = session

    def url(self, path: str, **query) -> str:
        query["t"] = self.session.token
        return path + "?" + urllib.parse.urlencode(query)

    def page(self, title: str, body: str, here: str = "", refresh: Optional[int] = None) -> str:
        nav = "".join(
            f'<a href="{self.url(path)}"{" class=here" if key == here else ""}>{label}</a>'
            for key, path, label in (("new", "/", "New run"), ("examples", "/examples", "Examples"),
                                     ("compare", "/compare", "Compare"), ("runs", "/runs", "Runs")))
        meta = f'<meta http-equiv="refresh" content="{refresh}">' if refresh else ""
        return (f'<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">'
                f'<meta name="viewport" content="width=device-width,initial-scale=1">{meta}'
                f"<title>{_e(title)} — AMR-ClonalShare</title><style>{_style()}</style></head><body>"
                f'<div class="sheet"><header><h1>AMR-ClonalShare<span>{_e(__version__)}</span></h1>'
                f"<nav>{nav}</nav></header>{body}"
                f'<div class="colophon">Working folder <code>{_e(self.session.root)}</code>. '
                f"The server answers this computer only; every run keeps its configuration, "
                f"record and reports in a folder of its own. "
                f'<button type="button" class="quiet small" onclick="navigator.clipboard.writeText(location.href)">Copy the address of this page</button>'
                + (f' <a href="{self.url("/manual")}">Manual</a>' if self.session.manual else "")
                + "</div></div></body></html>")

    # -- forms

    def start(self, message: str = "", repeat: str = "") -> str:
        error = f'<div class="error">{_e(message)}</div>' if message else ""
        if repeat:
            error += (f'<div class="status"><span class="kicker">Repeat with new tables</span><p>The column choices, the tested '
                      f'concentrations, the budgets and the seed of run <code>{_e(repeat)}</code> will be proposed for the tables '
                      f'chosen below.</p></div>')
        examples = self.session.example_runs()
        first = (f'<div class="status"><span class="kicker">First run</span><p>To see what a run produces before '
                 f'preparing your own tables, open <a href="{self.url("/examples")}">Examples</a>: the collections shipped '
                 f'with the package run with one click, the small one in seconds.</p></div>'
                 if examples else "")
        templates = ""
        if self.session.examples and (self.session.examples / "workflows" / "data").is_dir():
            links = ", ".join(f'<a href="{self.url("/examples/file", name=name)}">{_e(name)}</a>'
                              for name in ("metadata.csv", "calls.csv", "mic.csv"))
            templates = (f'<p class="lead">Small tables of the right shape, to open in a spreadsheet and fill with your own '
                         f'isolates: {links}.</p>')
        body = f"""{first}
<div class="status"><span class="kicker">One run from two or three tables</span>
<p>Export the tables from your laboratory system or spreadsheet as CSV or Excel workbooks (one sheet each).
Column names are free: the next page lists them and guesses which is which.</p></div>{error}
<table class="shapes"><thead><tr><th>Table</th><th>One row per</th><th>Columns</th></tr></thead><tbody>
<tr><td>Metadata</td><td>isolate</td><td>the isolate identifier, the lineage label (sequence type, cluster, serovar) and, if used, the year, laboratory, country, farm or collection of the isolate</td></tr>
<tr><td>Calls</td><td>isolate and antimicrobial (or one column per antimicrobial)</td><td>the isolate identifier, the antimicrobial and the call: S/I/R, wild-type/non-wild-type or 0/1</td></tr>
<tr><td>MICs</td><td>isolate and antimicrobial</td><td>the isolate identifier, the antimicrobial, the recorded concentration with its sign (<code>&lt;=0.5</code>, <code>8</code>, <code>&gt;16</code>) and, if there are several, the laboratory or panel</td></tr>
</tbody></table>{templates}
<form method="post" action="{self.url('/upload')}" enctype="multipart/form-data">
<input type="hidden" name="from" value="{_e(repeat)}">
<h2><span class="no">1</span>Tables</h2>
<div class="field"><label for="metadata">Metadata table</label><input id="metadata" type="file" name="metadata">
<div class="note">One row per isolate with its lineage label. Leave empty when the call or MIC table carries the lineage column itself: one sheet is then enough.</div></div>
<div class="field"><label for="calls">Call table</label><input id="calls" type="file" name="calls">
<div class="note">Susceptibility calls. Leave empty to analyse the MIC table alone.</div></div>
<div class="field"><label for="mic">MIC table</label><input id="mic" type="file" name="mic">
<div class="note">Recorded minimum inhibitory concentrations, one row per isolate and antimicrobial, or one column per antimicrobial with the sign inside the cell (<code>&lt;=0.5</code>). Leave empty to analyse the calls alone.</div></div>
<div class="field"><label for="name">Name of the run</label><input id="name" type="text" name="name" value="collection" maxlength="60">
<div class="note">Names the folder of the run and the record.</div></div>
<details class="more"><summary>Paste a table instead of choosing a file</summary>
<div class="field"><label for="paste_metadata">Metadata, pasted</label><textarea id="paste_metadata" name="paste_metadata" rows="3" placeholder="Copy the cells from the spreadsheet, header row included, and paste them here."></textarea>
<div class="note">Cells copied from a spreadsheet, with the header row; the pasted table is saved beside the run.</div></div>
<div class="field"><label for="paste_calls">Calls, pasted</label><textarea id="paste_calls" name="paste_calls" rows="3"></textarea></div>
<div class="field"><label for="paste_mic">MICs, pasted</label><textarea id="paste_mic" name="paste_mic" rows="3"></textarea></div>
</details>
<div class="actions"><button type="submit">Read the tables</button></div>
</form>"""
        return self.page("New run", body, "new")

    def mapping(self, folder: Path, files: Dict[str, str], info: Dict[str, Any],
                message: str = "") -> str:
        columns: Dict[str, List[str]] = info["columns"]
        guesses: Dict[str, Any] = info["guesses"]
        error = f'<div class="error">{_e(message)}</div>' if message else ""

        def options(table: str, chosen, none_label: str = "none") -> str:
            items = [f'<option value=""{" selected" if not chosen else ""}>— {none_label} —</option>']
            for column in columns.get(table, []):
                selected = " selected" if column == chosen else ""
                items.append(f'<option value="{_e(column)}"{selected}>{_e(column)}</option>')
            return "".join(items)

        def boxes(key: str, tables: Sequence[str], chosen: Sequence[str]) -> str:
            items = []
            for table in tables:
                for column in columns.get(table, []):
                    checked = " checked" if column in chosen else ""
                    items.append(f'<label class="box"><input type="checkbox" name="{key}" value="{_e(column)}"{checked}> {_e(column)}</label>')
            return '<div class="boxes">' + "".join(items) + "</div>"

        def select(key: str, label: str, table: str, note: str, required: bool = False) -> str:
            req = " required" if required else ""
            return (f'<div class="field"><label for="{key}" title="{_e(note)}">{_e(label)}</label>'
                    f'<select id="{key}" name="{key}"{req}>{options(table, guesses.get(key))}</select>'
                    f'<div class="note">{note}</div></div>')

        def choice(key: str, label: str, pairs, chosen, note: str) -> str:
            items = "".join(f'<option value="{_e(v)}"{" selected" if v == chosen else ""}>{_e(text)}</option>'
                            for v, text in pairs)
            return (f'<div class="field"><label for="{key}" title="{_e(note)}">{_e(label)}</label>'
                    f'<select id="{key}" name="{key}">{items}</select><div class="note">{note}</div></div>')

        def text(key: str, label: str, note: str, value: str = "") -> str:
            return (f'<div class="field"><label for="{key}" title="{_e(note)}">{_e(label)}</label>'
                    f'<input id="{key}" type="text" name="{key}" value="{_e(value)}">'
                    f'<div class="note">{note}</div></div>')

        def number(key: str, label: str, value, note: str, step: str = "1", minimum: str = "0") -> str:
            return (f'<div class="field"><label for="{key}" title="{_e(note)}">{_e(label)}</label>'
                    f'<input id="{key}" type="number" name="{key}" value="{_e(value)}" step="{step}" min="{minimum}">'
                    f'<div class="note">{note}</div></div>')

        def more(summary: str, inner: str) -> str:
            return f'<details class="more"><summary>{summary}</summary>{inner}</details>'

        def preview(table: str) -> str:
            name = files.get(table)
            return _preview(folder / "data" / name) if name else ""

        tables = "".join(f"<code>{_e(name)}</code> ({len(columns.get(table, []))} columns)"
                         + ("; " if i < len(files) - 1 else "")
                         for i, (table, name) in enumerate(files.items()))
        parts = [f"""
<div class="status"><span class="kicker">Columns of the run</span>
<p>Read from {tables}. Each choice below was guessed from the column names; the first rows of every table are
shown so that the guesses can be checked. The settings under "More" are optional and most runs leave them as they are.</p></div>{error}
<form method="post" action="{self.url('/run', run=folder.name)}">
<h2><span class="no">1</span>Isolates and lineages (metadata table)</h2>{preview("metadata")}"""]
        parts.append(select("strain_id_column", "Isolate identifier", "metadata",
                            "The column naming the isolate; the same identifiers appear in every table.", True))
        parts.append(select("lineage_column", "Lineage", "metadata",
                            "The label whose share is measured: sequence type, cluster, serovar.", True))
        parts.append(select("batch_column", "Year or intake", "metadata",
                            "Year or batch of arrival: the evidence across agents is then also followed intake by intake. Optional."))
        inner = select("phenotype_covariate_column", "Strata of the calls", "metadata",
                       "Country or laboratory: lineage labels are compared only within its levels, so a difference between countries is not read as a lineage effect. Optional.")
        inner += select("unit_column", "Sampling unit", "metadata",
                        "Farm or host: the permutation test exchanges lineage labels within units. Optional.")
        inner += select("stratify_by", "Repeat within", "metadata",
                        "Every analysis is repeated on the isolates of each level of this column, beside the whole collection. Optional.")
        inner += select("contrast_column", "Two collections to compare", "metadata",
                        "A column naming two collections, such as two periods: the difference in prevalence between them is decomposed into a change of lineage composition and a change within lineages. Optional.")
        inner += text("contrast_levels", "Their two values",
                      "The two values of that column, separated by a semicolon. Every difference is the first minus the second, so put the later collection first.", str(guesses.get("contrast_levels") or ""))
        parts.append(more("More: strata, sampling units, two collections", inner))
        if "calls" in files:
            parts.append(f'<h2><span class="no">2</span>Calls (call table)</h2>{preview("calls")}')
            parts.append(select("phenotype_id_column", "Isolate identifier", "calls", "The column naming the isolate.", True))
            parts.append(select("phenotype_antibiotic_column", "Antimicrobial", "calls",
                                "Long shape: the column naming the antimicrobial of each row."))
            parts.append(select("phenotype_call_column", "Call", "calls",
                                "Long shape: the column holding the call of each row."))
            chosen = guesses.get("phenotype_agent_columns")
            agents = [str(c) for c in chosen] if isinstance(chosen, (list, tuple)) else []
            parts.append(f'<div class="field"><span class="lab">Agent columns</span>{boxes("phenotype_agent_columns", ["calls"], agents)}'
                         f'<div class="note">Wide shape, one column per antimicrobial: tick the agent columns instead of naming an antimicrobial and a call column.</div></div>')
            parts.append(choice("phenotype_kind", "Kind of call", _KINDS, guesses.get("phenotype_kind") or "undeclared",
                                "What the values record. Declaring it fixes how they are read."))
            inner = choice("phenotype_intermediate", "Intermediate calls", _INTERMEDIATE, str(guesses.get("phenotype_intermediate") or ""),
                           "How an intermediate (I) call is counted; clinical calls only.")
            inner += choice("duplicate_policy", "Conflicting duplicates", _DUPLICATES, str(guesses.get("duplicate_policy") or "error"),
                            "What happens when one isolate has two different calls for one agent.")
            inner += text("phenotype_positive_definition", "Positive outcome",
                          "What a positive call means, written into the record and the report.", str(guesses.get("phenotype_positive_definition") or ""))
            inner += text("phenotype_source", "Interpretation source", "Where the calls were interpreted.", str(guesses.get("phenotype_source") or ""))
            inner += text("ast_standard", "AST standard", "CLSI, EUCAST, ...", str(guesses.get("ast_standard") or ""))
            inner += text("ast_version", "AST version", "", str(guesses.get("ast_version") or ""))
            parts.append(more("More: intermediate calls, duplicates, how the calls were made", inner))
        if "mic" in files:
            no = "3" if "calls" in files else "2"
            parts.append(f'<h2><span class="no">{no}</span>MIC readings (MIC table)</h2>{preview("mic")}')
            parts.append(select("mic_id_column", "Isolate identifier", "mic", "The column naming the isolate.", True))
            parts.append(select("mic_antibiotic_column", "Antimicrobial", "mic", "Long shape: the column naming the antimicrobial of each row."))
            parts.append(select("mic_value_column", "Measurement", "mic", "Long shape: the recorded concentration of each row, in mg/L."))
            parts.append(select("mic_operator_column", "Censoring sign", "mic",
                                "The column holding <, <=, > or >= where the export writes the sign apart from the value. A sign written inside the value is read as it stands. Optional."))
            chosen = guesses.get("mic_agent_columns")
            mic_agents = [str(c) for c in chosen] if isinstance(chosen, (list, tuple)) else []
            parts.append(f'<div class="field"><span class="lab">Agent columns</span>{boxes("mic_agent_columns", ["mic"], mic_agents)}'
                         f'<div class="note">Wide shape, one column per antimicrobial with the sign inside the cell: tick the agent columns instead of naming an antimicrobial and a measurement column.</div></div>')
            unified = info.get("agent_names") or {}
            if unified:
                parts.append('<div class="field"><span class="lab">Agent names</span><div>'
                             + ", ".join(f"<code>{_e(k)}</code> read as {_e(v)}" for k, v in unified.items())
                             + '</div><div class="note">Codes and abbreviations are read as the agents they stand for, so that the call and the MIC table meet; the record lists them.</div></div>')
            parts.append(select("mic_panel_column", "Laboratory or panel", "mic",
                                "Readings are ordered within its levels and the tested wells are found within each; needed when laboratories tested different ranges. Optional."))
            inner = select("mic_unit_column", "Unit", "mic",
                           "The column holding the unit of every reading; rows in a unit other than mg/L are set aside. Optional.")
            inner += (f'<div class="field"><label for="mic_wells">Tested concentrations</label>'
                      f'<textarea id="mic_wells" name="mic_wells" rows="4" placeholder="ciprofloxacin: 0.5, 1, 2&#10;gentamicin: 2, 4, 8">{_e(_wells_to_text(guesses.get("mic_wells")))}</textarea>'
                      f'<div class="note">One agent per line, with the concentrations of its dilution range, lowest to highest well. Without them the range is taken from the readings themselves. Optional.</div></div>')
            inner += _ranges_block(info.get("mic_ranges") or [])
            presets = [("", "none: the ranges above, or the readings themselves")] + [
                (name, preset["title"]) for name, preset in PANEL_PRESETS.items()]
            inner += choice("mic_panel_preset", "Panel preset", presets, guesses.get("mic_panel_preset") or "",
                            "A shipped panel whose published ranges fill the tested concentrations of every agent not written above (Commission Implementing Decision (EU) 2020/1729); the record names it, and the plate sheet of the laboratory remains the authority.")
            covariates = guesses.get("mic_covariate_columns")
            covariates = [str(c) for c in covariates] if isinstance(covariates, (list, tuple)) else []
            inner += (f'<div class="field"><span class="lab">Covariates of the reading</span>{boxes("mic_covariate_columns", ["mic", "metadata"], covariates)}'
                      f'<div class="note">Country or year: readings are ordered within their levels crossed with the panel. Optional.</div></div>')
            parts.append(more("More: unit, tested concentrations, panel preset, covariates", inner))
        no = str(2 + ("calls" in files) + ("mic" in files))
        budgets: Dict[str, Any] = dict(guesses["budgets"]) if isinstance(guesses.get("budgets"), dict) else {}
        inner = number("n_perm", "Permutations", budgets.get("n_perm", 999), "Permuted labellings for the control and the p-value; the smallest p-value is 1/(permutations + 1).", minimum="1")
        inner += number("n_boot", "Bootstrap draws", budgets.get("n_boot", 999), "Draws for the confidence interval; 0 reports the estimate only.")
        inner += number("folds", "Folds", budgets.get("folds", 5), "Folds of the cross-validated score.", minimum="2")
        inner += number("repeats", "Fold assignments", budgets.get("repeats", 20), "Random fold assignments averaged.", minimum="1")
        inner += number("surveillance_n_boot", "Decomposition draws", budgets.get("surveillance_n_boot", 2000), "Bootstrap draws of the decomposition of a prevalence difference.")
        inner += number("alpha", "Evidence level", budgets.get("alpha", 0.05), "Level of the e-value rule (declared at 1/level).", step="any", minimum="0.001")
        inner += choice("end_wells_censored", "End wells", (("on", "censored: a reading on an end well is an interval"), ("off", "exact")), str(budgets.get("end_wells_censored", "on")),
                        "How a reading on the lowest or highest tested well is read.")
        parts.append(f'<h2><span class="no">{no}</span>Seed and budgets</h2>')
        parts.append(number("seed", "Seed", budgets.get("seed", 42), "The run is a function of the data, the settings and the seed; the same seed repeats it."))
        parts.append(more("More: permutations, bootstrap draws, folds (the defaults are the validated settings)", inner))
        parts.append("""
<div class="actions"><button type="submit" name="action" value="run">Run the analysis</button>
<button type="submit" name="action" value="check" class="quiet">Check the input only</button></div>
<p class="note">The input check reads the tables, counts what every analysis would use and set aside, and writes no estimate; it takes seconds.</p>
</form>""")
        return self.page("Columns", "".join(parts), "new")

    def compare_form(self, message: str = "", folder: Optional[Path] = None,
                     file: str = "", columns: Sequence[str] = ()) -> str:
        error = f'<div class="error">{_e(message)}</div>' if message else ""
        if folder is None:
            body = f"""
<div class="status"><span class="kicker">Two lineage definitions on matched records</span>
<p>One table with one row per isolate: the isolate identifier, a 0/1 outcome and two lineage
columns. The comparison splits the difference between the two lineage shares into the records each
definition keeps, the singletons each sets aside, and the relabelling of the isolates both score, each with a
paired confidence interval.</p></div>{error}
<form method="post" action="{self.url('/compare/upload')}" enctype="multipart/form-data">
<h2><span class="no">1</span>Table</h2>
<div class="field"><label for="table">Table</label><input id="table" type="file" name="table" required>
<div class="note">CSV, one row per isolate.</div></div>
<div class="field"><label for="name">Name</label><input id="name" type="text" name="name" value="comparison" maxlength="60"></div>
<div class="actions"><button type="submit">Read the table</button></div></form>"""
            return self.page("Compare", body, "compare")

        def select(key: str, label: str, chosen: Optional[str], note: str) -> str:
            items = "".join(f'<option value="{_e(c)}"{" selected" if c == chosen else ""}>{_e(c)}</option>' for c in columns)
            return (f'<div class="field"><label for="{key}">{_e(label)}</label>'
                    f'<select id="{key}" name="{key}" required>{items}</select><div class="note">{note}</div></div>')

        cols = list(columns)
        guess_id = _match(cols, _PATTERNS["strain_id_column"])
        guess_outcome = _match(cols, ("outcome", "call", "resistant", "positive", "result"))
        lineage_like = [c for c in cols if _match([c], _PATTERNS["lineage_column"])]
        a = lineage_like[0] if lineage_like else None
        b = lineage_like[1] if len(lineage_like) > 1 else None
        body = f"""
<div class="status"><span class="kicker">Columns of the comparison</span>
<p>Read from <code>{_e(file)}</code> ({len(cols)} columns).</p></div>{error}
<form method="post" action="{self.url('/compare/run', run=folder.name)}">
<h2><span class="no">1</span>Columns</h2>
{select('id_column', 'Isolate identifier', guess_id, 'Column naming the isolate.')}
{select('outcome', 'Outcome', guess_outcome, 'Column holding the 0/1 outcome.')}
{select('lineage_a', 'First definition', a, 'The first lineage column, A.')}
{select('lineage_b', 'Second definition', b, 'The second lineage column, B.')}
<h2><span class="no">2</span>Budgets and seed</h2>
<div class="field"><label for="permutations">Permutations</label><input id="permutations" type="number" name="permutations" value="999" min="1"></div>
<div class="field"><label for="bootstraps">Bootstrap draws</label><input id="bootstraps" type="number" name="bootstraps" value="999" min="0"></div>
<div class="field"><label for="folds">Folds</label><input id="folds" type="number" name="folds" value="5" min="2"></div>
<div class="field"><label for="repeats">Fold assignments</label><input id="repeats" type="number" name="repeats" value="20" min="1"></div>
<div class="field"><label for="seed">Seed</label><input id="seed" type="number" name="seed" value="42" min="0"></div>
<div class="actions"><button type="submit">Run the comparison</button></div></form>"""
        return self.page("Compare", body, "compare")

    def examples(self, message: str = "") -> str:
        error = f'<div class="error">{_e(message)}</div>' if message else ""
        items = self.session.example_runs()
        parts = [f"""
<div class="status"><span class="kicker">Collections ready to run</span>
<p>The collections shipped with the package, each with its configuration file. "Run the analysis" copies
the tables and the configuration into a folder of the run and starts it: the fictional collection finishes
in seconds, the real collections in a few minutes. The page of the run shows the results and the report
as they are written, and the same run can be repeated from the command line with
<code>amr-clonalshare --config</code> on the folder's configuration file.</p></div>{error}"""]
        if not items:
            parts.append('<p class="lead">No examples folder was given to this session '
                         '(<code>amr-clonalshare-gui --examples &lt;folder&gt;</code>).</p>')
        for no, item in enumerate(items, 1):
            parts.append(f'<h2><span class="no">{no}</span>{_e(item["title"])}</h2>')
            if item["description"]:
                parts.append(f'<p class="lead">{_e(item["description"])}</p>')
            parts.append(f'<p class="note"><code>{_e(item["id"])}</code> in <code>{_e(self.session.examples)}</code></p>')
            if item["kind"] == "analysis":
                parts.append(f"""<form method="post" action="{self.url('/examples/run')}" class="row">
<input type="hidden" name="example" value="{_e(item['id'])}">
<button type="submit" name="action" value="run">Run the analysis</button>
<button type="submit" name="action" value="check" class="quiet">Check the input only</button>
<label class="inline" for="seed_{no}">Seed</label><input id="seed_{no}" type="number" name="seed" value="{item.get('seed', 42)}" min="0"></form>""")
            else:
                parts.append(f"""<form method="post" action="{self.url('/examples/compare')}" class="row">
<input type="hidden" name="example" value="{_e(item['id'])}">
<button type="submit">Read the table for the comparison</button></form>""")
        parts.append(f"""<h2><span class="no">{len(items) + 1}</span>A configuration file of your own</h2>
<p class="lead">The run the command line makes from a configuration file already on this computer:
the tables stay where the file names them, and the run gets a folder of its own for its configuration and results.</p>
<form method="post" action="{self.url('/config/run')}">
<div class="field"><label for="config_path">Configuration file</label>
<input id="config_path" type="text" name="config_path" placeholder="C:\\data\\collection\\config.yaml" size="60">
<div class="note">The full path of a <code>config.yaml</code>; its <code>data_dir</code> is read relative to the file, as the command line reads it.</div></div>
<div class="field"><label for="config_seed">Seed</label><input id="config_seed" type="number" name="seed" value="42" min="0"></div>
<div class="actions"><button type="submit" name="action" value="run">Run the analysis</button>
<button type="submit" name="action" value="check" class="quiet">Check the input only</button></div></form>""")
        return self.page("Examples", "".join(parts), "examples")

    # -- jobs

    def job(self, job: Job) -> str:
        running = job.status in ("queued", "running")
        kicker = {"queued": "Queued", "running": "Running", "done": "Completed", "failed": "Not completed"}[job.status]
        flag = " flag" if job.status == "failed" else ""
        elapsed = f"{job.seconds:.0f} s"
        head = (f'<div class="status"><span class="kicker{flag}">{kicker}</span>'
                f"<p>{_e(job.name)} — {_e(job.kind)}, {elapsed}. Folder <code>{_e(job.folder)}</code>.</p></div>")
        if job.error:
            head += (f'<div class="error">{_e(job.error)}<br><span class="note">The log below names the step that stopped; '
                     f'an input error is usually a column named wrongly or a value the tables cannot hold, and "Check the input only" '
                     f'lists what every table was read as.</span></div>')
        with job.lock:
            lines = list(job.lines[-40:])
        log = "<pre>" + _e("\n".join(lines) or "starting") + "</pre>"
        body = head
        no = 1
        if job.status == "done":
            body += self._results(job)
            no = 3 if job.kind == "analysis" and (job.results / "report.html").is_file() else 2
        body += f'<h2><span class="no">{no}</span>Log</h2>' + log
        if running:
            body += '<p class="lead">This page refreshes itself every three seconds; the results appear here when the run ends.</p>'
        return self.page(job.name, body, "runs", refresh=3 if running else None)

    def _results(self, job: Job) -> str:
        parts = ['<h2><span class="no">1</span>Results</h2>']
        record_path = job.results / "clonal_share_result.json"
        check_path = job.results / "input_qc.md"
        if job.kind == "analysis" and record_path.is_file():
            try:
                record = json.loads(record_path.read_text(encoding="utf-8"))
                parts.append(_headline(record))
                parts.append(_share_table(record))
                summary = _summary_paragraph(record, job.name)
                parts.append(f'<details class="more"><summary>A paragraph for a laboratory report</summary>'
                             f'<p class="lead" id="summary">{_e(summary)}</p>'
                             f'<div class="actions"><button type="button" class="quiet" data-text="{_e(summary)}" '
                             f'onclick="navigator.clipboard.writeText(this.dataset.text)">Copy the paragraph</button></div></details>')
            except (OSError, ValueError, KeyError):
                pass
        elif job.kind == "input check" and check_path.is_file():
            parts.append("<pre>" + _e(check_path.read_text(encoding="utf-8", errors="replace")) + "</pre>")
        parts.append('<div class="files">')
        labels = {"report.html": "Report", "report.md": "Report as Markdown", "results.csv": "Results table",
                  "strata_results.csv": "Results within strata", "clonal_share_result.json": "Record",
                  "input_qc.md": "Input check", "input_qc.json": "Input check as JSON",
                  "run_manifest.json": "Manifest of the files", "comparison.json": "Comparison record",
                  "arms.csv": "The six analyses"}
        order = list(labels)
        files = [p for p in job.results.iterdir() if p.is_file()] if job.results.is_dir() else []
        for path in sorted(files, key=lambda p: (order.index(p.name) if p.name in order else len(order), p.name)):
            size = path.stat().st_size
            parts.append(f'<a href="{self.url(f"/job/{job.id}/file/{path.name}")}">{_e(labels.get(path.name, path.name))}'
                         f"<small>{_e(path.name)}, {size:,} bytes</small></a>")
        parts.append(f'<a href="{self.url(f"/job/{job.id}/zip")}">Everything, as one zip archive</a>')
        parts.append("</div>")
        parts.append(f'<div class="actions"><a class="button" href="{self.url(f"/job/{job.id}/folder")}">Open the folder of this run</a>'
                     + (f'<a class="button" href="{self.url("/", **{"from": job.folder.name})}">Repeat with new tables</a>'
                        if job.kind == "analysis" and (job.folder / "config.yaml").is_file() else "")
                     + '</div>')
        if job.kind == "analysis":
            parts.append(f'<p class="lead">The same run from the command line: '
                         f'<code>amr-clonalshare --config "{_e(job.folder / "config.yaml")}" --results-dir &lt;folder&gt;'
                         + (f' --seed {job.seed}' if job.seed is not None else '') + '</code></p>')
            report = job.results / "report.html"
            if report.is_file():
                src = self.url(f"/job/{job.id}/file/report.html")
                parts.append(f'<h2><span class="no">2</span>Report</h2>'
                             f'<p class="lead">The report of the run, with its figures; <a href="{src}" target="_blank">open it in a tab of its own</a> to print or save it.</p>'
                             f'<iframe class="report" src="{src}" title="Report"></iframe>')
        return "".join(parts)

    def runs(self) -> str:
        items = self.session.previous_runs()
        rows = []
        for item in items:
            started = _dt.datetime.fromtimestamp(float(item.get("started") or 0)).strftime("%Y-%m-%d %H:%M")
            live = self.session.jobs.get(str(item.get("id")))
            link = (f'<a href="{self.url(f"/job/{live.id}")}">{_e(item.get("name"))}</a>' if live
                    else f'<a href="{self.url("/open", folder=item["folder"].name)}">{_e(item.get("name"))}</a>')
            remove = (f'<form method="post" action="{self.url("/runs/remove")}" class="inline">'
                      f'<input type="hidden" name="folder" value="{_e(item["folder"].name)}">'
                      f'<button type="submit" class="quiet small">Set aside</button></form>')
            rows.append(f"<tr><td>{started}</td><td>{link}</td><td>{_e(item.get('kind'))}</td>"
                        f"<td>{_e(item.get('status'))}</td><td><code>{_e(item['folder'].name)}</code></td><td>{remove}</td></tr>")
        for job in self.session.jobs.values():
            if job.status in ("queued", "running"):
                started = _dt.datetime.fromtimestamp(job.started).strftime("%Y-%m-%d %H:%M")
                rows.insert(0, f'<tr><td>{started}</td><td><a href="{self.url(f"/job/{job.id}")}">{_e(job.name)}</a></td>'
                               f"<td>{_e(job.kind)}</td><td>{_e(job.status)}</td><td><code>{_e(job.folder.name)}</code></td><td></td></tr>")
        table = ("<table><thead><tr><th>Started</th><th>Name</th><th>Kind</th><th>Status</th><th>Folder</th><th></th></tr></thead>"
                 f"<tbody>{''.join(rows)}</tbody></table>" if rows else '<p class="lead">No run yet.</p>')
        body = (f'<div class="status"><span class="kicker">Runs</span><p>Every run under '
                f"<code>{_e(self.session.runs)}</code>. \"Set aside\" moves a run into the <code>_removed</code> folder "
                f"beside the others; nothing is deleted.</p></div>{table}")
        return self.page("Runs", body, "runs")


def _headline(record: dict) -> str:
    """One line on what the run analysed."""
    diagnostics = record.get("metadata_diagnostics") or {}
    calls = diagnostics.get("clonal_share") or {}
    orders = ((diagnostics.get("censored_share") or {}).get("per_agent") or {})
    blocks = [b for b in calls.values() if isinstance(b, dict)]
    blocks += [(o.get("order") or {}) for o in orders.values() if isinstance(o, dict)]
    groups = max((b.get("n_groups") or 0 for b in blocks), default=0)
    singletons = max((b.get("n_singletons_set_aside") or 0 for b in blocks), default=0)
    agents = sorted(set(calls) | set(orders))
    computed = sum(1 for a in agents if (calls.get(a) or {}).get("estimable") or ((orders.get(a) or {}).get("order") or {}).get("estimable"))
    lineage = diagnostics.get("lineage_column") or "lineage"
    n_agents = len(agents)
    text = (f"{record.get('n_isolates', '?')} isolates, {groups} lineages ({_e(lineage)}), "
            f"{n_agents} antimicrobial{'s' if n_agents != 1 else ''}; a lineage share was computed for "
            f"{'all' if computed == n_agents and n_agents > 1 else computed}{' of them' if computed != n_agents or n_agents == 1 else ''}.")
    if singletons:
        text += f" {singletons} isolates in lineages seen once were set aside and counted."
    return f'<p class="lead">{text}</p>'


def _summary_paragraph(record: dict, name: str) -> str:
    """One paragraph a laboratory report can carry: the collection, the
    agents whose lineage structure is established, the three numbers of the
    strongest agent, the seed and the version."""
    diagnostics = record.get("metadata_diagnostics") or {}
    rows = summary_rows(record)
    blocks = [b for b in (diagnostics.get("clonal_share") or {}).values() if isinstance(b, dict)]
    blocks += [(o.get("order") or {}) for o in ((diagnostics.get("censored_share") or {}).get("per_agent") or {}).values()]
    groups = max((b.get("n_groups") or 0 for b in blocks), default=0)
    lineage = diagnostics.get("lineage_column") or "lineage"
    established = [r["agent"] for r in rows if any(
        (r.get(k) or {}).get("label") == "lineage structure established" for k in ("call_conclusion", "mic_conclusion"))]
    text = (f"{name}: {record.get('n_isolates', '?')} isolates in {groups} lineages ({lineage}), "
            f"{len(rows)} antimicrobial{'s' if len(rows) != 1 else ''} analysed with AMR-ClonalShare {__version__} "
            f"(seed {record.get('seed')}). ")
    if established:
        text += ("Lineage structure was established for " + ", ".join(established)
                 + f" ({len(established)} of {len(rows)}). ")
    else:
        text += "Lineage structure was established for no antimicrobial. "
    strongest = None
    for r in rows:
        for key in ("mic_share", "call_share"):
            try:
                value = float(r.get(key, "—"))
            except ValueError:
                continue
            if strongest is None or value > strongest[1]:
                strongest = (r, value, key)
    if strongest:
        r, value, key = strongest
        if key == "mic_share":
            text += (f"The largest lineage share of the MIC ordering was {r['mic_share']} for {r['agent']} "
                     f"(95 % interval {r['mic_interval']}, p {r['mic_p']}; lower bound on the latent ordering "
                     f"{r['lower_bound']}, lower confidence limit {r['lower_limit']}). ")
        else:
            text += (f"The largest lineage share of a call was {r['call_share']} for {r['agent']} "
                     f"(95 % interval {r['call_interval']}, p {r['call_p']}). ")
    text += "The shares describe the represented lineages of this collection and say nothing about transmission."
    return text


def _share_table(record: dict) -> str:
    """The lineage shares of the record, agent by agent, each with its conclusion."""
    rows = summary_rows(record)
    if not rows:
        return ""
    with_calls = any("call_share" in r for r in rows)
    with_orders = any("mic_share" in r for r in rows)
    head = ["<th>Antimicrobial</th>"]
    if with_calls:
        head += ['<th class="num">Lineage share of the call</th>', '<th class="num">95% CI</th>', '<th class="num">p</th>',
                 '<th class="num">q</th>', '<th class="num">Lower bound</th>', '<th class="num">Lower confidence limit</th>',
                 '<th>Conclusion for the call</th>']
    if with_orders:
        head += ['<th class="num">Lineage share of the MIC ordering</th>', '<th class="num">95% CI</th>', '<th class="num">p</th>',
                 '<th class="num">q</th>', '<th class="num">Lower bound</th>', '<th class="num">Lower confidence limit</th>',
                 '<th>Conclusion for the MICs</th>']

    def conclusion(cell) -> str:
        if not cell:
            return "<td></td>"
        return f'<td title="{_e(cell["text"])}">{_e(short_label(cell["label"]))}</td>'

    body = []
    for r in rows:
        cells = [f"<td>{_e(r['agent'])}</td>"]
        if with_calls:
            cells += [f'<td class="num">{_e(r.get("call_share", "—"))}</td>', f'<td class="num">{_e(r.get("call_interval", "—"))}</td>',
                      f'<td class="num">{_e(r.get("call_p", "—"))}</td>', f'<td class="num">{_e(r.get("call_q", "—"))}</td>',
                      f'<td class="num">{_e(r.get("call_lower_bound", "—"))}</td>',
                      f'<td class="num">{_e(r.get("call_lower_limit", "—"))}</td>', conclusion(r.get("call_conclusion"))]
        if with_orders:
            cells += [f'<td class="num">{_e(r.get("mic_share", "—"))}</td>', f'<td class="num">{_e(r.get("mic_interval", "—"))}</td>',
                      f'<td class="num">{_e(r.get("mic_p", "—"))}</td>', f'<td class="num">{_e(r.get("mic_q", "—"))}</td>',
                      f'<td class="num">{_e(r.get("lower_bound", "—"))}</td>',
                      f'<td class="num">{_e(r.get("lower_limit", "—"))}</td>', conclusion(r.get("mic_conclusion"))]
        body.append("<tr>" + "".join(cells) + "</tr>")
    sentences = "".join(f"<li><b>{_e(r['agent'])}</b>: " + " ".join(
        _e(r[k]["text"]) for k in ("call_conclusion", "mic_conclusion") if r.get(k)) + "</li>" for r in rows)
    return (f"<div class=\"wide\"><table><thead><tr>{''.join(head)}</tr></thead><tbody>{''.join(body)}</tbody></table></div>"
            '<p class="note">The lineage share is the part of the variation of the call, or of the recorded MIC ordering, '
            'that lies between lineages, with its 95% confidence interval, the p-value of the permutation test '
            '(≤ marks the smallest value the permutations can give) and the q-value, the p-value adjusted for the number '
            'of antimicrobials tested (Benjamini-Yekutieli). '
            'The lower bound is the smallest lineage share of the latent MIC ordering compatible with the readings, '
            'and the lower confidence limit extends it to further isolates of the same lineages. '
            f'A dash marks a result the collection could not give; the report states the reason. The conclusion of a row is '
            f'established (the interval lies above zero and q is at most 0.05, or the lower confidence limit is above zero), '
            f'not established, whole range (the interval spans 0 to 1) or not estimable; the sentence behind each conclusion '
            f'is folded below. {_e(VERDICT_NOTE)}</p>'
            f'<details class="more"><summary>The conclusion for every agent in full</summary><ul class="readings">{sentences}</ul></details>')


# ---------------------------------------------------------------- requests


def _parse_multipart(content_type: str, body: bytes) -> Tuple[Dict[str, List[str]], Dict[str, Tuple[str, bytes]]]:
    """Fields and files of a multipart/form-data body, by field name."""
    message = email.message_from_bytes(
        b"Content-Type: " + content_type.encode("ascii", "replace") + b"\r\nMIME-Version: 1.0\r\n\r\n" + body,
        policy=email.policy.HTTP)
    fields: Dict[str, List[str]] = {}
    files: Dict[str, Tuple[str, bytes]] = {}
    for part in message.iter_parts():
        name = part.get_param("name", header="content-disposition")
        if not name:
            continue
        payload = part.get_payload(decode=True)
        data = payload if isinstance(payload, bytes) else b""
        filename = part.get_filename()
        if filename:
            if data:
                files[str(name)] = (Path(filename).name, data)
        else:
            fields.setdefault(str(name), []).append(data.decode("utf-8", "replace"))
    return fields, files


class _LoopbackServer(ThreadingHTTPServer):
    """The server on 127.0.0.1. ``HTTPServer.server_bind`` looks up the fully
    qualified name of the bound address, a reverse DNS query that can take
    many seconds on a machine without a resolver for the loopback address
    (macOS runners), and the name is never used: the address printed and
    opened is the numeric one."""
    daemon_threads = True

    def server_bind(self) -> None:
        socketserver.TCPServer.server_bind(self)
        self.server_name = str(self.server_address[0])
        self.server_port = self.server_address[1]


def _serve(session: Session, port: int, open_browser: bool, ready: Optional[Callable[[str], None]] = None) -> None:
    pages = Pages(session)

    class Handler(BaseHTTPRequestHandler):
        server_version = f"amr-clonalshare/{__version__}"

        def log_message(self, *_):  # the console is the launcher's, not a request log
            return None

        # -- plumbing

        def _send(self, body: str, status: HTTPStatus = HTTPStatus.OK, content_type: str = "text/html; charset=utf-8") -> None:
            data = body.encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.end_headers()
            self.wfile.write(data)

        def _send_file(self, path: Path, download: Optional[str] = None) -> None:
            kind, _ = mimetypes.guess_type(path.name)
            if path.suffix in (".md", ".csv", ".json", ".yaml", ".log", ".txt"):
                kind = "text/plain; charset=utf-8"
            data = path.read_bytes()
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", kind or "application/octet-stream")
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "no-store")
            if download:
                self.send_header("Content-Disposition", f'attachment; filename="{download}"')
            self.end_headers()
            self.wfile.write(data)

        def _authorised(self, query: Dict[str, List[str]]) -> bool:
            return secrets.compare_digest(query.get("t", [""])[0], session.token)

        def _refuse(self) -> None:
            self._send(pages.page("Not this address", '<div class="status"><span class="kicker flag">Not this address</span>'
                                  "<p>Open the address the launcher printed, which carries the token of this session.</p></div>"),
                       HTTPStatus.FORBIDDEN)

        def _body(self) -> Tuple[Dict[str, List[str]], Dict[str, Tuple[str, bytes]]]:
            length = int(self.headers.get("Content-Length") or 0)
            content_type = self.headers.get("Content-Type") or ""
            body = self.rfile.read(length) if length else b""
            if content_type.startswith("multipart/form-data"):
                return _parse_multipart(content_type, body)
            return urllib.parse.parse_qs(body.decode("utf-8", "replace"), keep_blank_values=True), {}

        # -- routes

        def do_GET(self) -> None:
            url = urllib.parse.urlsplit(self.path)
            query = urllib.parse.parse_qs(url.query)
            if not self._authorised(query):
                return self._refuse()
            path = url.path
            if path == "/":
                return self._send(pages.start(repeat=query.get("from", [""])[0]))
            if path == "/compare":
                return self._send(pages.compare_form())
            if path == "/examples":
                return self._send(pages.examples())
            if path == "/runs":
                return self._send(pages.runs())
            if path == "/open":
                return self._open(query.get("folder", [""])[0])
            if path == "/manual" or path.startswith("/manual/"):
                return self._manual(path[len("/manual"):].lstrip("/"))
            if path == "/examples/file":
                return self._example_file(query.get("name", [""])[0])
            match = re.fullmatch(r"/job/([0-9a-f]+)(?:/(file/([^/]+)|zip|folder))?", path)
            if match:
                job = session.jobs.get(match.group(1))
                if job is None:
                    return self._send(pages.page("Unknown run", '<div class="error">No such run in this session.</div>'), HTTPStatus.NOT_FOUND)
                if match.group(2) is None:
                    return self._send(pages.job(job))
                if match.group(2) == "zip":
                    return self._zip(job)
                if match.group(2) == "folder":
                    return self._folder(job)
                name = urllib.parse.unquote(match.group(3))
                target = (job.results / name)
                if not _safe_name(name) or not target.is_file() or target.parent != job.results:
                    return self._send(pages.page("Unknown file", '<div class="error">No such file.</div>'), HTTPStatus.NOT_FOUND)
                return self._send_file(target)
            return self._send(pages.page("Unknown page", '<div class="error">No such page.</div>'), HTTPStatus.NOT_FOUND)

        def do_POST(self) -> None:
            url = urllib.parse.urlsplit(self.path)
            query = urllib.parse.parse_qs(url.query)
            if not self._authorised(query):
                return self._refuse()
            try:
                fields, files = self._body()
            except (ValueError, UnicodeDecodeError) as error:
                return self._send(pages.start(f"The form could not be read: {error}"), HTTPStatus.BAD_REQUEST)
            if url.path == "/upload":
                return self._upload(fields, files)
            if url.path == "/run":
                return self._run(query.get("run", [""])[0], fields)
            if url.path == "/compare/upload":
                return self._compare_upload(fields, files)
            if url.path == "/compare/run":
                return self._compare_run(query.get("run", [""])[0], fields)
            if url.path == "/examples/run":
                return self._example_run(fields)
            if url.path == "/examples/compare":
                return self._example_compare(fields)
            if url.path == "/config/run":
                return self._config_run(fields)
            if url.path == "/runs/remove":
                return self._remove(fields)
            return self._send(pages.page("Unknown page", '<div class="error">No such page.</div>'), HTTPStatus.NOT_FOUND)

        # -- the analysis

        def _upload(self, fields, files) -> None:
            files = dict(files)
            for table in ("metadata", "calls", "mic"):
                pasted = (fields.get(f"paste_{table}") or [""])[0]
                if table not in files and pasted.strip():
                    # cells copied from a spreadsheet arrive tab-separated
                    files[table] = (f"pasted_{table}.csv", pasted.strip().encode("utf-8") + b"\n")
            if "calls" not in files and "mic" not in files:
                return self._send(pages.start("Give a call table, a MIC table, or both."), HTTPStatus.BAD_REQUEST)
            name = (fields.get("name") or ["collection"])[0].strip() or "collection"
            folder = session.new_folder("analysis", name)
            saved: Dict[str, str] = {}
            paths: Dict[str, Path] = {}
            for table in ("metadata", "calls", "mic"):
                if table in files:
                    filename, payload = files[table]
                    target = folder / "data" / filename
                    if target.exists():  # two tables uploaded under one file name
                        target = folder / "data" / f"{table}_{filename}"
                    target.write_bytes(payload)
                    saved[table] = target.name
                    paths[table] = target
            if "metadata" not in paths:
                # One sheet: the table that carries the lineage column is the
                # metadata too.
                try:
                    carrier = next((t for t in ("calls", "mic")
                                    if t in paths and _match(_columns(paths[t]), _PATTERNS["lineage_column"])), None)
                except (OSError, ValueError, UnicodeDecodeError) as error:
                    return self._send(pages.start(f"A table could not be read: {error}"), HTTPStatus.BAD_REQUEST)
                if carrier is None:
                    return self._send(pages.start("Give a metadata table, or a call or MIC table that carries the "
                                                  "lineage column (sequence type, cluster, serovar) itself."),
                                      HTTPStatus.BAD_REQUEST)
                saved["metadata"] = saved[carrier]
                paths["metadata"] = paths[carrier]
            try:
                info = _table_guesses(paths)
            except (OSError, ValueError, UnicodeDecodeError) as error:
                return self._send(pages.start(f"A table could not be read: {error}"), HTTPStatus.BAD_REQUEST)
            earlier = (fields.get("from") or [""])[0]
            if earlier and _safe_name(earlier) and (session.runs / earlier / "config.yaml").is_file():
                _carry_settings(session.runs / earlier / "config.yaml", info)
            (folder / "tables.json").write_text(json.dumps({"name": name, "files": saved}, indent=1), encoding="utf-8")
            (folder / "columns.json").write_text(json.dumps(info["columns"], indent=1), encoding="utf-8")
            (folder / "ranges.json").write_text(json.dumps({"mic_ranges": info.get("mic_ranges") or [],
                                                            "agent_names": info.get("agent_names") or {}},
                                                           indent=1, default=str), encoding="utf-8")
            return self._send(pages.mapping(folder, saved, info))

        def _run(self, run: str, fields) -> None:
            folder = session.runs / run
            tables = folder / "tables.json"
            if not _safe_name(run) or not tables.is_file():
                return self._send(pages.start("The tables of this run were not found; upload them again."), HTTPStatus.NOT_FOUND)
            if session.live_job(folder):
                return self._send(pages.start("This run is still in progress; its page shows when it ends."), HTTPStatus.CONFLICT)
            meta = json.loads(tables.read_text(encoding="utf-8"))
            files: Dict[str, str] = meta["files"]
            try:
                seed = _seed(fields)
                raw = _config_from_form(fields, files)
                (folder / "config.yaml").write_text(
                    "# Written by amr-clonalshare-gui from the form; the same run from the command line:\n"
                    f"#   amr-clonalshare --config config.yaml --results-dir <folder> --seed {seed}\n"
                    + yaml.safe_dump(raw, sort_keys=False, allow_unicode=True), encoding="utf-8")
                load_config(folder / "config.yaml")
            except ConfigError as error:
                return self._back_to_form(folder, files, fields, f"The configuration was refused: {error}")
            check_only = (fields.get("action") or ["run"])[0] == "check"
            if not check_only:
                # The tables are read once before the run starts, so that a
                # column named wrongly comes back to this page with its
                # reason instead of ending the run on its page.
                try:
                    from .io import load_dataset
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        load_dataset(load_config(folder / "config.yaml"))
                except ConfigError as error:
                    return self._back_to_form(folder, files, fields, f"The tables could not be read as named: {error}")
                except (OSError, ValueError) as error:
                    return self._back_to_form(folder, files, fields, f"The tables could not be read: {error}")
            job = Job("input check" if check_only else "analysis", meta["name"], folder)
            job.seed = seed
            session.start(job, _analysis_target(session, seed, check_only))
            self.send_response(HTTPStatus.SEE_OTHER)
            self.send_header("Location", pages.url(f"/job/{job.id}"))
            self.send_header("Content-Length", "0")
            self.end_headers()

        def _back_to_form(self, folder: Path, files: Dict[str, str], fields, message: str) -> None:
            columns = json.loads((folder / "columns.json").read_text(encoding="utf-8"))
            info = {"columns": columns, "guesses": {k: (v[0] if v else None) for k, v in fields.items()}}
            info["guesses"]["phenotype_agent_columns"] = fields.get("phenotype_agent_columns", [])
            info["guesses"]["mic_agent_columns"] = fields.get("mic_agent_columns", [])
            info["guesses"]["mic_covariate_columns"] = fields.get("mic_covariate_columns", [])
            # the settings the reader entered come back as entered
            info["guesses"]["budgets"] = {k: fields[k][0] for k in (
                "n_perm", "n_boot", "folds", "repeats", "surveillance_n_boot", "alpha",
                "end_wells_censored", "seed") if fields.get(k) and fields[k][0].strip()}
            extra = folder / "ranges.json"
            if extra.is_file():
                info.update(json.loads(extra.read_text(encoding="utf-8")))
            return self._send(pages.mapping(folder, files, info, message), HTTPStatus.BAD_REQUEST)

        # -- the examples, and a configuration file of the reader's own

        def _start_analysis(self, folder: Path, name: str, fields, seed: int) -> None:
            check_only = (fields.get("action") or ["run"])[0] == "check"
            job = Job("input check" if check_only else "analysis", name, folder)
            job.seed = seed
            session.start(job, _analysis_target(session, seed, check_only))
            self.send_response(HTTPStatus.SEE_OTHER)
            self.send_header("Location", pages.url(f"/job/{job.id}"))
            self.send_header("Content-Length", "0")
            self.end_headers()

        def _example_run(self, fields) -> None:
            item = session.example((fields.get("example") or [""])[0])
            if item is None or item["kind"] != "analysis":
                return self._send(pages.examples("No such example."), HTTPStatus.NOT_FOUND)
            name = f"{item['path'].parent.name}_{item['path'].stem}"
            try:
                seed = _seed(fields)
                folder = _run_folder_from_config(session, item["path"], name, copy_data=True, seed=seed)
            except (ConfigError, OSError, yaml.YAMLError) as error:
                return self._send(pages.examples(f"The example could not be prepared: {error}"), HTTPStatus.BAD_REQUEST)
            return self._start_analysis(folder, name, fields, seed)

        def _example_compare(self, fields) -> None:
            item = session.example((fields.get("example") or [""])[0])
            if item is None or item["kind"] != "comparison":
                return self._send(pages.examples("No such example."), HTTPStatus.NOT_FOUND)
            name = f"{item['path'].parent.name}_comparison"
            folder = session.new_folder("comparison", name)
            target = folder / "data" / item["path"].name
            target.write_bytes(item["path"].read_bytes())
            try:
                columns = _columns(target)
            except (OSError, ValueError, UnicodeDecodeError) as error:
                return self._send(pages.examples(f"The table could not be read: {error}"), HTTPStatus.BAD_REQUEST)
            (folder / "tables.json").write_text(json.dumps({"name": name, "files": {"table": target.name}}, indent=1),
                                                encoding="utf-8")
            return self._send(pages.compare_form("", folder, target.name, columns))

        def _config_run(self, fields) -> None:
            given = (fields.get("config_path") or [""])[0].strip().strip('"')
            source = Path(given).expanduser()
            if not given or not source.is_file():
                return self._send(pages.examples(f"No configuration file at {given or '(empty)'}."), HTTPStatus.BAD_REQUEST)
            name = source.parent.name or source.stem
            try:
                seed = _seed(fields)
                folder = _run_folder_from_config(session, source.resolve(), name, copy_data=False, seed=seed)
            except (ConfigError, OSError, yaml.YAMLError) as error:
                return self._send(pages.examples(f"The configuration was refused: {error}"), HTTPStatus.BAD_REQUEST)
            return self._start_analysis(folder, name, fields, seed)

        # -- the comparison

        def _compare_upload(self, fields, files) -> None:
            if "table" not in files:
                return self._send(pages.compare_form("The table is required."), HTTPStatus.BAD_REQUEST)
            name = (fields.get("name") or ["comparison"])[0].strip() or "comparison"
            folder = session.new_folder("comparison", name)
            filename, payload = files["table"]
            target = folder / "data" / filename
            target.write_bytes(payload)
            try:
                columns = _columns(target)
            except (OSError, ValueError, UnicodeDecodeError) as error:
                return self._send(pages.compare_form(f"The table could not be read: {error}"), HTTPStatus.BAD_REQUEST)
            (folder / "tables.json").write_text(json.dumps({"name": name, "files": {"table": target.name}}, indent=1),
                                                encoding="utf-8")
            return self._send(pages.compare_form("", folder, target.name, columns))

        def _compare_run(self, run: str, fields) -> None:
            folder = session.runs / run
            tables = folder / "tables.json"
            if not _safe_name(run) or not tables.is_file():
                return self._send(pages.compare_form("The table of this comparison was not found; upload it again."),
                                  HTTPStatus.NOT_FOUND)
            if session.live_job(folder):
                return self._send(pages.compare_form("This comparison is still in progress; its page shows when it ends."),
                                  HTTPStatus.CONFLICT)
            meta = json.loads(tables.read_text(encoding="utf-8"))
            table = folder / "data" / meta["files"]["table"]
            arguments = ["--input", str(table)]
            for key, option in (("id_column", "--id-column"), ("outcome", "--outcome"),
                                ("lineage_a", "--lineage-a"), ("lineage_b", "--lineage-b")):
                value = (fields.get(key) or [""])[0].strip()
                if not value:
                    return self._send(pages.compare_form(f"{key} is required.", folder, table.name, _columns(table)),
                                      HTTPStatus.BAD_REQUEST)
                arguments += [option, value]
            for key, option in (("seed", "--seed"), ("folds", "--folds"), ("repeats", "--repeats"),
                                ("permutations", "--permutations"), ("bootstraps", "--bootstraps")):
                value = (fields.get(key) or [""])[0].strip()
                if value:
                    arguments += [option, value]
            job = Job("comparison", meta["name"], folder)
            # a second comparison of the same table keeps the first one's files
            number = 2
            while job.results.exists() and any(job.results.iterdir()):
                job.results = folder / f"results_{number}"
                number += 1
            session.start(job, _comparison_target(session, arguments))
            self.send_response(HTTPStatus.SEE_OTHER)
            self.send_header("Location", pages.url(f"/job/{job.id}"))
            self.send_header("Content-Length", "0")
            self.end_headers()

        # -- files of a run

        def _zip(self, job: Job) -> None:
            buffer = io.BytesIO()
            with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
                for path in sorted(job.folder.rglob("*")):
                    if path.is_file():
                        archive.write(path, path.relative_to(job.folder).as_posix())
            data = buffer.getvalue()
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", "application/zip")
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Content-Disposition", f'attachment; filename="{job.folder.name}.zip"')
            self.end_headers()
            self.wfile.write(data)

        def _example_file(self, name: str) -> None:
            """One of the three tables of the fictional collection, as a template."""
            base = (session.examples / "workflows" / "data") if session.examples else None
            target = (base / name) if base and name else None
            if target is None or not _safe_name(name) or not target.is_file() or target.parent != base:
                return self._send(pages.page("Unknown file", '<div class="error">No such file.</div>'), HTTPStatus.NOT_FOUND)
            return self._send_file(target, download=name)

        def _folder(self, job: Job) -> None:
            """Open the folder of the run in the file manager of this
            computer, and come back to the page of the run."""
            opened = _open_folder(job.folder)
            if opened:
                self.send_response(HTTPStatus.SEE_OTHER)
                self.send_header("Location", pages.url(f"/job/{job.id}"))
                self.send_header("Content-Length", "0")
                self.end_headers()
                return None
            return self._send(pages.page("Folder", f'<div class="status"><span class="kicker">Folder of the run</span>'
                                         f'<p>The file manager could not be opened from here; the folder is '
                                         f'<code>{_e(job.folder)}</code>.</p></div>'))

        def _manual(self, name: str) -> None:
            """A page of the manual shipped beside the examples, when there is one."""
            folder = session.manual
            if folder is None:
                return self._send(pages.page("Manual", '<div class="error">No manual folder was given to this session; '
                                             'the manual is at https://maciejkochanowski.github.io/amr-clonalshare/.</div>'),
                                  HTTPStatus.NOT_FOUND)
            target = (folder / (name or "index.html")).resolve()
            if target.is_dir():
                target = target / "index.html"
            if not target.is_file() or folder.resolve() not in target.parents:
                return self._send(pages.page("Manual", '<div class="error">No such page of the manual.</div>'), HTTPStatus.NOT_FOUND)
            return self._send_file(target)

        def _remove(self, fields) -> None:
            """Set a run aside: its folder is moved under _removed, and no
            file is deleted."""
            name = (fields.get("folder") or [""])[0]
            folder = session.runs / name
            if not _safe_name(name) or not folder.is_dir():
                return self._send(pages.runs(), HTTPStatus.NOT_FOUND)
            removed = session.runs / "_removed"
            removed.mkdir(exist_ok=True)
            target = removed / name
            if target.exists():
                target = removed / f"{name}_{secrets.token_hex(3)}"
            try:
                folder.rename(target)
            except OSError as error:
                return self._send(pages.page("Runs", f'<div class="error">The run could not be set aside: {error}</div>'),
                                  HTTPStatus.CONFLICT)
            for job_id, job in list(session.jobs.items()):
                if job.folder == folder:
                    del session.jobs[job_id]
            self.send_response(HTTPStatus.SEE_OTHER)
            self.send_header("Location", pages.url("/runs"))
            self.send_header("Content-Length", "0")
            self.end_headers()

        def _open(self, name: str) -> None:
            """A run of an earlier session, re-listed from its folder."""
            folder = session.runs / name
            record = folder / "job.json"
            if not _safe_name(name) or not record.is_file():
                return self._send(pages.runs(), HTTPStatus.NOT_FOUND)
            item = json.loads(record.read_text(encoding="utf-8"))
            job = Job(item.get("kind", "analysis"), item.get("name", name), folder)
            if _safe_name(str(item.get("results") or "")):
                job.results = folder / item["results"]
            if isinstance(item.get("seed"), int):
                job.seed = item["seed"]
            if _safe_name(str(item.get("id") or "")) and item["id"] not in session.jobs:
                job.id = item["id"]
            job.status = item.get("status", "done")
            job.error = item.get("error")
            job.started = float(item.get("started") or job.started)
            job.finished = float(item.get("finished") or job.started)
            job.summary = item.get("summary") or {}
            log = folder / "gui.log"
            if log.is_file():
                job.lines = log.read_text(encoding="utf-8", errors="replace").splitlines()
            session.jobs[job.id] = job
            self.send_response(HTTPStatus.SEE_OTHER)
            self.send_header("Location", pages.url(f"/job/{job.id}"))
            self.send_header("Content-Length", "0")
            self.end_headers()

    server = _LoopbackServer(("127.0.0.1", port), Handler)
    address = f"http://127.0.0.1:{server.server_address[1]}/?t={session.token}"
    # a session with example collections opens on them, so that the first run is one click away
    landing = (f"http://127.0.0.1:{server.server_address[1]}/examples?t={session.token}"
               if session.example_runs() else address)
    if ready is not None:
        ready(address)
    else:
        print(f"AMR-ClonalShare {__version__} is open at {address}", flush=True)
        print("If the browser did not open, or its tab was closed, paste that address into the browser.", flush=True)
        print(f"Runs are kept under {session.runs}", flush=True)
        print("Close this window, or press Ctrl+C, to stop.", flush=True)
    if open_browser:
        threading.Timer(0.6, lambda: webbrowser.open(landing)).start()
    asked = False
    while True:
        try:
            server.serve_forever()
            break
        except KeyboardInterrupt:
            running = [j.name for j in session.jobs.values() if j.status in ("queued", "running")]
            if running and not asked:
                asked = True
                print(f"A run is in progress ({', '.join(running)}); press Ctrl+C again to stop it and close.", flush=True)
                continue
            break
    server.server_close()


def serve(root: Path, port: int = 0, open_browser: bool = False, threads: Optional[int] = None,
          ready: Optional[Callable[[str], None]] = None, examples: Optional[Path] = None,
          manual: Optional[Path] = None) -> None:
    """Serve the form until interrupted; ``ready`` receives the address."""
    session = Session(root, threads, examples)
    session.manual = manual if manual is not None and (manual / "index.html").is_file() else None
    _serve(session, port, open_browser, ready)


def _open_folder(folder: Path) -> bool:
    """Open a folder in the file manager of this computer; False when no
    way to do so is known here."""
    import subprocess
    try:
        if sys.platform.startswith("win"):
            os.startfile(str(folder))  # type: ignore[attr-defined]
            return True
        opener = "open" if sys.platform == "darwin" else "xdg-open"
        subprocess.Popen([opener, str(folder)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except (OSError, AttributeError):
        return False


def default_root() -> Path:
    home = Path(os.environ.get("AMR_CLONALSHARE_HOME") or Path.home() / "AMR-ClonalShare")
    return home


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="amr-clonalshare-gui", description=__doc__.splitlines()[0])
    parser.add_argument("--root", default=None, help="working folder for uploads and runs (default: ~/AMR-ClonalShare)")
    parser.add_argument("--port", type=int, default=0, help="port on 127.0.0.1 (default: a free one)")
    parser.add_argument("--threads", type=int, default=None, help="processor threads of a run (default: all)")
    parser.add_argument("--no-browser", action="store_true", help="print the address and do not open the browser")
    parser.add_argument("--examples", default=None,
                        help="folder of the example collections listed on the Examples page "
                             "(default: the examples folder of a source checkout, if the package runs from one)")
    parser.add_argument("--manual", default=None,
                        help="folder of the manual as HTML pages, served under /manual (the Windows folder ships one)")
    parser.add_argument("--version", action="version", version=f"amr-clonalshare {__version__}")
    args = parser.parse_args(argv)
    if args.threads is not None and args.threads < 1:
        print("--threads must be at least 1", file=sys.stderr)
        return 2
    if args.threads:
        from .cli import _set_threads
        _set_threads(args.threads)
    root = Path(args.root).expanduser().resolve() if args.root else default_root()
    examples = examples_folder()
    if args.examples:
        examples = Path(args.examples).expanduser().resolve()
        if not examples.is_dir():
            print(f"--examples: {examples} is not a folder", file=sys.stderr)
            return 2
    try:
        manual = Path(args.manual).expanduser().resolve() if args.manual else None
        serve(root, port=args.port, open_browser=not args.no_browser, threads=args.threads, examples=examples,
              manual=manual)
    except OSError as error:
        print(f"the server could not start: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
