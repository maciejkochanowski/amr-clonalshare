"""Coverage of the calibrated MIC test when laboratories read different panels,
and on the design of the S. suis example.

``--mode multipanel``: 30 lineages of uneven size; three laboratories with
    lab0  a wide panel, cut points -4 to 4
    lab1  a narrow panel set low, cut points -2 to 1
    lab2  one cut point at 0 (a breakpoint-only reading)
Laboratory membership is aligned with lineage (85% of a lineage's isolates in
its home laboratory) or random; the laboratories are either not shifted or
shifted by -0.5, 0 and +0.5 total standard deviations, the shift entering the
analysis as a fixed effect. Every dataset is analysed twice: ``perlab`` reads
all records, each on its own panel; ``onelab`` reads the records of lab0 only,
the analysis available without per-laboratory panels. Coverage is the share of
datasets whose test at the generating rho retains it; with ``--bounds-every``
the full interval is also computed, for the widths.

``--mode plasmode``: the lineages, laboratories, inferred panels and countries
of one S. suis agent as recorded; MICs are simulated from the Gaussian model
fitted to that agent (country as a fixed effect), with rho at the fitted value,
and read on the laboratory's panel. This checks coverage on the exact design
of the example; it does not test the Gaussian assumption, which the model
check of ``shape_check.py`` addresses.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from amr_clonalshare._mic_panels import _record_offset, _validated_problem, observe_panels
from amr_clonalshare.mic_null_bootstrap import (_unrestricted_fit, null_bootstrap_interval,
                                                null_bootstrap_test)

ROOT_SEED = 20261001
PLASMODE_SEED = 20261002
PANELS = {"lab0": [float(x) for x in np.arange(-4., 5.)], "lab1": [-2., -1., 0., 1.], "lab2": [0.]}
SHIFT = {"lab0": -.5, "lab1": 0., "lab2": .5}
LABS = ("lab0", "lab1", "lab2")


def panel_design_grid():
    cells = []
    for rho in (.1, .5, .9):
        for alignment in ("aligned", "random"):
            for shift in (False, True):
                cells.append(dict(cell=len(cells), rho=rho, alignment=alignment, shift=shift))
    return cells


def generate(cell, replicate, root_seed=ROOT_SEED):
    rng = np.random.default_rng(np.random.SeedSequence([root_seed, cell["cell"], replicate]))
    groups = 30
    sizes = np.rint(np.exp(rng.normal(2.2, 1., groups))).astype(int).clip(2, 100)
    labels = np.repeat(np.arange(groups), sizes)
    home = rng.integers(0, 3, groups)
    if cell["alignment"] == "aligned":
        stay = rng.random(len(labels)) < .85
        lab = np.where(stay, home[labels], rng.integers(0, 3, len(labels)))
    else:
        lab = rng.integers(0, 3, len(labels))
    lab = np.asarray(LABS)[lab]
    rho = cell["rho"]
    y = rng.normal(0, np.sqrt(rho), groups)[labels] + rng.normal(0, np.sqrt(1 - rho), len(labels))
    if cell["shift"]:
        y = y + np.array([SHIFT[x] for x in lab])
    names = sorted(set(lab))
    a, b = observe_panels(y, ([np.asarray(PANELS[n]) for n in names], np.searchsorted(names, lab)))
    seed = int(np.random.SeedSequence([root_seed, cell["cell"], replicate, 991]).generate_state(1)[0])
    return a, b, labels, lab, seed


def read_agent(inputs: Path, agent: str):
    records = json.loads((inputs / "manifest.json").read_text())["records"]
    record = next(r for r in records if r["agent"] == agent)
    source = inputs / record["file"]
    if hashlib.sha256(source.read_bytes()).hexdigest() != record["sha256"]:
        raise ValueError(f"{source} does not match its manifest")
    with source.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    lo = np.array([float(r["lo"]) for r in rows])
    hi = np.array([float(r["hi"]) for r in rows])
    labels = np.array([r["lineage"] for r in rows])
    panel = np.array([r["panel"] for r in rows])
    country = np.array([r["country"] for r in rows])
    return record, lo, hi, labels, panel, country


def plasmode_state(inputs: Path, agent: str):
    """The recorded design of one agent and the Gaussian model fitted to it."""
    record, lo, hi, labels, panel, country = read_agent(inputs, agent)
    problem, panels, _ = _validated_problem(lo, hi, labels, record["panel_edges"], panel, country)
    fit = _unrestricted_fit(problem, 24)
    return dict(record=record, labels=labels, panel=panel, country=country, problem=problem,
                panels=panels, fit=fit, offset=_record_offset(problem, fit))


def plasmode_dataset(state, replicate):
    """Readings simulated from the fitted model on the recorded design, and the
    seed of their analysis."""
    record, problem, fit = state["record"], state["problem"], state["fit"]
    rng = np.random.default_rng(np.random.SeedSequence([PLASMODE_SEED, int(record["seed"]), replicate]))
    y = (fit.mean + state["offset"] + rng.normal(0, fit.total_sd * np.sqrt(fit.rho), problem.groups)[problem.code]
         + rng.normal(0, fit.total_sd * np.sqrt(1 - fit.rho), problem.n))
    a, b = observe_panels(y, state["panels"])
    seed = int(np.random.SeedSequence([PLASMODE_SEED, int(record["seed"]), replicate, 991]).generate_state(1)[0])
    return a, b, seed


def _run(arm, a, b, labels, panel_edges, panel, covariate, rho, seed, bounds):
    row = dict(arm=arm, n=int(len(a)))
    started = time.perf_counter()
    try:
        test = null_bootstrap_test(a, b, labels, rho=rho, panel_edges=panel_edges, panel=panel,
                                   covariate=covariate, n_boot=199, seed=seed)
        row.update(reportable=bool(test.reportable), covered=bool(test.accepted),
                   statistic=test.statistic, critical=test.critical,
                   failed=test.bootstrap_failed, redrawn=test.bootstrap_redrawn,
                   rho_hat=test.unrestricted_fit.rho)
        if bounds:
            interval = null_bootstrap_interval(a, b, labels, panel_edges=panel_edges, panel=panel,
                                               covariate=covariate, n_boot=199, seed=seed)
            row.update(low=interval.low, high=interval.high,
                       interval_reportable=bool(interval.reportable), status=interval.status)
    except (ValueError, ArithmeticError) as exc:
        # a dataset the analysis refuses (a laboratory censored on one side in
        # every reading) is counted as refused, not as covered
        row.update(reportable=False, covered=False, refusal=f"{type(exc).__name__}: {exc}")
    row["seconds"] = time.perf_counter() - started
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--mode", choices=["multipanel", "plasmode"], required=True)
    parser.add_argument("--cell", type=int, default=0)
    parser.add_argument("--agent")
    parser.add_argument("--inputs", type=Path, help="directory written by prepare_empirical.py")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--bounds-every", type=int, default=0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    rows = []
    if args.mode == "plasmode":
        state = plasmode_state(args.inputs, args.agent)
        record, labels, panel, country, fit = (state[k] for k in ("record", "labels", "panel", "country", "fit"))
    for rep in range(args.start, args.start + args.replicates):
        bounds = bool(args.bounds_every and rep % args.bounds_every == 0)
        if args.mode == "multipanel":
            cell = panel_design_grid()[args.cell]
            a, b, labels, lab, seed = generate(cell, rep)
            base = dict(mode="multipanel", cell=args.cell, replicate=rep, rho=cell["rho"],
                        alignment=cell["alignment"], shift=cell["shift"])
            rows.append({**base, **_run("perlab", a, b, labels, PANELS, lab,
                                        lab if cell["shift"] else None, cell["rho"], seed, bounds)})
            one = lab == "lab0"
            rows.append({**base, **_run("onelab", a[one], b[one], labels[one], PANELS["lab0"], None,
                                        None, cell["rho"], seed, bounds)})
        else:
            a, b, seed = plasmode_dataset(state, rep)
            base = dict(mode="plasmode", agent=args.agent, replicate=rep, rho=fit.rho)
            rows.append({**base, **_run("plasmode", a, b, labels, record["panel_edges"], panel, country,
                                        fit.rho, seed, bounds)})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(rows, indent=1, allow_nan=True, default=float))


if __name__ == "__main__":
    main()
