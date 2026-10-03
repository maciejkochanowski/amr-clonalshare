#!/usr/bin/env python3
"""Write the queue files of the confirmatory campaign (CONFIRMATORY_PROTOCOL.md).

    python benchmarks/campaign/make_confirmatory_queues.py "<command that sets up R>"

for instance ``"module load GCC/13.2.0 R/4.3.3"``; the R packages are those
``competitors/install_r_bench.R`` installs in ``$ROOT/rlib``. Writes, under
``benchmarks/campaign/commands``:

``confirmatory_main.txt``         Studies 1 to 3, a slow grid cell in parts over
                                  the same datasets;
``confirmatory_competitors.txt``  Study 4, which runs R;
``confirmatory_summaries.txt``    the merges, the summaries and last the
                                  verdicts of the protocol, in order.

Every line of a queue is ``task<TAB>cores<TAB>command``; the commands are dealt
to tasks longest first, each to the task with the least planned work, so that
every task fills one node for about TARGET seconds.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "benchmarks" / "conditional"))
from benchmarks.competitors.compare_order import confirmatory_designs
from benchmarks.estimator_benchmark import design as grid_design
from benchmarks.order_calibration import order_design_grid
import designs as conditional

OUT = ROOT / "benchmarks" / "campaign" / "commands"
#: Planning figures in core-seconds per replicate.
COST = {"grid": {5: .5, 10: .6, 30: .8, 100: 1.3, 300: 1.5, 1000: 3},
        "grid_calls": {5: .6, 10: .7, 30: 1, 100: 1.3, 300: 3, 1000: 15}, "order": 1.5,
        "observed": .6, "latent": .3, "units": 6, "comparison": 7, "competitors": 15}
#: Parts of a grid cell by its number of lineages; a part scores every
#: parts-th replicate of the cell's one stream.
GRID_PARTS = {30: 2, 100: 4, 300: 8, 1000: 32}
#: Replicates per command.
CHUNK = {"order": 250, "observed": 250, "latent": 500, "units": 100, "comparison": 50, "competitors": 50}
REPLICATES = {"observed": 10000, "latent": 10000, "units": 20000, "comparison": 5000}
CORES = 192
TARGET = 45 * 60


def item(cores, cost, command):
    return dict(cores=cores, cost=cost, command=command)


def chunks(n, size):
    return [(start, min(size, n - start)) for start in range(0, n, size)]


def deal(items):
    tasks = max(1, math.ceil(sum(i["cost"] for i in items) / (CORES * TARGET)))
    load = [0.] * tasks
    lines = []
    for it in sorted(items, key=lambda i: -i["cost"]):
        t = min(range(tasks), key=load.__getitem__)
        load[t] += it["cost"] / it["cores"]
        lines.append(f"{t}\t{it['cores']}\t{it['command']}")
    return tasks, sorted(lines, key=lambda line: int(line.split("\t", 1)[0]))


def grid():
    out = []
    for cell in grid_design():
        parts = GRID_PARTS.get(cell["n_groups"], 1)
        cost = COST["grid_calls" if cell["binary"] else "grid"][cell["n_groups"]] * cell["replicates"] / parts
        base = f"$PY benchmarks/estimator_benchmark.py --cell {cell['index']} --out $ROOT/outputs/estimator_grid"
        out += ([item(1, cost, base)] if parts == 1 else
                [item(1, cost, f"{base} --parts {parts} --part {k}") for k in range(parts)])
    return out


def order():
    return [item(1, COST["order"] * n,
                 f"$PY benchmarks/order_calibration.py --cell {cell['cell']} --start {start} --replicates {n} "
                 f"--output $ROOT/outputs/order/{cell['design']}_{start:05d}.json")
            for cell in order_design_grid() for start, n in chunks(cell["replicates"], CHUNK["order"])]


def conditional_study():
    out = []
    families = {"observed": conditional.observed_designs, "latent": conditional.latent_designs,
                "units": conditional.unit_designs, "comparison": conditional.comparison_designs}
    for family, make in families.items():
        for d in make():
            for start, n in chunks(REPLICATES[family], CHUNK[family]):
                out.append(item(1, COST[family] * n,
                                f"$PY benchmarks/conditional/run.py --family {family} --design {d.name} "
                                f"--start {start} --replicates {n} "
                                f"--output $ROOT/outputs/conditional/{family}/{d.name}_{start:05d}.json"))
    return out


def competitors(setup):
    return [item(1, COST["competitors"] * n,
                 f"{setup} >/dev/null 2>&1 && R_LIBS_BENCH=$ROOT/rlib $PY -m benchmarks.competitors.compare_order "
                 f"--design {name} --start {start} --replicates {n} "
                 f"--output $ROOT/outputs/competitors/{name}_{start:05d}.csv")
            for name, cell in sorted(confirmatory_designs().items())
            for start, n in chunks(cell["replicates"], CHUNK["competitors"])]


def summaries(setup):
    logged = "$PY benchmarks/campaign/run_logged.py --receipt {out}/RUN_RECEIPT.json --outputs {out} -- "
    merges = [f"$PY benchmarks/estimator_benchmark.py --cell {c['index']} --merge --out $ROOT/outputs/estimator_grid"
              for c in grid_design() if GRID_PARTS.get(c["n_groups"], 1) > 1]
    return merges + [
        logged.format(out="$ROOT/outputs/estimator_grid") + "$PY benchmarks/estimator_benchmark.py --aggregate "
        "--out $ROOT/outputs/estimator_grid",
        f"{setup} >/dev/null 2>&1 && " + logged.format(out="$ROOT/outputs/competitors_summary")
        + "env R_LIBS_BENCH=$ROOT/rlib $PY -m benchmarks.competitors.compare_order --summarize "
        "$ROOT/outputs/competitors --output $ROOT/outputs/competitors_summary",
        "mkdir -p $ROOT/outputs/verdicts && " + logged.format(out="$ROOT/outputs/verdicts")
        + "$PY benchmarks/confirmatory_verdicts.py --grid $ROOT/outputs/estimator_grid/estimator_benchmark.json "
        "--order $ROOT/outputs/order --conditional $ROOT/outputs/conditional --output $ROOT/outputs/verdicts",
    ]


def main(argv) -> int:
    if len(argv) != 1:
        print(__doc__, file=sys.stderr)
        return 2
    OUT.mkdir(parents=True, exist_ok=True)
    for name, items in (("confirmatory_main.txt", grid() + order() + conditional_study()),
                        ("confirmatory_competitors.txt", competitors(argv[0]))):
        tasks, text = deal(items)
        (OUT / name).write_text("\n".join(text) + "\n", encoding="utf-8")
        print(f"{name}: {len(items)} commands, {tasks} tasks, "
              f"{sum(i['cost'] for i in items) / 3600:.0f} planned core-hours")
    (OUT / "confirmatory_summaries.txt").write_text("\n".join(summaries(argv[0])) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
