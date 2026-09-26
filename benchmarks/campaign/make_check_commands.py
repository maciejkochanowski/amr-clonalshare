#!/usr/bin/env python3
"""Write the command files of the model-check, multi-panel and mixture arrays.

    python benchmarks/campaign/make_check_commands.py

As in ``make_population_commands.py``, replicates are split into chunks and one
array task runs PACK chunks at once, one per core; every replicate draws from
its own seed, so the split changes nothing in the results.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmarks.campaign.make_population_commands import OUT, pack
from benchmarks.mic_inference.calibrate_extension import extension_design_grid
from benchmarks.mic_inference.calibrate_mixture import mixture_design_grid
from benchmarks.mic_inference.calibrate_null import null_design_grid
from benchmarks.mic_inference.calibrate_panels import panel_design_grid

#: The S. suis agents of the empirical inputs, in the order of their manifest.
MANIFEST = ROOT / "benchmarks" / "results_mic_release" / "mic_empirical_inputs" / "manifest.json"
INPUTS = "$ROOT/outputs/empirical_inputs"
PY = "$PY benchmarks/mic_inference/"
SHAPE_REPLICATES, PANEL_REPLICATES, PLASMODE_REPLICATES, MIXTURE_REPLICATES = 400, 1000, 500, 500


def chunks(total, size):
    return [(start, min(size, total - start)) for start in range(0, total, size)]


def main() -> int:
    agents = [r["agent"] for r in json.loads(MANIFEST.read_text())["records"]]
    shape = []
    for source, grid in (("null", null_design_grid()), ("extension", extension_design_grid())):
        for cell in grid:
            if cell["reading"] == "exact" or cell["size_pattern"] == "all_singletons":
                continue
            for start, n in chunks(SHAPE_REPLICATES, 10):
                shape.append(f"{PY}shape_check.py --source {source} --cell {cell['cell']} --start {start} "
                             f"--replicates {n} --output $ROOT/outputs/mic_shape/{source}_{cell['cell']:03d}_{start:04d}.json")
    for agent in agents:
        for start, n in chunks(SHAPE_REPLICATES, 5):
            shape.append(f"{PY}shape_check.py --source plasmode --agent {agent} --inputs {INPUTS} --start {start} "
                         f"--replicates {n} --output $ROOT/outputs/mic_shape/plasmode_{agent}_{start:04d}.json")
    panels = []
    for cell in panel_design_grid():
        for start, n in chunks(PANEL_REPLICATES, 10):
            panels.append(f"{PY}calibrate_panels.py --mode multipanel --cell {cell['cell']} --start {start} "
                          f"--replicates {n} --bounds-every 50 "
                          f"--output $ROOT/outputs/mic_panels/multipanel_{cell['cell']:02d}_{start:04d}.json")
    for agent in agents:
        for start, n in chunks(PLASMODE_REPLICATES, 5):
            panels.append(f"{PY}calibrate_panels.py --mode plasmode --agent {agent} --inputs {INPUTS} "
                          f"--start {start} --replicates {n} "
                          f"--output $ROOT/outputs/mic_panels/plasmode_{agent}_{start:04d}.json")
    mixture = []
    for cell in mixture_design_grid():
        for start, n in chunks(MIXTURE_REPLICATES, 5):
            mixture.append(f"{PY}calibrate_mixture.py --cell {cell['cell']} --start {start} --replicates {n} "
                           "--critical-value $(cat $ROOT/outputs/population/critical_value.txt) "
                           f"--output $ROOT/outputs/mic_mixture/cell_{cell['cell']}_{start:04d}.json")
    files = {"mic_shape_packed.txt": shape, "mic_panels_packed.txt": panels, "mic_mixture_packed.txt": mixture}
    for name, commands in files.items():
        packed = pack(commands)
        (OUT / name).write_text("\n".join(packed) + "\n", encoding="utf-8")
        print(name, len(commands), "chunks in", len(packed), "tasks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
