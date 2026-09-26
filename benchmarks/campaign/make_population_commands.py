#!/usr/bin/env python3
"""Write the command files of the population-model arrays.

    python benchmarks/campaign/make_population_commands.py

Replicates are split into chunks and the chunks packed so that one array task
runs PACK of them at once, one per core; the split changes nothing in the
results, since every replicate draws from its own seed.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmarks.population_model.design import cells

OUT = ROOT / "benchmarks" / "campaign" / "commands"
PACK = 48
RUN = ("$PY -m benchmarks.population_model.run_cell --cell {c} --start {s} --replicates {r} "
       "--method {m}{extra} --output $ROOT/outputs/population/{d}/cell_{c:03d}_{s:04d}.csv")
CRITICAL = " --critical-value $(cat $ROOT/outputs/population/critical_value.txt)"


def lines(families, method, chunk, prefix, extra="", replicates=None):
    out = []
    for cell in cells():
        if cell["family"] not in families:
            continue
        total = replicates or cell["replicates"]
        directory = prefix if len(families) == 1 else f"{prefix}_{cell['family']}"
        for start in range(0, total, chunk):
            out.append(RUN.format(c=cell["cell"], s=start, r=min(chunk, total - start), m=method,
                                  extra=extra, d=directory))
    return out


def pack(commands):
    packed = []
    for i in range(0, len(commands), PACK):
        body = " ".join(f"{{ {c}; }} & P+=($!);" for c in commands[i:i + PACK])
        packed.append(f'P=(); F=0; {body} for p in "${{P[@]}}"; do wait $p || F=1; done; [ $F -eq 0 ]')
    return packed


def main() -> int:
    files = {
        "population_calibration_packed.txt": lines({"calibration"}, "lr", 50, "lr_calibration"),
        "population_validation_lr_packed.txt": lines({"validation", "robustness"}, "lr", 50, "lr", CRITICAL),
        "population_validation_general_packed.txt": lines({"validation", "robustness"}, "general", 5,
                                                          "general", replicates=500),
        "population_benefit_packed.txt": lines({"benefit"}, "benefit", 20, "benefit"),
        "population_mixing_packed.txt": lines({"validation", "robustness"}, "mixing", 10, "mixing",
                                              CRITICAL, replicates=500),
    }
    for name, commands in files.items():
        packed = pack(commands)
        (OUT / name).write_text("\n".join(packed) + "\n", encoding="utf-8")
        print(name, len(commands), "chunks in", len(packed), "tasks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
