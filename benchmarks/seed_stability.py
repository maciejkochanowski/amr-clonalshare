"""How much of the reported ordering is the seed?

The per-agent share is cross-validated and its control is a permutation, so
the reported number carries a Monte Carlo error of its own, beside the
sampling uncertainty the interval describes. This script runs the recorded
analysis under many seeds and writes, for every antimicrobial, the spread of
the estimate and how often it attains the highest and the lowest share. A
ranking that changes with the seed is a ranking the data do not resolve; one
that holds across every seed may be named.

    OMP_NUM_THREADS=1 python benchmarks/seed_stability.py 40 18 out.json

The thread limits are not decoration: a threaded BLAS makes the timing and,
through the order of its reductions, the last digits vary between runs.
"""
from __future__ import annotations

import json
import statistics
import sys
import tempfile
import warnings
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "examples" / "ssuis" / "config.yaml"


def one(seed: int):
    warnings.simplefilter("ignore")
    from amr_clonalshare import cli
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "out"
        status = cli.main(["--config", str(CONFIG), "--results-dir", str(out),
                           "--seed", str(seed), "--quiet"])
        assert status == 0, status
        record = json.loads((out / "clonal_share_result.json").read_text(encoding="utf-8"))
    share = record["metadata_diagnostics"]["clonal_share"]
    return seed, {agent: block["kappa_adj"] for agent, block in share.items()
                  if block.get("kappa_adj") is not None}


def summarise(runs: dict) -> dict:
    values = defaultdict(list)
    highest, lowest = Counter(), Counter()
    for share in runs.values():
        for agent, value in share.items():
            values[agent].append(value)
        highest[max(share, key=share.get)] += 1
        lowest[min(share, key=share.get)] += 1
    per_agent = {
        agent: {"median": statistics.median(v), "min": min(v), "max": max(v),
                "sd": statistics.pstdev(v)}
        for agent, v in sorted(values.items())}
    return {"n_seeds": len(runs), "per_agent": per_agent,
            "highest_share": dict(highest), "lowest_share": dict(lowest)}


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    n_seeds, workers, destination = int(argv[0]), int(argv[1]), Path(argv[2])
    runs = {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for seed, share in pool.map(one, range(1, n_seeds + 1)):
            runs[str(seed)] = share
    destination.write_text(json.dumps(runs, indent=1), encoding="utf-8")
    print(json.dumps(summarise(runs), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
