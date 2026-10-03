"""Where the time and the memory of one analysis go.

Total runtime says what to buy, not what to change. This records the cost by
phase -- reading the tables, the per-agent estimators, writing the artefacts --
and the twenty functions with the largest cumulative time, together with the
peak resident set and the largest single allocation. Run it before optimising
anything, and again afterwards.

    OMP_NUM_THREADS=1 python benchmarks/profile_run.py examples/ssuis/config.yaml out.json

The thread limits matter: with a threaded BLAS the profile attributes to
Python what is spent in a pool of worker threads.
"""
from __future__ import annotations

import cProfile
import json
import pstats
import resource
import sys
import tempfile
import time
import tracemalloc
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _phases(config: Path) -> dict:
    """Wall time of the three phases, measured by instrumenting the seams."""
    from amr_clonalshare import core, io, outputs
    marks: dict[str, float] = {}

    def timed(name, function):
        def wrapper(*args, **kwargs):
            start = time.perf_counter()
            try:
                return function(*args, **kwargs)
            finally:
                marks[name] = marks.get(name, 0.0) + time.perf_counter() - start
        return wrapper

    original_load, original_publish = io.load_dataset, outputs.publish_results
    core.load_dataset = timed("read_tables", original_load)
    outputs.publish_results = timed("write_artefacts", original_publish)
    try:
        with tempfile.TemporaryDirectory() as tmp:
            start = time.perf_counter()
            from amr_clonalshare.config import load_config
            core.run(load_config(config), results_dir=Path(tmp) / "out")
            marks["total"] = time.perf_counter() - start
    finally:
        core.load_dataset, outputs.publish_results = original_load, original_publish
    marks["estimators"] = (marks["total"] - marks.get("read_tables", 0.0)
                           - marks.get("write_artefacts", 0.0))
    return marks


def main(argv=None) -> int:
    warnings.simplefilter("ignore")
    argv = list(sys.argv[1:] if argv is None else argv)
    config = Path(argv[0])
    destination = Path(argv[1]) if len(argv) > 1 else None

    phases = _phases(config)

    tracemalloc.start()
    profiler = cProfile.Profile()
    with tempfile.TemporaryDirectory() as tmp:
        from amr_clonalshare import cli
        profiler.enable()
        cli.main(["--config", str(config), "--results-dir", str(Path(tmp) / "out"),
                  "--quiet"])
        profiler.disable()
    traced_peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()

    stats = pstats.Stats(profiler)
    stats.sort_stats("cumulative")
    rows = []
    for (path, line, name), (_calls, primitive, total, cumulative, _c) in list(
            stats.stats.items()):
        rows.append({"function": f"{Path(path).name}:{line}:{name}",
                     "calls": primitive, "tottime": round(total, 4),
                     "cumtime": round(cumulative, 4)})
    rows.sort(key=lambda r: -r["cumtime"])
    report = {
        "config": str(config.relative_to(ROOT) if config.is_absolute() else config),
        "phase_seconds": {k: round(v, 3) for k, v in sorted(phases.items())},
        "peak_rss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1),
        "peak_traced_python_mb": round(traced_peak / 2 ** 20, 1),
        "hot_functions": rows[:20],
        "own_time": sorted(rows, key=lambda r: -r["tottime"])[:15],
    }
    text = json.dumps(report, indent=1)
    if destination is not None:
        destination.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
