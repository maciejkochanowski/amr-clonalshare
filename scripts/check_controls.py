#!/usr/bin/env python3
"""Assert the planted control behaves as a control.

Two agents follow the lineage and one does not. The share must separate them,
the permuted-label control must manufacture no share on any of them, and the
run must write the record and both reports. The cohort is written to a
temporary directory, so the check needs nothing the repository does not ship.

    python scripts/check_controls.py [--out DIR]
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def write_cohort(root: Path, *, n_per_lineage=12, n_lineages=8, seed=11) -> Path:
    rng = np.random.default_rng(seed)
    ids, lineages = [], []
    for g in range(n_lineages):
        for i in range(n_per_lineage):
            ids.append(f"iso{g:02d}_{i:02d}")
            lineages.append(f"L{g:02d}")
    carrier = {f"L{g:02d}": int(g < n_lineages // 2) for g in range(n_lineages)}
    rows = []
    for iso, lin in zip(ids, lineages):
        for agent, clonal in (("agentA", True), ("agentB", True),
                              ("agentC", False)):
            p = (0.9 if carrier[lin] else 0.1) if clonal else 0.5
            call = "non-susceptible" if rng.random() < p else "susceptible"
            rows.append({"iid": iso, "antibiotic": agent, "call": call})
    data = root / "data"
    data.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(data / "calls.csv", index=False)
    pd.DataFrame({"iid": ids, "lineage": lineages,
                  "intake": [f"20{10 + (i % 3)}" for i in range(len(ids))]}
                 ).to_csv(data / "meta.csv", index=False)
    raw = {
        "dataset": {"name": "planted_control", "strain_id_column": "iid",
                    "data_dir": "data", "metadata": "meta.csv",
                    "lineage_column": "lineage", "batch_column": "intake",
                    "phenotype": "calls.csv", "phenotype_id_column": "iid",
                    "phenotype_antibiotic_column": "antibiotic",
                    "phenotype_call_column": "call"},
        "attribution": {"n_boot": 100, "n_perm": 100},
        "evidence": {"folds": 3, "repeats": 2},
    }
    cfg = root / "planted.yaml"
    cfg.write_text(yaml.safe_dump(raw, sort_keys=False))
    return cfg


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None,
                    help="results directory (default: a temporary one)")
    a = ap.parse_args(argv)
    work = Path(tempfile.mkdtemp(prefix="planted-"))
    out = Path(a.out) if a.out else work / "out"
    cfg = write_cohort(work)
    rc = subprocess.call([sys.executable, "-m", "amr_clonalshare.cli",
                          "--config", str(cfg), "--results-dir", str(out),
                          "--seed", "7", "--quiet", "--threads", "1"])
    if rc != 0:
        print(f"planted control: the run exited {rc}", file=sys.stderr)
        return 1
    record = json.loads((out / "clonal_share_result.json").read_text())
    share = record["metadata_diagnostics"]["clonal_share"]
    failures = []
    for name in ("agentA", "agentB"):
        if share[name]["kappa_adj"] - share[name]["null_mean"] <= 0.30:
            failures.append(f"{name}: share {share[name]['kappa_adj']:.3f} does "
                            f"not separate from its control {share[name]['null_mean']:.3f}")
    if share["agentC"]["kappa_adj"] >= 0.15:
        failures.append(f"agentC: share {share['agentC']['kappa_adj']:.3f} on an "
                        "agent that does not follow the lineage")
    for name, d in share.items():
        if not -0.15 < d["null_mean"] < 0.05:
            failures.append(f"{name}: permuted control {d['null_mean']:.3f}")
    for f in ("report.md", "report.html", "input_qc.json", "input_qc.md"):
        if not (out / f).is_file():
            failures.append(f"{f} was not written")
    if failures:
        print("planted control FAILED:", *failures, sep="\n  ", file=sys.stderr)
        return 1
    print("planted control OK: agentA %.3f, agentB %.3f follow the lineage; "
          "agentC %.3f does not; controls within band"
          % (share["agentA"]["kappa_adj"], share["agentB"]["kappa_adj"],
             share["agentC"]["kappa_adj"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
