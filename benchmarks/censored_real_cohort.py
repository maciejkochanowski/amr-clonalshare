#!/usr/bin/env python3
"""The coarsening ladder on the shipped Streptococcus suis cohort.

    python benchmarks/censored_real_cohort.py --out <dir>

For every antimicrobial with both a dichotomised non-wild-type call and a
recorded dilution, the same estimator runs twice: on the call, which is a
single cut point, and on the panel, which is a set of intervals. Both report
the latent variance share, so what separates them is the resolution of the
reading and nothing else. A comparison against the package's out-of-sample
predictive share would compare two estimands and is deliberately not made here.

Agents with a panel but no epidemiological cut-off are reported as well. They
have no call to compare against, which is the point: they are analysable at
all only from the dilution.
"""
from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from amr_clonalshare.censored import (censored_clonal_share,
                                         intervals_from_binary,
                                         intervals_from_mic, panel_geometry)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "examples" / "ssuis" / "data"
LINEAGE_COLUMN = "baps_cluster"
SEED = 0


def width(record: dict, high: str, low: str) -> float:
    a, b = record[high], record[low]
    return a - b if np.isfinite(a) and np.isfinite(b) else float("nan")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path,
                        default=ROOT / "benchmarks" / "results_censored")
    parser.add_argument("--n-boot", type=int, default=400)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(DATA / "metadata.csv", dtype=str).set_index("genome_id")
    calls = (pd.read_csv(DATA / "ribo.csv", dtype={"genome_id": str})
             .set_index("genome_id")
             .join(pd.read_csv(DATA / "cell.csv", dtype={"genome_id": str})
                   .set_index("genome_id"), how="inner"))
    mic = pd.read_csv(DATA / "mic_panel.csv", dtype={"genome_id": str})
    ids = calls.index
    lineage = meta.reindex(ids)[LINEAGE_COLUMN].to_numpy(dtype=object)

    rows = []
    for agent in sorted(mic["antibiotic"].astype(str).unique()):
        block = (mic[mic["antibiotic"] == agent]
                 .drop_duplicates(subset=["genome_id"])
                 .set_index("genome_id"))
        values = pd.to_numeric(block["measurement"],
                               errors="coerce").reindex(ids)
        present = values.notna().to_numpy()
        if present.sum() < 8:
            continue
        v = values.to_numpy(dtype=float)[present]
        lin = lineage[present]
        entry: dict = {"agent": agent, "n": int(present.sum()),
                       "has_call": bool(agent in calls.columns)}
        cut = float("nan")
        if entry["has_call"]:
            y = calls[agent].to_numpy(dtype=float)[present]
            entry["prevalence"] = float(y.mean())
            wild = v[y < 0.5]
            # The cut-off is read back from the shipped call rather than
            # assumed: the highest tested dilution the cohort calls wild type.
            cut = float(np.log2(wild.max())) if wild.size else float("nan")
            entry["cutoff_log2"] = cut
            lo, hi = intervals_from_binary(
                y, cutoff_log2=cut if np.isfinite(cut) else 0.0)
            entry["binary"] = censored_clonal_share(
                lo, hi, lin, n_boot=args.n_boot, seed=SEED).as_dict()
        lo, hi = intervals_from_mic(v)
        entry["mic"] = censored_clonal_share(
            lo, hi, lin, n_boot=args.n_boot, seed=SEED).as_dict()
        entry["panel"] = panel_geometry(
            v, cutoff=2.0 ** cut if np.isfinite(cut) else None).as_dict()
        rows.append(entry)

    by_agent = {r["agent"]: r for r in rows}
    paired = [r for r in rows if "binary" in r]
    kb = np.array([r["binary"]["kappa"] for r in paired])
    km = np.array([r["mic"]["kappa"] for r in paired])
    pb = np.array([width(r["binary"], "profile_high", "profile_low")
                   for r in paired])
    pm = np.array([width(r["mic"], "profile_high", "profile_low")
                   for r in paired])
    ib = np.array([width(r["binary"], "ci_high", "ci_low") for r in paired])
    im = np.array([width(r["mic"], "ci_high", "ci_low") for r in paired])
    # The same fit answers the other question: how much sits between the
    # thirty lineages this cohort holds, rather than between lineages drawn
    # afresh. At thirty lineages that is the difference between an interval a
    # laboratory can act on and one it cannot.
    rm_ = np.array([width(r["mic"], "realised_high", "realised_low")
                    for r in rows])
    im_all = np.array([width(r["mic"], "ci_high", "ci_low") for r in rows])
    usable = np.isfinite(rm_) & np.isfinite(im_all)
    bb = np.array([width(r["binary"], "boot_high", "boot_low") for r in paired])
    bm = np.array([width(r["mic"], "boot_high", "boot_low") for r in paired])
    summary = {
        "n_isolates": int(len(ids)),
        "n_lineages": int(pd.Series(lineage).nunique()),
        "lineage_column": LINEAGE_COLUMN,
        "n_agents_with_panel": len(rows),
        "n_agents_paired": len(paired),
        "n_agents_panel_only": len(rows) - len(paired),
        "n_agents_call_refused": sum(
            1 for r in paired if not r["binary"]["estimable"]),
        "n_agents_dilution_refused": sum(
            1 for r in rows if not r["mic"]["estimable"]),
        "n_agents_with_realised_interval": int(usable.sum()),
        "median_realised_width_mic": float(np.median(rm_[usable])),
        "median_superpopulation_width_mic": float(np.median(im_all[usable])),
        "realised_width_ratio_mic": float(np.median(im_all[usable])
                                          / np.median(rm_[usable])),
        "agents_with_narrower_realised_interval": int(
            np.sum(rm_[usable] < im_all[usable])),
        "pearson_kappa": float(np.corrcoef(kb, km)[0, 1]),
        "median_interval_width_binary": float(np.nanmedian(ib)),
        "median_interval_width_mic": float(np.nanmedian(im)),
        "agents_with_narrower_mic_interval": int(np.nansum(im < ib)),
        "median_profile_width_binary": float(np.nanmedian(pb)),
        "median_profile_width_mic": float(np.nanmedian(pm)),
        "agents_with_narrower_mic_profile": int(np.nansum(pm < pb)),
        "median_bootstrap_width_binary": float(np.nanmedian(bb)),
        "median_bootstrap_width_mic": float(np.nanmedian(bm)),
        "agents_with_narrower_mic_bootstrap": int(np.nansum(bm < bb)),
    }
    payload = {"generated": date.today().isoformat(),
               "bootstrap_draws": args.n_boot,
               "summary": summary, "per_agent": by_agent}
    target = args.out / "censored_real_cohort.json"
    target.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("written", target)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
