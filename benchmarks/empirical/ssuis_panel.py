"""What the S. suis MIC panel resolves that coarser panels and the call do not.

The 628 isolates read by LGC Fordham are read on one panel per agent. Every
agent is analysed on that panel, on the panels that keep every second and
every fourth of its cut points, and, for an agent with a cut-off on that
panel, on the one cut of its call. The cut points kept are chosen so that
the cut-off is among them, so the four panels are nested: every reading of a
coarser panel is the union of readings of the finer one. For each panel the
lineage (population cluster) share of the ordering of the readings, with its
interval for the represented lineages, and the bounds on the lineage share of
the latent MIC ordering with its lower confidence limit are computed
within the country of isolation, the covariate of the pipeline. On the same
isolates and strata the lower bound the readings establish can only fall as the
panel coarsens (Theorem 1), so the steps down measure what each coarsening
loses.

    python benchmarks/empirical/ssuis_panel.py <output dir>
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from amr_clonalshare.censored import intervals_by_panel
from amr_clonalshare.mic_order import mic_order_share

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "examples/ssuis/data"
SEED = 20260915
LAB = "LGC Fordham"
KEYS = ("kappa_adj", "observed_low", "observed_high", "p_value", "n_readings",
        "latent_order_lower", "latent_order_upper", "latent_order_upper_bound",
        "latent_order_lower_limit", "latent_order_midpoint", "panel_resolution", "share_end_wells")


def _cutoffs():
    """The cut-off of every agent with a call: the published value where
    there is one, the fitted value otherwise (DATA_PROVENANCE.md)."""
    out = {}
    for r in json.loads((DATA / "ecoff_derived.json").read_text()):
        out[r["agent"]] = r["published_value"] if r["published_value"] is not None else r["fit"]["ecoff"]
    return out


def coarsen(lo, hi, kept):
    """Each reading (lo, hi] widened to the panel whose cut points are
    ``kept``: from the largest kept cut at or below lo to the smallest kept
    cut at or above hi."""
    kept = np.sort(np.asarray(kept, dtype=float))
    ext = np.r_[-np.inf, kept, np.inf]
    new_lo = ext[np.searchsorted(ext, lo, side="right") - 1]
    new_hi = ext[np.searchsorted(ext, hi, side="left")]
    return new_lo, new_hi


def main(out: Path) -> int:
    started = time.perf_counter()
    out.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(DATA / "metadata.csv", dtype=str, keep_default_na=False).set_index("genome_id")
    mic = pd.read_csv(DATA / "mic_panel.csv", dtype=str, keep_default_na=False)
    mic = mic[mic.testing_laboratory.eq(LAB)]
    cutoffs = _cutoffs()
    rows = []
    for agent, frame in mic.groupby("antibiotic", sort=True):
        frame = frame.set_index("genome_id")
        ids = frame.index.to_numpy()
        lineage = meta.loc[ids, "baps_cluster"].to_numpy(dtype=object)
        country = meta.loc[ids, "isolation_country"].to_numpy(dtype=object)
        lo, hi = intervals_by_panel(frame.measurement.astype(float), frame.testing_laboratory,
                                    operators=frame.measurement_sign)
        cuts = np.unique(np.r_[lo[np.isfinite(lo)], hi[np.isfinite(hi)]])
        cut = np.log2(float(cutoffs[agent])) if agent in cutoffs else None
        on_panel = cut is not None and bool(np.any(np.isclose(cuts, cut)))
        anchor = int(np.argmin(np.abs(cuts - cut))) if on_panel else 0
        panels = {"panel": cuts}
        for step in (2, 4):
            panels[f"every {step}"] = cuts[(np.arange(cuts.size) - anchor) % step == 0]
        if on_panel:
            panels["call"] = cuts[[anchor]]
        for name, kept in panels.items():
            print(agent, name, flush=True)
            clo, chi = coarsen(lo, hi, kept)
            share = mic_order_share(clo, chi, lineage, strata=country, stratified_by="isolation_country",
                                    seed=SEED)
            rows.append(dict(agent=agent, panel=name, n_cuts=int(np.size(kept)), n=int(ids.size),
                             cutoff_log2=cut if on_panel else None, **{k: share[k] for k in KEYS}))
    table = pd.DataFrame(rows)
    table.to_csv(out / "ssuis_panel.csv", index=False)
    files = [DATA / "metadata.csv", DATA / "mic_panel.csv", DATA / "ecoff_derived.json", Path(__file__)]
    files += sorted((ROOT / "src/amr_clonalshare").glob("*.py"))
    receipt = dict(seed=SEED, laboratory=LAB, strata="isolation_country",
                   elapsed_seconds=time.perf_counter() - started,
                   source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in files})
    (out / "receipt.json").write_text(json.dumps(receipt, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1])))
