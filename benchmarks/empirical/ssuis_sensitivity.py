"""Sensitivity of the S. suis MIC and call results to the choices the analysis makes.

Every MIC agent is analysed under four settings of the lineage share of the
MIC ordering and of the bounds the readings place on the share of the
latent ordering:

``baseline``        scored and permuted within testing laboratory × country of
                    isolation, as the pipeline runs (ssuis_analysis.py);
``collection``      scored and permuted within source collection × clinical
                    category (Genome.Set × Pathogen of the source table), the
                    strata of the sampling design;
``without_rcuk_uk`` the baseline without the RCUK_UK collection, 124 isolates
                    all non-clinical and from one sampling campaign, the part
                    of the collection most likely to hold several isolates of a
                    farm;
``third_panel``     marbofloxacin only: the Lola_UK collection read on a panel
                    of its own (it alone records readings of 0.015 mg/L), so
                    that the end wells are found within that panel.

Every agent with a call is analysed at its cut-off and one dilution either
side: the share of the call with its interval for the represented lineages, and
the bounds the call places on the share of the latent ordering, pooled as
the pipeline computes them and within testing laboratory × country, the
strata of the MIC bounds, so that the two can be compared: a call is one cut
of the panel, so its bounds contain the panel's.

    python benchmarks/empirical/ssuis_sensitivity.py <output dir>
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from amr_clonalshare.attribution import clonal_share
from amr_clonalshare.censored import intervals_by_panel
from amr_clonalshare.latent_order import call_bounds
from amr_clonalshare.mic_order import mic_order_share

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "examples/ssuis/data"
SEED = 20260915
SETTINGS = dict(folds=5, repeats=20, n_boot=999, n_perm=999)
ORDER_KEYS = ("kappa_adj", "observed_low", "observed_high", "p_value", "n_readings", "n_groups_repeated",
              "latent_order_lower", "latent_order_upper", "latent_order_lower_limit",
              "latent_order_midpoint", "panel_resolution", "share_end_wells")


def _cutoffs():
    """The cut-off of every agent with a call: the published value where
    there is one, the fitted value otherwise (DATA_PROVENANCE.md)."""
    out = {}
    for r in json.loads((DATA / "ecoff_derived.json").read_text()):
        out[r["agent"]] = r["published_value"] if r["published_value"] is not None else r["fit"]["ecoff"]
    return out


def main(out: Path) -> int:
    started = time.perf_counter()
    out.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(DATA / "metadata.csv", dtype=str, keep_default_na=False).set_index("genome_id")
    source = pd.read_csv(DATA / "source_table_s1.csv", dtype=str, keep_default_na=False).set_index("Strain")
    meta = meta.join(source[["Genome.Set", "Pathogen"]], on="source_strain")
    assert meta["Genome.Set"].ne("").all() and len(meta) == 677
    mic = pd.read_csv(DATA / "mic_panel.csv", dtype=str, keep_default_na=False)
    calls = pd.read_csv(DATA / "calls_long.csv", dtype=str, keep_default_na=False)
    called = sorted(calls.antibiotic.unique())
    cutoffs = _cutoffs()
    lineage = meta.baps_cluster.to_numpy(dtype=object)
    lab_country = None
    order_rows, call_rows = [], []
    for agent, frame in mic.groupby("antibiotic", sort=True):
        frame = frame.set_index("genome_id").reindex(meta.index)
        value = frame.measurement.astype(float).to_numpy()
        lab = frame.testing_laboratory.astype(str).to_numpy()
        lab_country = (frame.testing_laboratory.astype(str) + " / " + meta.isolation_country.astype(str)).to_numpy()
        settings = {"baseline": (lab, lab_country, np.ones(len(meta), bool)),
                    "collection": (lab, (meta["Genome.Set"] + " / " + meta["Pathogen"]).to_numpy(),
                                   np.ones(len(meta), bool)),
                    "without_rcuk_uk": (lab, lab_country, meta["Genome.Set"].ne("RCUK_UK").to_numpy())}
        if agent == "marbofloxacin":
            third = np.where(meta["Genome.Set"].eq("Lola_UK").to_numpy(), lab + " (Lola_UK panel)", lab)
            settings["third_panel"] = (third, lab_country, np.ones(len(meta), bool))
        for name, (panel, strata, keep) in settings.items():
            print("MIC", agent, name, flush=True)
            lo, hi = intervals_by_panel(value[keep], panel[keep], operators=frame.measurement_sign.to_numpy()[keep])
            r = mic_order_share(lo, hi, lineage[keep], strata=strata[keep], seed=SEED, **SETTINGS)
            order_rows.append(dict(agent=agent, setting=name, n=int(keep.sum()),
                                   n_strata=len(r["strata"]), **{k: r[k] for k in ORDER_KEYS}))
        if agent not in called:
            continue
        cut = float(cutoffs[agent])
        # the same dilution is recorded at different precision (0.06 and
        # 0.0625), so readings and cut-offs are compared as log2 dilutions
        well = np.rint(np.log2(value))
        recorded = calls[calls.antibiotic == agent].set_index("genome_id").reindex(meta.index).call
        for step in (-1, 0, 1):
            c = cut * 2. ** step
            y = (well > np.rint(np.log2(c))).astype(float)
            if step == 0:
                assert (y == recorded.map({"1": 1., "0": 0.}).to_numpy()).all(), agent
            print("call", agent, c, flush=True)
            s = clonal_share(y, lineage, seed=SEED, **SETTINGS).as_dict()
            pooled = call_bounds(y, lineage, seed=SEED)
            within = call_bounds(y, lineage, strata=lab_country, seed=SEED)
            call_rows.append(dict(agent=agent, cutoff=c, dilutions_from_cutoff=step, prevalence=float(y.mean()),
                                  share=s["kappa_adj"], observed_low=s["observed_low"],
                                  observed_high=s["observed_high"], p_value=s["p_value"],
                                  pooled_lower=pooled["latent_order_lower"], pooled_upper=pooled["latent_order_upper"],
                                  pooled_limit=pooled["latent_order_lower_limit"],
                                  strata_lower=within["latent_order_lower"], strata_upper=within["latent_order_upper"],
                                  strata_limit=within["latent_order_lower_limit"]))
    pd.DataFrame(order_rows).to_csv(out / "ssuis_sensitivity_order.csv", index=False)
    pd.DataFrame(call_rows).to_csv(out / "ssuis_sensitivity_calls.csv", index=False)
    files = [DATA / "metadata.csv", DATA / "mic_panel.csv", DATA / "calls_long.csv", DATA / "ecoff_derived.json",
             DATA / "source_table_s1.csv", Path(__file__)] + sorted((ROOT / "src/amr_clonalshare").glob("*.py"))
    receipt = dict(seed=SEED, settings=SETTINGS, elapsed_seconds=time.perf_counter() - started,
                   source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in files})
    (out / "receipt.json").write_text(json.dumps(receipt, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1])))
