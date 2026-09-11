"""Shared fixtures. Cohorts are small so the fast lane stays under a minute."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import yaml

from amr_clonalshare.config import from_dict


def planted_cohort(tmp_path, n_per_lineage=12, n_lineages=8, seed=11, *,
                   intake=True, mic=False):
    """A planted cohort on disk: two agents that follow the lineage, one that
    does not, as a long susceptibility table and a metadata table, and with
    ``mic`` a dilution table on a doubling lattice beside them. Returns the raw
    configuration mapping that reads it."""
    rng = np.random.default_rng(seed)
    ids, lineages = [], []
    for g in range(n_lineages):
        for i in range(n_per_lineage):
            ids.append(f"iso{g:02d}_{i:02d}")
            lineages.append(f"L{g:02d}")
    carrier = {f"L{g:02d}": int(g < n_lineages // 2) for g in range(n_lineages)}
    rows, dilutions = [], []
    wells = [0.25, 0.5, 1, 2, 4, 8, 16, 32]
    for iso, lin in zip(ids, lineages):
        for agent, clonal in (("agentA", True), ("agentB", True),
                              ("agentC", False)):
            p = (0.9 if carrier[lin] else 0.1) if clonal else 0.5
            resistant = rng.random() < p
            call = "non-susceptible" if resistant else "susceptible"
            rows.append({"iid": iso, "antibiotic": agent, "call": call})
            well = (4 + rng.integers(0, 4)) if resistant else rng.integers(0, 4)
            dilutions.append({"iid": iso, "antibiotic": agent,
                              "measurement": wells[well]})
    d = tmp_path / "data"
    d.mkdir(exist_ok=True)
    pd.DataFrame(rows).to_csv(d / "calls.csv", index=False)
    if mic:
        pd.DataFrame(dilutions).to_csv(d / "mic.csv", index=False)
    meta = {"iid": ids, "lineage": lineages}
    if intake:
        meta["intake"] = [f"20{10 + (i % 3)}" for i in range(len(ids))]
    pd.DataFrame(meta).to_csv(d / "meta.csv", index=False)
    raw = {
        "dataset": {
            "name": "planted_share",
            "strain_id_column": "iid",
            "data_dir": "data",
            "metadata": "meta.csv",
            "lineage_column": "lineage",
            "phenotype": "calls.csv",
            "phenotype_id_column": "iid",
            "phenotype_antibiotic_column": "antibiotic",
            "phenotype_call_column": "call",
        },
        "attribution": {"n_boot": 60, "n_perm": 60},
        "evidence": {"folds": 3, "repeats": 1},
    }
    if intake:
        raw["dataset"]["batch_column"] = "intake"
    if mic:
        raw["dataset"].update({"mic": "mic.csv", "mic_id_column": "iid"})
        raw["censored"] = {"n_boot": 20}
    return raw


def write_config(tmp_path, raw):
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(yaml.safe_dump(raw, sort_keys=False))
    return cfg_path, from_dict(raw, config_path=cfg_path).validate()


@pytest.fixture
def share_cfg(tmp_path):
    return write_config(tmp_path, planted_cohort(tmp_path))


@pytest.fixture
def rng():
    return np.random.default_rng(0)
