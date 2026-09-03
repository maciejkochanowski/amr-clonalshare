"""Shared fixtures. Cohorts are small so the fast lane stays under a minute."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import yaml

from amr_clonalshare.config import from_dict
from amr_clonalshare.synthetic import synth_cluster_archetypes, synth_lineage_cohort


def _write_cohort(tmp_path, amr, vir, meta=None, kind_amr="wide_binary"):
    d = tmp_path / "data"
    d.mkdir(exist_ok=True)
    amr.rename_axis("isolate").to_csv(d / "amr.csv")
    vir.rename_axis("isolate").to_csv(d / "vir.csv")
    if meta is not None:
        meta.rename_axis("isolate").to_csv(d / "metadata.csv")
    raw = {
        "dataset": {"name": "t", "strain_id_column": "isolate", "data_dir": "data"},
        "files": {"amr": {"path": "amr.csv", "kind": kind_amr},
                  "vir": {"path": "vir.csv", "kind": "wide_binary"}},
        "trait_cluster": {"layers": ["amr", "vir"], "k_range": [2, 3, 4],
                          "prevalence_gate": {"lo": 0.0, "hi": 1.0, "min_count": 2}},
        "snf": {"K": 10, "T": 10},
        "tva": {"n_splits": 7},
        "influence": {"n_perm": 0},
    }
    if meta is not None:
        raw["dataset"].update({"metadata": "metadata.csv", "lineage_column": "lineage"})
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(yaml.safe_dump(raw, sort_keys=False))
    return from_dict(raw, config_path=cfg_path).validate()


@pytest.fixture
def planted_cfg(tmp_path):
    # Sized so that valid (correlation-blocked) feature splitting still has
    # units to work with. The previous 10 + 14 columns collapsed to six split
    # units, which is enough to detect structure but not enough for any feature
    # to be held out often enough to be testable - a property of a 24-column
    # fixture, not of the method (the shipped planted example, with 60 columns,
    # validates 48 defining features).
    amr, vir, labels = synth_cluster_archetypes(n=100, p_amr=24, p_vir=30,
                                                k_true=3, overlap=0.05, seed=3)
    meta = pd.DataFrame({"truth": labels,
                         "lineage": [f"L{i % 20}" for i in range(len(labels))]},
                        index=amr.index)
    cfg = _write_cohort(tmp_path, amr, vir, meta)
    return cfg, labels


@pytest.fixture
def null_cfg(tmp_path):
    amr, vir, labels = synth_cluster_archetypes(n=60, p_amr=10, p_vir=14,
                                                k_true=1, overlap=0.25, seed=4)
    return _write_cohort(tmp_path, amr, vir)


@pytest.fixture
def clonal_cfg(tmp_path):
    amr, vir, meta = synth_lineage_cohort(n=80, n_lineages=10, p_amr=10,
                                          p_vir=12, seed=5)
    return _write_cohort(tmp_path, amr, vir, meta)


@pytest.fixture
def rng():
    return np.random.default_rng(0)
