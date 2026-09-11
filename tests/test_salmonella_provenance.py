"""The shipped poultry-meat *Salmonella* cell is the cell the article reads.

The tables in ``examples/salmonella_poultry/data`` are cut from one pinned
NCBI Pathogen Detection release by ``build_cell.py``; the release itself is
not in the repository, so these tests check the shipped tables against the
counts the article and the cell receipt state, and that both configurations
load and pass the input check.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from amr_clonalshare.config import load_config
from amr_clonalshare.io import load_dataset

ROOT = Path(__file__).resolve().parents[1] / "examples" / "salmonella_poultry"
DATA = ROOT / "data"


@pytest.fixture(scope="module")
def tables():
    meta = pd.read_csv(DATA / "metadata.csv", dtype={"genome_id": str})
    calls = pd.read_csv(DATA / "calls_long.csv", dtype={"genome_id": str})
    receipt = json.loads((DATA / "cell_receipt.json").read_text())
    return meta, calls, receipt


def test_the_cell_is_the_one_the_receipt_and_the_article_describe(tables):
    meta, calls, receipt = tables
    assert receipt["release"] == "PDG000000002.4210"
    assert len(meta) == receipt["n_isolates"] == 7049
    assert meta["genome_id"].is_unique
    assert meta["serovar"].nunique() == receipt["n_serovars"]
    assert meta["pds_cluster"].nunique() == receipt["n_clusters"]
    assert set(calls["genome_id"]) <= set(meta["genome_id"])
    assert set(calls["call"]) == {"susceptible", "non-susceptible"}
    assert calls["antibiotic"].nunique() == receipt["n_antimicrobials"] == 22
    # The nalidixic-acid reading of Section 3.1: 6,915 isolates with a call,
    # 584 non-susceptible, 84 serovars, 496 clusters.
    nal = calls[calls["antibiotic"] == "nalidixic acid"].merge(meta, on="genome_id")
    assert len(nal) == 6915
    assert int((nal["call"] == "non-susceptible").sum()) == 584
    assert nal["serovar"].nunique() == 84
    assert nal["pds_cluster"].nunique() == 496
    assert (meta["period"] == "from_2016").eq(meta["collection_year"] >= 2016).all()


@pytest.mark.parametrize("name", ["config.yaml", "config_cluster.yaml"])
def test_both_configurations_load_and_pass_the_input_check(name):
    cfg = load_config(ROOT / name)
    ds = load_dataset(cfg)
    assert ds.input_qc is not None
    assert ds.input_qc["n_isolates"] == 7049
    groups = ds.input_qc["lineage"]
    assert groups["estimable"] is True
    assert groups["support"] > 0.95
