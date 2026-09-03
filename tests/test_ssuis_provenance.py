"""The shipped S. suis phenotype matrix must be derivable from its raw inputs.

A binary matrix of non-wild-type calls is two decisions away from the
measurements it came from: which cut-off, and which side of it. Both are
recorded, so both can be checked, and a test that rebuilds all 8,801 cells is
worth more than a sentence saying they were built correctly.

The identifier is read as a string throughout. `1307.7510` parsed as a float
loses its trailing zero and joins nothing, which is the trap that once disabled
the lineage diagnostic on this cohort while 458 labels sat in the file.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

DATA = Path(__file__).resolve().parents[1] / "examples" / "ssuis" / "data"
pytestmark = pytest.mark.skipif(not (DATA / "mic_long.csv").is_file(),
                                reason="S. suis example data not present")


@pytest.fixture(scope="module")
def cohort():
    mic = pd.read_csv(DATA / "mic_long.csv", dtype={"genome_id": str})
    calls = (pd.read_csv(DATA / "ribo.csv")
             .merge(pd.read_csv(DATA / "cell.csv"), on="genome_id")
             .set_index("genome_id"))
    cutoffs = {row["agent"]: (row["published_value"] or row["fit"]["ecoff"])
               for row in json.loads(
                   (DATA / "ecoff_derived.json").read_text())}
    return mic, calls, cutoffs


def test_every_call_rebuilds_from_the_raw_mic_and_the_recorded_cutoff(cohort):
    mic, calls, cutoffs = cohort
    wide = mic.assign(key="g" + mic["genome_id"]).pivot(
        index="key", columns="antibiotic", values="measurement")
    wide = wide.reindex(calls.index)
    assert wide[calls.columns].isna().sum().sum() == 0, (
        "a call has no MIC record; the identifier join has failed")
    rebuilt = (wide[calls.columns] > pd.Series(cutoffs)[calls.columns]).astype(int)
    mismatched = (rebuilt != calls).sum()
    assert mismatched.sum() == 0, (
        f"cells that do not rebuild: {mismatched[mismatched > 0].to_dict()}")
    assert calls.size == 8801


def test_the_panel_carries_no_censored_measurement(cohort):
    # Every cut-off comparison below is a strict inequality on an exact value.
    # One censored record would make that rule wrong rather than merely
    # incomplete, so its absence is asserted rather than assumed.
    mic, _, _ = cohort
    signs = mic["measurement_sign"].dropna().unique()
    assert len(signs) == 0, f"censored measurements present: {signs}"


def test_three_of_thirteen_cutoffs_are_published_and_ten_are_derived(cohort):
    _, calls, _ = cohort
    rows = {row["agent"]: row
            for row in json.loads((DATA / "ecoff_derived.json").read_text())}
    published = [a for a in calls.columns if rows[a]["published_value"]]
    derived = [a for a in calls.columns if not rows[a]["published_value"]]
    assert sorted(published) == ["doxycycline", "erythromycin", "tetracycline"]
    assert len(derived) == 10


def test_the_metadata_carries_no_date_coerced_serotype():
    from amr_clonalshare.io import metadata_quality
    metadata = pd.read_csv(DATA / "metadata.csv", dtype=str)
    assert "date_coerced_values" not in metadata_quality(metadata), (
        "run examples/ssuis/repair_metadata.py")


def test_every_isolate_carries_a_population_cluster():
    metadata = pd.read_csv(DATA / "metadata.csv", dtype=str)
    assert metadata["baps_cluster"].notna().all()
    assert metadata["baps_cluster"].nunique() == 30
    # the sequence type is the one that is missing, and informatively so
    assert metadata["mlst"].isna().sum() == 219


def test_the_lineage_join_receipt_records_an_exact_match():
    receipt = json.loads((DATA / "lineage_link_receipt.json").read_text())
    checks = receipt["validation"]
    assert checks["match_rate"] == 1.0
    assert checks["n_matched"] == checks["n_distinct_strains_claimed"] == 677
    # neither of these takes part in the join
    assert checks["mic_agreement"] == 1.0
    assert checks["n_mic_cells_compared"] == 10832
    assert checks["year_agreement"] == 1.0
    assert checks["adjusted_mutual_information_with_mlst"] > 0.5
