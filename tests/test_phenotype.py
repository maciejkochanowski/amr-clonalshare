"""Scoring a partition against measured phenotype, and against what."""
import numpy as np
import pandas as pd
import pytest

from amr_clonalshare.phenotype import phenotype_concordance, to_non_susceptible


def _long(rows):
    return pd.DataFrame(rows, columns=["Strain_ID", "antibiotic",
                                       "resistant_phenotype"])


def test_intermediate_policy_is_stated_not_assumed():
    df = _long([("a", "mero", "Resistant"), ("b", "mero", "Intermediate"),
                ("c", "mero", "Susceptible")])
    ns = to_non_susceptible(df, intermediate="non_susceptible")
    assert ns.loc["b", "mero"] == 1.0
    sus = to_non_susceptible(df, intermediate="susceptible")
    assert sus.loc["b", "mero"] == 0.0
    dropped = to_non_susceptible(df, intermediate="drop")
    assert "b" not in dropped.index
    with pytest.raises(ValueError):
        to_non_susceptible(df, intermediate="whatever")


def test_untested_stays_distinguishable_from_susceptible():
    df = _long([("a", "mero", "Resistant"), ("b", "cipro", "Susceptible")])
    ns = to_non_susceptible(df)
    assert np.isnan(ns.loc["a", "cipro"])          # untested
    assert ns.loc["b", "cipro"] == 0.0             # tested, susceptible


def test_repeat_testing_takes_the_non_susceptible_result():
    df = _long([("a", "mero", "Susceptible"), ("a", "mero", "Resistant")])
    assert to_non_susceptible(df).loc["a", "mero"] == 1.0


def test_a_constant_phenotype_is_reported_not_scored():
    ids = [f"i{i}" for i in range(40)]
    labels = [0] * 20 + [1] * 20
    pheno = pd.DataFrame({"amp": np.ones(40)}, index=ids)   # intrinsic resistance
    out = phenotype_concordance(labels, ids, pheno, min_tested=5)
    rec = out["per_antibiotic"][0]
    assert rec["status"] == "uninformative"
    assert out["n_fdr_significant"] == 0


def test_a_partition_that_loses_to_one_column_is_reported_as_losing():
    """The whole point: association is not the question, beating a column is."""
    rng = np.random.default_rng(0)
    n = 200
    ids = [f"i{i}" for i in range(n)]
    truth = rng.integers(0, 2, n)                  # what actually drives phenotype
    labels = np.where(rng.random(n) < 0.25, 1 - truth, truth)   # noisy copy
    pheno = pd.DataFrame(
        {f"ab{j}": np.where(rng.random(n) < 0.1, 1 - truth, truth).astype(float)
         for j in range(4)}, index=ids)
    out = phenotype_concordance(labels, ids, pheno,
                                competing_rules={"the_truth": truth},
                                min_tested=20)
    assert out["status"] == "ok"
    assert out["n_fdr_significant"] == 4
    ba = out["mean_balanced_accuracy"]
    assert ba["the_truth"] > ba["partition"]
    assert out["partition_beats"]["the_truth"] is False


def test_length_mismatch_is_an_error_not_a_silent_reindex():
    ids = ["a", "b", "c"]
    pheno = pd.DataFrame({"x": [1.0, 0.0, 1.0]}, index=ids)
    with pytest.raises(ValueError):
        phenotype_concordance([0, 1], ids, pheno)
    with pytest.raises(ValueError):
        phenotype_concordance([0, 1, 0], ids, pheno,
                              competing_rules={"bad": [0, 1]})
