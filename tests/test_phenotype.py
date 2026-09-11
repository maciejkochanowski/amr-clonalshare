"""Reading categorical susceptibility calls into a trait panel."""
import numpy as np
import pandas as pd
import pytest

from amr_clonalshare.phenotype import to_non_susceptible


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
