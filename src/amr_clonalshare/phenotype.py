"""phenotype.py — read categorical susceptibility calls into a trait panel.

Long antimicrobial-susceptibility records, one row per isolate and agent with
a susceptible / intermediate / resistant call, become the isolate-by-agent
matrix of non-wild-type indicators the estimators read. The handling of
"Intermediate" is stated rather than defaulted silently: ``"non_susceptible"``
(the surveillance convention, and the default) groups I with R, ``"drop"``
sets those results aside for that agent, ``"susceptible"`` groups I with S. The choice moves
prevalence, so it is echoed into the record by the loader.
"""
from __future__ import annotations

import pandas as pd

__all__ = ["to_non_susceptible"]

_SIR = {"resistant": 1, "intermediate": None, "susceptible": 0,
        "nonsusceptible": 1, "non-susceptible": 1}


def to_non_susceptible(long_df: pd.DataFrame, *, id_column: str = "Strain_ID",
                       antibiotic_column: str = "antibiotic",
                       call_column: str = "resistant_phenotype",
                       intermediate: str = "non_susceptible") -> pd.DataFrame:
    """Long AST records -> wide isolate x antibiotic matrix of 0/1/NaN.

    ``NaN`` means "not tested", which is the overwhelming majority of cells and
    must stay distinguishable from "tested and susceptible".
    """
    if intermediate not in ("non_susceptible", "drop", "susceptible"):
        raise ValueError(f"unknown intermediate policy {intermediate!r}")
    df = long_df[[id_column, antibiotic_column, call_column]].copy()
    df[call_column] = df[call_column].astype(str).str.strip().str.lower()
    mapped = df[call_column].map(_SIR)
    inter = df[call_column].eq("intermediate")
    if intermediate == "non_susceptible":
        mapped = mapped.where(~inter, 1.0)
    elif intermediate == "susceptible":
        mapped = mapped.where(~inter, 0.0)
    df["_ns"] = mapped
    df = df.dropna(subset=["_ns"])
    # An isolate can carry several records for one drug (repeat testing or
    # several typing methods). Any non-susceptible result wins, which is the
    # conservative reading for surveillance.
    wide = df.pivot_table(index=id_column, columns=antibiotic_column,
                          values="_ns", aggfunc="max")
    return wide
