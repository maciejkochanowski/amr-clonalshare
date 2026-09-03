"""Lineage concentration of a trait, and which way it departs from chance.

The index is an inverse Simpson number, so the cases that pin it down are the
ones where the answer is arithmetically forced: a trait confined to one lineage
must give one effective carrier, a trait spread evenly over L equally sized
lineages must give L, and a trait distributed exactly in proportion to lineage
abundance must show zero departure however uneven the lineages are.

The departure statistic is a magnitude, so the direction has to be pinned
separately. The case that matters is a trait *avoiding* the dominant lineage:
it produces a large, significant departure that a reader would misread as
clonal concentration unless the direction is reported with it.
"""
from __future__ import annotations

import numpy as np
import pytest

from amr_clonalshare.clonality import trait_concentration


def _cohort(spec):
    """Build (y, lineage) from {lineage: (n, n_positive)}."""
    y, lin = [], []
    for name, (n, k) in spec.items():
        y.extend([1] * k + [0] * (n - k))
        lin.extend([name] * n)
    return np.array(y, dtype=int), np.array(lin, dtype=object)


def test_a_trait_confined_to_one_lineage_has_one_effective_carrier():
    y, lin = _cohort({"A": (20, 20), "B": (20, 0), "C": (20, 0)})
    r = trait_concentration(y, lin, n_perm=200, rng=np.random.default_rng(0))
    assert r["status"] == "ok"
    assert r["effective_number_of_lineages"] == pytest.approx(1.0)
    assert r["n_lineages_carrying"] == 1
    assert r["herfindahl_index"] == pytest.approx(1.0)


def test_a_trait_spread_evenly_gives_the_number_of_lineages():
    y, lin = _cohort({"A": (20, 10), "B": (20, 10), "C": (20, 10),
                      "D": (20, 10)})
    r = trait_concentration(y, lin, n_perm=200, rng=np.random.default_rng(1))
    assert r["effective_number_of_lineages"] == pytest.approx(4.0)
    assert r["herfindahl_index"] == pytest.approx(0.25)


def test_proportional_spread_shows_no_excess_however_uneven_the_lineages():
    # lineage sizes 80/16/4, every lineage at the same rate: the resistant pool
    # mirrors the population exactly, so the excess is zero by construction
    # even though the concentration index itself is far from even.
    y, lin = _cohort({"BIG": (80, 40), "MID": (16, 8), "SMALL": (4, 2)})
    r = trait_concentration(y, lin, n_perm=200, rng=np.random.default_rng(2))
    assert r["departure_from_proportional_bits"] == pytest.approx(0.0, abs=1e-12)
    assert r["effective_number_of_lineages"] < 3.0


def test_excess_is_positive_when_a_small_lineage_carries_the_trait():
    y, lin = _cohort({"BIG": (90, 0), "SMALL": (10, 10)})
    r = trait_concentration(y, lin, n_perm=500, rng=np.random.default_rng(3))
    assert r["departure_from_proportional_bits"] > 1.0
    assert r["effective_number_of_lineages"] == pytest.approx(1.0)
    assert r["p_value"] < 0.05


def test_a_trait_scattered_at_random_is_not_flagged():
    rng = np.random.default_rng(4)
    lin = np.array([f"L{i % 20:02d}" for i in range(400)], dtype=object)
    y = (rng.random(400) < 0.3).astype(int)
    r = trait_concentration(y, lin, n_perm=500, rng=np.random.default_rng(5))
    assert r["p_value"] > 0.05
    assert abs(r["z"]) < 3.0


def test_permutation_null_preserves_prevalence():
    y, lin = _cohort({"A": (30, 15), "B": (30, 5), "C": (30, 25)})
    r = trait_concentration(y, lin, n_perm=300, rng=np.random.default_rng(6))
    # the null shuffles carriage across isolates, so the number carrying the
    # trait is fixed and only its lineage placement varies
    assert r["n_carriers"] == 45
    assert r["null_mean_effective_number"] > r["effective_number_of_lineages"]


def test_concentration_refuses_a_trait_nobody_carries():
    y, lin = _cohort({"A": (10, 0), "B": (10, 0)})
    r = trait_concentration(y, lin, n_perm=50, rng=np.random.default_rng(7))
    assert r["status"] == "skipped"
    assert "carries" in r["reason"]


def test_concentration_refuses_when_no_lineage_label_survives():
    r = trait_concentration(np.array([1, 0]),
                            np.array([None, np.nan], dtype=object))
    assert r["status"] == "skipped"


def test_missing_lineage_labels_are_dropped_and_counted():
    y = np.array([1, 1, 1, 0])
    lin = np.array(["A", "A", None, "B"], dtype=object)
    r = trait_concentration(y, lin, n_perm=50, rng=np.random.default_rng(8))
    assert r["n_dropped_missing_lineage"] == 1
    assert r["n_carriers"] == 2


# --- direction, which the departure magnitude cannot carry ------------------

def test_direction_is_concentrated_when_carriage_piles_into_one_lineage():
    # Ten equally sized lineages, all carriage in one. Chance would spread the
    # carriers over all ten, so the observed single carrier is unambiguous.
    spec = {f"L{i}": (20, 20 if i == 0 else 0) for i in range(10)}
    y, lin = _cohort(spec)
    r = trait_concentration(y, lin, n_perm=500, rng=np.random.default_rng(20))
    assert r["direction"] == "concentrated"
    assert r["effective_number_of_lineages"] < r["null_mean_effective_number"]


def test_direction_has_no_resolution_in_a_cohort_of_two_lineages():
    # One lineage holds 90% of the cohort and none of the carriage. The
    # departure fires, because carriage does not track abundance - but the null
    # is itself concentrated, so the *number* of carrying lineages is
    # indistinguishable from chance. The two fields answer different questions
    # (which lineages carry, versus how many) and are allowed to disagree.
    y, lin = _cohort({"BIG": (90, 0), "SMALL": (10, 10)})
    r = trait_concentration(y, lin, n_perm=500, rng=np.random.default_rng(23))
    assert r["p_value"] < 0.05
    assert r["direction"] == "proportional"


def test_direction_is_dispersed_when_the_trait_avoids_the_dominant_lineage():
    # The shipped S. suis panel behaves this way: the dominant sequence type is
    # nearly free of beta-lactam carriage, so the carrier pool spreads over more
    # lineages than chance allows. The departure is large and significant, and
    # calling that "concentrated" would invert the finding.
    spec = {"BIG": (90, 0)}
    spec.update({f"S{i}": (5, 5) for i in range(10)})
    y, lin = _cohort(spec)
    r = trait_concentration(y, lin, n_perm=500, rng=np.random.default_rng(21))
    assert r["direction"] == "dispersed"
    assert r["effective_number_of_lineages"] > r["null_mean_effective_number"]
    assert r["p_value"] < 0.05          # departure fires...
    assert r["direction"] != "concentrated"   # ...the other way


def test_direction_is_proportional_when_carriage_tracks_abundance():
    y, lin = _cohort({"BIG": (80, 40), "MID": (16, 8), "SMALL": (4, 2)})
    r = trait_concentration(y, lin, n_perm=500, rng=np.random.default_rng(22))
    assert r["direction"] == "proportional"


# --- Monte Carlo resolution --------------------------------------------------

def test_the_permutation_floor_is_reported_with_the_p_value():
    # A panel run at a small budget returns 1/(n_perm+1) for every trait that
    # clears it, which reads as many strong results and is one statement about
    # the budget. The floor is therefore reported next to the value.
    spec = {f"L{i}": (20, 20 if i == 0 else 0) for i in range(10)}
    y, lin = _cohort(spec)
    r = trait_concentration(y, lin, n_perm=100, rng=np.random.default_rng(30))
    assert r["p_value"] == pytest.approx(r["p_value_floor"])
    assert r["p_value_floor"] == pytest.approx(1 / 101)
    assert r["n_perm"] == 100


def test_the_tail_interval_is_exact_and_brackets_the_estimate():
    spec = {f"L{i}": (20, 20 if i == 0 else 0) for i in range(10)}
    y, lin = _cohort(spec)
    r = trait_concentration(y, lin, n_perm=500, rng=np.random.default_rng(31))
    interval = r["tail_probability_ci95"]
    assert interval["method"] == "Clopper-Pearson exact binomial"
    assert r["permutation_exceedances"] == 0
    assert interval["low"] == 0.0
    assert 0.0 < interval["high"] < 0.01
    assert r["resolved_at_alpha_0_05"] is True


def test_a_budget_too_small_to_resolve_the_decision_says_so():
    # Ten permutations cannot put the tail below 0.05 whatever they return: the
    # exact interval on 0 of 10 reaches 0.31.
    spec = {f"L{i}": (20, 20 if i == 0 else 0) for i in range(10)}
    y, lin = _cohort(spec)
    r = trait_concentration(y, lin, n_perm=10, rng=np.random.default_rng(32))
    assert r["permutation_exceedances"] == 0
    assert r["tail_probability_ci95"]["high"] > 0.05
    assert r["resolved_at_alpha_0_05"] is False


def test_a_trait_scattered_at_random_is_unresolved_rather_than_negative():
    rng = np.random.default_rng(33)
    lin = np.array([f"L{i % 20:02d}" for i in range(400)], dtype=object)
    y = (rng.random(400) < 0.3).astype(int)
    r = trait_concentration(y, lin, n_perm=2000, rng=np.random.default_rng(34))
    assert r["p_value"] > 0.05
    interval = r["tail_probability_ci95"]
    assert interval["low"] < interval["high"]
    assert interval["low"] > 0.05 or not r["resolved_at_alpha_0_05"]
