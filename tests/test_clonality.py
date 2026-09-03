"""Composition-versus-acquisition decomposition and lineage-resolved prevalence.

The two estimators are exact algebraic identities plus a resampling interval, so
most of what has to be tested is that the identity holds, that the two extreme
designs (mix-only, rate-only) are separated cleanly, and that the conventions
for unshared lineages and missing labels are the declared ones rather than
whatever the arithmetic happened to do.
"""
from __future__ import annotations

import numpy as np
import pytest

from amr_clonalshare.clonality import (
    decompose_panel,
    decompose_prevalence_difference,
    lineage_resolved_prevalence,
)


def _cohort(spec):
    """Build (y, lineage) from {lineage: (n, n_positive)}."""
    y, lin = [], []
    for name, (n, k) in spec.items():
        y.extend([1] * k + [0] * (n - k))
        lin.extend([name] * n)
    return np.array(y, dtype=int), np.array(lin, dtype=object)


# --- the identity -----------------------------------------------------------

def test_components_sum_to_the_observed_difference():
    rng = np.random.default_rng(0)
    lin_a = rng.choice(list("ABCDE"), 300)
    lin_b = rng.choice(list("CDEFG"), 200)
    ya = rng.integers(0, 2, 300)
    yb = rng.integers(0, 2, 200)
    r = decompose_prevalence_difference(ya, lin_a, yb, lin_b, n_boot=50,
                                        rng=np.random.default_rng(1))
    assert r["status"] == "ok"
    assert r["composition"] + r["within_lineage"] == pytest.approx(
        r["difference"], abs=1e-12)
    assert abs(r["identity_residual"]) < 1e-12


def test_identity_holds_when_the_two_groups_share_no_lineage():
    ya, lin_a = _cohort({"A": (40, 30), "B": (60, 12)})
    yb, lin_b = _cohort({"X": (50, 5), "Y": (50, 40)})
    r = decompose_prevalence_difference(ya, lin_a, yb, lin_b, n_boot=25,
                                        rng=np.random.default_rng(2))
    assert r["composition"] + r["within_lineage"] == pytest.approx(
        r["difference"], abs=1e-12)
    assert r["n_lineages_shared"] == 0
    assert r["turnover_share"] == pytest.approx(1.0)


# --- the two extreme designs ------------------------------------------------

def test_a_pure_mix_difference_loads_entirely_on_composition():
    # identical within-lineage rates, different lineage shares
    ya, lin_a = _cohort({"X": (80, 40), "Y": (20, 18)})   # 0.5 / 0.9
    yb, lin_b = _cohort({"X": (20, 10), "Y": (80, 72)})   # 0.5 / 0.9
    r = decompose_prevalence_difference(ya, lin_a, yb, lin_b, n_boot=25,
                                        rng=np.random.default_rng(3))
    assert r["within_lineage"] == pytest.approx(0.0, abs=1e-12)
    assert r["composition"] == pytest.approx(r["difference"], abs=1e-12)


def test_a_pure_rate_difference_loads_entirely_on_within_lineage():
    # identical lineage shares, different within-lineage rates
    ya, lin_a = _cohort({"X": (50, 10), "Y": (50, 20)})   # 0.2 / 0.4
    yb, lin_b = _cohort({"X": (50, 30), "Y": (50, 40)})   # 0.6 / 0.8
    r = decompose_prevalence_difference(ya, lin_a, yb, lin_b, n_boot=25,
                                        rng=np.random.default_rng(4))
    assert r["composition"] == pytest.approx(0.0, abs=1e-12)
    assert r["within_lineage"] == pytest.approx(r["difference"], abs=1e-12)


def test_a_lineage_present_in_one_group_only_is_charged_to_composition():
    # X is shared and identical; Z exists only in A. The declared convention
    # sends the whole of Z's contribution to composition, because "that lineage
    # is not there at all" is a statement about the mix, not about a rate.
    ya, lin_a = _cohort({"X": (50, 25), "Z": (50, 50)})
    yb, lin_b = _cohort({"X": (50, 25)})
    r = decompose_prevalence_difference(ya, lin_a, yb, lin_b, n_boot=25,
                                        rng=np.random.default_rng(5))
    assert r["within_lineage"] == pytest.approx(0.0, abs=1e-12)
    assert r["n_lineages_only_a"] == 1
    # Z moves 0.50 of the mass, the shared X moves 0.25 of it the other way.
    assert r["turnover_share"] == pytest.approx(2 / 3)


# --- coverage, missingness, refusal -----------------------------------------

def test_missing_lineage_labels_are_dropped_and_counted():
    ya = np.array([1, 0, 1, 0, 1, 1])
    lin_a = np.array(["X", "X", None, "Y", "Y", np.nan], dtype=object)
    yb = np.array([0, 1, 0, 1])
    lin_b = np.array(["X", "X", "Y", "Y"], dtype=object)
    r = decompose_prevalence_difference(ya, lin_a, yb, lin_b, n_boot=25,
                                        rng=np.random.default_rng(6))
    assert r["status"] == "ok"
    assert r["n_dropped_missing_lineage"] == 2
    assert r["n_a"] == 4


def test_decomposition_refuses_when_no_lineage_label_survives():
    ya = np.array([1, 0]); yb = np.array([0, 1])
    miss = np.array([None, np.nan], dtype=object)
    r = decompose_prevalence_difference(ya, miss, yb, miss)
    assert r["status"] == "skipped"
    assert "lineage" in r["reason"]


def test_decomposition_refuses_an_empty_group():
    ya, lin_a = _cohort({"X": (10, 5)})
    r = decompose_prevalence_difference(ya, lin_a, np.array([]),
                                        np.array([], dtype=object))
    assert r["status"] == "skipped"


def test_bootstrap_interval_brackets_the_point_estimate():
    rng = np.random.default_rng(7)
    lin_a = rng.choice(list("ABCD"), 400)
    lin_b = rng.choice(list("ABCD"), 400)
    ya = (rng.random(400) < 0.6).astype(int)
    yb = (rng.random(400) < 0.3).astype(int)
    r = decompose_prevalence_difference(ya, lin_a, yb, lin_b, n_boot=400,
                                        rng=np.random.default_rng(8))
    lo, hi = r["within_lineage_ci95"]
    assert lo <= r["within_lineage"] <= hi
    lo, hi = r["composition_ci95"]
    assert lo <= r["composition"] <= hi


# --- lineage-resolved prevalence --------------------------------------------

def test_per_lineage_equals_per_isolate_when_every_lineage_is_a_singleton():
    y = np.array([1, 0, 1, 1, 0])
    lin = np.array(list("ABCDE"), dtype=object)
    r = lineage_resolved_prevalence(y, lin, n_boot=50,
                                    rng=np.random.default_rng(9))
    assert r["status"] == "ok"
    assert r["prevalence_per_lineage"] == pytest.approx(
        r["prevalence_per_isolate"])


def test_a_dominant_susceptible_clone_pulls_the_two_estimands_apart():
    # 90 isolates of one fully susceptible lineage, 10 singletons all resistant
    spec = {"BIG": (90, 0)}
    spec.update({f"S{i}": (1, 1) for i in range(10)})
    y, lin = _cohort(spec)
    r = lineage_resolved_prevalence(y, lin, n_boot=200,
                                    rng=np.random.default_rng(10))
    assert r["prevalence_per_isolate"] == pytest.approx(10 / 100)
    assert r["prevalence_per_lineage"] == pytest.approx(10 / 11)
    assert r["difference_per_lineage_minus_per_isolate"] > 0.5
    assert r["largest_lineage_share"] == pytest.approx(0.9)


def test_fraction_of_lineages_carrying_the_trait_at_all():
    y, lin = _cohort({"X": (10, 3), "Y": (10, 0), "Z": (10, 10)})
    r = lineage_resolved_prevalence(y, lin, n_boot=50,
                                    rng=np.random.default_rng(11))
    assert r["fraction_of_lineages_with_any"] == pytest.approx(2 / 3)
    assert r["prevalence_per_lineage"] == pytest.approx((0.3 + 0.0 + 1.0) / 3)


def test_prevalence_interval_resamples_lineages_not_isolates():
    # One lineage carries everything. Resampling lineages must be able to
    # exclude it, so the interval has to reach down to zero; an interval built
    # by resampling isolates within lineages could not.
    y, lin = _cohort({"A": (10, 10), "B": (10, 0), "C": (10, 0)})
    r = lineage_resolved_prevalence(y, lin, n_boot=1000,
                                    rng=np.random.default_rng(12))
    lo, hi = r["prevalence_per_lineage_ci95"]
    assert lo == pytest.approx(0.0, abs=1e-9)
    assert hi > r["prevalence_per_lineage"]


def test_prevalence_refuses_when_no_lineage_label_survives():
    r = lineage_resolved_prevalence(np.array([1, 0]),
                                    np.array([None, np.nan], dtype=object))
    assert r["status"] == "skipped"


# --- the estimability gate ---------------------------------------------------

def test_within_lineage_is_not_estimable_when_the_collections_barely_overlap():
    # One shared lineage of three, holding a fifth of each collection. The
    # component is still reported, because a number a reader can see is better
    # than a silence, but it is marked not estimable and the reason is named.
    ya, lin_a = _cohort({"S": (20, 10), "A1": (40, 30), "A2": (40, 5)})
    yb, lin_b = _cohort({"S": (20, 12), "B1": (40, 8), "B2": (40, 32)})
    r = decompose_prevalence_difference(ya, lin_a, yb, lin_b, n_boot=200,
                                        rng=np.random.default_rng(30))
    assert r["shared_support_isolate_share"] == pytest.approx(0.2)
    assert r["within_lineage_estimable"] is False
    assert r["shared_support_margin"] < 0
    assert any("support" in reason for reason in r["not_estimable_because"])
    assert r["within_lineage"] is not None


def test_a_well_overlapped_contrast_is_estimable():
    ya, lin_a = _cohort({"X": (100, 40), "Y": (100, 60)})
    yb, lin_b = _cohort({"X": (100, 30), "Y": (100, 50)})
    r = decompose_prevalence_difference(ya, lin_a, yb, lin_b, n_boot=200,
                                        rng=np.random.default_rng(31))
    assert r["shared_support_isolate_share"] == pytest.approx(1.0)
    assert r["within_lineage_estimable"] is True
    assert r["not_estimable_because"] == []


def test_the_gate_threshold_is_a_parameter_and_reports_its_margin():
    ya, lin_a = _cohort({"S": (50, 20), "A": (50, 40)})
    yb, lin_b = _cohort({"S": (50, 25), "B": (50, 10)})
    strict = decompose_prevalence_difference(ya, lin_a, yb, lin_b, n_boot=100,
                                             rng=np.random.default_rng(32),
                                             min_shared_support=0.9)
    loose = decompose_prevalence_difference(ya, lin_a, yb, lin_b, n_boot=100,
                                            rng=np.random.default_rng(32),
                                            min_shared_support=0.1)
    assert strict["within_lineage_estimable"] is False
    assert loose["within_lineage_estimable"] is True
    assert strict["shared_support_margin"] < 0 < loose["shared_support_margin"]


# --- differential and informative label missingness --------------------------

def _labelled(n_labelled, k_labelled, n_missing, k_missing, tag):
    """A collection where some isolates carry no lineage label."""
    y = [1] * k_labelled + [0] * (n_labelled - k_labelled)
    lin = [f"{tag}{i % 4}" for i in range(n_labelled)]
    y += [1] * k_missing + [0] * (n_missing - k_missing)
    lin += [None] * n_missing
    return np.array(y, dtype=float), np.array(lin, dtype=object)


def test_differential_missingness_alone_does_not_fire_the_gate():
    # Coverage differs sharply, but labelling is unrelated to the trait in both
    # collections, so the difference between them is not distorted.
    ya, lin_a = _labelled(40, 20, 60, 30, "L")   # 40 % labelled, rate 0.5 both
    yb, lin_b = _labelled(90, 45, 10, 5, "L")    # 90 % labelled, rate 0.5 both
    r = decompose_prevalence_difference(ya, lin_a, yb, lin_b, n_boot=100,
                                        rng=np.random.default_rng(33))
    a = r["lineage_label_availability"]
    assert a["label_coverage_differs"] is True
    assert a["missingness_informative"] is False
    assert a["labels_representative"] is True


def test_informative_missingness_alone_does_not_fire_the_gate():
    # The unlabelled isolates are far more often positive, but both collections
    # are labelled to the same degree, so the selection applies equally.
    ya, lin_a = _labelled(50, 10, 50, 45, "L")
    yb, lin_b = _labelled(50, 12, 50, 44, "L")
    r = decompose_prevalence_difference(ya, lin_a, yb, lin_b, n_boot=100,
                                        rng=np.random.default_rng(34))
    a = r["lineage_label_availability"]
    assert a["label_coverage_differs"] is False
    assert a["missingness_informative"] is True
    assert a["labels_representative"] is True


def test_the_gate_fires_when_missingness_is_both_differential_and_informative():
    # This is the shipped S. suis period contrast in miniature: the later
    # collection is mostly unlabelled and its unlabelled isolates are the
    # resistant ones, so the labelled subsets are not comparable.
    ya, lin_a = _labelled(30, 6, 70, 63, "L")
    yb, lin_b = _labelled(90, 20, 10, 9, "L")
    r = decompose_prevalence_difference(ya, lin_a, yb, lin_b, n_boot=100,
                                        rng=np.random.default_rng(35))
    a = r["lineage_label_availability"]
    assert a["label_coverage_differs"] is True
    assert a["missingness_informative"] is True
    assert a["labels_representative"] is False
    assert r["within_lineage_estimable"] is False
    assert any("differentially missing" in reason
               for reason in r["not_estimable_because"])


def test_a_non_binary_trait_is_reported_as_not_assessable():
    y = np.array([0.3, 1.7, 2.2, 0.9])
    lin = np.array(["A", "A", None, "B"], dtype=object)
    r = decompose_prevalence_difference(y, lin, y, lin, n_boot=50,
                                        rng=np.random.default_rng(36))
    assert r["lineage_label_availability"]["assessable"] is False


# --- bootstrap p-values ------------------------------------------------------

def test_a_component_that_is_zero_by_construction_is_not_a_discovery():
    # The within-lineage rates are identical by construction, so the observed
    # component is exactly zero. Resampling perturbs the rates, so the p-value
    # is near one rather than exactly one, and that is the correct behaviour:
    # the bootstrap is asked how far the component could have moved, not
    # whether the arithmetic held.
    ya, lin_a = _cohort({"X": (80, 40), "Y": (20, 18)})
    yb, lin_b = _cohort({"X": (20, 10), "Y": (80, 72)})
    r = decompose_prevalence_difference(ya, lin_a, yb, lin_b, n_boot=500,
                                        rng=np.random.default_rng(37))
    assert r["within_lineage"] == pytest.approx(0.0, abs=1e-12)
    assert r["within_lineage_p"] > 0.9
    lo, hi = r["within_lineage_ci95"]
    assert lo <= 0.0 <= hi


def test_the_p_value_agrees_with_the_interval_it_inverts():
    rng = np.random.default_rng(38)
    lin_a = rng.choice(list("ABCD"), 400)
    lin_b = rng.choice(list("ABCD"), 400)
    ya = (rng.random(400) < 0.65).astype(int)
    yb = (rng.random(400) < 0.30).astype(int)
    r = decompose_prevalence_difference(ya, lin_a, yb, lin_b, n_boot=1000,
                                        rng=np.random.default_rng(39))
    lo, hi = r["within_lineage_ci95"]
    assert (lo * hi > 0) == (r["within_lineage_p"] < 0.05)
    assert r["p_value_floor"] == pytest.approx(2 / 1001)


# --- the panel ---------------------------------------------------------------

def test_a_panel_of_identical_columns_is_worth_one_independent_agent():
    import pandas as pd
    rng = np.random.default_rng(40)
    lin_a = rng.choice(list("ABCDE"), 200)
    lin_b = rng.choice(list("ABCDE"), 200)
    base_a = (rng.random(200) < 0.5).astype(float)
    base_b = (rng.random(200) < 0.3).astype(float)
    frame_a = pd.DataFrame({f"agent{i}": base_a for i in range(5)})
    frame_b = pd.DataFrame({f"agent{i}": base_b for i in range(5)})
    out = decompose_panel(frame_a, lin_a, frame_b, lin_b, n_boot=200,
                          rng=np.random.default_rng(41))
    assert out["family"]["n_agents"] == 5
    assert out["family"]["effective_independent_agents"] == pytest.approx(1.0)


def test_the_panel_controls_the_false_discovery_rate():
    import pandas as pd
    rng = np.random.default_rng(42)
    lin_a = rng.choice(list("ABCDE"), 300)
    lin_b = rng.choice(list("ABCDE"), 300)
    frame_a = pd.DataFrame({f"agent{i}": (rng.random(300) < 0.4).astype(float)
                            for i in range(10)})
    frame_b = pd.DataFrame({f"agent{i}": (rng.random(300) < 0.4).astype(float)
                            for i in range(10)})
    out = decompose_panel(frame_a, lin_a, frame_b, lin_b, n_boot=400,
                          rng=np.random.default_rng(43), q=0.05)
    family = out["family"]
    # nothing is real here, so control must not exceed the nominal count and
    # the adjusted values must dominate the raw ones
    assert family["n_within_lineage_discoveries"] <= family["n_within_lineage_nominal"]
    for record in out["per_agent"].values():
        assert record["within_lineage_q"] >= record["within_lineage_p"] - 1e-12
        assert 0.0 <= record["composition_q"] <= 1.0


def test_a_panel_refuses_mismatched_columns():
    import pandas as pd
    frame_a = pd.DataFrame({"x": [0.0, 1.0]})
    frame_b = pd.DataFrame({"y": [0.0, 1.0]})
    with pytest.raises(ValueError, match="same agent columns"):
        decompose_panel(frame_a, ["A", "B"], frame_b, ["A", "B"], n_boot=10)


def test_a_non_finite_trait_is_not_reported_as_a_missing_lineage_label():
    """The analysis mask is "label present **and** trait finite", and its
    complement was published as ``n_dropped_missing_lineage`` and as a fall in
    ``label_coverage_a`` -- claims SWX-057 and SWX-071. A NaN well is not an
    untyped isolate, and a genuinely untyped isolate must still be counted."""
    y = np.array([1.0, 0.0] * 10)
    y[3] = np.nan
    y[8] = np.nan
    lin = np.array([f"ST{i % 4}" for i in range(20)], dtype=object)

    per_lineage = lineage_resolved_prevalence(y, lin, n_boot=50,
                                              rng=np.random.default_rng(50))
    assert per_lineage["n_dropped_missing_lineage"] == 0
    assert per_lineage["n_dropped_non_finite_trait"] == 2
    assert per_lineage["n_isolates"] == 18

    both = decompose_prevalence_difference(y, lin, y, lin, n_boot=50,
                                           rng=np.random.default_rng(51))
    assert both["n_dropped_missing_lineage"] == 0
    assert both["n_dropped_non_finite_trait"] == 4
    assert both["lineage_label_availability"]["label_coverage_a"] == 1.0

    untyped = lin.copy()
    untyped[5] = None
    mixed = lineage_resolved_prevalence(y, untyped, n_boot=50,
                                        rng=np.random.default_rng(52))
    assert mixed["n_dropped_missing_lineage"] == 1
    assert mixed["n_dropped_non_finite_trait"] == 2
