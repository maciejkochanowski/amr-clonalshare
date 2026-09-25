"""How much of the answer is the data, and how much is the arithmetic.

Two implementations agreeing to twelve digits says nothing about whether the
quantity they compute is determined by the collection. These tests measure the
other two sources of movement: the Monte Carlo draw, which the seed fixes but
does not remove, and the single record, whose removal should not carry the
conclusion. Where the measurement is a property of the statistical problem
rather than of the code, the package must report it, not hide it.
"""
from __future__ import annotations

import numpy as np
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from amr_clonalshare.attribution import clonal_share
from amr_clonalshare.evalues import e_process
from amr_clonalshare.realised import realised_share

BUDGETS = dict(folds=3, repeats=5, n_boot=60, n_perm=99)


def _cohort(seed=5, lineages=9, per=6, prevalence=0.4):
    rng = np.random.default_rng(seed)
    lineage = np.repeat([f"L{i:02d}" for i in range(lineages)], per)
    return (rng.random(lineages * per) < prevalence).astype(float), lineage


@pytest.mark.slow
def test_the_draw_moves_the_estimate_far_less_than_the_interval_does():
    """The reported interval must dominate the Monte Carlo error.

    Measured on this fixture: a standard deviation of about 0.04 across twelve
    seeds against an interval 0.60 wide, that is one part in fourteen. The
    test allows a fifth, so it passes comfortably today and fails if a change
    to the resampling makes the draw a comparable source of movement -- at
    which point the interval no longer describes the uncertainty of the
    reported number.
    """
    y, lineage = _cohort()
    estimates = [clonal_share(y, lineage, seed=s, **BUDGETS).kappa_adj
                 for s in range(12)]
    reference = clonal_share(y, lineage, seed=0, **BUDGETS)
    width = reference.ci_high - reference.ci_low
    assert np.std(estimates) < width / 5.0, (
        f"Monte Carlo sd {np.std(estimates):.4f} against an interval {width:.4f} wide")


@pytest.mark.slow
def test_no_single_isolate_carries_the_estimate():
    """Leave one record out, one at a time.

    Measured on this fixture: the largest move is about 0.11, against an
    interval 0.60 wide. A collection in which one isolate moves the estimate
    by more than its whole interval is reporting that isolate, not the
    collection.
    """
    y, lineage = _cohort()
    reference = clonal_share(y, lineage, seed=0, **BUDGETS)
    width = reference.ci_high - reference.ci_low
    moves = []
    for i in range(y.size):
        keep = np.ones(y.size, dtype=bool)
        keep[i] = False
        moves.append(abs(clonal_share(y[keep], lineage[keep], seed=0, **BUDGETS).kappa_adj
                         - reference.kappa_adj))
    assert max(moves) < width, (
        f"dropping one isolate moved the estimate by {max(moves):.4f}, "
        f"more than the reported interval of {width:.4f}")


@pytest.mark.slow
def test_one_reading_read_one_dilution_out_does_not_overturn_the_reading():
    """Measurement error of one dilution on one isolate.

    A dilution panel is read to the nearest doubling, so one reading in the
    wrong well is the ordinary error of the method. The variance ratio may
    move; it may not move by more than its own interval.
    """
    rng = np.random.default_rng(31)
    lineage = np.repeat([f"L{i}" for i in range(8)], 7)
    y = rng.normal(size=8)[np.repeat(np.arange(8), 7)] + rng.normal(size=56) * 0.8
    reference = realised_share(y, lineage)
    width = reference.ci_high - reference.ci_low
    for i in (0, 13, 27, 55):
        moved = y.copy()
        moved[i] += 1.0        # one doubling on the log2 scale
        assert abs(realised_share(moved, lineage).kappa - reference.kappa) < width


# --------------------------------------------------------------------------- #
# Inputs outside the domain: a named refusal, never a quiet NaN
# --------------------------------------------------------------------------- #

_PATHOLOGICAL = st.lists(
    st.one_of(st.just(0.0), st.just(1.0), st.just(float("nan")),
              st.just(float("inf")), st.just(float("-inf")),
              st.floats(-1e300, 1e300), st.just(1e-308), st.just(0.5)),
    min_size=1, max_size=30)
_LABELS = st.lists(st.sampled_from(["A", "B", "C", "", None, "  ", "nan"]),
                   min_size=1, max_size=30)


@settings(max_examples=120, deadline=None,
          suppress_health_check=[HealthCheck.too_slow])
@given(_PATHOLOGICAL, _LABELS)
def test_a_pathological_cohort_is_refused_or_reported_as_not_estimable(values, labels):
    """Never a NaN presented as an estimate.

    Hypothesis supplies NaN, both infinities, denormals, values far outside
    the binary domain, blank and missing labels, and cohorts of one. The
    contract is that the call either raises with a sentence, or returns a
    result whose numbers are finite, or one that says it is not estimable and
    why. A NaN that travels on as if it were a number is the failure this
    looks for.
    """
    size = min(len(values), len(labels))
    y, lineage = values[:size], labels[:size]
    try:
        # An interval is asked for: with a budget of zero the package
        # documents a point-only route, and a missing interval would then be
        # the contract rather than a defect.
        result = clonal_share(y, lineage, folds=2, repeats=1, n_boot=20,
                              n_perm=9, seed=0)
    except ValueError as refusal:
        assert str(refusal), "a refusal must carry a sentence"
        return
    if not result.estimable:
        return
    for name in ("kappa", "kappa_adj", "prevalence"):
        value = getattr(result, name)
        assert value is None or np.isfinite(value), (
            f"{name} is {value!r} in a result that calls itself estimable")
    # A cohort can carry an estimate and no interval, when resampling whole
    # lineages has nothing to vary. That is allowed, and the record must then
    # say so rather than leave the missing limits unexplained.
    if not (np.isfinite(result.ci_low) and np.isfinite(result.ci_high)):
        assert result.reason, "a missing interval with no reason recorded"


@settings(max_examples=60, deadline=None,
          suppress_health_check=[HealthCheck.too_slow])
@given(_PATHOLOGICAL, _LABELS)
def test_the_evidence_is_refused_or_finite(values, labels):
    size = min(len(values), len(labels))
    try:
        result = e_process(values[:size], labels[:size], folds=2, repeats=1, seed=0)
    except ValueError as refusal:
        assert str(refusal)
        return
    if result.n_splits == 0:
        # Nothing could be scored, so there is no evidence to report. The
        # package's convention for an undefined quantity is a non-finite
        # value in memory, written as null in the record and reported as
        # unavailable; what it may not be is a number.
        assert np.isnan(result.e_value) and np.isnan(result.log_e)
        return
    assert result.e_value >= 0.0 and not np.isnan(result.e_value)
    assert not np.isnan(result.log_e)
