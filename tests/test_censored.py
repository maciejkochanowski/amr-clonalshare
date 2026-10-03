"""MIC readings as intervals on the log2 scale, and the panel they were read on."""
import numpy as np
import pytest

from amr_clonalshare.censored import (intervals_by_panel, intervals_from_binary,
                                      intervals_from_mic, panel_geometry)


def test_binary_constructor_brackets_the_cut():
    lo, hi = intervals_from_binary(np.array([0., 1., np.nan]), cutoff_log2=1.)
    assert not np.isfinite(lo[0]) and hi[0] == 1.
    assert lo[1] == 1. and not np.isfinite(hi[1])
    assert np.isnan(lo[2]) and np.isnan(hi[2])
    with pytest.raises(ValueError):
        intervals_from_binary(np.array([0., 2.]))
    with pytest.raises(ValueError):
        intervals_from_binary(np.array([0., 1.]), cutoff_log2=np.inf)


def test_panel_geometry_counts_the_readings_on_the_end_wells():
    values = np.array([128.0] * 40 + [1.0, 2.0, 4.0, 8.0, 16.0])
    g = panel_geometry(values)
    assert g.share_on_highest == pytest.approx(40 / 45)
    assert g.share_on_lowest == pytest.approx(1 / 45)
    assert panel_geometry([]).n_wells == 0


def test_intervals_from_mic_bracket_the_reading():
    values = np.array([1.0, 2.0, 4.0, 8.0, 2.0, 4.0])
    lo, hi = intervals_from_mic(values, treat_end_wells_as_censored=False)
    inner = values > values.min()
    assert np.all(lo[inner] < np.log2(values[inner]))
    assert np.all(hi[inner] == np.log2(values[inner]))


def test_end_wells_are_censored_by_default():
    values = np.array([1.0, 1.0, 2.0, 4.0, 8.0, 8.0])
    lo, hi = intervals_from_mic(values)
    assert np.all(~np.isfinite(lo[values == 1.0]))
    assert np.all(~np.isfinite(hi[values == 8.0]))


def test_recorded_operator_wins_over_the_end_well_heuristic():
    values = np.array([2.0, 4.0, 8.0, 4.0])
    ops = np.array(["", "", "", ">"], dtype=object)
    lo, hi = intervals_from_mic(values, operators=ops,
                                treat_end_wells_as_censored=True)
    assert lo[3] == pytest.approx(2.0)
    assert not np.isfinite(hi[3])


@pytest.mark.parametrize('sign', ['>', '>='])
def test_blank_sign_on_the_highest_well_is_exact_beside_a_censoring_sign(sign):
    """An export that writes the sign only on censored readings: a blank 32
    beside a >32 is the reading of the 32 well, not a right-censored one,
    and the readings of one panel then do not overlap."""
    from amr_clonalshare.latent_order import Readings
    values = np.array([8.0, 16.0, 32.0, 32.0, 32.0, 32.0])
    lo, hi = intervals_from_mic(values, operators=['', '', '', '', '', sign])
    assert lo[2] == pytest.approx(4.0) and hi[2] == pytest.approx(5.0)
    assert lo[3] == pytest.approx(4.0) and hi[3] == pytest.approx(5.0)
    assert lo[5] == pytest.approx(5.0 if sign == '>' else 4.0)
    assert not np.isfinite(hi[5])
    if sign == '>':
        readings = Readings(lo, hi, np.array([0, 0, 1, 1, 2, 2]), np.zeros(6, dtype=object))
        assert readings.sharp


def test_blank_signs_alone_keep_the_end_well_heuristic():
    """Without a censoring sign on the highest well a blank reading there is
    read as the heuristic reads it, whatever other signs the call carries."""
    values = np.array([8.0, 16.0, 32.0, 32.0, 2.0])
    lo, hi = intervals_from_mic(values, operators=['', '', '', '', '<='])
    assert np.all(lo[values == 32.0] == pytest.approx(4.0))
    assert np.all(~np.isfinite(hi[values == 32.0]))
    lo, hi = intervals_from_mic(values[:4], operators=['', '', '', ''])
    assert np.all(~np.isfinite(hi[2:]))


def test_blank_sign_on_the_lowest_well_stays_censored_beside_a_less_sign():
    """No growth at the lowest well is what <= records, so a blank reading
    there is the same reading and the two do not overlap."""
    from amr_clonalshare.latent_order import Readings
    values = np.array([0.5, 0.5, 1.0, 2.0, 0.5, 1.0])
    lo, hi = intervals_from_mic(values, operators=['<=', '', '', '', '', ''])
    assert not np.isfinite(lo[0]) and not np.isfinite(lo[1])
    assert hi[0] == pytest.approx(-1.0) and hi[1] == pytest.approx(-1.0)
    readings = Readings(lo, hi, np.array([0, 0, 1, 1, 2, 2]), np.zeros(6, dtype=object))
    assert readings.sharp


@pytest.mark.parametrize('sign', ['>=', '\u2265'])
def test_greater_or_equal_includes_the_recorded_well(sign):
    # >=4 on a doubling panel: the MIC is 4 or more, so above the 2 well.
    values = np.array([1.0, 2.0, 4.0, 8.0, 4.0])
    lo, hi = intervals_from_mic(values, operators=['', '', '', '', sign])
    assert lo[4] == pytest.approx(1.0)
    assert not np.isfinite(hi[4])


@pytest.mark.parametrize('sign', ['<', '<=', '\u2264'])
def test_less_than_reads_at_or_below_the_recorded_well(sign):
    values = np.array([1.0, 2.0, 4.0, 1.0])
    lo, hi = intervals_from_mic(values, operators=['', '', '', sign])
    assert not np.isfinite(lo[3])
    assert hi[3] == pytest.approx(0.0)


def test_panel_geometry_reports_a_missing_dilution():
    """A dilution no isolate landed on was still tested: the lattice holds
    five wells, four of them recorded, and the recorded ratio of 4 says where
    the unoccupied one sits."""
    values = np.array([1.0, 2.0, 8.0, 16.0])
    g = panel_geometry(values)
    assert g.doubling and g.lattice == "doubling"
    assert g.n_wells == 5 and g.n_wells_recorded == 4
    assert 4.0 in g.lattice_ratios
    lo, hi = intervals_from_mic(values, treat_end_wells_as_censored=False)
    assert lo[2] == pytest.approx(2.0) and hi[2] == pytest.approx(3.0)


def test_rounding_variants_of_one_dilution_share_a_well():
    """0.06 and 0.064 are one well recorded two ways, so a reading of 0.064
    is bounded below by the previous doubling and not by 0.06."""
    values = np.array([0.06, 0.064, 0.125, 0.25, 0.5])
    g = panel_geometry(values)
    assert g.n_wells == 4 and g.n_wells_recorded == 5
    lo, hi = intervals_from_mic(values, treat_end_wells_as_censored=False)
    assert hi[0] == hi[1] == pytest.approx(-4.0)
    assert lo[2] == pytest.approx(-4.0)


def test_recorded_wells_place_the_readings():
    values = np.array([1.0, 4.0, 4.0, 16.0])
    lo, hi = intervals_from_mic(values, wells=[0.5, 1, 2, 4, 8, 16, 32],
                                treat_end_wells_as_censored=False)
    assert lo[0] == pytest.approx(-1.0) and hi[0] == pytest.approx(0.0)
    assert lo[1] == pytest.approx(1.0) and hi[1] == pytest.approx(2.0)
    assert panel_geometry(values, wells=[0.5, 1, 2, 4, 8, 16, 32]).n_wells == 7


def test_readings_off_any_doubling_lattice_are_their_own_wells():
    values = np.array([1.0, 1.5, 2.0, 3.0])
    g = panel_geometry(values)
    assert g.lattice == "recorded" and not g.doubling and g.n_wells == 4


def test_each_laboratory_keeps_its_own_end_wells():
    # Laboratory B tested one dilution lower, so B's 2.0 is an interior well
    # while A's 2.0 is A's lowest well and left-censored.
    values = np.array([2.0, 4.0, 8.0, 1.0, 2.0, 8.0])
    panel = np.array(["A", "A", "A", "B", "B", "B"], dtype=object)
    lo, hi = intervals_by_panel(values, panel)
    pooled_lo, _ = intervals_from_mic(values)
    assert not np.isfinite(lo[0]) and np.isfinite(pooled_lo[0])
    assert lo[4] == pytest.approx(0.0) and hi[4] == pytest.approx(1.0)
    with pytest.raises(ValueError):
        intervals_by_panel(values, panel[:-1])
