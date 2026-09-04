"""Interval-censored clonal share: reductions, identification, calibration.

The claim this module has to earn is that a dichotomised call, a recorded
dilution and a censored reading are one likelihood at three interval widths.
Two tests hold that claim to account. The first fixes the estimator against a
closed form on exact data; the second fixes the coarsest rung of the ladder
against the constructor that produces it. The rest are the properties a
reviewer is entitled to check without reading the source.
"""
import numpy as np
import pytest

from amr_clonalshare.censored import (CENSORED_GROUP_LIMIT,
                                         SINGLE_CUT_PREVALENCE,
                                         censored_clonal_share,
                                         intervals_from_binary,
                                         intervals_from_mic, panel_geometry,
                                         profile_interval,
                                         scale_is_identified,
                                         sensitivity_endpoints)

TAU2 = 1.0
SIG2 = 1.0
TRUE_SHARE = TAU2 / (TAU2 + SIG2)


def _latent(rng, sizes):
    """A one-way random-effects cohort on the latent log2 scale."""
    code = np.repeat(np.arange(len(sizes)), sizes)
    mu = rng.normal(0.0, np.sqrt(TAU2), len(sizes))
    z = mu[code] + rng.normal(0.0, np.sqrt(SIG2), code.size)
    lineage = np.array(["L%d" % g for g in code], dtype=object)
    return z, code, lineage


def _closed_form_icc(z, code):
    """One-way variance-component ratio, textbook form. The estimator under
    test must reproduce this on a balanced exact design."""
    G = int(code.max()) + 1
    cnt = np.bincount(code, minlength=G).astype(float)
    group_mean = np.bincount(code, weights=z, minlength=G) / np.maximum(cnt, 1)
    grand = float(z.mean())
    n, gp = z.size, int((cnt > 0).sum())
    msw = float(((z - group_mean[code]) ** 2).sum() / (n - gp))
    msb = float((cnt * (group_mean - grand) ** 2).sum() / (gp - 1))
    n0 = float((n - (cnt ** 2).sum() / n) / (gp - 1))
    tau2 = max(0.0, (msb - msw) / n0)
    return tau2 / (tau2 + msw)


# --------------------------------------------------------------------------
# Reduction 1: zero-width intervals are exact readings
# --------------------------------------------------------------------------
def test_zero_width_intervals_reproduce_the_closed_form():
    """With every interval collapsed to a point the estimator must agree with
    the classical one-way variance-component ratio. A balanced design is used
    because that is where the two coincide by construction; on an unbalanced
    one the restricted-likelihood estimator is deliberately different, and
    better, which the next test measures."""
    rng = np.random.default_rng(11)
    z, code, lineage = _latent(rng, [25] * 30)
    got = censored_clonal_share(z, z, lineage, n_boot=0).kappa
    assert got == pytest.approx(_closed_form_icc(z, code), abs=1e-3)


def test_unbalanced_design_beats_the_moment_estimator():
    """The departure from the closed form on an unbalanced design is not a
    defect. Over repeated cohorts with group sizes from 2 to 90 the estimator
    must be closer to the truth than the moment estimator it departs from."""
    sizes = [2, 3, 4, 5, 8, 60, 12, 90, 6, 30, 7, 45, 9, 33, 11, 20, 15, 4, 55, 18]
    rng = np.random.default_rng(11)
    ours, moment = [], []
    for _ in range(40):
        z, code, lineage = _latent(rng, sizes)
        ours.append(censored_clonal_share(z, z, lineage, n_boot=0).kappa)
        moment.append(_closed_form_icc(z, code))
    ours, moment = np.asarray(ours), np.asarray(moment)
    assert abs(ours.mean() - TRUE_SHARE) < abs(moment.mean() - TRUE_SHARE)
    assert ours.std(ddof=1) < moment.std(ddof=1)


# --------------------------------------------------------------------------
# Reduction 2: a single cut point is the coarsest rung of the same ladder
# --------------------------------------------------------------------------
def test_dichotomising_a_mic_gives_the_binary_constructor_intervals():
    """Coarsening a recorded MIC at a cut-off and taking the dichotomised call
    to the binary constructor must produce the identical intervals, and
    therefore the identical share. This is the join in the ladder: if the two
    entry points disagreed, the three input modes would not be one estimator."""
    rng = np.random.default_rng(3)
    z, _code, lineage = _latent(rng, [20] * 25)
    cut = 1.0
    from_call = intervals_from_binary((z > cut).astype(float), cutoff_log2=cut)
    by_hand = (np.where(z > cut, cut, -np.inf), np.where(z > cut, np.inf, cut))
    assert np.array_equal(from_call[0], by_hand[0])
    assert np.array_equal(from_call[1], by_hand[1])
    a = censored_clonal_share(*from_call, lineage, n_boot=0).kappa
    b = censored_clonal_share(*by_hand, lineage, n_boot=0).kappa
    assert a == pytest.approx(b, abs=1e-12)


def test_the_ladder_recovers_one_latent_share_from_three_readings():
    """Point, interval and single-cut readings of the same latent values must
    land on the same number. The bias is allowed to grow as the reading
    coarsens; it is not allowed to change the answer."""
    wells = np.arange(-4.0, 6.0)
    rng = np.random.default_rng(11)
    point, interval, binary = [], [], []
    for _ in range(12):
        z, _code, lineage = _latent(rng, [25] * 30)
        point.append(censored_clonal_share(z, z, lineage, n_boot=0).kappa)
        idx = np.searchsorted(wells, z, side="left")
        lo = np.where(idx > 0, wells[np.clip(idx - 1, 0, wells.size - 1)], -np.inf)
        hi = np.where(idx < wells.size, wells[np.clip(idx, 0, wells.size - 1)],
                      np.inf)
        interval.append(censored_clonal_share(lo, hi, lineage, n_boot=0).kappa)
        blo, bhi = intervals_from_binary((z > 1.0).astype(float), cutoff_log2=1.0)
        binary.append(censored_clonal_share(blo, bhi, lineage, n_boot=0).kappa)
    for name, got, tol in (("point", point, 0.03), ("interval", interval, 0.04),
                           ("binary", binary, 0.06)):
        assert abs(float(np.mean(got)) - TRUE_SHARE) < tol, name


# --------------------------------------------------------------------------
# Properties
# --------------------------------------------------------------------------
def test_share_is_near_zero_when_lineage_carries_nothing():
    rng = np.random.default_rng(5)
    code = np.repeat(np.arange(30), 25)
    z = rng.normal(0.0, 1.0, code.size)
    lineage = np.array(["L%d" % g for g in code], dtype=object)
    assert censored_clonal_share(z, z, lineage, n_boot=0).kappa < 0.06


def test_scale_is_not_identified_from_a_single_cut_point():
    y = np.array([0, 1, 1, 0, 1, 0, 0, 1], dtype=float)
    lo, hi = intervals_from_binary(y)
    assert not scale_is_identified(lo, hi)
    lo2, hi2 = intervals_from_mic(np.array([1.0, 2, 4, 8, 1, 2, 4, 8]))
    assert scale_is_identified(lo2, hi2)


def test_binary_result_says_the_scale_is_a_convention():
    rng = np.random.default_rng(7)
    z, _code, lineage = _latent(rng, [20] * 20)
    lo, hi = intervals_from_binary((z > 0.0).astype(float))
    res = censored_clonal_share(lo, hi, lineage, n_boot=0)
    assert "not identified" in res.reason
    assert res.sigma_within == pytest.approx(1.0)


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


def test_panel_geometry_refuses_the_point_mode_on_a_piled_end_well():
    values = np.array([128.0] * 40 + [1.0, 2.0, 4.0, 8.0, 16.0])
    g = panel_geometry(values)
    assert "point" not in g.admissible_modes
    assert "interval" in g.admissible_modes
    assert g.share_on_highest > 0.5


def test_panel_geometry_reports_a_missing_dilution():
    values = np.array([1.0, 2.0, 8.0, 16.0])
    g = panel_geometry(values)
    assert not g.doubling
    assert 4.0 in g.lattice_ratios


def test_lineage_length_mismatch_is_refused():
    lo, hi = intervals_from_binary(np.zeros(10))
    with pytest.raises(ValueError):
        censored_clonal_share(lo, hi, ["A"] * 9, n_boot=0)


def test_missing_readings_are_dropped_with_their_lineage():
    rng = np.random.default_rng(9)
    z, _code, lineage = _latent(rng, [20] * 20)
    lo, hi = z.copy(), z.copy()
    lo[:7] = np.nan
    hi[:7] = np.nan
    res = censored_clonal_share(lo, hi, lineage, n_boot=0)
    assert res.n == z.size - 7


def test_sensitivity_reports_both_readings_of_the_end_wells():
    rng = np.random.default_rng(13)
    z, _code, lineage = _latent(rng, [20] * 20)
    mic = 2.0 ** np.clip(np.round(z), -3, 3)
    out = sensitivity_endpoints(mic, lineage, n_boot=40, seed=1)
    assert set(out) >= {"end_wells_censored", "end_wells_exact",
                        "kappa_bracket_width", "interval_width",
                        "assumption_dominates"}
    assert out["kappa_bracket_width"] >= 0.0


def test_too_few_lineages_is_refused_rather_than_estimated():
    rng = np.random.default_rng(2)
    z = rng.normal(size=40)
    res = censored_clonal_share(z, z, ["A"] * 40, n_boot=0)
    assert not res.estimable
    assert np.isnan(res.kappa)


def test_the_reported_interval_brackets_the_point_estimate():
    rng = np.random.default_rng(17)
    z, _code, lineage = _latent(rng, [20] * 25)
    res = censored_clonal_share(z, z, lineage, n_boot=200, seed=4)
    assert res.ci_low <= res.kappa <= res.ci_high
    assert 0.0 <= res.ci_low <= res.ci_high <= 1.0
    assert res.boot_low <= res.kappa <= res.boot_high


def test_the_reported_interval_covers_a_known_truth():
    """The interval the package reports is the variance-ratio one, chosen
    because the cluster bootstrap was measured and found to cover a median
    0.61 of the time. Twenty-five cohorts is a coarse check, not the
    calibration, so the bar is set where a nominal 0.95 procedure fails it
    about once in a thousand runs."""
    rng = np.random.default_rng(31)
    hits = 0
    for _ in range(25):
        z, _code, lineage = _latent(rng, [25] * 30)
        res = censored_clonal_share(z, z, lineage, n_boot=0, profile=False)
        hits += int(res.ci_low <= TRUE_SHARE <= res.ci_high)
    assert hits >= 20


def test_the_interval_widens_as_the_reading_coarsens():
    """More lineages, same isolates: an exact reading pins the share more
    tightly than a single cut point on the same cohort."""
    rng = np.random.default_rng(41)
    z, _code, lineage = _latent(rng, [25] * 30)
    exact = censored_clonal_share(z, z, lineage, n_boot=0, profile=False)
    lo, hi = intervals_from_binary((z > 0.0).astype(float))
    call = censored_clonal_share(lo, hi, lineage, n_boot=0, profile=False)
    assert (exact.ci_high - exact.ci_low) < (call.ci_high - call.ci_low)


def test_a_cut_in_a_tail_is_refused_and_the_dilution_is_not():
    """A single cut point with almost every isolate on one side leaves no
    contrast to divide. The recorded concentration on the same cohort does,
    and that difference is the practical argument for reading the panel."""
    rng = np.random.default_rng(51)
    z, _code, lineage = _latent(rng, [25] * 30)
    cut = float(np.quantile(z, 0.97))
    lo, hi = intervals_from_binary((z > cut).astype(float), cutoff_log2=cut)
    call = censored_clonal_share(lo, hi, lineage, n_boot=0, profile=False)
    assert not call.estimable
    assert "calibrated window" in call.reason
    wells = np.arange(-4.0, 6.0)
    idx = np.searchsorted(wells, z, side="left")
    dlo = np.where(idx > 0, wells[np.clip(idx - 1, 0, wells.size - 1)], -np.inf)
    dhi = np.where(idx < wells.size, wells[np.clip(idx, 0, wells.size - 1)],
                   np.inf)
    dilution = censored_clonal_share(dlo, dhi, lineage, n_boot=0, profile=False)
    assert dilution.estimable


def test_a_cut_inside_the_window_is_not_refused_for_prevalence():
    rng = np.random.default_rng(53)
    z, _code, lineage = _latent(rng, [25] * 30)
    lo, hi = intervals_from_binary((z > 0.0).astype(float))
    res = censored_clonal_share(lo, hi, lineage, n_boot=0, profile=False)
    assert SINGLE_CUT_PREVALENCE[0] < 0.5 < SINGLE_CUT_PREVALENCE[1]
    assert res.estimable
    assert "calibrated window" not in res.reason


# --------------------------------------------------------------------------
# The estimability gate reads the isolate share, not the lineage count
# --------------------------------------------------------------------------
def _left_censored(rng, cut, sizes=(25,) * 30):
    z, code, lineage = _latent(rng, list(sizes))
    below = z <= cut
    return np.where(below, -np.inf, cut), np.where(below, cut, np.inf), lineage


def test_a_few_wholly_censored_lineages_do_not_refuse_the_agent():
    """On a susceptible agent some lineages sit entirely below the lowest
    tested dilution. Refusing the agent for that would refuse exactly the
    result a veterinary laboratory wants, and the calibration says the
    shrinkage absorbs it: the bias stays under 0.005 while a third of the
    isolates are in such lineages."""
    rng = np.random.default_rng(101)
    lo, hi, lineage = _left_censored(rng, cut=1.0)
    res = censored_clonal_share(lo, hi, lineage, n_boot=0)
    assert res.n_groups_fully_censored > 0
    assert res.share_in_censored_groups < CENSORED_GROUP_LIMIT
    assert res.estimable
    assert abs(res.kappa - TRUE_SHARE) < 0.12


def test_a_cohort_mostly_beyond_the_panel_is_refused():
    rng = np.random.default_rng(101)
    lo, hi, lineage = _left_censored(rng, cut=3.0)
    res = censored_clonal_share(lo, hi, lineage, n_boot=0)
    assert res.share_in_censored_groups > CENSORED_GROUP_LIMIT
    assert not res.estimable
    assert "not estimable" in res.reason


def test_the_gate_names_the_quantity_it_read():
    rng = np.random.default_rng(101)
    lo, hi, lineage = _left_censored(rng, cut=3.0)
    reason = censored_clonal_share(lo, hi, lineage, n_boot=0).reason
    assert "%" in reason and "wholly beyond the panel" in reason


def test_the_profile_width_answers_to_the_level_it_was_asked_for():
    """``profile_interval`` accepts ``alpha`` and must use it. The cut was a
    hard-coded 0.95 chi-squared quantile, so the width did not move with the
    argument. The default still reproduces that quantile bit for bit, so no
    recorded number changes; what changes is that a level other than 0.05 now
    reaches the endpoints.

    The design can fail rather than merely pass: on an inestimable cohort every
    arm returns NaN and the comparison would be vacuous, so the peak is
    required to be finite before the widths are compared. The grid is refined
    because the endpoints are grid points and a coarse grid can tie."""
    rng = np.random.default_rng(7)
    z, _, lineage = _latent(rng, [8] * 30)
    width = {}
    for alpha in (0.5, 0.05, 0.01):
        lower, upper, peak = profile_interval(z, z, lineage, alpha=alpha,
                                              grid=201)
        assert np.isfinite(peak)
        width[alpha] = upper - lower
    assert width[0.5] < width[0.05] < width[0.01]
