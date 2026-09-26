"""Checks of the models behind the intervals, the numerics of the MIC fit, and
the panel selection of one analysis."""
from __future__ import annotations

import numpy as np
import pytest

from amr_clonalshare import _mic_shape
from amr_clonalshare._mic_likelihood import GaussianMICLikelihood
from amr_clonalshare._mic_panels import _validated_problem, observe_panel
from amr_clonalshare.mic_null_bootstrap import _likelihood_ratio, _unrestricted_fit, calibrate_null


def _readings(kind, seed=3, groups=25):
    rng = np.random.default_rng(seed)
    sizes = rng.integers(3, 25, groups)
    lab = np.repeat(np.arange(groups), sizes)
    if kind == "gaussian":
        y = rng.normal(0, np.sqrt(.4), groups)[lab] + rng.normal(0, np.sqrt(.6), len(lab))
    else:   # a wild-type and a non-wild-type mode, membership driven by lineage
        z = (rng.normal(0, .8, groups)[lab] + rng.normal(0, .6, len(lab))) > .3
        y = -2 + 4 * z + rng.normal(0, .5, len(lab))
    edges = np.arange(-4., 5.)
    a, b = observe_panel(y, edges)
    return a, b, lab, edges


@pytest.mark.parametrize("kind,rejected", [("gaussian", False), ("bimodal", True)])
def test_shape_check_reads_the_readings_against_the_fitted_model(kind, rejected):
    a, b, lab, edges = _readings(kind)
    problem, panels, exact = _validated_problem(a, b, lab, edges)
    fit = _unrestricted_fit(problem, 24)
    test = calibrate_null(problem, panels, exact, rho=fit.rho, n_boot=49, seed=1,
                          full_fit=fit, shape_check=True)
    check = test.shape_check
    assert check["simulated_datasets"] == 49
    assert check["rejected"] is rejected
    assert check["p_value"] == pytest.approx(min(1., 2 * min(check["p_marginal"], check["p_conditional"])))


def test_shape_check_is_not_run_on_exact_readings_or_away_from_the_estimate():
    rng = np.random.default_rng(5)
    lab = np.repeat(np.arange(8), 6)
    y = rng.normal(size=48)
    problem, panels, exact = _validated_problem(y, y, lab, np.array([-1., 1.]))
    fit = _unrestricted_fit(problem, 24)
    assert _mic_shape.statistics(problem, fit, panels, exact) is None
    a, b, lab, edges = _readings("gaussian")
    problem, panels, exact = _validated_problem(a, b, lab, edges)
    test = calibrate_null(problem, panels, exact, rho=.3, n_boot=19, seed=1)
    assert test.shape_check is None


def test_refused_simulated_datasets_are_redrawn_not_counted_as_exceedances():
    # Ten readings, one above the lower of two cut points: many simulated
    # datasets put every reading in the lowest interval, which the analysis
    # refuses as a single cut point.
    lab = np.repeat(np.arange(4), [2, 2, 3, 3])
    y = np.array([-2., -1.5, -2.2, -1.8, -2.5, -1.9, -2.1, .5, -1.7, -2.3])
    a, b = observe_panel(y, np.array([0., 1.]))
    problem, panels, exact = _validated_problem(a, b, lab, np.array([0., 1.]))
    test = calibrate_null(problem, panels, exact, rho=.5, n_boot=19, seed=2)
    assert test.bootstrap_redrawn > 0
    assert len(test.bootstrap_statistics) == 19


def test_piecewise_rule_gives_the_fixed_effect_scores_analytically():
    a, b, lab, edges = _readings("gaussian", seed=7)
    level = np.where(np.arange(len(a)) % 3 == 0, "x", "y")
    p = GaussianMICLikelihood(a, b, lab, covariate=level)
    ll, grad = p.evaluate(.1, 1.1, .995, order=0, coefficients=[.3])
    h = 1e-5
    up = p.evaluate(.1, 1.1, .995, order=0, coefficients=[.3 + h], beta_gradient=False)[0]
    down = p.evaluate(.1, 1.1, .995, order=0, coefficients=[.3 - h], beta_gradient=False)[0]
    assert grad[3] == pytest.approx((up - down) / (2 * h), rel=1e-5, abs=1e-5)


def test_a_null_maximum_above_the_unrestricted_one_gives_a_zero_statistic():
    a, b, lab, edges = _readings("gaussian", seed=9)
    p = GaussianMICLikelihood(a, b, lab)
    fit = _unrestricted_fit(p, 24)
    _, _, statistic = _likelihood_ratio(p, fit.rho, 24, fit)
    assert statistic == pytest.approx(0., abs=1e-6)


def test_population_mixing_check_separates_a_clone_from_a_gaussian_law():
    from amr_clonalshare.population_probit import population_probit_icc
    rng = np.random.default_rng(11)
    sizes = np.full(50, 10)
    gauss = rng.binomial(sizes, 1 / (1 + np.exp(-(rng.normal(0, 1.2, 50) - 1.))))
    clone = np.where(rng.random(50) < .2, rng.binomial(sizes, .95), rng.binomial(sizes, .03))
    assert population_probit_icc(clone, sizes).mixing_check["rejected"] is True
    kept = population_probit_icc(gauss, sizes).mixing_check
    assert kept["rejected"] is False and kept["simulated_datasets"] >= 90


def test_single_analysis_selects_agents_by_the_permutation_p_values():
    from amr_clonalshare.core import _lineage_selection
    shares = {name: {"p_value": p, "p_floor": .001} for name, p in
              (("a", .001), ("b", .001), ("c", .001), ("d", .4), ("e", .9))}
    sel = _lineage_selection(shares, .05)
    assert sel["method"] == "benjamini_yekutieli"
    assert sel["rejected_features"] == ["a", "b", "c"]
    assert _lineage_selection({"x": {"p_value": float("nan")}}, .05)["n_tested"] == 0
