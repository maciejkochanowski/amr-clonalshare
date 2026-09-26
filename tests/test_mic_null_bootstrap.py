"""The null model, not the unconstrained estimate, generates calibration data."""
import importlib
import numpy as np
import pytest


def api():
    spec = importlib.util.find_spec('amr_clonalshare.mic_null_bootstrap')
    assert spec is not None, 'Null-wise bootstrap inference is not implemented'
    return importlib.import_module('amr_clonalshare.mic_null_bootstrap')


def data():
    rng = np.random.default_rng(612)
    labels = np.repeat(np.arange(8), [3, 4, 8, 9, 4, 3, 8, 9])
    y = rng.normal(size=8)[labels] + rng.normal(size=len(labels))
    return y, labels


def test_calibration_uses_tested_rho_and_retains_every_draw():
    y, labels = data()
    r = api().null_bootstrap_test(y, y, labels, rho=.2,
          panel_edges=[-1., 0., 1.], n_boot=19, seed=122)
    assert r.reportable
    assert r.null_fit.rho == .2
    assert len(r.bootstrap_statistics) == 19
    expected = (1 + np.sum(np.array(r.bootstrap_statistics) >= r.statistic))/20
    assert r.pvalue == expected
    assert r.accepted == (r.pvalue > .05)


def test_exact_null_bootstrap_is_location_scale_equivariant():
    y, labels = data()
    first = api().null_bootstrap_test(y, y, labels, rho=.4,
        panel_edges=[-1., 0., 1.], n_boot=19, seed=77)
    z = 3*y + 7
    second = api().null_bootstrap_test(z, z, labels, rho=.4,
        panel_edges=[4., 7., 10.], n_boot=19, seed=77)
    assert first.pvalue == second.pvalue
    np.testing.assert_allclose(first.bootstrap_statistics,
                               second.bootstrap_statistics, atol=2e-6)
    assert first.statistic == pytest.approx(second.statistic, abs=2e-6)


@pytest.mark.parametrize('rho', [-.1, 1., np.nan, True])
def test_invalid_null_is_rejected(rho):
    y, labels = data()
    with pytest.raises(ValueError):
        api().null_bootstrap_test(y, y, labels, rho=rho,
            panel_edges=[-1., 0., 1.], n_boot=19, seed=2)


def test_test_inversion_contains_the_likelihood_maximum():
    y, labels = data()
    result = api().null_bootstrap_interval(y, y, labels,
        panel_edges=[-1., 0., 1.], n_boot=19, seed=441, tolerance=.03)
    assert result.reportable
    assert 0 <= result.low <= result.fit.rho <= result.high <= 1
    assert len(result.tests) > 2
    assert all(test.null_fit.rho == test.rho for test in result.tests)
    assert result.method == 'nullwise_parametric_bootstrap_LR_inversion'


def test_unresolved_observed_fit_has_no_executed_bootstrap(monkeypatch):
    from dataclasses import replace
    from amr_clonalshare._mic_likelihood import GaussianMICLikelihood
    y, labels = data()
    problem = GaussianMICLikelihood(y, y, labels)
    fit = GaussianMICLikelihood.fit
    monkeypatch.setattr(GaussianMICLikelihood, 'fit', lambda self, **kw: replace(
        fit(self, **kw), converged=False, message='deliberate test'))
    result = api().calibrate_null(problem, np.array([-1., 0., 1.]),
        np.ones(len(y), dtype=bool), rho=.3, n_boot=19, seed=5)
    assert not result.reportable
    assert result.bootstrap_attempted == 0
    assert result.pvalue == 1.
    assert result.accepted


def test_parallel_bootstrap_fits_give_identical_results():
    y, labels = data()
    from amr_clonalshare._mic_panels import observe_panel
    edges = [-3., -2., -1., 0., 1., 2., 3.]
    a, b = observe_panel(y, edges)
    kwargs = dict(panel_edges=edges, n_boot=19, seed=31)
    serial = api().null_bootstrap_test(a, b, labels, rho=.3, **kwargs)
    parallel = api().null_bootstrap_test(a, b, labels, rho=.3, workers=2, **kwargs)
    assert serial.bootstrap_statistics == parallel.bootstrap_statistics
    assert serial.pvalue == parallel.pvalue


@pytest.mark.parametrize('workers', [-1, 1.5, True])
def test_invalid_worker_count_is_rejected(workers):
    y, labels = data()
    with pytest.raises(ValueError):
        api().null_bootstrap_test(y, y, labels, rho=.2,
            panel_edges=[-1., 0., 1.], n_boot=19, seed=2, workers=workers)


def test_zero_workers_means_every_available_cpu():
    from amr_clonalshare.mic_null_bootstrap import _workers
    assert _workers(0) >= 1
    assert _workers(3) == 3
    with pytest.raises(ValueError):
        _workers(-1)


def test_records_read_on_different_panels_are_simulated_on_their_own_panels():
    y, labels = data()
    from amr_clonalshare._mic_panels import observe_panels, _panels
    # one laboratory read a single cut point, a binary reading of its records
    panel = np.where(np.arange(y.size) % 2 == 0, 'narrow', 'wide')
    edges = {'narrow': [0.], 'wide': [-3., -2., -1., 0., 1., 2., 3.]}
    a, b = observe_panels(y, _panels(edges, panel, y.size))
    result = api().null_bootstrap_test(a, b, labels, rho=.3, panel_edges=edges,
                                       panel=panel, n_boot=19, seed=31)
    assert result.reportable and result.bootstrap_attempted == 19
    with pytest.raises(ValueError):
        api().null_bootstrap_test(a, b, labels, rho=.3, panel_edges=edges['wide'],
                                  n_boot=19, seed=31)
    with pytest.raises(ValueError):
        api().null_bootstrap_test(a, b, labels, rho=.3, panel_edges={'narrow': edges['narrow']},
                                  panel=panel, n_boot=19, seed=31)
    with pytest.raises(ValueError):
        api().null_bootstrap_test(a, b, labels, rho=.3, panel_edges={'narrow': [0.], 'wide': [1.]},
                                  panel=panel, n_boot=19, seed=31)


def test_covariate_is_refitted_under_the_null_and_in_every_simulated_dataset():
    rng=np.random.default_rng(23)
    labels=np.repeat(np.arange(20),8); n=len(labels)
    laboratory=np.where(rng.random(n)<.85,labels%2,1-labels%2)
    y=rng.normal(0,np.sqrt(.3),20)[labels]+rng.normal(0,np.sqrt(.7),n)+1.*laboratory
    from amr_clonalshare._mic_panels import observe_panel
    edges=[-3.,-2.,-1.,0.,1.,2.,3.]
    a,b=observe_panel(y,edges)
    result=api().null_bootstrap_test(a,b,labels,rho=.3,panel_edges=edges,n_boot=19,seed=5,covariate=laboratory)
    assert result.reportable and result.bootstrap_attempted == 19
    assert len(result.unrestricted_fit.coefficients) == 1 and len(result.null_fit.coefficients) == 1
    assert result.null_fit.rho == .3
    plain=api().null_bootstrap_test(a,b,labels,rho=.3,panel_edges=edges,n_boot=19,seed=5)
    assert plain.unrestricted_fit.coefficients == ()
    with pytest.raises(ValueError):
        api().null_bootstrap_test(a,b,labels,rho=.3,panel_edges=edges,n_boot=19,seed=5,covariate=laboratory[:-1])


@pytest.mark.parametrize('edges', [[-1., 0., 1.], [-3., -2., -1., 0., 1., 2., 3.]])
@pytest.mark.parametrize('rho', [0., .3, .9])
def test_bootstrap_computes_the_observed_statistic_on_identical_data(edges, rho):
    # T* must be the same function of the data as T; otherwise the bootstrap
    # calibrates a different statistic from the one it is compared with.
    from amr_clonalshare._mic_panels import observe_panel
    y, labels = data()
    a, b = observe_panel(y, edges)
    observed = api().null_bootstrap_test(a, b, labels, rho=rho,
        panel_edges=edges, n_boot=19, seed=3).statistic
    replicate = api()._bootstrap_statistic((a, b, labels, rho, 24, None))
    assert replicate == observed


def test_failed_simulations_retain_the_tested_value(monkeypatch):
    # A simulated dataset without a statistic counts as an exceedance: with
    # every simulation failed the value cannot be rejected, and the observed
    # fit is still a result, so the inversion goes on rather than stopping.
    y, labels = data()
    module = api()
    monkeypatch.setattr(module, '_bootstrap_statistic', lambda job: float('inf'))
    r = module.null_bootstrap_test(y, y, labels, rho=.2, panel_edges=[-1., 0., 1.], n_boot=19, seed=4)
    assert r.reportable and r.accepted
    assert r.pvalue == 1. and r.critical == float('inf')
    assert r.bootstrap_failed == 19
    assert r.status.endswith('failed_simulations_retain_value')


def test_a_simulated_level_censored_on_one_side_still_gives_a_statistic():
    # The observed data must identify every covariate level; a simulated
    # dataset in which one level fell wholly beyond the panel is scored by
    # the same box-bounded likelihood ratio instead of being lost.
    from amr_clonalshare._mic_likelihood import GaussianMICLikelihood
    y, labels = data()
    edges = [-1., 0., 1.]
    from amr_clonalshare._mic_panels import observe_panel
    a, b = observe_panel(y, edges)
    level = np.where(labels == 7, 'high', 'rest')
    a[labels == 7], b[labels == 7] = 1., np.inf
    with pytest.raises(ValueError, match='censored on one side'):
        GaussianMICLikelihood(a, b, labels, covariate=level)
    statistic = api()._bootstrap_statistic((a, b, labels, .3, 24, [level]))
    assert np.isfinite(statistic)


def test_every_candidate_of_a_box_fit_is_fitted_over_the_box(monkeypatch):
    # The rho = 0 candidate of an unrestricted fit is a fit of its own; if it
    # were not told that maxima on the box count, a dataset whose best fit
    # is that candidate with a coefficient on the box came back unresolved.
    from amr_clonalshare._mic_likelihood import GaussianMICLikelihood
    y, labels = data()
    problem = GaussianMICLikelihood(y, y + 1., labels)
    calls = []
    fit = GaussianMICLikelihood.fit
    def spy(self, **kwargs):
        calls.append(kwargs)
        return fit(self, **kwargs)
    monkeypatch.setattr(GaussianMICLikelihood, 'fit', spy)
    problem.fit(box_maximum=True)
    zero = [k for k in calls if k.get('rho') == 0.]
    assert zero and all(k.get('box_maximum') for k in zero)


def test_an_interval_names_values_kept_only_by_failed_simulations(monkeypatch):
    y, labels = data()
    module = api()
    real = module._bootstrap_statistic
    def fail_high(job):
        return float('inf') if job[3] > .95 else real(job)
    monkeypatch.setattr(module, '_bootstrap_statistic', fail_high)
    r = module.null_bootstrap_interval(y, y, labels, panel_edges=[-1., 0., 1.], n_boot=19, seed=441,
                                       tolerance=.03)
    assert r.reportable and r.high == 1.
    assert 'retained_by_failed_simulations: 0.97, 0.995' in r.status
