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


def test_unresolved_observed_fit_has_no_executed_bootstrap():
    from dataclasses import replace
    from amr_clonalshare._mic_likelihood import GaussianMICLikelihood
    y, labels = data()
    problem = GaussianMICLikelihood(y, y, labels)
    failed = replace(problem.fit(), converged=False, message='deliberate test')
    result = api().calibrate_null(problem, np.array([-1., 0., 1.]),
        np.ones(len(y), dtype=bool), rho=failed.rho, n_boot=19, seed=5,
        full_fit=failed)
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
