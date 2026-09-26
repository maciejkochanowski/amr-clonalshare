"""Convergence checks for a Gaussian grouped moment estimator."""
import numpy as np
from amr_clonalshare import censored as c


def test_variance_change_prevents_false_convergence():
    y = np.tile([-1., 1., -1., 1.], 5)
    code = np.repeat(np.arange(5), 4)
    diagnostics = {}
    c._fit_moment_iteration(y, y, code, 5, iters=2, fixed_scale=1., diagnostics=diagnostics)
    assert diagnostics['iterations'] == 2
    assert not diagnostics['converged']
    assert diagnostics['max_scaled_change'] > 1e-10


def test_result_records_numerical_convergence():
    y = np.tile([-1., 1., -1., 1.], 5)
    result = c.censored_clonal_share(y, y, np.repeat(np.arange(5), 4))
    assert isinstance(result.fit_converged, bool)
    assert result.fit_iterations > 0
    assert np.isfinite(result.fit_max_scaled_change)
