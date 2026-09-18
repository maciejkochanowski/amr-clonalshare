"""Null-wise calibration of a fully profiled Gaussian MIC likelihood.

Every tested rho has its own nuisance fit and bootstrap distribution. The
panel and lineage sizes are retained. Nuisance parameters are still plugged
in: this is not a finite-sample exact procedure for censored observations.
Failed bootstrap fits count as infinite statistics, not deleted replicates.
"""
from __future__ import annotations
from concurrent.futures import ProcessPoolExecutor
from contextlib import nullcontext
from dataclasses import asdict, dataclass
from numbers import Integral
import os
import numpy as np
from ._mic_likelihood import GaussianMICLikelihood, MICFit
from ._mic_panels import _as_panels, _record_covariate, _record_offset, _validated_problem, observe_panels
from .mic_inference import _alpha


@dataclass(frozen=True)
class NullBootstrapTest:
    rho: float
    pvalue: float
    statistic: float
    critical: float
    alpha: float
    accepted: bool
    unrestricted_fit: MICFit
    null_fit: MICFit
    bootstrap_statistics: tuple[float, ...]
    bootstrap_failed: int
    reportable: bool
    status: str
    bootstrap_attempted: int = 0

    def as_dict(self):
        return asdict(self)


def _bootstrap_statistic(job):
    a, b, code, rho, start, order, covariate = job
    try:
        sample = GaussianMICLikelihood(a, b, code, covariate=covariate)
        box = dict(order=order, box_maximum=True)
        unrestricted = sample.fit(start=start, multistart=False, **box)
        restricted = sample.fit(rho=rho, start=unrestricted, multistart=False, **box)
        if not unrestricted.converged or unrestricted.loglik < restricted.loglik-1e-5:
            unrestricted = sample.fit(start=restricted, multistart=True, **box)
        if (not unrestricted.converged or not restricted.converged
                or unrestricted.loglik < restricted.loglik-1e-5):
            return float('inf')
        return max(0., 2*(unrestricted.loglik-restricted.loglik))
    except (ValueError, ArithmeticError, FloatingPointError):
        return float('inf')


def _workers(workers):
    if isinstance(workers, bool) or not isinstance(workers, Integral) or workers < 0:
        raise ValueError('workers must be a nonnegative integer; 0 means every CPU the process may use')
    if workers == 0:
        return len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else (os.cpu_count() or 1)
    return int(workers)


def _one_thread():
    # Each process fits one small model at a time; nested BLAS threads only
    # oversubscribe the machine. The limit does not change any result.
    try:
        from threadpoolctl import threadpool_limits
    except ImportError:
        return
    threadpool_limits(1)


def _pool(workers):
    return ProcessPoolExecutor(workers, initializer=_one_thread) if workers > 1 else nullcontext()


def _settings(rho, n_boot, seed, alpha):
    alpha = _alpha(alpha)
    if isinstance(rho, bool) or not np.isscalar(rho) or not np.isfinite(rho) or not 0 <= rho < 1:
        raise ValueError('tested rho must be finite and lie in [0,1)')
    if isinstance(n_boot, bool) or not isinstance(n_boot, Integral) or n_boot < 19:
        raise ValueError('n_boot must be an integer at least 19')
    if isinstance(seed, bool) or not isinstance(seed, Integral) or seed < 0:
        raise ValueError('an explicit nonnegative integer seed is required')
    return alpha


def calibrate_null(problem, panels, exact, *, rho, n_boot=199, seed,
                   alpha=.05, order=24, full_fit=None, executor=None):
    """Calibrate the LR test at rho, profiling the mean, the total scale and
    the fixed effects. ``panels`` is the pair ``_panels`` returns or, for
    records read on one panel, its edge array."""
    alpha = _settings(rho, n_boot, seed, alpha)
    panels = _as_panels(panels, problem.n)
    fit = problem.fit(order=order) if full_fit is None else full_fit
    null = problem.fit(rho=float(rho), start=fit, order=order)
    if null.loglik > fit.loglik + 1e-5:
        fit = problem.fit(start=null, order=order, multistart=True)
    statistics = np.full(n_boot, np.inf)
    resolved = fit.converged and null.converged and null.loglik <= fit.loglik + 1e-5
    observed = max(0., 2*(fit.loglik-null.loglik)) if resolved else 0.
    rng = np.random.default_rng(seed)
    if resolved:
        # Draw every dataset first, in the original order, so the result does
        # not depend on how the fits are distributed across processes.
        jobs = []
        for _ in range(n_boot):
            effects = rng.normal(0, null.total_sd*np.sqrt(rho), problem.groups)
            latent = null.mean + _record_offset(problem, null) + effects[problem.code] + rng.normal(
                0, null.total_sd*np.sqrt(1-rho), problem.n)
            a, b = observe_panels(latent, panels)
            a[exact] = latent[exact]
            b[exact] = latent[exact]
            jobs.append((a, b, problem.code, float(rho), null, order, _record_covariate(problem)))
        if executor is None:
            values = map(_bootstrap_statistic, jobs)
        else:
            values = executor.map(_bootstrap_statistic, jobs,
                                  chunksize=max(1, n_boot//(4*getattr(executor, '_max_workers', 1))))
        statistics = np.array(list(values), dtype=float)
    rank = int(np.ceil((n_boot+1)*(1-alpha)))
    critical = float(np.sort(statistics)[rank-1]) if rank <= n_boot else float('inf')
    pvalue = float((1+np.count_nonzero(statistics >= observed))/(n_boot+1))
    reportable = bool(resolved and np.isfinite(critical))
    status = 'nullwise_parametric_bootstrap; plug_in_nuisance; empirical_calibration_required'
    if not reportable:
        status = 'unresolved_null_calibration; no_rejection'
        pvalue = 1.
    return NullBootstrapTest(float(rho), pvalue, float(observed), critical, alpha,
        bool(pvalue > alpha), fit, null, tuple(float(x) for x in statistics),
        int(np.isinf(statistics).sum()) if resolved else 0, reportable, status,
        int(n_boot) if resolved else 0)


def null_bootstrap_test(lo, hi, lineage, *, rho, panel_edges, n_boot=199,
                        seed, alpha=.05, order=24, workers=1, panel=None, covariate=None):
    """Test one population variance fraction with a declared observation panel.

    Reject when pvalue <= alpha. The (b+1)/(B+1) calculation includes ties and
    failed simulations conservatively. The nuisance parameters (mean, total
    scale and any fixed effects) are fitted under the tested rho; simulation is
    not performed at the unrestricted rho.
    ``workers`` > 1 distributes the bootstrap fits; results are unchanged.
    """
    _settings(rho, n_boot, seed, alpha)
    workers = _workers(workers)
    problem, edges, exact = _validated_problem(lo, hi, lineage, panel_edges, panel, covariate)
    with _pool(workers) as executor:
        return calibrate_null(problem, edges, exact, rho=rho, n_boot=n_boot,
                              seed=seed, alpha=alpha, order=order, executor=executor)


@dataclass(frozen=True)
class NullBootstrapInterval:
    low: float
    high: float
    confidence: float
    fit: MICFit
    tests: tuple[NullBootstrapTest, ...]
    method: str
    status: str
    reportable: bool
    tolerance: float

    def as_dict(self):
        return asdict(self)


def null_bootstrap_interval(lo, hi, lineage, *, panel_edges, n_boot=199, seed,
                            alpha=.05, order=24, tolerance=.002, workers=1, panel=None,
                            covariate=None):
    """Numerical hull of null-wise LR test inversion, with outer brackets.

    A fixed grid plus the MLE locates the outer crossings. Bisection retains
    the rejected outer brackets and uses common random numbers across rho.
    This is numerical test inversion, not a proof of profile connectedness.
    Every tested null and its bootstrap sample are retained in the result.
    ``workers`` > 1 distributes the bootstrap fits; results are unchanged.
    ``panel`` names the panel of every record when the records were read on
    several; ``panel_edges`` then maps each name to that panel's cut points,
    and every simulated reading is made on the panel of the record it stands
    for. ``covariate`` names, per record, the level of one categorical
    covariate whose fixed effects are estimated with the other parameters,
    re-estimated under every candidate rho and in every simulated dataset.
    """
    alpha = _settings(0., n_boot, seed, alpha)
    workers = _workers(workers)
    with _pool(workers) as executor:
        return _invert(lo, hi, lineage, panel_edges, n_boot, seed, alpha,
                       order, tolerance, executor, panel, covariate)


def _invert(lo, hi, lineage, panel_edges, n_boot, seed, alpha, order, tolerance, executor,
            panel=None, covariate=None):
    if isinstance(tolerance, bool) or not np.isfinite(tolerance) or not 1e-5 <= tolerance <= .05:
        raise ValueError('tolerance must be between 1e-5 and .05')
    problem, edges, exact = _validated_problem(lo, hi, lineage, panel_edges, panel, covariate)
    fit = problem.fit(order=order)
    cache = {}

    def result(low=0., high=1., reportable=False, status='unresolved_inversion; full_range'):
        return NullBootstrapInterval(float(low), float(high), 1-alpha, fit,
            tuple(cache[k] for k in sorted(cache)),
            'nullwise_parametric_bootstrap_LR_inversion', status, reportable,
            float(tolerance))

    def accepted(rho):
        # A certified limiting maximum includes rho=1 in the closure.
        # Do not simulate from zero or infinite nuisance-scale limits.
        if rho == 1. and fit.rho == 1.:
            return True
        key = float(rho)
        if key not in cache:
            cache[key] = calibrate_null(problem, edges, exact, rho=key,
                n_boot=n_boot, seed=seed, alpha=alpha, order=order, full_fit=fit,
                executor=executor)
        test = cache[key]
        if not test.reportable or test.unrestricted_fit.loglik > fit.loglik+1e-4:
            raise ArithmeticError('unresolved null or improved unrestricted maximum')
        return test.accepted

    if not fit.converged:
        return result()
    grid = np.unique(np.r_[0., .05, .1, .25, .5, .75, .9, .97, .995, fit.rho])
    try:
        inside = np.flatnonzero([accepted(r) for r in grid])
        if not len(inside):
            return result(status='no_accepted_grid_point; full_range')
        first, last = int(inside[0]), int(inside[-1])
        low = 0.
        if first > 0:
            left, right = float(grid[first-1]), float(grid[first])
            while right-left > tolerance:
                midpoint = (left+right)/2
                if accepted(midpoint):
                    right = midpoint
                else:
                    left = midpoint
            low = left
        high = 1.
        if last < len(grid)-1:
            left, right = float(grid[last]), float(grid[last+1])
            while right-left > tolerance:
                midpoint = (left+right)/2
                if accepted(midpoint):
                    left = midpoint
                else:
                    right = midpoint
            high = right
        return result(low, high, True,
            'nullwise_bootstrap_numerical_hull; plug_in_nuisance; Gaussian_model')
    except (ValueError, ArithmeticError, FloatingPointError):
        return result()
