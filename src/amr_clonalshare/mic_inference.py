"""Population MIC variance-fraction inference, separate from the moment heuristic.

The exact-reading interval inverts a Gaussian F pivot for arbitrary fixed
lineage sizes. It does not use shrunken group means or effective group sizes.
The interval-reading likelihood/bootstrap route is implemented separately.
Neither route corrects informative sampling or unrecorded host dependence.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from scipy.optimize import brentq
from scipy.stats import f

from .attribution import _codes, _is_untyped


@dataclass(frozen=True)
class MICInterval:
    low: float
    high: float
    confidence: float
    n: int
    groups: int
    method: str
    status: str
    raw_set_empty: bool = False

    def as_dict(self):
        return asdict(self)


def _alpha(alpha):
    if isinstance(alpha, bool) or not np.isscalar(alpha) or not np.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError('alpha must be a finite number strictly between zero and one')
    return float(alpha)


def _exact_summary(values, lineage):
    y = np.asarray(values, dtype=float)
    labels = np.asarray(lineage, dtype=object)
    if y.ndim != 1 or labels.ndim != 1 or y.shape != labels.shape or not np.isfinite(y).all():
        raise ValueError('finite exact values and lineage labels must be matching one-dimensional arrays')
    if any(_is_untyped(label) for label in labels):
        raise ValueError('untyped lineage labels must be handled explicitly before inference')
    if y.size < 3:
        raise ValueError('at least three exact observations are required')
    code = _codes(labels)
    counts = np.bincount(code).astype(float)
    groups = counts.size
    if groups < 2 or y.size <= groups:
        raise ValueError('at least two lineages and within-lineage replication are required')
    # Center first to avoid loss of significance for large location shifts.
    centered = y - y.mean()
    means = np.bincount(code, weights=centered) / counts
    sse = float(np.sum((centered - means[code])**2))
    if not np.isfinite(sse) or sse <= 0:
        raise ValueError('positive within-lineage variation is required by the Gaussian F pivot')
    return counts, means, sse, int(y.size)


def _pivot(counts, means, sse, n, rho):
    if rho == 1:
        return 0.0
    weights = counts*(1-rho)/(1-rho+counts*rho)
    mu = float(np.dot(weights, means)/weights.sum())
    q = float(np.dot(weights, (means-mu)**2))
    return q/(counts.size-1)/(sse/(n-counts.size))


def exact_f_statistic(values, lineage, rho: float) -> float:
    """F pivot at a candidate rho, with nuisance mean and scale eliminated.

    Under rho in [0,1), independent Gaussian random intercepts/residuals and
    fixed group sizes, this follows F(G-1,N-G). At rho=1 this returns the
    limiting statistic, not a claim for a singular residual model.
    """
    if isinstance(rho, bool) or not np.isfinite(rho) or not 0 <= rho <= 1:
        raise ValueError('rho must lie between zero and one')
    counts, means, sse, n = _exact_summary(values, lineage)
    return _pivot(counts, means, sse, n, float(rho))


def exact_gaussian_interval(values, lineage, *, alpha: float = .05) -> MICInterval:
    """Invert the equal-tail generalized F pivot for population rho.

    For unequal sizes, the weighting depends on the candidate rho. Both
    limits are roots of the pivot; replacing the sizes by a mean is avoided.
    If the unrestricted confidence set lies below rho=0, return [0,0] and
    raw_set_empty=True. This closure is conservative at the boundary and
    does not change coverage for interior rho. It is not an absence-of-effect
    claim or an identification interval for non-Gaussian data.
    """
    alpha = _alpha(alpha)
    counts, means, sse, n = _exact_summary(values, lineage)
    quantiles = f.ppf([alpha/2, 1-alpha/2], counts.size-1, n-counts.size)
    start = _pivot(counts, means, sse, n, 0.)

    def root(target):
        if start <= target:
            return 0.0
        return float(brentq(lambda rho: _pivot(counts, means, sse, n, rho)-target,
                            0., 1., xtol=1e-12, rtol=1e-12))

    return MICInterval(root(quantiles[1]), root(quantiles[0]), 1-alpha, n,
                       int(counts.size), 'exact_generalized_F',
                       'Gaussian_model_exact; fixed_noninformative_group_sizes',
                       bool(start < quantiles[0]))
