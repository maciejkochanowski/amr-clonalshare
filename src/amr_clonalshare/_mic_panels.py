"""Dilution panels as the MIC likelihood sees them.

A latent log2 value read on a panel becomes the interval between the two cut
points it fell between, open at either end of the panel. Records may have been
read on different panels, so a panel index is carried per record. The checked
problem built here is shared by the exact pivot and the null-wise bootstrap.
"""
from __future__ import annotations

import numpy as np

from ._mic_likelihood import GaussianMICLikelihood


def _panel(edges, minimum=2):
    if edges is None:
        raise ValueError('the tested panel must be supplied explicitly')
    e=np.asarray(edges,dtype=float)
    if e.ndim!=1 or e.size<minimum or not np.isfinite(e).all() or np.any(np.diff(e)<=0):
        raise ValueError(f'panel edges must be at least {minimum} strictly increasing finite cut points')
    return e.copy()


def observe_panel(values, edges):
    """Return intervals (-infinity,e0], (e0,e1], ..., (elast,infinity)."""
    return _observe(values, _panel(edges))


def _observe(values, e):
    y=np.asarray(values,dtype=float)
    if not np.isfinite(y).all():
        raise ValueError('latent observations must be finite')
    k=np.searchsorted(e,y,side='left')
    return np.r_[-np.inf,e][k], np.r_[e,np.inf][k]


def _panels(panel_edges, panel, n):
    """One edge array per panel and the panel index of every record.

    ``panel`` is ``None`` when every record was read on the one panel
    ``panel_edges`` describes, and otherwise names the panel of each record,
    with ``panel_edges`` mapping every name to that panel's cut points.
    """
    if panel is None:
        return [_panel(panel_edges)], np.zeros(n, dtype=int)
    labels = np.asarray([str(x) for x in np.asarray(panel, dtype=object).ravel()])
    if labels.size != n:
        raise ValueError('panel must name the panel of every record')
    if not isinstance(panel_edges, dict):
        raise ValueError('with per-record panels, panel_edges must map each panel name to its cut points')
    declared = {str(k): v for k, v in panel_edges.items()}
    names = sorted(set(labels))
    missing = [name for name in names if name not in declared]
    if missing:
        raise ValueError(f'no cut points declared for panel(s) {missing}')
    # One panel may be a single cut point, a binary reading of its records,
    # as long as another panel carries the dilution steps that identify the scale.
    edges = [_panel(declared[name], minimum=1) for name in names]
    if max(e.size for e in edges) < 2:
        raise ValueError('at least one panel needs two or more cut points')
    index = np.searchsorted(names, labels)
    return edges, index


def _as_panels(panels, n):
    """Accept either the pair ``_panels`` returns or one bare edge array."""
    if isinstance(panels, tuple) and len(panels) == 2:
        return panels
    return _panels(panels, None, n)


def observe_panels(values, panels):
    """``observe_panel`` with each record read on its own panel."""
    edges, index = _as_panels(panels, len(values))
    y = np.asarray(values, dtype=float)
    a, b = np.empty_like(y), np.empty_like(y)
    for k, e in enumerate(edges):
        m = index == k
        a[m], b[m] = _observe(y[m], e)
    return a, b


def _validated_problem(lo,hi,lineage,panel_edges,panel=None,covariate=None):
    problem=GaussianMICLikelihood(lo,hi,lineage,covariate=covariate)
    panels=_panels(panel_edges,panel,problem.n)
    exact=np.isfinite(problem.lo)&(problem.lo==problem.hi)
    for k,e in enumerate(panels[0]):
        allowed=set(zip(np.r_[-np.inf,e],np.r_[e,np.inf]))
        here=~exact&(panels[1]==k)
        if any((a,b) not in allowed for a,b in zip(problem.lo[here],problem.hi[here])):
            raise ValueError('observed intervals do not match the declared panel')
    return problem,panels,exact


def _record_offset(problem, fit):
    """The fixed effect of every reading's covariate level under ``fit``."""
    if not problem.n_coefficients:
        return 0.
    betas = problem._split_coefficients(fit.coefficients)
    return np.sum([beta[problem.record_level[:, k]] for k, beta in enumerate(betas)], axis=0)


def _record_covariate(problem):
    """The covariate levels of every reading, by their original names so that a
    simulated dataset is coded in the same order as the observed one; None
    without a covariate."""
    if not problem.n_coefficients:
        return None
    return [np.asarray(problem.levels[k], dtype=object)[problem.record_level[:, k]]
            for k in range(problem.record_level.shape[1])]
