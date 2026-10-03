"""Coverage and level of the lineage share of the MIC ordering.

Every replicate draws lineages, their readings and the laboratory of every
reading, reads them on the laboratory's panel, and estimates the share with
:func:`amr_clonalshare.mic_order.mic_order_share` as the pipeline does
(scored and permuted within laboratory, the pipeline's resampling settings).
The target is the share of the same functional for the represented lineages: every
lineage is read again ``TRUTH_DRAWS`` times, each laboratory's distribution
function is the isolate-weighted one of those readings at its cut points, and
the share is the isolate-weighted between-lineage variance of the scores over
their total variance. The truth is written here from that definition; it
shares no code with the package's scoring (the nonparametric
maximum-likelihood distribution of the readings, which the estimate uses).

Every replicate also records the interval for the represented lineages (the
studentized bootstrap of attribution._smoothed_bootstrap) against the same
target, and the bounds on the lineage share of the latent MIC ordering
(latent_order.py) against two targets. The first is the share of the ordering
of the latent log2 MICs of the sampled isolates themselves, each placed by the
share of the repeated isolates of its laboratory below it: every ordering
consistent with the readings lies in the set the bounds are drawn from
(Theorem 1), so that share must lie within the bounds in every replicate,
without tolerance for sampling. The second is the share of the represented lineages:
the same lineages' further draws of the latent log2 MIC, each placed by the
isolate-weighted distribution of those draws within its laboratory, and the
isolate-weighted between-lineage variance of those positions times 12, the
target of the one-sided lower confidence limit. Both are computed here
without the readings, the panels or any of the package's code.

Designs
-------
``gauss_*``        Gaussian lineage effects and residuals on the wide panel,
                   share of the latent variance 0.1, 0.5 and 0.9;
``null``           no lineage effect;
``heavy_censoring`` a three-well panel with most readings on an end well;
``few``            eight lineages; ``singleton_rich``: half the lineages
                   singletons;
``clone``          a two-point distribution of lineage effects (one lineage in five
                   carries the effect);
``t4``, ``contaminated``, ``two_component``: residuals with heavy tails, a 5 %
                   outlying component, and two modes within every lineage;
``wt_nwt``         readings in two modes, the upper one reached by lineage
                   membership (the *S. suis* lineage sizes);
``informative``    lineage size ordered as the lineage effect;
``two_labs``       two laboratories on different panels, 1 log2 apart, with
                   85 % of a lineage's readings in its home laboratory.
``pairs_*``        half the lineages of two isolates, the rest uneven, share
                   of the latent variance 0.5 and 0.9;
``all_pairs_*``    every lineage of two isolates, 10 and 30 lineages, share
                   0.7 (CONFIRMATORY_PROTOCOL.md, Study 2);
``null_lab_*``     no lineage effect, laboratories that differ in offset,
                   panel or spread, lineages aligned with laboratories (85 % or
                   all readings in the home laboratory): the level of the
                   permutation test within laboratories, beside the test that
                   ignores the laboratory.

One JSON list per run, one row per replicate.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from amr_clonalshare.mic_order import mic_order_share

ROOT_SEED = 20261002207
TRUTH_DRAWS = 16000
SSUIS = np.array((161, 60, 57, 50, 44, 42, 40, 36, 31, 22, 20, 17, 13, 11, 11, 10, 8, 7, 6, 5,
                  5, 4, 4, 3, 3, 2, 2, 1, 1, 1))
WIDE, NARROW, SHIFTED = np.arange(-4., 5.), np.arange(-1., 2.), np.arange(-2., 7.)


def _design(name, rho=.5, groups=30, sizes="uneven", effects="gaussian", residual="normal",
            labs=((0., WIDE, 1.),), alignment=1., mean=0., replicates=5000):
    return dict(design=name, rho=rho, groups=groups, sizes=sizes, effects=effects,
                residual=residual, labs=labs, alignment=alignment, mean=mean,
                replicates=replicates)


def order_design_grid():
    two = ((-.5, WIDE, 1.), (.5, NARROW, 1.))
    grid = [
        _design("gauss_0.1", rho=.1), _design("gauss_0.5"), _design("gauss_0.9", rho=.9),
        _design("null", rho=0., replicates=20000),
        _design("heavy_censoring", labs=((0., NARROW, 1.),), mean=-1.3),
        _design("few", groups=8), _design("singleton_rich", sizes="singleton_rich"),
        _design("clone", effects="two_point"), _design("t4", residual="t4"),
        _design("contaminated", residual="contaminated"),
        _design("two_component", residual="two_component"),
        _design("wt_nwt", rho=.7, sizes="ssuis", residual="wt_nwt"),
        _design("informative", sizes="informative"),
        _design("two_labs", labs=two, alignment=.85),
        _design("null_lab_offset_85", rho=0., labs=two, alignment=.85, replicates=20000),
        _design("null_lab_offset_100", rho=0., labs=two, replicates=20000),
        _design("null_lab_panel_100", rho=0., labs=((0., WIDE, 1.), (1.5, SHIFTED, 1.)), replicates=20000),
        _design("null_lab_spread_85", rho=0., labs=((0., WIDE, .6), (0., WIDE, 1.8)), alignment=.85,
                replicates=20000),
        _design("null_lab_few_100", rho=0., groups=8, labs=two, replicates=20000),
        _design("pairs_0.5", sizes="pairs"), _design("pairs_0.9", rho=.9, sizes="pairs"),
        _design("all_pairs_10", rho=.7, groups=10, sizes="two"),
        _design("all_pairs_30", rho=.7, sizes="two"),
    ]
    for k, cell in enumerate(grid):
        cell["cell"] = k
    return grid


def _sizes(cell, rng):
    g = cell["groups"]
    if cell["sizes"] == "ssuis":
        return SSUIS.copy()
    if cell["sizes"] == "singleton_rich":
        return np.r_[np.ones(g // 2, dtype=int), np.full(g - g // 2, 20)]
    if cell["sizes"] == "three":
        return np.full(g, 3)
    if cell["sizes"] == "two":
        return np.full(g, 2)
    if cell["sizes"] == "pairs":
        rest = np.rint(np.exp(rng.normal(2.2, 1., g - g // 2))).astype(int).clip(3, 100)
        return np.r_[np.full(g // 2, 2), rest]
    return np.rint(np.exp(rng.normal(2.2, 1., g))).astype(int).clip(2, 100)


def _effects(cell, rng, g):
    rho = cell["rho"]
    if cell["effects"] == "two_point":
        z = (rng.random(g) < .2).astype(float)
        return (z - .2) / np.sqrt(.16) * np.sqrt(rho)
    return rng.normal(0, np.sqrt(rho), g)


def _residual(cell, rng, n):
    kind = cell["residual"]
    if kind == "t4":
        return rng.standard_t(4, n) / np.sqrt(2.)
    if kind == "contaminated":
        return rng.normal(size=n) * np.where(rng.random(n) < .05, 5., 1.) / np.sqrt(2.2)
    if kind == "two_component":
        shifted = rng.random(n) < .3
        return (rng.normal(size=n) + 2.5 * shifted - .75) / np.sqrt(1. + .21 * 6.25)
    return rng.normal(size=n)


def _readings(cell, g, n, effects, home, extra, rng):
    """Latent log2 MIC and laboratory of n readings of lineage g."""
    rho = cell["rho"]
    lab = np.where(rng.random(n) < cell["alignment"], home[g], 1 - home[g]) if len(cell["labs"]) > 1 \
        else np.zeros(n, dtype=int)
    offset = np.array([lb[0] for lb in cell["labs"]])[lab]
    spread = np.array([lb[2] for lb in cell["labs"]])[lab]
    if cell["residual"] == "wt_nwt":
        # a lower and an upper mode; membership of the upper mode follows the
        # lineage liability, the position within a mode a lineage shift
        upper = effects[g] + rng.normal(0, np.sqrt(1 - rho), n) > .5244
        return cell["mean"] - 1. + 5. * upper + extra[g] + rng.normal(0, .8, n), lab
    y = cell["mean"] + offset + spread * (effects[g] + np.sqrt(1 - rho) * _residual(cell, rng, n))
    return y, lab


def _observe(y, edges):
    """Readings (lo, hi] on a panel with the given cut points, written here."""
    k = np.searchsorted(edges, y, side="left")
    ext = np.r_[-np.inf, edges, np.inf]
    return ext[k], ext[k + 1]


def generate(cell, replicate):
    rng = np.random.default_rng(np.random.SeedSequence([ROOT_SEED, cell["cell"], replicate]))
    sizes = _sizes(cell, rng)
    g = sizes.size
    effects = _effects(cell, rng, g)
    extra = rng.normal(0, .3, g)
    if cell["sizes"] == "informative":
        sizes = np.sort(sizes)[np.argsort(np.argsort(effects))]
    home = rng.integers(0, 2, g)
    parts = [_readings(cell, k, sizes[k], effects, home, extra, rng) for k in range(g)]
    y = np.concatenate([p[0] for p in parts])
    lab = np.concatenate([p[1] for p in parts])
    lineage = np.repeat(np.arange(g), sizes)
    lo, hi = np.empty(y.size), np.empty(y.size)
    for k, (_, edges, _) in enumerate(cell["labs"]):
        m = lab == k
        lo[m], hi[m] = _observe(y[m], edges)
    return dict(lo=lo, hi=hi, lineage=lineage, lab=lab, sizes=sizes, effects=effects, extra=extra,
                home=home, latent=y)


def sample_latent_share(data):
    """12 sum_g w_g (mean position of lineage g - 1/2)^2 over the repeated
    lineages, a position being the share of the repeated isolates of the
    same laboratory below the latent value plus half the share at it."""
    lineage, lab, y = data["lineage"], data["lab"], data["latent"]
    keep = np.bincount(lineage)[lineage] >= 2
    lineage, lab, y = lineage[keep], lab[keep], y[keep]
    u = np.empty(y.size)
    for k in np.unique(lab):
        m = np.flatnonzero(lab == k)
        order = np.argsort(y[m], kind="stable")
        rank = np.empty(m.size)
        rank[order] = np.arange(m.size)
        u[m] = (rank + .5) / m.size
    code = np.unique(lineage, return_inverse=True)[1]
    size = np.bincount(code).astype(float)
    means = np.bincount(code, weights=u) / size
    return float(12. * np.sum(size / size.sum() * (means - .5) ** 2))


def truth(cell, data, draws=TRUTH_DRAWS):
    """The share of the functional for the represented lineages, from ``draws``
    readings of every repeated lineage. Without a lineage effect every
    lineage's readings in a laboratory follow that laboratory's distribution, so every
    lineage's mean score is one half and the share is zero exactly."""
    if cell["rho"] == 0.:
        return 0.
    rng = np.random.default_rng(np.random.SeedSequence([ROOT_SEED, cell["cell"], 10**6]))
    sizes = data["sizes"]
    repeated = np.flatnonzero(sizes >= 2)
    w = sizes[repeated] / sizes[repeated].sum()
    ys, labs = [], []
    for g in repeated:
        y, lab = _readings(cell, g, draws, data["effects"], data["home"], data["extra"], rng)
        ys.append(y)
        labs.append(lab)
    y, lab = np.concatenate(ys), np.concatenate(labs)
    weight = np.repeat(w, draws)
    score = np.empty(y.size)
    for k, (_, edges, _) in enumerate(cell["labs"]):
        m = lab == k
        if not m.any():
            continue
        cdf = np.r_[0., [np.sum(weight[m] * (y[m] <= e)) / weight[m].sum() for e in edges], 1.]
        idx = np.searchsorted(edges, y[m], side="left")
        score[m] = (cdf[idx] + cdf[idx + 1]) / 2.
    per = score.reshape(len(repeated), draws)
    means, variances = per.mean(1), per.var(1)
    centre = np.sum(w * means)
    between = np.sum(w * (means - centre) ** 2)
    total = between + np.sum(w * variances)
    return float(between / total) if total > 0 else 0.


def latent_truth(cell, data, draws=TRUTH_DRAWS):
    """12 times the between-lineage variance of U = F_lab(latent MIC) for the
    represented lineages of two or more isolates, from ``draws`` latent draws of each, with F_lab
    the isolate-weighted distribution of those draws in the laboratory."""
    if cell["rho"] == 0.:
        return 0.
    rng = np.random.default_rng(np.random.SeedSequence([ROOT_SEED, cell["cell"], 10**6, 17]))
    sizes = data["sizes"]
    repeated = np.flatnonzero(sizes >= 2)
    w = sizes[repeated] / sizes[repeated].sum()
    ys, labs = [], []
    for g in repeated:
        y, lab = _readings(cell, g, draws, data["effects"], data["home"], data["extra"], rng)
        ys.append(y)
        labs.append(lab)
    y, lab = np.concatenate(ys), np.concatenate(labs)
    weight = np.repeat(w, draws)
    u = np.empty(y.size)
    for k in np.unique(lab):
        m = lab == k
        order = np.argsort(y[m], kind="stable")
        cum = np.cumsum(weight[m][order])
        pos = np.empty(int(m.sum()))
        pos[order] = (cum - weight[m][order] / 2.) / cum[-1]
        u[m] = pos
    means = u.reshape(repeated.size, draws).mean(1)
    return float(12. * np.sum(w * (means - .5) ** 2))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cell", type=int, required=True)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    grid = order_design_grid()
    if not 0 <= args.cell < len(grid) or args.start < 0 or args.replicates < 1:
        parser.error("invalid simulation arguments")
    if args.output.exists():
        raise FileExistsError(args.output)
    cell = grid[args.cell]
    rows = []
    for rep in range(args.start, args.start + args.replicates):
        started = time.perf_counter()
        data = generate(cell, rep)
        strata = data["lab"] if len(cell["labs"]) > 1 else None
        seed = int(np.random.SeedSequence([ROOT_SEED, cell["cell"], rep, 7]).generate_state(1)[0])
        r = mic_order_share(data["lo"], data["hi"], data["lineage"], strata=strata,
                            stratified_by="laboratory" if strata is not None else None, seed=seed)
        target = truth(cell, data)
        row = dict(cell=cell["cell"], design=cell["design"], replicate=rep, rho=cell["rho"],
                   n=int(data["lineage"].size), groups=int(data["sizes"].size), truth=target,
                   estimate=r["kappa_adj"], p_value=r["p_value"], estimable=bool(r["estimable"]),
                   end_wells=r["share_end_wells"],
                   observed_low=r["observed_low"], observed_high=r["observed_high"],
                   observed_covered=(bool(r["observed_low"] <= target <= r["observed_high"])
                                     if r["estimable"] and np.isfinite(r["observed_low"]) else None))
        latent = latent_truth(cell, data)
        sample = sample_latent_share(data)
        row.update(latent_truth=latent, latent_sample=sample, latent_lower=r["latent_order_lower"],
                   latent_upper=r["latent_order_upper"], latent_upper_bound=r["latent_order_upper_bound"],
                   latent_upper_exact=bool(r["latent_order_upper_exact"]), latent_sharp=bool(r["latent_order_sharp"]),
                   latent_limit=r["latent_order_lower_limit"],
                   latent_midpoint=r["latent_order_midpoint"], resolution=r["panel_resolution"],
                   sample_within_bounds=bool(r["latent_order_lower"] - 1e-9 <= sample
                                             <= r["latent_order_upper_bound"] + 1e-9),
                   latent_limit_below=bool(r["latent_order_lower_limit"] <= latent))
        if strata is not None:
            pooled = mic_order_share(data["lo"], data["hi"], data["lineage"], seed=seed)
            row.update(pooled_estimate=pooled["kappa_adj"], pooled_p_value=pooled["p_value"])
        row["seconds"] = time.perf_counter() - started
        rows.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output.with_suffix(".partial")
    partial.write_text(json.dumps(rows, indent=1, allow_nan=True, default=float))
    partial.replace(args.output)


if __name__ == "__main__":
    main()
