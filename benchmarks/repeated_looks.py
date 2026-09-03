"""Operating characteristics of multiplicity rules read repeatedly over time.

WHY THIS EXISTS. The manuscript asserts, from theory alone, that the e-value
carries its guarantee through an unbounded number of inspections, so a
surveillance panel re-read every year keeps its error control, whereas a false
discovery rate procedure recomputed at unplanned looks controls nothing. Both
halves of that sentence are claims about operating characteristics, and neither
has been measured on this code. This module measures them. It simulates a panel
of antimicrobials followed for a fixed number of yearly collections, applies
four decision rules at every look, and records what a surveillance programme
actually experiences, which is not the error rate of one look but the chance
that some agent with no lineage signal was declared at some point during the
programme.

WHAT IS COMPARED. The first rule is the Benjamini-Hochberg step-up recomputed
on the whole accumulated panel at every look. The second is the same step-up
under the Benjamini-Yekutieli arbitrary-dependence correction, which a reviewer
will name because cross-resistance makes the agents dependent with unknown
sign. The third is the LOND rule of Javanmard and Montanari, an online
procedure designed for exactly this sequential setting, which is the fairest
comparator available and the one whose omission would be the first objection.
The fourth is the e-BH procedure applied to the e-process as the package ships
it, which recomputes the e-value on the accumulated cohort at every look. A
fifth, registered as a secondary arm, applies e-BH to the running product of
per-year e-values, which is the object the anytime-valid theory actually
describes, because a product of e-values over disjoint batches is a test
martingale and a value recomputed on overlapping accumulated data is not.

WHY THE GENERATOR IS NOT THE ESTIMATOR. The estimator is a one-way Bernoulli
model with empirical-Bayes shrinkage of group probabilities, scored by a split
likelihood ratio. The generator here is a latent liability threshold model with
random lineage effects, a shared isolate-level factor standing for co-carriage
of determinants on mobile elements, and unequal lineage sizes drawn from a
Dirichlet. No part of the estimator's likelihood appears in it. The lineage
effects and the noise are drawn from a Gaussian in the calibration arm and from
a t distribution on four degrees of freedom in the validation arm, standardised
by theoretical moments rather than by the sample, so that the two arms differ in
shape and not in scale and the between-sample variability the estimator exists
to quantify is left intact.
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import sys
import time
from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy import stats as sps

from amr_clonalshare import stats as bstats
from amr_clonalshare.evalues import e_process, e_bh, combine_independent

LAWS = ("gaussian", "t4")
RULES = ("bh_accumulated", "by_accumulated", "lond_batch", "ebh_accumulated",
         "ebh_product", "lond_finite_horizon")


def standardised_draws(rng: np.random.Generator, shape, law: str) -> np.ndarray:
    """Draws with mean zero and variance one, scaled by theoretical moments.

    Scaling a draw by its own sample mean and standard deviation would force
    every simulated panel to a sample mean of exactly zero and remove most of
    the randomness the procedures are being asked to control. The theoretical
    standard deviation of a t variate on nu degrees of freedom is the square
    root of nu over nu minus two, which is used here so that the heavy-tailed
    arm and the Gaussian arm differ in shape alone.
    """
    if law == "gaussian":
        return rng.standard_normal(shape)
    if law == "t4":
        nu = 4.0
        return rng.standard_t(nu, size=shape) / np.sqrt(nu / (nu - 2.0))
    raise ValueError(f"law must be one of {LAWS}, not {law!r}")


def simulate_panel(rng: np.random.Generator, *, n_agents: int,
                   null_fraction: float, n_lineages: int, n_per_year: int,
                   n_years: int, law: str, prevalence: float, rho: float,
                   shares: Sequence[float]) -> Dict[str, np.ndarray]:
    """One surveillance panel followed for ``n_years`` yearly collections.

    Every agent is a binary resistance trait measured on the same isolates, so
    the agents are dependent through the lineage composition they share and
    through an isolate-level factor that stands for determinants travelling
    together on a mobile element. An agent whose lineage share is zero is null
    by construction, because its liability then involves no term indexed by
    lineage, and this remains true whatever the shared factor does. The
    liability is standardised to unit variance for every share, so that
    changing the share moves the lineage signal and not the prevalence.
    """
    n_total = int(n_per_year) * int(n_years)
    freq = rng.dirichlet(np.full(int(n_lineages), 1.0))
    lineage = rng.choice(int(n_lineages), size=n_total, p=freq)
    year = np.repeat(np.arange(1, int(n_years) + 1), int(n_per_year))
    n_null = int(round(float(null_fraction) * int(n_agents)))
    share = np.zeros(int(n_agents), dtype=float)
    tail = np.arange(n_null, int(n_agents))
    if tail.size:
        share[tail] = np.asarray(shares, dtype=float)[
            np.arange(tail.size) % len(shares)]
    order = rng.permutation(int(n_agents))
    share = share[order]
    is_null = share <= 0.0
    u = standardised_draws(rng, (int(n_lineages), int(n_agents)), law)
    c = standardised_draws(rng, (n_total, 1), law)
    z = standardised_draws(rng, (n_total, int(n_agents)), law)
    mu = float(sps.norm.ppf(float(prevalence)))
    liab = (mu + np.sqrt(share) * u[lineage]
            + np.sqrt((1.0 - share) * float(rho)) * c
            + np.sqrt((1.0 - share) * (1.0 - float(rho))) * z)
    return {"Y": (liab > 0.0).astype(float), "lineage": lineage,
            "year": year, "is_null": is_null, "share": share}


def permutation_pvalues(Y: np.ndarray, lineage: np.ndarray, *, n_perm: int,
                        rng: np.random.Generator) -> np.ndarray:
    """Exact permutation p-values for lineage association, one per agent.

    The statistic is the between-lineage sum of squares of the trait, which is
    a monotone transform of the Pearson statistic for the lineage by trait
    table and needs no asymptotics. That matters here because a panel of forty
    lineages read after one year of thirty isolates has almost every cell
    empty, and an asymptotic p-value would fail at a single look for reasons
    that have nothing to do with repeated looking. The permutation p-value is
    exactly valid at every single look by construction, so any excess false
    discovery rate the study reports can only come from the looking. The
    correction of Phipson and Smyth is applied through the shipped estimator so
    that no p-value is reported as zero.

    Permuting the trait against the lineage labels is the same as randomly
    partitioning the isolates into groups of the observed sizes, so each
    permutation is one shuffle followed by a cumulative sum read at the group
    boundaries. All agents share the same shuffle, which is correct because the
    agents are measured on the same isolates and their joint dependence is part
    of what the multiplicity rules must handle.
    """
    Y = np.ascontiguousarray(np.asarray(Y, dtype=float))
    n, A = Y.shape
    codes, counts = np.unique(np.asarray(lineage), return_counts=True)
    code = np.searchsorted(codes, np.asarray(lineage))
    G = codes.size
    cnt = counts.astype(float)
    obs = np.empty((G, A))
    for a in range(A):
        obs[:, a] = np.bincount(code, weights=Y[:, a], minlength=G)
    stat_obs = (obs * obs / cnt[:, None]).sum(axis=0)
    edges = np.concatenate(([0], np.cumsum(counts)))
    csum = np.zeros((n + 1, A))
    null = np.empty((int(n_perm), A))
    for b in range(int(n_perm)):
        idx = rng.permutation(n)
        np.cumsum(Y[idx], axis=0, out=csum[1:])
        S = csum[edges[1:]] - csum[edges[:-1]]
        null[b] = (S * S / cnt[:, None]).sum(axis=0)
    return np.array([bstats.permutation_pvalue(null[:, a], float(stat_obs[a]),
                                               tail="greater")
                     for a in range(A)], dtype=float)


_LOND_TOTAL: Dict[int, float] = {}


def _lond_series_total(cut: int = 10 ** 7) -> float:
    """The sum of the LOND weight sequence over all hypotheses ever to come.

    Javanmard and Montanari require the sequence to sum to the nominal level
    over the infinite horizon, which is what makes the rule valid when the
    analyst does not know in advance how many hypotheses will arrive. The
    weight decays like one over ell times log squared ell, so the series
    converges but only at the rate of one over log, and truncating it at the
    length of this particular study would spend a budget the rule is not
    entitled to spend. The partial sum is taken to ten million terms and the
    remainder is added from the exact integral of the same envelope.
    """
    if cut in _LOND_TOTAL:
        return _LOND_TOTAL[cut]
    total = 0.0
    step = 10 ** 6
    for start in range(1, cut + 1, step):
        ell = np.arange(start, min(start + step, cut + 1), dtype=float)
        total += float((1.0 / (ell * np.log(np.maximum(ell, 2.0)) ** 2)).sum())
    total += 1.0 / np.log(float(cut))
    _LOND_TOTAL[cut] = total
    return total


def lond_levels(n_hyp: int, alpha: float, *, horizon: str) -> np.ndarray:
    """The LOND weight sequence, normalised over one of two horizons.

    The rule tests hypothesis j at level beta_j times one plus the number of
    discoveries already made, so the sequence beta fixes the whole procedure.
    The horizon named infinite is the one the theorem requires and is reported
    as the primary comparator. The horizon named finite normalises the same
    weights over the hypotheses this study will actually see, which is not a
    valid online rule because it presumes the stopping point, but it is
    reported as a registered secondary arm so that the primary result cannot be
    dismissed as an artefact of handicapping the comparator.
    """
    ell = np.arange(1, int(n_hyp) + 1, dtype=float)
    w = 1.0 / (ell * np.log(np.maximum(ell, 2.0)) ** 2)
    if horizon == "infinite":
        total = _lond_series_total()
    elif horizon == "finite":
        total = float(w.sum())
    else:
        raise ValueError(f"horizon must be 'infinite' or 'finite', "
                         f"not {horizon!r}")
    return float(alpha) * w / total


def lond_declare(p_batch: np.ndarray, alpha: float, *,
                 horizon: str) -> np.ndarray:
    """Apply LOND to the yearly batches in arrival order.

    Each agent at each look is one new hypothesis, tested on the p-value of
    that year's own isolates rather than on the accumulated cohort, because
    LOND is proved for a stream of hypotheses whose p-values are independent
    and reusing the earlier years would break that in the comparator's
    disfavour for a reason unrelated to the question being asked.
    """
    T, A = p_batch.shape
    beta = lond_levels(T * A, alpha, horizon=horizon)
    out = np.zeros((T, A), dtype=bool)
    discoveries = 0
    j = 0
    for t in range(T):
        for a in range(A):
            level = beta[j] * (discoveries + 1)
            if p_batch[t, a] <= level:
                out[t, a] = True
                discoveries += 1
            j += 1
    return out


def _safe_e(value: float) -> float:
    """An e-value that could not be computed is read as no evidence.

    A batch in which the trait is constant leaves the split likelihood ratio
    undefined rather than large, and the honest reading of an undefined ratio
    is one, which is the value that never contributes a rejection.
    """
    return 1.0 if not np.isfinite(value) or value < 0.0 else float(value)


def run_replicate(seed_seq: np.random.SeedSequence, *, n_agents: int,
                  null_fraction: float, n_lineages: int, n_per_year: int,
                  n_years: int, law: str, prevalence: float, rho: float,
                  shares: Sequence[float], alpha: float, n_perm: int,
                  folds: int, repeats: int) -> Dict[str, object]:
    """Simulate one panel, apply every rule at every look, return declarations.

    Two p-values and two e-values are formed at each look. The accumulated ones
    use every isolate collected so far and feed the two step-up procedures and
    the e-BH arm as the package ships it. The per-year ones use that year's
    isolates alone and feed the online comparator and the running product,
    which is the only combination of the shipped functions that forms a test
    martingale over the looks.
    """
    child = seed_seq.spawn(4)
    panel = simulate_panel(np.random.Generator(np.random.PCG64(child[0])),
                           n_agents=n_agents, null_fraction=null_fraction,
                           n_lineages=n_lineages, n_per_year=n_per_year,
                           n_years=n_years, law=law, prevalence=prevalence,
                           rho=rho, shares=shares)
    Y, lineage = panel["Y"], panel["lineage"]
    rng_p = np.random.Generator(np.random.PCG64(child[1]))
    ss_e = child[2]
    T, A = int(n_years), int(n_agents)
    p_acc = np.ones((T, A))
    p_bat = np.ones((T, A))
    e_acc = np.ones((T, A))
    e_bat = np.ones((T, A))
    for t in range(T):
        hi = (t + 1) * int(n_per_year)
        lo = t * int(n_per_year)
        p_acc[t] = permutation_pvalues(Y[:hi], lineage[:hi], n_perm=n_perm,
                                       rng=rng_p)
        p_bat[t] = permutation_pvalues(Y[lo:hi], lineage[lo:hi],
                                       n_perm=n_perm, rng=rng_p)
        seeds = ss_e.spawn(2 * A)
        for a in range(A):
            e_acc[t, a] = _safe_e(e_process(Y[:hi, a], lineage[:hi],
                                            folds=folds, repeats=repeats,
                                            seed=seeds[a]).e_value)
            e_bat[t, a] = _safe_e(e_process(Y[lo:hi, a], lineage[lo:hi],
                                            folds=folds, repeats=repeats,
                                            seed=seeds[A + a]).e_value)
    return _apply_rules(p_acc, p_bat, e_acc, e_bat, alpha=alpha,
                        is_null=panel["is_null"], prevalence_seen=
                        float(Y.mean()))


def _apply_rules(p_acc, p_bat, e_acc, e_bat, *, alpha, is_null,
                 prevalence_seen) -> Dict[str, object]:
    """Turn the four streams of evidence into per-look declaration matrices.

    The two step-up procedures and the two e-BH arms are recomputed from
    scratch at every look, so their declarations need not be nested across
    looks, and what a programme experiences is the union over looks. LOND
    accumulates by construction. Both readings are recorded so that the
    per-look rate and the rate over the whole programme can be reported side by
    side rather than one standing in for the other.
    """
    T, A = p_acc.shape
    decl: Dict[str, np.ndarray] = {}
    for name, pv, dep in (("bh_accumulated", p_acc, "independent"),
                          ("by_accumulated", p_acc, "arbitrary")):
        out = np.zeros((T, A), dtype=bool)
        for t in range(T):
            out[t] = bstats.benjamini_hochberg(pv[t], alpha,
                                               nan_policy="raise",
                                               dependence=dep)[1]
        decl[name] = out
    decl["lond_batch"] = lond_declare(p_bat, alpha, horizon="infinite")
    decl["lond_finite_horizon"] = lond_declare(p_bat, alpha, horizon="finite")
    out = np.zeros((T, A), dtype=bool)
    for t in range(T):
        out[t, e_bh(e_acc[t], alpha=alpha)["rejected"]] = True
    decl["ebh_accumulated"] = out
    out = np.zeros((T, A), dtype=bool)
    running = np.ones(A)
    for t in range(T):
        running = np.array([combine_independent((running[a], e_bat[t, a]))
                            for a in range(A)])
        out[t, e_bh(running, alpha=alpha)["rejected"]] = True
    decl["ebh_product"] = out
    n_null = int(is_null.sum())
    n_alt = int((~is_null).sum())
    rec: Dict[str, object] = {"prevalence_seen": prevalence_seen,
                              "n_null": n_null, "n_alt": n_alt}
    for name, d in decl.items():
        ever = d.any(axis=0)
        n_ever = int(ever.sum())
        per_look = np.array([
            (d[t] & is_null).sum() / max(1, int(d[t].sum())) for t in range(T)])
        rec[name] = {
            "fdp_ever": float((ever & is_null).sum() / max(1, n_ever)),
            "fdp_per_look": float(per_look.mean()),
            "fdp_final_look": float(per_look[-1]),
            "power_ever": float((ever & ~is_null).sum() / n_alt)
            if n_alt else float("nan"),
            "power_final_look": float((d[-1] & ~is_null).sum() / n_alt)
            if n_alt else float("nan"),
            "n_declared_ever": n_ever,
        }
    return rec


METRICS = ("fdp_ever", "fdp_per_look", "fdp_final_look", "power_ever",
           "power_final_look", "n_declared_ever")


def grid(null_fractions, lineages, per_year, laws) -> List[Dict[str, object]]:
    """The registered configuration grid, in a fixed and reproducible order.

    The order is fixed because each cell takes its seed from its position, so a
    reordering would silently change every number in the table.
    """
    cells = []
    for law in laws:
        for nf in null_fractions:
            for L in lineages:
                for m in per_year:
                    cells.append({"law": law, "null_fraction": float(nf),
                                  "n_lineages": int(L), "n_per_year": int(m)})
    return cells


def _worker(args):
    idx, cell, common = args
    ss = np.random.SeedSequence(entropy=common["master_seed"],
                                spawn_key=(int(cell["cell_index"]), int(idx)))
    return run_replicate(ss, null_fraction=cell["null_fraction"],
                         n_lineages=cell["n_lineages"],
                         n_per_year=cell["n_per_year"], law=cell["law"],
                         **{k: common[k] for k in
                            ("n_agents", "n_years", "prevalence", "rho",
                             "shares", "alpha", "n_perm", "folds", "repeats")})


def summarise(records: List[Dict[str, object]]) -> Dict[str, object]:
    """Monte-Carlo mean and standard error of every metric for every rule.

    The standard error is the sample standard deviation of the per-replicate
    value divided by the square root of the replicate count, which is the
    right quantity for a false discovery proportion because its denominator is
    itself random and a binomial formula would understate it.
    """
    R = len(records)
    out: Dict[str, object] = {"replicates": R,
                              "prevalence_seen": float(np.mean(
                                  [r["prevalence_seen"] for r in records]))}
    for rule in RULES:
        block = {}
        for metric in METRICS:
            v = np.array([float(r[rule][metric]) for r in records], dtype=float)
            good = v[np.isfinite(v)]
            block[metric] = {
                "mean": float(good.mean()) if good.size else float("nan"),
                "mc_se": float(good.std(ddof=1) / np.sqrt(good.size))
                if good.size > 1 else float("nan"),
                "n_finite": int(good.size)}
        out[rule] = block
    return out


def environment() -> Dict[str, object]:
    """Everything needed to say which machine and which library made a number."""
    import scipy
    import amr_clonalshare as btc
    return {"python": sys.version.split()[0], "numpy": np.__version__,
            "scipy": scipy.__version__,
            "amr_clonalshare": getattr(btc, "__version__", "unknown"),
            "package_file": btc.__file__, "platform": platform.platform(),
            "node": platform.node()}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cell-index", type=int, default=None)
    ap.add_argument("--reps", type=int, default=6000)
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--n-agents", type=int, default=13)
    ap.add_argument("--n-years", type=int, default=8)
    ap.add_argument("--n-perm", type=int, default=999)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--repeats", type=int, default=4)
    ap.add_argument("--prevalence", type=float, default=0.35)
    ap.add_argument("--rho", type=float, default=0.3)
    ap.add_argument("--master-seed", type=int, default=20260902)
    ap.add_argument("--out", type=str, required=True)
    ap.add_argument("--list-cells", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args(argv)

    cells = grid((1.0, 0.5, 0.8), (15, 40), (30, 60), LAWS)
    for i, c in enumerate(cells):
        c["cell_index"] = i
    if a.list_cells:
        print(json.dumps(cells, indent=1))
        return 0
    common = {"n_agents": a.n_agents, "n_years": a.n_years,
              "prevalence": a.prevalence, "rho": a.rho, "shares": (0.1, 0.3),
              "alpha": a.alpha, "n_perm": a.n_perm, "folds": a.folds,
              "repeats": a.repeats, "master_seed": a.master_seed}
    todo = cells if a.cell_index is None else [cells[a.cell_index]]
    results = []
    for cell in todo:
        t0 = time.time()
        payload = [(i, cell, common) for i in range(a.reps)]
        if a.procs > 1:
            import multiprocessing as mp
            with mp.get_context("fork").Pool(a.procs) as pool:
                recs = pool.map(_worker, payload, chunksize=4)
        else:
            recs = [_worker(x) for x in payload]
        s = summarise(recs)
        s.update({"cell": cell, "seconds": time.time() - t0,
                  "seconds_per_replicate": (time.time() - t0) / max(1, a.reps),
                  "config": common, "environment": environment()})
        results.append(s)
        print(f"cell {cell['cell_index']} done in {s['seconds']:.1f}s", flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump(results, fh, indent=1, sort_keys=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
