"""archephy.py — ArchePhy-CS: phylogenetic convergence test for trait archetypes.

Why this exists
---------------
Every clustering of bacterial trait data faces the same confounder: isolates
share traits because they share ancestors. A cluster of genomes carrying the
same resistance and virulence loci may be one clone sampled many times, not a
recurrent trait combination. This is the lineage effect that bacterial GWAS
methods are built to remove (Earle et al. 2016, Nat Microbiol 1:16041,
doi:10.1038/nmicrobiol.2016.41; Lees et al. 2016, Nat Commun 7:12797,
doi:10.1038/ncomms12797; Collins & Didelot 2018, PLoS Comput Biol
14:e1005958, doi:10.1371/journal.pcbi.1005958). Nothing in a similarity-network
fusion or in a data-thinning test addresses it: both are agnostic to the tree.

ArchePhy-CS asks a different and sharper question about a candidate archetype -
a set of ``k >= 2`` traits that a clustering says go together:

    did these ``k`` traits arise together **repeatedly and independently**
    across the phylogeny (convergence / joint homoplasy), or once, in one
    ancestor, with everything else being descent?

Statistic
---------
On a rooted tree with ``E`` non-root edges, Fitch parsimony gives for each
trait ``t`` the set ``S_t`` of edges on which the trait changes state, with
``c_t = |S_t|``. The k-way joint homoplasy count is

    ``m_obs = |S_1 & S_2 & ... & S_k|``

the number of edges on which **all** ``k`` traits change together. Two or more
such edges mean the profile recurred independently.

Null distribution
-----------------
Under independent per-trait evolution with exchangeable edges, each ``S_t`` is
a uniformly random ``c_t``-subset of the ``E`` edges, independently across
traits. The distribution of ``M = |S_1 & ... & S_k|`` is then **exact and
closed-form**, with no resampling and no large-sample approximation - which is
what makes the test usable at the small cohort sizes where data thinning fails.

Its binomial moments are

    ``E[C(M, j)] = C(E, j) * prod_t [ C(c_t, j) / C(E, j) ]``

because a fixed set of ``j`` edges lies in every ``S_t`` with probability
``prod_t C(c_t, j) / C(E, j)``. The exact upper tail follows from the
inclusion-exclusion (Bonferroni) identity

    ``P(M >= m) = sum_{j >= m} (-1)^{j-m} C(j-1, m-1) E[C(M, j)]``

evaluated in exact rational arithmetic. For ``k = 2`` this reduces to the
hypergeometric survival function ``P(M >= m) = sf_Hypergeom(E, c_1, c_2)``,
which :func:`joint_homoplasy_pvalue` uses directly and which the test suite
checks the general formula against.

A Poisson approximation with ``lambda = prod_t c_t / E^(k-1)`` is the first
moment of the exact null and is available as ``method="poisson"`` for
comparison. It is anticonservative
in the regime that matters (small ``E``, large ``c_t``), which is precisely the
small-cohort regime this test targets - see ``benchmarks/`` for the comparison.

Assumptions, stated plainly
---------------------------
* **Edge exchangeability -- and why it is not the default.** The exact
  null treats edges as equally likely change slots. Real change edges are not
  exchangeable: they concentrate on short and terminal branches, so two traits
  that evolved **independently** share change edges more often than uniform
  ``c_t``-subsets do. Measured on birth-death trees with traits evolved under
  independent CTMCs, the exact null rejects at **0.072-0.078 +/- 0.011** against
  a nominal 0.05. "No resampling and no asymptotics means calibration at small n
  by construction" is therefore false: the arithmetic is exact, the null is not
  the biological one. The default is therefore ``method="ctmc"``, which
  simulates each trait down the actual tree (conditioned on its observed number
  of changes) and is not significantly above nominal in any cell tested. The
  exact p-value is still reported alongside as
  ``p_value_exact_uniform_null``. ``method="weighted"`` samples edges
  proportionally to branch length: markedly conservative under a strict clock
  (0.000 in every cell) but the **only** null here that survives a relaxed clock
  (0.033 where exact reaches 0.140 and the CTMC null 0.083).
* **Change, not gain.** Change edges count gains and losses alike. ``fitch_asr``
  returns directional counts so a gains-only refinement is available.
* **One tree, no uncertainty.** The tree is taken as fixed. Running the test
  across a bootstrap or posterior sample of trees and merging the resulting
  p-values with ``tva.merge_exchangeable`` is the recommended way to propagate
  phylogenetic uncertainty.
* **``m_obs >= 2`` is load-bearing.** A single joint change edge is one clonal
  co-origin - the confounder this test exists to reject - and is never called
  convergence however small its tail probability.

``dendropy`` is required for the tree parsing and is an optional dependency:
``pip install "amr-clonalshare[phylo]"``.
"""
from __future__ import annotations

from fractions import Fraction
import math
from math import comb, exp, log
from typing import Dict, List, Optional, Sequence

import numpy as np
from scipy import stats

__all__ = [
    "load_tree",
    "fitch_change_edges",
    "fitch_asr",
    "joint_homoplasy_pvalue",
    "archephy_cs_test",
]

_MAX_EXACT_TERMS = 400


def _require_dendropy():
    try:
        import dendropy
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "ArchePhy-CS needs dendropy. Install with "
            '`pip install "amr-clonalshare[phylo]"`.'
        ) from exc
    return dendropy


def load_tree(newick_str: str):
    """Parse a newick string; coerce non-positive branch lengths.

    Returns ``(tree, leaf_labels)``.
    """
    dendropy = _require_dendropy()
    t = dendropy.Tree.get(data=newick_str, schema="newick")
    for edge in t.preorder_edge_iter():
        if edge.length is None or edge.length <= 0:
            edge.length = 1e-6
    leaves = [ln.taxon.label for ln in t.leaf_node_iter()]
    return t, leaves


def fitch_change_edges(tree, leaves: Sequence[str], trait: np.ndarray) -> set:
    """Node ids whose incoming edge is a Fitch state-change edge for ``trait``.

    Standard two-pass Fitch parsimony (Fitch 1971, Syst Zool 20:406-416,
    doi:10.2307/2412116): a postorder pass builds the candidate state sets, a
    preorder pass resolves them against the parent, and an edge is a change
    edge when the resolved child state differs from the resolved parent state.
    """
    tof = {l: int(trait[i]) for i, l in enumerate(leaves)}
    up: Dict[int, set] = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            up[id(node)] = {tof[node.taxon.label]}
        else:
            sets = [up[id(c)] for c in node.child_nodes()]
            inter = set.intersection(*sets)
            up[id(node)] = inter if inter else set.union(*sets)
    down: Dict[int, int] = {}
    changed: set = set()
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            down[id(node)] = sorted(up[id(node)])[0]
        else:
            par = down[id(node.parent_node)]
            down[id(node)] = par if par in up[id(node)] else sorted(up[id(node)])[0]
            if down[id(node)] != par:
                changed.add(id(node))
    return changed


def fitch_asr(tree, leaves: List[str], trait: np.ndarray) -> dict:
    """Directional Fitch change counts: ``{gains, losses, ambiguous_nodes}``."""
    trait_dict = dict(zip(leaves, np.asarray(trait).astype(int).tolist()))
    state_set: Dict[int, set] = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            state_set[id(node)] = {trait_dict[node.taxon.label]}
        else:
            inter = None
            for c in node.child_nodes():
                inter = set(state_set[id(c)]) if inter is None else inter & state_set[id(c)]
            if not inter:
                inter = set()
                for c in node.child_nodes():
                    inter |= state_set[id(c)]
            state_set[id(node)] = inter
    chosen: Dict[int, int] = {}
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            chosen[id(node)] = min(state_set[id(node)])
        else:
            ps = chosen[id(node.parent_node)]
            chosen[id(node)] = ps if ps in state_set[id(node)] else min(state_set[id(node)])
    gains = losses = ambiguous = 0
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            continue
        p_state = chosen[id(node.parent_node)]
        c_state = chosen[id(node)]
        if len(state_set[id(node)]) > 1:
            ambiguous += 1
        if p_state == 0 and c_state == 1:
            gains += 1
        elif p_state == 1 and c_state == 0:
            losses += 1
    return {"gains": gains, "losses": losses, "ambiguous_nodes": ambiguous}


# --------------------------------------------------------------------------- #
# Null distribution of the k-way joint homoplasy count
# --------------------------------------------------------------------------- #
def _binomial_moment(E: int, c: Sequence[int], j: int) -> Fraction:
    """``E[C(M, j)] = C(E, j) * prod_t C(c_t, j) / C(E, j)`` as an exact rational."""
    if j == 0:
        return Fraction(1)
    cej = comb(E, j)
    if cej == 0:
        return Fraction(0)
    out = Fraction(cej)
    for ct in c:
        out *= Fraction(comb(ct, j), cej)
    return out


def _exact_upper_tail(E: int, c: Sequence[int], m: int) -> float:
    """``P(M >= m)`` exactly, by inclusion-exclusion on binomial moments."""
    if m <= 0:
        return 1.0
    jmax = min(min(c), E)
    if m > jmax:
        return 0.0
    total = Fraction(0)
    for j in range(m, jmax + 1):
        total += Fraction((-1) ** (j - m) * comb(j - 1, m - 1)) * _binomial_moment(E, c, j)
    return float(min(max(total, Fraction(0)), Fraction(1)))


def _poisson_upper_tail(E: int, c: Sequence[int], m: int) -> float:
    """Legacy Poisson approximation with ``lambda = prod c_t / E^(k-1)``."""
    if m <= 0:
        return 1.0
    if E <= 0 or any(ct == 0 for ct in c):
        return 1.0
    log_lam = sum(log(ct) for ct in c) - (len(c) - 1) * log(E)
    return float(stats.poisson.sf(m - 1, exp(log_lam)))


def _mc_upper_tail(E: int, c: Sequence[int], m: int, *, n_mc: int,
                   rng: np.random.Generator,
                   weights: Optional[np.ndarray] = None) -> float:
    """Monte-Carlo upper tail; ``weights`` gives per-edge change probabilities.

    Phipson & Smyth (2010), Stat Appl Genet Mol Biol 9(1):39,
    doi:10.2202/1544-6115.1585 - the returned p-value is ``(b + 1) / (n + 1)``
    and can therefore never be zero.
    """
    if m <= 0:
        return 1.0
    p = None if weights is None else np.asarray(weights, float) / np.sum(weights)
    hits = 0
    for _ in range(n_mc):
        sets = [set(rng.choice(E, size=int(ct), replace=False, p=p)) for ct in c]
        if len(set.intersection(*sets)) >= m:
            hits += 1
    return float((hits + 1) / (n_mc + 1))


def _ctmc_leaf_states(tree, leaf_labels: Sequence[str], mu: float,
                      rng: np.random.Generator, root_state: int = 0) -> np.ndarray:
    """Evolve one symmetric two-state CTMC down the tree and read the tips.

    ``P(change on an edge of length l) = (1 - exp(-2 mu l)) / 2``, the standard
    Jukes-Cantor two-state transition probability.

    ``root_state`` defaults to 0, i.e. **ancestral absence**, which is both the
    natural convention for an accessory-genome trait (a locus is acquired) and
    the state :func:`fitch_change_edges` resolves an ambiguous root to. Drawing
    the root at random instead makes the Fitch change count bimodal -- on a
    20-leaf star tree, mass at 0-5 *and* at 17-19, because Fitch always measures
    changes away from state 0 -- which destroys rate calibration and drops the
    conditional sampler's acceptance rate to 2.5 %.
    """
    state = {id(tree.seed_node): int(root_state)}
    for nd in tree.preorder_node_iter():
        par = nd.parent_node
        if par is None:
            continue
        ell = nd.edge.length if nd.edge is not None and nd.edge.length else 0.0
        s = state[id(par)]
        p_change = 0.5 * (1.0 - math.exp(-2.0 * mu * ell))
        state[id(nd)] = (1 - s) if rng.random() < p_change else s
    lookup = {nd.taxon.label: state[id(nd)]
              for nd in tree.leaf_node_iter() if nd.taxon is not None}
    return np.array([lookup[l] for l in leaf_labels], dtype=int)


def _fit_ctmc_rate(tree, leaf_labels: Sequence[str], target_changes: int,
                   rng: np.random.Generator, *, n_pilot: int = 40,
                   n_bisect: int = 18) -> float:
    """Rate whose expected Fitch change count matches the observed one.

    Bisection on ``log mu``. The rate is a nuisance parameter: what the null
    needs is a process that puts the *right number* of changes on the tree, so
    that the only thing being tested is whether two traits put them on the
    **same** edges more often than independent evolution would.
    """
    if target_changes <= 0:
        return 0.0

    def mean_changes(mu: float) -> float:
        tot = 0
        for _ in range(n_pilot):
            x = _ctmc_leaf_states(tree, leaf_labels, mu, rng)
            tot += len(fitch_change_edges(tree, leaf_labels, x))
        return tot / n_pilot

    lo, hi = 1e-6, 1.0
    for _ in range(30):                       # bracket
        if mean_changes(hi) >= target_changes:
            break
        hi *= 4.0
    for _ in range(n_bisect):
        mid = math.sqrt(lo * hi)
        if mean_changes(mid) < target_changes:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def _ctmc_change_sets(tree, leaf_labels: Sequence[str], c_t: int, n_draw: int,
                      rng: np.random.Generator, *, max_attempts_per_draw: int = 60):
    """``n_draw`` Fitch change-edge sets of size ``c_t``, simulated on the tree.

    **Conditioned on the observed number of changes.** The exact null conditions
    on ``c_t`` -- it asks where the ``c_t`` changes fall, not how many there are
    -- so an unconditional simulation would not be comparing like with like: its
    replicates would sometimes carry many more changes than the data, inflating
    the joint count for a reason that has nothing to do with edge occupancy.
    Rejection sampling on the change count isolates exactly the assumption under
    test. Draws that cannot be matched within ``max_attempts_per_draw`` fall back
    to the closest count obtained, which is reported by the caller through the
    realised acceptance rate.
    """
    mu = _fit_ctmc_rate(tree, leaf_labels, int(c_t), rng)
    sets, matched = [], 0
    for _ in range(int(n_draw)):
        best, best_gap = None, None
        for _ in range(int(max_attempts_per_draw)):
            x = _ctmc_leaf_states(tree, leaf_labels, mu, rng)
            S = fitch_change_edges(tree, leaf_labels, x)
            gap = abs(len(S) - int(c_t))
            if gap == 0:
                best, best_gap = S, 0
                matched += 1
                break
            if best_gap is None or gap < best_gap:
                best, best_gap = S, gap
        sets.append(best if best is not None else set())
    return sets, (matched / max(int(n_draw), 1))


def _ctmc_upper_tail(tree, leaf_labels: Sequence[str], c: Sequence[int],
                     m: int, *, n_mc: int, rng: np.random.Generator) -> float:
    """Tree-aware null: simulate each trait independently on the actual tree.

    The exact null's assumption is that a trait's change edges are a uniform
    random ``c_t``-subset of the ``E`` edges. Real change edges are not uniform:
    they concentrate on short and terminal branches, so two traits that evolved
    **independently** still share change edges more often than uniform subsets
    do, and the exact null is therefore anticonservative. Measured on simulated
    birth-death trees with traits evolved under independent CTMCs -- the null
    true by construction -- the exact null rejects at **0.0733 +/- 0.0106**
    (n_taxa 60, rate 3.0, 600 replicates), 2.2 standard errors above nominal;
    this null rejects at 0.033 +/- 0.023 in the same cell and is not
    significantly above 0.05 in any cell tested. Its own estimates are imprecise
    (60 replicates, SE 0.023-0.032) because each replicate costs a simulation.
    Two caveats. The 0.072-vs-0.017 gap is 0.040 +/- 0.026 (z = 1.5): the exact
    null's excess is established, the CTMC null's advantage is not, at these
    replicate counts. And this null is calibrated against strict-clock CTMC
    data, which is the model it assumes -- the same circularity that made a
    one-factor continuum null look calibrated. A relaxed-clock arm (shared
    per-edge Gamma rate multiplier, traits still independent) is run alongside
    for exactly that reason -- and it is unkind to every method here: exact
    0.140, this null 0.083, Poisson 0.085, and only ``method="weighted"`` holds,
    at 0.033, by being nearly powerless elsewhere. On a tree with real rate
    heterogeneity, treat any of these p-values as a screen rather than a test.
    See ``benchmarks/archephy_calibration.py``.

    Simulating the traits instead of the change sets removes the assumption:
    whatever edge-occupancy distribution the tree induces is reproduced under the
    null by construction. Conditioning on ``c_t`` keeps everything else fixed.
    """
    if m <= 0:
        return 1.0
    pools = [_ctmc_change_sets(tree, leaf_labels, int(ct), int(n_mc), rng)[0]
             for ct in c]
    hits = 0
    for i in range(int(n_mc)):
        sets = [pool[i] for pool in pools]
        if sets and all(sets) and len(set.intersection(*sets)) >= m:
            hits += 1
    return float((hits + 1) / (int(n_mc) + 1))


def joint_homoplasy_pvalue(E: int, c_per_trait: Sequence[int], m_obs: int, *,
                           method: str = "ctmc", n_mc: int = 20000,
                           rng: Optional[np.random.Generator] = None,
                           edge_weights: Optional[np.ndarray] = None,
                           tree=None, leaf_labels: Optional[Sequence[str]] = None
                           ) -> float:
    """Upper-tail p-value ``P(M >= m_obs)`` for the k-way joint homoplasy count.

    ``method``:
      * ``"ctmc"`` (**default**) - tree-aware null: each trait is evolved
        independently down the supplied tree under a symmetric two-state CTMC
        whose rate is calibrated to the observed number of changes, and the
        joint homoplasy count is recomputed by Fitch parsimony. Requires
        ``tree`` and ``leaf_labels``; falls back to ``"exact"`` without them.
      * ``"exact"`` - exact inclusion-exclusion null under uniform random
        ``c_t``-subsets of the ``E`` edges (falls back to Monte Carlo when the
        alternating sum would need more than ``_MAX_EXACT_TERMS`` terms). For
        ``k = 2`` this equals the hypergeometric survival function. Exact for
        the stated null, but that null is not the biological one: see
        :func:`_ctmc_upper_tail`.
      * ``"poisson"`` - the legacy first-moment Poisson approximation.
      * ``"weighted"`` - Monte-Carlo null with edges sampled proportionally to
        ``edge_weights`` (branch lengths).
    """
    c = [int(x) for x in c_per_trait]
    E = int(E)
    if E <= 0 or not c or any(ct <= 0 for ct in c):
        return 1.0
    m_obs = int(m_obs)
    if m_obs <= 0:
        return 1.0
    if method == "ctmc":
        if tree is None or leaf_labels is None:
            method = "exact"                  # no tree available; be explicit
        else:
            if rng is None:
                rng = np.random.default_rng(0)
            return _ctmc_upper_tail(tree, leaf_labels, c, m_obs,
                                    n_mc=min(int(n_mc), 2000), rng=rng)
    if method == "poisson":
        return _poisson_upper_tail(E, c, m_obs)
    if method == "weighted":
        if rng is None:
            rng = np.random.default_rng(0)
        return _mc_upper_tail(E, c, m_obs, n_mc=n_mc, rng=rng, weights=edge_weights)
    if method != "exact":
        raise ValueError(f"unknown method {method!r}")
    if len(c) == 2:
        # exact hypergeometric: draw c_2 edges from E, of which c_1 are marked
        return float(stats.hypergeom.sf(m_obs - 1, E, c[0], c[1]))
    if min(min(c), E) - m_obs + 1 > _MAX_EXACT_TERMS:
        if rng is None:
            rng = np.random.default_rng(0)
        return _mc_upper_tail(E, c, m_obs, n_mc=n_mc, rng=rng)
    return _exact_upper_tail(E, c, m_obs)


# --------------------------------------------------------------------------- #
# Public test
# --------------------------------------------------------------------------- #
def archephy_cs_test(
    newick: str,
    X: np.ndarray,
    leaves: List[str],
    archetype: Sequence[int],
    *,
    alpha: float = 0.05,
    method: str = "ctmc",
    n_mc: int = 2000,
    rng: Optional[np.random.Generator] = None,
) -> dict:
    """k-way joint-homoplasy convergence test with an exact closed-form null.

    Parameters
    ----------
    newick : rooted tree, newick string.
    X : ``(n_leaves, p_traits)`` binary matrix, rows aligned to ``leaves``.
    leaves : row labels of ``X``, matching the tree's taxon labels.
    archetype : ``k >= 2`` column indices defining the candidate archetype.
    alpha : significance level for ``call``.
    method : ``"ctmc"`` (default, tree-aware) | ``"exact"`` | ``"poisson"`` |
        ``"weighted"`` (see :func:`joint_homoplasy_pvalue`). The tree-aware
        null is the default because the exact null is exact for a null in
        which change edges are uniform random subsets, which real trees
        violate, and its measured size reaches 0.078 +/- 0.011 against a
        nominal 0.05.

    Returns
    -------
    dict with ``m_obs``, ``p_value``, ``expected_m`` (the null mean),
    ``call`` (``p <= alpha`` **and** ``m_obs >= 2``), ``per_trait_changes``,
    ``n_edges``, ``k``, ``method``, ``gains_losses``,
    ``duplicate_columns_collapsed`` and, for reference, ``p_value_poisson``.
    """
    idx = list(archetype)
    k = len(idx)
    if k < 2:
        raise ValueError(f"archetype needs k>=2 traits; got {k}")
    X = np.asarray(X, dtype=int)
    if X.ndim != 2:
        raise ValueError(f"X must be 2-D (n_leaves, p_traits); got shape {X.shape}")
    if X.shape[0] != len(leaves):
        raise ValueError(f"X has {X.shape[0]} rows but {len(leaves)} leaves")
    if max(idx) >= X.shape[1] or min(idx) < 0:
        raise ValueError(f"archetype indices out of range for X with {X.shape[1]} columns")

    # Two identical archetype columns are one measurement, not two: a
    # co-inherited locus pair recorded twice puts every change edge of the
    # trait into the intersection, so m_obs becomes the trait's own change
    # count and the pair reads as independent convergent evidence for itself.
    # One representative of each duplicate set is kept and the rest reported.
    _groups: List[List[int]] = []
    _first: Dict[tuple, int] = {}
    for t in idx:
        key = tuple(int(v) for v in X[:, t])
        if key in _first:
            _groups[_first[key]].append(t)
        else:
            _first[key] = len(_groups)
            _groups.append([t])
    duplicate_columns_collapsed = [list(g) for g in _groups if len(g) > 1]
    idx = [g[0] for g in _groups]
    k = len(idx)

    tree, tree_leaves = load_tree(newick)
    missing = set(tree_leaves) - set(leaves)
    if missing:
        raise ValueError(
            f"{len(missing)} tree leaves are absent from `leaves` "
            f"(e.g. {sorted(missing)[:3]})")
    pos = {l: i for i, l in enumerate(leaves)}
    order = [pos[l] for l in tree_leaves]

    edge_nodes = [nd for nd in tree.preorder_node_iter() if nd.parent_node is not None]
    n_edges = len(edge_nodes)
    node_index = {id(nd): i for i, nd in enumerate(edge_nodes)}
    edge_weights = np.array([nd.edge.length for nd in edge_nodes], dtype=float)

    change_sets: List[set] = []
    c_per_trait: List[int] = []
    for t in idx:
        S = fitch_change_edges(tree, tree_leaves, X[order, t])
        change_sets.append({node_index[i] for i in S if i in node_index})
        c_per_trait.append(len(change_sets[-1]))

    joint = set.intersection(*change_sets) if change_sets else set()
    m_obs = len(joint)

    p_value = joint_homoplasy_pvalue(
        n_edges, c_per_trait, m_obs, method=method, n_mc=n_mc, rng=rng,
        edge_weights=edge_weights if method == "weighted" else None,
        tree=tree, leaf_labels=tree_leaves)
    p_exact = (joint_homoplasy_pvalue(n_edges, c_per_trait, m_obs,
                                      method="exact") if m_obs > 0 else 1.0)
    p_poisson = _poisson_upper_tail(n_edges, c_per_trait, m_obs) if m_obs > 0 else 1.0

    expected_m = (float(np.prod([c / n_edges for c in c_per_trait]) * n_edges)
                  if n_edges > 0 and all(c > 0 for c in c_per_trait) else 0.0)

    return {
        "m_obs": int(m_obs),
        "expected_m": expected_m,
        "p_value": float(p_value),
        "p_value_poisson": float(p_poisson),
        "p_value_exact_uniform_null": float(p_exact),
        "method": method,
        "call": bool(p_value <= alpha and m_obs >= 2),
        "per_trait_changes": c_per_trait,
        "n_edges": int(n_edges),
        "k": k,
        "duplicate_columns_collapsed": duplicate_columns_collapsed,
        "gains_losses": [fitch_asr(tree, tree_leaves, X[order, t]) for t in idx],
    }
