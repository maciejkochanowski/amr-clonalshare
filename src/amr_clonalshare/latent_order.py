"""Bounds on the lineage share of the latent MIC ordering.

What is bounded
---------------
Within a stratum s, normally a testing laboratory, let F_s be the distribution
of the latent MIC and U = F_s(MIC) its probability integral transform,
which is uniform on (0, 1) in every stratum. The lineage share of the latent
MIC ordering is

    rho = Var(E[U | lineage]) / Var(U) = 12 * sum_g w_g (E[U | g] - 1/2)^2,

the rank intraclass correlation of the latent MIC within strata, with the
lineages weighted by their share w_g of the isolates. The lineages and their
counts in every stratum are those observed; F_s is the mixture of the
lineages' distributions in stratum s at those counts. U is pooled over strata, so the
share measures lineage differences that hold across strata: opposite lineage
effects in two strata cancel. ``latent_order_lower_within_strata`` reports
the lower bound of every stratum on its own beside the pooled one.

A dilution panel does not observe U. It observes the reading (a, b], and
within its stratum U then lies in [F_s(a), F_s(b)]. How the isolates of
different lineages are ordered inside one reading is not observed, so rho is
not identified from the readings.

What the readings identify (Theorem 1)
--------------------------------------
The readings identify, for every reading r of a stratum, its range
[L_r, H_r] = [F_s(a), F_s(b)], its share P(r) of all isolates and the share
pi_gr of its isolates that lineage g holds. Inside reading r the lineages'
conditional distributions of U can be anything that mixes back to the uniform distribution on
[L_r, H_r]. By Strassen's theorem the conditional means m_gr they can have are
exactly those whose distribution, with weights pi_gr, is below that uniform distribution in the
convex order, which for y_gr = pi_gr m_gr reads

    sum_{g in A} y_gr <= f_r(A) = pi_r(A) H_r - (H_r - L_r) pi_r(A)^2 / 2

for every set A of lineages, with equality for all of them. f_r is
submodular and the feasible y form its base polytope. The lineage means enter
rho through x_g = w_g E[U | g] = sum_r P(r) y_gr, and a sum of base
polytopes is the base polytope of the sum of the functions, so the x the
readings allow are exactly the base polytope of F = sum_r P(r) f_r, and
rho = 12 (sum_g x_g^2 / w_g - 1/4) ranges over an interval [rho_min, rho_max]
as x ranges over it.

* rho_min, the lower bound the readings establish, minimises a convex function over
  the polytope. It is computed by the minimum-norm-point algorithm of Wolfe
  (1976) in the norm sum x_g^2 / w_g, whose linear step is the greedy vertex
  of Edmonds (1970); the returned gap bounds the distance to the minimum,
  and the reported lower end is certified not to exceed it.
* rho_max is attained at a vertex: one order of the lineages, the same in
  every reading. At the vertex of an order, x_g = u_g - sum_{h above g} c_gh
  with u_g = F({g}) and c_gh = sum_r P(r) (H_r - L_r) pi_gr pi_hr, so the
  largest sum x_g^2 / w_g over orders is a maximum over the sets of lineages
  placed above each lineage. An insertion search returns the largest share
  it attains, which the readings allow, and a branch-and-bound over the
  orders certifies an upper bound on rho_max; where the two meet, the upper
  end is exact. Where they do not and there are at most EXACT_LINEAGES
  lineages, dynamic programming over the sets finds rho_max exactly (with
  at most DP_FIRST lineages it is used directly); with more, the gap between
  the attained share and the certified bound is reported. Every share
  between rho_min and the attained share is allowed by the readings.

The dual certificate (Theorem 2)
--------------------------------
For every vector z of lineage scores,

    rho >= 12 G(z)_+^2 / Q(z),
    G(z) = sum_s sum_{g, h: z_g > z_h} omega_ghs (z_g - z_h) (P(Y_h < Y_g) - 1/2),

with omega_ghs = n_gs n_hs / (n n_s), Q(z) = sum_g w_g (z_g - zbar)^2 and
P(Y_h < Y_g) the probability that a reading of lineage h in stratum s lies
wholly below a reading of lineage g there. It follows from E[U | g, s] - 1/2 =
sum_h (n_hs / n_s) (P(U_h < U_g) - 1/2), from P(U_h < U_g) >= P(Y_h < Y_g)
and from the Cauchy-Schwarz inequality, and it needs only that every reading
holds its MIC: not one panel per stratum, not coarsening at random. When every
stratum is read on one panel, the largest value over z is rho_min, attained at
z = E[U | g] of the minimum-norm point. G(z) is linear in the probabilities
P(Y_h < Y_g), and the readings of different isolates are independent, so for
a fixed z it has an unbiased estimate, a U-statistic over pairs of isolates.

The lower confidence limit (split certificate)
----------------------------------------------
Within every lineage-by-stratum cell the isolates are split at random into
two halves, and the split is used both ways: in each direction one half
chooses z, the minimum-norm point of its own readings, and the other
estimates G(z) without bias. Scaled by sqrt(Q(z)), every estimate has an
expectation T(z) = G(z) / sqrt(Q(z)) <= sqrt(rho_min / 12) <= sqrt(rho / 12),
whatever z the other half chose. The limit averages the scaled estimates over
the two directions of SPLITS random splits: the average has an expectation
below sqrt(rho_min / 12), and its standard deviation is at most the average
of their standard deviations (Minkowski's inequality), whatever the
dependence between the splits and the directions. With Student's t at the
smallest of their degrees of freedom,

    L = 12 max(0, mean T - t_(1 - alpha) mean sd)^2 <= rho_min <= rho

with probability at least 1 - alpha, to the accuracy of the normal
approximation to the U-statistics. Pairs of lineages one of which has no
isolate of a stratum in the estimating half are counted at the worst case,
P = 0. The variance of an estimate is the variance of its first-order
projection, estimated from all isolates of every cell under the cell's distribution
smoothed by Perks' prior of total mass one over the readings of its stratum
(as the bootstrap of the share smooths it), so that a cell whose isolates
share one reading is not given a variance of zero, plus the variance of the
pair terms, bounded above; its degrees of freedom are Satterthwaite's. The
limit is reported at most at the lower end of the observed readings: lowering
a lower confidence limit keeps it one.

Assumptions: each reading holds the MIC it stands for, and the latent MIC
distribution of each stratum is continuous; the lineage share of the
latent ordering is then bounded by the readings. The bounds are sharp
(``latent_order_sharp``) when the distinct readings of every stratum are
disjoint intervals, as they are when a stratum is read on one panel. Where
readings of one stratum overlap, every run of overlapping readings is taken
as one cell: the isolates of the run may then stand in any order inside it,
which allows every arrangement the readings allow and some they do not, so
the lower end stays a lower bound and the upper end an upper bound, and the
lower confidence limit stays valid (its certificate counts only readings that
lie wholly below one another). The limit treats the
readings of different isolates as independent given the represented lineages
(and any sampling units, such as farms), and conditions on the number of
isolates of every lineage in every stratum. Singleton lineages are set
aside, as for the share.
"""
from __future__ import annotations

from typing import Any, Dict

import numpy as np

from .mic_order import _cell_ranges

__all__ = ["Readings", "order_bounds", "call_bounds", "bounds_from_bins"]

#: Stopping rule of the minimum-norm-point algorithm, on the scale of rho:
#: the gap it returns bounds the distance of the lower end to the minimum.
GAP_TOLERANCE = 1e-10
MAX_MAJOR = 20000
#: Major cycles without a new smallest sum x^2 / w after which only rounding
#: is taken to be left.
PATIENCE = 10
#: Lineages up to which the upper end is found exactly by dynamic programming
#: when the branch-and-bound leaves a gap, and up to which the dynamic
#: programme is tried first, because it is then the quicker.
EXACT_LINEAGES = 20
DP_FIRST = 12
#: Nodes the branch-and-bound over orders may open when certifying the upper
#: end with more lineages; a count, not a time, so the result does not depend
#: on the machine.
UPPER_NODES = 500
#: Random splits, each used both ways, averaged into the lower confidence limit.
SPLITS = 25


class _Polytope:
    """Lineage-by-reading shares and reading ranges of one data set."""

    def __init__(self, counts: np.ndarray, low: np.ndarray, high: np.ndarray):
        counts = np.asarray(counts, dtype=float)
        self.counts = counts
        self.p = counts / counts.sum()
        self.Pr = self.p.sum(0)
        self.w = self.p.sum(1)
        self.L = np.asarray(low, dtype=float)
        self.H = np.asarray(high, dtype=float)
        self.D = self.H - self.L
        with np.errstate(invalid="ignore", divide="ignore"):
            self.pi = np.where(self.Pr > 0, self.p / np.where(self.Pr > 0, self.Pr, 1.), 0.)

    def vertex(self, order) -> np.ndarray:
        """x at the vertex for one order of the lineages, top slice first."""
        order = np.asarray(order, dtype=int)
        P = self.pi[order]
        T = np.cumsum(P, axis=0) - P
        m = self.H - self.D * (T + P / 2.)
        x = np.empty(order.size)
        x[order] = np.sum(self.p[order] * m, axis=1)
        return x

    def rho(self, x) -> float:
        return float(12. * (np.sum(np.asarray(x) ** 2 / self.w) - 0.25))

    def midpoint(self) -> np.ndarray:
        return self.p @ ((self.L + self.H) / 2.)

    def resolution(self) -> float:
        return float(1. - np.sum(self.Pr * self.D ** 2))

    def pieces(self):
        """u_g = F({g}), the coordinate of lineage g at the top of the order;
        c_gh = sum_r P(r) D_r pi_gr pi_hr, by which h above g lowers it; and
        l_g = F(V) - F(V - g), its coordinate at the bottom."""
        u = self.p @ self.H - np.sum(self.p * self.D * self.pi, axis=1) / 2.
        C = (self.pi * (self.Pr * self.D)) @ self.pi.T
        low = u - (C.sum(1) - np.diag(C))
        return u, C, low

    def merged(self):
        """The same readings with the lineages of one profile merged, and the
        group of every lineage.

        Lineages whose isolates spread over the readings in the same
        proportions can pool their distributions inside every reading without leaving
        the polytope, which moves their means to the weighted average and
        cannot raise sum x^2 / w. At the minimum they share one mean, so
        merging them leaves the minimum unchanged, and it removes the ties
        that otherwise hold the minimum-norm point on a face of many alike
        lineages, as with calls, where lineages of one rate are alike."""
        profile = self.counts / self.counts.sum(1, keepdims=True)
        _, group = np.unique(profile, axis=0, return_inverse=True)
        group = group.ravel()
        counts = np.zeros((int(group.max()) + 1, self.counts.shape[1]))
        np.add.at(counts, group, self.counts)
        return _Polytope(counts, self.L, self.H), group


def _minimum_norm(poly: _Polytope):
    """Minimum of sum x^2 / w over the base polytope (Wolfe's algorithm).

    Returns a certified lower end, the gap 12 * 2 <x, x - q>, which bounds
    how far rho at the point x can be above its minimum (q is the greedy
    vertex at x), and E[U | g] of every lineage at the point. By convexity
    rho(x) - gap is below the minimum at any x, so every point visited
    certifies a lower end; the best one is returned: rho at its point less
    the larger of its gap and GAP_TOLERANCE, and never below zero. It cannot
    exceed the minimum through the stopping rule or through rounding, so
    that a share of zero is reported as zero and not as a positive value in
    the last digits.

    Lineages of one profile are merged first, which leaves the minimum
    unchanged (_Polytope.merged). The affine step solves its least-squares
    problem on the corral's differences rather than on their Gram matrix.
    The algorithm stops when the gap falls to GAP_TOLERANCE, when PATIENCE
    major cycles in a row have not lowered sum x^2 / w below its smallest
    value, which in exact arithmetic every cycle does, so that only rounding
    is left, or after MAX_MAJOR cycles; the gap then says how sharp the lower
    end is."""
    poly, group = poly.merged()
    w = poly.w
    root = np.sqrt(w)

    def ip(a, b):
        return float(np.sum(a * b / w))

    x = poly.vertex(np.argsort(poly.midpoint() / w, kind="stable"))
    corral = [x.copy()]
    lam = np.array([1.0])
    best, best_gap, best_x, smallest, since = -np.inf, np.inf, x.copy(), np.inf, 0
    for _ in range(MAX_MAJOR):
        q = poly.vertex(np.argsort(x / w, kind="stable"))
        norm = ip(x, x)
        gap = max(24. * (norm - ip(x, q)), 0.)
        if poly.rho(x) - max(gap, GAP_TOLERANCE) > best:
            best, best_gap, best_x = poly.rho(x) - max(gap, GAP_TOLERANCE), gap, x.copy()
        smallest, since = (norm, 0) if norm < smallest else (smallest, since + 1)
        if gap <= GAP_TOLERANCE or since >= PATIENCE:
            break
        corral.append(q)
        lam = np.append(lam, 0.)
        while True:
            A = np.asarray(corral)
            scaled = A / root
            mu = np.linalg.lstsq((scaled[1:] - scaled[0]).T, -scaled[0], rcond=None)[0]
            a = np.r_[1. - mu.sum(), mu]
            if np.all(a > 1e-14):
                lam, x = a, a @ A
                break
            neg = a <= 1e-14
            theta = float(np.min(lam[neg] / (lam[neg] - a[neg])))
            lam = lam + theta * (a - lam)
            keep = lam > 1e-14
            corral = [s for s, kept in zip(corral, keep) if kept]
            lam = lam[keep] / lam[keep].sum()
            x = lam @ np.asarray(corral)
    return max(0., best), best_gap, (best_x / w)[group]


def _lower(poly: _Polytope):
    """The certified lower end and its gap (see :func:`_minimum_norm`)."""
    lower, gap, _ = _minimum_norm(poly)
    return lower, gap


def _upper_search(poly: _Polytope, starts: int, rng: np.random.Generator):
    """Largest sum x^2 / w found over vertices by insertion search from
    several orders, and the order that attains it.

    Every lineage in turn is taken out of the order and put back where the
    sum is largest. Taking lineage g out raises the coordinate of every
    lineage below it by c_l = sum_r p_lr D_r pi_gr; putting it back at
    position j lowers those below j by the same c_l, so with prefix sums over
    the order the value at every position j comes at once, in O(G B) for one
    lineage. A pass that improves the sum is followed by another."""
    G = poly.w.size
    w = poly.w
    pD = poly.p * poly.D
    top = poly.p @ poly.H - np.sum(pD * poly.pi, axis=1) / 2.     # x_g at the top of the order
    e = poly.midpoint() / w
    orders = [np.argsort(-e, kind="stable")] + [rng.permutation(G) for _ in range(starts)]

    def coordinates(order):
        P = poly.pi[order]
        return top[order] - np.sum(pD[order] * (np.cumsum(P, axis=0) - P), axis=1)

    best, best_order = -np.inf, orders[0]
    for order in orders:
        order = np.asarray(order, dtype=np.int64)
        x = coordinates(order)
        value = float(np.sum(x ** 2 / w[order]))
        improved = True
        while improved:
            improved = False
            for i in range(G):
                g = order[i]
                rest = np.delete(order, i)
                c = pD[rest] @ poly.pi[g]
                x0 = np.delete(x, i)
                x0[i:] += c[i:]
                wr = w[rest]
                kept_above = np.r_[0., np.cumsum(x0 ** 2 / wr)]
                pushed_down = np.r_[0., np.cumsum((x0 - c) ** 2 / wr)]
                above = np.vstack([np.zeros(poly.pi.shape[1]), np.cumsum(poly.pi[rest], axis=0)])
                xg = top[g] - above @ pD[g]
                values = kept_above + (pushed_down[-1] - pushed_down) + xg ** 2 / w[g]
                j = int(np.argmax(values))
                if values[j] > value + 1e-12 / 12.:
                    value, improved = float(values[j]), True
                    order = np.insert(rest, j, g)
                    x = coordinates(order)
        if value > best:
            best, best_order = value, order
    return best, best_order


def _upper_exact(poly: _Polytope) -> float:
    """Largest sum x^2 / w over the vertices, by dynamic programming.

    At the vertex of an order the coordinate of lineage g depends only on
    the set S of lineages above it, x_g = u_g - sum_{h in S} c_gh, so the
    largest sum over the orders of a set T placed at the top is
    V(T) = max_{g in T} V(T - g) + (u_g - sum_{h in T - g} c_gh)^2 / w_g,
    g the lowest of T. V over all lineages is the maximum; the sets are
    visited by size, O(2^G G^2) in all."""
    from itertools import combinations
    u, C, _ = poly.pieces()
    w = poly.w
    G = w.size
    V = np.full(1 << G, -np.inf)
    V[0] = 0.
    weights = 1 << np.arange(G)
    for k in range(1, G + 1):
        members = np.array(list(combinations(range(G), k)), dtype=np.int64)
        masks = weights[members].sum(1)
        inside = np.zeros((masks.size, G))
        np.put_along_axis(inside, members, 1., axis=1)
        pulled = inside @ C
        best = np.full(masks.size, -np.inf)
        for j in range(k):
            g = members[:, j]
            x = u[g] - (pulled[np.arange(masks.size), g] - C[g, g])
            best = np.maximum(best, V[masks ^ weights[g]] + x * x / w[g])
        V[masks] = best
    return float(V[-1])


def _upper_certificate(poly: _Polytope, incumbent: float):
    """An upper bound on the largest sum x^2 / w, by branch and bound over
    orders built from the top, and the largest sum it attains.

    A node is a set S placed at the top with the best sum found for it. The
    lineages below it lie in the base polytope of F(S + T) - F(S), where
    lineage g ranges from l_g (at the bottom) to u_g - sum_{h in S} c_gh
    (directly below S). Two bounds on their sum hold there: the chord bound,
    x^2 <= (l + u) x - l u on [l, u], whose linear right side is largest at
    the greedy vertex, where it is sum_g a_g (u_g - sum_S c_gh) less the sum
    over pairs of min(a_g, a_h) c_gh, a_g = (l_g + u_g) / w_g; and the box
    bound, max(l^2, u^2) per lineage. A node enters the queue under the bound
    of its parent, which holds for it too, and its own bound is computed when
    it is taken out: a node whose bound cannot beat the incumbent is closed,
    one whose bound fell is put back under it, and one whose bound stands is
    opened. A set reached again with a smaller sum is closed. Lineages with
    the same counts in every reading are interchangeable, so of those not yet
    placed only the first is placed next. At most UPPER_NODES bounds are
    computed; the certified bound is the largest of the incumbent and the
    bounds left in the queue, each of which holds for every order below its
    node."""
    import heapq
    u, C, low = poly.pieces()
    w = poly.w
    G = w.size
    _, kind = np.unique(poly.counts, axis=0, return_inverse=True)
    kind = kind.ravel()

    def bound(pulled, rest):
        if rest.size == 0:
            return 0.
        up = u[rest] - pulled[rest]
        lo = low[rest]
        a = (lo + up) / w[rest]
        pairs = np.minimum.outer(a, a) * C[np.ix_(rest, rest)]
        linear = float(np.sum(a * up) - (pairs.sum() - np.trace(pairs)) / 2.)
        chord = linear - float(np.sum(lo * up / w[rest]))
        box = float(np.sum(np.maximum(lo * lo, up * up) / w[rest]))
        return min(chord, box)

    everyone = np.arange(G)
    best = incumbent
    seen: Dict[int, float] = {}
    # (-ceiling, whether the ceiling is the node's own, sum, set, tie-break)
    heap = [(-np.inf, False, 0., 0, 0)]
    computed = 0
    ticket = 1
    left = -np.inf
    while heap:
        negative, own, value, mask, _ = heap[0]
        if -negative <= best:
            break
        if not own and computed >= UPPER_NODES:
            left = -negative
            break
        heapq.heappop(heap)
        inside = np.array([(mask >> g) & 1 for g in range(G)], dtype=bool)
        rest = everyone[~inside]
        pulled = C[:, inside].sum(1)
        if not own:
            computed += 1
            ceiling = value + bound(pulled, rest)
            if ceiling > best:
                heapq.heappush(heap, (-min(ceiling, -negative), True, value, mask, ticket))
                ticket += 1
            continue
        placed_kinds = set()
        for g in rest:
            if kind[g] in placed_kinds:
                continue
            placed_kinds.add(kind[g])
            x = u[g] - pulled[g]
            total = value + x * x / w[g]
            child = mask | (1 << int(g))
            if seen.get(child, -np.inf) >= total:
                continue
            seen[child] = total
            if rest.size == 1:
                best = max(best, total)
                continue
            heapq.heappush(heap, (negative, False, total, child, ticket))
            ticket += 1
    return max(best, left), best


def _upper(poly: _Polytope, starts: int, rng: np.random.Generator):
    """(attained, bound, exact) on the scale of rho: the largest share at a
    vertex found, an upper bound on the largest share the readings allow,
    and whether the two agree."""
    to_rho = lambda s: float(12. * (s - 0.25))
    if poly.w.size <= DP_FIRST:
        value = _upper_exact(poly)
        return to_rho(value), to_rho(value), True
    found, _ = _upper_search(poly, starts, rng)
    bound, attained = _upper_certificate(poly, found)
    attained = max(found, attained)
    if bound - attained > 1e-12 and poly.w.size <= EXACT_LINEAGES:
        value = _upper_exact(poly)
        return to_rho(value), to_rho(value), True
    exact = bound - attained <= 1e-12
    return to_rho(attained), to_rho(max(bound, attained)), bool(exact)


class Readings:
    """Readings of the repeated lineages, grouped by stratum and reading.

    ``bins`` indexes the distinct readings of every stratum and ``cell`` the
    disjoint cell every reading falls in: the reading itself where the
    readings of the stratum are disjoint, otherwise the run of overlapping
    readings it belongs to. The cells of a stratum are ordered, so a
    subsample recomputes their ranges from its own counts. ``below`` holds,
    per stratum, which reading lies wholly below which."""

    def __init__(self, lo, hi, code, strata):
        self.code = np.asarray(code, dtype=int)
        self.G = int(self.code.max()) + 1 if self.code.size else 0
        strata = np.asarray(strata, dtype=object)
        self.stratum = np.empty(self.code.size, dtype=int)
        self.bins = np.empty(self.code.size, dtype=int)
        self.groups = []            # per stratum: (bin ids, first/last cells of each reading, disjoint)
        self.cells = []             # per stratum: the ordered cell ids of the stratum
        self.below = []             # per stratum: below[a, b] = reading a wholly below reading b
        self.names = sorted(set(strata.tolist()), key=str)
        start = 0
        cell_start = 0
        cell_of_bin = []
        for s, name in enumerate(self.names):
            idx = np.flatnonzero(strata == name)
            self.stratum[idx] = s
            first, last, _ = _cell_ranges(lo[idx], hi[idx])
            pairs, inverse = np.unique(np.column_stack([first, last]), axis=0, return_inverse=True)
            inverse = inverse.ravel()
            order = np.lexsort((pairs[:, 1], pairs[:, 0]))
            rank = np.empty(order.size, dtype=int)
            rank[order] = np.arange(order.size)
            pairs = pairs[order]
            self.bins[idx] = start + rank[inverse]
            # readings sorted by their first cell: a reading opens a new cell
            # when it starts above everything before it, otherwise it joins
            # the run of the readings it overlaps
            run = np.zeros(len(pairs), dtype=int)
            reach = pairs[0, 1]
            for b in range(1, len(pairs)):
                if pairs[b, 0] > reach:
                    run[b] = run[b - 1] + 1
                else:
                    run[b] = run[b - 1]
                reach = max(reach, pairs[b, 1])
            disjoint = bool(run[-1] == len(pairs) - 1)
            cell_of_bin.extend(cell_start + run)
            self.groups.append((np.arange(start, start + len(pairs)), pairs, disjoint))
            self.cells.append(np.arange(cell_start, cell_start + run[-1] + 1))
            self.below.append((pairs[:, 1][:, None] < pairs[:, 0][None, :]).astype(float))
            start += len(pairs)
            cell_start += run[-1] + 1
        self.n_bins = start
        self.cell = np.asarray(cell_of_bin, dtype=int)
        self.n_cells = cell_start

    @property
    def sharp(self) -> bool:
        return all(group[2] for group in self.groups)

    def ranges(self, cell_counts: np.ndarray):
        """[L, H] of every cell from the counts of one (sub)sample: the cells
        of a stratum are disjoint and ordered, so each takes its share of the
        stratum's isolates, in order."""
        low = np.zeros(self.n_cells)
        high = np.zeros(self.n_cells)
        for ids in self.cells:
            c = cell_counts[ids].astype(float)
            if not (c > 0).any():
                continue
            cum = np.r_[0., np.cumsum(c / c.sum())]
            low[ids], high[ids] = cum[:-1], cum[1:]
        return low, high

    def polytope(self, rows=None, lineages=None) -> _Polytope:
        """The polytope of the rows, over the cells; ``lineages`` keeps those
        rows (lineages without a row would have no weight)."""
        code = self.code if rows is None else self.code[rows]
        cells = self.cell[self.bins] if rows is None else self.cell[self.bins[rows]]
        G = self.G if lineages is None else lineages.size
        if lineages is not None:
            index = -np.ones(self.G, dtype=int)
            index[lineages] = np.arange(lineages.size)
            code = index[code]
        counts = np.zeros((G, self.n_cells))
        np.add.at(counts, (code, cells), 1.)
        low, high = self.ranges(counts.sum(0))
        keep = counts.sum(0) > 0
        return _Polytope(counts[:, keep], low[keep], high[keep])

    def cell_counts(self, rows):
        """Per stratum, the counts of the rows by lineage and reading."""
        out = []
        for s, (ids, _, _) in enumerate(self.groups):
            counts = np.zeros((self.G, ids.size))
            sel = rows[self.stratum[rows] == s]
            np.add.at(counts, (self.code[sel], self.bins[sel] - ids[0]), 1.)
            out.append(counts)
        return out


def _certificate(readings: Readings, z: np.ndarray, part, full, n_cell):
    """Estimate of G(z) from the counts ``part`` of the evaluation half,
    with the variance of the estimate and its degrees of freedom.

    ``full`` holds the counts of all isolates, from which the variances are
    estimated, and ``n_cell`` the fixed counts per lineage and stratum that
    weight the pairs. See the module docstring."""
    n = n_cell.sum()
    estimate = 0.
    terms, dfs = [], []
    for s, (B, F, D) in enumerate(zip(part, full, readings.below)):
        ns = n_cell[:, s].sum()
        if ns == 0:
            continue
        mB, mF = B.sum(1), F.sum(1)
        present = n_cell[:, s] > 0
        omega = np.outer(n_cell[:, s], n_cell[:, s]) / (n * ns)
        dz = z[:, None] - z[None, :]
        A = np.where((dz > 0) & present[:, None] & present[None, :], omega * dz, 0.)   # g above h
        seen = mB > 0
        MB = np.where(seen[:, None], B / np.maximum(mB, 1)[:, None], 0.)
        both = seen[:, None] & seen[None, :]
        below = MB @ D @ MB.T                   # [h, g]: a reading of h below one of g
        estimate += float(np.sum(np.where(both, A * (below.T - .5), -.5 * A)))
        # the first-order projection of every reading of every cell, from
        # the readings of all isolates
        MF = np.where((mF > 0)[:, None], F / np.maximum(mF, 1)[:, None], 0.)
        Ab = np.where(both, A, 0.)
        phi = Ab @ (MF @ D) + Ab.T @ (MF @ D.T)
        K = F.shape[1]
        for g in np.flatnonzero(seen):
            counts, f = F[g], phi[g]
            size = counts.sum()
            # the variance of the projection under the cell's distribution smoothed
            # as the bootstrap smooths it, Perks' prior of total mass one
            # over the readings of the stratum, so that a cell whose isolates
            # share one reading is not given a variance of zero
            law = (counts + 1. / K) / (size + 1.)
            spread = float(law @ (f - law @ f) ** 2) * size / max(size - 1., 1.)
            terms.append(spread / mB[g])
            dfs.append(size)
        # the pair terms: Var(k) <= theta (1 - theta) for a kernel in {0, 1}
        whole = MF @ D @ MF.T
        pairs = np.outer(mF, mF)
        theta = (whole.T * pairs + .5) / (pairs + 1.)
        second = float(np.sum(np.where(both, Ab ** 2 * theta * (1. - theta)
                                       / np.maximum(np.outer(mB, mB), 1.), 0.)))
        if second > 0:
            terms.append(second)
            dfs.append(np.inf)
    terms_a = np.asarray(terms)
    variance = float(terms_a.sum())
    spread_df = float(np.sum(terms_a ** 2 / np.asarray(dfs))) if terms else 0.
    df = variance * variance / spread_df if spread_df > 0 else np.inf
    return estimate, variance, df


def _split(readings: Readings, rng: np.random.Generator):
    """Rows of the two halves, each cell (lineage, stratum) split at random,
    an odd cell's extra isolate falling to either half alike."""
    key = readings.code * (len(readings.names) + 1) + readings.stratum
    order = np.argsort(key, kind="stable")
    edges = np.flatnonzero(np.diff(np.r_[-1, key[order], -1]))
    first, second = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        rows = rng.permutation(order[a:b])
        k = rows.size / 2.
        k = int(np.floor(k)) + int(rng.random() < k - np.floor(k))
        first.append(rows[:k])
        second.append(rows[k:])
    return np.concatenate(first), np.concatenate(second)


def _direction(readings: Readings, rows: np.ndarray, w: np.ndarray):
    """z from the minimum-norm point of the rows (lineages without a row take
    the weighted mean), centred, with Q(z); None if it has no spread."""
    chosen = np.unique(readings.code[rows])
    if chosen.size < 2:
        return None
    _, _, means = _minimum_norm(readings.polytope(rows, lineages=chosen))
    z = np.zeros(readings.G)
    z[chosen] = means - .5
    others = np.setdiff1d(np.arange(readings.G), chosen)
    z[others] = float(np.sum(w[chosen] * z[chosen]) / w[chosen].sum())
    z = z - float(np.sum(w * z))
    Q = float(np.sum(w * z * z))
    return (z, Q) if Q > 1e-15 else None


def _certified_terms(readings: Readings, rng: np.random.Generator, full, n_cell):
    """One random split used both ways: in each direction one half chooses
    z and the other estimates G(z); returns the scaled estimates
    G / sqrt(Q), their standard deviations and degrees of freedom."""
    w = n_cell.sum(1) / n_cell.sum()
    halves = _split(readings, rng)
    out = []
    for choose, evaluate in (halves, halves[::-1]):
        found = _direction(readings, choose, w)
        if found is None:
            out.append((0., 0., np.inf))
            continue
        z, Q = found
        estimate, variance, df = _certificate(readings, z, readings.cell_counts(evaluate), full, n_cell)
        out.append((estimate / np.sqrt(Q), np.sqrt(variance / Q), df))
    return out


def _lower_limit(readings: Readings, alpha: float, rng: np.random.Generator) -> float:
    """The one-sided (1 - alpha) lower confidence limit of rho_min and rho
    from SPLITS random splits, each used both ways (module docstring)."""
    from scipy.stats import norm, t as student
    n_cell = np.zeros((readings.G, len(readings.names)))
    np.add.at(n_cell, (readings.code, readings.stratum), 1.)
    full = readings.cell_counts(np.arange(readings.code.size))
    terms = [term for _ in range(SPLITS) for term in _certified_terms(readings, rng, full, n_cell)]
    estimate = float(np.mean([t for t, _, _ in terms]))
    sd = float(np.mean([s for _, s, _ in terms]))
    df = float(min(d for _, _, d in terms))
    q = float(student.ppf(1. - alpha, df)) if np.isfinite(df) else float(norm.ppf(1. - alpha))
    return 12. * max(estimate - q * sd, 0.) ** 2


def bounds_from_bins(readings: Readings, *, alpha: float = 0.05, starts: int = 10,
                     seed=0) -> Dict[str, Any]:
    """The interval of rho the readings allow, with a lower confidence limit."""
    rng = np.random.default_rng(seed)
    poly = readings.polytope()
    lower, gap = _lower(poly)
    attained, bound, exact = _upper(poly, starts, rng)
    lower = float(max(lower, 0.))
    # a lower confidence limit of rho_min stays one when it is lowered, so it
    # is never reported above the lower end of the observed readings
    limit = min(_lower_limit(readings, alpha, np.random.default_rng(rng.integers(2 ** 63))), lower)
    record: Dict[str, Any] = dict(
        latent_order_lower=lower,
        latent_order_upper=float(min(attained, 1.)),
        latent_order_upper_bound=float(min(bound, 1.)),
        latent_order_upper_exact=bool(exact),
        latent_order_lower_limit=float(limit),
        latent_order_midpoint=float(poly.rho(poly.midpoint())),
        panel_resolution=float(poly.resolution()),
        latent_order_gap=float(gap),
        latent_order_sharp=readings.sharp,
        latent_order_readings=int(poly.Pr.size),
        latent_order_alpha=float(alpha),
        latent_order_splits=SPLITS,
        latent_order_upper_starts=int(starts + 1))
    if len(readings.names) > 1:
        within = {}
        for s, name in enumerate(readings.names):
            rows = np.flatnonzero(readings.stratum == s)
            lineages, counts = np.unique(readings.code[rows], return_counts=True)
            keep = np.isin(readings.code[rows], lineages[counts >= 2])
            chosen = lineages[counts >= 2]
            within[str(name)] = (float(_lower(readings.polytope(rows[keep], lineages=chosen))[0])
                                 if chosen.size >= 2 else None)
        record["latent_order_lower_within_strata"] = within
    return record


def _empty(alpha: float, starts: int) -> Dict[str, Any]:
    nan = float("nan")
    return dict(latent_order_lower=nan, latent_order_upper=nan, latent_order_upper_bound=nan,
                latent_order_upper_exact=False, latent_order_lower_limit=nan,
                latent_order_midpoint=nan, panel_resolution=nan, latent_order_gap=nan,
                latent_order_sharp=False, latent_order_readings=0,
                latent_order_alpha=float(alpha), latent_order_splits=SPLITS,
                latent_order_upper_starts=int(starts + 1))


def _repeated(lineage) -> np.ndarray:
    labels = np.asarray([str(v) for v in lineage], dtype=object)
    _, code, sizes = np.unique(labels, return_inverse=True, return_counts=True)
    return sizes[code.ravel()] >= 2


def order_bounds(lo, hi, lineage, *, strata=None, alpha: float = 0.05,
                 starts: int = 10, seed=0) -> Dict[str, Any]:
    """Bounds on the lineage share of the latent MIC ordering from readings.

    ``lo``, ``hi`` and ``strata`` are as for
    :func:`amr_clonalshare.mic_order.mic_order_share`. Readings without a
    value, unbounded readings and untyped isolates are left out, and
    singleton lineages are set aside, as they are for the share."""
    from .attribution import _is_untyped
    if not 0. < alpha < 1.:
        raise ValueError("alpha must satisfy 0 < alpha < 1")
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    lineage = np.asarray(lineage, dtype=object)
    strata = np.zeros(lo.size, dtype=object) if strata is None else \
        np.asarray([str(s) for s in strata], dtype=object)
    keep = ~(np.isnan(lo) | np.isnan(hi)) & ~(np.isneginf(lo) & np.isposinf(hi)) \
        & np.array([not _is_untyped(v) for v in lineage], dtype=bool)
    lo, hi, lineage, strata = lo[keep], hi[keep], lineage[keep], strata[keep]
    rep = _repeated(lineage)
    lo, hi, lineage, strata = lo[rep], hi[rep], lineage[rep], strata[rep]
    if lineage.size == 0 or len(set(map(str, lineage))) < 2:
        return _empty(alpha, starts)
    _, code = np.unique(np.asarray([str(v) for v in lineage], dtype=object), return_inverse=True)
    return bounds_from_bins(Readings(lo, hi, code.ravel(), strata), alpha=alpha,
                            starts=starts, seed=seed)


def call_bounds(y, lineage, *, strata=None, alpha: float = 0.05,
                starts: int = 10, seed=0) -> Dict[str, Any]:
    """The same bounds from a binary call, a reading with one cut point.

    A call of 1 is the reading above the cut point and 0 the reading below
    it, so the bounds are those of a panel with one cut point, and they
    contain the bounds of any panel the call was cut from."""
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(y)
    if np.any(ok & ~np.isin(y, (0., 1.))):
        raise ValueError("calls must be 0, 1 or missing")
    lo = np.where(y == 1., 0., -np.inf)
    hi = np.where(y == 1., np.inf, 0.)
    lo[~ok] = np.nan
    hi[~ok] = np.nan
    return order_bounds(lo, hi, lineage, strata=strata, alpha=alpha,
                        starts=starts, seed=seed)
