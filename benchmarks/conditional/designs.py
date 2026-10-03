"""Designs with the lineages, their sizes and their distributions held fixed.

Every design here fixes the distribution of the readings of every lineage, so the
target of an interval for the represented lineages, the share of the lineages
held fixed, is known exactly, and a replicate draws only the isolates again.
Four families:

``observed``    distributions over the wells of a panel, per lineage and laboratory,
                with the laboratory of every isolate drawn with its reading
                in the proportions of its lineage: the interval for the share
                of the represented lineages (``mic_order_share``; a design of two
                wells is a call, which the package reads identically).
``latent``      a latent log2 MIC per lineage (normal, t with four degrees
                of freedom, or a contaminated normal), read on the panel of a
                laboratory: the bounds on the share of the latent
                ordering and their lower confidence limit.
``units``       farms that shift the readings of their isolates, lineages
                that sit mostly or wholly on one farm, and no lineage effect
                within farms: the level of the test within sampling units.
``comparison``  records with two lineage labellings and outcomes drawn from
                the distribution of the records sharing both labels: the paired
                intervals of the comparison of two lineage definitions.

Nothing here is used by the package; the truths are written from their
definitions and share no code with the estimates.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import stats


# --------------------------------------------------------------------------- #
# observed: fixed distributions over wells
# --------------------------------------------------------------------------- #

@dataclass
class ObservedDesign:
    """laws[g][s]: probabilities over the K_s wells of laboratory s;
    sizes[g][s]: expected isolates of lineage g in laboratory s. The
    laboratory of an isolate is drawn with its reading, in the lineage's
    proportions (random membership), the lineage's total size fixed."""
    name: str
    laws: list
    sizes: np.ndarray
    wells: list
    eta: float = field(init=False)

    def __post_init__(self):
        self.sizes = np.asarray(self.sizes, dtype=int)
        if self.sizes.ndim == 1:
            self.sizes = self.sizes[:, None]
        self.G, self.S = self.sizes.shape
        self.laws = [[np.asarray(self.laws[g][s], float) / np.sum(self.laws[g][s]) for s in range(self.S)]
                     for g in range(self.G)]
        self.eta = self._truth()

    def _truth(self) -> float:
        """The share of the mid-distribution positions within laboratory,
        the laboratory distribution being the mixture of the lineages' distributions at the
        expected counts, over the lineages held fixed."""
        n = self.sizes.sum()
        mu = np.zeros(self.G)
        m2 = np.zeros(self.G)
        for s in range(self.S):
            ns = self.sizes[:, s].sum()
            if ns == 0:
                continue
            P = sum(self.sizes[g, s] * self.laws[g][s] for g in range(self.G)) / ns
            score = np.r_[0, np.cumsum(P)][:-1] + P / 2
            for g in range(self.G):
                f = self.sizes[g, s] / self.sizes[g].sum()
                mu[g] += f * self.laws[g][s] @ score
                m2[g] += f * self.laws[g][s] @ score ** 2
        w = self.sizes.sum(1) / n
        # lineages of one distribution have one mean: the share is zero exactly, not
        # the rounding of their weighted spread
        between = 0. if np.ptp(mu) == 0 else w @ (mu - w @ mu) ** 2
        return float(between / (between + w @ (m2 - mu ** 2)))

    def draw(self, rng):
        lo, hi, lin, lab = [], [], [], []
        for g in range(self.G):
            total = self.sizes[g].sum()
            counts = rng.multinomial(total, self.sizes[g] / total) if self.S > 1 else self.sizes[g]
            for s in range(self.S):
                m = counts[s]
                if m == 0:
                    continue
                K = self.wells[s]
                r = rng.choice(K, size=m, p=self.laws[g][s])
                lo.append(np.where(r == 0, -np.inf, r - 1.))
                hi.append(np.where(r == K - 1, np.inf, r * 1.))
                lin += [g] * m
                lab += [s] * m
        return np.concatenate(lo), np.concatenate(hi), np.array(lin), np.array(lab)


def _dominant(K, k, p):
    v = np.full(K, (1 - p) / (K - 1))
    v[k] = p
    return v


def _normal_wells(shift, K, spread=1.):
    cuts = stats.norm.ppf(np.arange(1, K) / K) * spread
    return np.diff(stats.norm.cdf(np.r_[-np.inf, cuts, np.inf] - shift))


def observed_designs() -> list:
    D = []
    add = lambda *a: D.append(ObservedDesign(*a))
    # near-deterministic lineages on few wells, the audit's two-well design
    add('audit_10x4_two_wells', [[_dominant(2, g % 2, .95)] for g in range(10)], [4] * 10, [2])
    add('k3_10x4_alternate', [[_dominant(3, 2 * (g % 2), .9)] for g in range(10)], [4] * 10, [3])
    add('k5_10x4_spread', [[_dominant(5, g % 5, .95)] for g in range(10)], [4] * 10, [5])
    add('k10_30x3_spread', [[_dominant(10, g % 10, .9)] for g in range(30)], [3] * 30, [10])
    add('k5_6x2_spread', [[_dominant(5, g % 5, .99)] for g in range(6)], [2] * 6, [5])
    add('k5_20x2_spread', [[_dominant(5, g % 5, .95)] for g in range(20)], [2] * 20, [5])
    add('k8_4x6_top_bottom', [[_dominant(8, 7 * (g % 2), .95)] for g in range(4)], [6] * 4, [8])
    add('k3_5x3_spread', [[_dominant(3, g % 3, .97)] for g in range(5)], [3] * 5, [3])
    # two laboratories on different panels, random membership
    add('two_labs', [[_dominant(4, g % 4, .9), _dominant(6, (g * 2) % 6, .9)] for g in range(12)],
        [[3, 2]] * 12, [4, 6])
    # no or weak lineage effect
    add('no_effect_one_well', [[_dominant(6, 2, .9)] for g in range(15)], [4] * 15, [6])
    add('null_15x5', [[_normal_wells(0., 8)] for g in range(15)], [5] * 15, [8])
    add('weak_20x6', [[_normal_wells(m, 8)] for m in np.linspace(-.4, .4, 20)], [6] * 20, [8])
    add('moderate_8x25', [[_normal_wells(m, 8)] for m in np.linspace(-1, 1, 8)], [25] * 8, [8])
    # uneven sizes, pairs
    add('mixed_sizes', [[_dominant(6, g % 6, .8)] for g in range(25)], [20] * 5 + [2] * 20, [6])
    add('pairs_5', [[_normal_wells(m, 12, 1.5)] for m in np.linspace(-1.5, 1.5, 5)], [2] * 5, [12])
    add('pairs_10', [[_normal_wells(m, 12, 1.5)] for m in np.linspace(-1.5, 1.5, 10)], [2] * 10, [12])
    # pairs with a weak effect, on twelve wells and on a practically
    # continuous panel of a thousand: the within-lineage variance rests on as
    # many degrees of freedom as there are pairs
    for G in (5, 10, 30):
        add(f'weak_pairs_{G}_fine', [[_normal_wells(m, 1000)] for m in np.linspace(-.45, .45, G)], [2] * G, [1000])
    for G in (10, 30):
        add(f'weak_pairs_{G}_12_wells', [[_normal_wells(m, 12)] for m in np.linspace(-.45, .45, G)], [2] * G, [12])
    add('pairs_100_fine', [[_normal_wells(m, 1000)] for m in np.linspace(-.9, .9, 100)], [2] * 100, [1000])
    # a rare call carried by one or three lineages (the sparse regime)
    for m in (5, 10, 20):
        for p0 in (.02, .05):
            add(f'rare_call_one_m{m}_p{p0}', [[np.array([1 - p0, p0])]] * 19 + [[np.array([.5, .5])]],
                [m] * 20, [2])
    for m in (5, 10):
        add(f'rare_call_three_m{m}', [[np.array([.98, .02])]] * 17 + [[np.array([.5, .5])]] * 3, [m] * 20, [2])
    # a shifted lineage on a panel
    base, shifted = _normal_wells(-1.2, 8), _normal_wells(1.2, 8)
    for m in (5, 10):
        add(f'one_shifted_m{m}', [[base]] * 19 + [[shifted]], [m] * 20, [8])
    return D


# --------------------------------------------------------------------------- #
# latent: a latent value per isolate, read on a panel
# --------------------------------------------------------------------------- #

LAWS = {
    'normal': (lambda rng, n: rng.normal(size=n), stats.norm.cdf),
    't4': (lambda rng, n: rng.standard_t(4, n) / np.sqrt(2.), lambda x: stats.t.cdf(x * np.sqrt(2.), 4)),
    'contaminated': (lambda rng, n: rng.normal(size=n) * np.where(rng.random(n) < .1, 4., 1.) / np.sqrt(2.5),
                     lambda x: .9 * stats.norm.cdf(x * np.sqrt(2.5)) + .1 * stats.norm.cdf(x * np.sqrt(2.5) / 4.)),
}


@dataclass
class LatentDesign:
    """mu[g]: latent shift of lineage g; sizes[g][s]: fixed isolates of g in
    laboratory s (fixed membership); cuts[s]: the panel of laboratory s;
    offsets[s]: its shift; distribution: the residual distribution; overlap: laboratory 0
    reads every second isolate on its panel moved half a well, so that one
    laboratory holds readings that overlap."""
    name: str
    mu: np.ndarray
    sizes: np.ndarray
    cuts: list
    law: str = 'normal'
    offsets: tuple = ()
    overlap: bool = False

    def __post_init__(self):
        self.mu = np.asarray(self.mu, float)
        self.sizes = np.asarray(self.sizes, int)
        if self.sizes.ndim == 1:
            self.sizes = self.sizes[:, None]
        self.G, self.S = self.sizes.shape
        self.cuts = [np.asarray(c, float) for c in self.cuts]
        self.off = np.zeros(self.S) if not self.offsets else np.asarray(self.offsets, float)
        self.rho, self.rho_min = self._truth()

    def _cdf(self, x):
        return LAWS[self.law][1](x)

    def _truth(self):
        """rho by quadrature of P(X_h < X_g), and rho_min from the polytope
        of the population reading probabilities (one panel per laboratory;
        with overlapping readings rho_min is not computed here)."""
        from amr_clonalshare.latent_order import _Polytope, _lower
        n = self.sizes.sum()
        w = self.sizes.sum(1) / n
        x = np.linspace(-12, 12, 24001)
        dens = np.gradient(self._cdf(x), x)
        EU = np.zeros(self.G)
        for s in range(self.S):
            ns = self.sizes[:, s].sum()
            if ns == 0:
                continue
            # P(X_h < X_g) = integral F(t - mu_h) f(t - mu_g) dt
            P = np.array([[np.trapezoid(self._cdf(x + self.mu[g] - self.mu[h]) * dens, x)
                           for h in range(self.G)] for g in range(self.G)])
            EU += self.sizes[:, s] / n * (P @ (self.sizes[:, s] / ns))
        # lineages of one shift share one distribution: rho is zero exactly, not the
        # error of the quadrature
        rho = 0. if np.ptp(self.mu) == 0 else float(12 * np.sum(w * (EU / w - .5) ** 2))
        if self.overlap:
            return rho, float('nan')
        blocks, L, H = [], [], []
        for s in range(self.S):
            edges = np.r_[-np.inf, self.cuts[s], np.inf]
            q = np.diff(self._cdf(edges[None, :] - self.mu[:, None] - self.off[s]), axis=1)
            counts = self.sizes[:, s][:, None] * q
            blocks.append(counts)
            P = counts.sum(0)
            cum = np.r_[0., np.cumsum(P / P.sum())] if P.sum() > 0 else np.zeros(P.size + 1)
            L.append(cum[:-1]); H.append(cum[1:])
        counts, L, H = np.hstack(blocks), np.concatenate(L), np.concatenate(H)
        keep = counts.sum(0) > 1e-300
        return rho, float(_lower(_Polytope(counts[:, keep], L[keep], H[keep]))[0])

    def draw(self, rng):
        """Readings (lo, hi], lineage, laboratory, and the latent values."""
        sample = LAWS[self.law][0]
        lo, hi, lin, lab, lat = [], [], [], [], []
        for g in range(self.G):
            for s in range(self.S):
                m = self.sizes[g, s]
                if m == 0:
                    continue
                x = self.mu[g] + self.off[s] + sample(rng, m)
                cuts = np.repeat(self.cuts[s][None, :], m, axis=0)
                if self.overlap and s == 0:
                    cuts[1::2] += .5 * np.diff(self.cuts[s]).mean()
                k = (x[:, None] > cuts).sum(1)
                ext = np.hstack([np.full((m, 1), -np.inf), cuts, np.full((m, 1), np.inf)])
                lo.append(ext[np.arange(m), k]); hi.append(ext[np.arange(m), k + 1])
                lin += [g] * m; lab += [s] * m; lat.append(x)
        return (np.concatenate(lo), np.concatenate(hi), np.array(lin), np.array(lab),
                np.concatenate(lat))


def latent_designs() -> list:
    u = lambda K: stats.norm.ppf(np.arange(1, K) / K) * np.sqrt(2.)
    D = []
    add = lambda *a, **k: D.append(LatentDesign(*a, **k))
    # the audit's nulls: a median call on two lineages of two, and many pairs
    # on a coarse and on a practically continuous panel
    add('null_2x2_median_call', np.zeros(2), [2, 2], [u(2)])
    add('null_30x2_9_wells', np.zeros(30), [2] * 30, [u(9)])
    add('null_30x2_1000_wells', np.zeros(30), [2] * 30, [u(1000)])
    add('null_100x2_1000_wells', np.zeros(100), [2] * 100, [u(1000)])
    add('null_10x10_12_wells', np.zeros(10), [10] * 10, [u(12)])
    add('near_null_30x4', np.linspace(-.3, .3, 30), [4] * 30, [u(12)])
    add('effect_8x10', np.linspace(-1, 1, 8), [10] * 8, [u(10)])
    add('effect_3x20', np.array([-1., 0, 1]), [20] * 3, [u(6)])
    add('call_10x4', np.tile([-1.645, 1.645], 5), [4] * 10, [np.array([0.])])
    add('two_labs_confined', np.linspace(-.8, .8, 12), np.array([[4, 0]] * 6 + [[0, 4]] * 6),
        [u(8), u(5)], offsets=(-.5, .5))
    for m in (5, 10, 20):
        add(f'rare_call_one_m{m}', [-2.054] * 19 + [0.], [m] * 20, [np.array([0.])])
    add('rare_call_three_m5', [-2.054] * 17 + [0.] * 3, [5] * 20, [np.array([0.])])
    add('rare_one_m10_8_wells', [-2.054] * 19 + [0.], [10] * 20, [np.linspace(-3, 1, 7)])
    for law in ('t4', 'contaminated'):
        add(f'effect_8x10_{law}', np.linspace(-1, 1, 8), [10] * 8, [u(10)], law=law)
        add(f'null_30x2_9_wells_{law}', np.zeros(30), [2] * 30, [u(9)], law=law)
    add('overlapping_panels_effect', np.linspace(-1, 1, 8), [10] * 8, [np.arange(-3., 3.5, 1.)], overlap=True)
    add('overlapping_panels_null', np.zeros(15), [4] * 15, [np.arange(-3., 3.5, 1.)], overlap=True)
    return D


def sample_latent_share(latent, lineage, lab):
    """The share of the ordering of the latent values of the sample: 12 sum_g w_g
    (mean position - 1/2)^2 over the lineages of two or more isolates, a
    position being the share of those isolates of the same laboratory below
    the value plus half the share at it."""
    keep = np.bincount(lineage)[lineage] >= 2
    lineage, lab, latent = lineage[keep], lab[keep], latent[keep]
    u = np.empty(latent.size)
    for k in np.unique(lab):
        m = np.flatnonzero(lab == k)
        order = np.argsort(latent[m], kind='stable')
        rank = np.empty(m.size)
        rank[order] = np.arange(m.size)
        u[m] = (rank + .5) / m.size
    code = np.unique(lineage, return_inverse=True)[1]
    size = np.bincount(code).astype(float)
    means = np.bincount(code, weights=u) / size
    return float(12. * np.sum(size / size.sum() * (means - .5) ** 2))


# --------------------------------------------------------------------------- #
# units: farms, no lineage effect within farms
# --------------------------------------------------------------------------- #

@dataclass
class UnitDesign:
    """``farms`` sampling units with effects of standard deviation
    ``farm_sd`` on the latent value; ``lineages`` lineages of ``size``
    isolates each (lognormal sizes if ``size`` is None), a share ``home`` of
    every lineage's isolates on its home farm and the rest spread; no lineage
    effect within farms. ``call`` cuts the latent value at its median into a
    call; otherwise it is read on a panel of eight wells. ``labs``
    laboratories, crossed with the farms, shift the readings of their
    isolates by one log2 step."""
    name: str
    farms: int
    lineages: int
    size: int | None
    home: float
    farm_sd: float = 1.
    call: bool = False
    labs: int = 1

    def draw(self, rng):
        if self.size is None:
            sizes = np.rint(np.exp(rng.normal(1.6, .8, self.lineages))).astype(int).clip(2, 60)
        else:
            sizes = np.full(self.lineages, self.size)
        home = rng.integers(0, self.farms, self.lineages)
        effect = rng.normal(0., self.farm_sd, self.farms)
        lineage = np.repeat(np.arange(self.lineages), sizes)
        farm = np.where(rng.random(lineage.size) < self.home, home[lineage],
                        rng.integers(0, self.farms, lineage.size))
        lab = rng.integers(0, self.labs, lineage.size)
        x = effect[farm] + 1. * lab + rng.normal(size=lineage.size)
        if self.call:
            cut = np.median(x)
            lo, hi = np.where(x > cut, 0., -np.inf), np.where(x > cut, np.inf, 0.)
        else:
            cuts = np.arange(-3., 4.)
            k = np.searchsorted(cuts, x)
            ext = np.r_[-np.inf, cuts, np.inf]
            lo, hi = ext[k], ext[k + 1]
        return lo, hi, lineage, farm, lab


def unit_designs() -> list:
    return [UnitDesign('farms10_home80_panel', 10, 20, 6, .8),
            UnitDesign('farms10_home80_call', 10, 20, 6, .8, call=True),
            UnitDesign('farms4_nested_uneven', 4, 40, None, 1., farm_sd=1.5),
            UnitDesign('farms30_home60_pairs', 30, 60, 2, .6),
            UnitDesign('farms6_home90_two_labs', 6, 24, 5, .9, labs=2)]


# --------------------------------------------------------------------------- #
# comparison: two labellings of the same records
# --------------------------------------------------------------------------- #

@dataclass
class ComparisonDesign:
    """Records of ``lineages_a`` lineages of ``size`` isolates under the first
    definition; the second merges consecutive lineages in groups of
    ``merge`` (split > 1 splits every lineage into ``split`` parts instead);
    a share ``missing_b`` of records has no second label. The rate of a
    record is set by its first-definition lineage: normal liabilities of
    share ``share`` at ``prevalence``, fixed for the design (one draw from
    ``seed``); outcomes are drawn again in every replicate."""
    name: str
    lineages_a: int
    size: int
    share: float
    prevalence: float = .25
    merge: int = 1
    split: int = 1
    missing_b: float = .1
    seed: int = 1

    def __post_init__(self):
        rng = np.random.default_rng(self.seed)
        n = self.lineages_a * self.size
        self.a = np.repeat(np.arange(self.lineages_a), self.size).astype(object)
        if self.split > 1:
            b = [f'{g}.{i % self.split}' for g, i in zip(np.repeat(np.arange(self.lineages_a), self.size),
                                                      np.tile(np.arange(self.size), self.lineages_a))]
        else:
            b = [f'{g // self.merge}' for g in np.repeat(np.arange(self.lineages_a), self.size)]
        self.b = np.asarray(b, dtype=object)
        self.b[rng.random(n) < self.missing_b] = None
        tau = np.sqrt(self.share / (1. - self.share)) if self.share > 0 else 0.
        effect = rng.normal(0., tau, self.lineages_a)
        c = stats.norm.ppf(self.prevalence) * np.sqrt(1. + tau * tau)
        self.p = stats.norm.cdf(c + effect)[np.repeat(np.arange(self.lineages_a), self.size)]
        self.ids = np.array([f'R{i:05d}' for i in range(n)])

    def draw(self, rng):
        return (rng.random(self.p.size) < self.p).astype(float)


def comparison_designs() -> list:
    return [ComparisonDesign('merge_pairs_share04', 60, 5, .4, merge=2),
            ComparisonDesign('merge_fives_share06', 60, 5, .6, merge=5),
            ComparisonDesign('split_in_two_share04', 40, 6, .4, split=2),
            ComparisonDesign('rare_merge_pairs', 60, 5, .5, prevalence=.08, merge=2),
            ComparisonDesign('null_merge_pairs', 60, 5, 0., merge=2),
            # added after the protocol, from the Salmonella example: many
            # lineages of three, most of them without a positive, merged in
            # fives, a rare outcome (CONFIRMATORY_PROTOCOL.md, Section 5)
            ComparisonDesign('small_rare_merge_fives', 200, 3, .9, prevalence=.08, merge=5)]
