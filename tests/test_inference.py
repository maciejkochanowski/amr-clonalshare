"""Validity of the post-clustering inference: merging rules, feature splitting,
count splitting, and the calibration claims the manuscript relies on."""
import numpy as np
import pandas as pd
import pytest
from sklearn.cluster import KMeans

from amr_clonalshare.inference import (binomial_thin, feature_split_report,
                                          feature_split_test, merge_exchangeable,
                                          merge_pvalues, nb_thin,
                                          ruger_merge, thinning_dependence)


def _km(block, k, seed):
    return KMeans(n_clusters=k, n_init=5, random_state=seed).fit(
        np.asarray(block, dtype=float)).labels_


# ----------------------------------------------------------------- merging --
def test_ruger_matches_its_definition():
    p = [0.02, 0.30, 0.44, 0.61, 0.90]
    # K = 5, k = ceil(5/2) = 3 -> (5/3) * p_(3) = (5/3) * 0.44
    assert ruger_merge(p) == pytest.approx(min(1.0, (5 / 3) * 0.44))


def test_exchangeable_ruger_dominates_classical_ruger():
    rng = np.random.default_rng(0)
    for _ in range(300):
        p = rng.random(9)
        assert merge_exchangeable(p) <= ruger_merge(p) + 1e-12


def test_exchangeable_ruger_formula_on_a_worked_example():
    p = [0.4, 0.1, 0.8, 0.2, 0.6]
    K, k = 5, 3
    best = min(np.sort(p[:l])[int(np.ceil(l * k / K)) - 1] for l in range(1, K + 1))
    assert merge_exchangeable(p) == pytest.approx(min(1.0, (K / k) * best))


@pytest.mark.slow
@pytest.mark.parametrize("rho", [0.0, 0.5, 0.9])
def test_merging_rules_hold_their_level_under_exchangeable_dependence(rho):
    """Every merging rule must hold its level on exchangeable p-values.

    The null is an equicorrelated Gaussian copula: ``p_m = Phi(sqrt(rho) W +
    sqrt(1-rho) e_m)`` with ``W, e_m ~ N(0,1)``. Each ``p_m`` is then exactly
    uniform (which a merging rule requires) while the ``p_m`` are exchangeable
    and, at ``rho = 0.9``, strongly dependent - the regime the rules are chosen
    for, since repeated splits of one fixed dataset are dependent through the
    data.
    """
    from scipy.stats import norm
    rng = np.random.default_rng(1)
    M, B = 9, 3000
    W = rng.standard_normal(B)[:, None]
    E = rng.standard_normal((B, M))
    P = norm.cdf(np.sqrt(rho) * W + np.sqrt(1 - rho) * E)
    for name in ("exchangeable_ruger", "ruger", "twice_mean", "mmb"):
        rate = np.mean([merge_pvalues(P[b], name) < 0.05 for b in range(B)])
        assert rate <= 0.065, f"{name} at rho={rho} rejected {rate:.3f}"


def test_merging_rules_are_marginally_valid_at_rho_zero():
    """Fast smoke version of the calibration test above."""
    from scipy.stats import norm
    rng = np.random.default_rng(11)
    P = norm.cdf(rng.standard_normal((600, 9)))
    for name in ("exchangeable_ruger", "ruger", "twice_mean", "mmb"):
        rate = np.mean([merge_pvalues(row, name) < 0.05 for row in P])
        assert rate <= 0.08, f"{name} rejected {rate:.3f}"


def test_merge_rejects_an_unknown_rule():
    with pytest.raises(ValueError):
        merge_pvalues([0.1, 0.2], "minimum")


# ---------------------------------------------------------- count splitting --
def test_binomial_hypergeometric_split_is_independent():
    rng = np.random.default_rng(2)
    X = rng.binomial(11, 0.3, size=(1500, 4))
    X1, X2 = binomial_thin(X, m=11, rng=rng)
    assert np.array_equal(X1 + X2, X)
    dep = thinning_dependence(X1, X2)
    assert np.nanmax(np.abs(dep)) < 0.06


def test_nb_split_is_independent_with_the_true_r_and_dependent_without_it():
    rng = np.random.default_rng(3)
    r_true = 5.0
    X = rng.negative_binomial(r_true, r_true / (r_true + 4.0), size=(1500, 4))
    _, dep_ok = (lambda a: (a, thinning_dependence(*a)))(
        nb_thin(X, r=r_true, rng=rng))
    assert np.nanmax(np.abs(dep_ok)) < 0.06
    # Neufeld et al. 2024 Prop. 11: a wrong r induces a non-zero covariance
    X1, X2 = nb_thin(X, r=0.4, rng=rng)
    dep_bad = thinning_dependence(X1, X2)
    assert np.nanmax(np.abs(dep_bad)) > 0.15


def test_nb_thin_validates_its_arguments():
    X = np.ones((5, 3), dtype=int)
    with pytest.raises(ValueError):
        nb_thin(X, r=1.0, eps=0.0)
    with pytest.raises(ValueError):
        nb_thin(X, r=-1.0)
    with pytest.raises(ValueError):
        nb_thin(X, r=[1.0, 2.0])          # wrong length


# --------------------------------------------------------- feature splitting --
@pytest.mark.slow
def test_feature_split_holds_its_level_under_a_binary_null():
    """No clusters: the merged p-value must reject at most 5% of the time."""
    rng = np.random.default_rng(4)
    rej = 0
    B = 120
    for _ in range(B):
        prev = rng.uniform(0.1, 0.6, size=24)
        X = (rng.random((120, 24)) < prev[None, :]).astype(int)
        r = feature_split_test(X, 2, cluster_fn=_km, rng=rng, n_splits=9)
        rej += r["p_value"] < 0.05
    assert rej / B <= 0.08, f"observed rejection rate {rej / B:.3f}"


def test_feature_split_has_power_against_real_structure():
    rng = np.random.default_rng(5)
    n, p = 150, 24
    g = rng.integers(0, 2, size=n)
    prev = np.full((n, p), 0.2)
    prev[g == 1, :p // 2] = 0.8
    X = (rng.random((n, p)) < prev).astype(int)
    r = feature_split_test(X, 2, cluster_fn=_km, rng=rng, n_splits=15)
    assert r["p_value"] < 1e-3


def test_feature_split_report_returns_a_usable_table():
    rng = np.random.default_rng(6)
    n, p = 120, 20
    g = rng.integers(0, 2, size=n)
    prev = np.full((n, p), 0.2)
    prev[g == 1, :6] = 0.85
    X = pd.DataFrame((rng.random((n, p)) < prev).astype(int),
                     columns=[f"f{j}" for j in range(p)])
    out = feature_split_report(X, 2, cluster_fn=_km, rng=rng, n_splits=15)
    assert out["status"] == "ok"
    tab = pd.DataFrame(out["per_feature"])
    assert len(tab) == p
    assert set(tab.columns) >= {"feature", "p_value", "p_value_bh", "is_defining",
                                "n_splits_tested"}
    # the six shifted features should dominate the ranking
    top6 = set(tab.nsmallest(6, "p_value")["feature"])
    assert len(top6 & {f"f{j}" for j in range(6)}) >= 5


def test_feature_split_declines_when_there_are_too_few_features():
    rng = np.random.default_rng(7)
    X = rng.integers(0, 2, size=(50, 3))
    assert feature_split_test(X, 2, cluster_fn=_km, rng=rng)["status"] == "skipped"


def test_single_split_is_far_more_variable_than_the_merged_verdict():
    """The reason multi-split exists: one draw is not reproducible."""
    rng = np.random.default_rng(8)
    prev = rng.uniform(0.1, 0.6, size=24)
    X = (rng.random((150, 24)) < prev[None, :]).astype(int)
    singles, merged = [], []
    for s in range(10):
        r = feature_split_test(X, 2, cluster_fn=_km,
                               rng=np.random.default_rng(100 + s), n_splits=9)
        singles.append(r["p_per_split"][0])
        merged.append(r["p_value"])
    assert np.ptp(singles) > np.ptp(merged)


# ----------------------------------------------------- discreteness vs gradient --
def _planted_discrete(rng, n=200, p=24, k=3, noise=0.05):
    sizes = np.full(k, n // k)
    sizes[-1] += n - sizes.sum()
    truth = np.concatenate([np.full(s, c) for c, s in enumerate(sizes)])
    X = np.zeros((n, p), dtype=int)
    per = p // k
    for c in range(k):
        X[truth == c, c * per:(c + 1) * per] = 1
    return np.clip(X + rng.binomial(1, noise, X.shape), 0, 1)


def _one_factor_data(rng, n=200, p=24):
    a = rng.normal(-0.5, 1.0, p)
    b = rng.normal(1.2, 0.4, p)
    z = rng.standard_normal(n)
    P = 1 / (1 + np.exp(-(a[None, :] + z[:, None] * b[None, :])))
    return (rng.random((n, p)) < P).astype(int)


def test_one_factor_loglik_is_a_proper_marginal_likelihood():
    from amr_clonalshare.inference import fit_one_factor, one_factor_loglik
    rng = np.random.default_rng(20)
    X = _one_factor_data(rng)
    a, b, _ = fit_one_factor(X)
    ll = one_factor_loglik(X, a, b)
    assert np.isfinite(ll) and ll < 0
    # a model with no discrimination is strictly worse on data that has some
    ll_flat = one_factor_loglik(X, a, np.zeros_like(b))
    assert ll > ll_flat


def test_fit_one_factor_bounds_its_coefficients():
    """Unbounded slopes made the bootstrap null useless (sd 2e6)."""
    from amr_clonalshare.inference import fit_one_factor
    rng = np.random.default_rng(21)
    X = _planted_discrete(rng, n=120, p=18, k=2, noise=0.0)   # separable
    a, b, z = fit_one_factor(X)
    assert np.abs(b).max() <= 6.0 + 1e-9
    assert np.abs(a).max() <= 6.0 + 1e-9
    assert np.isfinite(z).all()


def test_continuum_null_detects_discrete_structure():
    from amr_clonalshare.inference import continuum_null_test
    rng = np.random.default_rng(22)
    out = continuum_null_test(_planted_discrete(rng), 3, rng=rng, n_boot=49)
    assert out["status"] == "ok"
    assert out["discrete_beyond_a_gradient"] is True
    assert out["observed"] > out["null_mean"]


def test_continuum_null_does_not_fire_on_a_gradient():
    """The whole point: a smooth latent gradient must not be called discrete."""
    from amr_clonalshare.inference import continuum_null_test
    rng = np.random.default_rng(23)
    out = continuum_null_test(_one_factor_data(rng, n=250), 3, rng=rng, n_boot=49)
    assert out["status"] == "ok"
    assert out["discrete_beyond_a_gradient"] is False


def test_continuum_null_reports_its_p_value_floor():
    from amr_clonalshare.inference import continuum_null_test
    rng = np.random.default_rng(24)
    out = continuum_null_test(_planted_discrete(rng, n=100, p=18, k=2), 2,
                              rng=rng, n_boot=19)
    assert out["p_value"] >= out["p_value_floor"] == pytest.approx(1 / 20)
    assert out["bootstrap_exceedances"] >= 0
    ci = out["bootstrap_tail_probability_ci95"]
    assert ci["method"] == "Clopper-Pearson exact binomial"
    assert 0 <= ci["low"] <= ci["high"] <= 1


def test_continuum_zero_exceedances_reports_monte_carlo_resolution(monkeypatch):
    """A floor p-value must expose finite-bootstrap uncertainty."""
    import amr_clonalshare.inference as inference
    import amr_clonalshare.baselines as baselines

    monkeypatch.setattr(
        inference,
        "select_latent_dimension",
        lambda X, q_max: {"q_selected": 1, "bic_by_q": {"1": 0.0},
                          "q_max": q_max, "at_boundary": False},
    )
    calls = iter([10.0] + [0.0] * 99)
    monkeypatch.setattr(
        inference, "fit_latent_trait",
        lambda X, q: (np.zeros(X.shape[1]), np.zeros((X.shape[1], q)),
                      np.zeros((len(X), q))))
    monkeypatch.setattr(inference, "latent_trait_loglik", lambda *a, **k: 0.0)
    monkeypatch.setattr(
        baselines, "bernoulli_mixture",
        lambda *a, **k: {"bic": -next(calls)},
    )
    out = inference.continuum_null_test(
        np.zeros((20, 3)), 2, rng=np.random.default_rng(242), n_boot=99, q=1)
    assert out["bootstrap_exceedances"] == 0
    assert out["p_value"] == pytest.approx(0.01)
    assert out["bootstrap_tail_probability_ci95"]["high"] == pytest.approx(
        0.03657574498347894)
    assert out["resolved_at_alpha_0_05"] is True


def test_continuum_boundary_is_withheld_before_bootstrap(monkeypatch):
    import amr_clonalshare.inference as inference

    monkeypatch.setattr(
        inference,
        "select_latent_dimension",
        lambda X, q_max: {"q_selected": q_max, "bic_by_q": {str(q_max): 1.0},
                          "q_max": q_max, "at_boundary": True},
    )
    out = inference.continuum_null_test(
        _planted_discrete(np.random.default_rng(240), n=80, p=12, k=2),
        2, rng=np.random.default_rng(241), n_boot=99, q_max=4)
    assert out["status"] == "withheld_under_dimensioned"
    assert out["discrete_beyond_a_gradient"] is None
    assert out["p_value"] is None
    assert out["n_boot"] == 0


def test_latent_dimension_extends_only_from_a_boundary(monkeypatch):
    import amr_clonalshare.inference as inference

    fitted = []

    def fake_fit(X, q, n_iter=30):
        fitted.append(q)
        return np.zeros(X.shape[1]), np.full((X.shape[1], q), q), np.zeros((len(X), q))

    # BIC is minimal at q=2, so q=4 must not be fitted.
    monkeypatch.setattr(inference, "fit_latent_trait", fake_fit)
    monkeypatch.setattr(
        inference,
        "latent_trait_loglik",
        lambda X, a, B, n_nodes=None: {1: -140.0, 2: -100.0, 3: -110.0, 4: -90.0}[B.shape[1]],
    )
    out = inference.select_latent_dimension(np.zeros((100, 5)), q_max=4)
    assert fitted == [1, 2, 3]
    assert out["q_selected"] == 2
    assert out["q_evaluated_max"] == 3
    assert out["adaptive_extension_used"] is False

    # If q=3 is provisionally best, q=4 is required and can close the curve.
    fitted.clear()
    monkeypatch.setattr(
        inference,
        "latent_trait_loglik",
        lambda X, a, B, n_nodes=None: {1: -160.0, 2: -130.0, 3: -100.0, 4: -110.0}[B.shape[1]],
    )
    out = inference.select_latent_dimension(np.zeros((100, 5)), q_max=4)
    assert fitted == [1, 2, 3, 4]
    assert out["q_selected"] == 3
    assert out["q_evaluated_max"] == 4
    assert out["adaptive_extension_used"] is True
    assert out["at_boundary"] is False


@pytest.mark.slow
def test_block_aware_split_is_calibrated_where_the_naive_one_is_not():
    """The headline calibration result, at a size the fast lane can afford."""
    from amr_clonalshare.inference import feature_split_test
    rng = np.random.default_rng(25)
    n, p, bs, rho = 150, 30, 5, 0.9
    naive, aware = [], []
    groups = np.repeat(np.arange(p // bs), bs)
    for _ in range(25):
        prev = rng.uniform(0.15, 0.5, size=p)
        X = np.zeros((n, p), dtype=int)
        for start in range(0, p, bs):
            drv = rng.random(n) < prev[start]
            for j in range(start, min(start + bs, p)):
                flip = rng.random(n) < (1 - rho)
                X[:, j] = np.where(flip, rng.random(n) < prev[j], drv)
        naive.append(feature_split_test(X, 2, cluster_fn=_km, rng=rng,
                                        n_splits=9)["p_value"])
        aware.append(feature_split_test(X, 2, cluster_fn=_km, rng=rng,
                                        n_splits=9, groups=groups)["p_value"])
    naive_rate = float(np.mean(np.asarray(naive) < 0.05))
    aware_rate = float(np.mean(np.asarray(aware) < 0.05))
    assert naive_rate > 0.5, f"expected the naive split to break; got {naive_rate}"
    assert aware_rate <= 0.12, f"block-aware split rejected {aware_rate}"


def test_merge_exchangeable_refuses_inputs_that_are_not_probabilities():
    """This is the default merge rule feeding Benjamini-Hochberg.
    ``[-1.0, 0.1, 0.1]`` merged to -1.5, a negative p-value, and
    ``[0.01, nan, 0.3]`` merged to 0.015 with the NaN sorted to the end of
    every prefix and silently dropped from the minimum. Both copies of the
    rule -- the one in `inference` and the one in `tva` -- must refuse."""
    from amr_clonalshare.tva import merge_exchangeable as tva_merge

    for merge in (merge_exchangeable, tva_merge):
        for bad in ([-1.0, 0.1, 0.1], [0.01, float("nan"), 0.3],
                    [1.4, 0.1, 0.1]):
            with pytest.raises(ValueError, match="not a probability"):
                merge(bad)
        assert 0.0 <= merge([0.01, 0.2, 0.3]) <= 1.0
