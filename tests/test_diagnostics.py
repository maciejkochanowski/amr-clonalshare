"""Layer influence, lineage confounding, baselines: the diagnostics that decide
whether a partition means anything."""
import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import adjusted_rand_score

from amr_clonalshare.baselines import (bernoulli_mixture,
                                          bernoulli_mixture_select_k,
                                          external_agreement)
from amr_clonalshare.core import spectral_from_similarity
from amr_clonalshare.fusion import snf_fuse, snf_kernel
from amr_clonalshare.influence import (effective_n_layers, hill_number,
                                          layer_influence)
from amr_clonalshare.lineage import (cluster_composition, dereplicate_index,
                                        lineage_concordance)
from amr_clonalshare.stats import binary_distance_matrix


def _kernel(M):
    return snf_kernel(binary_distance_matrix(np.asarray(M)), K=10)


def _fuse(Ws):
    return snf_fuse(Ws, K=10, T=10)


def _blocks(n, p, k, rng, noise=0.05):
    """n isolates, k planted groups owning disjoint blocks of p features."""
    sizes = np.full(k, n // k)
    sizes[-1] += n - sizes.sum()
    truth = np.concatenate([np.full(s, c) for c, s in enumerate(sizes)])
    X = np.zeros((n, p), dtype=int)
    per = p // k
    for c in range(k):
        X[truth == c, c * per:(c + 1) * per] = 1
    return np.clip(X + rng.binomial(1, noise, X.shape), 0, 1), truth


# ------------------------------------------------------------- Hill numbers --
def test_hill_number_endpoints():
    assert hill_number([1.0, 0.0, 0.0]) == pytest.approx(1.0)
    assert hill_number([1 / 3] * 3) == pytest.approx(3.0)
    assert hill_number([0.5, 0.5], order=2) == pytest.approx(2.0)


def test_effective_n_layers_is_uniform_when_no_layer_moves_the_answer():
    out = effective_n_layers([0.0, 0.0, 0.0])
    assert out["n_eff"] == pytest.approx(3.0)


# ---------------------------------------------------------- layer influence --
def test_influence_detects_a_pure_noise_layer():
    rng = np.random.default_rng(0)
    signal, truth = _blocks(90, 12, 3, rng)
    noise = rng.integers(0, 2, size=(90, 12))
    out = layer_influence([signal, noise], ["signal", "noise"],
                          kernel_fn=_kernel, fuse_fn=_fuse,
                          cluster_fn=spectral_from_similarity, k=3, rng=rng)
    infl = {d["layer"]: d["delta_loo"] for d in out["per_layer"]}
    assert infl["signal"] > infl["noise"]
    assert out["dominant_layer"] == "signal"


def test_influence_reports_full_diversity_for_two_equally_informative_layers():
    rng = np.random.default_rng(1)
    a, truth = _blocks(90, 12, 3, rng)
    b, _ = _blocks(90, 12, 3, np.random.default_rng(2))
    out = layer_influence([a, b], ["a", "b"], kernel_fn=_kernel, fuse_fn=_fuse,
                          cluster_fn=spectral_from_similarity, k=3, rng=rng)
    assert out["n_eff"] >= 1.5
    assert out["collapse"] is False


def test_influence_permutation_test_runs_and_returns_valid_pvalues():
    rng = np.random.default_rng(3)
    a, _ = _blocks(60, 9, 3, rng)
    b = rng.integers(0, 2, size=(60, 9))
    out = layer_influence([a, b], ["a", "b"], kernel_fn=_kernel, fuse_fn=_fuse,
                          cluster_fn=spectral_from_similarity, k=3, rng=rng,
                          n_perm=5)
    for d in out["per_layer"]:
        assert 0 < d["p_contribution"] <= 1.0
        assert d["p_contribution"] >= 1 / 6      # Phipson-Smyth floor at B=5


def test_influence_handles_a_single_layer():
    rng = np.random.default_rng(4)
    a, _ = _blocks(40, 9, 3, rng)
    out = layer_influence([a], ["only"], kernel_fn=_kernel, fuse_fn=_fuse,
                          cluster_fn=spectral_from_similarity, k=3, rng=rng)
    assert out["n_layers"] == 1
    assert out["collapse"] is False


# ------------------------------------------------------------------ lineage --
def test_lineage_concordance_detects_a_pure_lineage_relabelling():
    lineage = [f"L{i // 5}" for i in range(100)]
    labels = [i // 25 for i in range(100)]      # nested inside lineages
    out = lineage_concordance(labels, lineage, n_perm=200,
                              rng=np.random.default_rng(0))
    assert out["status"] == "ok"
    assert out["concordance_observed"] == pytest.approx(1.0)
    assert out["z"] > 5
    assert out["p_value"] < 0.01


def test_lineage_concordance_is_null_when_labels_are_independent_of_lineage():
    rng = np.random.default_rng(1)
    lineage = [f"L{i % 20}" for i in range(200)]
    labels = rng.integers(0, 3, size=200)
    out = lineage_concordance(labels, lineage, n_perm=300, rng=rng)
    assert out["p_value"] > 0.05


def test_lineage_concordance_skips_when_every_lineage_is_a_singleton():
    out = lineage_concordance([0, 1, 2], ["a", "b", "c"], n_perm=10)
    assert out["status"] == "skipped"


def test_dereplicate_index_keeps_one_per_lineage():
    lineage = ["A", "A", "B", "B", "B", "C"]
    idx = dereplicate_index(lineage)
    assert len(idx) == 3
    assert sorted(np.asarray(lineage)[idx]) == ["A", "B", "C"]


def test_cluster_composition_reports_the_dominant_category():
    meta = pd.DataFrame({"ST": ["a", "a", "a", "b"]})
    out = cluster_composition([0, 0, 0, 1], meta, ["ST"])
    assert out["0"]["ST"]["max_share"] == pytest.approx(1.0)


# ---------------------------------------------------------------- baselines --
def test_bernoulli_mixture_recovers_planted_classes():
    rng = np.random.default_rng(5)
    X, truth = _blocks(150, 15, 3, rng, noise=0.03)
    fit = bernoulli_mixture(X, 3, rng=rng, n_init=8)
    assert adjusted_rand_score(truth, fit["labels"]) > 0.9
    assert fit["responsibilities"].shape == (150, 3)


def test_bic_selects_the_planted_k():
    rng = np.random.default_rng(6)
    X, _ = _blocks(200, 15, 3, rng, noise=0.03)
    sel = bernoulli_mixture_select_k(X, [1, 2, 3, 4, 5], rng=rng, n_init=6)
    assert sel["k_selected"] == 3


def test_bic_selects_one_cluster_on_noise():
    rng = np.random.default_rng(7)
    X = rng.binomial(1, 0.3, size=(300, 12))
    sel = bernoulli_mixture_select_k(X, [1, 2, 3], rng=rng, n_init=6)
    assert sel["k_selected"] == 1


def test_external_agreement_ignores_missing_and_constant_labels():
    out = external_agreement([0, 0, 1, 1],
                             {"good": ["x", "x", "y", "y"],
                              "constant": ["z"] * 4,
                              "missing": [np.nan] * 4})
    assert out["good"]["ari"] == pytest.approx(1.0)
    assert "constant" not in out and "missing" not in out


def test_influence_regime_separates_the_four_degenerate_cases():
    """n_eff alone gives the same verdict to the best and worst cases.

    Three concordant informative layers and three pure-noise layers both drive
    n_eff towards m. The `regime` field, which reads the solo agreements too,
    has to tell them apart - and a fusion of noise must set the collapse gate.
    """
    rng = np.random.default_rng(30)

    def informative(seed):
        return _blocks(200, 20, 2, np.random.default_rng(seed), noise=0.05)[0]

    cases = {
        "unstructured": [rng.integers(0, 2, size=(200, 20)) for _ in range(3)],
        "redundant": [informative(1)] * 2,
        "collapsed": [informative(2), rng.integers(0, 2, size=(200, 20))],
    }
    got = {}
    for name, mats in cases.items():
        out = layer_influence(mats, [f"L{i}" for i in range(len(mats))],
                              kernel_fn=_kernel, fuse_fn=_fuse,
                              cluster_fn=spectral_from_similarity, k=2, rng=rng)
        got[name] = out
    assert got["unstructured"]["regime"] == "unstructured"
    assert got["unstructured"]["collapse"] is True          # must set the gate
    assert got["redundant"]["regime"] == "redundant"
    assert got["collapsed"]["regime"] == "collapsed"
    assert got["collapsed"]["collapse"] is True
    # and n_eff really is uninformative on its own: noise and redundancy both
    # push it towards the number of layers
    assert got["unstructured"]["n_eff"] > 2.0
    assert got["redundant"]["n_eff"] == pytest.approx(2.0)
