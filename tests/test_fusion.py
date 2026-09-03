"""SNF kernel fidelity, order invariance, single-layer handling."""
import numpy as np
import pytest

from amr_clonalshare.core import spectral_from_similarity
from amr_clonalshare.fusion import (_cross_diffuse, _renormalise,
                                       _rowstochastic, default_K,
                                       knn_tie_diagnostics, snf_fuse,
                                       snf_kernel, sparsify)
from amr_clonalshare.stats import binary_distance_matrix
from sklearn.metrics import adjusted_rand_score


def _reference_kernel(D, mu, K):
    """Wang et al. 2014 Eq. 1, written out independently of the implementation."""
    n = D.shape[0]
    mean_d = np.empty(n)
    for i in range(n):
        others = np.delete(D[i], i)
        mean_d[i] = np.sort(others)[:K].mean()
    eps = (mean_d[:, None] + mean_d[None, :] + D) / 3.0
    W = np.exp(-(D ** 2) / (mu * np.maximum(eps, 1e-12)))
    np.fill_diagonal(W, 0.0)
    return (W + W.T) / 2


def test_kernel_uses_knn_mean_not_global_mean():
    rng = np.random.default_rng(0)
    D = binary_distance_matrix(rng.integers(0, 2, size=(40, 12)))
    got = snf_kernel(D, mu=0.5, K=8)
    assert np.allclose(got, _reference_kernel(D, 0.5, 8), atol=1e-12)
    # the global-mean variant is a different matrix
    mean_all = D.mean(axis=1)
    eps = (mean_all[:, None] + mean_all[None, :] + D) / 3.0
    old = np.exp(-(D ** 2) / (0.5 * np.maximum(eps, 1e-12)))
    np.fill_diagonal(old, 0.0)
    assert not np.allclose(got, (old + old.T) / 2, atol=1e-6)


def test_single_layer_fusion_does_not_raise():
    rng = np.random.default_rng(1)
    W = snf_kernel(binary_distance_matrix(rng.integers(0, 2, size=(20, 6))), K=5)
    out = snf_fuse([W], K=5, T=5)
    assert out.shape == (20, 20)
    assert np.allclose(out, out.T)


def test_inclusive_sparsification_is_row_order_invariant():
    rng = np.random.default_rng(2)
    X = rng.integers(0, 2, size=(60, 8))
    W = snf_kernel(binary_distance_matrix(X), K=10)
    perm = rng.permutation(60)
    inv = np.argsort(perm)
    Wp = snf_kernel(binary_distance_matrix(X[perm]), K=10)
    S = sparsify(W, 10, tie_policy="inclusive")
    Sp = sparsify(Wp, 10, tie_policy="inclusive")[np.ix_(inv, inv)]
    assert np.allclose(S, Sp, atol=1e-12)


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_pipeline_is_permutation_equivariant(seed):
    """Relabelling the samples must not change the partition.

    This is the regression test for the defect that makes a partition
    a function of the input file's row order.
    """
    rng = np.random.default_rng(seed)
    n = 90
    truth = np.repeat([0, 1, 2], n // 3)
    X1 = np.zeros((n, 9), dtype=int)
    X2 = np.zeros((n, 12), dtype=int)
    for c in range(3):
        X1[truth == c, c * 3:(c + 1) * 3] = 1
        X2[truth == c, c * 4:(c + 1) * 4] = 1
    X1 = np.clip(X1 + rng.binomial(1, 0.05, X1.shape), 0, 1)
    X2 = np.clip(X2 + rng.binomial(1, 0.05, X2.shape), 0, 1)

    def labels_of(mats):
        Ws = [snf_kernel(binary_distance_matrix(M), K=12) for M in mats]
        return spectral_from_similarity(snf_fuse(Ws, K=12, T=10), 3)

    base = labels_of([X1, X2])
    perm = rng.permutation(n)
    inv = np.argsort(perm)
    permuted = labels_of([X1[perm], X2[perm]])[inv]
    assert adjusted_rand_score(base, permuted) == pytest.approx(1.0)


def test_tie_diagnostics_detect_a_degenerate_graph():
    W = np.ones((30, 30))          # every affinity tied
    np.fill_diagonal(W, 0.0)
    d = knn_tie_diagnostics(W, K=5)
    assert d["tie_inflation"] > 5
    assert d["frac_rows_tied"] == 1.0


def test_default_K_matches_the_documented_heuristic():
    assert default_K(50) == 10
    assert default_K(1500) == 30
    assert default_K(3) == 2


def test_add_identity_update_reproduces_snfpy():
    """`update="add_identity"` is the 2014 MATLAB/SNFtool<=2.2/snfpy update.

    Verified against snfpy itself rather than against a transcription of it.
    Skipped when snfpy is not installed; `benchmarks/snf_update_benchmark.py`
    runs the same check and reports the number.
    """
    snfc = pytest.importorskip("snf.compute")
    orig = snfc.check_array

    def shim(a, **kw):                      # sklearn >= 1.6 renamed the kwarg
        if "force_all_finite" in kw:
            try:
                return orig(a, **kw)
            except TypeError:
                kw["ensure_all_finite"] = kw.pop("force_all_finite")
        return orig(a, **kw)

    snfc.check_array = shim
    try:
        n, m, K, T = 40, 3, 10, 15
        rng = np.random.default_rng(0)
        affs = []
        for _ in range(m):
            A = rng.random((n, n))
            A = (A + A.T) / 2.0
            np.fill_diagonal(A, 0.0)
            affs.append(A)
        theirs = snfc.snf(*[a.copy() for a in affs], K=K, t=T, alpha=1.0)

        idx = np.arange(n)
        P = [_rowstochastic(a) for a in affs]
        P = [(p + p.T) / 2.0 for p in P]
        S = [snfc._find_dominate_set(p, K) for p in P]
        P = _cross_diffuse(P, S, T=T, alpha=1.0, update="add_identity", idx=idx)
        ours = _rowstochastic(sum(P) / m)
        ours = (ours + ours.T + np.eye(n)) / 2.0
        assert np.abs(ours - theirs).max() < 1e-10
    finally:
        snfc.check_array = orig


def test_renormalise_matches_snftool_normalize_semantics():
    """SNFtool >= 2.2.1 `normalize`: off-diagonal mass 1-alpha, diagonal alpha."""
    rng = np.random.default_rng(1)
    W = rng.random((12, 12))
    W = (W + W.T) / 2.0
    np.fill_diagonal(W, 0.0)
    P = _renormalise(W, 0.5, np.arange(12))
    assert np.allclose(np.diag(P), 0.5)
    off = P - np.diag(np.diag(P))
    assert np.allclose(off.sum(axis=1), 0.5)


def test_the_two_updates_agree_on_planted_structure():
    """The choice of update is documented as immaterial; hold it to that."""
    rng = np.random.default_rng(7)
    n = 90
    truth = np.repeat(np.arange(3), n // 3)
    mats = []
    for p in (12, 20):
        base = rng.binomial(1, 0.5, size=(3, p))
        X = base[truth]
        X = np.abs(X - rng.binomial(1, 0.05, size=X.shape))
        mats.append(X)
    Ws = [snf_kernel(binary_distance_matrix(M), K=12) for M in mats]
    a = spectral_from_similarity(snf_fuse(Ws, K=12, T=15, alpha=0.5,
                                          update="renormalise"), 3)
    b = spectral_from_similarity(snf_fuse(Ws, K=12, T=15, alpha=1.0,
                                          update="add_identity"), 3)
    assert adjusted_rand_score(truth, a) == pytest.approx(1.0)
    assert adjusted_rand_score(a, b) == pytest.approx(1.0)
