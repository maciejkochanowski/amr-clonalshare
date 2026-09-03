"""Clonal share and partition attribution: properties, planted truth, gates."""
import numpy as np
import pytest
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

from amr_clonalshare.attribution import (SUPPORT_THRESHOLD, attribute_partition,
                                            clonal_share, concordance_z,
                                            layer_clonal_share)


def _two_lineage_cohort(n=400, sep=1.0, seed=0):
    """Half the isolates in lineage A, half in B, prevalence separated by
    ``sep`` (0 = lineage carries nothing, 1 = lineage determines the trait)."""
    rng = np.random.default_rng(seed)
    lin = np.repeat(["A", "B"], n // 2)
    q = 0.5 + sep * np.where(lin == "A", -0.45, 0.45)
    y = (rng.random(n) < q).astype(float)
    return y, lin


def test_clonal_share_is_near_zero_when_lineage_carries_nothing():
    y, lin = _two_lineage_cohort(sep=0.0, seed=1)
    r = clonal_share(y, lin, n_boot=100, n_perm=100, seed=1)
    assert abs(r.kappa_adj) < 0.08
    assert r.ci_low <= 0.0 <= r.ci_high
    assert r.p_value > 0.05


def test_clonal_share_is_near_one_when_lineage_determines_the_trait():
    y, lin = _two_lineage_cohort(sep=1.0, seed=2)
    r = clonal_share(y, lin, n_boot=100, n_perm=100, seed=2)
    assert r.kappa_adj > 0.75
    assert r.p_value <= 0.05


def test_debias_removes_the_penalty_the_permuted_run_measures():
    """The raw score is pushed down by the cost of estimating group means; the
    corrected one is not. With many small lineages the gap must be visible."""
    rng = np.random.default_rng(3)
    lin = rng.integers(0, 60, 300)
    y = (rng.random(300) < 0.4).astype(float)
    r = clonal_share(y, lin, n_boot=60, n_perm=120, seed=3)
    assert r.null_mean < -0.02
    assert r.kappa < r.kappa_adj
    assert abs(r.kappa_adj) < 0.12


def test_support_gate_fires_on_a_singleton_heavy_lineage_variable():
    rng = np.random.default_rng(4)
    lin = np.arange(300)                     # every isolate its own lineage
    y = (rng.random(300) < 0.4).astype(float)
    r = clonal_share(y, lin, n_boot=40, n_perm=40, seed=4)
    assert r.support == 0.0
    assert r.estimable is False

    lin2 = np.repeat(np.arange(30), 10)      # ten isolates per lineage
    r2 = clonal_share(y, lin2, n_boot=40, n_perm=40, seed=4)
    assert r2.support == 1.0
    assert r2.estimable is True
    assert SUPPORT_THRESHOLD <= 1.0


def test_missing_lineage_labels_become_their_own_level_not_dropped():
    y, lin = _two_lineage_cohort(sep=0.6, seed=5)
    lin = np.array(lin, dtype=object)
    lin[:50] = None
    r = clonal_share(y, lin, n_boot=40, n_perm=40, seed=5)
    assert r.n == len(y)                    # nothing silently removed
    assert r.n_groups == 3                  # A, B and __missing__


def test_layer_share_lies_between_its_agents():
    rng = np.random.default_rng(6)
    lin = np.repeat(np.arange(20), 20)
    clonal = (rng.random(400) < np.where(lin < 10, 0.1, 0.9)).astype(float)
    free = (rng.random(400) < 0.5).astype(float)
    k_clonal = clonal_share(clonal, lin, n_boot=60, n_perm=60, seed=6).kappa_adj
    k_free = clonal_share(free, lin, n_boot=60, n_perm=60, seed=6).kappa_adj
    k_layer = layer_clonal_share(np.column_stack([clonal, free]), lin,
                                 n_boot=60, n_perm=60, seed=6).kappa_adj
    assert k_free < k_layer < k_clonal


def _planted(kind, n=400, seed=7):
    """Traits driven either by a lineage-independent element or by the lineage."""
    rng = np.random.default_rng(seed)
    lin = np.repeat(np.arange(20), n // 20)
    driver = (rng.random(n) < 0.4) if kind == "mobile" else (lin % 2 == 0)
    block = np.column_stack([np.where(driver, rng.random(n) < 0.95,
                                      rng.random(n) < 0.05) for _ in range(6)])
    X = np.hstack([block.astype(float), (rng.random((n, 3)) < 0.3).astype(float)])
    lab = (X[:, :6].mean(1) > 0.5).astype(int)
    return X, lab, lin


def test_attribution_separates_a_mobile_archetype_from_a_clonal_one():
    Xm, labm, linm = _planted("mobile", seed=7)
    Xc, labc, linc = _planted("clonal", seed=7)
    am = attribute_partition(Xm, labm, linm, n_boot=80, repeats=6, seed=7)
    ac = attribute_partition(Xc, labc, linc, n_boot=80, repeats=6, seed=7)
    assert am.lam < 0.15, "a lineage-independent element must not read as clonal"
    assert ac.lam > 0.85, "a lineage-driven partition must read as clonal"
    assert am.lam < ac.lam


def test_shapley_shares_add_up_to_the_joint_variance_explained():
    X, lab, lin = _planted("clonal", seed=8)
    a = attribute_partition(X, lab, lin, n_boot=40, repeats=6, seed=8)
    assert a.shapley_lineage + a.shapley_partition == pytest.approx(a.r2_joint,
                                                                    abs=1e-9)


@pytest.mark.parametrize("kind", ["mobile", "clonal"])
def test_joint_term_is_never_below_either_factor_alone(kind):
    """Out-of-sample variance explained is not monotone in model nesting: an
    additive backfit that carries 20 useless lineage levels can predict worse
    than the 2-level partition alone. Letting the fitted value stand would turn
    that estimation noise into apparent shared variance, so the joint term is
    the best member of the additive family on the same folds."""
    X, lab, lin = _planted(kind, seed=9)
    a = attribute_partition(X, lab, lin, n_boot=40, repeats=6, seed=9)
    assert a.r2_joint >= a.r2_lineage - 1e-9
    assert a.r2_joint >= a.r2_partition - 1e-9
    assert 0.0 <= a.lam <= 1.0


def test_attribution_does_not_depend_on_how_finely_the_cohort_was_typed():
    """Splitting each lineage in two leaves the biology unchanged and doubles
    the levels. An undebiased split would move; this one must not."""
    X, lab, lin = _planted("clonal", n=600, seed=13)
    rng = np.random.default_rng(13)
    fine = lin * 2 + rng.integers(0, 2, len(lin))
    coarse = attribute_partition(X, lab, lin, n_boot=0, repeats=8, seed=13)
    split = attribute_partition(X, lab, fine, n_boot=0, repeats=8, seed=13)
    assert abs(coarse.lam - split.lam) < 0.20


def _weakly_linked(n, seed=10):
    """A mostly mobile element whose carriage is mildly higher in half the
    lineages. Small, fixed, non-zero lineage effect -- the real-cohort case."""
    rng = np.random.default_rng(seed)
    lin = np.repeat(np.arange(20), n // 20)
    p = np.where(lin % 2 == 0, 0.35, 0.55)
    driver = rng.random(n) < p
    block = np.column_stack([np.where(driver, rng.random(n) < 0.95,
                                      rng.random(n) < 0.05) for _ in range(6)])
    X = np.hstack([block.astype(float), (rng.random((n, 3)) < 0.3).astype(float)])
    return X, (X[:, :6].mean(1) > 0.5).astype(int), lin


def test_gate_statistic_grows_with_n_while_the_attribution_does_not():
    """The reason the module exists. One generating process with a small fixed
    lineage effect, three cohort sizes: the permutation z climbs with the number
    of pairs, the variance share does not. A gate on the z therefore refuses
    larger cohorts for being larger."""
    zs, lams = [], []
    for n in (200, 800, 3200):
        X, lab, lin = _weakly_linked(n)
        lams.append(attribute_partition(X, lab, lin, n_boot=0, repeats=4,
                                        seed=10).lam)
        zs.append(concordance_z(lab, lin, n_perm=200, seed=10)[3])
    assert zs[-1] > 2.5 * zs[0], "z must grow with cohort size"
    assert max(lams) - min(lams) < 0.20, "the attribution must be stable in n"


def test_cluster_bootstrap_interval_contains_the_point_estimate():
    """Resampling isolates *within* a drawn lineage puts duplicate rows into
    training and held-out folds at once, which leaks and lifts every replicate
    above the estimate. Lineages are therefore taken whole."""
    rng = np.random.default_rng(11)
    lin = np.repeat(np.arange(40), 8)
    y = (rng.random(320) < np.where(lin < 20, 0.25, 0.75)).astype(float)
    r = clonal_share(y, lin, n_boot=200, n_perm=100, seed=11)
    assert r.ci_low <= r.kappa_adj <= r.ci_high


def test_constant_trait_returns_nan_rather_than_a_number():
    lin = np.repeat(np.arange(10), 20)
    r = clonal_share(np.ones(200), lin, n_boot=20, n_perm=20, seed=12)
    assert np.isnan(r.kappa) or np.isnan(r.kappa_adj)


def test_a_non_finite_trait_value_is_dropped_rather_than_pinning_the_p_value_to_its_floor():
    """A missing well made every ``null >= kappa`` comparison false, so the
    permutation test counted no exceedances and reported 1 / (n_perm + 1) --
    the most significant value obtainable -- beside a NaN kappa_adj and
    estimable true, on a cohort where the lineage carries almost nothing."""
    rng = np.random.default_rng(21)
    lin = np.array([f"ST{i % 6}" for i in range(60)])
    y = (rng.random(60) < 0.4).astype(float)
    y[7] = np.nan
    r = clonal_share(y, lin, n_boot=50, n_perm=200, repeats=5, seed=1)
    assert r.n_dropped_non_finite == 1
    assert r.n == 59
    assert r.p_value > 1.0 / (200 + 1)
    assert np.isfinite(r.kappa_adj)
    assert r.as_dict()["n_dropped_non_finite"] == 1


def test_a_length_mismatch_is_refused_rather_than_transposed_into_a_cohort_of_one():
    """Transposing on any mismatch turned 60 traits against 59 labels into
    n = 1 and reported a number for it. A transpose is the answer only when
    the transpose is what resolves the mismatch."""
    rng = np.random.default_rng(22)
    lin = np.array([f"ST{i % 6}" for i in range(60)])
    y = (rng.random(60) < 0.4).astype(float)
    with pytest.raises(ValueError, match="60 rows and lineage has 59"):
        clonal_share(y, lin[:59], n_boot=20, n_perm=20, repeats=2, seed=1)
    with pytest.raises(ValueError, match="60 rows and labels has 59"):
        attribute_partition(np.column_stack([y, y]),
                            np.array(["A"] * 30 + ["B"] * 29), lin[:59],
                            n_boot=20, n_perm=3, repeats=2, seed=1)
    # a genuine (p, n) block still has to be recognised and turned round
    assert layer_clonal_share(np.column_stack([y, y]).T, lin, n_boot=20,
                              n_perm=20, repeats=2, seed=1).n == 60


def test_fewer_than_two_folds_is_refused_rather_than_read_as_no_lineage_signal():
    """One fold holds nothing out, so every score was taken on an empty
    held-out set and the run returned kappa 0, a zero-width interval, p = 1 and
    estimable true: a clean acquittal of the lineage that nothing was tested
    for. Zero folds reached a modulo by zero first."""
    rng = np.random.default_rng(23)
    lin = np.array([f"ST{i % 6}" for i in range(60)])
    y = (rng.random(60) < 0.4).astype(float)
    for folds in (1, 0):
        with pytest.raises(ValueError, match="at least two folds"):
            clonal_share(y, lin, folds=folds, n_boot=20, n_perm=20, repeats=2,
                         seed=1)
    with pytest.raises(ValueError, match="at least two folds"):
        attribute_partition(np.column_stack([y, y]),
                            np.array(["A"] * 30 + ["B"] * 30), lin, folds=1,
                            n_boot=10, n_perm=3, repeats=2, seed=1)


# --------------------------------------------------- fold-internal refitting
def _mixed_cohort(n=300, seed=12):
    """A cohort where the lineage and the partition each carry part of six
    binary traits and neither carries all of it, so every term of the split
    and both ends of the interval sit away from a boundary. Recorded numbers
    taken on a boundary would pass the test below without testing anything."""
    rng = np.random.default_rng(seed)
    lin = rng.integers(0, 6, n)
    lab = np.where(rng.random(n) < 0.7, lin % 3, rng.integers(0, 3, n))
    q = 0.15 + 0.20 * (lin % 2) + 0.35 * (lab == 1)
    X = (rng.random((n, 6)) < q[:, None]).astype(float)
    return X, lab, lin


def _separable_cohort(n=240, seed=5):
    """Three archetypes far enough apart that dropping a fifth of the isolates
    cannot move them, and a lineage drawn independently of all three."""
    rng = np.random.default_rng(seed)
    lab = np.repeat([0, 1, 2], n // 3)
    centres = np.array([[1, 1, 1, 0, 0, 0],
                        [0, 0, 0, 1, 1, 1],
                        [1, 0, 1, 0, 1, 0]], dtype=float)
    X = np.abs(centres[lab] - (rng.random((n, 6)) < 0.02))
    return X, lab, rng.integers(0, 8, n)


def _ward_refit(X, k):
    """A refit callable of the shape attribute_partition asks for: rebuild the
    partition from the training rows and put each held-out row on the nearest
    training centroid. Ward rather than the pipeline's own fusion, because a
    test of the plumbing should not carry a dependency on the clustering."""
    def refit(train, held_out, labels=None):
        if labels is None:
            d = pdist(X[train], "euclidean")
            lab_train = fcluster(linkage(d, "ward"), k, "maxclust") - 1
        else:
            lab_train = np.asarray(labels)[train]
        present = np.unique(lab_train)
        centres = np.array([X[train][lab_train == c].mean(axis=0)
                            for c in present])
        far = ((X[held_out][:, None, :] - centres[None, :, :]) ** 2).sum(axis=2)
        out = np.empty(len(train) + len(held_out), dtype=int)
        out[train] = lab_train
        out[held_out] = present[far.argmin(axis=1)]
        return out
    return refit


def test_refit_none_reproduces_the_recorded_numbers_exactly():
    """The gate on the whole refit change. These are the values this function
    returned before a refit path existed, recorded to full precision, and the
    comparison is exact rather than approximate on purpose: the refit work
    reorganises the scorers, and anything that perturbs a caller who did not
    ask for a refit is a defect however small the perturbation looks."""
    X, lab, lin = _mixed_cohort()
    a = attribute_partition(X, lab, lin, folds=5, repeats=6, n_boot=40,
                            n_perm=30, seed=17)
    assert a.r2_lineage == 0.10119220358989915
    assert a.r2_partition == 0.13611294066043805
    assert a.r2_joint == 0.1693841803620396
    assert a.shared == 0.06792096388829762
    assert a.lam == 0.4989681503517153
    assert a.shapley_lineage == 0.06723172164575035
    assert a.shapley_partition == 0.10215245871628924
    assert a.ci_low == 0.09714590544629859
    assert a.ci_high == 0.7039208339917482
    assert a.lam_clipped is False
    assert a.partition_refit is False
    assert np.isnan(a.r2_partition_assignment_control)
    assert np.isnan(a.lam_assignment_control)


def test_a_refit_that_returns_the_fixed_partition_agrees_with_no_refit():
    """The refit path and the fixed path are the same estimator whenever the
    callable hands back the partition it was going to be given anyway. They do
    not draw their folds in the same order, so the agreement is to Monte Carlo
    error rather than exact, but a plumbing mistake in the per-fold codes would
    not land within it."""
    X, lab, lin = _mixed_cohort()
    fixed = attribute_partition(X, lab, lin, folds=5, repeats=12, n_boot=0,
                                n_perm=40, seed=31)
    refit = attribute_partition(X, lab, lin, folds=5, repeats=12, n_boot=0,
                                n_perm=40, seed=31,
                                refit=lambda train, held_out, labels=None: lab)
    assert refit.partition_refit is True
    assert abs(refit.r2_partition - fixed.r2_partition) < 0.02
    assert abs(refit.r2_lineage - fixed.r2_lineage) < 0.02
    assert abs(refit.lam - fixed.lam) < 0.05
    assert np.isfinite(refit.lam_cv_sd)
    assert np.isfinite(refit.r2_partition_assignment_control)


def test_a_separable_partition_is_unchanged_by_refitting_inside_the_fold():
    """Where the archetypes are far apart the training rows rebuild the same
    partition and place the held-out rows back where they were, so refitting
    has nothing to correct. The share stays near zero because the lineage was
    drawn independently. A refit that moved this number would be measuring its
    own instability rather than leakage."""
    X, lab, lin = _separable_cohort()
    fixed = attribute_partition(X, lab, lin, folds=5, repeats=8, n_boot=0,
                                n_perm=30, seed=41)
    refit = attribute_partition(X, lab, lin, folds=5, repeats=8, n_boot=0,
                                n_perm=30, seed=41, refit=_ward_refit(X, 3))
    assert fixed.lam < 0.15
    assert refit.lam < 0.15
    assert abs(refit.r2_partition - fixed.r2_partition) < 0.05


def test_the_lineage_bootstrap_is_refused_together_with_a_refit():
    """A lineage drawn twice puts the same isolates in the cohort twice, and a
    partition rebuilt on that cohort trains on rows that are also held out. The
    interval that would come back is narrower than the truth for a reason the
    caller cannot see, so the combination is refused rather than served."""
    X, lab, lin = _mixed_cohort()
    with pytest.raises(ValueError, match="n_boot=10 together with refit"):
        attribute_partition(X, lab, lin, folds=5, repeats=4, n_boot=10,
                            n_perm=10, seed=17,
                            refit=lambda train, held_out, labels=None: lab)
