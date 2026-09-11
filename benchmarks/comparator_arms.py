#!/usr/bin/env python3
"""Two comparator arms for the estimator benchmark, registered before running.

    python comparator_arms.py --folds-identity --out <dir>
    python comparator_arms.py --anchor --out <dir>
    python comparator_arms.py --synthetic s1 --cell <i> --out <dir>
    python comparator_arms.py --synthetic s2 --cell <i> --out <dir>
    python comparator_arms.py --real --raw <raw> --atlas <vet_atlas.json> --out <dir>

WHY THIS EXISTS. Two substitutes will be named by a reviewer. The first is
population-structure-only classification accuracy: Moradigaravand and
colleagues report a mean accuracy of 0.79 across eleven compounds for a
predictor that sees population structure and nothing else, and the objection is
that the clonal share re-expresses that number. The second is mixed-model
heritability, which the benchmark already runs under the name reml_lmm without
saying what it is. Neither objection is answered by argument, so both are
turned into measurements here.

ARM B, WHAT IT IS AND WHAT IT IS NOT. The arm predicts the binary phenotype
from lineage alone, out of fold, by the training-fold majority class inside the
isolate's lineage, with the overall training majority as the fallback for a
lineage the training fold never saw. That is the most favourable honest form of
the population-structure-only predictor on a cohort whose only structural
variable is a cluster label. It is scored three ways, because a single accuracy
is not interpretable: plain accuracy, balanced accuracy, and the majority-class
baseline max(p, 1 - p), which is what predicting the majority for every isolate
already scores with no lineage information whatsoever. At prevalence 0.90 that
baseline is 0.90, which is above the accuracy the objection is built on, and
the point of reporting all three together is that the reader can see how much
of an accuracy is prevalence and how much is population structure.

WHY THE FOLDS MUST BE THE SHIPPED ONES. A correlation between an accuracy and a
variance share means nothing if the two were computed on different partitions
of the cohort. The arm therefore imports the package's own _folds and rebuilds
the generator state that layer_clonal_share creates, so the fold assignments
are the same objects and not merely the same recipe. That identity is not
assumed: --folds-identity reconstructs the estimator's own kappa from the fold
vectors this module draws and refuses to continue unless it matches the
shipped kappa to machine precision.

ARM A, WHAT IS VERIFIED AND WHAT IS REFUSED. The shipped reml_lmm profiles the
restricted likelihood of the one-way random-effects model on the variance
ratio. That is the mixed-model heritability estimator applied to a discrete
lineage factor. The claim is verified against an anchor written from the model
statement in the general marginal form, with a dense covariance matrix, a
generalised-least-squares mean and explicit log determinants, sharing no line
of algebra with the profiled scalar form. What is refused is a genomic
relatedness arm. Bugwas and pyseer build kinship from genome-wide variants, and
these cohorts carry a cluster accession and a susceptibility call and no
variant matrix, so such an arm could only be faked. It is not run, and what
would be required is written in the report instead.

WHAT WOULD MEAN THIS MEASURED THE WRONG THING. The prereg names it: if every
cell of a grid sits near prevalence one half, the majority baseline is a
constant and accuracy minus baseline is accuracy shifted, so every correlation
agrees by construction. leverage_check.py computes that verdict before any
number here is produced, and the shipped binary grid is reported as the narrow
case it is rather than as a confirmation.
"""
from __future__ import annotations

import argparse
import json
import platform
import time
from datetime import date, timezone, datetime
from pathlib import Path

import numpy as np
import scipy
from scipy import linalg, optimize

from amr_clonalshare.attribution import (_codes, _folds, _rng, _skill,
                                            clonal_share)
from estimator_benchmark import (ALPHA, SEED, cohort, design, in_sample_r2,
                                 reml_component)
from leverage_check import SWEEP_PREVALENCES_EXTENDED, sweep_design

#: The settings the shipped benchmark passes to clonal_share, so that the
#: comparator is scored against the published estimate and not a variant of it.
BENCH_FOLDS, BENCH_REPEATS, BENCH_BOOT, BENCH_PERM = 5, 10, 200, 100
#: The settings benchmarks/vet_atlas.py passed for the published atlas.
VET_SEED, VET_FOLDS, VET_REPEATS = 42, 5, 10
PREREG_SHA256 = ("4d54749cf39ed17db74ad997fec31dd617"
                 "38a5500c77961b8809627780438df9")


# ------------------------------------------------------------ the fold scheme
def shipped_folds(n: int, seed, folds: int, repeats: int) -> list[np.ndarray]:
    """The fold assignments layer_clonal_share itself would draw.

    layer_clonal_share builds its generator with _rng(seed) and _repeated_skill
    is the first thing to draw from it, taking one _folds(n, folds, rng) per
    repeat. Rebuilding the generator the same way and drawing the same calls in
    the same order therefore reproduces the estimator's own partitions rather
    than merely imitating the recipe. Nothing here re-implements the scheme:
    _folds and _rng are the package's own.
    """
    rng = _rng(seed)
    return [_folds(n, folds, rng) for _ in range(repeats)]


def folds_identity(y, lineage, seed, folds: int, repeats: int) -> dict:
    """Prove the partitions are the estimator's, rather than assume it.

    The raw kappa the estimator reports is the mean over repeats of the skill
    on each fold vector. Recomputing that mean from the fold vectors this
    module draws and comparing it with the shipped kappa is a check that fails
    loudly if the generator state has moved, which is the failure the prereg
    declares void: two numbers computed on different partitions of a cohort are
    not comparable and no correlation between them means anything.
    """
    X = np.asarray(y, dtype=float).reshape(-1, 1)
    code = _codes(np.asarray(lineage, dtype=object).ravel())
    vectors = shipped_folds(X.shape[0], seed, folds, repeats)
    rebuilt = float(np.nanmean([_skill(X, code, f) for f in vectors]))
    shipped = clonal_share(y, lineage, folds=folds, repeats=repeats,
                           n_boot=20, n_perm=5, seed=seed).kappa
    return {"kappa_rebuilt_from_our_folds": rebuilt,
            "kappa_from_shipped_estimator": float(shipped),
            "absolute_difference": abs(rebuilt - float(shipped)),
            "identical": bool(abs(rebuilt - float(shipped)) < 1e-12)}


# ------------------------------------------------- arm B, the accuracy arm
def _majority(values: np.ndarray, fallback: int) -> int:
    """Majority class of a training block, ties resolved without randomness.

    A tie inside a lineage carries no information about which class that
    lineage favours, so it defers to the fallback the caller supplies, which is
    the overall training majority. Breaking such a tie at random would put
    Monte-Carlo noise into a quantity the reader will compare across cells.
    """
    if values.size == 0:
        return fallback
    m = float(values.mean())
    if m > 0.5:
        return 1
    if m < 0.5:
        return 0
    return fallback


def structure_only_accuracy(y, lineage, *, seed, folds: int,
                            repeats: int) -> dict:
    """Out-of-fold accuracy of the lineage-majority predictor, and its baseline.

    Three numbers are produced for every repeat and then averaged. The first is
    plain accuracy, which is the quantity the population-structure-only
    literature reports. The second is balanced accuracy, the mean of
    sensitivity and specificity, which is the same prediction scored so that
    the rare class counts as much as the common one. The third is the accuracy
    of the same fold scheme with the lineage term removed, so the baseline pays
    the identical cross-validation cost and the difference between the two is
    attributable to lineage rather than to the design.

    The marginal baseline max(p, 1 - p) is reported beside them because it is
    the number a reader can compute from the prevalence alone, and because the
    whole force of the comparison is that at prevalence 0.90 it is already 0.90.
    """
    y = np.asarray(y, dtype=float).ravel()
    code = _codes(np.asarray(lineage, dtype=object).ravel())
    n = y.size
    vectors = shipped_folds(n, seed, folds, repeats)
    acc, bal, base_oof = [], [], []
    for fold in vectors:
        pred = np.empty(n, dtype=float)
        base = np.empty(n, dtype=float)
        for f in np.unique(fold):
            tr, te = fold != f, fold == f
            if not tr.any() or not te.any():
                pred[te], base[te] = np.nan, np.nan
                continue
            grand = _majority(y[tr], 1)
            table = {}
            for c in np.unique(code[tr]):
                table[int(c)] = _majority(y[tr][code[tr] == c], grand)
            pred[te] = [table.get(int(c), grand) for c in code[te]]
            base[te] = grand
        scored = np.isfinite(pred)
        acc.append(float((pred[scored] == y[scored]).mean()))
        base_oof.append(float((base[scored] == y[scored]).mean()))
        pos, neg = scored & (y == 1), scored & (y == 0)
        bal.append(float(0.5 * ((pred[pos] == 1).mean() + (pred[neg] == 0).mean()))
                   if pos.any() and neg.any() else float("nan"))
    prevalence = float(y.mean())
    baseline = float(max(prevalence, 1.0 - prevalence))
    out = {"n": int(n), "n_lineages": int(code.max()) + 1,
           "prevalence": prevalence, "baseline_marginal": baseline,
           "repeats": int(repeats), "folds": int(folds)}
    for name, series in (("accuracy", acc), ("balanced_accuracy", bal),
                         ("baseline_out_of_fold", base_oof)):
        v = np.asarray(series, dtype=float)
        out[name] = float(np.nanmean(v))
        out[f"{name}_mc_se"] = (float(np.nanstd(v, ddof=1) / np.sqrt(v.size))
                                if v.size > 1 else float("nan"))
    out["lift"] = out["accuracy"] - baseline
    out["lift_over_out_of_fold_baseline"] = (out["accuracy"]
                                             - out["baseline_out_of_fold"])
    out["binomial_se_of_accuracy"] = float(
        np.sqrt(max(out["accuracy"] * (1.0 - out["accuracy"]), 0.0) / n))
    return out


# ---------------------------------------- arm A, the mixed-model identification
def anchor_reml(y, lineage) -> dict:
    """REML for the one-way random-effects model, derived from the model.

    Written from the model statement rather than from the shipped code. The
    marginal law of the data is y normal with mean X mu and covariance
    sigma2 (I + r Z Z transpose), where Z is the lineage incidence matrix and
    r is the ratio of the between-lineage variance to the within-lineage one.
    The restricted log likelihood is formed here in that general form, with a
    dense covariance, a Cholesky factorisation, a generalised-least-squares
    mean, an explicit log determinant of the covariance and an explicit log
    determinant of the information for the fixed effect. The scale is profiled
    out analytically because it enters the restricted likelihood only through
    the residual quadratic form.

    None of the algebra that makes the shipped estimator fast appears here: no
    eigenvalues written down from the block structure, no Woodbury identity, no
    lineage-mean weights. That is the point. Two implementations that share
    those steps would agree even if the steps were wrong, and the question this
    answers is whether the shipped arm is the mixed model it is about to be
    called in the manuscript.
    """
    y = np.asarray(y, dtype=float).ravel()
    _labels, code = np.unique(np.asarray(lineage, dtype=object).ravel(),
                              return_inverse=True)
    n, groups = y.size, int(code.max()) + 1
    Z = np.zeros((n, groups))
    Z[np.arange(n), code] = 1.0
    K = Z @ Z.T
    X = np.ones((n, 1))
    identity = np.eye(n)
    dfe = float(n - 1)

    def negative_restricted(log_ratio: float) -> float:
        ratio = float(np.exp(log_ratio))
        try:
            factor = linalg.cho_factor(identity + ratio * K, lower=True)
        except linalg.LinAlgError:
            return float("inf")
        solved_x = linalg.cho_solve(factor, X)
        solved_y = linalg.cho_solve(factor, y)
        information = float((X.T @ solved_x).item())
        mu = float((X.T @ solved_y).item()) / information
        residual = y - mu
        quadratic = float(residual @ linalg.cho_solve(factor, residual))
        if quadratic <= 0.0:
            return float("inf")
        log_det = 2.0 * float(np.sum(np.log(np.diag(factor[0]))))
        return 0.5 * (dfe * np.log(quadratic / dfe) + log_det
                      + np.log(information) + dfe)

    grid = np.linspace(-12.0, 6.0, 41)
    values = np.array([negative_restricted(v) for v in grid])
    start = float(grid[int(np.argmin(values))])
    fit = optimize.minimize_scalar(negative_restricted,
                                   bounds=(start - 0.6, start + 0.6),
                                   method="bounded")
    ratio = float(np.exp(float(fit.x)))
    return {"ratio": ratio, "point": ratio / (1.0 + ratio),
            "negative_restricted_likelihood": float(fit.fun),
            "n": int(n), "n_groups": int(groups)}


def run_anchor() -> dict:
    """Compare the anchor with the shipped reml_lmm across the Gaussian grid.

    One replicate is drawn from each of the thirty Gaussian cells of the
    shipped design, using that cell's own generator and its own seed, so the
    verification spans every lineage count, every cohort size and every true
    share the benchmark reports on. The registered tolerance is an absolute
    difference of 0.01 on the share scale. A larger difference does not mean
    the anchor is wrong, it means the identification of the shipped arm as the
    mixed model is not established, and the manuscript sentence is then not
    licensed.
    """
    rows = []
    for cell in [c for c in design() if not c["binary"]]:
        rng = np.random.default_rng(SEED + 977 * cell["index"])
        y, lineage, truth = cohort(rng, cell["n_groups"], cell["group_size"],
                                   cell["share"], binary=False,
                                   unbalance=cell["unbalance"])
        t0 = time.perf_counter()
        anchor = anchor_reml(y, lineage)
        seconds = time.perf_counter() - t0
        shipped = reml_component(y, lineage)
        rows.append({"index": cell["index"], "n_groups": cell["n_groups"],
                     "group_size": cell["group_size"], "share": cell["share"],
                     "unbalance": cell["unbalance"],
                     "truth_superpopulation": truth["superpopulation"],
                     "truth_realised": truth["realised"],
                     "anchor_point": anchor["point"],
                     "shipped_reml_point": shipped["point"],
                     "absolute_difference": abs(anchor["point"]
                                                - shipped["point"]),
                     "in_sample_r2_point": in_sample_r2(y, lineage)["point"],
                     "anchor_seconds": seconds})
    differences = np.array([r["absolute_difference"] for r in rows])
    return {"tolerance": 0.01, "cohorts": rows,
            "max_absolute_difference": float(differences.max()),
            "mean_absolute_difference": float(differences.mean()),
            "agrees_within_tolerance": bool(differences.max() <= 0.01)}


# ------------------------------------------------------------ synthetic grids
def _mean_and_mc_se(values) -> tuple[float, float]:
    """A replicate mean and the Monte-Carlo error actually achieved on it."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan"), float("nan")
    return (float(v.mean()),
            float(v.std(ddof=1) / np.sqrt(v.size)) if v.size > 1
            else float("nan"))


def run_synthetic_cell(cell: dict) -> dict:
    """One synthetic cell: the accuracy arm and the clonal share, paired.

    Both quantities are computed on the same replicate, from the same fold
    vectors, so a rank correlation between their cell means is a statement
    about the two estimators rather than about two different draws.
    """
    rng = np.random.default_rng(SEED + 977 * cell["index"])
    kwargs = {"binary": True, "unbalance": cell["unbalance"]}
    if "prevalence" in cell:
        kwargs["prevalence"] = cell["prevalence"]
    series: dict[str, list] = {k: [] for k in
                               ("accuracy", "balanced_accuracy",
                                "baseline_marginal", "baseline_out_of_fold",
                                "lift", "lift_over_out_of_fold_baseline",
                                "prevalence", "kappa_adj", "kappa_raw",
                                "truth", "n_lineages_used")}
    identity = None
    for replicate in range(cell["replicates"]):
        y, lineage, truth = cohort(rng, cell["n_groups"], cell["group_size"],
                                   cell["share"], **kwargs)
        if identity is None:
            identity = folds_identity(y, lineage, SEED, BENCH_FOLDS,
                                      BENCH_REPEATS)
            if not identity["identical"]:
                raise RuntimeError(
                    "the fold vectors this module draws are not the ones the "
                    "shipped estimator uses, so the two arms are not "
                    f"comparable: {identity}")
        arm = structure_only_accuracy(y, lineage, seed=SEED,
                                      folds=BENCH_FOLDS, repeats=BENCH_REPEATS)
        share = clonal_share(y, lineage, folds=BENCH_FOLDS,
                             repeats=BENCH_REPEATS, n_boot=BENCH_BOOT,
                             n_perm=BENCH_PERM, seed=SEED).as_dict()
        for key in ("accuracy", "balanced_accuracy", "baseline_marginal",
                    "baseline_out_of_fold", "lift",
                    "lift_over_out_of_fold_baseline", "prevalence"):
            series[key].append(arm[key])
        series["kappa_adj"].append(share["kappa_adj"])
        series["kappa_raw"].append(share["kappa"])
        series["truth"].append(truth["superpopulation"])
        series["n_lineages_used"].append(share["n_groups"])
    out = dict(cell)
    out["folds_identity"] = identity
    for key, values in series.items():
        mean, mc_se = _mean_and_mc_se(values)
        out[key] = mean
        out[f"{key}_mc_se"] = mc_se
    # The loop above writes the realised marginal prevalence into the field the
    # design used for the target prevalence, and the two are not the same
    # number: a lineage effect on the probit scale pulls the marginal towards
    # one half. The target is what indexes the sweep, so it is kept under its
    # own name rather than left to be recovered from the cell index.
    out["prevalence_target"] = float(cell.get("prevalence", 0.25))
    out["prevalence_realised"] = out["prevalence"]
    return out


def synthetic_grid(name: str) -> list[dict]:
    """The registered grids. s3 is the addendum extension of s2.

    s2 and s3 share the generator, the geometry and the seed rule and differ
    only in how finely the prevalence axis is sampled, so a cell that appears
    in both is the same cohort drawn from the same stream only when its index
    coincides, which it does not in general. Both are reported.
    """
    if name == "s1":
        return [c for c in design() if c["binary"]]
    if name == "s2":
        return sweep_design()
    if name == "s3":
        return sweep_design(SWEEP_PREVALENCES_EXTENDED)
    raise ValueError(f"grid={name}; the registered grids are s1, s2 and s3")


# ------------------------------------------------------ the veterinary cells
def run_real(raw: Path, atlas: Path) -> dict:
    """Arm B on the published veterinary cells, against the published share.

    The cohort is not rebuilt from a rewritten loader. benchmarks/vet_atlas.py
    is imported and its own load, select and parse_ast are called, so the host
    cut, the matrix cut, the treatment of an undetermined susceptibility call
    and the lineage variable are the ones the published atlas used. The clonal
    share is read back from the artefact and never recomputed, so what the
    accuracy is compared against is the published number.

    The isolate count of every cell is checked against the artefact before the
    cell is kept. A mismatch means the arm is describing a different cohort
    than the share it is paired with, which the prereg declares void for that
    cell, and the mismatched cells are reported rather than dropped silently.
    """
    import vet_atlas

    payload = json.loads(atlas.read_text())
    wanted: dict[str, list[dict]] = {}
    for record in payload["cohorts"]:
        wanted.setdefault(record["organism"], []).append(record)
    rows, mismatched, identity_checked = [], [], None
    for organism in sorted(wanted):
        frame = vet_atlas.load(raw, organism)
        if frame is None:
            mismatched.append({"organism": organism,
                               "reason": "source tables absent"})
            continue
        for record in wanted[organism]:
            cell = {"organism": organism, "host_group": record["host_group"],
                    "matrix": record["matrix"]}
            sub = vet_atlas.select(frame, cell)
            if len(sub) != record["n_isolates"]:
                mismatched.append(dict(cell, rebuilt=int(len(sub)),
                                       published=int(record["n_isolates"]),
                                       reason="cohort size differs"))
                continue
            parsed = [vet_atlas.parse_ast(s) for s in sub["AST_phenotypes"]]
            lineage = sub["PDS_acc"].tolist()
            for agent, published in record["agents"].items():
                if "kappa" not in published:
                    continue
                y = np.array([r.get(agent, np.nan) for r in parsed],
                             dtype=float)
                ok = np.isfinite(y)
                if int(ok.sum()) != int(published["n"]):
                    mismatched.append(dict(cell, agent=agent,
                                           rebuilt=int(ok.sum()),
                                           published=int(published["n"]),
                                           reason="agent isolate count differs"))
                    continue
                yy = y[ok]
                ll = [lineage[i] for i in np.where(ok)[0]]
                if identity_checked is None:
                    identity_checked = folds_identity(yy, ll, VET_SEED,
                                                      VET_FOLDS, VET_REPEATS)
                    if not identity_checked["identical"]:
                        raise RuntimeError(
                            "the fold vectors this module draws are not the "
                            "ones the published atlas used, so the arms are "
                            f"not comparable: {identity_checked}")
                arm = structure_only_accuracy(yy, ll, seed=VET_SEED,
                                              folds=VET_FOLDS,
                                              repeats=VET_REPEATS)
                rows.append(dict(cell, agent=agent, **arm,
                                 kappa=float(published["kappa"]),
                                 kappa_ci=[published["ci_low"],
                                           published["ci_high"]],
                                 kappa_permuted=float(
                                     published["kappa_permuted"]),
                                 support=float(published["support"]),
                                 clonal_estimable=bool(published["estimable"]),
                                 published_prevalence=float(
                                     published["prevalence"])))
            print(f"{organism} {record['host_group']} {record['matrix']}: "
                  f"{len(rows)} cells so far", flush=True)
    return {"cells": rows, "mismatched": mismatched,
            "folds_identity": identity_checked}


def provenance() -> dict:
    return {"generated": date.today().isoformat(),
            "generated_utc": datetime.now(timezone.utc).strftime(
                "%Y-%m-%dT%H:%M:%SZ"),
            "prereg_sha256": PREREG_SHA256,
            "seed": SEED, "vet_seed": VET_SEED, "alpha": ALPHA,
            "estimator_settings": {"folds": BENCH_FOLDS,
                                   "repeats": BENCH_REPEATS,
                                   "bootstrap": BENCH_BOOT,
                                   "permutations": BENCH_PERM},
            "python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--synthetic", choices=("s1", "s2", "s3"))
    ap.add_argument("--cell", type=int)
    ap.add_argument("--plan", action="store_true")
    ap.add_argument("--anchor", action="store_true")
    ap.add_argument("--folds-identity", action="store_true")
    ap.add_argument("--real", action="store_true")
    ap.add_argument("--raw", type=Path)
    ap.add_argument("--atlas", type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    if args.plan:
        for grid in ("s1", "s2", "s3"):
            print(f"{grid} cells: {len(synthetic_grid(grid))}")
        return 0

    if args.folds_identity:
        rng = np.random.default_rng(SEED)
        checks = []
        for spec in ((10, 5, 0.3), (30, 20, 0.5), (100, 5, 0.1)):
            y, lineage, _ = cohort(rng, spec[0], spec[1], spec[2], binary=True)
            checks.append({"n_groups": spec[0], "group_size": spec[1],
                           "share": spec[2],
                           "bench_settings": folds_identity(
                               y, lineage, SEED, BENCH_FOLDS, BENCH_REPEATS),
                           "vet_settings": folds_identity(
                               y, lineage, VET_SEED, VET_FOLDS, VET_REPEATS)})
        payload = {"provenance": provenance(), "checks": checks,
                   "all_identical": all(c[k]["identical"] for c in checks
                                        for k in ("bench_settings",
                                                  "vet_settings"))}
        (args.out / "folds_identity.json").write_text(
            json.dumps(payload, indent=1) + "\n", encoding="utf-8")
        print(json.dumps(payload, indent=1))
        return 0

    if args.anchor:
        payload = {"provenance": provenance(), "anchor": run_anchor()}
        (args.out / "arm_a_anchor.json").write_text(
            json.dumps(payload, indent=1) + "\n", encoding="utf-8")
        print(json.dumps({k: v for k, v in payload["anchor"].items()
                          if k != "cohorts"}, indent=1))
        return 0

    if args.real:
        if args.raw is None or args.atlas is None:
            raise SystemExit("--real needs --raw and --atlas")
        payload = {"provenance": provenance(),
                   **run_real(args.raw, args.atlas)}
        (args.out / "arm_b_real.json").write_text(
            json.dumps(payload, indent=1) + "\n", encoding="utf-8")
        print(f"cells: {len(payload['cells'])} "
              f"mismatched: {len(payload['mismatched'])}")
        return 0

    if args.synthetic is None or args.cell is None:
        raise SystemExit("choose --plan, --folds-identity, --anchor, --real, "
                         "or --synthetic with --cell")
    cells = synthetic_grid(args.synthetic)
    cell = cells[args.cell]
    print(f"[{args.synthetic}:{args.cell}] G={cell['n_groups']} "
          f"m={cell['group_size']} share={cell['share']} "
          f"prevalence={cell.get('prevalence', 0.25)}", flush=True)
    result = run_synthetic_cell(cell)
    target = args.out / f"cells_{args.synthetic}"
    target.mkdir(exist_ok=True)
    (target / f"{args.cell:03d}.json").write_text(
        json.dumps(result, indent=1) + "\n", encoding="utf-8")
    print(f"  accuracy={result['accuracy']:.4f}"
          f"+-{result['accuracy_mc_se']:.4f} "
          f"baseline={result['baseline_marginal']:.4f} "
          f"lift={result['lift']:+.4f} "
          f"kappa={result['kappa_adj']:+.4f}"
          f"+-{result['kappa_adj_mc_se']:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
