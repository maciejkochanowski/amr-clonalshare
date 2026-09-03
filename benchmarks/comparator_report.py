#!/usr/bin/env python3
"""Aggregate the comparator arms and score them against what was registered.

    python comparator_report.py --evidence <dir>

WHAT THIS DOES AND DOES NOT DECIDE. Every threshold and every predicted
direction read here comes from prereg_comparators.json, whose hash is carried
in the output. Nothing is chosen after seeing a number. An outcome that failed
is written as failed and stays in the table, because a registration that only
ever records confirmations is a record of nothing.

WHY THE CORRELATIONS ARE RANK CORRELATIONS. The clonal share and a
classification accuracy are on different scales and neither is a linear
function of the other under any model anyone proposes, so the question the
manuscript has to answer is whether the two order the cells the same way. That
is a rank question, and Spearman answers it. The interval on each coefficient
is a bootstrap over cells rather than the usual normal approximation, because
the cells are not independent draws from one population: many share an organism
and a host, and the resampling makes that dependence visible in the width
rather than hiding it in a formula.

WHY BOTH ACCURACY AND ACCURACY MINUS BASELINE. Raw accuracy contains the
prevalence. Two cells with the same lineage structure and different prevalence
have different accuracies, and two cells with the same accuracy can have
opposite amounts of lineage structure. Subtracting the majority baseline
removes the part of the accuracy that any predictor gets for free. If the
correlation with the share is stronger after that subtraction, the part that
was removed was noise with respect to the share, which is the whole claim.
"""
from __future__ import annotations

import argparse
import json
import platform
from datetime import date, datetime, timezone
from pathlib import Path

import numpy as np
import scipy
from scipy import stats

PREREG_SHA256 = ("4d54749cf39ed17db74ad997fec31dd617"
                 "38a5500c77961b8809627780438df9")
ADDENDUM_SHA256 = ("27b653a0e23aafd8b6575f87fdd1b973c"
                   "75e6515c03a30bed249b64a7a22b360")
#: The parent registration set a floor of thirty cells for a rank correlation
#: and in the same document registered a sweep of twenty five cells, so the
#: sweep fell below the parent's own floor. The addendum repairs that with a
#: forty five cell grid, and this constant is what the repair is scored on.
CELL_FLOOR = 30
BOOTSTRAP = 4000
BOOTSTRAP_SEED = 20260902


def _spearman(x, y, *, label: str, seed: int = BOOTSTRAP_SEED) -> dict:
    """Spearman rho with a bootstrap interval over cells."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 3:
        return {"label": label, "n_cells": int(x.size), "rho": float("nan"),
                "p_value": float("nan"), "ci_low": float("nan"),
                "ci_high": float("nan"), "reportable": False}
    result = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(BOOTSTRAP):
        idx = rng.integers(0, x.size, x.size)
        if np.unique(x[idx]).size < 3 or np.unique(y[idx]).size < 3:
            continue
        draws.append(float(stats.spearmanr(x[idx], y[idx]).statistic))
    draws = np.asarray(draws, dtype=float)
    return {"label": label, "n_cells": int(x.size),
            "rho": float(result.statistic), "p_value": float(result.pvalue),
            "ci_low": float(np.percentile(draws, 2.5)),
            "ci_high": float(np.percentile(draws, 97.5)),
            "bootstrap_se": float(draws.std(ddof=1)),
            "reportable": bool(x.size >= 30)}


def _paired_difference(x, a, b, *, label: str) -> dict:
    """Bootstrap the difference of two rank correlations on the same cells.

    The registered outcome is a comparison between two coefficients computed on
    one set of cells, so the two are not independent and their intervals cannot
    be compared by eye. The difference is therefore resampled directly.
    """
    x, a, b = (np.asarray(v, dtype=float) for v in (x, a, b))
    ok = np.isfinite(x) & np.isfinite(a) & np.isfinite(b)
    x, a, b = x[ok], a[ok], b[ok]
    rng = np.random.default_rng(BOOTSTRAP_SEED + 1)
    observed = (abs(float(stats.spearmanr(x, b).statistic))
                - abs(float(stats.spearmanr(x, a).statistic)))
    draws = []
    for _ in range(BOOTSTRAP):
        idx = rng.integers(0, x.size, x.size)
        if np.unique(x[idx]).size < 3:
            continue
        draws.append(abs(float(stats.spearmanr(x[idx], b[idx]).statistic))
                     - abs(float(stats.spearmanr(x[idx], a[idx]).statistic)))
    draws = np.asarray(draws, dtype=float)
    return {"label": label, "n_cells": int(x.size),
            "difference_of_absolute_rho": observed,
            "ci_low": float(np.percentile(draws, 2.5)),
            "ci_high": float(np.percentile(draws, 97.5)),
            "share_of_draws_positive": float((draws > 0).mean())}


def _rate(values, mc_ses=None) -> dict:
    """A mean over cells, with the two errors that apply to it kept apart.

    The dispersion across cells is a property of the cells and would not shrink
    if the computation were run longer. The Monte-Carlo error is the part that
    would. Reporting one in place of the other is how a campaign concludes that
    an effect is absent when it only failed to resolve it.
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    out = {"n": int(v.size), "mean": float(v.mean()),
           "median": float(np.median(v)),
           "quartiles": [float(x) for x in np.percentile(v, [25, 75])],
           "min": float(v.min()), "max": float(v.max()),
           "between_cell_se_of_the_mean": float(v.std(ddof=1)
                                                / np.sqrt(v.size))}
    if mc_ses is not None:
        m = np.asarray(mc_ses, dtype=float)
        m = m[np.isfinite(m)]
        out["monte_carlo_se_of_the_mean"] = float(
            np.sqrt((m ** 2).sum()) / m.size) if m.size else float("nan")
        out["median_per_cell_monte_carlo_se"] = (float(np.median(m))
                                                 if m.size else float("nan"))
    return out


def _arm_b_block(cells: list[dict], *, name: str) -> dict:
    """Every Arm B quantity for one grid, plus the correlations."""
    kappa = [c["kappa"] for c in cells]
    accuracy = [c["accuracy"] for c in cells]
    balanced = [c["balanced_accuracy"] for c in cells]
    baseline = [c["baseline_marginal"] for c in cells]
    lift = [c["lift"] for c in cells]
    lift_oof = [c["lift_over_out_of_fold_baseline"] for c in cells]
    high = [c for c in cells if c["baseline_marginal"] >= 0.85]
    block = {
        "grid": name, "n_cells": len(cells),
        "accuracy": _rate(accuracy, [c["accuracy_mc_se"] for c in cells]),
        "balanced_accuracy": _rate(balanced,
                                   [c["balanced_accuracy_mc_se"]
                                    for c in cells]),
        "baseline_marginal": _rate(baseline),
        "baseline_out_of_fold": _rate([c["baseline_out_of_fold"]
                                       for c in cells],
                                      [c["baseline_out_of_fold_mc_se"]
                                       for c in cells]),
        "lift_over_marginal_baseline": _rate(lift),
        "lift_over_out_of_fold_baseline": _rate(lift_oof),
        "clonal_share": _rate(kappa),
        "share_of_cells_where_accuracy_beats_out_of_fold_baseline": float(
            np.mean([c["accuracy"] > c["baseline_out_of_fold"]
                     for c in cells])),
        "share_of_cells_where_accuracy_beats_marginal_baseline": float(
            np.mean([c["accuracy"] > c["baseline_marginal"] for c in cells])),
        "high_baseline_subset": {
            "definition": "cells whose majority baseline is at least 0.85",
            "n_cells": len(high),
            "median_lift": (float(np.median([c["lift"] for c in high]))
                            if high else float("nan")),
            "median_accuracy": (float(np.median([c["accuracy"]
                                                 for c in high]))
                                if high else float("nan")),
            "median_baseline": (float(np.median([c["baseline_marginal"]
                                                 for c in high]))
                                if high else float("nan")),
            "median_accuracy_minus_balanced_accuracy": (
                float(np.median([c["accuracy"] - c["balanced_accuracy"]
                                 for c in high])) if high else float("nan")),
        },
        "rank_correlations": {
            "share_vs_accuracy": _spearman(kappa, accuracy,
                                           label="clonal share against raw "
                                                 "accuracy"),
            "share_vs_lift": _spearman(kappa, lift,
                                       label="clonal share against accuracy "
                                             "minus its majority baseline"),
            "share_vs_lift_out_of_fold": _spearman(
                kappa, lift_oof,
                label="clonal share against accuracy minus the out-of-fold "
                      "baseline"),
            "baseline_vs_accuracy": _spearman(
                baseline, accuracy,
                label="majority baseline against raw accuracy"),
            "share_vs_balanced_accuracy_exploratory": _spearman(
                kappa, balanced,
                label="clonal share against balanced accuracy, exploratory"),
        },
        "paired_difference": _paired_difference(
            kappa, accuracy, lift,
            label="absolute rho for accuracy minus baseline, minus absolute "
                  "rho for raw accuracy"),
    }
    return block


def _load_synthetic(evidence: Path, grid: str) -> list[dict]:
    """Synthetic cells, with the estimate named the way the real cells name it.

    The synthetic run carries the bias-corrected clonal share under its own
    field name kappa_adj, and the published atlas carries the same quantity
    under kappa. Renaming here rather than in either producer keeps both
    artefacts as they were written.
    """
    rows = []
    for path in sorted((evidence / f"cells_{grid}").glob("*.json")):
        cell = json.loads(path.read_text())
        cell["kappa"] = cell["kappa_adj"]
        rows.append(cell)
    return rows


def _verdict(identifier: str, statement: str, met: bool, value) -> dict:
    return {"outcome": identifier, "registered": statement,
            "met": bool(met), "measured": value}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--evidence", type=Path, required=True)
    args = ap.parse_args()
    evidence = args.evidence

    leverage = json.loads((evidence / "leverage_check.json").read_text())
    identity = json.loads((evidence / "folds_identity.json").read_text())
    anchor = json.loads((evidence / "arm_a_anchor.json").read_text())["anchor"]
    real = json.loads((evidence / "arm_b_real.json").read_text())
    s1 = _load_synthetic(evidence, "s1")
    s2 = _load_synthetic(evidence, "s2")
    s3 = _load_synthetic(evidence, "s3")

    # The published atlas carries human cohorts beside the animal ones, because
    # the cut by host is the thing it measures. The registered outcomes are
    # scored on all of its cells. The animal-only split is reported beside them
    # as an exploratory robustness check with no registered direction, since a
    # veterinary reader will ask whether the human cells carried the result.
    animal = [c for c in real["cells"] if c["host_group"] != "human"]
    blocks = {"R1_real": _arm_b_block(real["cells"], name="R1_real"),
              "R1_real_animal_only_exploratory": _arm_b_block(
                  animal, name="R1_real_animal_only_exploratory"),
              "S1_shipped_binary": _arm_b_block(s1, name="S1_shipped_binary"),
              "S2_prevalence_sweep": _arm_b_block(s2,
                                                  name="S2_prevalence_sweep"),
              "S3_prevalence_sweep_extended": _arm_b_block(
                  s3, name="S3_prevalence_sweep_extended")}

    r1, sweep = blocks["R1_real"], blocks["S2_prevalence_sweep"]
    r1_rho = r1["rank_correlations"]
    s1_rho = blocks["S1_shipped_binary"]["rank_correlations"]
    s2_rho = sweep["rank_correlations"]
    s3_rho = blocks["S3_prevalence_sweep_extended"]["rank_correlations"]
    null_high = [c for c in s2 if c["prevalence_target"] == 0.90
                 and c["share"] == 0.0]
    lift_by = {}
    for cell in s2:
        lift_by[(cell["prevalence_target"], cell["share"])] = cell["lift"]
    monotone = [lift_by[(0.90, s)] < lift_by[(0.50, s)]
                for s in (0.1, 0.3, 0.5, 0.7)]

    verdicts = [
        _verdict("B-001", "at least 0.20 of the real cells carry a majority "
                          "baseline of 0.80 or more",
                 leverage["R1_real"]["leverage"]["share_baseline_at_least_0.80"]
                 >= 0.20,
                 leverage["R1_real"]["leverage"]["share_baseline_at_least_0.80"]),
        _verdict("B-002", "in real cells with a baseline of 0.85 or more the "
                          "median lift is below 0.05",
                 r1["high_baseline_subset"]["median_lift"] < 0.05,
                 r1["high_baseline_subset"]["median_lift"]),
        _verdict("B-003", "on the real cells the share correlates with raw "
                          "accuracy more weakly than with accuracy minus "
                          "baseline",
                 abs(r1_rho["share_vs_accuracy"]["rho"])
                 < abs(r1_rho["share_vs_lift"]["rho"]),
                 [r1_rho["share_vs_accuracy"]["rho"],
                  r1_rho["share_vs_lift"]["rho"]]),
        _verdict("B-004", "on the real cells the majority baseline and raw "
                          "accuracy correlate at 0.5 or more",
                 r1_rho["baseline_vs_accuracy"]["rho"] >= 0.5,
                 r1_rho["baseline_vs_accuracy"]["rho"]),
        _verdict("B-005", "on the real cells the share and raw accuracy "
                          "correlate below 0.7 in absolute value",
                 abs(r1_rho["share_vs_accuracy"]["rho"]) < 0.7,
                 r1_rho["share_vs_accuracy"]["rho"]),
        _verdict("B-006", "in real cells with a baseline of 0.85 or more, "
                          "plain accuracy exceeds balanced accuracy by a "
                          "median of more than 0.05",
                 r1["high_baseline_subset"][
                     "median_accuracy_minus_balanced_accuracy"] > 0.05,
                 r1["high_baseline_subset"][
                     "median_accuracy_minus_balanced_accuracy"]),
        _verdict("S-001", "the shipped binary grid spans a baseline range "
                          "below 0.30",
                 leverage["S1_shipped_binary"]["leverage"]["baseline_range"]
                 < 0.30,
                 leverage["S1_shipped_binary"]["leverage"]["baseline_range"]),
        _verdict("S-002", "on the shipped binary grid the share and raw "
                          "accuracy correlate at 0.8 or more in absolute value",
                 abs(s1_rho["share_vs_accuracy"]["rho"]) >= 0.8,
                 s1_rho["share_vs_accuracy"]["rho"]),
        _verdict("S-003", "on the prevalence sweep the share correlates with "
                          "raw accuracy more weakly than with accuracy minus "
                          "baseline",
                 abs(s2_rho["share_vs_accuracy"]["rho"])
                 < abs(s2_rho["share_vs_lift"]["rho"]),
                 [s2_rho["share_vs_accuracy"]["rho"],
                  s2_rho["share_vs_lift"]["rho"]]),
        _verdict("S-004", "at prevalence 0.90 and true share zero, accuracy "
                          "is within 0.02 of the baseline and the share is "
                          "within 0.05 of zero",
                 bool(null_high) and all(
                     abs(c["accuracy"] - c["baseline_marginal"]) <= 0.02
                     and abs(c["kappa"]) <= 0.05 for c in null_high),
                 [{"accuracy": c["accuracy"],
                   "baseline": c["baseline_marginal"],
                   "kappa": c["kappa"]} for c in null_high]),
        _verdict("S-005", "on the prevalence sweep the lift at prevalence "
                          "0.90 is below the lift at prevalence 0.50 at every "
                          "non-zero true share",
                 all(monotone),
                 {str(s): [lift_by[(0.90, s)], lift_by[(0.50, s)]]
                  for s in (0.1, 0.3, 0.5, 0.7)}),
        _verdict("A-001", "the anchor REML agrees with the shipped reml_lmm "
                          "to 0.01 on the share scale on every verification "
                          "cohort",
                 anchor["agrees_within_tolerance"],
                 anchor["max_absolute_difference"]),
        _verdict("S-006", "ADDENDUM. on the extended sweep, which clears the "
                          "thirty cell floor the parent registration set, the "
                          "share correlates with raw accuracy more weakly "
                          "than with accuracy minus baseline",
                 abs(s3_rho["share_vs_accuracy"]["rho"])
                 < abs(s3_rho["share_vs_lift"]["rho"]),
                 [s3_rho["share_vs_accuracy"]["rho"],
                  s3_rho["share_vs_lift"]["rho"]]),
        _verdict("S-007", "ADDENDUM. on the extended sweep the majority "
                          "baseline and raw accuracy correlate at 0.5 or more "
                          "in absolute value",
                 abs(s3_rho["baseline_vs_accuracy"]["rho"]) >= 0.5,
                 s3_rho["baseline_vs_accuracy"]["rho"]),
        _verdict("S-008", "ADDENDUM. the extended sweep spans a baseline "
                          "range of at least 0.30",
                 (blocks["S3_prevalence_sweep_extended"]["baseline_marginal"]
                  ["max"]
                  - blocks["S3_prevalence_sweep_extended"]["baseline_marginal"]
                  ["min"]) >= 0.30,
                 blocks["S3_prevalence_sweep_extended"]["baseline_marginal"]
                 ["max"]
                 - blocks["S3_prevalence_sweep_extended"]["baseline_marginal"]
                 ["min"]),
    ]


    void = {
        "prevalence_leverage": {
            grid: leverage[grid]["leverage"]["void_for_this_comparison"]
            for grid in ("R1_real", "S1_shipped_binary",
                         "S2_prevalence_sweep")},
        "fold_scheme_identical": identity["all_identical"],
        "fold_scheme_identical_on_the_real_cells": (
            real["folds_identity"]["identical"]
            if real.get("folds_identity") else None),
        "isolate_sets_matched": {
            "cells_kept": len(real["cells"]),
            "cells_voided_for_a_count_mismatch": len(real["mismatched"]),
            "detail": real["mismatched"]},
        "predictor_has_leverage": {
            grid: blocks[grid][
                "share_of_cells_where_accuracy_beats_out_of_fold_baseline"]
            for grid in blocks},
        "enough_cells_for_a_rank_correlation": {
            grid: bool(blocks[grid]["n_cells"] >= CELL_FLOOR)
            for grid in blocks},
        "registration_defect": (
            "The parent registration set a floor of thirty cells for a rank "
            "correlation and in the same document registered the prevalence "
            "sweep S2 at twenty five cells, so by its own terms the S2 rank "
            "correlations are below the floor and are reported as description "
            "rather than as evidence. The addendum registers the extended "
            "sweep S3 at forty five cells before it was run, and S-006 is the "
            "version of that comparison that may be quoted. Both grids are "
            "reported and neither is discarded."),
        "anchor_reproduced_the_shipped_estimator":
            anchor["agrees_within_tolerance"],
    }

    payload = {
        "provenance": {
            "generated": date.today().isoformat(),
            "generated_utc": datetime.now(timezone.utc).strftime(
                "%Y-%m-%dT%H:%M:%SZ"),
            "prereg_sha256": PREREG_SHA256,
            "prereg_addendum_sha256": ADDENDUM_SHA256,
            "python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__, "bootstrap_replicates": BOOTSTRAP,
            "arm_provenance": {"real": real["provenance"],
                               "anchor": json.loads(
                                   (evidence / "arm_a_anchor.json").read_text()
                               )["provenance"]}},
        "void_conditions": void,
        "leverage": {grid: leverage[grid]["leverage"] for grid in
                     ("R1_real", "S1_shipped_binary", "S2_prevalence_sweep")},
        "arm_a_mixed_model": {
            "what_the_shipped_arm_fits": (
                "restricted maximum likelihood for the one-way random-effects "
                "model with a lineage random intercept, profiled on the ratio "
                "of the between-lineage variance to the within-lineage one, "
                "reported as ratio over one plus ratio"),
            "anchor_tolerance": anchor["tolerance"],
            "max_absolute_difference": anchor["max_absolute_difference"],
            "mean_absolute_difference": anchor["mean_absolute_difference"],
            "agrees_within_tolerance": anchor["agrees_within_tolerance"],
            "n_verification_cohorts": len(anchor["cohorts"]),
            "what_it_is_not": (
                "the genomic-relatedness version of the same estimator, in "
                "which the kinship is built from genome-wide variants rather "
                "than from a discrete cluster label")},
        "arm_b_accuracy": blocks,
        "registered_outcomes": verdicts,
        "outcomes_met": sum(1 for v in verdicts if v["met"]),
        "outcomes_not_met": sum(1 for v in verdicts if not v["met"]),
    }
    (evidence / "comparator_report.json").write_text(
        json.dumps(payload, indent=1) + "\n", encoding="utf-8")

    print(json.dumps(void, indent=1))
    for grid in blocks:
        block = blocks[grid]
        print(f"\n{grid}: {block['n_cells']} cells")
        print(f"  accuracy        {block['accuracy']['mean']:.4f} "
              f"mc_se {block['accuracy']['monte_carlo_se_of_the_mean']:.5f} "
              f"between-cell se {block['accuracy']['between_cell_se_of_the_mean']:.4f}")
        print(f"  balanced        {block['balanced_accuracy']['mean']:.4f}")
        print(f"  baseline        {block['baseline_marginal']['mean']:.4f}")
        print(f"  lift            {block['lift_over_marginal_baseline']['mean']:+.4f} "
              f"median {block['lift_over_marginal_baseline']['median']:+.4f}")
        for key, rho in block["rank_correlations"].items():
            print(f"  rho {key:<42} {rho['rho']:+.4f} "
                  f"[{rho['ci_low']:+.4f}, {rho['ci_high']:+.4f}] "
                  f"n={rho['n_cells']}")
    print("\nregistered outcomes")
    for v in verdicts:
        print(f"  {v['outcome']}  {'MET    ' if v['met'] else 'NOT MET'}  "
              f"{v['measured']}")
    print(f"\nmet {payload['outcomes_met']} of {len(verdicts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
