"""cli.py — command-line entry point.

    amr-clonalshare --config CONFIG [--results-dir DIR] [--validate]
                       [--seed N] [--k-select {auto,mdl,prediction_strength,gap,bic_mixture}]
                       [--threads N] [--quiet] [--no-check-files]

Exit codes
----------
0   analysis completed
1   ``--validate`` ran and the planted-truth recovery threshold was not met
2   configuration error
3   the analysis completed but a diagnostic gate failed: the fusion collapsed
    onto a single layer, the partition is the trait-absent stratum, or the
    post-clustering inference could not run. The results are still written; the
    non-zero status exists so a pipeline does not silently treat such a run as
    a finding.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Any, Dict

from . import __version__
from .config import ConfigError, load_config


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="amr-clonalshare",
        description="Multi-layer trait clustering for bacterial genotype data, "
                    "with artefact diagnostics and valid post-clustering inference.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--version", action="version",
                   version=f"amr-clonalshare {__version__}")
    p.add_argument("--config", required=True, help="path to a dataset YAML config")
    p.add_argument("--results-dir", default=None,
                   help="directory for cluster_result.json and archetype_profiles.tsv")
    p.add_argument("--validate", action="store_true",
                   help="run the planted-truth recovery check instead of the analysis")
    p.add_argument("--seed", type=int, default=42,
                   help="master seed; every stochastic stage is spawned from it")
    p.add_argument("--k-select", default="auto",
                   choices=["auto", "mdl", "prediction_strength", "gap", "bic_mixture"],
                   help="number-of-archetypes criterion (auto uses the config)")
    p.add_argument("--threads", type=int, default=None,
                   help="limit BLAS/OpenMP threads (sets OMP_NUM_THREADS et al.)")
    p.add_argument("--check-input", action="store_true",
                   help="load and check the input data, print the input check "
                        "in plain language, write input_qc.json and "
                        "input_qc.md to --results-dir if given, and stop "
                        "before any estimate")
    p.add_argument("--quiet", action="store_true", help="suppress the stdout summary")
    p.add_argument("--no-check-files", action="store_true",
                   help="skip the existence check on configured data files")
    return p


def _set_threads(n: int) -> None:
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[var] = str(n)


def _excludes_zero(ci) -> bool:
    return bool(ci) and len(ci) == 2 and ci[0] * ci[1] > 0


def _surveillance_summary(meta: dict) -> dict:
    """Compact reading of the lineage-aware surveillance block.

    The headline the full record buries is the *offsetting* case: an agent
    whose reported prevalence barely moved while both components moved a long
    way in opposite directions. A report quoting prevalence alone calls that
    agent stable, which is the one reading the decomposition exists to prevent.
    """
    out: Dict[str, Any] = {}
    lrp = meta.get("lineage_resolved_prevalence") or {}
    gaps = {f: r["difference_per_lineage_minus_per_isolate"]
            for f, r in lrp.items() if r.get("status") == "ok"}
    if gaps:
        widest = max(gaps, key=lambda f: abs(gaps[f]))
        out["lineage_prevalence_widest_gap_feature"] = widest
        out["lineage_prevalence_widest_gap"] = round(gaps[widest], 4)

    conc = meta.get("trait_concentration") or {}
    ok = [r for r in conc.values() if r.get("status") == "ok"]
    if ok:
        counts: Dict[str, int] = {}
        for r in ok:
            counts[r["direction"]] = counts.get(r["direction"], 0) + 1
        out["carriage_direction_counts"] = counts
        out["n_features_departing_from_proportional_carriage"] = sum(
            1 for r in ok if (r.get("p_value") or 1.0) < 0.05)

    dec = meta.get("prevalence_decomposition") or {}
    per = dec.get("per_feature") or {}
    rows = {f: r for f, r in per.items() if r.get("status") == "ok"}
    if rows:
        comp = {f for f, r in rows.items()
                if _excludes_zero(r.get("composition_ci95"))}
        within = {f for f, r in rows.items()
                  if _excludes_zero(r.get("within_lineage_ci95"))}
        offsetting = sorted(
            f for f in comp & within
            if rows[f]["composition"] * rows[f]["within_lineage"] < 0
            and abs(rows[f]["difference"]) < min(abs(rows[f]["composition"]),
                                                 abs(rows[f]["within_lineage"])))
        out["decomposition"] = {
            "contrast": f"{dec.get('contrast_column')}: "
                        f"{' vs '.join(dec.get('levels') or [])}",
            "n_features": len(rows),
            "n_composition_significant": len(comp),
            "n_within_lineage_significant": len(within),
            "n_offsetting": len(offsetting),
            "offsetting_features": offsetting,
            "note": "offsetting: both components significant, opposite in "
                    "sign, and each larger than the difference they produce",
        }
    return out


def _resolved(att: dict):
    """Does the 95 % interval put the attributable share on one side of the gate?

    A precision statement, not a verdict. ``True`` means the interval lies
    wholly below or wholly at-or-above the threshold, so a wider cohort would
    not change the reading; ``False`` means the cohort cannot yet separate the
    two. The gate itself is decided by the point estimate, for the reason given
    where ``lineage_triggered`` is computed.
    """
    lo, hi, gate = att.get("ci_low"), att.get("ci_high"), att.get("lambda_gate", 0.5)
    if lo is None or hi is None or lo != lo or hi != hi:
        return None
    return bool(hi < gate or lo >= gate)


def _summary(result: dict) -> dict:
    infer = result.get("post_clustering_inference") or {}
    disc = (infer.get("discreteness") or {})
    art = result.get("artifact_diagnostics") or {}
    infl = result.get("layer_influence") or {}
    lc = ((result.get("metadata_diagnostics") or {}).get("lineage_concordance")
          or {})
    att = ((result.get("metadata_diagnostics") or {}).get("lineage_attribution")
           or {})
    exact_p = infer.get("diagnostic_p_value_structure", infer.get("p_value_structure"))
    order_sensitivity = infer.get("p_value_structure_order_sensitivity") or {}
    max_order_p = order_sensitivity.get("max_over_orders")
    if infer.get("status") == "withheld_inadequate_split_design":
        p_report = "withheld: inadequate split design"
    elif max_order_p is not None and max_order_p < 1e-100:
        p_report = "p < 1e-100 across all checked split orderings"
    elif exact_p is not None and exact_p < 0.001:
        p_report = "p < 0.001"
    elif exact_p is not None:
        p_report = "p >= 0.001"
    else:
        p_report = None
    return {
        "selected_k": result.get("selected_k"),
        "n_isolates": result.get("n_isolates"),
        "no_structure": (result.get("k_selection") or {}).get("no_structure"),
        "inference_status": infer.get("status"),
        "p_value_structure_report": p_report,
        "diagnostic_p_value_structure_exact": exact_p,
        "structure_detected": infer.get("structure_detected"),
        "discrete_beyond_a_gradient": disc.get("discrete_beyond_a_gradient"),
        "discreteness_status": disc.get("status"),
        "discreteness_verdict": disc.get("verdict"),
        "discreteness_p_value": disc.get("p_value"),
        "continuum_bootstrap_exceedances": disc.get("bootstrap_exceedances"),
        "continuum_tail_probability_ci95": disc.get(
            "bootstrap_tail_probability_ci95"),
        "continuum_decision_resolved_at_alpha_0_05": disc.get(
            "resolved_at_alpha_0_05"),
        # Without these two, a small p-value next to a False verdict reads as a
        # bug rather than as the withholding it is.
        "continuum_latent_dimension": (disc.get("latent_dimension") or {}).get(
            "q_selected"),
        "continuum_null_under_dimensioned": (disc.get("latent_dimension") or {}).get(
            "at_boundary"),
        "split_design_adequate": infer.get("split_design_adequate"),
        "phenotype_n_fdr_significant": (result.get("phenotype_validation") or {}
                                        ).get("n_fdr_significant"),
        "phenotype_partition_beats_every_single_column_rule": (
            None if not (result.get("phenotype_validation") or {}).get("partition_beats")
            else all((result["phenotype_validation"]["partition_beats"]).values())),
        "n_defining_validated": infer.get("n_defining"),
        "n_defining_descriptive": result.get("n_defining_features_descriptive"),
        "fusion_n_eff_layers": infl.get("n_eff"),
        "fusion_collapse": infl.get("collapse"),
        "max_empty_stratum_ari": art.get("max_empty_stratum_ari"),
        "lineage_concordance_excess": (
            None if not lc or lc.get("status") != "ok"
            else round(lc["concordance_observed"] - lc["concordance_null_mean"], 4)),
        "lineage_concordance_z": lc.get("z"),
        # The magnitude, beside the significance test it replaces as the gate.
        "lineage_attributable_share": att.get("lam"),
        "lineage_attributable_share_ci95": (
            None if not att else [att.get("ci_low"), att.get("ci_high")]),
        "lineage_attribution_resolved_at_95ci": _resolved(att) if att else None,
        "r2_lineage_out_of_sample": att.get("r2_lineage"),
        "r2_partition_out_of_sample": att.get("r2_partition"),
        "shapley_lineage": att.get("shapley_lineage"),
        "shapley_partition": att.get("shapley_partition"),
        "clonal_share_all_features": (
            (result.get("metadata_diagnostics") or {})
            .get("clonal_share_by_layer", {})
            .get("all_clustering_features", {}).get("kappa_adj")),
        "surveillance": _surveillance_summary(
            result.get("metadata_diagnostics") or {}) or None,
        "lineage_confounded": (
            bool(att.get("lam", 0.0) >= float(att.get("lambda_gate", 0.5)))
            if att else
            None if not lc or lc.get("status") != "ok"
            else bool(lc.get("p_value", 1.0) < 0.05 and lc.get("z", 0) > 2)),
        "mdl_p_value": (result.get("mdl_calibration") or {}).get("p_MDL"),
        "mdl_gain_fraction_of_null_code": (
            (result.get("k_selection") or {}).get("mdl_gain_fraction")),
        "claim_level": (result.get("interpretation") or {}).get("claim_level"),
        "claim_status": (result.get("interpretation") or {}).get("claim_status"),
        "active_gate_codes": (result.get("interpretation") or {}).get(
            "active_gate_codes", []),
    }


def _apply_release_gates(result: dict) -> None:
    """Withhold inferential fields that fail a mandatory design gate.

    Exact numerical values remain available under explicitly diagnostic keys,
    but no public result artifact can simultaneously say that the split design
    is inadequate and release an ``ok`` status, a positive verdict, and a
    precise p-value.
    """
    infer = result.get("post_clustering_inference") or {}
    if infer.get("status") == "ok" and infer.get("split_design_adequate") is False:
        infer["diagnostic_status"] = infer.get("status")
        infer["diagnostic_p_value_structure"] = infer.get("p_value_structure")
        infer["diagnostic_structure_detected"] = infer.get("structure_detected")
        infer["status"] = "withheld_inadequate_split_design"
        infer["p_value_structure"] = None
        infer["structure_detected"] = False


def _gates(result: dict) -> list:
    """Diagnostics that make the reported partition uninterpretable if tripped."""
    problems = []
    infl = result.get("layer_influence") or {}
    if infl.get("collapse"):
        problems.append(
            f"fusion collapsed onto layer {infl.get('dominant_layer')!r} "
            f"(effective number of contributing layers "
            f"{infl.get('n_eff'):.2f} of {infl.get('n_layers')}); the result is "
            f"a single-layer clustering, not an integrated one")
    art = result.get("artifact_diagnostics") or {}
    if art.get("empty_stratum_warning"):
        problems.append(
            f"the partition largely reproduces the 'carries no feature of layer "
            f"{art.get('empty_stratum_layer')!r}' indicator "
            f"(ARI {art.get('max_empty_stratum_ari'):.2f}); this is an artefact "
            f"of the empty-union distance convention, not trait structure")
    infer = result.get("post_clustering_inference") or {}
    if infer.get("status") not in ("ok", "skipped", None):
        problems.append(f"post-clustering inference failed: {infer.get('status')}")
    md_diag = result.get("metadata_diagnostics") or {}
    lc = md_diag.get("lineage_concordance") or {}
    att = md_diag.get("lineage_attribution") or {}
    if att:
        # The gate is a magnitude. The concordance z is reported with it and no
        # longer carries the threshold: at unchanged structure it runs from 1.5
        # at n = 67 to 25.5 at n = 677 on the shipped S. suis cohort, because
        # its null standard deviation shrinks with the number of pairs while
        # the share of the partition that lineage explains does not move.
        gate = float(att.get("lambda_gate", 0.5))
        lam = att.get("lam")

        # A record read back from JSON carries null where the run held NaN,
        # so a missing number is formatted from None as well as from NaN.
        def num(key):
            v = att.get(key)
            return float("nan") if v is None else float(v)

        if lam is not None and lam >= gate:
            problems.append(
                f"the partition is a lineage relabelling: {lam:.2f} of what it "
                f"explains about the traits (95 % CI "
                f"{num('ci_low'):.2f} to {num('ci_high'):.2f}) is attributable to "
                f"the lineage label, at or above the registered threshold of "
                f"{gate:.2f}; out-of-sample variance explained is "
                f"{num('r2_lineage'):.3f} for lineage alone, "
                f"{num('r2_partition'):.3f} for the partition "
                f"alone and {num('r2_joint'):.3f} for both")
        # An interval that straddles the threshold is not a refusal. The point
        # estimate decides; whether the cohort is large enough to put the whole
        # interval on one side is a precision statement, carried by
        # `lineage_attribution_resolved_at_95ci` in the summary and by
        # `resolved_at_95ci` in the gate record. Making it block would restore
        # exactly the sample-size dependence this gate removes.
    elif lc.get("status") == "ok" and lc.get("p_value", 1.0) < 0.05 and lc.get("z", 0) > 2:
        problems.append(
            f"the partition is lineage-confounded: same-lineage isolates are "
            f"co-clustered {lc['concordance_observed']:.3f} of the time against a "
            f"permutation null of {lc['concordance_null_mean']:.3f} "
            f"(z = {lc['z']:.1f}); the clusters may be a lineage relabelling "
            f"rather than recurrent trait combinations. No attribution was "
            f"computed, so this significance test is all that is available and "
            f"it carries no effect size")
    if infer.get("status") == "ok" and infer.get("split_design_adequate") is False:
        problems.append(
            "the feature-split design is not usable: "
            + "; ".join(infer.get("split_design_problems") or ["unspecified"]))

    ph = result.get("phenotype_validation") or {}
    if ph.get("status") == "ok":
        lost = [n for n, won in (ph.get("partition_beats") or {}).items() if not won]
        if lost:
            ba = ph.get("mean_balanced_accuracy_over_fdr_family") or {}
            boot = ph.get("head_to_head_bootstrap") or {}
            best = max(lost, key=lambda n: ba.get(n, 0.0))
            b = boot.get(best, {})
            floor = (ph.get("orientation_null_floor") or {}).get("partition")
            problems.append(
                f"the partition is significantly beaten on measured phenotype by "
                f"the single rule '{best}': mean balanced accuracy "
                f"{ba.get('partition', float('nan')):.3f} against "
                f"{ba.get(best, float('nan')):.3f} over the "
                f"{ph.get('n_in_fdr_family')} antibiotics with enough tested "
                f"isolates (permutation floor {floor:.3f}), difference "
                f"{b.get('mean_difference_minus_partition', float('nan')):+.3f} "
                f"with bootstrap 95 % CI "
                f"[{b.get('ci95_low', float('nan')):+.3f}, "
                f"{b.get('ci95_high', float('nan')):+.3f}]. Whatever the fusion "
                f"found, a laboratory could have got it from one column")

    disc = infer.get("discreteness") or {}
    if disc.get("status") == "withheld_under_dimensioned":
        dim = disc.get("latent_dimension") or {}
        problems.append(
            f"the discreteness verdict is WITHHELD, not negative: BIC selects "
            f"a {dim.get('q_selected')}-dimensional latent continuum at "
            f"q_max = {dim.get('q_max')}, so the continuum null is known to be "
            "under-dimensioned; no bootstrap p-value was released")
    elif disc.get("status") == "ok" and disc.get("discrete_beyond_a_gradient") is False:
        dim = disc.get("latent_dimension") or {}
        problems.append(
            f"structure was detected but is not distinguishable from a "
            f"continuous gradient (continuum-null p = {disc.get('p_value')}, "
            f"latent dimension {dim.get('q_selected')}); reporting discrete "
            f"archetypes is not supported")
    return problems


def _interpretation_contract(result: dict) -> dict:
    """Return a machine-readable claim ladder and threshold ledger.

    The ladder prevents downstream code from promoting a partition directly to
    a biological archetype. ``archetype_candidate`` is deliberately narrower
    than biological confirmation: it means that all diagnostics applicable to
    this run passed, not that causality, clinical utility or transportability
    has been established.
    """
    infer = result.get("post_clustering_inference") or {}
    disc = infer.get("discreteness") or {}
    infl = result.get("layer_influence") or {}
    art = result.get("artifact_diagnostics") or {}
    lc = ((result.get("metadata_diagnostics") or {}).get("lineage_concordance")
          or {})
    att = ((result.get("metadata_diagnostics") or {}).get("lineage_attribution")
           or {})
    ph = result.get("phenotype_validation") or {}
    k = int(result.get("selected_k") or 1)

    # The gate is the magnitude where one is available. The permutation z is
    # still recorded, but it grows with the number of pairs at fixed structure
    # and a threshold on it refuses cohorts for being large.
    #
    # The point estimate decides and the interval does not. Requiring the whole
    # interval to clear the threshold would make the verdict a function of how
    # wide the bootstrap is, and the bootstrap is wide when the cohort is small
    # -- which is the sample-size dependence this gate exists to remove, in
    # mirror image. Whether the interval excludes the threshold is a separate
    # statement about precision, reported as `resolved_at_95ci`, and the
    # Kitagawa gate on the other half of the package draws the same line.
    if att and att.get("lam") is not None:
        lineage_triggered = bool(att["lam"] >= float(att.get("lambda_gate", 0.5)))
    else:
        lineage_triggered = bool(
            lc.get("status") == "ok"
            and lc.get("p_value", 1.0) < 0.05
            and lc.get("z", 0.0) > 2.0
        )
    phenotype_lost = [
        name for name, won in (ph.get("partition_beats") or {}).items()
        if not won
    ] if ph.get("status") == "ok" else []
    best_bootstrap = None
    if phenotype_lost:
        scores = ph.get("mean_balanced_accuracy_over_fdr_family") or {}
        best = max(phenotype_lost, key=lambda name: scores.get(name, 0.0))
        best_bootstrap = (ph.get("head_to_head_bootstrap") or {}).get(best) or {}

    records = [
        {
            "code": "fusion_collapse",
            "applicable": bool(infl),
            "triggered": bool(infl.get("collapse")),
            "value": infl.get("n_eff"),
            "threshold": infl.get("collapse_threshold", 1.5),
            "margin_to_failure": (
                None if infl.get("n_eff") is None else
                float(infl["n_eff"] - infl.get("collapse_threshold", 1.5))
            ),
        },
        {
            "code": "empty_stratum",
            "applicable": art.get("max_empty_stratum_ari") is not None,
            "triggered": bool(art.get("empty_stratum_warning")),
            "value": art.get("max_empty_stratum_ari"),
            "threshold": 0.5,
            "margin_to_failure": (
                None if art.get("max_empty_stratum_ari") is None else
                float(0.5 - art["max_empty_stratum_ari"])
            ),
        },
        {
            "code": "post_clustering_inference",
            "applicable": k > 1,
            "triggered": bool(k > 1 and infer.get("status") != "ok"),
            "value": infer.get("status"),
            "threshold": "status == ok",
            "margin_to_failure": None,
        },
        {
            "code": "split_design",
            "applicable": k > 1,
            "triggered": bool(k > 1 and infer.get("split_design_adequate") is False),
            "value": infer.get("split_design_adequate"),
            "threshold": "adequate == true",
            "margin_to_failure": None,
        },
        {
            "code": "lineage_confounding",
            "applicable": bool(att) or lc.get("status") == "ok",
            "triggered": lineage_triggered,
            # When an attribution is available the gate is the share, and the
            # concordance test is carried alongside as a descriptor. It stays
            # in the record because it is the statistic the literature reports
            # and a reader comparing runs needs to see both.
            "value": (
                {"lineage_attributable_share": att.get("lam"),
                 "ci95": [att.get("ci_low"), att.get("ci_high")],
                 "r2_lineage": att.get("r2_lineage"),
                 "r2_partition": att.get("r2_partition"),
                 "r2_joint": att.get("r2_joint"),
                 "shapley_lineage": att.get("shapley_lineage"),
                 "shapley_partition": att.get("shapley_partition"),
                 "resolved_at_95ci": _resolved(att),
                 "concordance_z_not_gated": lc.get("z")}
                if att else {"p_value": lc.get("p_value"), "z": lc.get("z")}
            ),
            "threshold": (
                {"lineage_attributable_share_below": att.get("lambda_gate", 0.5)}
                if att else {"p_value_below": 0.05, "z_above": 2.0}
            ),
            "margin_to_failure": (
                float(att.get("lambda_gate", 0.5) - att["lam"])
                if att and att.get("lam") is not None else
                None if lc.get("status") != "ok" or lc.get("z") is None else
                float(2.0 - lc["z"])
            ),
        },
        {
            "code": "phenotype_superiority",
            "applicable": ph.get("status") == "ok",
            "triggered": bool(phenotype_lost),
            "value": {
                "rules_beating_partition": phenotype_lost,
                "worst_ci95_low_minus_partition": (
                    None if not best_bootstrap else
                    best_bootstrap.get("ci95_low")
                ),
            },
            "threshold": "no tested simple rule has bootstrap CI95 low > 0",
            "margin_to_failure": (
                None if not best_bootstrap else
                -float(best_bootstrap.get("ci95_low", 0.0))
            ),
        },
        {
            "code": "continuum_discreteness",
            "applicable": k > 1 and bool(infer.get("structure_detected")),
            "triggered": bool(
                k > 1 and bool(infer.get("structure_detected"))
                and not (disc.get("status") == "ok"
                         and disc.get("discrete_beyond_a_gradient") is True)
            ),
            "value": {
                "status": disc.get("status"),
                "p_value": disc.get("p_value"),
                "tail_probability_ci95": disc.get(
                    "bootstrap_tail_probability_ci95"),
            },
            "threshold": "status == ok and p_value < 0.05",
            "margin_to_failure": (
                None if disc.get("p_value") is None else
                float(0.05 - disc["p_value"])
            ),
        },
    ]
    active = [r["code"] for r in records if r["applicable"] and r["triggered"]]
    no_structure = bool((result.get("k_selection") or {}).get("no_structure"))
    if no_structure or k <= 1:
        level, status = 0, "no_structure"
    elif infer.get("status") != "ok" or not infer.get("structure_detected"):
        level, status = 1, "descriptive_partition"
    elif disc.get("status") != "ok" or not disc.get("discrete_beyond_a_gradient"):
        level, status = 2, "validated_structure"
    elif active:
        level, status = 3, "statistically_discrete_partition"
    else:
        level, status = 4, "archetype_candidate"
    return {
        "contract_version": "1.0",
        "claim_level": level,
        "claim_status": status,
        "active_gate_codes": active,
        "gate_ledger": records,
        "scope_note": (
            "Highest claim supported by diagnostics applicable to this run; "
            "archetype_candidate is not causal, clinical or cross-cohort "
            "biological confirmation."
        ),
    }


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.threads:
        _set_threads(args.threads)

    # NumPy/SciPy read BLAS/OpenMP limits while their extension modules are
    # imported. Importing ``core`` at module load made ``--threads`` cosmetic:
    # by the time the variables were set, the thread pools already existed.
    # Delay the numerical stack until after the execution contract is applied.
    from . import core
    from .jsonio import dumps, write_json

    try:
        cfg = load_config(args.config, check_files_exist=not args.no_check_files)
    except ConfigError as exc:
        print(f"config error: {exc}", file=sys.stderr)
        return 2

    if args.validate and not args.check_input:
        res = core.validate(cfg, seed=args.seed)
        print(dumps(res))
        return 0 if res.get("passed", False) else 1

    rd = Path(args.results_dir).expanduser().resolve() if args.results_dir else None
    if rd is not None:
        rd.mkdir(parents=True, exist_ok=True)

    # The input check runs before any estimate: a layer with empty cells, a
    # value that is not 0 or 1, or a lineage column that cannot support the
    # estimator is reported here in plain language, and a refusal leaves with
    # the config exit code rather than a traceback.
    from .io import load_dataset
    from .qc import render_markdown
    try:
        ds = load_dataset(cfg)
    except ConfigError as exc:
        print(f"input error: {exc}", file=sys.stderr)
        return 2
    assert ds.input_qc is not None
    qc_text = render_markdown(ds.input_qc)
    if rd is not None:
        write_json(ds.input_qc, rd / "input_qc.json")
        (rd / "input_qc.md").write_text(qc_text, encoding="utf-8")
    if args.check_input:
        print(qc_text)
        return 0
    input_record = ds.input_qc
    del ds

    result = core.run(cfg, results_dir=rd, seed=args.seed, k_select=args.k_select)
    _apply_release_gates(result)
    result["interpretation"] = _interpretation_contract(result)

    problems = _gates(result)
    if rd is not None:
        from .report import render_report
        from .report_html import render_html_report
        write_json(result, rd / "cluster_result.json")
        write_json(_summary(result), rd / "summary.json")
        (rd / "report.md").write_text(
            render_report(result, _summary(result), input_qc=input_record,
                          problems=problems), encoding="utf-8")
        # The same run, drawn. The prose form states the result; the drawn
        # form is where a reader sees that every per-trait interval crosses
        # zero, or that one gate held by 0.04 while the rest held by 0.5.
        (rd / "report.html").write_text(
            render_html_report(result, _summary(result)), encoding="utf-8")

    if not args.quiet:
        print("[amr-clonalshare] summary:")
        print(dumps(_summary(result)))
        if rd is not None:
            print(f"[amr-clonalshare] outputs written to {rd}")
        for p in problems:
            print(f"[amr-clonalshare] DIAGNOSTIC FAILURE: {p}", file=sys.stderr)
    return 3 if problems else 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
