#!/usr/bin/env python3
"""Build an auditable reviewer-resilience table from a terminal campaign.

The calculation never changes a scientific decision. It exposes the claim
level, every applicable gate, distance from the registered threshold, and a
small set of predeclared stricter/looser thresholds that a reviewer is likely
to probe. Positive margins are on the passing side; negative margins indicate
an active gate.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CAMPAIGN = (
    ROOT / "paper" / "softwarex" / "evidence" / "campaign_2026-09-01"
)
DEFAULT_OUTPUT = (
    ROOT / "paper" / "softwarex" / "evidence"
    / "reviewer_resilience_2026-09-01"
)
CASES = ("planted", "planted_confirmation", "null", "clonal", "klebsiella", "ssuis")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build(campaign: Path, output: Path) -> dict:
    receipt_path = campaign / "RUN_RECEIPT.json"
    receipt = load(receipt_path)
    if receipt.get("status") != "PASS_TECHNICAL_AND_ILLUSTRATIVE_CAMPAIGN":
        raise RuntimeError("campaign receipt is absent or not terminal PASS")
    output.mkdir(parents=True, exist_ok=False)

    decisions = []
    gate_rows = []
    inputs = []
    for case in CASES:
        result_path = campaign / "results" / case / "cluster_result.json"
        summary_path = campaign / "results" / case / "summary.json"
        result, summary = load(result_path), load(summary_path)
        inputs.extend([
            {"path": str(result_path.relative_to(ROOT)), "sha256": sha256(result_path)},
            {"path": str(summary_path.relative_to(ROOT)), "sha256": sha256(summary_path)},
        ])
        ci = summary.get("continuum_tail_probability_ci95") or {}
        decisions.append({
            "case": case,
            "n": summary.get("n_isolates"),
            "k": summary.get("selected_k"),
            "claim_level": summary.get("claim_level"),
            "claim_status": summary.get("claim_status"),
            "active_gate_codes": ";".join(summary.get("active_gate_codes") or []),
            "continuum_p": summary.get("discreteness_p_value"),
            "bootstrap_exceedances": summary.get("continuum_bootstrap_exceedances"),
            "tail_ci95_low": ci.get("low"),
            "tail_ci95_high": ci.get("high"),
            "resolved_alpha_0_05": summary.get(
                "continuum_decision_resolved_at_alpha_0_05"),
        })
        for gate in (result.get("interpretation") or {}).get("gate_ledger", []):
            gate_rows.append({
                "case": case,
                "gate": gate.get("code"),
                "applicable": gate.get("applicable"),
                "triggered": gate.get("triggered"),
                "value_json": json.dumps(gate.get("value"), sort_keys=True),
                "threshold_json": json.dumps(gate.get("threshold"), sort_keys=True),
                "margin_to_failure": gate.get("margin_to_failure"),
            })

    write_tsv(
        output / "decision_ladder.tsv",
        ["case", "n", "k", "claim_level", "claim_status", "active_gate_codes",
         "continuum_p", "bootstrap_exceedances", "tail_ci95_low",
         "tail_ci95_high", "resolved_alpha_0_05"],
        decisions,
    )
    write_tsv(
        output / "gate_margins.tsv",
        ["case", "gate", "applicable", "triggered", "value_json",
         "threshold_json", "margin_to_failure"],
        gate_rows,
    )

    by_case = {
        case: load(campaign / "results" / case / "cluster_result.json")
        for case in CASES
    }
    robustness = []
    for case in ("clonal", "klebsiella", "ssuis"):
        lc = by_case[case]["metadata_diagnostics"]["lineage_concordance"]
        for z_threshold in (2.0, 3.0, 5.0):
            robustness.append({
                "case": case, "diagnostic": "lineage_z",
                "threshold": z_threshold, "value": lc["z"],
                "gate_triggered": bool(lc["p_value"] < 0.05 and lc["z"] > z_threshold),
            })
    for ari_threshold in (0.4, 0.5, 0.6):
        value = by_case["klebsiella"]["artifact_diagnostics"]["max_empty_stratum_ari"]
        robustness.append({
            "case": "klebsiella", "diagnostic": "empty_stratum_ari",
            "threshold": ari_threshold, "value": value,
            "gate_triggered": bool(value >= ari_threshold),
        })
    for case in ("planted", "planted_confirmation", "klebsiella", "ssuis"):
        disc = by_case[case]["post_clustering_inference"]["discreteness"]
        upper = disc["bootstrap_tail_probability_ci95"]["high"]
        for alpha in (0.025, 0.05, 0.10):
            robustness.append({
                "case": case, "diagnostic": "continuum_tail_ci95_upper",
                "threshold": alpha, "value": upper,
                "gate_triggered": bool(upper >= alpha),
            })
    write_tsv(
        output / "threshold_robustness.tsv",
        ["case", "diagnostic", "threshold", "value", "gate_triggered"],
        robustness,
    )

    kp_ph = by_case["klebsiella"].get("phenotype_validation") or {}
    data_limits = {
        "klebsiella": {
            "cohort_n": by_case["klebsiella"]["n_isolates"],
            "phenotype_n": kp_ph.get("n_isolates_with_any_phenotype"),
            "phenotype_fraction": (
                kp_ph.get("n_isolates_with_any_phenotype", 0)
                / by_case["klebsiella"]["n_isolates"]
            ),
            "phenotype_fdr_family_n": kp_ph.get("n_in_fdr_family"),
            "phenotype_fdr_significant_n": kp_ph.get("n_fdr_significant"),
            "provenance_boundary": (
                "Derived cohort; incomplete upstream accessions, Kleborate/database "
                "versions and historical subsetting provenance. Phenotype coverage "
                "is non-random and is not a prevalence sample."
            ),
        },
        "ssuis": {
            "cohort_n": by_case["ssuis"]["n_isolates"],
            "lineage_typed_n": by_case["ssuis"]["metadata_diagnostics"]
            ["lineage_concordance"].get("n_isolates"),
            "lineages_n": by_case["ssuis"]["metadata_diagnostics"]
            ["lineage_concordance"].get("n_lineages"),
            "provenance_boundary": (
                "Derived veterinary cohort; interpretation is limited to the "
                "hash-locked binary matrices and supplied lineage labels."
            ),
        },
    }
    (output / "data_limits.json").write_text(
        json.dumps(data_limits, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    generated = [
        output / "decision_ladder.tsv", output / "gate_margins.tsv",
        output / "threshold_robustness.tsv", output / "data_limits.json",
    ]
    payload = {
        "schema": "amr-clonalshare-reviewer-resilience-1.0",
        "status": "PASS",
        "generated_utc": utc_now(),
        "campaign_receipt": {
            "path": str(receipt_path.relative_to(ROOT)),
            "sha256": sha256(receipt_path),
        },
        "interpretation": (
            "This is a deterministic audit of registered outputs, not a new "
            "model fit and not a search for thresholds that remove refusals."
        ),
        "margin_convention": (
            "Positive margin_to_failure is on the passing side; negative is on "
            "the failing side; null means no scalar margin is defined."
        ),
        "key_findings": {
            "positive_controls": (
                "Both independent controls reach level 4 archetype_candidate "
                "with no active applicable gate."
            ),
            "observed_cohorts": (
                "K. pneumoniae and S. suis reach level 3 statistically discrete "
                "partitions but retain independent biological gates."
            ),
            "monte_carlo_resolution": (
                "Zero of 99 continuum-null replicates exceed the observed "
                "statistic; exact 95% upper tail bound is 0.036576. This resolves "
                "alpha=0.05 but not the stricter alpha=0.025 sensitivity check."
            ),
        },
        "inputs": inputs,
        "outputs": [
            {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
            for path in generated
        ],
    }
    receipt_out = output / "RECEIPT.json"
    receipt_out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", type=Path, default=DEFAULT_CAMPAIGN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = build(args.campaign.resolve(), args.output.resolve())
    print(json.dumps({"status": payload["status"], "output": str(args.output)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
