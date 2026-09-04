#!/usr/bin/env python3
"""Run and receipt the SoftwareX 1.0.0 evidence campaign.

Each invocation runs one named stage into a new campaign directory. The
``finalize`` stage fails closed unless all required commands have terminal
records, expected exit codes, required artifacts, schema 1.0, and the declared
positive/negative/adversarial/real-cohort decisions.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CAMPAIGN = (
    ROOT / "paper" / "softwarex" / "evidence" / "campaign_2026-09-01"
)
REQUIRED_OUTPUTS = (
    "summary.json",
    "cluster_result.json",
    "archetype_profiles.tsv",
    "assignment.tsv",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path) -> dict:
    return {
        "path": path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=ROOT, check=True, capture_output=True, text=True
    ).stdout.strip()


def campaign_inputs() -> list[Path]:
    paths = [
        ROOT / "pyproject.toml",
        ROOT / "uv.lock",
        ROOT / "CITATION.cff",
        ROOT / "LICENSE",
        ROOT / "Licence.txt",
    ]
    for directory in (ROOT / "src", ROOT / "tests"):
        paths.extend(p for p in directory.rglob("*.py") if p.is_file())
    for directory in (
        ROOT / "examples" / "synthetic" / "data",
        ROOT / "examples" / "klebsiella" / "data",
        ROOT / "examples" / "ssuis" / "data",
    ):
        paths.extend(p for p in directory.rglob("*") if p.is_file())
    paths.extend(
        [
            ROOT / "examples" / "synthetic" / "planted.yaml",
            ROOT / "examples" / "synthetic" / "planted_confirmation.yaml",
            ROOT / "examples" / "synthetic" / "null.yaml",
            ROOT / "examples" / "synthetic" / "clonal.yaml",
            ROOT / "examples" / "synthetic" / "make_data.py",
            ROOT / "examples" / "klebsiella" / "config.yaml",
            ROOT / "examples" / "ssuis" / "config.yaml",
        ]
    )
    return sorted(set(paths))


def command_specs(campaign: Path) -> dict[str, dict]:
    py = sys.executable
    result_root = campaign / "results"
    return {
        "pytest": {
            "command": [
                py,
                "-m",
                "pytest",
                "-q",
                "--junitxml",
                str(campaign / "pytest-junit.xml"),
            ],
            "allowed_exit_codes": [0],
        },
        "make_synthetic": {
            "command": [py, "examples/synthetic/make_data.py"],
            "allowed_exit_codes": [0],
        },
        "validate_planted": {
            "command": [
                py,
                "-m",
                "amr_clonalshare.cli",
                "--config",
                "examples/synthetic/planted.yaml",
                "--validate",
                "--seed",
                "42",
            ],
            "allowed_exit_codes": [0],
        },
        "planted": {
            "command": [
                py,
                "-m",
                "amr_clonalshare.cli",
                "--config",
                "examples/synthetic/planted.yaml",
                "--results-dir",
                str(result_root / "planted"),
                "--seed",
                "42",
                "--threads",
                "1",
            ],
            "allowed_exit_codes": [0],
        },
        "planted_confirmation": {
            "command": [
                py,
                "-m",
                "amr_clonalshare.cli",
                "--config",
                "examples/synthetic/planted_confirmation.yaml",
                "--results-dir",
                str(result_root / "planted_confirmation"),
                "--seed",
                "42",
                "--threads",
                "1",
            ],
            "allowed_exit_codes": [0],
        },
        "null": {
            "command": [
                py,
                "-m",
                "amr_clonalshare.cli",
                "--config",
                "examples/synthetic/null.yaml",
                "--results-dir",
                str(result_root / "null"),
                "--seed",
                "42",
                "--threads",
                "1",
            ],
            "allowed_exit_codes": [0],
        },
        "clonal": {
            "command": [
                py,
                "-m",
                "amr_clonalshare.cli",
                "--config",
                "examples/synthetic/clonal.yaml",
                "--results-dir",
                str(result_root / "clonal"),
                "--seed",
                "42",
                "--threads",
                "1",
            ],
            "allowed_exit_codes": [3],
        },
        "klebsiella": {
            "command": [
                py,
                "-m",
                "amr_clonalshare.cli",
                "--config",
                "examples/klebsiella/config.yaml",
                "--results-dir",
                str(result_root / "klebsiella"),
                "--seed",
                "42",
                "--threads",
                "1",
            ],
            "allowed_exit_codes": [3],
        },
        "ssuis": {
            "command": [
                py,
                "-m",
                "amr_clonalshare.cli",
                "--config",
                "examples/ssuis/config.yaml",
                "--results-dir",
                str(result_root / "ssuis"),
                "--seed",
                "42",
                "--threads",
                "1",
            ],
            # A permutation gate on same-lineage pairs would fire on this
            # cohort as it would on any cohort of 677 isolates, and the run
            # would exit 3. Under the attribution gate the lineage question
            # resolves below threshold and the run completes clean, so the
            # expected code is 0 and the campaign fails closed if it is not.
            "allowed_exit_codes": [0],
        },
        # The prevalence decomposition, run four ways. The two hierBAPS runs
        # are the result; the two MLST runs are the demonstration that the
        # estimability gate fires, and they are part of the campaign for the
        # same reason the null control is: a gate that is never seen to fire
        # is not evidence that it would.
        "decompose_period": {
            "command": [py, "examples/ssuis/decompose_trend.py",
                        "--contrast", "period", "--lineage", "baps_cluster",
                        "--n-boot", "4000", "--seed", "23",
                        "--json", str(result_root / "decompose_period.json")],
            "allowed_exit_codes": [0],
        },
        "decompose_country": {
            "command": [py, "examples/ssuis/decompose_trend.py",
                        "--contrast", "country", "--lineage", "baps_cluster",
                        "--n-boot", "4000", "--seed", "23",
                        "--json", str(result_root / "decompose_country.json")],
            "allowed_exit_codes": [0],
        },
        "decompose_period_mlst": {
            "command": [py, "examples/ssuis/decompose_trend.py",
                        "--contrast", "period", "--lineage", "mlst",
                        "--n-boot", "4000", "--seed", "23",
                        "--json", str(result_root / "decompose_period_mlst.json")],
            "allowed_exit_codes": [0],
        },
        "decompose_country_mlst": {
            "command": [py, "examples/ssuis/decompose_trend.py",
                        "--contrast", "country", "--lineage", "mlst",
                        "--n-boot", "4000", "--seed", "23",
                        "--json", str(result_root / "decompose_country_mlst.json")],
            "allowed_exit_codes": [0],
        },
        # The sensitivity analysis on the three agents with a published EUCAST
        # cut-off. It is a campaign stage rather than a note because the
        # manuscript states what it found, including that one of the two
        # conclusions cannot be tested this way.
        "published_cutoff_country": {
            "command": [py, "examples/ssuis/decompose_trend.py",
                        "--contrast", "country", "--agents", "published_cutoff",
                        "--n-boot", "4000", "--seed", "23",
                        "--json", str(result_root / "published_cutoff_country.json")],
            "allowed_exit_codes": [0],
        },
        "published_cutoff_period": {
            "command": [py, "examples/ssuis/decompose_trend.py",
                        "--contrast", "period", "--agents", "published_cutoff",
                        "--n-boot", "4000", "--seed", "23",
                        "--json", str(result_root / "published_cutoff_period.json")],
            "allowed_exit_codes": [0],
        },
        "genotype_check": {
            "command": [py, "examples/ssuis/validate_with_determinants.py",
                        "--contrast", "country", "--n-boot", "4000",
                        "--seed", "7",
                        "--json", str(result_root / "genotype_check.json")],
            "allowed_exit_codes": [0],
        },
        "provenance": {
            "command": [py, "examples/ssuis/link_source_lineages.py", "--check"],
            "allowed_exit_codes": [0],
        },
        "metadata_repair": {
            "command": [py, "examples/ssuis/repair_metadata.py", "--check"],
            "allowed_exit_codes": [0],
        },
    }


def initialize(campaign: Path) -> None:
    campaign.mkdir(parents=True, exist_ok=True)
    (campaign / "logs").mkdir(exist_ok=True)
    (campaign / "results").mkdir(exist_ok=True)
    dependency_names = (
        "numpy",
        "pandas",
        "scipy",
        "scikit-learn",
        "pyyaml",
        "pandera",
        "pytest",
        "dendropy",
    )
    environment = {
        "created_utc": utc_now(),
        "campaign": "amr-clonalshare-1.0.0-softwarex",
        "product_version": "1.0.0",
        "python_executable": sys.executable,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "git_head": git("rev-parse", "HEAD"),
        "git_branch": git("branch", "--show-current"),
        "git_status_scoped": git("status", "--short", "--", ".").splitlines(),
        "dependencies": {},
    }
    for name in dependency_names:
        try:
            environment["dependencies"][name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            environment["dependencies"][name] = None
    (campaign / "environment.json").write_text(
        json.dumps(environment, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    inputs = [file_record(path) for path in campaign_inputs()]
    (campaign / "input_manifest.json").write_text(
        json.dumps(inputs, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    records = campaign / "command_records.json"
    if not records.exists():
        records.write_text("{}\n", encoding="utf-8")


def run_stage(campaign: Path, stage: str) -> int:
    initialize(campaign)
    specs = command_specs(campaign)
    if stage not in specs:
        raise SystemExit(f"unknown stage: {stage}")
    records_path = campaign / "command_records.json"
    records = json.loads(records_path.read_text(encoding="utf-8"))
    spec = specs[stage]
    command = spec["command"]
    start_utc = utc_now()
    start = time.monotonic()
    print(f"START {stage} {start_utc}", flush=True)
    proc = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    end_utc = utc_now()
    stdout_path = campaign / "logs" / f"{stage}.stdout.log"
    stderr_path = campaign / "logs" / f"{stage}.stderr.log"
    stdout_path.write_text(proc.stdout, encoding="utf-8")
    stderr_path.write_text(proc.stderr, encoding="utf-8")
    output_dir = campaign / "results" / stage
    outputs = []
    if output_dir.exists():
        outputs = [file_record(path) for path in sorted(output_dir.rglob("*")) if path.is_file()]
    records[stage] = {
        "stage": stage,
        "command": command,
        "command_display": " ".join(command),
        "cwd": str(ROOT),
        "start_utc": start_utc,
        "end_utc": end_utc,
        "elapsed_seconds": round(time.monotonic() - start, 3),
        "exit_code": proc.returncode,
        "allowed_exit_codes": spec["allowed_exit_codes"],
        "terminal": True,
        "stdout": file_record(stdout_path),
        "stderr": file_record(stderr_path),
        "outputs": outputs,
    }
    records_path.write_text(
        json.dumps(records, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    ok = proc.returncode in spec["allowed_exit_codes"]
    print(
        f"END {stage} {end_utc} exit={proc.returncode} "
        f"elapsed={records[stage]['elapsed_seconds']}s expected={ok}",
        flush=True,
    )
    return 0 if ok else 1


def load_summary(campaign: Path, case: str) -> dict:
    return json.loads(
        (campaign / "results" / case / "summary.json").read_text(encoding="utf-8")
    )


def load_gate(campaign: Path, case: str):
    """The lineage-attribution threshold the run itself applied.

    Read back from the gate ledger in cluster_result.json rather than written
    into this script, so a check on the margin cannot drift from the threshold
    the run was actually judged against.
    """
    path = campaign / "results" / case / "cluster_result.json"
    if not path.is_file():
        return None
    result = json.loads(path.read_text(encoding="utf-8"))
    for record in (result.get("interpretation") or {}).get("gate_ledger", []):
        if record.get("code") == "lineage_confounding":
            return (record.get("threshold") or {}).get(
                "lineage_attributable_share_below")
    return None


def acceptance_checks(campaign: Path, records: dict) -> list[dict]:
    """Every declared decision the campaign must reproduce, with its evidence.

    The two real-cohort checks read the attributable share and its margin to
    the gate rather than whether a bootstrap interval resolved. The partition
    is now rebuilt inside every cross-validation fold, and under that refit the
    lineage bootstrap cannot be run: it draws lineages whole, so a lineage
    drawn twice puts the same isolates on both sides of a fold, and a partition
    rebuilt on that cohort trains on rows it is also scored against. The run
    therefore reports the spread of the share across repeated fold draws as
    lam_cv_sd and leaves the interval empty, and the gate has always been
    decided by the point estimate, so the margin is the quantity to assert.
    """
    checks: list[dict] = []

    def add(check_id: str, passed: bool, evidence: object) -> None:
        checks.append({"id": check_id, "passed": bool(passed), "evidence": evidence})

    required_stages = tuple(command_specs(campaign))
    add(
        "all_stages_terminal",
        all(records.get(name, {}).get("terminal") for name in required_stages),
        required_stages,
    )
    add(
        "all_exit_codes_expected",
        all(
            records.get(name, {}).get("exit_code")
            in records.get(name, {}).get("allowed_exit_codes", [])
            for name in required_stages
        ),
        {name: records.get(name, {}).get("exit_code") for name in required_stages},
    )
    junit = campaign / "pytest-junit.xml"
    add("terminal_junit_exists", junit.is_file() and junit.stat().st_size > 0, str(junit))

    for case in (
        "planted",
        "planted_confirmation",
        "null",
        "clonal",
        "klebsiella",
        "ssuis",
    ):
        directory = campaign / "results" / case
        missing = [name for name in REQUIRED_OUTPUTS if not (directory / name).is_file()]
        add(f"{case}_four_public_outputs", not missing, missing)
        result_path = directory / "cluster_result.json"
        if result_path.is_file():
            result = json.loads(result_path.read_text(encoding="utf-8"))
            add(f"{case}_schema_1_0", result.get("schema_version") == "1.0", result.get("schema_version"))

    for case in ("planted", "planted_confirmation"):
        result_path = campaign / "results" / case / "cluster_result.json"
        if result_path.is_file():
            summary = load_summary(campaign, case)
            result = json.loads(result_path.read_text(encoding="utf-8"))
            external = (
                (result.get("metadata_diagnostics") or {})
                .get("external_agreement", {})
                .get("planted_cluster", {})
            )
            evidence = {
                "selected_k": summary.get("selected_k"),
                "ari_vs_truth": external.get("ari"),
                "inference_status": summary.get("inference_status"),
                "split_design_adequate": summary.get("split_design_adequate"),
                "discrete_beyond_a_gradient": summary.get(
                    "discrete_beyond_a_gradient"
                ),
                "lineage_confounded": summary.get("lineage_confounded"),
                "claim_status": summary.get("claim_status"),
                "active_gate_codes": summary.get("active_gate_codes"),
                "continuum_bootstrap_exceedances": summary.get(
                    "continuum_bootstrap_exceedances"),
                "continuum_tail_probability_ci95": summary.get(
                    "continuum_tail_probability_ci95"),
                "exit_code": records.get(case, {}).get("exit_code"),
            }
            add(
                f"{case}_positive_archetype_acceptance",
                summary.get("selected_k") == 3
                and (external.get("ari") or 0.0) >= 0.95
                and summary.get("inference_status") == "ok"
                and summary.get("split_design_adequate") is True
                and summary.get("discrete_beyond_a_gradient") is True
                and summary.get("lineage_confounded") is False
                and summary.get("claim_status") == "archetype_candidate"
                and summary.get("active_gate_codes") == []
                and summary.get("continuum_bootstrap_exceedances") == 0
                and (summary.get("continuum_tail_probability_ci95") or {}).get(
                    "high", 1.0) < 0.05
                and records.get(case, {}).get("exit_code") == 0,
                evidence,
            )
    if (campaign / "results" / "null" / "summary.json").is_file():
        null = load_summary(campaign, "null")
        add(
            "null_refuses_clusters",
            null.get("selected_k") == 1 and null.get("no_structure") is True,
            {"selected_k": null.get("selected_k"), "no_structure": null.get("no_structure")},
        )
    if (campaign / "results" / "clonal" / "summary.json").is_file():
        clonal = load_summary(campaign, "clonal")
        add("clonal_lineage_gate", clonal.get("lineage_confounded") is True, clonal.get("lineage_confounded"))
    if (campaign / "results" / "klebsiella" / "summary.json").is_file():
        kp = load_summary(campaign, "klebsiella")
        kp_share = kp.get("lineage_attributable_share")
        kp_gate = load_gate(campaign, "klebsiella")
        add(
            "klebsiella_discrete_signal_but_other_gates_block",
            kp.get("discrete_beyond_a_gradient") is True
            and kp.get("continuum_null_under_dimensioned") is False
            and kp.get("claim_status") == "statistically_discrete_partition"
            # `lineage_confounding` is deliberately not required here any more.
            # Under 1.0 it was active on this cohort because the permutation
            # test fires on any 1500 isolates, not because the lineage label
            # was the cohort's problem. Its actual problems are the two gates
            # below, and the attribution is asserted separately.
            and {"empty_stratum", "phenotype_superiority"}
            .issubset(set(kp.get("active_gate_codes") or []))
            and kp_share is not None
            and kp_gate is not None
            and kp_share < kp_gate
            and kp.get("continuum_bootstrap_exceedances") == 0
            and (kp.get("continuum_tail_probability_ci95") or {}).get(
                "high", 1.0) < 0.05
            and records.get("klebsiella", {}).get("exit_code") == 3,
            {
                "discrete_beyond_a_gradient": kp.get("discrete_beyond_a_gradient"),
                "continuum_null_under_dimensioned": kp.get(
                    "continuum_null_under_dimensioned"
                ),
                "claim_status": kp.get("claim_status"),
                "active_gate_codes": kp.get("active_gate_codes"),
                "lineage_attributable_share": kp_share,
                "lambda_gate": kp_gate,
                "margin_to_gate": (
                    None if kp_share is None or kp_gate is None
                    else kp_gate - kp_share),
                "continuum_tail_probability_ci95": kp.get(
                    "continuum_tail_probability_ci95"),
                "exit_code": records.get("klebsiella", {}).get("exit_code"),
            },
        )
    if (campaign / "results" / "ssuis" / "summary.json").is_file():
        ssuis = load_summary(campaign, "ssuis")
        share = ssuis.get("lineage_attributable_share")
        ci = ssuis.get("lineage_attributable_share_ci95") or [None, None]
        gate = load_gate(campaign, "ssuis")
        # Under 1.0 this check asserted that the veterinary cohort was refused.
        # That was not a scientific expectation, it was the arithmetic of a
        # permutation test over 229,000 pairs, and it was true of every cohort
        # large enough to analyse. What is asserted now is that the lineage
        # question resolves: the attributable share sits below the threshold
        # the run recorded by a margin of at least 0.05, so the gate is not
        # active and the partition is released. The concordance z is recorded
        # beside it and is expected to be large, which is the point.
        add(
            "ssuis_lineage_attribution_resolves_below_gate",
            share is not None
            and gate is not None
            and share < gate
            and gate - share >= 0.05
            and ssuis.get("lineage_confounded") is False
            and "lineage_confounding" not in (ssuis.get("active_gate_codes") or [])
            and ssuis.get("continuum_bootstrap_exceedances") == 0
            and (ssuis.get("continuum_tail_probability_ci95") or {}).get(
                "high", 1.0) < 0.05
            and records.get("ssuis", {}).get("exit_code") == 0,
            {
                "lineage_attributable_share": share,
                "lineage_attributable_share_ci95": ci,
                "lambda_gate": gate,
                "margin_to_gate": (
                    None if share is None or gate is None else gate - share),
                "lineage_concordance_z_not_gated": ssuis.get(
                    "lineage_concordance_z"),
                "clonal_share_all_features": ssuis.get(
                    "clonal_share_all_features"),
                "lineage_confounded": ssuis.get("lineage_confounded"),
                "claim_status": ssuis.get("claim_status"),
                "active_gate_codes": ssuis.get("active_gate_codes"),
                "continuum_tail_probability_ci95": ssuis.get(
                    "continuum_tail_probability_ci95"),
                "exit_code": records.get("ssuis", {}).get("exit_code"),
            },
        )
    decomposition_checks(campaign, add)
    add("licence_files_identical", (ROOT / "LICENSE").read_bytes() == (ROOT / "Licence.txt").read_bytes(), None)
    return checks


def _panel(campaign: Path, name: str) -> dict:
    path = campaign / "results" / f"{name}.json"
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}


def decomposition_checks(campaign: Path, add) -> None:
    """The decomposition must give both answers and refuse in both bad cases.

    A method that only ever returns "composition" would be indistinguishable
    from one that measures lineage turnover and calls it a decomposition, so
    the campaign requires one contrast that loads on composition and one that
    loads on the within-lineage rate. It requires the same panel to be refused
    under the sequence-type labels, which is what shows the gate is a gate
    rather than a field that is always true.
    """
    period, country = _panel(campaign, "decompose_period"), _panel(campaign, "decompose_country")
    if period and country:
        pf, cf = period["family"], country["family"]
        add("decomposition_period_loads_on_composition",
            pf["n_composition_discoveries"] >= 1
            and pf["n_within_lineage_discoveries"] == 0,
            {"composition": pf["n_composition_discoveries"],
             "within_lineage": pf["n_within_lineage_discoveries"]})
        add("decomposition_country_loads_on_within_lineage",
            cf["n_within_lineage_discoveries"] >= 1
            and cf["n_composition_discoveries"] == 0,
            {"composition": cf["n_composition_discoveries"],
             "within_lineage": cf["n_within_lineage_discoveries"]})
        residual = max(abs(r["identity_residual"])
                       for panel in (period, country)
                       for r in panel["per_agent"].values()
                       if r.get("status") == "ok")
        add("decomposition_identity_exact", residual < 1e-12,
            {"largest_absolute_residual": residual})
        add("decomposition_estimable_under_population_clusters",
            all(r["within_lineage_estimable"]
                for panel in (period, country)
                for r in panel["per_agent"].values() if r.get("status") == "ok"),
            None)
    for name, reason in (("decompose_period_mlst", "differentially missing"),
                         ("decompose_country_mlst", "support")):
        panel = _panel(campaign, name)
        if not panel:
            continue
        records = [r for r in panel["per_agent"].values() if r.get("status") == "ok"]
        add(f"gate_fires_for_{name}",
            bool(records) and not any(r["within_lineage_estimable"] for r in records)
            and any(reason in cause for r in records
                    for cause in r["not_estimable_because"]),
            {"reasons": sorted({cause for r in records
                                for cause in r["not_estimable_because"]})})
    anchored = _panel(campaign, "published_cutoff_country")
    if anchored:
        add("country_conclusion_survives_published_cutoffs_only",
            anchored["family"]["n_within_lineage_discoveries"]
            == anchored["family"]["n_agents_decomposed"],
            {"agents": anchored["family"]["n_agents_decomposed"],
             "within_lineage_discoveries":
                 anchored["family"]["n_within_lineage_discoveries"]})
    genotype = _panel(campaign, "genotype_check")
    if genotype:
        agreement = genotype["agreement"]["within_lineage_correlation"]
        add("genotype_layer_reproduces_the_split", agreement > 0.8,
            {"within_lineage_correlation": agreement})


def finalize(campaign: Path) -> int:
    initialize(campaign)
    records = json.loads((campaign / "command_records.json").read_text(encoding="utf-8"))
    checks = acceptance_checks(campaign, records)
    passed = all(check["passed"] for check in checks)
    output_files = [
        file_record(path)
        for path in sorted(campaign.rglob("*"))
        if path.is_file() and path.name not in {"RUN_RECEIPT.json", "SHA256SUMS"}
    ]
    receipt = {
        "schema": "amr-clonalshare-softwarex-campaign-receipt-1.0",
        "status": "PASS_TECHNICAL_AND_ILLUSTRATIVE_CAMPAIGN" if passed else "FAIL_CLOSED",
        "product_version": "1.0.0",
        "finalized_utc": utc_now(),
        "git_head": git("rev-parse", "HEAD"),
        "scientific_boundary": (
            "A passed receipt confirms execution, artifact integrity and declared control/gate "
            "behaviour. It does not authorize discrete archetype claims when a gate fires."
        ),
        "checks": checks,
        "commands": records,
        "files": output_files,
    }
    receipt_path = campaign / "RUN_RECEIPT.json"
    receipt_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest_files = [p for p in sorted(campaign.rglob("*")) if p.is_file() and p.name != "SHA256SUMS"]
    lines = [f"{sha256(path)}  {path.relative_to(campaign).as_posix()}" for path in manifest_files]
    (campaign / "SHA256SUMS").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"FINAL {receipt['status']} checks={len(checks)} receipt={receipt_path}", flush=True)
    return 0 if passed else 1


def sync_fixtures(campaign: Path) -> int:
    mapping = {
        "planted": ROOT / "examples" / "synthetic" / "expected_planted",
        "planted_confirmation": (
            ROOT / "examples" / "synthetic" / "expected_planted_confirmation"
        ),
        "null": ROOT / "examples" / "synthetic" / "expected_null",
        "clonal": ROOT / "examples" / "synthetic" / "expected_clonal",
        "klebsiella": ROOT / "examples" / "klebsiella" / "expected",
    }
    for case, destination in mapping.items():
        source = campaign / "results" / case
        missing = [name for name in REQUIRED_OUTPUTS if not (source / name).is_file()]
        if missing:
            raise SystemExit(f"cannot sync {case}; missing {missing}")
        destination.mkdir(parents=True, exist_ok=True)
        for name in REQUIRED_OUTPUTS:
            shutil.copy2(source / name, destination / name)
    print("SYNCED fresh campaign outputs into tracked example fixtures", flush=True)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "stage",
        choices=[*command_specs(DEFAULT_CAMPAIGN).keys(), "finalize", "sync-fixtures"],
    )
    parser.add_argument("--campaign", type=Path, default=DEFAULT_CAMPAIGN)
    args = parser.parse_args()
    campaign = args.campaign.resolve()
    if args.stage == "finalize":
        return finalize(campaign)
    if args.stage == "sync-fixtures":
        return sync_fixtures(campaign)
    return run_stage(campaign, args.stage)


if __name__ == "__main__":
    raise SystemExit(main())
