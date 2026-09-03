import pytest

from amr_clonalshare.cli import (
    _apply_release_gates,
    _interpretation_contract,
    _summary,
)


def test_inadequate_split_design_is_fail_closed_in_public_artifacts():
    result = {
        "post_clustering_inference": {
            "status": "ok",
            "split_design_adequate": False,
            "split_design_problems": ["only 3 split units"],
            "p_value_structure": 1e-30,
            "structure_detected": True,
        }
    }

    _apply_release_gates(result)
    infer = result["post_clustering_inference"]
    public = _summary(result)

    assert infer["status"] == "withheld_inadequate_split_design"
    assert infer["p_value_structure"] is None
    assert infer["structure_detected"] is False
    assert infer["diagnostic_p_value_structure"] == 1e-30
    assert public["p_value_structure_report"].startswith("withheld")
    assert public["diagnostic_p_value_structure_exact"] == 1e-30


def _candidate_result(*, lineage_z=0.0, empty_warning=False):
    return {
        "selected_k": 3,
        "k_selection": {"no_structure": False},
        "layer_influence": {"collapse": False, "n_eff": 2.0,
                            "collapse_threshold": 1.5},
        "artifact_diagnostics": {
            "empty_stratum_warning": empty_warning,
            "max_empty_stratum_ari": 0.6 if empty_warning else 0.1,
        },
        "metadata_diagnostics": {"lineage_concordance": {
            "status": "ok", "p_value": 0.001 if lineage_z > 2 else 0.8,
            "z": lineage_z, "concordance_observed": 0.5,
            "concordance_null_mean": 0.5,
        }},
        "post_clustering_inference": {
            "status": "ok", "structure_detected": True,
            "split_design_adequate": True,
            "discreteness": {
                "status": "ok", "discrete_beyond_a_gradient": True,
                "p_value": 0.01,
                "bootstrap_tail_probability_ci95": {"low": 0.0, "high": 0.037},
            },
        },
    }


def test_claim_ladder_accepts_only_a_gate_clean_candidate():
    result = _candidate_result()
    result["interpretation"] = _interpretation_contract(result)
    public = _summary(result)

    assert public["claim_level"] == 4
    assert public["claim_status"] == "archetype_candidate"
    assert public["active_gate_codes"] == []


def test_claim_ladder_preserves_discreteness_but_blocks_archetype_promotion():
    result = _candidate_result(lineage_z=8.0, empty_warning=True)
    interpretation = _interpretation_contract(result)

    assert interpretation["claim_level"] == 3
    assert interpretation["claim_status"] == "statistically_discrete_partition"
    assert interpretation["active_gate_codes"] == [
        "empty_stratum", "lineage_confounding"]
    ledger = {row["code"]: row for row in interpretation["gate_ledger"]}
    assert ledger["empty_stratum"]["margin_to_failure"] == pytest.approx(-0.1)
    assert ledger["lineage_confounding"]["margin_to_failure"] == pytest.approx(-6.0)
