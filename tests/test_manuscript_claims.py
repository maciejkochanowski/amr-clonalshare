"""Guard the current evidence hierarchy while preserving the frozen long paper.

``paper/manuscript.md`` is a historical methodology manuscript and must remain
byte-identical. Current numerical assertions come from the regenerated 1.0.0
fixtures: discreteness is estimable and detected in the real cohort, while
independent artifact, lineage and phenotype gates still prevent biological
archetype interpretation.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
#: The retired methodology manuscript belongs to the article bundle, not to the
#: software. Where it is present its digest is pinned here as well as in
#: `paper/softwarex/scripts/qa_manuscript.py`; where it is not, the two checks
#: that read it skip, so that a distribution of the software carries no part of
#: the article.
_FROZEN = ROOT / "paper" / "manuscript.md"
MANUSCRIPT = _FROZEN.read_text(encoding="utf-8") if _FROZEN.is_file() else None
needs_article_bundle = pytest.mark.skipif(
    MANUSCRIPT is None,
    reason="paper/manuscript.md ships in the article bundle, not here")
README = (ROOT / "README.md").read_text(encoding="utf-8")
README_FLAT = " ".join(README.split())
SUMMARY = json.loads(
    (ROOT / "examples" / "klebsiella" / "expected" / "summary.json").read_text()
)
RESULT = json.loads(
    (ROOT / "examples" / "klebsiella" / "expected" / "cluster_result.json").read_text()
)
#: The retired methodological manuscript, pinned so that work on the SoftwareX
#: article cannot alter it silently. Re-pinned once, on 2026-08-31, when five
#: Polish placeholder tokens and one place name were rendered in English; the
#: previous digest was
#: 577650c30a5d731c5e2db311bcd2f560b631f608a0ac5daa373ec239f3f0e9f6 and the
#: change is recorded in paper/README.md. No sentence of its content changed.
FROZEN_MANUSCRIPT_SHA256 = (
    "f35f6d06767191b62918758d6105a8a9e561e2c39b83fcdb2a4f084323e4f7b0"
)


@needs_article_bundle
def test_frozen_methodology_manuscript_is_unchanged():
    assert hashlib.sha256(MANUSCRIPT.encode("utf-8")).hexdigest() == (
        FROZEN_MANUSCRIPT_SHA256
    )
    assert MANUSCRIPT.startswith(
        '---\ntitle: "A resistance-load gradient, not discrete archetypes:'
    )


def test_current_narrative_matches_the_shipped_discreteness_result():
    assert SUMMARY["structure_detected"] is True
    assert SUMMARY["discrete_beyond_a_gradient"] is True
    assert SUMMARY["discreteness_status"] == "ok"
    assert SUMMARY["discreteness_verdict"] == "discrete"
    assert SUMMARY["continuum_null_under_dimensioned"] is False
    assert SUMMARY["claim_level"] == 3
    assert SUMMARY["claim_status"] == "statistically_discrete_partition"
    assert set(SUMMARY["active_gate_codes"]) == {
        "empty_stratum", "lineage_confounding", "phenotype_superiority"
    }
    assert SUMMARY["continuum_bootstrap_exceedances"] == 0
    assert SUMMARY["continuum_tail_probability_ci95"]["high"] < 0.05
    assert SUMMARY["continuum_decision_resolved_at_alpha_0_05"] is True
    assert "discrete signal is detected" in README_FLAT
    assert "does not establish biological archetypes" in README_FLAT
    assert "claim level 4" in README_FLAT
    assert "0-0.0366" in README_FLAT


def test_public_p_value_reporting_is_thresholded_and_order_stable():
    assert SUMMARY["p_value_structure_report"] == (
        "p < 1e-100 across all checked split orderings"
    )
    assert SUMMARY["diagnostic_p_value_structure_exact"] < 1e-100
    assert "p_value_structure" not in SUMMARY


@needs_article_bundle
def test_the_exact_p_value_is_not_headlined_in_the_frozen_manuscript():
    assert "1.4 × 10⁻¹⁹³" not in MANUSCRIPT


def test_gradient_proxy_numbers_are_derived_from_the_result_artifact():
    ph = RESULT["phenotype_validation"]
    ba = ph["mean_balanced_accuracy_over_fdr_family"]
    boot = ph["head_to_head_bootstrap"]["external_resistance_score_above_median"]

    assert (
        f"{ba['external_resistance_score_above_median']:.3f} versus "
        f"{ba['partition']:.3f}" in README
    )
    assert boot["significantly_better_than_partition"] is True
    assert f"{boot['mean_difference_minus_partition']:+.3f}" in README
    assert f"[{boot['ci95_low']:+.3f}, {boot['ci95_high']:+.3f}]" in README
