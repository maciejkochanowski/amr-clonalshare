"""Relations between two runs, for the results no oracle can state.

Nobody can say what the clonal share of a fictional table ought to be, so a
test cannot compare it with a known answer. It can compare two runs: change
the input in a way whose effect on the answer is known in advance, and check
that the answer changed in exactly that way, or did not change at all. A
violated relation is a defect even when neither number can be called wrong on
its own.
"""
from __future__ import annotations

import csv
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

from amr_clonalshare import clonality, comparison, evalues, missingness
from amr_clonalshare.attribution import clonal_share
from amr_clonalshare.realised import realised_share
from amr_clonalshare.stats import benjamini_hochberg

RELATIVE = dict(rel=1e-9, abs=1e-12)


def _cohort(seed=5, n_lineages=9, per=6, prevalence=0.4):
    rng = np.random.default_rng(seed)
    lineage = np.repeat([f"L{i:02d}" for i in range(n_lineages)], per)
    y = (rng.random(n_lineages * per) < prevalence).astype(float)
    return y, lineage


# --------------------------------------------------------------------------- #
# Relations of the estimators
# --------------------------------------------------------------------------- #

def test_calling_the_other_outcome_positive_leaves_the_variance_ratio_alone():
    """y -> 1-y is an affine map; a ratio of mean squares cannot see it."""
    y, lineage = _cohort()
    here, flipped = realised_share(y, lineage), realised_share(1.0 - y, lineage)
    assert flipped.f_ratio == pytest.approx(here.f_ratio, **RELATIVE)
    assert flipped.kappa == pytest.approx(here.kappa, **RELATIVE)
    assert flipped.ci_low == pytest.approx(here.ci_low, **RELATIVE)
    assert flipped.ci_high == pytest.approx(here.ci_high, **RELATIVE)


@pytest.mark.parametrize("scale, shift", [(2.0, 0.0), (1.0, 7.5), (-3.0, 1.25),
                                          (0.125, -40.0)])
def test_a_change_of_unit_or_origin_leaves_the_variance_ratio_alone(scale, shift):
    """A dilution scale read in log2 or in log10, or from another zero."""
    rng = np.random.default_rng(11)
    lineage = np.repeat([f"L{i}" for i in range(7)], 8)
    y = rng.normal(size=lineage.size) + np.repeat(rng.normal(size=7), 8)
    here = realised_share(y, lineage)
    moved = realised_share(scale * y + shift, lineage)
    assert moved.f_ratio == pytest.approx(here.f_ratio, **RELATIVE)
    assert moved.kappa == pytest.approx(here.kappa, **RELATIVE)
    assert moved.ci_high == pytest.approx(here.ci_high, **RELATIVE)


def test_calling_the_other_outcome_positive_leaves_the_evidence_alone():
    """The likelihood ratio of a Bernoulli fit is the same under 1-y."""
    y, lineage = _cohort(seed=9)
    kw = dict(folds=3, repeats=3, seed=7)
    assert evalues.e_process(1.0 - y, lineage, **kw).e_value == pytest.approx(
        evalues.e_process(y, lineage, **kw).e_value, **RELATIVE)


def test_calling_the_other_outcome_positive_leaves_the_attributed_share_alone():
    """Predicting 1-y with 1-p is the same Brier score, and the same skill."""
    y, lineage = _cohort(seed=13)
    kw = dict(folds=3, repeats=3, n_boot=40, n_perm=39, seed=4)
    here, flipped = clonal_share(y, lineage, **kw), clonal_share(1.0 - y, lineage, **kw)
    assert flipped.kappa_adj == pytest.approx(here.kappa_adj, **RELATIVE)
    assert flipped.ci_low == pytest.approx(here.ci_low, **RELATIVE)
    assert flipped.ci_high == pytest.approx(here.ci_high, **RELATIVE)
    assert flipped.p_value == here.p_value


def test_exchanging_the_two_collections_turns_the_decomposition_around():
    """Every component of A - B is the negative of the component of B - A."""
    ya, la = [1, 0, 1, 1, 0, 0, 1, 0], list("AABBCCDD")
    yb, lb = [0, 0, 1, 0, 1, 1, 1, 1], list("AABBCCDD")
    here = clonality.decompose_prevalence_difference(ya, la, yb, lb, n_boot=0)
    back = clonality.decompose_prevalence_difference(yb, lb, ya, la, n_boot=0)
    for key in ("difference", "composition", "within_lineage", "typed_difference",
                "observed_difference"):
        assert back[key] == pytest.approx(-here[key], **RELATIVE)
    assert back["identity_residual"] == pytest.approx(here["identity_residual"],
                                                      **RELATIVE)


def test_exchanging_the_two_definitions_turns_the_comparison_around():
    y, lineage = _cohort(seed=21)
    other = np.array([f"C{i % 5}" for i in range(lineage.size)], dtype=object)
    ids = [f"iso_{i:03d}" for i in range(lineage.size)]
    kw = dict(ids=ids, seed=3, folds=3, repeats=3, n_perm=19)
    here = comparison.compare_lineage_definitions(y, lineage, other, name_a="v1",
                                                  name_b="v2", **kw)
    back = comparison.compare_lineage_definitions(y, other, lineage, name_a="v2",
                                                  name_b="v1", **kw)
    for key, value in here["difference_decomposition"].items():
        mirrored = back["difference_decomposition"][key]
        if value is None or mirrored is None:
            assert value is mirrored is None
        else:
            assert mirrored == pytest.approx(-value, **RELATIVE)


def test_reordering_a_family_reorders_its_rejections_and_nothing_else():
    family = [0.001, 0.2, 0.04, 0.6, 0.011, 0.9]
    order = [4, 0, 5, 2, 1, 3]
    adjusted, rejected = benjamini_hochberg(family, 0.1)
    moved_adjusted, moved_rejected = benjamini_hochberg(
        [family[i] for i in order], 0.1)
    assert list(moved_rejected) == [rejected[i] for i in order]
    assert list(moved_adjusted) == [adjusted[i] for i in order]


def test_counting_every_record_twice_leaves_the_identification_bounds_alone():
    """The bounds are shares of the frame; replicating the frame is invisible."""
    values = [1, 0, None, 1, 0, None, 1]
    here = missingness.finite_collection_bounds(values, total_count=9)
    twice = missingness.finite_collection_bounds(values + values, total_count=18)
    assert Fraction(twice["lower_bound"]) == Fraction(here["lower_bound"])
    assert Fraction(twice["upper_bound"]) == Fraction(here["upper_bound"])
    assert twice["observed_prevalence"] == pytest.approx(here["observed_prevalence"],
                                                         **RELATIVE)


def test_exchanging_the_two_collections_turns_the_difference_bounds_around():
    a = {"lower_bound": 0.2, "upper_bound": 0.5}
    b = {"lower_bound": 0.1, "upper_bound": 0.9}
    here = missingness.difference_bounds(a, b)
    back = missingness.difference_bounds(b, a)
    assert back["lower_bound"] == pytest.approx(-here["upper_bound"], **RELATIVE)
    assert back["upper_bound"] == pytest.approx(-here["lower_bound"], **RELATIVE)


# --------------------------------------------------------------------------- #
# Relations of a whole run: the record is a function of the data, not of the
# order the data happen to be written in
# --------------------------------------------------------------------------- #

_LINEAGES = [f"L{i // 4 + 1:02d}" for i in range(40)]
_PATTERN = {"ag_alpha": "RSRSS", "ag_beta": "RRSSS", "ag_gamma": "SRRSR"}
_CONFIG = """\
dataset:
  name: metamorphic_fixture
  data_dir: data
  metadata: metadata.csv
  strain_id_column: isolate_id
  lineage_column: lineage
  phenotype: calls.csv
  phenotype_id_column: isolate_id
  phenotype_antibiotic_column: agent
  phenotype_call_column: call
  phenotype_kind: clinical_sir
attribution: {folds: 2, repeats: 2, n_boot: 20, n_perm: 19}
surveillance: {n_boot: 20}
evidence: {folds: 2, repeats: 2}
"""


def _write_case(root: Path, rows, ids=None, config=None):
    """Write one fixture; `rows` is the call table in the order it is given."""
    data = root / "data"
    data.mkdir(parents=True, exist_ok=True)
    names = ids or [f"iso_{i:02d}" for i in range(40)]
    with (data / "metadata.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["isolate_id", "lineage", "period"])
        writer.writerows([[name, lineage, "earlier" if i % 2 else "later"]
                          for i, (name, lineage) in enumerate(zip(names, _LINEAGES))])
    with (data / "calls.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["isolate_id", "agent", "call"])
        writer.writerows(rows)
    (root / "config.yaml").write_text(config or _CONFIG, encoding="utf-8")
    return root / "config.yaml"


def _rows(agents=("ag_alpha", "ag_beta", "ag_gamma"), ids=None):
    names = ids or [f"iso_{i:02d}" for i in range(40)]
    return [[names[i], agent, _PATTERN[agent][i % 5]]
            for agent in agents for i in range(40)]


def _run(config: Path, out: Path) -> dict:
    from amr_clonalshare import cli
    assert cli.main(["--config", str(config), "--results-dir", str(out),
                     "--quiet"]) == 0
    return json.loads((out / "clonal_share_result.json").read_text(encoding="utf-8"))


_CONTRAST_CONFIG = _CONFIG.replace(
    "  phenotype_kind: clinical_sir\n",
    "  phenotype_kind: clinical_sir\n  contrast_column: period\n"
    "  contrast_levels: [earlier, later]\n")


def _numbers(record: dict) -> dict:
    """The reported quantities, keyed by agent, with no file path in sight."""
    diagnostics = record["metadata_diagnostics"]
    out = {}
    for agent, block in diagnostics["clonal_share"].items():
        out[f"share/{agent}"] = {k: v for k, v in block.items() if k != "reason"}
    for agent, block in diagnostics["realised_share"].items():
        out[f"realised/{agent}"] = {k: v for k, v in block.items() if k != "reason"}
    for agent, block in diagnostics["lineage_evidence"]["per_feature"].items():
        out[f"evidence/{agent}"] = dict(block)
    decomposition = diagnostics.get("prevalence_decomposition")
    if decomposition:
        for agent, block in decomposition.get("per_feature", {}).items():
            out[f"decomposition/{agent}"] = {k: v for k, v in block.items()
                                             if k != "reason"}
    return out


def test_the_row_order_of_the_call_table_does_not_change_one_number(tmp_path):
    """A table sorted differently is the same table."""
    straight = _rows()
    shuffled = list(straight)
    np.random.default_rng(17).shuffle(shuffled)
    assert shuffled != straight
    first = _numbers(_run(_write_case(tmp_path / "a", straight), tmp_path / "out_a"))
    second = _numbers(_run(_write_case(tmp_path / "b", shuffled), tmp_path / "out_b"))
    assert first == second


def test_the_order_of_the_antimicrobials_does_not_change_one_number(tmp_path):
    """Reading the drugs in another order is not another analysis."""
    first = _numbers(_run(_write_case(tmp_path / "a", _rows()), tmp_path / "out_a"))
    reversed_rows = _rows(agents=("ag_gamma", "ag_beta", "ag_alpha"))
    second = _numbers(_run(_write_case(tmp_path / "b", reversed_rows), tmp_path / "out_b"))
    assert first == second


def test_adding_an_antimicrobial_does_not_move_the_others(tmp_path):
    """A per-agent estimate is a property of that agent's calls."""
    two = ("ag_alpha", "ag_beta")
    first = _numbers(_run(_write_case(tmp_path / "a", _rows(agents=two)),
                          tmp_path / "out_a"))
    three = _rows(agents=("ag_alpha", "ag_beta", "ag_gamma"))
    second = _numbers(_run(_write_case(tmp_path / "b", three), tmp_path / "out_b"))
    for key, value in first.items():
        assert key in second, key
        assert second[key] == value, key


def test_renaming_the_isolates_in_order_does_not_change_one_number(tmp_path):
    """The identifier labels a record; it is not data about it."""
    plain = [f"iso_{i:02d}" for i in range(40)]
    renamed = [f"strain-{i:04d}" for i in range(40)]
    first = _numbers(_run(_write_case(tmp_path / "a", _rows(ids=plain), ids=plain),
                          tmp_path / "out_a"))
    second = _numbers(_run(_write_case(tmp_path / "b", _rows(ids=renamed), ids=renamed),
                           tmp_path / "out_b"))
    assert first == second


def test_adding_an_antimicrobial_does_not_move_the_decomposition_of_the_others(tmp_path):
    """The bootstrap of one antimicrobial is that antimicrobial's own."""
    two = ("ag_alpha", "ag_beta")
    first = _numbers(_run(_write_case(tmp_path / "a", _rows(agents=two),
                                      config=_CONTRAST_CONFIG), tmp_path / "out_a"))
    three = _rows(agents=("ag_alpha", "ag_beta", "ag_gamma"))
    second = _numbers(_run(_write_case(tmp_path / "b", three,
                                       config=_CONTRAST_CONFIG), tmp_path / "out_b"))
    moved = [key for key, value in first.items()
             if key.startswith("decomposition/") and second.get(key) != value]
    assert not moved, moved
    assert any(key.startswith("decomposition/") for key in first), (
        "the fixture produced no decomposition to compare")


def test_renaming_the_antimicrobials_does_not_change_one_number(tmp_path):
    """The name labels a drug; the analysis reads its calls."""
    plain = _rows()
    renamed = [[row[0], {"ag_alpha": "drug_x", "ag_beta": "drug_y",
                         "ag_gamma": "drug_z"}[row[1]], row[2]] for row in plain]
    first = _numbers(_run(_write_case(tmp_path / "a", plain), tmp_path / "out_a"))
    second = _numbers(_run(_write_case(tmp_path / "b", renamed), tmp_path / "out_b"))
    back = {"drug_x": "ag_alpha", "drug_y": "ag_beta", "drug_z": "ag_gamma"}
    mapped = {f"{kind}/{back[agent]}": value
              for kind, agent, value in ((k.split("/")[0], k.split("/")[1], v)
                                         for k, v in second.items())}
    assert mapped == first
