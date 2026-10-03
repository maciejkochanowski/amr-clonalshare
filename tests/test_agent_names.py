"""Agent codes and abbreviations are read as the agents they stand for."""
import pandas as pd
import pytest

from amr_clonalshare.agents import canonical_agent, unify_agents
from amr_clonalshare.config import from_dict
from amr_clonalshare.io import load_dataset


@pytest.mark.parametrize("written,read", [
    ("CIP", "ciprofloxacin"), ("cip", "ciprofloxacin"), ("CIP_NM", "ciprofloxacin"), ("CIP_ND5", "ciprofloxacin"),
    ("GEN (mg/L)", "gentamicin"), ("SXT", "trimethoprim-sulfamethoxazole"),
    ("trimethoprim/sulfamethoxazole", "trimethoprim-sulfamethoxazole"),
    ("AMC", "amoxicillin-clavulanic acid"), ("Pen G", "penicillin"),
    ("ciprofloxacin", "ciprofloxacin"), ("amoxicillin-clavulanic acid", "amoxicillin-clavulanic acid"),
    ("demo_agent", "demo_agent"), ("ST", "ST"), ("  tetracycline ", "tetracycline"),
])
def test_codes_are_read_as_names_and_names_are_kept(written, read):
    assert canonical_agent(written) == read


def test_the_names_the_shipped_collections_use_are_their_own():
    shipped = ["amikacin", "amoxicillin-clavulanic acid", "ampicillin", "azithromycin", "cefoxitin",
               "ceftazidime", "ceftiofur", "ceftriaxone", "chloramphenicol", "ciprofloxacin", "colistin",
               "gentamicin", "imipenem", "kanamycin", "meropenem", "nalidixic acid", "spectinomycin",
               "streptomycin", "sulfamethoxazole", "sulfisoxazole", "tetracycline",
               "trimethoprim-sulfamethoxazole", "amoxicillin", "cefquinome", "doxycycline", "enrofloxacin",
               "erythromycin", "florfenicol", "lincomycin", "marbofloxacin", "penicillin", "tiamulin",
               "tilmicosin", "trimethoprim", "tylosin", "cefotaxime", "demo_agent", "agent_a"]
    mapping, changed = unify_agents(shipped)
    assert changed == {}
    assert list(mapping) == shipped


def _tables(tmp_path, agents):
    tmp_path.mkdir(exist_ok=True)
    (tmp_path / "metadata.csv").write_text("id,st\nA1,ST1\nA2,ST1\nA3,ST2\nA4,ST2\nA5,ST3\nA6,ST3\n", encoding="utf-8")
    calls = ["id,agent,call"]
    for i, isolate in enumerate(("A1", "A2", "A3", "A4", "A5", "A6")):
        for j, agent in enumerate(agents):
            calls.append(f"{isolate},{agent},{'R' if (i + j) % 2 else 'S'}")
    (tmp_path / "calls.csv").write_text("\n".join(calls) + "\n", encoding="utf-8")
    return from_dict({"dataset": {"name": "t", "data_dir": str(tmp_path), "metadata": "metadata.csv",
                                  "strain_id_column": "id", "lineage_column": "st", "phenotype": "calls.csv",
                                  "phenotype_id_column": "id", "phenotype_antibiotic_column": "agent",
                                  "phenotype_call_column": "call", "phenotype_kind": "clinical_sir"}})


def test_a_call_table_with_codes_gives_the_panel_of_the_same_table_with_names(tmp_path):
    coded = load_dataset(_tables(tmp_path / "coded", ["CIP", "GEN", "SXT"]))
    named = load_dataset(_tables(tmp_path / "named", ["ciprofloxacin", "gentamicin", "trimethoprim-sulfamethoxazole"]))
    pd.testing.assert_frame_equal(coded.panel, named.panel)
    assert coded.input_qc["phenotype_reading"]["agent_names_unified"] == {
        "CIP": "ciprofloxacin", "GEN": "gentamicin", "SXT": "trimethoprim-sulfamethoxazole"}
    assert "agent_names_unified" not in named.input_qc["phenotype_reading"]
