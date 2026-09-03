"""The configuration is a contract, so a key it does not read is refused.

A misspelt key is the case that matters. `lineage_colum` silently disables the
lineage diagnostic, and the run then reports no population-structure gate for
the same reason it would report none on a cohort with no lineage labels at all.
The two outcomes are indistinguishable in the output, which is why the key has
to be refused rather than ignored.
"""
from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

from amr_clonalshare.config import ConfigError, from_dict

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


@pytest.fixture(scope="module")
def raw():
    path = EXAMPLES / "ssuis" / "config.yaml"
    if not path.is_file():
        pytest.skip("S. suis example config not present")
    return yaml.safe_load(path.read_text(encoding="utf-8"))


@pytest.mark.parametrize("path", sorted(EXAMPLES.rglob("*.yaml")))
def test_every_shipped_config_still_loads(path):
    from amr_clonalshare.config import load_config
    load_config(path, check_files_exist=False)


def test_an_unknown_top_level_section_is_refused(raw):
    bad = copy.deepcopy(raw)
    bad["nonsense_section"] = {"x": 1}
    with pytest.raises(ConfigError, match="unknown key 'nonsense_section'"):
        from_dict(bad)


def test_an_unknown_key_inside_a_section_is_refused(raw):
    bad = copy.deepcopy(raw)
    bad["dataset"]["nonsense_key"] = 1
    with pytest.raises(ConfigError, match="unknown key 'nonsense_key'"):
        from_dict(bad)


@pytest.mark.parametrize("section,typo,intended", [
    ("dataset", "lineage_colum", "lineage_column"),
    ("dataset", "external_column", "external_columns"),
    ("snf", "alpah", "alpha"),
    ("surveillance", "n_bootstrap", "n_boot"),
])
def test_a_misspelt_key_is_refused_and_the_intended_one_is_named(
        raw, section, typo, intended):
    bad = copy.deepcopy(raw)
    bad.setdefault(section, {})[typo] = 1
    with pytest.raises(ConfigError) as caught:
        from_dict(bad)
    assert typo in str(caught.value)
    assert intended in str(caught.value)


def test_a_nested_section_is_checked_too(raw):
    bad = copy.deepcopy(raw)
    bad["trait_cluster"].setdefault("prevalence_gate", {})["hii"] = 1
    with pytest.raises(ConfigError, match="prevalence_gate"):
        from_dict(bad)


def test_the_surveillance_section_validates_its_own_ranges(raw):
    for key, value in (("q_fdr", 1.5), ("min_shared_support", 2.0),
                       ("label_alpha", 0.0), ("n_boot", -1)):
        bad = copy.deepcopy(raw)
        bad["surveillance"] = {key: value}
        with pytest.raises(ConfigError, match="surveillance"):
            from_dict(bad).validate(check_files_exist=False)


def test_surveillance_defaults_match_the_calibrated_gate(raw):
    from amr_clonalshare.clonality import DEFAULT_MIN_SHARED_SUPPORT
    cfg = from_dict(copy.deepcopy(raw))
    assert cfg.surveillance.min_shared_support == DEFAULT_MIN_SHARED_SUPPORT
    assert cfg.surveillance.n_boot >= 2000
    assert cfg.surveillance.n_perm >= 2000
