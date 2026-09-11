"""Public release-contract checks.

The public product version is exactly 1.0.0. The output schema is versioned
separately and is 1.0.
"""
from __future__ import annotations

import importlib.metadata
import re
from pathlib import Path

from amr_clonalshare import __version__, core


ROOT = Path(__file__).resolve().parents[1]


#: The one place the public product version is written down for the tests.
PUBLIC_VERSION = "1.0.0"

#: The Zenodo concept DOI of the deposit. It always resolves to the version
#: record deposited from this tag, so every file that cites the deposit must
#: cite this one string and not a version DOI of its own.
CONCEPT_DOI = "10.5281/zenodo.22306353"


def test_all_current_public_version_sources_agree():
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    zenodo = (ROOT / ".zenodo.json").read_text(encoding="utf-8")
    v = re.escape(PUBLIC_VERSION)

    assert re.search(rf'^version = "{v}"$', pyproject, flags=re.MULTILINE)
    assert __version__ == PUBLIC_VERSION
    assert importlib.metadata.version("amr-clonalshare") == PUBLIC_VERSION
    assert f'version: "{PUBLIC_VERSION}"' in citation
    assert f'"version": "{PUBLIC_VERSION}"' in zenodo


def test_first_public_output_schema_is_1_0(share_cfg):
    _, cfg = share_cfg
    result = core.run(cfg, seed=42)
    assert result["schema_version"] == "1.0"


def test_every_file_that_cites_the_deposit_cites_the_same_concept_doi():
    for name in ("README.md", "CITATION.cff", "REPRODUCIBILITY.md",
                 "docs/CODE_METADATA.md"):
        assert CONCEPT_DOI in (ROOT / name).read_text(encoding="utf-8"), name


def test_softwarex_licence_copy_is_byte_identical():
    assert (ROOT / "Licence.txt").read_bytes() == (ROOT / "LICENSE").read_bytes()
