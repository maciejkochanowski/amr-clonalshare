"""The artefacts a run writes, compared whole against a committed copy.

A test that reads one field of a result checks that field. These tests read
what the user is actually given -- the record, the two input diagnostics, the
result table and both reports -- and compare them in full. A renamed key, a
reworded heading, a changed rounding or a moved number therefore fails here,
and the diff of the golden file is the description of what changed.

Regenerate after an intended change with `python scripts/update_golden.py`,
read the diff, and commit it together with the change that caused it.
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
GOLDEN = Path(__file__).resolve().parent / "golden"

# The two parts of an artefact that a second run, or another machine, does not
# reproduce: the hour it was issued, and the digest of a record that names the
# dependency versions it was produced with. The record's own content is
# compared separately and in full, and the manifest digests are recomputed
# from the files they describe, so nothing is lost by normalising them here.
_ISSUED = re.compile(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}Z")
_DIGEST = re.compile(r"sha256:[0-9a-f]{6,}")

PIPELINE_CASES = ("calls", "contrasts", "mic", "rough")
PIPELINE_FILES = ("clonal_share_result.json", "input_qc.json", "input_qc.md",
                  "report.html", "report.md", "results.csv")
COMPARE_FILES = ("comparison.json", "report.md", "arms.csv")
DRAFT_CASES = {
    "one_table": ["tests/golden/draft_data/one_table.csv"],
    "wide": ["tests/golden/draft_data/wide.csv"],
    "unnamed": ["tests/golden/draft_data/unnamed.csv"],
    "two_tables": ["examples/workflows/data/metadata.csv",
                   "examples/workflows/data/calls.csv"],
}
CASES = PIPELINE_CASES + ("compare", "draft")


def normalise(text: str) -> str:
    return _DIGEST.sub("sha256:<digest>", _ISSUED.sub("<issued>", text))


def run_case(name: str, out_dir: Path) -> None:
    """Write the artefacts of one case into `out_dir`, from the repository root."""
    import os

    out_dir.mkdir(parents=True, exist_ok=True)
    here = Path.cwd()
    os.chdir(ROOT)
    try:
        if name in PIPELINE_CASES:
            from amr_clonalshare import cli
            status = cli.main(["--config", str(GOLDEN / f"{name}.yaml"),
                               "--results-dir", str(out_dir), "--quiet"])
            assert status == 0, f"{name} exited {status}"
        elif name == "compare":
            from amr_clonalshare import comparison
            status = comparison.main([
                "--input", str(ROOT / "examples/matched_lineages/input.csv"),
                "--id-column", "isolate_id", "--outcome", "positive",
                "--lineage-a", "lineage_a", "--lineage-b", "lineage_b",
                "--output", str(out_dir)])
            assert status == 0, f"compare exited {status}"
        elif name == "draft":
            from amr_clonalshare.draft import draft_config
            for case, paths in DRAFT_CASES.items():
                # The draft names its data folder by its absolute path, which
                # differs between checkouts; the golden copy names it relative
                # to the repository.
                text = draft_config([Path(p) for p in paths]).replace(
                    json.dumps(str(ROOT)).strip('"'), "<repository>")
                (out_dir / f"{case}.yaml").write_text(text, encoding="utf-8")
        else:  # pragma: no cover - the case table is closed
            raise AssertionError(name)
    finally:
        os.chdir(here)


def case_files(name: str) -> tuple[str, ...]:
    if name in PIPELINE_CASES:
        return PIPELINE_FILES
    if name == "compare":
        return COMPARE_FILES
    return tuple(sorted(f"{case}.yaml" for case in DRAFT_CASES))


def _numbers_agree(got, want) -> bool:
    if isinstance(got, bool) or isinstance(want, bool):
        return got is want
    if not isinstance(got, (int, float)) or not isinstance(want, (int, float)):
        return False
    if got != got and want != want:  # both nan
        return True
    return abs(got - want) <= 1e-9 * max(1.0, abs(want))


def compare_json(got, want, where: str = "") -> list[str]:
    """Keys and strings exactly; numbers to the precision a rerun reproduces."""
    if isinstance(want, dict):
        if not isinstance(got, dict):
            return [f"{where}: expected an object, got {type(got).__name__}"]
        problems = []
        for key in sorted(set(want) | set(got)):
            if key not in got:
                problems.append(f"{where}/{key}: missing from the run")
            elif key not in want:
                problems.append(f"{where}/{key}: not in the golden record")
            else:
                problems += compare_json(got[key], want[key], f"{where}/{key}")
        return problems
    if isinstance(want, list):
        if not isinstance(got, list) or len(got) != len(want):
            return [f"{where}: expected {len(want) if isinstance(want, list) else want} "
                    f"entries, got {len(got) if isinstance(got, list) else got}"]
        problems = []
        for i, (g, w) in enumerate(zip(got, want)):
            problems += compare_json(g, w, f"{where}[{i}]")
        return problems
    if isinstance(want, (int, float)) and not isinstance(want, bool):
        if not _numbers_agree(got, want):
            return [f"{where}: {got!r} is not {want!r}"]
        return []
    if got != want:
        return [f"{where}: {got!r} is not {want!r}"]
    return []


def compare_csv(got_text: str, want_text: str) -> list[str]:
    got = list(csv.reader(got_text.splitlines()))
    want = list(csv.reader(want_text.splitlines()))
    if len(got) != len(want):
        return [f"{len(got)} rows, golden has {len(want)}"]
    problems = []
    for r, (grow, wrow) in enumerate(zip(got, want)):
        if len(grow) != len(wrow):
            problems.append(f"row {r}: {len(grow)} fields, golden has {len(wrow)}")
            continue
        for c, (g, w) in enumerate(zip(grow, wrow)):
            if g == w:
                continue
            try:
                if _numbers_agree(float(g), float(w)):
                    continue
            except ValueError:
                pass
            problems.append(f"row {r} field {c} ({want[0][c] if r else 'header'}): "
                            f"{g!r} is not {w!r}")
    return problems


# Provenance, not result: the record names the versions it was produced with,
# and two machines with different point releases write different strings while
# writing the same numbers. Checked separately, below.
VOLATILE = ("versions",)


def _without_provenance(record):
    if not isinstance(record, dict):
        return record
    return {k: v for k, v in record.items() if k not in VOLATILE}


def compare_artefact(name: str, got_text: str, want_text: str) -> list[str]:
    if name.endswith(".json"):
        return compare_json(_without_provenance(json.loads(got_text)),
                            _without_provenance(json.loads(want_text)))
    if name.endswith(".csv"):
        return compare_csv(got_text, want_text)
    got, want = normalise(got_text), normalise(want_text)
    if got == want:
        return []
    got_lines, want_lines = got.splitlines(), want.splitlines()
    for i, (g, w) in enumerate(zip(got_lines, want_lines), start=1):
        if g != w:
            return [f"line {i}: {g!r}\n   golden: {w!r}"]
    return [f"{len(got_lines)} lines, golden has {len(want_lines)}"]


@pytest.mark.parametrize("case", CASES)
def test_written_artefacts_match_the_committed_copy(case, tmp_path):
    run_case(case, tmp_path)
    for name in case_files(case):
        golden = GOLDEN / case / name
        assert golden.exists(), (
            f"{golden} is missing; run python scripts/update_golden.py")
        problems = compare_artefact(name, (tmp_path / name).read_text(encoding="utf-8"),
                                    golden.read_text(encoding="utf-8"))
        assert not problems, (
            f"{case}/{name} differs from the committed artefact:\n  "
            + "\n  ".join(problems[:12]))


@pytest.mark.parametrize("case", PIPELINE_CASES)
def test_the_record_names_the_versions_it_was_produced_with(case, tmp_path):
    """The one part of the record that two machines may write differently."""
    run_case(case, tmp_path)
    versions = json.loads((tmp_path / "clonal_share_result.json")
                          .read_text(encoding="utf-8"))["versions"]
    assert {"python", "numpy", "pandas", "scipy"} <= set(versions)
    assert all(isinstance(v, str) and v for v in versions.values())


@pytest.mark.parametrize("case", PIPELINE_CASES)
def test_the_manifest_lists_every_artefact_with_its_digest(case, tmp_path):
    import hashlib

    run_case(case, tmp_path)
    manifest = json.loads((tmp_path / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "complete"
    assert sorted(manifest["files"]) == sorted(PIPELINE_FILES)
    for name, digest in manifest["files"].items():
        written = hashlib.sha256((tmp_path / name).read_bytes()).hexdigest()
        assert written == digest, f"{name}: the manifest digest is not the file"


def test_a_reported_number_is_never_written_at_full_precision_in_prose():
    """Text artefacts are compared as text, which only works while they round."""
    long_decimal = re.compile(r"\d\.\d{6,}")
    for case in PIPELINE_CASES + ("compare",):
        for name in case_files(case):
            if not name.endswith((".md", ".html")):
                continue
            text = (GOLDEN / case / name).read_text(encoding="utf-8")
            found = long_decimal.findall(text)
            assert not found, f"{case}/{name} prints {found[:3]} unrounded"
