#!/usr/bin/env python3
"""Place the outputs of the re-issue runs in the results folders and check them.

    python benchmarks/campaign/collect_reissue.py CAMPAIGN_ROOT [--apply]

Without --apply nothing is written: the script reads CAMPAIGN_ROOT/outputs,
checks that every receipt names the package the checked-out source is
(the digest of every Python file of src/amr_clonalshare must equal the
receipt's), that every step exited 0, and reports what would change against
the results folders in the repository: the example records (the MIC ordering
must be unchanged field by field; every changed field of the call analyses is
listed), the decomposition gate, the null-uniformity summary, the repeated-looks
rules and the test logs. With --apply it copies the outputs over the results
folders (benchmarks/results_empirical, examples/*/expected*,
benchmarks/results_decomposition_calibration,
benchmarks/results_decomposition_vs_regression,
benchmarks/results_null_uniformity, benchmarks/results_repeated_looks,
benchmarks/results_confirmatory/campaign/pytest_<commit>.log) and, when a
mutation run is present, leaves its survivors file beside the record for the
classification of any new survivor (scripts/check_mutation_survivors.py).
"""
from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "benchmarks" / "campaign"))
from run_logged import source_digest

EMPIRICAL = {
    # output folder under outputs/empirical -> results folder, example folder
    "ssuis_run": ("benchmarks/results_empirical/ssuis_run", "examples/ssuis/expected"),
    "ssuis_contrast": ("benchmarks/results_empirical/ssuis_contrast", None),
    "salmonella_serovar": ("benchmarks/results_empirical/salmonella_serovar", "examples/salmonella_poultry/expected_serovar"),
    "salmonella_cluster": ("benchmarks/results_empirical/salmonella_cluster", "examples/salmonella_poultry/expected_cluster"),
    "ecoli": ("benchmarks/results_empirical/ecoli", "examples/ecoli_swine/expected"),
    "ecoli_two_categories": ("benchmarks/results_empirical/ecoli_two_categories", "examples/ecoli_swine/expected_two_categories"),
    "ecoli_without_st410": ("benchmarks/results_empirical/ecoli_without_st410", "examples/ecoli_swine/expected_without_st410"),
    "ecoli_periods": ("benchmarks/results_empirical/ecoli_periods", "examples/ecoli_swine/expected_periods"),
    "ecoli_comparison": ("benchmarks/results_empirical/ecoli_comparison", "examples/ecoli_swine/expected_comparison"),
    "ecoli_comparison_clermont": ("benchmarks/results_empirical/ecoli_comparison_clermont",
                                  "examples/ecoli_swine/expected_comparison_clermont"),
    "ssuis": ("benchmarks/results_empirical/ssuis", None),
    "ssuis_panel": ("benchmarks/results_empirical/ssuis_panel", None),
    "ssuis_sensitivity": ("benchmarks/results_empirical/ssuis_sensitivity", None),
    "salmonella": ("benchmarks/results_empirical/salmonella", None),
    "seed_stability": ("benchmarks/results_empirical/seed_stability", None),
    "profile": ("benchmarks/results_empirical/profile", None),
}
# the two empirical scripts write their receipt under the name of the script
RECEIPT_NAMES = {"ssuis": "ssuis_analysis", "salmonella": "salmonella_analysis"}
# what a second run does not reproduce (tests/test_golden_artefacts.py masks the same)
_ISSUED = re.compile(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}Z")
_DIGEST = re.compile(r"sha256:[0-9a-f]{6,}")
AUXILIARY = {
    "decomposition_calibration": "benchmarks/results_decomposition_calibration",
    "decomposition_vs_regression": "benchmarks/results_decomposition_vs_regression",
    "null_uniformity": "benchmarks/results_null_uniformity",
    "repeated_looks": "benchmarks/results_repeated_looks",
}

problems: list[str] = []
notes: list[str] = []


def only_volatile(patch: Path) -> bool:
    """True when every changed line of a diff differs only in the issue time or a digest."""
    removed, added = [], []
    for line in patch.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("---", "+++")):
            continue
        if line[:1] in "-+":
            text = _DIGEST.sub("sha256:<digest>", _ISSUED.sub("<issued>", line[1:]))
            (removed if line[0] == "-" else added).append(text)
    return sorted(removed) == sorted(added)


def flatten(d, prefix=""):
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(flatten(v, f"{prefix}/{k}"))
    elif isinstance(d, list):
        for i, v in enumerate(d):
            out.update(flatten(v, f"{prefix}[{i}]"))
    else:
        out[prefix] = d
    return out


def same(a, b) -> bool:
    if isinstance(a, float) and isinstance(b, float):
        return (math.isnan(a) and math.isnan(b)) or a == b
    return a == b


def check_receipt(path: Path, digests: dict) -> None:
    if not path.exists():
        problems.append(f"missing receipt {path}")
        return
    receipt = json.loads(path.read_text(encoding="utf-8"))
    if receipt.get("exit_status") != 0:
        problems.append(f"{path}: exit status {receipt.get('exit_status')}")
    files = receipt.get("source", {}).get("files", {})
    differing = sorted(k for k in set(files) | set(digests) if files.get(k) != digests.get(k))
    if differing:
        problems.append(f"{path}: package digests differ from src/amr_clonalshare in {', '.join(differing)}")
    stack = receipt.get("stack", {})
    notes.append(f"{path.name}: commit {str(receipt.get('commit'))[:7]}, numpy {stack.get('numpy')}, "
                 f"pandas {stack.get('pandas')}, scipy {stack.get('scipy')}, {receipt.get('wall_seconds')} s")


def number(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def compare_records(old_path: Path, new_path: Path, label: str) -> None:
    if not new_path.exists():
        problems.append(f"{label}: record missing ({new_path})")
        return
    if not old_path.exists():
        notes.append(f"{label}: a new results folder, nothing to compare against")
        return
    old = flatten(json.loads(old_path.read_text(encoding="utf-8")))
    new = flatten(json.loads(new_path.read_text(encoding="utf-8")))
    changed = [k for k in old if k in new and not same(old[k], new[k])]
    # a number of the MIC ordering must not move; its wording may
    mic = [k for k in changed if "/censored_share/" in k and not k.endswith("/mic_units") and number(old[k])]
    if mic:
        problems.append(f"{label}: {len(mic)} MIC-ordering numbers changed, e.g. {mic[:5]}")
    else:
        notes.append(f"{label}: MIC ordering unchanged number by number ({len(changed)} fields changed, "
                     f"{sum(1 for k in changed if number(old[k]))} of them numbers; "
                     f"{len([k for k in new if k not in old])} added)")
    calls = sorted({k.split('/')[-1] for k in changed if '/clonal_share/' in k})
    if calls:
        notes.append(f"{label}: call-analysis fields changed: {', '.join(calls)}")


def compare_scalar(label: str, old, new) -> None:
    if same(old, new):
        notes.append(f"{label}: unchanged ({old})")
    else:
        notes.append(f"{label}: {old} -> {new}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("campaign_root", type=Path)
    parser.add_argument("--apply", action="store_true", help="copy the outputs over the results folders")
    args = parser.parse_args(argv)
    out = args.campaign_root / "outputs"
    digests = source_digest()["files"]
    copies: list[tuple[Path, Path]] = []

    # empirical queue
    for name, (results, example) in EMPIRICAL.items():
        src = out / "empirical" / name
        receipt = out / "empirical" / "receipts" / f"{RECEIPT_NAMES.get(name, name)}.json"
        if not src.exists():
            problems.append(f"missing output folder {src}")
            continue
        check_receipt(receipt, digests)
        copies.append((receipt, ROOT / "benchmarks" / "results_empirical" / "receipts" / receipt.name))
        for path in sorted(src.rglob("*")):
            if path.is_file():
                copies.append((path, ROOT / results / path.relative_to(src)))
                if example:
                    copies.append((path, ROOT / example / path.relative_to(src)))
        if (src / "clonal_share_result.json").exists():
            compare_records(ROOT / results / "clonal_share_result.json", src / "clonal_share_result.json", name)

    # auxiliary queue
    for name, results in AUXILIARY.items():
        src = out / "aux" / name
        if not src.exists():
            problems.append(f"missing output folder {src}")
            continue
        if name == "repeated_looks":
            for i in range(24):
                cell = src / f"cell_{i}"
                check_receipt(cell / "RUN_RECEIPT.json", digests)
                copies.append((cell / f"cell_{i}.json", ROOT / results / f"cell_{i}.json"))
                copies.append((cell / "RUN_RECEIPT.json", ROOT / results / f"RUN_RECEIPT_cell_{i}.json"))
                old_cell = ROOT / results / f"cell_{i}.json"
                if old_cell.exists() and (cell / f"cell_{i}.json").exists():
                    old = json.loads(old_cell.read_text())[0]
                    new = json.loads((cell / f"cell_{i}.json").read_text())[0]
                    for rule in ("bh_accumulated", "by_accumulated", "lond_batch", "ebh_accumulated", "ebh_product"):
                        if rule in old and rule in new:
                            compare_scalar(f"repeated_looks cell {i} {rule} fdp_ever",
                                           round(old[rule]["fdp_ever"]["mean"], 4), round(new[rule]["fdp_ever"]["mean"], 4))
            continue
        for receipt_name in ("RUN_RECEIPT.json", "GATE_RECEIPT.json"):
            if (src / receipt_name).exists() or receipt_name == "RUN_RECEIPT.json":
                check_receipt(src / receipt_name, digests)
        for path in sorted(src.rglob("*")):
            if path.is_file():
                copies.append((path, ROOT / results / path.relative_to(src)))
        if name == "decomposition_calibration":
            old = ROOT / results / "decomposition_gate.json"
            new = src / "decomposition_gate.json"
            if old.exists() and new.exists():
                compare_scalar("decomposition gate", json.loads(old.read_text())["gate"], json.loads(new.read_text())["gate"])
        if name == "null_uniformity":
            old = ROOT / results / "null_uniformity.json"
            new = src / "null_uniformity.json"
            if old.exists() and new.exists():
                o, n = json.loads(old.read_text()), json.loads(new.read_text())
                for key in ("null", "power_at_0.05"):
                    compare_scalar(f"null_uniformity {key}", json.dumps(o.get(key), sort_keys=True)[:200],
                                   json.dumps(n.get(key), sort_keys=True)[:200])

    # checks: the test logs, and the golden artefacts regenerated on the campaign stack
    checks = out / "checks"
    goldens_patch = None
    if checks.exists():
        for log in sorted(checks.glob("pytest_*.log")):
            tail = log.read_text(encoding="utf-8", errors="replace").strip().splitlines()[-1:]
            notes.append(f"{log.name}: {tail[0] if tail else '(empty)'}")
            copies.append((log, ROOT / "benchmarks" / "results_confirmatory" / "campaign" / log.name))
        goldens_patch = checks / "goldens.patch"
        if goldens_patch.exists() and goldens_patch.stat().st_size and not only_volatile(goldens_patch):
            stat = (checks / "goldens_stat.txt").read_text(errors="replace").strip().splitlines()
            notes.append("goldens regenerated on the campaign stack: " + (stat[-1] if stat else "(no stat)"))
            values = [line for line in stat if "tests/golden/values" in line]
            if values:
                notes.append(f"goldens: {len(values)} value record(s) differ from the committed ones; read the patch "
                             "before applying it (a value record should not move without a change of the code)")
        else:
            goldens_patch = None
            notes.append("goldens: no difference from the committed artefacts beyond the issue time and digests")
    else:
        problems.append("missing outputs/checks (the test suite and the oracles were not run)")

    # mutation
    mutation = out / "mutation"
    if mutation.exists():
        check_receipt(mutation / "RUN_RECEIPT.json", digests)
        # survivors_final.txt exists when mutants whose test process crashed in
        # the first pass were run again alone; it carries their verdicts too
        survivors = mutation / "survivors_final.txt"
        if not survivors.exists():
            survivors = mutation / "survivors.txt"
        if survivors.exists():
            commit = json.loads((mutation / "RUN_RECEIPT.json").read_text())["commit"][:7]
            copies.append((survivors, ROOT / "benchmarks" / "results_confirmatory" / "campaign" / f"mutation_{commit}.txt"))
            log = mutation / "check_mutation_survivors.log"
            notes.append("mutation: " + (log.read_text(errors="replace").strip().splitlines()[-1] if log.exists() else "no check log"))
            notes.append("mutation: classify every new survivor in tests/mutation_survivors.json and update its 'run' block "
                         "from mutmut_run.log before the release")

    for line in notes:
        print("note  " + line)
    for line in problems:
        print("PROBLEM " + line)
    if args.apply and not problems:
        for src, dst in copies:
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
        print(f"{len(copies)} files copied into the repository")
        if goldens_patch is not None:
            import subprocess
            done = subprocess.run(["git", "-C", str(ROOT), "apply", str(goldens_patch)], capture_output=True, text=True)
            print("goldens patch applied" if done.returncode == 0 else f"goldens patch NOT applied: {done.stderr.strip()}")
    elif args.apply:
        print("nothing copied: fix the problems above first")
    else:
        print(f"{len(copies)} files would be copied; rerun with --apply")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
