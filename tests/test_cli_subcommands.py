"""The subcommands spell out the options of the main command; the old form stays valid."""
import json
from pathlib import Path

import pytest
import yaml

from amr_clonalshare import __version__
from amr_clonalshare.cli import EXAMPLES, completion_script, doctor_report, examples_folder, main

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_DIR = ROOT / "examples"
needs_examples = pytest.mark.skipif(not (EXAMPLES_DIR / "workflows").is_dir(), reason="examples folder not shipped")


def test_version_doctor_and_completion(capsys):
    assert main(["version"]) == 0
    assert capsys.readouterr().out.strip() == f"amr-clonalshare {__version__}"
    assert main(["doctor"]) == 0
    out = capsys.readouterr().out
    assert "numpy" in out and "Processor threads" in out and "Excel workbooks" in out
    for shell in ("bash", "powershell"):
        assert main(["completion", shell]) == 0
        script = capsys.readouterr().out
        assert "--results-dir" in script and "example" in script and "ssuis" in script
    assert main(["completion", "fish"]) == 2
    assert "complete -o default" in completion_script("bash")
    assert "Register-ArgumentCompleter" in completion_script("powershell")
    assert doctor_report().startswith("amr-clonalshare")


@needs_examples
def test_run_check_init_and_example_are_the_options_of_the_main_command(tmp_path, capsys):
    config = EXAMPLES_DIR / "workflows" / "mic.yaml"
    assert main(["check", "--config", str(config), "--results-dir", str(tmp_path / "check")]) == 0
    assert capsys.readouterr().out.startswith("# Input check")
    assert (tmp_path / "check" / "input_qc.md").is_file()
    assert main(["init", str(EXAMPLES_DIR / "workflows" / "data" / "calls.csv")]) == 0
    assert "phenotype_call_column" in capsys.readouterr().out
    assert main(["run", "--config", str(config), "--results-dir", str(tmp_path / "run"), "--seed", "3",
                 "--permutations", "7", "--bootstrap", "5", "--alpha", "0.1", "--quiet"]) == 0
    record = json.loads((tmp_path / "run" / "clonal_share_result.json").read_text(encoding="utf-8"))
    assert record["config"]["attribution"]["n_perm"] == 7
    assert record["config"]["attribution"]["n_boot"] == 5
    assert record["config"]["evidence"]["alpha"] == 0.1
    assert (tmp_path / "run" / "run.log").read_text(encoding="utf-8").count("\n") >= 3
    assert "run.log" in json.loads((tmp_path / "run" / "run_manifest.json").read_text(encoding="utf-8"))["files"]
    assert main(["example"]) == 0
    assert "fictional" in capsys.readouterr().out
    assert examples_folder() == EXAMPLES_DIR
    assert main(["example", "fictional", "--results-dir", str(tmp_path / "example"), "--quiet"]) == 0
    record = json.loads((tmp_path / "example" / "clonal_share_result.json").read_text(encoding="utf-8"))
    assert record["seed"] == EXAMPLES["fictional"][1]


@needs_examples
def test_dry_run_prints_the_configuration_with_every_default_and_writes_nothing(tmp_path, capsys):
    config = EXAMPLES_DIR / "workflows" / "mic.yaml"
    assert main(["run", "--config", str(config), "--dry-run", "--results-dir", str(tmp_path / "none"), "--alpha", "0.2"]) == 0
    printed = yaml.safe_load(capsys.readouterr().out)
    assert printed["dataset"]["name"] == "fictional_workflow"
    assert printed["evidence"]["alpha"] == 0.2
    assert printed["surveillance"]["min_shared_support"] == 0.9
    assert not (tmp_path / "none").exists()


def test_overrides_are_checked_like_the_file_values(tmp_path, capsys):
    config = EXAMPLES_DIR / "workflows" / "mic.yaml"
    if not config.is_file():
        pytest.skip("examples folder not shipped")
    assert main(["run", "--config", str(config), "--permutations", "0"]) == 2
    assert "--permutations" in capsys.readouterr().err
    assert main(["run", "--config", str(config), "--alpha", "1.5"]) == 2


@needs_examples
def test_a_failed_run_leaves_its_log_with_the_error(tmp_path, capsys):
    raw = yaml.safe_load((EXAMPLES_DIR / "workflows" / "mic.yaml").read_text(encoding="utf-8"))
    raw["dataset"]["data_dir"] = str(EXAMPLES_DIR / "workflows" / "data")
    raw["dataset"]["lineage_column"] = "no_such_column"
    (tmp_path / "bad.yaml").write_text(yaml.safe_dump(raw), encoding="utf-8")
    assert main(["run", "--config", str(tmp_path / "bad.yaml"), "--results-dir", str(tmp_path / "out")]) == 2
    log = (tmp_path / "out" / "run.log").read_text(encoding="utf-8")
    assert "Reading input tables" in log and "error: input error" in log and "no_such_column" in log
    # a later run replaces the log rather than refusing the folder
    assert main(["run", "--config", str(EXAMPLES_DIR / "workflows" / "mic.yaml"), "--results-dir", str(tmp_path / "out"),
                 "--overwrite", "--quiet"]) == 0
    assert (tmp_path / "out" / "clonal_share_result.json").is_file()
    # a failed run over that finished one leaves it as it was, its log beside it
    finished = (tmp_path / "out" / "run.log").read_text(encoding="utf-8")
    assert main(["run", "--config", str(tmp_path / "bad.yaml"), "--results-dir", str(tmp_path / "out"),
                 "--overwrite"]) == 2
    assert (tmp_path / "out" / "run.log").read_text(encoding="utf-8") == finished
    assert "error: input error" in (tmp_path / "out" / "run_failed.log").read_text(encoding="utf-8")
