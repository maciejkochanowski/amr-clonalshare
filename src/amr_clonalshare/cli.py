"""Command-line entry point.

    amr-clonalshare run --config CONFIG --results-dir DIR [--seed N]
    amr-clonalshare init TABLE [TABLE ...]
    amr-clonalshare check --config CONFIG [--results-dir DIR]
    amr-clonalshare example ecoli|ssuis|salmonella --results-dir DIR
    amr-clonalshare compare ... | gui ... | doctor | version | completion SHELL

The same options without a subcommand (``amr-clonalshare --config ...``,
``--init``, ``--check-input``) remain valid, as do ``amr-clonalshare-compare``
and ``amr-clonalshare-gui``.

Exit codes
----------
0   the run completed, including a run whose estimate was refused; the
    record then names the condition that failed
2   configuration or input error
"""
from __future__ import annotations

import argparse
import os
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, Optional

from . import __version__
from .config import ConfigError, load_config


SUBCOMMANDS = ("run", "init", "check", "example", "compare", "gui", "doctor", "version", "completion")
#: The shipped collections ``example`` runs, with the configuration file of
#: each relative to the examples folder and the seed of the published record.
EXAMPLES = {
    "ecoli": ("ecoli_swine/config.yaml", 20261001),
    "ecoli-two-categories": ("ecoli_swine/config_two_categories.yaml", 20261001),
    "ecoli-without-st410": ("ecoli_swine/config_without_st410.yaml", 20261001),
    "ecoli-periods": ("ecoli_swine/config_periods.yaml", 20261001),
    "ssuis": ("ssuis/config.yaml", 42),
    "ssuis-contrast": ("ssuis/config_contrast.yaml", 42),
    "salmonella": ("salmonella_poultry/config.yaml", 42),
    "salmonella-clusters": ("salmonella_poultry/config_cluster.yaml", 42),
    "fictional": ("workflows/mic.yaml", 42),
    "fictional-calls": ("workflows/calls.yaml", 42),
    "fictional-contrasts": ("workflows/contrasts.yaml", 42),
}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="amr-clonalshare",
        description="How much of a susceptibility result do the lineages of a collection account for.\n"
                    "One run reads a lineage label per isolate with calls, MICs or both, and writes a\n"
                    "record, a results table and a report.\n"
                    "\n"
                    "  amr-clonalshare run --config analysis.yaml --results-dir out/run\n"
                    "  amr-clonalshare init data/metadata.csv data/mic.csv > analysis.yaml\n"
                    "  amr-clonalshare example ecoli --results-dir out/ecoli\n"
                    "\n"
                    "Subcommands: run, init, check, example, compare, gui, doctor, version, completion;\n"
                    "the options below also work without a subcommand.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--version", action="version",
                   version=f"amr-clonalshare {__version__}")
    p.add_argument("--config", help="path to a dataset YAML config")
    p.add_argument("--init", nargs="+", metavar="CSV", default=None,
                   help="read the column names of these CSV tables, print a first "
                        "configuration with every guess marked, and stop; it writes "
                        "no file and runs no analysis")
    p.add_argument("--results-dir", default=None,
                   help="directory for clonal_share_result.json and the two reports")
    p.add_argument("--seed", type=int, default=42,
                   help="master seed; every stochastic stage is spawned from it")
    p.add_argument("--threads", type=int, default=None,
                   help="limit BLAS/OpenMP threads (sets OMP_NUM_THREADS et al.)")
    p.add_argument("--check-input", action="store_true",
                   help="load and check the input data, print the input check "
                        "in plain language, write input_qc.json and "
                        "input_qc.md to --results-dir if given, and stop "
                        "before any estimate")
    p.add_argument("--json", action="store_true",
                   help="print the machine-readable summary of the run instead of the results table")
    p.add_argument("--quiet", action="store_true",
                   help="suppress the stdout summary of a run; the input check is still printed")
    p.add_argument("--overwrite", action="store_true",
                   help="replace an existing results directory after successful completion; "
                        "files in it that the tool did not write are kept")
    p.add_argument("--no-check-files", action="store_true",
                   help="skip the existence check on configured data files")
    p.add_argument("--permutations", type=int, default=None, metavar="N",
                   help="permutations of the lineage labels, replacing attribution.n_perm of the configuration")
    p.add_argument("--bootstrap", "--bootstraps", type=int, default=None, metavar="N",
                   help="bootstrap draws of the interval, replacing attribution.n_boot")
    p.add_argument("--alpha", type=float, default=None,
                   help="level of the selection and the e-value rule, replacing evidence.alpha")
    p.add_argument("--dry-run", action="store_true",
                   help="print the configuration as the run would read it, every default filled in, and stop")
    return p


def _set_threads(n: int) -> None:
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[var] = str(n)


def _surveillance_summary(meta: dict) -> dict:
    """Compact reading of the lineage-aware surveillance block.

    The headline the full record buries is the *offsetting* case: an agent
    whose reported prevalence barely moved while both components moved a long
    way in opposite directions. A report quoting prevalence alone calls that
    agent stable, which is the one reading the decomposition exists to prevent.
    """
    out: Dict[str, Any] = {}
    dec = meta.get("prevalence_decomposition") or {}
    per = dec.get("per_feature") or {}
    rows = {f: r for f, r in per.items() if r.get("status") == "ok"}
    if rows:
        # Discoveries are the step-up's, within each component family, and a
        # within-lineage component the shared-support gate refused is not
        # one whatever its interval says.
        comp = {f for f, r in rows.items()
                if r.get("composition_discovery")
                and r.get("composition_estimable", True)}
        within = {f for f, r in rows.items()
                  if r.get("within_lineage_discovery")
                  and r.get("within_lineage_estimable")}
        refused = sorted(f for f, r in rows.items()
                         if not r.get("within_lineage_estimable"))
        refused_comp = sorted(f for f, r in rows.items()
                              if not r.get("composition_estimable", True))
        offsetting = sorted(
            f for f in comp & within
            if rows[f]["composition"] * rows[f]["within_lineage"] < 0
            and abs(rows[f]["difference"]) < min(abs(rows[f]["composition"]),
                                                 abs(rows[f]["within_lineage"])))
        fam = dec.get("family") or {}
        out["decomposition"] = {
            "contrast": f"{dec.get('contrast_column')}: "
                        f"{' minus '.join(dec.get('levels') or [])}",
            "n_features": len(rows),
            "n_composition_discoveries": len(comp),
            "n_within_lineage_discoveries": len(within),
            "n_within_lineage_refused": len(refused),
            "within_lineage_refused_features": refused,
            "n_composition_refused": len(refused_comp),
            "composition_refused_features": refused_comp,
            "n_offsetting": len(offsetting),
            "offsetting_features": offsetting,
            "method": fam.get("method"),
            "smallest_attainable_q": fam.get("smallest_attainable_q"),
            "warning": fam.get("warning"),
            "note": "offsetting: both components discoveries, opposite in "
                    "sign, and each larger than the difference they produce",
        }
    return out


def _evidence_summary(meta: dict) -> Optional[dict]:
    le = meta.get("lineage_evidence") or {}
    ebh = le.get("e_bh") or {}
    per = le.get("per_feature") or {}
    if not per:
        return None
    return {
        "n_agents": len(per),
        "alpha": ebh.get("alpha"),
        "threshold": ebh.get("threshold"),
        "n_rejected": ebh.get("n_rejected"),
        "rejected_features": list(ebh.get("rejected_features") or []),
        "note": le.get("note"),
    }


def _share_summary(share: dict) -> dict:
    """The headline of a run: one row per antimicrobial, and the
    conditions that refused one."""
    md = (share.get("metadata_diagnostics") or {})
    per = (md.get("clonal_share") or {})
    rows = {}
    for agent, d in per.items():
        rows[agent] = {
            "share": d.get("kappa_adj"),
            "ci95": [d.get("observed_low"), d.get("observed_high")],
            "latent_bounds": [d.get("latent_order_lower"), d.get("latent_order_upper_bound")],
            "latent_lower_limit": d.get("latent_order_lower_limit"),
            "permuted": d.get("null_mean"),
            "p_value": d.get("p_value"),
            "support": d.get("support"),
            "estimable": d.get("estimable"),
        }
    ev = (md.get("lineage_evidence") or {}).get("e_bh") or {}
    seq = (md.get("lineage_evidence") or {}).get("sequential") or {}
    return {
        "n_isolates": share.get("n_isolates"),
        "lineage_column": md.get("lineage_column"),
        "n_traits": share.get("n_traits"),
        "n_estimable": sum(1 for r in rows.values() if r["estimable"]),
        "per_agent": rows,
        "selected": ((md.get("lineage_selection") or {}).get("rejected_features")),
        "e_bh_rejected": ev.get("rejected_features"),
        "sequential_e_bh_rejected": ((seq.get("e_bh") or {})
                                     .get("rejected_features")),
    }


def _summary(record: dict) -> dict:
    """The digest of a run: the headline rows and the readings the report needs."""
    md = record.get("metadata_diagnostics") or {}
    out = _share_summary(record)
    out["evidence"] = _evidence_summary(md)
    out["surveillance"] = _surveillance_summary(md) or None
    from .outputs import analysis_summary
    out['analyses'] = analysis_summary(record)
    return out


def _warning_line(message, category, filename, lineno, line=None) -> str:
    """A warning in the voice of the command line, not of the interpreter.

    The default format prints the installed file and the source line that
    raised the warning, which says nothing to the reader of a run.
    """
    del category, filename, lineno, line
    return f"[amr-clonalshare] warning: {message}\n"


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in SUBCOMMANDS:
        return _subcommand(argv[0], argv[1:])
    args = build_parser().parse_args(argv)
    warnings.formatwarning = _warning_line
    if args.init:
        from .draft import draft_config
        try:
            print(draft_config([Path(p) for p in args.init]), end="")
        except (OSError, ValueError) as exc:
            print(f"init error: {exc}", file=sys.stderr)
            return 2
        return 0
    if not args.config:
        print("config error: --config is required unless --init is given",
              file=sys.stderr)
        return 2
    if args.seed is not None and args.seed < 0:
        print(f"config error: --seed must be a non-negative integer, got {args.seed}",
              file=sys.stderr)
        return 2
    if args.threads is not None and args.threads < 1:
        print(f"config error: --threads must be at least 1, got {args.threads}",
              file=sys.stderr)
        return 2
    if args.threads:
        _set_threads(args.threads)

    # NumPy/SciPy read BLAS/OpenMP limits while their extension modules are
    # imported. Importing ``core`` at module load made ``--threads`` cosmetic:
    # by the time the variables were set, the thread pools already existed.
    # Delay the numerical stack until after the execution contract is applied.
    from . import core
    from .jsonio import dumps
    from .outputs import validate_destination, publish_input_check

    try:
        cfg = load_config(args.config, check_files_exist=not args.no_check_files)
        cfg = _with_overrides(cfg, args)
    except ConfigError as exc:
        print(f"config error: {exc}", file=sys.stderr)
        return 2
    if args.dry_run:
        import yaml
        print(yaml.safe_dump(core._config_record(cfg), sort_keys=False, allow_unicode=True), end="")
        return 0

    rd = Path(args.results_dir).expanduser().resolve() if args.results_dir else None
    if rd is not None:
        try:
            validate_destination(rd, overwrite=args.overwrite)
        except OSError as exc:
            print(f"output error: {exc}", file=sys.stderr)
            return 2
    def progress(message):
        if not args.quiet:
            print(f'[amr-clonalshare] {message}', file=sys.stderr, flush=True)
    log: list = []
    try:
        if args.check_input:
            from .io import load_dataset
            from .qc import render_markdown
            ds = load_dataset(cfg)
            assert ds.input_qc is not None
            if rd is not None:
                publish_input_check(ds.input_qc, rd, overwrite=args.overwrite)
            print(render_markdown(ds.input_qc))
            return 0
        with core.log_warnings(log):
            record = core.run(cfg, results_dir=rd, seed=args.seed,
                              overwrite=args.overwrite, progress=progress,
                              check_files_exist=not args.no_check_files, log=log)
        summary = _summary(record)
    except ConfigError as exc:
        print(f'input error: {exc}', file=sys.stderr)
        _failure_log(rd, log, f'input error: {exc}')
        return 2
    except (OSError, ValueError) as exc:
        print(f'run error: {exc}', file=sys.stderr)
        _failure_log(rd, log, f'run error: {exc}')
        return 2
    except Exception as exc:
        print(f'run error: {type(exc).__name__}: {exc}', file=sys.stderr)
        _failure_log(rd, log, f'run error: {type(exc).__name__}: {exc}')
        return 2
    if not args.quiet:
        if args.json:
            print(dumps(summary))
        else:
            print(results_table(record))
        if rd is not None:
            print(f"Report: {rd / 'report.html'}")
            print(f"Results: {rd / 'results.csv'}; record: {rd / 'clonal_share_result.json'}")
        else:
            print("Nothing was written: add --results-dir <folder> to keep the record, the results table and the reports.")
    return 0


def _with_overrides(cfg, args):
    """The configuration with the budgets and the level the command line
    replaces, each checked as the file's value would be."""
    from dataclasses import replace
    changes = {}
    if args.permutations is not None:
        if args.permutations < 1:
            raise ConfigError(f"--permutations must be at least 1, got {args.permutations}")
        changes["attribution"] = replace(cfg.attribution, n_perm=args.permutations)
    if args.bootstrap is not None:
        if args.bootstrap < 0:
            raise ConfigError(f"--bootstrap must be 0 or more, got {args.bootstrap}")
        changes["attribution"] = replace(changes.get("attribution", cfg.attribution), n_boot=args.bootstrap)
    if args.alpha is not None:
        if not 0 < args.alpha < 1:
            raise ConfigError(f"--alpha must lie between 0 and 1, got {args.alpha}")
        changes["evidence"] = replace(cfg.evidence, alpha=args.alpha)
    return replace(cfg, **changes) if changes else cfg


def _subcommand(name: str, rest: list) -> int:
    """The subcommands, each spelt out as the options of the main command
    or handed to the command it stands for."""
    if name == "run":
        return main(rest)
    if name == "init":
        return main(["--init", *rest]) if rest else main(["--init"])
    if name == "check":
        return main(["--check-input", *rest])
    if name == "compare":
        from .comparison import main as compare_main
        return compare_main(rest)
    if name == "gui":
        from .gui import main as gui_main
        return gui_main(rest)
    if name == "version":
        print(f"amr-clonalshare {__version__}")
        return 0
    if name == "doctor":
        print(doctor_report())
        return 0
    if name == "completion":
        shell = rest[0] if rest else ""
        if shell not in ("bash", "powershell"):
            print("completion: name the shell, bash or powershell", file=sys.stderr)
            return 2
        print(completion_script(shell), end="")
        return 0
    return _example(rest)


def examples_folder() -> Optional[Path]:
    """The examples folder of a source checkout or of the Windows folder,
    found beside the package or named by AMR_CLONALSHARE_EXAMPLES."""
    named = os.environ.get("AMR_CLONALSHARE_EXAMPLES")
    candidates = [Path(named)] if named else []
    here = Path(__file__).resolve()
    candidates += [here.parents[2] / "examples", here.parents[3] / "examples", Path.cwd() / "examples"]
    for folder in candidates:
        if (folder / "workflows" / "mic.yaml").is_file():
            return folder
    return None


def _example(rest: list) -> int:
    parser = argparse.ArgumentParser(prog="amr-clonalshare example",
                                     description="Run a shipped collection with the seed of its published record.")
    parser.add_argument("name", nargs="?", choices=sorted(EXAMPLES), help="the collection")
    parser.add_argument("--results-dir", required=False)
    parser.add_argument("--seed", type=int, default=None, help="default: the seed of the published record")
    parser.add_argument("--threads", type=int, default=None)
    parser.add_argument("--check-input", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--json", action="store_true")
    args, passed = parser.parse_known_args(rest)
    folder = examples_folder()
    if args.name is None:
        print("Shipped collections: " + ", ".join(sorted(EXAMPLES))
              + (f"; found under {folder}" if folder else "; no examples folder was found (set AMR_CLONALSHARE_EXAMPLES)"))
        return 0
    if folder is None:
        print("example error: no examples folder was found beside the package; set AMR_CLONALSHARE_EXAMPLES "
              "to the examples folder of a source checkout or of the Windows folder", file=sys.stderr)
        return 2
    relative, seed = EXAMPLES[args.name]
    argv = ["--config", str(folder / relative), "--seed", str(args.seed if args.seed is not None else seed)]
    if args.results_dir:
        argv += ["--results-dir", args.results_dir]
    if args.threads:
        argv += ["--threads", str(args.threads)]
    for flag in ("check_input", "overwrite", "quiet", "json"):
        if getattr(args, flag):
            argv.append("--" + flag.replace("_", "-"))
    return main(argv + passed)


def doctor_report() -> str:
    """What this computer offers the package: the versions it runs on, the
    threads, whether workbooks can be read, and where the examples are."""
    import platform
    lines = [f"amr-clonalshare {__version__}", f"Python {platform.python_version()} on {platform.platform()}"]
    for module in ("numpy", "pandas", "scipy", "yaml"):
        try:
            mod = __import__(module)
            lines.append(f"{module} {getattr(mod, '__version__', 'present')}")
        except ImportError:
            lines.append(f"{module} missing")
    try:
        import openpyxl  # noqa: F401
        lines.append("Excel workbooks: readable (openpyxl present)")
    except ImportError:
        lines.append("Excel workbooks: not readable; pip install 'amr-clonalshare[excel]' or save the sheet as CSV")
    lines.append(f"Processor threads: {os.cpu_count() or 'unknown'}"
                 + (f" (OMP_NUM_THREADS={os.environ['OMP_NUM_THREADS']})" if os.environ.get("OMP_NUM_THREADS") else ""))
    folder = examples_folder()
    lines.append(f"Examples: {folder}" if folder else "Examples: not found beside the package")
    return "\n".join(lines)


def completion_script(shell: str) -> str:
    """A completion script for bash or PowerShell, listing the subcommands,
    the options of the main command and the shipped collections."""
    options = sorted(a for action in build_parser()._actions for a in action.option_strings)
    words = " ".join([*SUBCOMMANDS, *options, *sorted(EXAMPLES)])
    if shell == "bash":
        return ("_amr_clonalshare() {\n"
                "    local cur=${COMP_WORDS[COMP_CWORD]}\n"
                f"    COMPREPLY=( $(compgen -W \"{words}\" -- \"$cur\") )\n"
                "}\n"
                "complete -o default -F _amr_clonalshare amr-clonalshare\n")
    return ("Register-ArgumentCompleter -Native -CommandName amr-clonalshare -ScriptBlock {\n"
            "    param($wordToComplete, $commandAst, $cursorPosition)\n"
            f"    '{words}'.Split(' ') | Where-Object {{ $_ -like \"$wordToComplete*\" }} |\n"
            "        ForEach-Object { [System.Management.Automation.CompletionResult]::new($_, $_, 'ParameterValue', $_) }\n"
            "}\n")


def _failure_log(rd, log, error: str) -> None:
    """After a failed run, the log with its error is left in the results
    folder (and nowhere else), so that one file says what happened."""
    if rd is not None:
        import traceback
        from .outputs import write_failure_log
        write_failure_log(rd, [*log, traceback.format_exc().rstrip()], error)


def results_table(record) -> str:
    """The results of a run as the form shows them: one line per agent with
    the share, its interval, the p- and q-values, the bounds and the
    conclusion."""
    from .verdict import VERDICT_NOTE, short_label, summary_rows
    rows = summary_rows(record)
    if not rows:
        return "No lineage share was computed; the input check and the report say why."
    with_calls = any("call_share" in r for r in rows)
    with_orders = any("mic_share" in r for r in rows)
    columns = [("Antimicrobial", "agent")]
    if with_calls:
        columns += [("Lineage share of the call", "call_share"), ("95% CI", "call_interval"), ("p", "call_p"),
                    ("q", "call_q"), ("Lower bound", "call_lower_bound"),
                    ("Lower confidence limit", "call_lower_limit"), ("Conclusion", "call_conclusion")]
    if with_orders:
        columns += [("Lineage share of the MIC ordering", "mic_share"), ("95% CI", "mic_interval"), ("p", "mic_p"),
                    ("q", "mic_q"), ("Lower bound", "lower_bound"), ("Lower confidence limit", "lower_limit"),
                    ("Conclusion", "mic_conclusion")]

    def cell(row, key):
        value = row.get(key, "")
        return short_label(value["label"]) if isinstance(value, dict) else str(value or "")

    table = [[head for head, _ in columns]] + [[cell(r, key) for _, key in columns] for r in rows]
    widths = [max(len(line[i]) for line in table) for i in range(len(columns))]
    lines = ["  ".join(value.ljust(widths[i]) if i == 0 or columns[i][1].endswith("conclusion") else value.rjust(widths[i])
                       for i, value in enumerate(line)).rstrip() for line in table]
    lines.insert(1, "  ".join("-" * w for w in widths))
    lines.append("")
    for r in rows:
        for key in ("call_conclusion", "mic_conclusion"):
            if r.get(key):
                lines.append(f"{r['agent']}: {r[key]['text']}")
    lines.append("")
    lines.append(VERDICT_NOTE)
    return "\n".join(lines)


if __name__ == "__main__":
    raise SystemExit(main())
