"""cli.py — command-line entry point.

    amr-clonalshare --config CONFIG [--results-dir DIR] [--seed N]
                    [--threads N] [--check-input] [--quiet] [--overwrite]
                    [--no-check-files]

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
from pathlib import Path
from typing import Any, Dict, Optional

from . import __version__
from .config import ConfigError, load_config


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="amr-clonalshare",
        description="Reads a phenotype and a lineage label per isolate and returns "
                    "the share of resistance variance that lineage membership "
                    "carries: a realised, a lineage-membership and a "
                    "superpopulation share with measured intervals, a gate before "
                    "each that refuses where the collection cannot identify the "
                    "quantity, a permuted-label control, an e-value per agent with "
                    "a sequential e-process over intakes, a Kitagawa decomposition "
                    "of a prevalence difference and an interval-censored reading of "
                    "dilution panels.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--version", action="version",
                   version=f"amr-clonalshare {__version__}")
    p.add_argument("--config", required=True, help="path to a dataset YAML config")
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
    p.add_argument("--quiet", action="store_true",
                   help="suppress the stdout summary of a run; the input check is still printed")
    p.add_argument("--overwrite", action="store_true",
                   help="replace an existing output bundle after successful completion; retain the previous directory as a backup")
    p.add_argument("--no-check-files", action="store_true",
                   help="skip the existence check on configured data files")
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
                        f"{' vs '.join(dec.get('levels') or [])}",
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
            "ci95": [d.get("ci_low"), d.get("ci_high")],
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


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
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
    except ConfigError as exc:
        print(f"config error: {exc}", file=sys.stderr)
        return 2

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
        record = core.run(cfg, results_dir=rd, seed=args.seed,
                          overwrite=args.overwrite, progress=progress,
                          check_files_exist=not args.no_check_files)
        summary = _summary(record)
    except ConfigError as exc:
        print(f'input error: {exc}', file=sys.stderr)
        return 2
    except (OSError, ValueError) as exc:
        print(f'run error: {exc}', file=sys.stderr)
        return 2
    if not args.quiet:
        print(dumps(summary))
        if rd is not None:
            print(f"[amr-clonalshare] outputs written to {rd}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
