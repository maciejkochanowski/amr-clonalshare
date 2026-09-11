"""cli.py — command-line entry point.

    amr-clonalshare --config CONFIG [--results-dir DIR] [--seed N]
                    [--threads N] [--check-input] [--quiet] [--no-check-files]

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
    p.add_argument("--quiet", action="store_true", help="suppress the stdout summary")
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
    lrp = meta.get("lineage_resolved_prevalence") or {}
    gaps = {f: r["difference_per_lineage_minus_per_isolate"]
            for f, r in lrp.items() if r.get("status") == "ok"}
    if gaps:
        widest = max(gaps, key=lambda f: abs(gaps[f]))
        out["lineage_prevalence_widest_gap_feature"] = widest
        out["lineage_prevalence_widest_gap"] = round(gaps[widest], 4)

    conc = meta.get("trait_concentration") or {}
    ok = [r for r in conc.values() if r.get("status") == "ok"]
    if ok:
        counts: Dict[str, int] = {}
        for r in ok:
            counts[r["direction"]] = counts.get(r["direction"], 0) + 1
        out["carriage_direction_counts"] = counts
        out["n_features_with_carriage_reading"] = len(ok)
        # `direction` compares the effective number of carrying lineages with
        # the null band; the p-value tests the size of the departure in bits.
        # They are different criteria and can disagree on the same trait, so
        # the direction is reported for the traits the p-value selects rather
        # than beside a count taken over all of them.
        departing = [r for r in ok if (r.get("p_value") or 1.0) < 0.05]
        out["n_features_departing_from_proportional_carriage"] = len(departing)
        dep_counts: Dict[str, int] = {}
        for r in departing:
            dep_counts[r["direction"]] = dep_counts.get(r["direction"], 0) + 1
        out["departing_carriage_direction_counts"] = dep_counts

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
    from .jsonio import dumps, write_json

    try:
        cfg = load_config(args.config, check_files_exist=not args.no_check_files)
    except ConfigError as exc:
        print(f"config error: {exc}", file=sys.stderr)
        return 2

    rd = Path(args.results_dir).expanduser().resolve() if args.results_dir else None
    if rd is not None:
        if rd.exists() and not rd.is_dir():
            print(f"config error: --results-dir {rd} exists and is not a directory",
                  file=sys.stderr)
            return 2
        try:
            rd.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            print(f"config error: cannot create --results-dir {rd}: {exc}",
                  file=sys.stderr)
            return 2

    # The input check runs before any estimate: a table that does not join,
    # or a lineage column that cannot support the estimator, is reported here
    # in plain language, and a refusal leaves with the config exit code rather
    # than a traceback.
    from .io import load_dataset
    from .qc import render_markdown
    try:
        ds = load_dataset(cfg)
    except ConfigError as exc:
        print(f"input error: {exc}", file=sys.stderr)
        return 2
    assert ds.input_qc is not None
    qc_text = render_markdown(ds.input_qc)
    if rd is not None:
        write_json(ds.input_qc, rd / "input_qc.json")
        (rd / "input_qc.md").write_text(qc_text, encoding="utf-8")
    if args.check_input:
        print(qc_text)
        return 0
    input_record = ds.input_qc
    del ds

    record = core.run(cfg, results_dir=rd, seed=args.seed)
    summary = _summary(record)
    if rd is not None:
        from .report import render_report
        from .report_html import render_html_report
        # Both reports are built from one structure, and both carry the
        # digest of the record file they describe, so a report and a record
        # can be matched without trusting either.
        record_bytes = (rd / "clonal_share_result.json").read_bytes()
        (rd / "report.md").write_text(
            render_report(record, summary, input_qc=input_record,
                          record_bytes=record_bytes),
            encoding="utf-8")
        (rd / "report.html").write_text(
            render_html_report(record, summary, input_qc=input_record,
                               record_bytes=record_bytes),
            encoding="utf-8")
    if not args.quiet:
        print(dumps(summary))
        if rd is not None:
            print(f"[amr-clonalshare] outputs written to {rd}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
