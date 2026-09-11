"""core.py — the amr-clonalshare run.

The run, in order:

1. load the phenotype table, the metadata and, where one is named, the
   dilution table, and align them on the isolates the phenotype table holds;
2. write the input check: the lineage groups and their sizes, the support,
   and the count of the rarer outcome per antimicrobial;
3. per antimicrobial, the lineage-membership share against its permuted
   control (:mod:`.attribution`), the e-value and, where the intake column is
   declared, the sequential e-process (:mod:`.evalues`);
4. the reading of a recorded dilution as an interval (:mod:`.censored`);
5. prevalence per isolate and per lineage, the concentration of carriage, and,
   where two collections are declared, the decomposition of their prevalence
   difference (:mod:`.clonality`);
6. the record, ``clonal_share_result.json``, from which the command line
   renders both reports.

Every estimator sits behind a gate that refuses where the collection cannot
identify the quantity, and the record carries the condition that failed.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd

from . import attribution as _attribution
from . import realised as _realised
from . import censored as _censored
from . import clonality as _clonality
from . import evalues as _evalues
from .config import Config
from .io import load_dataset

__all__ = ["run"]

# The share draws from the seventh of eight generators spawned off the master
# seed; the position is part of the record format, so a record reproduces at
# the same seed.
_LINEAGE_STREAM = 6


def _lineage_rng(seed: int) -> np.random.Generator:
    children = np.random.SeedSequence(seed).spawn(8)
    return np.random.default_rng(children[_LINEAGE_STREAM])


def run(cfg: Config, *, results_dir=None, seed: int = 42) -> dict:
    """Read a cohort for the share its lineages carry.

    Per antimicrobial, the lineage-membership share against its permuted
    control, the e-value and, where the intake column is declared, the
    sequential e-process, the reading of a recorded dilution, and the
    decomposition of a prevalence difference, each behind its own gate. The
    record is written to ``results_dir/clonal_share_result.json`` when a
    directory is given.
    """
    rng = _lineage_rng(seed)
    ds = load_dataset(cfg)
    strain_ids = ds.strain_ids
    n_isolates = 0 if strain_ids is None else len(strain_ids)
    if n_isolates < 4:
        raise ValueError(f"need at least 4 aligned isolates; got {n_isolates}")

    assert ds.panel is not None
    X_df = ds.panel.reindex(strain_ids).dropna(axis=1, how="all")

    meta_block = _metadata_diagnostics(cfg, ds, strain_ids, rng,
                                       X_df=X_df, seed=seed)
    out = {
        "schema_version": "1.0",
        "seed": int(seed),
        "config": _config_record(cfg),
        "versions": _versions(),
        "n_isolates": int(n_isolates),
        "n_traits": int(X_df.shape[1]),
        "traits": [str(c) for c in X_df.columns],
        "metadata_diagnostics": meta_block,
        "input_qc": ds.input_qc,
    }
    if results_dir is not None:
        out_dir = Path(results_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        from .jsonio import write_json
        write_json(out, out_dir / "clonal_share_result.json")
    return out


def _config_record(cfg: Config) -> dict:
    """The configuration as run, every section, so that the record alone
    says what produced it. The record carries no run identifier of its own;
    the report derives one by digesting this block."""
    from dataclasses import asdict, is_dataclass

    def plain(obj):
        if is_dataclass(obj):
            return {k: plain(v) for k, v in asdict(obj).items()}
        if isinstance(obj, Path):
            return str(obj)
        if isinstance(obj, (list, tuple)):
            return [plain(v) for v in obj]
        return obj

    record = {}
    for section in ("dataset", "attribution", "surveillance", "censored",
                    "evidence"):
        value = getattr(cfg, section, None)
        if value is not None:
            record[section] = plain(value)
    record["config_file"] = (cfg.config_path.name
                             if cfg.config_path is not None else None)
    return record


def _versions() -> dict:
    import platform

    import pandas
    import scipy

    from . import __version__
    return {"amr_clonalshare": __version__, "python": platform.python_version(),
            "numpy": np.__version__, "pandas": pandas.__version__,
            "scipy": scipy.__version__}


def _batch_order(levels) -> list:
    """Intake labels in arrival order.

    Labels that all read as numbers, years or intake counters, are ordered
    numerically, so that intake 10 follows intake 9 rather than intake 1;
    anything else is ordered as text. The order used is written into the
    record beside the product it produced.
    """
    levels = list(levels)
    try:
        return sorted(levels, key=lambda v: (float(v), v))
    except ValueError:
        return sorted(levels)


def _metadata_diagnostics(cfg, ds, strain_ids, rng,
                          X_df=None, seed: int = 0) -> Optional[dict]:
    md = getattr(ds, "metadata", None)
    if md is None or cfg.dataset.lineage_column is None:
        return None
    md = md.reindex(strain_ids)
    out: Dict[str, Any] = {}
    col = cfg.dataset.lineage_column
    if col in md.columns:
        out["lineage_column"] = col
        raw = md[col].tolist()
        att = getattr(cfg, "attribution", None)
        if (att is not None and getattr(att, "enabled", True)
                and X_df is not None and len(X_df.columns)):
            # The per-agent share is a property of the label alone: the
            # lineage means scored on isolates the estimator did not see,
            # against the same score of a permuted labelling.
            out["clonal_share"] = {
                str(feat): _attribution.clonal_share(
                    X_df[feat].to_numpy(dtype=float), raw, folds=att.folds,
                    repeats=att.repeats, n_boot=att.n_boot,
                    n_perm=att.n_perm, seed=rng).as_dict()
                for feat in X_df.columns
            }
            # The realised share of the same call, with the exact interval its
            # kurtosis gate admits or the reason it refuses: on a binary call
            # the gate closes wherever the prevalence is extreme, and the
            # record says so rather than leaving the estimator out.
            out["realised_share"] = {
                str(feat): _realised.realised_share(
                    X_df[feat].to_numpy(dtype=float), raw).as_dict()
                for feat in X_df.columns
            }
        # The share says how much. Whether there is anything there at all is a
        # separate question, and surveillance asks it again every year on the
        # same panel. A p-value recomputed at each of an unplanned number of
        # looks controls nothing; an e-value may be read as often as wanted.
        ev = getattr(cfg, "evidence", None)
        if (ev is not None and getattr(ev, "enabled", True)
                and X_df is not None and len(X_df.columns)):
            per_feature = {
                str(feat): _evalues.e_process(
                    X_df[feat].to_numpy(dtype=float), raw, folds=ev.folds,
                    repeats=ev.repeats, seed=rng).as_dict()
                for feat in X_df.columns
            }
            names = list(per_feature)
            decision = _evalues.e_bh(
                [per_feature[k]["e_value"] for k in names], alpha=ev.alpha)
            out["lineage_evidence"] = {
                "per_feature": per_feature,
                "e_bh": dict(decision, rejected_features=[
                    names[i] for i in decision["rejected"]]),
                "note": ("e-value in the betting sense of Vovk and Wang, not "
                         "the BLAST expectation value and not the E-value of "
                         "VanderWeele and Ding"),
            }
            # One look is one e-value. A programme of intakes is a product of
            # them, and only the product carries a guarantee that survives
            # being read after every intake. The batches are the levels of
            # dataset.batch_column in ascending order, which is why that
            # column has to sort into arrival order; an isolate with no
            # intake recorded is set aside and counted rather than assigned
            # to one.
            bc = getattr(cfg.dataset, "batch_column", None)
            if bc and bc in md.columns:
                stamp = md[bc]
                known = stamp.notna().to_numpy()
                levels = _batch_order({str(v) for v in stamp[known]})
                blocks = [(np.asarray(stamp.astype(str) == lv) & known)
                          for lv in levels]
                seq = {}
                for feat in X_df.columns:
                    y = X_df[feat].to_numpy(dtype=float)
                    seq[str(feat)] = _evalues.sequential_e_process(
                        [y[b] for b in blocks],
                        [np.asarray(raw, dtype=object)[b] for b in blocks],
                    ).as_dict()
                sq_names = list(seq)
                sq_decision = _evalues.e_bh_log(
                    [seq[k]["log_e"][-1] if seq[k]["log_e"] else float("-inf")
                     for k in sq_names], alpha=ev.alpha)
                out["lineage_evidence"]["sequential"] = {
                    "batch_column": bc,
                    "batches": levels,
                    "n_per_batch": [int(b.sum()) for b in blocks],
                    "n_without_batch": int((~known).sum()),
                    "per_feature": seq,
                    "e_bh": dict(sq_decision, rejected_features=[
                        sq_names[i] for i in sq_decision["rejected"]]),
                    "note": ("the running product over intakes is a test "
                             "supermartingale, so this decision holds at "
                             "whatever intake the programme is read at; the "
                             "panel decision is taken on the log scale, "
                             "because a product over many intakes overflows "
                             "the e-value long before it overflows its "
                             "logarithm"),
                }
        cen = getattr(cfg, "censored", None)
        if cen is not None and getattr(cen, "enabled", True):
            found = _censored_diagnostics(cfg, ds, strain_ids, raw, cen, seed)
            if found:
                out["censored_share"] = found
        # Prevalence per isolate and per lineage are different estimands and a
        # surveillance report quotes only the first. They separate exactly when
        # sampling across lineages is uneven, which is the normal state of a
        # submission-driven collection, so both are emitted per feature.
        surv = getattr(cfg, "surveillance", None)
        n_boot_s = getattr(surv, "n_boot", 2000) if surv else 2000
        n_perm_s = getattr(surv, "n_perm", 2000) if surv else 2000
        surv_on = getattr(surv, "enabled", True) if surv else True
        if surv_on and X_df is not None and len(X_df.columns):
            out["lineage_resolved_prevalence"] = {
                str(feat): _clonality.lineage_resolved_prevalence(
                    X_df[feat].to_numpy(), raw, n_boot=n_boot_s, rng=rng)
                for feat in X_df.columns
            }
            # How many lineages effectively carry each trait, and whether that
            # number departs from what the cohort's own lineage abundances
            # would produce by chance. Direction is reported because the
            # departure statistic fires on dispersion as readily as on
            # clonality.
            out["trait_concentration"] = {
                str(feat): _clonality.trait_concentration(
                    X_df[feat].to_numpy(), raw, n_perm=n_perm_s, rng=rng)
                for feat in X_df.columns
            }
        cc = cfg.dataset.contrast_column
        if (surv_on and X_df is not None and len(X_df.columns) and cc
                and cc in md.columns
                and len(cfg.dataset.contrast_levels) == 2):
            # A prevalence difference between two collections splits into a
            # change in the lineage mix and a change in the within-lineage
            # rate. The two imply opposite interventions, and the reported
            # prevalence cannot separate them.
            lv_a, lv_b = (str(v) for v in cfg.dataset.contrast_levels)
            side = md[cc].astype(str)
            in_a = (side == lv_a).to_numpy()
            in_b = (side == lv_b).to_numpy()
            lin_arr = np.asarray(raw, dtype=object)
            panel = _clonality.decompose_panel(
                X_df.loc[in_a], lin_arr[in_a],
                X_df.loc[in_b], lin_arr[in_b],
                n_boot=n_boot_s, rng=rng,
                q=getattr(surv, "q_fdr", 0.05) if surv else 0.05,
                min_shared_support=(getattr(surv, "min_shared_support", 0.8)
                                    if surv else 0.8),
                label_alpha=(getattr(surv, "label_alpha", 0.05)
                             if surv else 0.05))
            out["prevalence_decomposition"] = {
                "contrast_column": cc,
                "levels": [lv_a, lv_b],
                "n_a": int(in_a.sum()),
                "n_b": int(in_b.sum()),
                "family": panel["family"],
                "per_feature": panel["per_agent"],
            }
    return out or None


def _censored_diagnostics(cfg, ds, strain_ids, lineage, cen,
                          seed: int = 0) -> Optional[dict]:
    """The clonal share of each antimicrobial read from the dilution panel.

    A dichotomised call and a recorded MIC are the same likelihood at two
    interval widths, so this reports the same quantity as the binary share on
    a scale that has not thrown away where in the panel each isolate fell. The
    panel geometry is reported first and decides which modes the data support:
    a panel with fewer than three tested wells cannot carry an interval, and a
    panel with a fifth of its readings piled on an end well cannot be read as
    a set of point values.

    Every agent's cluster bootstrap is given the same integer seed, so the
    per-agent interval widths are common random numbers rather than
    independent draws, and two agents on the same cohort move together more
    than two independent runs would.
    """
    mic = getattr(ds, "mic", None)
    if mic is None or not len(mic):
        return None
    idc = cfg.dataset.mic_id_column
    abc = cfg.dataset.mic_antibiotic_column
    vc = cfg.dataset.mic_value_column
    opc = cfg.dataset.mic_operator_column
    if not opc and "_operator_from_value" in mic.columns:
        # the loader read a censoring sign off the value column
        opc = "_operator_from_value"
    ids = pd.Index(strain_ids).astype(str)
    lineage = pd.Series(list(lineage), index=ids)

    join = getattr(ds, "mic_join", None)
    duplicates: Dict[str, int] = {}
    out: Dict[str, Any] = {"join": dict(join) if join is not None else None,
                              "per_agent": {}}
    for agent, block in mic.groupby(mic[abc].astype(str), sort=True):
        rows_read = len(block)
        block = block.drop_duplicates(subset=[idc], keep="first")
        # A repeated (isolate, antimicrobial) reading is a property of the
        # supplied table, not of the cohort, so the second row is dropped and
        # the count is reported rather than the drop being silent.
        if rows_read > len(block):
            duplicates[str(agent)] = int(rows_read - len(block))
        block = block.set_index(block[idc].astype(str))
        values = pd.to_numeric(block[vc], errors="coerce").reindex(ids)
        present = values.notna().to_numpy()
        if present.sum() < 8:
            # Every agent that reaches the estimator leaves a row: an agent
            # that vanished from per_agent could not be told from one the
            # panel never carried.
            out["per_agent"][str(agent)] = {
                "n": int(present.sum()),
                "share": {
                    "estimable": False,
                    "reason": f"{int(present.sum())} of {len(ids)} aligned "
                              f"isolate(s) carry a reading for this agent, "
                              f"against the 8 the interval likelihood needs "
                              f"before a between-lineage component is "
                              f"identified",
                },
            }
            continue
        v = values.to_numpy(dtype=float)[present]
        lin = lineage.to_numpy(dtype=object)[present]
        ops = (block[opc].reindex(ids).to_numpy(dtype=object)[present]
               if opc else None)
        geometry = _censored.panel_geometry(v)
        lo, hi = _censored.intervals_from_mic(
            v, operators=ops,
            treat_end_wells_as_censored=cen.end_wells_censored)
        share = _censored.censored_clonal_share(lo, hi, lin,
                                                n_boot=cen.n_boot, seed=seed)
        entry: Dict[str, Any] = {
            "n": int(present.sum()),
            "panel": geometry.as_dict(),
            "share": share.as_dict(),
        }
        if cen.sensitivity and "point" not in geometry.admissible_modes:
            entry["end_well_sensitivity"] = _censored.sensitivity_endpoints(
                v, lin, operators=ops, n_boot=max(cen.n_boot // 4, 25), seed=seed)
        out["per_agent"][str(agent)] = entry
    if out["join"] is not None:
        out["join"]["n_duplicate_rows"] = int(sum(duplicates.values()))
        out["join"]["duplicate_rows_by_agent"] = duplicates
    return out if out["per_agent"] else None
