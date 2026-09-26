"""core.py — the amr-clonalshare run.

The run, in order:

1. load the phenotype table, the metadata and, where one is named, the
   dilution table, keeping the union of recorded phenotype and MIC isolates;
2. write the input check: the lineage groups and their sizes, the support,
   and the count of the rarer outcome per antimicrobial;
3. per antimicrobial, the lineage-membership share against its permuted
   control (:mod:`.attribution`), the e-value and, where the intake column is
   declared, the sequential e-process (:mod:`.evalues`);
4. the reading of a recorded dilution as an interval (:mod:`.censored`);
5. prevalence per isolate and per lineage, the concentration of carriage, and,
   where two collections are declared, the decomposition of their prevalence
   difference (:mod:`.clonality`);
6. one versioned record, from which the API and command line render both
   reports, CSV and a completion manifest in a protected output directory.

Every estimator sits behind a gate that refuses where the collection cannot
identify the quantity, and the record carries the condition that failed.
"""
from __future__ import annotations

import math
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
from .io import _METADATA_MISSING, load_dataset
from ._seeding import estimator_rng as _estimator_rng
from .phenotype import _token

__all__ = ["run"]

# The share draws from the seventh of eight generators spawned off the master
# seed; the position is part of the record format, so a record reproduces at
# the same seed.
_LINEAGE_STREAM = 6
_EVIDENCE_STREAM = 7


def _lineage_rng(seed: int) -> np.random.Generator:
    children = np.random.SeedSequence(seed).spawn(8)
    return np.random.default_rng(children[_LINEAGE_STREAM])




def run(cfg: Config, *, results_dir=None, seed: int = 42, overwrite: bool = False,
        progress=None, check_files_exist: bool = True) -> dict:
    """Read a cohort for the share its lineages carry.

    Per antimicrobial, the lineage-membership share against its permuted
    control, the e-value and, where the intake column is declared, the
    sequential e-process, the reading of a recorded dilution, and the
    decomposition of a prevalence difference, each behind its own gate. A complete JSON, CSV, QC and report bundle is
    written when ``results_dir`` is given. An existing nonempty destination requires explicit overwrite.
    """
    cfg.validate(check_files_exist=check_files_exist)
    from .outputs import result_rows, validate_destination, publish_results
    from .missingness import finite_collection_bounds
    if results_dir is not None:
        validate_destination(results_dir, overwrite=overwrite)
    notify = progress or (lambda _message: None)
    notify('Reading input tables and checking data')
    rng = _lineage_rng(seed)
    ds = load_dataset(cfg)
    strain_ids = ds.strain_ids
    n_isolates = 0 if strain_ids is None else len(strain_ids)
    assert ds.panel is not None
    raw_panel = ds.panel.reindex(strain_ids)
    X_df = raw_panel.dropna(axis=1, how="all")
    notify('Computing analyses on their observed and typed subsets')

    meta_block = _metadata_diagnostics(cfg, ds, strain_ids, rng,
                                       X_df=X_df, seed=seed, progress=notify)
    out = {
        "schema_version": "1.0",
        "seed": int(seed),
        "config": _config_record(cfg),
        "versions": _versions(),
        "n_isolates": int(n_isolates),
        "n_traits": int(raw_panel.shape[1]),
        "traits": [str(c) for c in raw_panel.columns],
        "metadata_diagnostics": meta_block,
        "input_qc": ds.input_qc,
        "collection_bounds": {str(agent): finite_collection_bounds(raw_panel[agent].to_numpy(dtype=float))
                              for agent in raw_panel.columns},
    }
    out['results'] = result_rows(out)
    if cfg.dataset.stratify_by:
        notify('Repeating the analyses within each stratum')
        out['strata'] = _strata(cfg, ds, strain_ids, raw_panel, seed, notify)
    if results_dir is not None:
        notify('Writing reports, CSV and completion manifest')
        publish_results(out, results_dir, overwrite=overwrite)
    notify('Analysis complete')
    return out


def _strata(cfg, ds, strain_ids, raw_panel, seed, notify) -> dict:
    """Every analysis repeated on the isolates of each level of ``stratify_by``."""
    from .outputs import result_rows
    from .missingness import finite_collection_bounds
    column = cfg.dataset.stratify_by
    labels = ds.metadata[column].reindex(pd.Index(strain_ids).astype(str)) if ds.metadata is not None else None
    if labels is None:
        raise ValueError(f"stratify_by needs the metadata column {column!r}")
    labels = labels.astype(object).where(labels.notna(), None)
    levels = sorted({str(v) for v in labels if v is not None and str(v).strip() != ""})
    strata: Dict[str, Any] = {"column": column, "levels": {},
                              "n_without_level": int(sum(1 for v in labels if v is None or str(v).strip() == ""))}
    for level in levels:
        keep = np.asarray([v is not None and str(v) == level for v in labels], dtype=bool)
        ids = pd.Index(strain_ids)[keep]
        panel = raw_panel.loc[ids]
        notify(f'Stratum {column} = {level}: {len(ids)} isolates')
        block = _metadata_diagnostics(cfg, ds, ids, _lineage_rng(seed),
                                      X_df=panel.dropna(axis=1, how="all"), seed=seed)
        # The stratum record carries what result_rows needs to list every
        # requested analysis, including the ones this stratum could not run.
        record = {"schema_version": "1.0", "config": _config_record(cfg),
                  "n_isolates": int(len(ids)), "n_traits": int(panel.shape[1]),
                  "traits": [str(c) for c in panel.columns],
                  "metadata_diagnostics": block,
                  "collection_bounds": {str(agent): finite_collection_bounds(panel[agent].to_numpy(dtype=float))
                                        for agent in panel.columns}}
        record["results"] = result_rows(record)
        strata["levels"][level] = record
    return strata


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
                    "evidence", "population_probit"):
        value = getattr(cfg, section, None)
        if value is not None:
            record[section] = plain(value)
            if section == "population_probit" and record[section].get("interval_method") != "general":
                for name in ("memory_budget_mb", "table_cache_mb", "cache_dir"):
                    record[section].pop(name, None)
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
        keys = [float(v) for v in levels]
    except ValueError:
        return sorted(levels)
    if not all(math.isfinite(k) for k in keys):
        return sorted(levels)
    return sorted(levels, key=lambda v: (float(v), v))


def _lineage_selection(shares: Dict[str, Any], alpha: float) -> Dict[str, Any]:
    """Which antimicrobials show a lineage effect, from one analysis.

    The permutation p-values of the clonal share, one per antimicrobial,
    enter the Benjamini-Yekutieli step-up, which controls the false discovery
    rate whatever the dependence between agents (cross-resistance and
    trade-offs give it either sign). On one look this selects far more true
    effects than e-BH on the e-values at the same guarantee; the e-values stay
    in the record for a programme read after every intake, where only they
    keep their guarantee (``lineage_evidence``).
    """
    from .stats import benjamini_hochberg
    names = [k for k, v in shares.items() if np.isfinite(v.get("p_value", float("nan")))]
    if not names:
        return {"method": "benjamini_yekutieli", "alpha": float(alpha), "rejected_features": [],
                "n_tested": 0}
    p = np.array([shares[k]["p_value"] for k in names], dtype=float)
    adjusted, reject = benjamini_hochberg(p, float(alpha), dependence="arbitrary")
    floor = min(float(shares[k].get("p_floor", float("nan"))) for k in names)
    return {"method": "benjamini_yekutieli", "alpha": float(alpha), "n_tested": len(names),
            "q_values": {k: float(q) for k, q in zip(names, adjusted)},
            "rejected_features": [k for k, r in zip(names, reject) if r],
            "p_value_floor": floor,
            "note": ("permutation p-values of the clonal share; the smallest attainable p-value "
                     "is 1/(n_perm + 1), so a panel of m agents needs n_perm large enough that "
                     "this floor is below alpha/m for a single agent to be selectable")}


def _population_profiles(panel, raw, interval_method="fixed_cutoff", *, seed=42, compute=None, progress=None) -> dict:
    """Group each exact tested/typed subset without consuming a classical random stream."""
    fit: Any
    unavailable: Any
    if interval_method == "general":
        from .general_probit import general_probit_icc, _unavailable_general_result
        fit, unavailable = general_probit_icc, _unavailable_general_result
    elif interval_method == "fixed_cutoff":
        from .population_probit import population_probit_icc, _unavailable_result
        fit, unavailable = population_probit_icc, _unavailable_result
    else:
        raise ValueError("Unknown population interval method")

    labels = np.asarray(raw, dtype=object)
    typed = np.array([not _attribution._is_untyped(v) for v in labels], dtype=bool)
    output = {}
    for name in panel:
        if progress is not None:
            progress(f'Population model: {name}')
        y = panel[name].to_numpy(dtype=float)
        finite = np.isfinite(y)
        keep = finite & typed
        values = y[keep]
        codes = _attribution._codes(labels[keep]) if keep.any() else np.zeros(0, dtype=int)
        sizes = np.bincount(codes)
        n = int(keep.sum())
        repeated = sizes >= 2
        support = float(sizes[repeated].sum()/n) if n else 0.
        if not n:
            result = unavailable(0, 0, 0, 0., "no_retained_calls",
                                         "No finite result with a recorded lineage")
        elif not np.isin(values, (0., 1.)).all():
            result = unavailable(n, len(sizes), int(repeated.sum()), support,
                                         "invalid_binary_input", "Retained calls must be zero or one", sizes=sizes)
        else:
            counts = np.bincount(codes, weights=values).astype(int)
            try:
                if interval_method == "general":
                    result = fit(counts, sizes, seed=seed, case_key=str(name), compute=compute)
                else:
                    result = fit(counts, sizes)
            except ValueError as error:
                result = unavailable(n, len(sizes), int(repeated.sum()), support,
                                             "invalid_grouped_input", str(error), sizes=sizes)
        row = result.as_dict()
        row.update(n_dropped_non_finite=int((~finite).sum()),
                   n_dropped_untyped=int((finite & ~typed).sum()))
        output[str(name)] = row
    return output


def _metadata_diagnostics(cfg, ds, strain_ids, rng,
                          X_df=None, seed: int = 0, progress=None) -> Optional[dict]:
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
                    repeats=att.repeats, n_boot=att.n_boot, n_perm=att.n_perm,
                    seed=_estimator_rng(seed, _LINEAGE_STREAM)).as_dict()
                for feat in X_df.columns
            }
            out["lineage_selection"] = _lineage_selection(
                out["clonal_share"], getattr(getattr(cfg, "evidence", None), "alpha", 0.05))
            # The realised share of the same call, with the exact interval its
            # kurtosis gate admits or the reason it refuses: on a binary call
            # the gate closes wherever the prevalence is extreme, and the
            # record says so rather than leaving the estimator out.
            out["realised_share"] = {
                str(feat): _realised.realised_share(
                    X_df[feat].to_numpy(dtype=float), raw).as_dict()
                for feat in X_df.columns
            }
        profile = getattr(cfg, "population_probit", None)
        if (profile is not None and profile.enabled and X_df is not None):
            if profile.interval_method == "general":
                from .general_probit import GeneralComputeOptions, GeneralComputeSession
                compute = GeneralComputeSession(GeneralComputeOptions(
                    memory_budget_mb=profile.memory_budget_mb,
                    table_cache_mb=profile.table_cache_mb, cache_dir=profile.cache_dir))
                out["population_probit_icc"] = _population_profiles(
                    X_df, raw, profile.interval_method, seed=seed, compute=compute, progress=progress)
                compute.clear()
            else:
                out["population_probit_icc"] = _population_profiles(X_df, raw, profile.interval_method, progress=progress)
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
                    repeats=ev.repeats,
                    seed=_estimator_rng(seed, _EVIDENCE_STREAM)).as_dict()
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
                    "note": ("Each running product requires a conditionally independent common-probability "
                             "Bernoulli null within each prespecified batch. The final-look panel decision "
                             "is taken on the log scale. A common stopping-time panel claim additionally "
                             "requires validity in the joint panel filtration; repeated unions of "
                             "rejection sets are not controlled."),
                }
        cen = getattr(cfg, "censored", None)
        if cen is not None and getattr(cen, "enabled", True):
            found = _censored_diagnostics(cfg, ds, strain_ids, raw, cen, seed)
            if found:
                out["censored_share"] = found
        surv = getattr(cfg, "surveillance", None)
        n_boot_s = getattr(surv, "n_boot", 2000) if surv else 2000
        surv_on = getattr(surv, "enabled", True) if surv else True
        cc = cfg.dataset.contrast_column
        if (surv_on and X_df is not None and len(X_df.columns) and cc
                and cc in md.columns
                and len(cfg.dataset.contrast_levels) == 2):
            # A prevalence difference between two collections splits into a
            # change in the lineage mix and a change in the within-lineage
            # rate. These are descriptive components; the marginal prevalence
            # alone cannot separate them or identify intervention effects.
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
                min_shared_support=(getattr(surv, "min_shared_support", 0.9)
                                    if surv else 0.9),
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


def _panel_cut_points(lo, hi, wells, end_wells_censored):
    """The log2 cut points of the dilution panel behind interval readings.

    With recorded wells every tested dilution is a cut point except the
    highest, whose reading is right-censored under the end-well assumption;
    in the sensitivity arm the lowest well is one doubling wide, so its lower
    edge is a cut point too. Without recorded wells the cut points are the
    finite endpoints of the readings, and the unobserved wells between them
    are filled in when every endpoint is a whole log2 step.
    """
    exact = lo == hi
    ends = np.unique(np.r_[lo[~exact], hi[~exact]])
    ends = ends[np.isfinite(ends)]
    if wells is not None:
        lattice = np.log2(np.sort(np.asarray(wells, dtype=float)))
        cuts = lattice[:-1] if end_wells_censored else np.r_[lattice[0] - 1., lattice]
        return np.unique(np.r_[cuts, ends])
    if ends.size and np.allclose(ends, np.round(ends)):
        return np.arange(ends.min(), ends.max() + .5, 1.)
    return ends


def _calibrated_mic_interval(lo, hi, lineage, cen, seed, wells=None, panel=None,
                             covariate=None, n_boot=None) -> Dict[str, Any]:
    """The validated population interval for one agent, or the reason it is absent.

    Exact readings use the exact generalized F pivot. Interval readings use
    the null-wise parametric bootstrap on the panel: the configured wells when
    the laboratory recorded them, otherwise the finite endpoints of the
    readings, with the unobserved wells between them filled in when every
    endpoint is a whole log2 step. With a ``panel`` label per reading the cut
    points are found within each panel and every simulated reading is made on
    the panel of the record it stands for. With a ``covariate`` level per
    reading, its fixed effects are part of the model; exact readings then use
    the bootstrap as well, since Wald's pivot has no fixed effects.
    """
    import os
    from .mic_inference import exact_gaussian_interval
    n_boot = cen.calibration_n_boot if n_boot is None else int(n_boot)
    from .mic_null_bootstrap import null_bootstrap_interval
    from ._mic_likelihood import GaussianMICLikelihood
    keep = ~(np.isnan(lo) | np.isnan(hi)) & ~(np.isneginf(lo) & np.isposinf(hi))
    lo, hi, lineage = lo[keep], hi[keep], lineage[keep]
    if panel is not None:
        panel = np.asarray([str(x) for x in panel])[keep]
    if covariate is not None:
        covariate = [np.asarray([str(x) for x in column])[keep] for column in covariate]
    try:
        if np.array_equal(lo, hi) and covariate is None:
            r = exact_gaussian_interval(lo, lineage)
            estimate = GaussianMICLikelihood(lo, hi, lineage).fit().rho
            return {"estimate": float(estimate), "low": r.low, "high": r.high,
                    "confidence": r.confidence, "method": r.method,
                    "status": r.status, "estimable": True, "n": r.n,
                    "n_groups": r.groups, "raw_set_empty": r.raw_set_empty}
        ends: Any
        if np.array_equal(lo, hi):
            # exact readings with a covariate: the panel plays no part in the
            # simulation, which reproduces the exact readings, but must be declared
            span = np.array([float(lo.min()), float(lo.max())])
            ends = span if panel is None else {str(name): span for name in np.unique(panel)}
        elif panel is None:
            ends = _panel_cut_points(lo, hi, wells, cen.end_wells_censored)
        else:
            ends = {str(name): _panel_cut_points(lo[panel == name], hi[panel == name], wells,
                                                 cen.end_wells_censored)
                    for name in np.unique(panel)}
        if cen.workers:
            workers = cen.workers
        elif hasattr(os, "sched_getaffinity"):
            workers = len(os.sched_getaffinity(0))  # CPUs this process may use
        else:
            workers = os.cpu_count() or 1
        r = null_bootstrap_interval(lo, hi, lineage, panel_edges=ends,
                                    n_boot=n_boot,
                                    seed=int(seed), workers=workers, panel=panel,
                                    covariate=covariate)
        adjustment: Dict[str, Any] = {}
        if covariate is not None:
            adjustment = {"adjustment": []}
            start = 0
            for column in covariate:
                levels = sorted(set(column))
                adjustment["adjustment"].append(
                    {"levels": levels, "reference_level": levels[0],
                     "coefficients": {name: float(c) for name, c in
                                      zip(levels[1:], r.fit.coefficients[start:start+len(levels)-1])}})
                start += len(levels)-1
        return {"estimate": r.fit.rho, "low": r.low, "high": r.high, **adjustment,
                "confidence": r.confidence, "method": r.method,
                "status": r.status, "estimable": bool(r.reportable),
                "n": int(lo.size), "n_groups": int(np.unique(lineage).size),
                "n_boot": n_boot, "seed": int(seed),
                "panel_edges": ([float(e) for e in ends] if panel is None else
                                {name: [float(e) for e in edges] for name, edges in ends.items()}),
                "tested_values": len(r.tests),
                "bootstrap_failed": int(sum(t.bootstrap_failed for t in r.tests)),
                "bootstrap_redrawn": int(sum(t.bootstrap_redrawn for t in r.tests)),
                "shape_check": r.shape_check,
                "reason": "" if r.reportable else r.status}
    except (ValueError, ArithmeticError) as exc:
        return {"estimable": False, "status": "not_computed", "reason": str(exc)}


def _dilution_table(lo, hi, lineage) -> Dict[str, Any]:
    """Counts of readings per lineage and dilution interval, for the report's map.

    Each column is one interval of the log2 scale as read, so an end well
    appears as an open interval. The columns are ordered by the edge that
    names them: a well by its upper edge, a left-open interval just before
    the well that shares its edge, a right-open interval just after the well
    it starts from, so that the open intervals of two laboratories with
    different end wells each sit where their edge is.
    """
    def place(pair):
        a, b = pair
        if np.isneginf(a):
            return (b, 0)
        if np.isposinf(b):
            return (a + .5, 1)
        return (b, 1)
    keep = ~(np.isnan(lo) | np.isnan(hi))
    pairs = sorted({(float(a), float(b)) for a, b in zip(lo[keep], hi[keep])}, key=place)
    index = {ab: k for k, ab in enumerate(pairs)}
    names = sorted({str(x) for x in lineage[keep]})
    counts = np.zeros((len(names), len(pairs)), dtype=int)
    row = {name: k for k, name in enumerate(names)}
    for a, b, name in zip(lo[keep], hi[keep], lineage[keep]):
        counts[row[str(name)], index[(float(a), float(b))]] += 1
    return {"lineages": names,
            "intervals": [[None if not np.isfinite(a) else a, None if not np.isfinite(b) else b] for a, b in pairs],
            "counts": counts.tolist()}


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
    """
    mic = getattr(ds, "mic", None)
    if mic is None:
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
    pc = cfg.dataset.mic_panel_column
    ccs = cfg.dataset.mic_covariate_columns
    by_isolate: Dict[str, pd.Series] = {}
    for column in (pc, *ccs):
        if column and column not in mic.columns:
            # a property of the isolate's testing, recorded once per isolate
            series = ds.metadata[column].astype(object)
            series.index = series.index.astype(str)
            by_isolate[column] = series

    join = getattr(ds, "mic_join", None)
    if join is not None and join.get("strains_aligned") != len(ids):
        # a stratum: the join report describes the whole collection, so only
        # the counts that hold for these isolates are kept
        with_mic = int(pd.Index(ids).isin(mic[idc].astype(str)).sum())
        join = {"path": join.get("path"), "antimicrobials": join.get("antimicrobials"),
                "strains_with_mic": with_mic, "strains_aligned": len(ids),
                "join_rate": with_mic / len(ids) if len(ids) else 0.,
                "mic_wells": join.get("mic_wells"), "mic_units": join.get("mic_units"),
                "scope": "stratum; duplicate and parsing counts are in the pooled record"}
    out: Dict[str, Any] = {"join": dict(join) if join is not None else None,
                              "per_agent": {}}
    agents = sorted(set(mic[abc].astype(str)) | set((join or {}).get('antimicrobials') or []))
    for agent in agents:
        block = mic.loc[mic[abc].astype(str).eq(agent)].copy()
        # Duplicate resolution and its provenance belong to the input boundary.
        # Never silently choose another record or replace the loader's counts.
        if block[idc].duplicated().any():
            raise ValueError('MIC input boundary left unresolved duplicate records')
        block = block.set_index(block[idc].astype(str))
        values = pd.to_numeric(block[vc], errors="coerce").reindex(ids)
        present = values.notna().to_numpy()
        typed = np.array([not _attribution._is_untyped(v) for v in lineage], dtype=bool)
        n_untyped = int((present & ~typed).sum())
        present = present & typed
        wells = cfg.dataset.mic_wells.get(str(agent))
        units = cfg.dataset.mic_units.get(str(agent))
        if present.sum() < 8:
            # Every agent that reaches the estimator leaves a row: an agent
            # that vanished from per_agent could not be told from one the
            # panel never carried.
            out["per_agent"][str(agent)] = {
                "n": int(present.sum()),
                "n_dropped_untyped": n_untyped,
                "mic_units": units,
                "panel_source": "configured" if wells is not None else "inferred",
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
        def per_reading(column, what):
            if not column:
                return None
            labels = (by_isolate[column] if column in by_isolate else block[column]).reindex(ids)
            missing = labels.isna() | labels.map(_token).isin(_METADATA_MISSING)
            if (missing.to_numpy() & present).any():
                raise ValueError(f"MIC readings of {agent!r} without a value in the "
                                 f"{what} column {column!r}; every reading needs one")
            return labels.astype(str).to_numpy(dtype=object)[present]
        pan = per_reading(pc, "panel")
        cov = [per_reading(c, "covariate") for c in ccs] or None
        if pan is None:
            geometry = _censored.panel_geometry(v, wells=wells).as_dict()
            lo, hi = _censored.intervals_from_mic(
                v, operators=ops,
                treat_end_wells_as_censored=cen.end_wells_censored, wells=wells)
        else:
            geometry = {"per_panel": {str(name): _censored.panel_geometry(v[pan == name], wells=wells).as_dict()
                                      for name in np.unique(pan)},
                        "panel_column": pc}
            lo, hi = _censored.intervals_by_panel(
                v, pan, operators=ops,
                treat_end_wells_as_censored=cen.end_wells_censored, wells=wells)
        try:
            share_record = _censored.censored_clonal_share(lo, hi, lin).as_dict()
        except ArithmeticError as exc:
            # A numerical failure for one agent is recorded for that agent and
            # does not end the run.
            share_record = {"estimable": False, "status": "not_computed",
                            "reason": f"numerical failure: {exc}"}
        entry: Dict[str, Any] = {
            "n": int(present.sum()),
            "n_dropped_untyped": n_untyped,
            "mic_units": units,
            "panel_source": "configured" if wells is not None else "inferred",
            "panel": geometry,
            "dilution_table": _dilution_table(lo, hi, lin),
            "share": share_record,
        }
        if cen.calibrated_interval:
            entry["calibrated"] = _calibrated_mic_interval(lo, hi, lin, cen, seed, wells=wells, panel=pan,
                                                           covariate=cov, n_boot=cen.draws_for(agent))
            if ccs:
                entry["calibrated"]["covariate_columns"] = list(ccs)
                for column, block_ in zip(ccs, entry["calibrated"].get("adjustment") or []):
                    block_["column"] = column
        out["per_agent"][str(agent)] = entry
    return out if out["per_agent"] else None
