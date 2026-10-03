"""The amr-clonalshare run.

The run, in order:

1. load the phenotype table, the metadata and, where one is named, the
   dilution table, keeping the union of recorded phenotype and MIC isolates;
2. write the input check: the lineage groups and their sizes, the support,
   and the count of the rarer outcome per antimicrobial;
3. per antimicrobial, the lineage share of the call against its permuted
   control (:mod:`.attribution`), the e-value and, where the intake column is
   declared, the sequential e-process (:mod:`.evalues`);
4. the reading of a recorded dilution as an interval (:mod:`.censored`);
5. prevalence per isolate and per lineage, the concentration of carriage, and,
   where two collections are declared, the decomposition of their prevalence
   difference (:mod:`.clonality`);
6. one versioned record, from which the API and command line render both
   reports, CSV and a completion manifest in a protected output directory.

Every estimator has reporting conditions and refuses where the collection
cannot identify the quantity; the record names the condition that failed.
"""
from __future__ import annotations

import math
import time
import warnings
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd

from . import attribution as _attribution
from . import censored as _censored
from . import clonality as _clonality
from . import evalues as _evalues
from .config import Config
from .io import _METADATA_MISSING, load_dataset
from ._seeding import estimator_rng as _estimator_rng
from .phenotype import _token

__all__ = ["run"]

# Every estimator draws from its own numbered stream of the master seed
# (``estimator_rng``): the lineage share and the decomposition of a difference
# between two collections stream 6, the e-values 7, the MIC ordering 8 and the
# bounds of a call 9. The numbers are part of the record format, so a record
# reproduces at the same seed; ``_lineage_rng`` gives the draws of stream 6.
_LINEAGE_STREAM = 6
_EVIDENCE_STREAM = 7
_ORDER_STREAM = 8
_CALL_BOUNDS_STREAM = 9


def _lineage_rng(seed: int) -> np.random.Generator:
    children = np.random.SeedSequence(seed).spawn(8)
    return np.random.default_rng(children[_LINEAGE_STREAM])




def run(cfg: Config, *, results_dir=None, seed: int = 42, overwrite: bool = False,
        progress=None, check_files_exist: bool = True, log=None) -> dict:
    """Read a collection for the share its lineages account for.

    Per antimicrobial, the lineage share of the call against its permuted
    control, the e-value and, where the intake column is declared, the
    sequential e-process, the reading of a recorded dilution, and the
    decomposition of a prevalence difference, each with its own reporting conditions. A complete JSON, CSV, QC and report bundle is
    written when ``results_dir`` is given. An existing nonempty destination requires explicit overwrite.
    """
    cfg.validate(check_files_exist=check_files_exist)
    from .outputs import result_rows, validate_destination, publish_results
    from .missingness import finite_collection_bounds
    if results_dir is not None:
        validate_destination(results_dir, overwrite=overwrite)
    log = [] if log is None else log
    notify = _logged(progress, log)
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
        publish_results(out, results_dir, overwrite=overwrite, log=log)
    notify('Analysis complete')
    return out


def _logged(progress, log: list):
    """A progress callback that also keeps every line, with the time, for
    the run.log of the results folder; the warnings the run raises are kept
    in the same list through :func:`log_warnings`."""
    forward = progress or (lambda _message: None)

    def notify(message: str) -> None:
        log.append(f"{time.strftime('%Y-%m-%d %H:%M:%S')}  {message}")
        forward(message)
    return notify


class log_warnings:
    """Keep the warnings raised inside the block in ``log`` and let them
    reach the console as before."""

    def __init__(self, log: list) -> None:
        self.log = log
        self.previous = None

    def __enter__(self):
        self.previous = warnings.showwarning

        def record(message, category, filename, lineno, file=None, line=None):
            self.log.append(f"{time.strftime('%Y-%m-%d %H:%M:%S')}  warning: {message}")
            self.previous(message, category, filename, lineno, file, line)
        warnings.showwarning = record
        return self.log

    def __exit__(self, *exc):
        warnings.showwarning = self.previous
        return False


class _Pace:
    """Progress lines for a loop over agents: which agent, how many are
    done, and, once the first has shown how long one takes, about how long
    the rest will take."""

    def __init__(self, notify, label: str, n: int) -> None:
        self.notify = notify or (lambda _message: None)
        self.label, self.n, self.done = label, n, 0
        self.started = time.monotonic()

    def step(self, agent: str) -> None:
        self.done += 1
        line = f"{self.label}: {agent} ({self.done} of {self.n})"
        if self.done > 1:
            elapsed = time.monotonic() - self.started
            left = elapsed / (self.done - 1) * (self.n - self.done + 1)
            line += f"; {elapsed:.0f} s so far, about {max(left, 1):.0f} s to go"
        self.notify(line)


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
    for section in ("dataset", "attribution", "surveillance", "censored", "evidence"):
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
        keys = [float(v) for v in levels]
    except ValueError:
        return sorted(levels)
    if not all(math.isfinite(k) for k in keys):
        return sorted(levels)
    return sorted(levels, key=lambda v: (float(v), v))


def _lineage_selection(shares: Dict[str, Any], alpha: float) -> Dict[str, Any]:
    """Which antimicrobials show a lineage effect, from one analysis.

    The permutation p-values of the lineage share, one per antimicrobial,
    enter the Benjamini-Yekutieli step-up, which controls the false discovery
    rate whatever the dependence between agents (cross-resistance and
    trade-offs give it either sign). On one look this selects far more true
    effects than e-BH on the e-values at the same error rate; the e-values stay
    in the record for a programme read after every intake, where only they
    keep their error control (``lineage_evidence``).
    """
    from .stats import benjamini_hochberg
    names = [k for k, v in shares.items() if np.isfinite(v.get("p_value", float("nan")))]
    if not names:
        return {"method": "benjamini_yekutieli", "alpha": float(alpha), "rejected_features": [],
                "n_tested": 0}
    p = np.array([shares[k]["p_value"] for k in names], dtype=float)
    adjusted, reject = benjamini_hochberg(p, float(alpha), dependence="arbitrary")
    floor = min(float(shares[k].get("p_floor", float("nan"))) for k in names)
    out = {"method": "benjamini_yekutieli", "alpha": float(alpha), "n_tested": len(names),
           "q_values": {k: float(q) for k, q in zip(names, adjusted)},
           "rejected_features": [k for k, r in zip(names, reject) if r],
           "p_value_floor": floor,
           "note": ("permutation p-values of the lineage share; the smallest attainable p-value "
                    "is 1/(n_perm + 1), so a panel of m agents needs n_perm large enough that "
                    "this floor is below alpha/m for a single agent to be selectable")}
    needed = permutations_for_one_selection(len(names), float(alpha))
    out["permutations_for_one_selection"] = needed
    if np.isfinite(floor) and floor > 0 and round(1. / floor) - 1 < needed:
        out["resolution_note"] = (f"with {int(round(1. / floor)) - 1} permutations one agent alone cannot be "
                                  f"selected among {len(names)} at alpha = {alpha:g}; at least {needed} "
                                  "permutations would let a single very small p-value pass the step-up")
    return out


def permutations_for_one_selection(m: int, alpha: float) -> int:
    """The smallest number of permutations at which one agent of a panel of
    m can be selected by the Benjamini-Yekutieli step-up at level alpha when
    every other p-value is large: the attainable minimum 1/(B + 1) must not
    exceed alpha / (m H_m), so B >= ceil(m H_m / alpha) - 1."""
    if m < 1:
        return 0
    harmonic = sum(1. / k for k in range(1, m + 1))
    return int(math.ceil(m * harmonic / alpha)) - 1


#: Prefix of the unit given to an isolate whose unit is not recorded: a unit
#: of its own.
_OWN_UNIT = "__isolate__"


def _units(cfg, md) -> Optional[np.ndarray]:
    """The sampling unit of every isolate, or None when none is declared; an
    isolate without a recorded unit is a unit of its own."""
    column = getattr(cfg.dataset, "unit_column", None)
    if not column:
        return None
    if column not in md.columns:
        raise ValueError(f"unit_column needs the metadata column {column!r}")
    values = md[column]
    missing = values.isna() | values.map(_token).isin(_METADATA_MISSING)
    return np.asarray([f"{_OWN_UNIT}{i}" if m else str(v)
                       for i, (v, m) in enumerate(zip(values.tolist(), missing.tolist()))], dtype=object)


def _call_strata(cfg, md, X_df) -> Optional[np.ndarray]:
    """The stratum of every isolate for the analyses of a call, from
    ``dataset.phenotype_covariate_column``, or None when none is declared.
    Every isolate with a readable call needs a value; an isolate without a
    call may lack one, since no analysis of a call reads it."""
    column = getattr(cfg.dataset, "phenotype_covariate_column", None)
    if not column:
        return None
    if column not in md.columns:
        raise ValueError(f"phenotype_covariate_column needs the metadata column {column!r}")
    values = md[column]
    missing = (values.isna() | values.map(_token).isin(_METADATA_MISSING)).to_numpy()
    with_call = (X_df.notna().any(axis=1).to_numpy() if X_df is not None and len(X_df.columns)
                 else np.zeros(len(md), dtype=bool))
    if (missing & with_call).any():
        raise ValueError(f"{int((missing & with_call).sum())} isolate(s) with a call have no value in the "
                         f"phenotype covariate column {column!r}; every isolate with a call needs one")
    return np.asarray([str(v) for v in values.tolist()], dtype=object)


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
        units = _units(cfg, md)
        if units is not None:
            out["sampling_units"] = {"column": cfg.dataset.unit_column,
                                     "n_units": int(len(set(units.tolist()))),
                                     "n_isolates_without_unit": int(sum(
                                         str(u).startswith(_OWN_UNIT) for u in units))}
        # The strata of the analyses of a call: the permutation test exchanges
        # lineage labels within them, the interval draws each isolate's
        # stratum with its call, and the bounds a call places on the
        # latent ordering are the bounds within them. The e-values, which
        # test one common probability, take no strata.
        call_strata = _call_strata(cfg, md, X_df)
        if call_strata is not None:
            with_call = X_df.notna().any(axis=1).to_numpy() if X_df is not None and len(X_df.columns) else np.zeros(len(md), dtype=bool)
            out["call_strata"] = {"column": cfg.dataset.phenotype_covariate_column,
                                  "n_strata": int(len(set(call_strata[with_call].tolist())))}
        att = getattr(cfg, "attribution", None)
        if (att is not None and getattr(att, "enabled", True)
                and X_df is not None and len(X_df.columns)):
            # The per-agent share is a property of the label alone: the
            # lineage means scored on isolates the estimator did not see,
            # against the same score of a permuted labelling.
            notify = progress or (lambda _message: None)
            out["clonal_share"] = {}
            pace = _Pace(notify, "Lineage share of the call", len(X_df.columns))
            for feat in X_df.columns:
                pace.step(str(feat))
                out["clonal_share"][str(feat)] = _attribution.clonal_share(
                    X_df[feat].to_numpy(dtype=float), raw, folds=att.folds,
                    repeats=att.repeats, n_boot=att.n_boot, n_perm=att.n_perm,
                    strata=call_strata, units=units,
                    seed=_estimator_rng(seed, _LINEAGE_STREAM)).as_dict()
            out["lineage_selection"] = _lineage_selection(
                out["clonal_share"], getattr(getattr(cfg, "evidence", None), "alpha", 0.05))
            # A call is a reading with one cut point, so it bounds the lineage
            # share of the latent ordering it was cut from
            # (latent_order.py); these bounds contain those of any panel.
            from .latent_order import call_bounds
            for feat in X_df.columns:
                out["clonal_share"][str(feat)].update(call_bounds(
                    X_df[feat].to_numpy(dtype=float), raw, strata=call_strata,
                    seed=_estimator_rng(seed, _CALL_BOUNDS_STREAM)))
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
            # them, and only the product keeps error control that survives
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
                # The one rule here whose error is controlled over the whole
                # programme of looks: an agent is declared once its running
                # product reaches m / alpha at any intake.
                programme = _evalues.anytime_bonferroni([seq[k]["log_e"] for k in sq_names],
                                                        alpha=ev.alpha)
                programme["rejected_features"] = [sq_names[i] for i in programme["rejected"]]
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
                    "anytime_bonferroni": programme,
                    "note": ("Each running product requires a conditionally independent common-probability "
                             "Bernoulli null within each prespecified batch. The final-look panel decision "
                             "is taken on the log scale. A common stopping-time panel claim additionally "
                             "requires validity in the joint panel filtration; repeated unions of "
                             "rejection sets are not controlled."),
                }
        cen = getattr(cfg, "censored", None)
        if cen is not None and getattr(cen, "enabled", True):
            found = _censored_diagnostics(cfg, ds, strain_ids, raw, cen, seed, units=units, progress=progress)
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
                units_a=None if units is None else units[in_a],
                units_b=None if units is None else units[in_b],
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


def _panel_source(cfg: Config, agent: str, wells) -> str:
    """Where the tested concentrations of an agent came from: the
    configuration, a shipped preset (which the plate sheet should confirm),
    or the readings themselves."""
    if wells is None:
        return "inferred"
    preset = cfg.dataset.mic_panel_preset
    if preset:
        from .panels import PANEL_PRESETS
        if list(PANEL_PRESETS[preset]["wells"].get(agent, ())) == list(wells):
            return f"preset {preset}; verify against the plate sheet"
    return "configured"


def _censored_diagnostics(cfg, ds, strain_ids, lineage, cen,
                          seed: int = 0, units=None, progress=None) -> Optional[dict]:
    """The lineage share of the ordering of each antimicrobial's MIC readings.

    A call and a recorded MIC are one reading at two widths, so this reports
    the share of the call on a scale that keeps where in the panel each
    isolate fell: the share of the within-stratum mid-distribution scores
    (``mic_order.mic_order_share``) and the bounds on the share of the
    latent ordering (``latent_order``). The panel of every laboratory is
    reported first, with the share of readings on an end well.
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
    pace = _Pace(progress, "Lineage share of the MIC ordering", len(agents))
    for agent in agents:
        pace.step(agent)
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
        panel_source = _panel_source(cfg, str(agent), wells)
        # the configured unit, else the unit the loader read from the unit
        # column of the table
        concentration_units = cfg.dataset.mic_units.get(str(agent)) or ((join or {}).get("mic_units") or {}).get(str(agent))
        if present.sum() < 2:
            # Every agent that reaches the estimator leaves a row: an agent
            # that vanished from per_agent could not be told from one the
            # panel never carried.
            out["per_agent"][str(agent)] = {
                "n": int(present.sum()),
                "n_dropped_untyped": n_untyped,
                "mic_units": concentration_units,
                "panel_source": panel_source,
                "order": {
                    "estimable": False,
                    "reason": f"{int(present.sum())} of {len(ids)} aligned "
                              f"isolate(s) carry a typed reading for this agent; "
                              f"an ordering needs two",
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
        entry: Dict[str, Any] = {
            "n": int(present.sum()),
            "n_dropped_untyped": n_untyped,
            "mic_units": concentration_units,
            "panel_source": panel_source,
            "panel": geometry,
            "dilution_table": _dilution_table(lo, hi, lin),
        }
        att = getattr(cfg, "attribution", None)
        try:
            # The lineage share of the recorded ordering, scored within the
            # testing laboratory: it assumes no form for the readings, and is
            # the reading of the dilution scale the report leads with. The
            # declared covariates join the laboratory in the strata, so that a
            # shift between their levels is not read as a lineage difference.
            from .mic_order import mic_order_share
            columns = [(c, labels) for c, labels in zip((pc, *ccs), (pan, *(cov or []))) if labels is not None]
            strata = (np.array([" / ".join(str(v) for v in levels)
                                for levels in zip(*(labels for _, labels in columns))], dtype=object)
                      if columns else None)
            entry["order"] = mic_order_share(
                lo, hi, lin, strata=strata,
                units=None if units is None else np.asarray(units, dtype=object)[present],
                stratified_by=" × ".join(c for c, _ in columns) if columns else None,
                folds=getattr(att, "folds", 5), repeats=getattr(att, "repeats", 20),
                n_boot=getattr(att, "n_boot", 999), n_perm=getattr(att, "n_perm", 999),
                seed=_estimator_rng(seed, _ORDER_STREAM))
        except (ValueError, ArithmeticError) as exc:
            entry["order"] = {"estimable": False, "reason": str(exc)}
        out["per_agent"][str(agent)] = entry
    if out["per_agent"]:
        out["order_selection"] = _lineage_selection(
            {name: e["order"] for name, e in out["per_agent"].items() if "order" in e},
            getattr(getattr(cfg, "evidence", None), "alpha", 0.05))
        if out["order_selection"].get("resolution_note") and progress is not None:
            progress("Selection across agents: " + out["order_selection"]["resolution_note"])
    return out if out["per_agent"] else None
