"""The run report: one fixed structure, plain language, every number read from the record.

``cluster_result.json`` is the record. This module renders it for a reader who
will act on the result without reading the statistics: a laboratory head, a
surveillance officer, a reviewer of a veterinary programme. The section order
never changes, each section states what was computed and what the number
means, and no figure appears here that is not in the record it was read from.
"""
from __future__ import annotations

from typing import Dict, List, Optional

__all__ = ["render_report"]

_CAP = 12


def _f(x, nd=3) -> str:
    if x is None:
        return "not computed"
    try:
        if x != x:
            return "not computed"
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return str(x)


def _pct(x) -> str:
    return "not computed" if x is None or x != x else f"{100 * float(x):.1f} %"


def _ci(lo, hi, nd=3) -> str:
    if lo is None or hi is None or lo != lo or hi != hi:
        return ""
    return f" (95 % interval {_f(lo, nd)} to {_f(hi, nd)})"


def _structure(summary: dict) -> List[str]:
    k = summary.get("selected_k")
    n = summary.get("n_isolates")
    out = [f"The cohort holds {n} isolates. The procedure that chooses the number "
           f"of profile groups selected {k}."]
    if summary.get("no_structure") or k == 1:
        out.append("One group means the resistance profiles do not fall into "
                   "separable types: the differences between isolates are no "
                   "larger than a single population without subgroups would "
                   "produce. This is a result, not a failure.")
        return out
    status = summary.get("inference_status")
    if status == "withheld_inadequate_split_design":
        out.append("Whether these groups are real rather than an artefact of "
                   "the clustering could not be tested: the panel does not hold "
                   "enough independent features to split into a training and a "
                   "test half. The groups are descriptive only.")
    else:
        p = summary.get("p_value_structure_report")
        detected = summary.get("structure_detected")
        out.append(
            f"Tested on features held out of the clustering, the groups "
            f"{'are' if detected else 'are not'} reproducible"
            + (f" ({p})" if p else "") + ".")
        verdict = summary.get("discreteness_verdict")
        if verdict is not None:
            out.append(
                f"Discreteness verdict: {verdict}. "
                + ("The groups are separated by gaps; isolates do not shade "
                   "continuously from one type into the next."
                   if summary.get("discrete_beyond_a_gradient") else
                   "The groups sit on a gradient; they are convenient divisions "
                   "of a continuum rather than distinct types."))
    return out


def _clone(summary: dict, meta: dict) -> List[str]:
    lam = summary.get("lineage_attributable_share")
    if lam is None:
        return ["No lineage column was supplied, so the question of how much of "
                "the structure is inherited along lineages was not asked."]
    ci = summary.get("lineage_attributable_share_ci95") or [None, None]
    out = [f"Share of the profile groups explained by lineage: {_f(lam)}"
           f"{_ci(*ci)}. A value near 1 means the groups are the lineages under "
           "another name; a value near 0 means the groups cut across lineages."]
    att = meta.get("lineage_attribution") or {}
    if not _ci(*ci) and att.get("partition_refit"):
        out.append(
            "The groups were re-fitted inside every cross-validation fold, so a "
            "lineage bootstrap interval is not given for this share (a lineage "
            "drawn twice would train a group on rows that are also held out). "
            f"The spread of the share across repeated cross-validation is "
            f"{_f(att.get('lam_cv_sd'))}.")
    if summary.get("lineage_confounded"):
        out.append("This cohort is lineage-confounded: the groups should be read "
                   "as lineages, and any trait attributed to a group is first of "
                   "all a trait of that lineage.")
    resolved = summary.get("lineage_attribution_resolved_at_95ci")
    if resolved is False:
        out.append("The interval spans the decision threshold, so a larger cohort "
                   "could still change this reading.")
    share = summary.get("clonal_share_all_features")
    if share is not None:
        out.append(f"Clonal share of the whole panel, estimated out of sample: "
                   f"{_f(share)}. This is the fraction of the variation in the "
                   "panel that knowing an isolate's lineage predicts for an "
                   "isolate the model has not seen.")
    per = meta.get("clonal_share") or {}
    rows = [(f, r) for f, r in per.items() if isinstance(r, dict)]
    if rows:
        out.append("")
        out.append("Per trait, the clonal share with its interval, the share of "
                   "isolates in lineages large enough to inform it (support), and "
                   "whether the estimator accepts the cell:")
        out.append("")
        out.append("| trait | clonal share | 95 % interval | support | estimable |")
        out.append("|---|---|---|---|---|")
        rows.sort(key=lambda fr: -(fr[1].get("kappa_adj") or 0.0))
        for f, r in rows[:_CAP]:
            out.append(f"| {f} | {_f(r.get('kappa_adj'))} | "
                       f"{_f(r.get('ci_low'))} to {_f(r.get('ci_high'))} | "
                       f"{_pct(r.get('support'))} | "
                       f"{'yes' if r.get('estimable') else 'no'} |")
        if len(rows) > _CAP:
            out.append(f"| and {len(rows) - _CAP} more | | | | |")
    return out


def _surveillance(summary: dict) -> List[str]:
    s = summary.get("surveillance") or {}
    if not s:
        return []
    out = []
    dec = s.get("decomposition")
    if dec:
        out.append(
            f"Contrast {dec['contrast']}: of {dec['n_features']} traits, "
            f"{dec['n_composition_significant']} differ because the two "
            f"collections hold different lineages (a change in mix), and "
            f"{dec['n_within_lineage_significant']} differ because resistance "
            f"changed inside lineages (a change in rate).")
        if dec.get("n_offsetting"):
            out.append(
                f"{dec['n_offsetting']} trait(s) show a mix change and a rate "
                f"change of opposite sign that cancel: "
                f"{', '.join(dec['offsetting_features'][:_CAP])}. A prevalence "
                "table calls these stable. They are not.")
    widest = s.get("lineage_prevalence_widest_gap_feature")
    if widest:
        out.append(
            f"Largest gap between per-isolate and per-lineage prevalence: "
            f"{widest} ({_f(s.get('lineage_prevalence_widest_gap'))}). When the "
            "two differ, a few large lineages carry most of the resistance.")
    return out


def _gates(summary: dict, problems: List[str]) -> List[str]:
    out = []
    level = summary.get("claim_level")
    if level:
        out.append(f"Highest claim the diagnostics support: **{level}**.")
    codes = summary.get("active_gate_codes") or []
    if codes:
        out.append(f"Gates active on this run: {', '.join(str(c) for c in codes)}.")
    for p in problems:
        out.append(f"Diagnostic failure: {p}")
    if not codes and not problems:
        out.append("No diagnostic gate tripped.")
    out.append("A gate is not an error. It marks a reading the data cannot "
               "support, and the record keeps the numbers behind it so the "
               "limit can be examined.")
    return out


def render_report(result: dict, summary: dict, *, input_qc: Optional[dict] = None,
                  problems: Optional[List[str]] = None) -> str:
    """Markdown report from the run record. Fixed section order."""
    meta = result.get("metadata_diagnostics") or {}
    prov = result.get("provenance") or {}
    sections: List[tuple] = [
        ("1. What was analysed",
         [f"Layers: {', '.join(result.get('layers') or [])}. "
          f"Isolates: {result.get('n_isolates')}. Seed: {result.get('seed')}."
          + (f" Version: {prov.get('version')}." if prov.get("version") else "")]),
    ]
    if input_qc:
        m = input_qc.get("missing") or {}
        total = sum(v.get("cells", 0) for v in (m.get("layers") or {}).values())
        lin = input_qc.get("lineage")
        qc_lines = [
            f"Empty cells in the binary layers: {total}"
            + (f", handled by the policy `{m.get('policy')}`" if total else "")
            + ". Nothing was filled in."]
        if lin:
            qc_lines.append(
                f"Lineage groups: {lin['n_groups']}, of which {lin['n_singletons']} "
                f"hold a single isolate. Support {_pct(lin['support'])} against "
                f"the {_pct(lin['support_threshold'])} the estimator needs: "
                + ("accepted." if lin["estimable"] else "not accepted at this "
                   "typing resolution."))
        sections.append(("2. Input check", qc_lines))
    else:
        sections.append(("2. Input check", ["No input record was attached."]))
    sections.append(("3. Are there resistance-profile groups?", _structure(summary)))
    sections.append(("4. How much of it is the clone?", _clone(summary, meta)))
    surv = _surveillance(summary)
    if surv:
        sections.append(("5. Mix or rate: the surveillance reading", surv))
    sections.append((f"{6 if surv else 5}. What may be concluded",
                     _gates(summary, problems or [])))
    lines = ["# amr-clonalshare run report", ""]
    lines.append("Every number below is read from `cluster_result.json` and "
                 "`summary.json`, written by the same run.")
    lines.append("")
    for title, body in sections:
        lines.append(f"## {title}")
        lines.append("")
        lines.extend(body)
        lines.append("")
    return "\n".join(lines)
