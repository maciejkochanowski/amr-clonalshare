"""The run report as a structure, before it is a page.

Why one model
-------------
The package writes the report twice, as prose that diffs cleanly between runs
and as a page with the diagnostics drawn. Two writers reading the same record
drift apart: one grew an input-check section the other never had. Building the
report once, as sections of typed blocks, and rendering that structure into
each form means the two cannot disagree about what a run found.

What the structure is shaped by
-------------------------------
A measurement report, not a pipeline log. The first section states the quantity
measured, the result with its interval, the conditions the estimate required
and whether they held, and what the run did not evaluate. Everything after it
is the evidence, in the order a reader asks for it: is the input admissible,
what is the share and against what null, what survives re-reading, what the
recorded resolution adds, whether a change between two collections is a change
of lineages or a change within them; then provenance.

What is never done here
-----------------------
Nothing is estimated. Every value is read from the record or from the summary
the command line derived from it, and every rule this module applies to those
values is stated in the text it produces. Interpretation is a separate block,
labelled as such, generated from fixed conditions on quantities the record
already holds; it never introduces a threshold the record does not carry. """
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Sequence, Tuple

from . import __version__
from .verdict import VERDICT_NOTE, agent_verdict, format_p, short_label

__all__ = ["Report", "Section", "Block", "build_report"]


# ----------------------------------------------------------------- structure
@dataclass
class Block:
    """One unit of a section.

    ``kind`` is one of ``p`` (a paragraph; ``text``), ``kv`` (label and value
    pairs; ``rows``), ``table`` (``caption``, ``headers`` as (label, align)
    pairs with align in {lab, num, flag}, ``rows``), ``figure`` (``figure`` is
    a name the renderer knows, ``data`` its input, ``caption``), ``callout``
    (``title`` and ``lines``). Text carries two inline marks only:
    ``**bold**`` and ```code```.
    """
    kind: str
    text: str = ""
    title: str = ""
    caption: str = ""
    headers: List[Tuple[str, str]] = field(default_factory=list)
    rows: List[Sequence[Any]] = field(default_factory=list)
    figure: str = ""
    data: Any = None
    lines: List[str] = field(default_factory=list)
    label: str = ""


@dataclass
class Section:
    key: str
    title: str
    blocks: List[Block] = field(default_factory=list)

    def p(self, text: str) -> None:
        self.blocks.append(Block("p", text=text))

    def kv(self, rows) -> None:
        self.blocks.append(Block("kv", rows=list(rows)))

    def table(self, label, caption, headers, rows) -> None:
        self.blocks.append(Block("table", label=label, caption=caption,
                                 headers=list(headers), rows=list(rows)))

    def figure(self, label, name, data, caption) -> None:
        self.blocks.append(Block("figure", label=label, figure=name, data=data,
                                 caption=caption))

    def callout(self, title, lines) -> None:
        self.blocks.append(Block("callout", title=title, lines=list(lines)))


@dataclass
class Report:
    ident: Dict[str, str]
    status: Dict[str, Any]
    sections: List[Section]
    symbols: List[Tuple[str, str]]
    colophon: List[str]


# ------------------------------------------------------------------ helpers
def finite(x) -> bool:
    return (isinstance(x, (int, float)) and not isinstance(x, bool)
            and math.isfinite(float(x)))


def num(x, nd: int = 3) -> str:
    if isinstance(x, bool):
        return "yes" if x else "no"
    if not finite(x):
        return "not computed"
    s = ("%%.%df" % nd) % float(x)
    # a value that rounds to zero is printed without a sign
    return s[1:] if s.startswith("-") and float(s) == 0 else s


def text(x) -> str:
    """A recorded value as text; a value the record does not hold is said so."""
    return "not recorded" if x is None else str(x)


def pct(x, nd: int = 1) -> str:
    return "not computed" if not finite(x) else num(100.0 * float(x), nd) + " %"


_SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


def evalue(x) -> str:
    """An e-value, compact: plain below a thousand, mantissa and power above."""
    if not finite(x):
        return "not computed"
    x = float(x)
    if x < 1000:
        return "%.1f" % x
    exp = int(math.floor(math.log10(x)))
    return "%.1f × 10%s" % (x / 10 ** exp, str(exp).translate(_SUP))


def interval(lo, hi, nd: int = 3) -> str:
    if not (finite(lo) and finite(hi)):
        return "not computed"
    return "%s to %s" % (num(lo, nd), num(hi, nd))


def _bounds(lo, hi, limit, exact=True) -> str:
    """The bounds the readings allow (an upper end that is only certified is
    marked ≤) and the lower confidence limit of the lower end."""
    if not (finite(lo) and finite(hi)):
        return "not computed"
    return "%s to %s%s; ≥ %s" % (num(lo), "" if exact else "≤ ", num(hi), num(limit))


def _join(items: Sequence[str]) -> str:
    items = [i for i in items if i]
    if not items:
        return ""
    if len(items) == 1:
        return items[0]
    return ", ".join(items[:-1]) + " and " + items[-1]


# ------------------------------------------------------------------- build
def build_report(record: dict, summary: dict, *,
                 input_qc: Optional[dict] = None,
                 record_bytes: Optional[bytes] = None) -> Report:
    """Assemble the report from a record and the summary derived from it.

    ``input_qc`` is the input check the command line writes beside the record.
    ``record_bytes`` is the serialised record the report describes; when given
    its SHA-256 is printed in the identification band so that a report and a
    record can be matched without trusting either.
    """
    md = record.get("metadata_diagnostics") or {}
    cfg = record.get("config") or {}

    digest = hashlib.sha256(
        json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:12]
    run_id = hashlib.sha256(
        (digest + str(record.get("seed")) + str(record.get("n_isolates")))
        .encode()).hexdigest()[:8].upper()
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ")
    record_sha = (hashlib.sha256(record_bytes).hexdigest()
                  if record_bytes is not None else None)

    # the per-agent share rows, once, for every section that reads them
    rows: List[dict] = []
    point_only: List[tuple] = []
    unread: List[tuple] = []
    for name, v in (md.get("clonal_share") or {}).items():
        if finite(v.get("kappa_adj")) and finite(v.get("observed_low")):
            rows.append({"name": name, "pt": v["kappa_adj"], "lo": v["observed_low"],
                         "hi": v["observed_high"], "null": v.get("null_mean"),
                         "p": v.get("p_value"), "p_floor": v.get("p_floor"),
                         "support": v.get("support"),
                         "estimable": v.get("estimable"),
                         "bound_lo": v.get("latent_order_lower"),
                         "bound_hi": v.get("latent_order_upper_bound"),
                         "bound_exact": v.get("latent_order_upper_exact", True),
                         "bound_limit": v.get("latent_order_lower_limit"),
                         "repeated": v.get("n_groups_repeated"),
                         "set_aside": v.get("n_singletons_set_aside"),
                         "effective": v.get("effective_groups")})
        elif finite(v.get("kappa_adj")):
            point_only.append((name, v))
        else:
            # a trait the estimator could not score at all: one lineage, a
            # constant call, or too few isolates; it is listed with the
            # reason rather than dropped from the report
            n_groups = v.get("n_groups")
            if n_groups == 0:
                why = "no readable outcome with a recorded lineage remains"
            elif n_groups == 1:
                why = "the isolates with a result sit in a single lineage, so there is no contrast between lineages"
            elif finite(v.get("prevalence")) and v["prevalence"] in (0.0, 1.0):
                why = "every isolate carries the same call, so there is no variance to attribute"
            elif isinstance(v.get("n"), int) and v["n"] < 2:
                why = "fewer than two isolates carry a result"
            else:
                why = "the share could not be scored on this collection"
            unread.append((name, why))
    rows.sort(key=lambda r: -r["pt"])
    n_cross_all = sum(1 for r in rows if r["lo"] <= 0.0 <= r["hi"])

    evidence = summary.get("evidence") or {}
    rejected = set(evidence.get("rejected_features") or [])
    # The panel decision of one analysis: the Benjamini-Yekutieli step-up on
    # the permutation p-values; the e-BH selection stays with the e-values.
    selection = md.get("lineage_selection") or {}
    lineage_selected = set(selection.get("rejected_features") or [])
    per_e = (md.get("lineage_evidence") or {}).get("per_feature") or {}
    seq_e = (md.get("lineage_evidence") or {}).get("sequential") or {}

    shares = [v for v in (md.get("clonal_share") or {}).values()
              if finite(v.get("support"))]
    groups_used = min((v.get("n_groups") for v in shares
                       if isinstance(v.get("n_groups"), int)), default=None)

    cz = md.get("censored_share") or {}
    cz_agents = cz.get("per_agent") or {}
    cz_join = cz.get("join") or {}

    call_alpha = selection.get("alpha") or 0.05
    call_q = selection.get("q_values") or {}

    n_refused = (sum(1 for r in rows if r.get("estimable") is False)
                 + sum(1 for _, v in point_only if v.get("estimable") is False))
    ident = {
        "run": run_id, "issued": stamp,
        "software": "amr-clonalshare %s" % __version__,
        "record": "clonal_share_result.json",
        "isolates": text(record.get("n_isolates")),
        "lineage": str(md["lineage_column"]) if md.get("lineage_column")
        else "none",
        "seed": text(record.get("seed")),
        "configuration": "sha256:%s" % digest,
        "record_sha256": record_sha or "",
    }
    n_scored = len(rows) + len(point_only)
    n_all = n_scored + len(unread)
    status = {
        "n_traits": n_all, "n_refused": n_refused + len(unread),
        "lead": ("%d of %d antimicrobials estimable" % (n_scored - n_refused,
                                                         n_all))
        if n_scored else ("no antimicrobial could be scored" if unread
                      else "no antimicrobial read"),
        "gates_line": (("The estimator withheld %d of %d; the condition that "
                        "failed is named per antimicrobial." % (n_refused, n_all)
                        if n_refused else "")
                       + ((" " if n_refused else "")
                          + "%d of %d could not be scored, for a reason named "
                          "below." % (len(unread), n_all) if unread else ""))
        if (n_refused or unread) else ("Every antimicrobial read was scored."
                                       if n_scored else
                                       "No antimicrobial carried a readable call."),
        "flag": bool(n_refused or unread),
    }
    mic_only = bool(cz_agents) and not n_all
    if mic_only:
        orders = [e.get("order") or {} for e in cz_agents.values()]
        mic_groups = [int(o["n_groups"]) for o in orders if isinstance(o.get("n_groups"), int)]
        groups_used = min(mic_groups) if mic_groups else None
    canonical = []
    if record.get('schema_version') == '1.0':
        from .outputs import result_rows, status_guidance
        canonical = result_rows(record)
        agents = {r['agent'] for r in canonical}
        n_computed = sum(r['status'] == 'computed' for r in canonical)
        n_full = sum(r['status'] == 'full_range' for r in canonical)
        n_missing = len(canonical) - n_computed - n_full
        status.update(n_traits=len(agents), n_refused=n_missing,
                      lead=f'{len(agents)} antimicrobials; {len(canonical)} analysis outcomes',
                      gates_line=f'{n_computed} computed, {n_full} full ranges, '
                                 + f'{n_missing} unavailable. Each result retains its own target and data subset.',
                      flag=bool(n_missing or n_full))
    elif mic_only:
        status.update(n_traits=len(cz_agents), lead=f'{len(cz_agents)} antimicrobials with MIC readings',
                      gates_line='MIC analysis statuses are reported per antimicrobial.')
    display_groups: Any = groups_used if groups_used is not None else "an unrecorded number of"

    S: List[Section] = []
    fig_no = [0]

    def fig_label() -> str:
        fig_no[0] += 1
        return "Figure %d" % fig_no[0]

    # ------------------------------------------------ 0 measurement summary
    s = Section("summary", "Measurement summary")
    if mic_only:
        s.p('**Quantity measured.** The share of the ordering of the recorded MIC readings within '
            'each stratum that lineage membership accounts for, and the lower bound on the lineage '
            'share of the latent MIC ordering established by the readings.')
    else:
        s.p("**Quantity measured.** The share of the variation in the recorded binary outcome "
            "across this collection that lineage membership accounts for, read "
            "from a lineage label and an interpreted result and scored on "
            "isolates the estimator did not see.")
    s.kv([
        ("Collection", "%s isolates in %s lineages%s"
         % (record.get("n_isolates"), display_groups,
            (", typed by `%s`" % md["lineage_column"])
            if md.get("lineage_column") else "; no lineage column was given")),
        ("Antimicrobials read", str(status['n_traits'])),
    ] + ([("Membership shares estimable", "%d of %d" % (n_scored - n_refused, n_all))]
         if n_scored else []))
    call_strata = md.get("call_strata") or {}
    if call_strata:
        s.p("**Strata of the calls.** The analyses of every call were read within the %s levels "
            "of `%s`: lineage labels are exchanged only within a level, the interval draws each "
            "isolate's level with its call, and the bounds a call places on the latent "
            "ordering are the bounds within levels, so a shift between levels is not read as a "
            "lineage difference. The e-values take no strata."
            % (call_strata.get("n_strata"), call_strata.get("column")))
    units = md.get("sampling_units") or {}
    if units:
        s.p("**Sampling units.** Isolates were grouped by `%s` into %s units (%s isolates "
            "without a recorded unit counted as units of their own). The permutation tests "
            "exchange lineage labels within units; the confidence intervals of the shares and "
            "the lower confidence limits treat the isolates of a lineage as independent, as "
            "stated with every result, and the decomposition resamples whole units."
            % (units.get("column"), units.get("n_units"), units.get("n_isolates_without_unit", 0)))
    s.p("**Reporting conditions.** %s An unavailable value failed a declared reporting condition "
        "of its estimator; the recorded reason distinguishes input limitations "
        "from numerical failure."
        % status["gates_line"])
    if canonical:
        s.table('Table 0', 'Available analyses and their data subsets. A full range is a completed, uninformative result; an unavailable value is not zero.',
                [('Agent', 'lab'), ('Analysis', 'lab'), ('N', 'num'), ('Status', 'lab'), ('Data scope', 'lab')],
                [(r['agent'], r['analysis'], r['n'], r['status'], r['analysis_scope']) for r in canonical])
        attention = [r for r in canonical if r['status'] != 'computed' or r.get('reason')]
        if attention:
            s.table('Table 0b', 'Reasons and next checks for results requiring interpretation. These checks do not change the recorded status or justify changing reporting conditions.',
                    [('Agent', 'lab'), ('Analysis / status', 'lab'), ('Reason', 'lab'), ('Next check', 'lab')],
                    [(r['agent'], r['analysis'] + ' / ' + r['status'], *status_guidance(r)) for r in attention])
    bounds = record.get('collection_bounds') or {}
    if bounds:
        s.table('Table 0a', 'Prevalence in the recorded collection: bounds allow every missing outcome to be either negative or positive. These are identification bounds, not confidence intervals. Observed prevalence uses observed outcomes only.',
                [('Agent', 'lab'), ('Recorded', 'num'), ('Observed', 'num'), ('Missing', 'num'), ('Observed prevalence', 'num'), ('Collection bounds', 'lab')],
                [(agent, b['total_count'], b['observed_count'], b['missing_count'],
                  pct(b['observed_prevalence']), interval(b['lower_bound'], b['upper_bound'])) for agent, b in bounds.items()])
    if (md.get('prevalence_decomposition') or {}).get('per_feature'):
        s.p('**Comparison scope.** Decomposition components describe the observed, labelled subsets. '
            'Equal typing coverage or a nonsignificant missingness test does not establish representativeness. '
            'Full observed prevalence and subset prevalence are recorded separately; missing data prevent automatic collection-wide generalization.')
    if rows:
        s.table("Table 1", "Every trait, ordered by share, with the 95 %% interval, "
                "the permutation p-value and the Benjamini-Yekutieli q-value across "
                "traits; selected is the Benjamini-Yekutieli selection at a "
                "false-discovery level of %s, which holds whatever the dependence "
                "between traits. Control is the share on shuffled lineage labels. "
                "A p-value written with ≤ is the smallest the permutations can give. "
                "The last column gives the conclusion of the results table of the form: "
                "established, lineage structure established, when the interval lies "
                "above zero and the trait is selected, or the lower confidence limit "
                "is above zero; not established; whole range, when the interval spans "
                "it; or not estimable. %s"
                % (num(call_alpha, 2), VERDICT_NOTE),
                [("Trait", "lab"), ("Share", "num"), ("95 % interval", "num"),
                 ("p", "num"), ("q", "num"), ("Control", "num"), ("e-value", "num"),
                 ("Selected", "lab"), ("Conclusion", "lab")],
                [(r["name"], num(r["pt"]), interval(r["lo"], r["hi"]),
                  format_p(r["p"], r["p_floor"]),
                  num(call_q.get(r["name"])), num(r["null"], 2),
                  evalue((per_e.get(r["name"]) or {}).get("e_value")),
                  "yes" if r["name"] in lineage_selected else "no",
                  short_label(agent_verdict(md["clonal_share"][r["name"]], alpha=call_alpha,
                                            q=call_q.get(r["name"]))["label"]))
                 for r in rows])
        s.table("Table 1b", "Lineage structure behind each share: the lineages "
                "with at least two isolates on which the share is scored, the "
                "isolates of singleton lineages set aside, the effective number of "
                "scored lineages (inverse of the sum of squared lineage shares) and "
                "the share of isolates in lineages with at least two members "
                "(support). Few effective lineages mean the share rests on a few "
                "lineages.",
                [("Trait", "lab"), ("Lineages scored", "num"), ("Singletons set aside", "num"),
                 ("Effective", "num"), ("Support", "num")],
                [(r["name"], num(r.get("repeated"), 0), num(r.get("set_aside"), 0),
                  num(r.get("effective"), 1), pct(r.get("support")))
                 for r in rows])
    if point_only:
        s.table("Table 1a", "Point estimates without a computed interval. "
                "The interval was not computed or no finite bootstrap interval "
                "was available; no interval-based interpretation is made.",
                [("Trait", "lab"), ("Share", "num"), ("Interval", "lab"),
                 ("Conclusion", "lab")],
                [(name, num(v["kappa_adj"]), "interval not computed",
                  "refused" if v.get("estimable") is False else "point estimate only")
                 for name, v in point_only])
    not_done = _not_evaluated(record, cz_agents)
    not_done = ["`%s`: not scored; %s." % (name, why) for name, why in unread] + not_done
    if not_done:
        s.callout("Not evaluated on this run, and why", not_done)
    inter = _interpretation(rows, lineage_selected)
    if inter:
        s.callout("Interpretation (generated from Table 1 by fixed rules)", inter)
    S.append(s)

    # --------------------------------------------- 1 admissibility of input
    s = Section("admissibility", "Admissibility of the input")
    prov = []
    if md.get("lineage_column"):
        prov.append(("Lineage column", "`%s`" % md["lineage_column"]))
    phenotype = (input_qc or record.get('input_qc') or {}).get('phenotype_reading') or {}
    if phenotype:
        prov.extend([('Phenotype interpretation', str(phenotype.get('phenotype_kind') or 'not supplied')),
                     ('Positive outcome', str(phenotype.get('positive_definition') or 'not supplied')),
                     ('Applied coding', str(phenotype.get('positive_coding') or 'not supplied')),
                     ('Intermediate policy', phenotype.get('intermediate_policy') or 'not applicable'),
                     ('Interpretation source', phenotype.get('source') or 'not supplied'),
                     ('AST standard/version', ' / '.join(str(phenotype.get(k) or 'not supplied') for k in ('ast_standard', 'ast_version')))])
    prov.append(("Susceptibility calls", "%s antimicrobials on %s isolates"
                 % (record.get("n_traits", len(rows)), record.get("n_isolates"))))
    if cz_join:
        prov.append(("Recorded dilutions",
                     "%s rows read, %s joined; %s of %s isolates carry a "
                     "recorded dilution" % (
                         cz_join.get("rows_read"), cz_join.get("rows_joined"),
                         cz_join.get("strains_with_mic"),
                         cz_join.get("strains_aligned"))))
    att = cfg.get("attribution") or {}
    if att:
        prov.append(("Resampling", "%s folds, %s repeats, %s bootstrap draws, "
                     "%s permutations per antimicrobial"
                     % (att.get("folds"), att.get("repeats"),
                        att.get("n_boot"), att.get("n_perm"))))
    if not cfg:
        prov.append(("Settings", "no configuration is recorded in this record"))
    s.kv(prov)
    if input_qc:
        lin = input_qc.get("lineage") or {}
        traits = input_qc.get("traits") or {}
        weak = sum(1 for t in traits.values()
                   if not t.get("adequate") and not t.get("constant"))
        if traits:
            s.p("Antimicrobials with fewer than %s isolates of the rarer "
                "outcome: %d of %d. Each method uses its own reporting conditions; "
                "sparse outcomes may yield wide intervals or an unavailable result."
                % (input_qc.get("min_minor_count"), weak, len(traits)))
        if lin:
            s.p("Lineage groups in the input: %s, of which %s hold a single "
                "isolate. Support, the share of isolates in lineages of at least two, "
                "is %s; the share is scored on those isolates and the singletons "
                "are set aside. %s"
                % (lin.get("n_groups"), lin.get("n_singletons"), pct(lin.get("support")),
                   "At least two lineages hold two or more isolates."
                   if lin.get("estimable")
                   else "Fewer than two lineages hold two or more isolates, so "
                   "there is no between-lineage comparison."))
            sizes = lin.get("group_sizes") or {}
            if isinstance(sizes, dict) and sizes:
                ordered = sorted(((str(k), int(v)) for k, v in sizes.items()),
                                 key=lambda kv: (-kv[1], kv[0]))
                n_show = min(len(ordered), 40)
                s.figure(fig_label(), "lineages", {
                    "sizes": ordered[:n_show], "n_total": len(ordered),
                    "n_isolates": lin.get("n_typed", record.get("n_isolates")),
                    "support": lin.get("support"),
                    "min_group": input_qc.get("min_group_size", 2)},
                    "Isolates per lineage, largest first%s. A lineage of one "
                    "isolate, drawn in orange, cannot be predicted out of "
                    "sample and is set aside; support is the share of isolates "
                    "in the other lineages, %s here. The largest lineage holds %s "
                    "of the isolates, which sets how much one lineage can "
                    "weigh in the share."
                    % ((", the %d largest of %d" % (n_show, len(ordered)))
                       if n_show < len(ordered) else "",
                       pct(lin.get("support")),
                       pct(ordered[0][1] / float(lin.get("n_typed") or
                                                  record.get("n_isolates") or 1))))
    else:
        s.p("No input check was attached to this record.")
    # the lineages each share is scored on and the support, over the calls and
    # the MIC orderings read on this run
    scored_blocks = [v for v in (md.get("clonal_share") or {}).values() if isinstance(v, dict)]
    scored_blocks += [(e.get("order") or {}) for e in cz_agents.values()]
    repeated_used = min((b["n_groups_repeated"] for b in scored_blocks
                         if isinstance(b.get("n_groups_repeated"), int)), default=None)
    support_used = min((b["support"] for b in scored_blocks if finite(b.get("support"))), default=None)
    s.table("Table 2", "Conditions the estimator requires before any result "
            "is reported.",
            [("Condition", "lab"), ("Observed", "num"), ("Required", "num"),
             ("Verdict", "lab")],
            [("Lineages with at least two isolates, fewest over the antimicrobials read",
              str(repeated_used) if repeated_used is not None else "not recorded",
              "≥ 2",
              "accepted" if repeated_used is not None and repeated_used >= 2
              else "refused" if repeated_used is not None else "not evaluated"),
             ("Support, lowest over the antimicrobials read",
              pct(support_used) if support_used is not None else "not recorded",
              "reported",
              "singletons set aside" if support_used is not None else "not evaluated")])
    feasibility = (input_qc or {}).get("agent_feasibility") or {}
    if feasibility:
        from .qc import FEASIBILITY_HEADERS, FEASIBILITY_SCOPE, feasibility_rows
        s.p(FEASIBILITY_SCOPE)
        s.table("Table 2a", "Per-agent input feasibility.",
                [(name, "lab") for name in FEASIBILITY_HEADERS],
                feasibility_rows(feasibility))
    S.append(s)

    # ------------------------------------------------ 2 the share in detail
    s = Section("share", "The lineage share of the call, trait by trait")
    if rows:
        n_show = min(len(rows), 18)
        shown = rows[:n_show]
        n_cross = sum(1 for r in shown if r["lo"] <= 0.0 <= r["hi"])
        s.p("The interval is for the share of the lineages in this collection: "
            "the lineages and their sizes are held fixed and the isolates of each "
            "lineage are drawn again from the smoothed distribution of its calls. For "
            "**%d of the %d traits** shown the interval reaches zero, so no "
            "lineage effect is distinguishable from none for that trait%s."
            % (n_cross, len(shown),
               ("; over the whole panel this holds for **%d of %d**"
                % (n_cross_all, len(rows))) if len(rows) > n_show else ""))
        s.p("The control column is the same estimator run on shuffled lineage "
            "labels; it should sit near zero, and a share is read against it "
            "rather than against zero.")
        s.figure(fig_label(), "intervals", shown,
                 "Lineage share of the call by trait, point estimate with 95 %% "
                 "interval for the represented lineages, %d largest of %d. Traits "
                 "whose interval reaches zero are drawn in grey and marked ‡."
                 % (n_show, len(rows)))
        with_bounds = any(finite(r.get("bound_lo")) for r in rows)
        if with_bounds:
            s.p("The last column asks what the call says about the ordering "
                "it was cut from: the lineage share of the latent MIC, or "
                "of any continuous value the call thresholds, the rank "
                "intraclass correlation of that value. A call does not identify "
                "that share; it bounds it. The bracket gives the smallest share "
                "any ordering consistent with the calls allows, the lower "
                "bound the calls establish, and the largest; an upper end marked ≤ is a "
                "certified bound rather than a share attained. The figure after "
                "it is the one-sided 95 %% lower confidence limit of the "
                "lower bound, from splitting the isolates of every lineage "
                "at random, choosing a direction on one half and testing it on "
                "the other: a share of the latent ordering that the lineages of "
                "this collection are shown to account for, with no model for the latent values beyond "
                "their agreeing with the readings. A MIC reading on a panel that "
                "contains the cut-off refines the call%s, and refining the "
                "readings can only narrow the bounds they allow." % (" (Table 5)" if cz_agents else ""))
        cols = [("Trait", "lab"), ("Share", "num"), ("95 % interval, represented lineages", "num")]
        if with_bounds:
            cols.append(("Latent ordering: bounds; 95 % lower limit", "num"))
        s.table("Table 3", "The share of every trait with its interval%s."
                % (", and the bounds the call places on the share of the "
                   "latent ordering" if with_bounds else ""),
                cols,
                [(r["name"], num(r["pt"]), interval(r["lo"], r["hi"]))
                 + (((_bounds(r.get("bound_lo"), r.get("bound_hi"), r.get("bound_limit"),
                              r.get("bound_exact", True)),)
                    if with_bounds else ()))
                 for r in rows])
    elif point_only:
        s.p("The per-trait point estimates are listed in the measurement summary; "
            "an interval not computed cannot support an interval-based reading.")
    else:
        s.p("No per-trait share is reported on this run.")
    S.append(s)

    # ------------------------------------- 3 evidence that survives re-reading
    s = Section("evidence", "Evidence that survives re-reading")
    if per_e:
        s.p("A p-value is a statement about one look at the data. A "
            "surveillance panel is looked at again every year, and a p-value "
            "recomputed each time loses its error control. The e-value is evidence "
            "on a scale made for that: this run's e-value is a statement about "
            "this collection, and a programme that adds an intake each year "
            "multiplies the e-value of each new intake, scored against the "
            "lineage rates learned from the earlier ones, into a running "
            "product whose error control holds at whatever intake it is read "
            "(sequential_e_process). The e-BH procedure controls the "
            "false-discovery rate across the panel whatever the dependence "
            "between traits. Larger is stronger; 1 is no evidence.")
        s.table("Table 4", "e-value per trait and the e-BH selection at level "
                "%s. %s" % (num(evidence.get("alpha"), 2),
                            ("The selection threshold on this run is %s; a "
                             "trait at or above it is selected."
                             % num(evidence.get("threshold"), 1))
                            if finite(evidence.get("threshold")) else
                            "No trait reached the e-BH threshold on this run."),
                [("Trait", "lab"), ("e-value", "num"), ("natural log", "num"),
                 ("Selected", "lab")],
                [(name, evalue(v.get("e_value")), num(v.get("log_e"), 2),
                  "yes" if name in rejected else "no")
                 for name, v in sorted(per_e.items(),
                                       key=lambda kv: -(kv[1].get("e_value")
                                                        or 0.0))])
        s.p("%s of %s traits are selected. %s"
            % (num(evidence.get("n_rejected"), 0), num(evidence.get("n_agents"), 0),
               ("Note: %s." % evidence["note"]) if evidence.get("note") else ""))
        e_rows = [{"name": name, "log_e": v.get("log_e"),
                   "selected": name in rejected}
                  for name, v in per_e.items() if finite(v.get("log_e"))]
        e_rows.sort(key=lambda r: -r["log_e"])
        if e_rows:
            s.figure(fig_label(), "evalues", {
                "rows": e_rows[:18], "alpha": evidence.get("alpha"),
                "threshold": evidence.get("threshold")},
                "Evidence per trait on the natural-log scale, %d largest of "
                "%d. The dashed rule is 1/α, the evidence one trait alone "
                "needs at level %s%s. Traits the e-BH procedure selected are "
                "drawn in blue; the scale is logarithmic, so equal steps are "
                "equal factors of evidence."
                % (min(len(e_rows), 18), len(e_rows),
                   num(evidence.get("alpha"), 2),
                   ("; the solid rule is the e-BH selection threshold on this "
                    "run, %s, which rises with the number of traits read "
                    "together, and a trait at or beyond it is selected"
                    % num(evidence.get("threshold"), 1))
                   if finite(evidence.get("threshold")) else
                   "; no trait reached the e-BH threshold on this run"))
        if seq_e:
            sizes = sorted(seq_e.get("n_per_batch") or [])
            s.p("This collection is also divided into intakes, %s of them read in the "
                "order of %s, from %s to %s isolates each%s. Each intake is "
                "scored against the lineage rates learned from the intakes "
                "before it, and the running product is the sequential "
                "e-value: it may be read after any intake without spending "
                "its error control. %s of %s traits are selected on the product, "
                "against %s at this single look, which is the price of "
                "error control that holds at whatever intake the programme is "
                "read at."
                % (len(seq_e.get("batches") or []),
                   seq_e.get("batch_column"),
                   sizes[0] if sizes else 0, sizes[-1] if sizes else 0,
                   (", with %s isolates carrying no intake and set aside"
                    % seq_e["n_without_batch"])
                   if seq_e.get("n_without_batch") else "",
                   (seq_e.get("e_bh") or {}).get("n_rejected"),
                   (seq_e.get("e_bh") or {}).get("m"),
                   evidence.get("n_rejected")))
            s.table("Table 4b", "The sequential e-value per trait after the "
                    "last intake, and the e-BH selection taken on the log "
                    "scale, where a product over intakes does not overflow.",
                    [("Trait", "lab"), ("natural log of the product", "num"),
                     ("Selected", "lab")],
                    [(name, num(v.get("log_e")[-1] if v.get("log_e") else None,
                                2),
                      "yes" if name in set(
                          (seq_e.get("e_bh") or {}).get(
                              "rejected_features") or []) else "no")
                     for name, v in sorted(
                         (seq_e.get("per_feature") or {}).items(),
                         key=lambda kv: -kv[1]["log_e"][-1]
                         if kv[1].get("log_e") and finite(kv[1]["log_e"][-1])
                         else float("inf"))])
            seq_sel = set((seq_e.get("e_bh") or {}).get("rejected_features")
                          or [])
            traj = [{"name": name, "log_e": [x for x in v.get("log_e") or []],
                     "selected": name in seq_sel}
                    for name, v in (seq_e.get("per_feature") or {}).items()
                    if v.get("log_e")]
            if traj and seq_e.get("batches"):
                s.figure(fig_label(), "sequential", {
                    "batches": [str(b) for b in seq_e.get("batches") or []],
                    "n_per_batch": list(seq_e.get("n_per_batch") or []),
                    "rows": traj, "alpha": evidence.get("alpha"),
                    "threshold": (seq_e.get("e_bh") or {}).get("threshold")},
                    "The running product of e-values over the %d intakes, one "
                    "line per trait, on the natural-log scale. Each intake is "
                    "scored against the lineage rates learned from the intakes "
                    "before it, so the first intakes carry no evidence and a "
                    "line may fall as well as rise. The dashed rule is 1/α; "
                    "%s"
                    "Traits selected on the product are drawn in blue and "
                    "named at the right. Stopping-time panel control requires "
                    "validity in the joint panel filtration; repeated rejection "
                    "unions are not controlled. The frame is cut at −25: a line "
                    "below it gives little or no evidence against lineage independence. "
                    "It does not establish absence of a lineage effect."
                    % (len(seq_e.get("batches") or []),
                       ("the solid rule is the e-BH threshold on the product, "
                        "%s. " % num((seq_e.get("e_bh") or {}).get("threshold"), 1))
                       if finite((seq_e.get("e_bh") or {}).get("threshold"))
                       else "no trait reached the e-BH threshold on the product. "))
    else:
        s.p("No e-values were computed on this run.")
    S.append(s)

    # ------------------------------------- 4 reading at recorded resolution
    s = Section("censored", "Reading at the recorded resolution")
    if cz_agents:
        order_selection = cz.get("order_selection") or {}
        order_selected = set(order_selection.get("rejected_features") or [])
        order_q = order_selection.get("q_values") or {}
        by = sorted({str((e.get("order") or {}).get("stratified_by")) for e in cz_agents.values()
                     if (e.get("order") or {}).get("stratified_by")})
        s.p("Where a dilution was recorded, the lineage question is also read "
            "from the dilution itself rather than from the call derived from "
            "it. Every reading is placed by its position among the readings "
            "of %s: the share of those readings below it plus half the share "
            "at the same reading, with a censored reading placed by the "
            "nonparametric maximum-likelihood distribution of the readings. "
            "The share of these positions that lineage accounts for is estimated as "
            "the share of a call is, and its control shuffles the lineage "
            "labels %s. No distribution is assumed for the readings, a change "
            "of the MIC scale such as mg/L to log2 leaves the share unchanged, "
            "and %s. The share of the call and this share are read on "
            "independently retained readable and typed subsets; their "
            "denominators may differ."
            % (("its own stratum, a level of `%s`" % _join(by)) if by else "the collection",
               "only within strata" if by else "across the collection",
               "a constant offset or a different panel between strata drops out"
               if by else "the readings of one panel are compared with one another"))
        crow = []
        within_notes = []
        for name, e in sorted(cz_agents.items()):
            o = e.get("order") or {}
            if o.get("estimable") and finite(o.get("kappa_adj")):
                crow.append((name, num(o.get("n_readings"), 0), num(o.get("kappa_adj")),
                             interval(o.get("observed_low"), o.get("observed_high")),
                             format_p(o.get("p_value"), o.get("p_floor")),
                             num(order_q.get(name)), num(o.get("null_mean"), 2),
                             "yes" if name in order_selected else "no",
                             _bounds(o.get("latent_order_lower"), o.get("latent_order_upper_bound"),
                                     o.get("latent_order_lower_limit"),
                                     o.get("latent_order_upper_exact", True)),
                             num(o.get("panel_resolution")),
                             pct(o.get("share_end_wells")), str(len(o.get("strata") or {}) or ""),
                             short_label(agent_verdict(o, alpha=order_selection.get("alpha") or 0.05,
                                                       q=order_q.get(name))["label"])))
                within = o.get("latent_order_lower_within_strata") or {}
                strong = {k: v for k, v in within.items() if finite(v) and finite(o.get("latent_order_lower"))
                          and v > o["latent_order_lower"] + 0.1}
                if strong:
                    within_notes.append("%s (pooled %s; %s)" % (
                        name, num(o.get("latent_order_lower")),
                        ", ".join("%s %s" % (k, num(v)) for k, v in sorted(strong.items()))))
            else:
                # the agent did not reach the estimator, or the share could
                # not be scored; the reason takes the place of the estimate
                if o.get("reason"):
                    why = o["reason"]
                elif (o.get("n_groups_repeated") or 0) < 2:
                    why = "fewer than two lineages hold repeated readings"
                else:
                    why = "the share could not be scored on these readings"
                crow.append((name, num(o.get("n_readings", e.get("n")), 0), "not estimable", why, "", "", "", "",
                             "", "", pct(o.get("share_end_wells")), str(len(o.get("strata") or {}) or ""),
                             "not estimable"))
        sharp = all(bool((e.get("order") or {}).get("latent_order_sharp", True)) for e in cz_agents.values())
        s.p("The readings also bound the lineage share of the MIC ordering "
            "they were read from, the ordering of the latent MICs "
            "within %s: the rank intraclass correlation of the MIC itself. "
            "A reading says only that the MIC lies in its range, so the "
            "readings do not identify that share; they bound it. The lower "
            "bound is the smallest share that any arrangement of the latent "
            "MICs within their readings allows, the lower bound the readings "
            "establish, computed with a certified error; the upper end is the "
            "largest share, attained by an arrangement, or, marked ≤, a certified "
            "bound on it. Every share between them is allowed. %sThe figure "
            "after the bounds is the one-sided 95 %% lower confidence limit of "
            "the lower bound, and so of the share: the isolates of every "
            "lineage are split at random, one half chooses the lineage scores "
            "on which the lower bound rests and the other half tests them, and "
            "then the other way round; the tests of 25 such splits are averaged, "
            "and their average spread sets the limit. A "
            "finer panel never lowers the lower bound, and a call, one cut "
            "of the panel, never raises it. Resolution is the share of the "
            "ordering the panel resolves, one less the expected squared width "
            "of a reading on the scale of positions; it is the factor by which "
            "the share of the midpoint arrangement falls short of the plug-in share "
            "of the scores."
            % (("a level of `%s`" % _join(by)) if by else "the collection",
               "" if sharp else "Where readings of one stratum overlap, each run of "
               "overlapping readings is treated as one cell within which the "
               "isolates may stand in any order, so the bounds hold for every "
               "arrangement the readings allow but may not be attained. "))
        if within_notes:
            s.p("The share of the latent ordering is pooled over strata: a "
                "lineage placed high in one stratum and low in another does not "
                "count as ordered. For %s the readings of some stratum on its own "
                "establish a larger lower bound than the pooled one, so lineage effects "
                "differ between strata: %s." % ("one agent" if len(within_notes) == 1
                                                 else "%d agents" % len(within_notes),
                                                 "; ".join(within_notes)))
        s.table("Table 5", "Lineage share of the MIC ordering per agent, with the "
                "95 %% interval for the lineages of this collection and the "
                "permutation p-value and the Benjamini-Yekutieli q-value across "
                "agents; selected is the Benjamini-Yekutieli selection at level %s. "
                "Control is the share on shuffled lineage labels. "
                "The bounds are those the "
                "readings place on the share of the latent MIC ordering, "
                "with the 95 %% lower confidence limit of the lower bound. "
                "Readings on an end well are tied at the panel's edge and carry "
                "less of the ordering; their share is given beside the estimate. "
                "A p-value written with ≤ is the smallest the permutations can give. "
                "The last column gives the conclusion: established, lineage structure "
                "established, when the interval lies above zero and the agent is "
                "selected, or the lower confidence limit is above zero; not "
                "established; whole range, when the interval spans it; or not "
                "estimable. %s"
                % (num(order_selection.get("alpha"), 2), VERDICT_NOTE),
                [("Agent", "lab"), ("Readings", "num"), ("Share", "num"),
                 ("95 % interval, represented lineages", "num"), ("p", "num"), ("q", "num"),
                 ("Control", "num"), ("Selected", "lab"),
                 ("Latent ordering: bounds; 95 % lower limit", "num"), ("Resolution", "num"),
                 ("End wells", "num"), ("Strata", "num"), ("Conclusion", "lab")], crow)
        if order_selection.get("resolution_note"):
            s.p("Resolution of the selection: %s." % order_selection["resolution_note"])
        s.table("Table 5b", "Lineage structure behind each share of the MIC "
                "ordering: the lineages with at least two readings on which the "
                "share is scored, the readings of singleton lineages set aside, the "
                "effective number of scored lineages (inverse of the sum of squared "
                "lineage shares) and the share of readings in lineages with at least "
                "two members (support). Few effective lineages mean the share rests "
                "on a few lineages.",
                [("Agent", "lab"), ("Lineages scored", "num"), ("Singletons set aside", "num"),
                 ("Effective", "num"), ("Support", "num")],
                [(name, num((e.get("order") or {}).get("n_groups_repeated"), 0),
                  num((e.get("order") or {}).get("n_singletons_set_aside"), 0),
                  num((e.get("order") or {}).get("effective_groups"), 1),
                  pct((e.get("order") or {}).get("support")))
                 for name, e in sorted(cz_agents.items())])
        call_pt = {r["name"]: r["pt"] for r in rows}
        drows = []
        for name, e in sorted(cz_agents.items()):
            o = e.get("order") or {}
            if o.get("estimable") and finite(o.get("kappa_adj")):
                drows.append({"name": name, "pt": o["kappa_adj"],
                              "lo": o.get("observed_low"), "hi": o.get("observed_high"),
                              "censored": o.get("share_end_wells"),
                              "selected": name in order_selected,
                              "call": call_pt.get(name)})
        drows.sort(key=lambda r: -r["pt"])
        if drows:
            s.figure(fig_label(), "dilution", drows[:18],
                     "Lineage share of the MIC ordering per agent, %d largest of "
                     "%d: the point estimate with its 95 %% interval; agents "
                     "selected across the panel are drawn in blue. The hollow "
                     "diamond is the share of the binary call from Table 1 for "
                     "the same agent, a different quantity whose retained "
                     "isolates may differ, placed here so that the two readings "
                     "can be seen side by side. The figure in parentheses at the "
                     "right is the share of readings on an end well."
                     % (min(len(drows), 18), len(drows)))
        _panel_detail(s, cz_agents, fig_label)
    else:
        s.p("No recorded dilutions were supplied on this run; the shares "
            "above are read from binary calls only.")
    S.append(s)

    # ------------------------------------------------------------ 4b strata
    strata = record.get("strata") or {}
    if strata.get("levels"):
        s = Section("strata", "The same analyses within each stratum")
        levels = strata["levels"]
        s.p("The run was repeated on the isolates of each level of `%s` "
            "(%s; %d isolates without a level were left out). Each stratum is "
            "an analysis of its own, on its own lineages and its own panel, and "
            "the shares below are not adjusted for one another. A share that "
            "holds within every stratum is a property of the lineages; a share "
            "seen only in the pooled run and in no stratum is accounted for by "
            "the differences between strata."
            % (strata.get("column"),
               _join(["%s (n = %d)" % (lv, rec.get("n_isolates", 0)) for lv, rec in levels.items()]),
               strata.get("n_without_level", 0)))
        st_rows = []
        for lv, rec in levels.items():
            mdd = rec.get("metadata_diagnostics") or {}
            calls = mdd.get("clonal_share") or {}
            mics = (mdd.get("censored_share") or {}).get("per_agent") or {}
            for name in sorted(set(calls) | set(mics)):
                c = calls.get(name) or {}
                o = (mics.get(name) or {}).get("order") or {}
                st_rows.append((str(lv), name, num(c.get("n"), 0) if finite(c.get("n")) else "",
                                num(c.get("kappa_adj")) if finite(c.get("kappa_adj"))
                                else ("not computed: %s" % c["reason"] if c.get("reason") else "not computed"),
                                interval(c.get("observed_low"), c.get("observed_high"))
                                if finite(c.get("observed_low")) else "",
                                format_p(c.get("p_value"), c.get("p_floor")) if finite(c.get("p_value")) else "",
                                ("%s (%s)" % (num(o.get("kappa_adj")),
                                              interval(o.get("observed_low"), o.get("observed_high"))))
                                if o.get("estimable") and finite(o.get("kappa_adj"))
                                else (("not estimable: %s" % o["reason"] if o.get("reason") else "not estimable")
                                      if o else ""),
                                format_p(o.get("p_value"), o.get("p_floor"))
                                if o.get("estimable") and finite(o.get("p_value")) else ""))
        s.table("Table 5d", "Lineage share per stratum and agent: the share of the "
                "binary call with its interval and p-value, and the lineage share "
                "of the MIC ordering with its interval and p-value where a MIC "
                "table was supplied; a share not computed carries its reason. "
                "strata_results.csv holds these and the other analyses of every "
                "stratum.",
                [("Stratum", "lab"), ("Agent", "lab"), ("n", "num"),
                 ("Call share", "num"), ("95 % interval", "num"), ("p", "num"),
                 ("MIC ordering share (95 % interval)", "num"), ("p", "num")], st_rows)
        S.append(s)

    # ------------------------------------------ 5 composition or rate
    s = Section("decomposition", "A change of lineages or a change within them")
    sv = summary.get("surveillance") or {}
    dec = sv.get("decomposition")
    if dec:
        n_ref = int(dec.get("n_within_lineage_refused") or 0)
        n_feat = int(dec.get("n_features") or 0)
        s.p("Contrast %s: of %s traits, %s have a composition component (a change "
            "in lineage mix) and %s a within-lineage component (a change in rate) "
            "selected by the Benjamini-Yekutieli step-up within each component "
            "family; the intervals and the step-up rest on nominal percentile-bootstrap "
            "tail probabilities (their calibration is reported with the package)."
            % (dec.get("contrast"), n_feat,
               dec.get("n_composition_discoveries"),
               dec.get("n_within_lineage_discoveries")))
        n_ref_comp = int(dec.get("n_composition_refused") or 0)
        if n_ref_comp:
            s.p("Both components were refused for %s: the lineage labels are "
                "missing unevenly between the two collections and in "
                "association with the trait, so the labelled isolates are not "
                "the collections and neither component is counted above."
                % ("every trait" if n_ref_comp == n_feat else
                   "%d of the %d traits" % (n_ref_comp, n_feat)))
        if n_ref > n_ref_comp:
            s.p("The within-lineage component was refused for %s: too few "
                "isolates sit in lineages the two collections share, so a "
                "change inside lineages is not identified there and is not "
                "counted above."
                % ("every trait" if n_ref - n_ref_comp == n_feat else
                   "%d of the %d traits" % (n_ref - n_ref_comp, n_feat)))
        if dec.get("warning"):
            s.p(dec["warning"][0].upper() + dec["warning"][1:] + ".")
        if dec.get("n_offsetting"):
            s.p("%s trait(s) show a mix change and a rate change of opposite "
                "sign that cancel: %s. A prevalence table shows these traits "
                "as unchanged although both components moved."
                % (dec["n_offsetting"],
                   _join(list(dec.get("offsetting_features") or [])[:12])))
    per_feature = (md.get("prevalence_decomposition") or {}).get("per_feature") or {}
    if per_feature:
        def selected(block):
            chosen = [name for name, key in
                      (("composition", "composition_discovery"),
                       ("within lineage", "within_lineage_discovery"))
                      if block.get(key)]
            return _join(chosen) if chosen else "neither"

        def limits(block, key):
            pair = block.get(key) or [None, None]
            return interval(pair[0], pair[1])

        levels = (md.get("prevalence_decomposition") or {}).get("levels") or ["first", "second"]
        s.table("Table 6", "The prevalence difference of each trait, %s minus %s, split "
                "into a change in lineage composition and a change in rate within "
                "lineages. The two components sum to the difference. The "
                "limits are bootstrap percentiles; p is the bootstrap p-value of "
                "each component and q its Benjamini-Yekutieli q-value within the "
                "component family; the last column names the components the "
                "false-discovery procedure selected." % (levels[0], levels[-1]),
                [("Trait", "lab"), ("Difference", "num"), ("Composition", "num"),
                 ("95 % interval", "num"), ("p", "num"), ("q", "num"), ("Within lineage", "num"),
                 ("95 % interval", "num"), ("p", "num"), ("q", "num"), ("Selected", "lab")],
                [(name, num(block.get("difference")),
                  num(block.get("composition")), limits(block, "composition_ci95"),
                  num(block.get("composition_p")), num(block.get("composition_q")),
                  num(block.get("within_lineage")),
                  limits(block, "within_lineage_ci95"),
                  num(block.get("within_lineage_p")), num(block.get("within_lineage_q")),
                  selected(block) if block.get("status") == "ok"
                  else str(block.get("status")))
                 for name, block in sorted(per_feature.items())])
        s.table("Table 6b", "The two collections behind each difference: the "
                "prevalence and the isolates of each, the lineages both hold and "
                "those seen in one only, the share of isolates in lineages both "
                "hold (shared support, the mean of the two collections; the "
                "within-lineage component needs at least %s) and the turnover "
                "share, the part of the difference carried by lineages seen in one "
                "collection only."
                % num(next((b.get("shared_support_threshold") for b in per_feature.values()
                            if finite(b.get("shared_support_threshold"))), None), 2),
                [("Trait", "lab"), ("Prevalence, %s" % levels[0], "num"), ("n", "num"),
                 ("Prevalence, %s" % levels[-1], "num"), ("n", "num"),
                 ("Lineages shared; only first; only second", "num"),
                 ("Shared support", "num"), ("Turnover share", "num")],
                [(name, pct(block.get("prevalence_a")), num(block.get("n_a"), 0),
                  pct(block.get("prevalence_b")), num(block.get("n_b"), 0),
                  "%s; %s; %s" % (num(block.get("n_lineages_shared"), 0), num(block.get("n_lineages_only_a"), 0),
                                  num(block.get("n_lineages_only_b"), 0)),
                  pct(block.get("shared_support_isolate_share")), num(block.get("turnover_share")))
                 for name, block in sorted(per_feature.items())])
    else:
        s.p("A prevalence difference between two collections can be split "
            "into a change in lineage composition and a change in rate within "
            "lineages. That decomposition needs two collections and was not "
            "run here: this record holds one.")
    S.append(s)

    # ---------------------------------------------------- provenance, terms
    s = Section("provenance", "Provenance and terms")
    s.kv([("Software", "amr-clonalshare %s, record schema %s"
           % (__version__, record.get("schema_version", ""))),
          ("Seed", str(record.get("seed"))),
          ("Configuration", "`sha256:%s`" % digest),
          ("Record", ("`sha256:%s`" % record_sha) if record_sha
           else "digest not available: the report was built from an object, "
                "not from the file")])
    interval_terms = [
        "**Interval for the represented lineages.** The 95 % interval of a share holds "
        "the lineages of this collection and their sizes fixed and draws the "
        "isolates of every lineage again from the smoothed distribution of its readings, "
        "scoring every draw as the data were scored and studentizing it by its "
        "own standard error (the bootstrap-t interval). It does not describe a "
        "fresh draw of lineages from the species."]
    model_terms = []
    if cz_agents:
        model_terms.append(
            "**MIC ordering share.** The share read from the dilution: how much "
            "of the ordering of the readings within a laboratory (or another "
            "declared stratum) is associated with lineage. It assumes no "
            "distribution for the readings.")
    if cz_agents or any(finite(r.get("bound_lo")) for r in rows):
        model_terms.append(
            "**Bounds on the latent ordering.** The smallest and the largest "
            "lineage share of the ordering of the latent values (the MIC "
            "behind a reading or a call) that the readings allow; the smallest "
            "is the lower bound the readings establish. They are identification "
            "bounds, not a confidence interval; the one-sided lower confidence "
            "limit beside them adds the sampling of isolates.")
    s.callout("Terms used in this report", [
        "**Share.** How much of the variation of the call between isolates "
        "is associated with their recorded lineage. Near 1: a strong "
        "association on the measured scale. Near 0: little association on "
        "that scale. This does not identify transmission or its mechanism.",
    ] + model_terms + interval_terms + [
        "**Control.** The same calculation on shuffled lineage labels, which "
        "should return roughly zero; a share is read against it.",
        "**Support.** The share of isolates in lineages with at least two "
        "members. A singleton lineage cannot be predicted out of sample, so the "
        "share is scored on the other isolates and the singletons are set "
        "aside; the share then speaks for singleton lineages only if they are "
        "drawn from the same population of lineages as the repeated ones.",
        "**Reporting condition.** A condition fixed before the run. If it does not hold, the "
        "software withholds the estimate and names the condition that failed.",
        "**e-value.** A fixed-look e-value and the running evidence over intakes have "
        "different uses. Repeated-look validity requires the sequential "
        "construction and its stated null-model assumptions; a fixed-look "
        "value alone does not justify repeated inspection. Larger values "
        "are stronger evidence against the specified null.",
    ])
    S.append(s)

    symbols = [
        ("‡", "The 95 % interval includes zero. No lineage effect is "
              "distinguishable from none for this trait, and the point "
              "estimate must not be read on its own."),
        ("§", "Derived arithmetically from a quantity recorded elsewhere in "
              "the record rather than estimated in this run."),
    ]
    colophon = [
        "amr-clonalshare %s · record schema %s · seed %s · configuration "
        "sha256:%s" % (__version__, record.get("schema_version", ""),
                       record.get("seed"), digest),
        "Cite this run as: “amr-clonalshare %s, run %s, record sha256:%s.”"
        % (__version__, run_id, (record_sha or "")[:12] or "not available"),
        "This report supersedes any earlier report bearing the same run "
        "identifier. It is regenerated from the record and holds no value "
        "that the record does not. The symbols carry the same wording in every "
        "run of this software; no symbol against a value means only that none "
        "of the listed conditions fired.",
    ]
    return Report(ident=ident, status=status, sections=S, symbols=symbols,
                  colophon=colophon)


#: lineages shown per agent in the heat map, the largest first
HEATMAP_LINEAGES = 14


def _panel_detail(s: "Section", cz_agents: dict, fig_label) -> None:
    """The panel each laboratory tested, followed by the lineage by dilution
    counts as a heat map."""
    geo_rows = []
    sources = set()
    for name, e in sorted(cz_agents.items()):
        g = e.get("panel") or {}
        per = g.get("per_panel") or ({"all readings": g} if g else {})
        sources.add(str(e.get("panel_source") or ""))
        for lab, geo in per.items():
            geo_rows.append((name, str(lab), str(geo.get("n_wells", "")),
                             "%s to %s" % (num(geo.get("lowest"), 4), num(geo.get("highest"), 4)),
                             "yes" if geo.get("doubling") else "no",
                             pct(geo.get("share_on_lowest")), pct(geo.get("share_on_highest"))))
    if geo_rows:
        s.p("Each laboratory is read on the wells that laboratory tested: an "
            "end-well reading is censored at the panel edge of its own "
            "laboratory, not at the widest edge in the collection. The table "
            "below gives, per agent and laboratory, the panel the readings were "
            "taken to come from.")
        presets = sorted(p for p in sources if p.startswith("preset "))
        if presets:
            s.p("The tested concentrations of %s were taken from a shipped panel preset (%s): "
                "a preset transcribes the published range of the panel, so the plate sheet "
                "of the laboratory should be checked against it." % (
                    "every agent" if len(sources) == 1 else "some agents",
                    "; ".join(p[len("preset "):].split(";")[0] for p in presets)))
        s.table("Table 5a", "Panel geometry per agent and testing laboratory: the "
                "wells taken as tested, their range, whether they form a doubling "
                "series, and the shares of readings on the lowest and the highest well.",
                [("Agent", "lab"), ("Laboratory", "lab"), ("Wells", "num"), ("Range", "num"),
                 ("Doubling", "lab"), ("On lowest well", "num"), ("On highest well", "num")],
                geo_rows)
    heat = [(name, e["dilution_table"]) for name, e in sorted(cz_agents.items())
            if isinstance(e.get("dilution_table"), dict) and e["dilution_table"].get("counts")]
    if heat:
        most = max(len(d.get("lineages") or []) for _, d in heat)
        s.figure(fig_label(), "heatmap", heat[:6],
                 "Readings per lineage and dilution interval for %d of %d agents%s: "
                 "each cell counts the isolates of one lineage whose reading fell "
                 "in one interval of the panel, darker for more. Open intervals "
                 "at the panel edges are the censored readings. A lineage whose "
                 "readings sit in one or two adjacent cells contributes to the "
                 "lineage share; readings spread along a row do not."
                 % (min(len(heat), 6), len(heat),
                    (", the %d largest of %d lineages" % (HEATMAP_LINEAGES, most))
                    if most > HEATMAP_LINEAGES else ""))


# ------------------------------------------------------------- the rules
def _not_evaluated(record: dict, cz_agents: dict) -> List[str]:
    """Every instrument the run did not use, with the reason the record gives."""
    out: List[str] = []
    md = record.get("metadata_diagnostics") or {}
    if not (md.get("prevalence_decomposition") or {}).get("per_feature"):
        out.append("Decomposition of a prevalence difference into lineage "
                   "composition and within-lineage rate: needs two collections; "
                   "this record holds one.")
    if not cz_agents:
        out.append("Reading at the recorded dilution: no dilutions were "
                   "supplied; shares are read from binary calls.")
    if not (md.get("lineage_evidence") or {}).get("per_feature"):
        out.append("e-values: not computed on this run.")
    return out


def _interpretation(rows: List[dict], selected: set) -> List[str]:
    """Three sentences fixed by the readings in Table 1. Empty sets drop out."""
    if not rows:
        return []
    admitted = [r for r in rows if r.get("estimable") is not False]
    with_evidence = [r["name"] for r in admitted
                     if not (r["lo"] <= 0.0 <= r["hi"]) and r["name"] in selected
                     and r["pt"] >= 0.5]
    small_effect = [r["name"] for r in admitted
                    if not (r["lo"] <= 0.0 <= r["hi"]) and r["name"] in selected
                    and r["pt"] < 0.5]
    unresolved = [r["name"] for r in admitted
                  if (r["lo"] <= 0.0 <= r["hi"]) and r["name"] in selected]
    none_detectable = [r["name"] for r in admitted
                       if (r["lo"] <= 0.0 <= r["hi"]) and r["name"] not in selected]
    interval_only = [r["name"] for r in admitted
                     if not (r["lo"] <= 0.0 <= r["hi"]) and r["name"] not in selected]
    refused = [r["name"] for r in rows if r.get("estimable") is False]
    out: List[str] = []
    if with_evidence:
        out.append("For %s, lineage membership is associated with the recorded outcome: "
                   "the estimated share is at or above one half in this "
                   "collection at the chosen typing resolution."
                   % _join(with_evidence))
    if small_effect:
        out.append("For %s, a lineage association is detected and the estimated "
                   "share is below one half. Most variation is not explained "
                   "by the lineage labels at this typing resolution."
                   % _join(small_effect))
    if unresolved:
        out.append("For %s, the panel selection finds a lineage effect while the "
                   "interval for its size still reaches zero: there is "
                   "evidence that positive outcomes are not spread evenly across the "
                   "lineages, and this collection is too small, or too "
                   "uneven across its lineages, to say how much of it the "
                   "lineages account for. Read the selection as the finding and the "
                   "share as not yet resolved." % _join(unresolved))
    if interval_only:
        out.append("For %s, the interval for the share excludes zero while the "
                   "panel selection does not find a lineage effect at its "
                   "false-discovery level: the two readings disagree, and the "
                   "selection, which carries the multiplicity correction, is the "
                   "one to read; treat the share as suggestive."
                   % _join(interval_only))
    if none_detectable:
        out.append("For %s, no lineage effect is distinguishable from none in "
                   "this collection at the chosen resolution. This does not "
                   "establish independence from lineage or identify the "
                   "reason for a prevalence change." % _join(none_detectable))
    if refused:
        out.append("For %s, the estimator refused because a declared reporting "
                   "condition failed. The result record names that condition; "
                   "no admitted estimate follows." % _join(refused))
    if out:
        out.append("These readings describe association in the sampled "
                   "collection. The analysis does not identify transmission, "
                   "horizontal transfer, selection, or the effect of an "
                   "intervention.")
    return out
