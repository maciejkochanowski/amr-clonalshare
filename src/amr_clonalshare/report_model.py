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
from importlib import resources
from typing import Any, Dict, List, Optional, Sequence, Tuple

from . import __version__

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
    return ("%%.%df" % nd) % float(x)


def text(x) -> str:
    """A recorded value as text; a value the record does not hold is said so."""
    return "not recorded" if x is None else str(x)


def pct(x, nd: int = 1) -> str:
    return "not computed" if not finite(x) else ("%%.%df %%%%" % nd) % (100.0 * float(x))


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


def _shape_reading(check) -> str:
    """The model check of a calibrated interval, as a table cell."""
    if not check:
        return "not run"
    return "%s (p = %s)" % ("rejected" if check.get("rejected") else "not rejected",
                             num(check.get("p_value")))


def interval(lo, hi, nd: int = 3) -> str:
    if not (finite(lo) and finite(hi)):
        return "not computed"
    return "%s to %s" % (num(lo, nd), num(hi, nd))


def _join(items: Sequence[str]) -> str:
    items = [i for i in items if i]
    if not items:
        return ""
    if len(items) == 1:
        return items[0]
    return ", ".join(items[:-1]) + " and " + items[-1]


def _load_validation_grid() -> Optional[dict]:
    try:
        text = (resources.files(__package__) / "validation_grid.json").read_text(
            encoding="utf-8")
        return json.loads(text)
    except Exception:
        return None


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
        if finite(v.get("kappa_adj")) and finite(v.get("ci_low")):
            rows.append({"name": name, "pt": v["kappa_adj"], "lo": v["ci_low"],
                         "hi": v["ci_high"], "null": v.get("null_mean"),
                         "sp_lo": v.get("superpopulation_low"),
                         "sp_pt": v.get("superpopulation_share"),
                         "sp_hi": v.get("superpopulation_high"),
                         "p": v.get("p_value"), "support": v.get("support"),
                         "estimable": v.get("estimable"),
                         "latent": v.get("latent_share"),
                         "latent_lo": v.get("latent_low"),
                         "latent_hi": v.get("latent_high"),
                         "dominant": v.get("dominant_lineage_share"),
                         "few": bool(v.get("few_lineages")),
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
    support = min((v["support"] for v in shares), default=None)
    groups_used = min((v.get("n_groups") for v in shares
                       if isinstance(v.get("n_groups"), int)), default=None)

    cz = md.get("censored_share") or {}
    cz_agents = cz.get("per_agent") or {}
    cz_join = cz.get("join") or {}

    def reading(r) -> str:
        crosses = r["lo"] <= 0.0 <= r["hi"]
        if r.get("estimable") is False:
            return "refused"
        if crosses and r["name"] in lineage_selected:
            return "lineage effect selected; its size is not resolved"
        if crosses:
            return "no detectable lineage effect"
        if r["name"] in lineage_selected:
            return "evidence of a lineage effect"
        return "interval excludes zero; not selected at the panel's level"

    def marks(r) -> str:
        m = ""
        if r["lo"] <= 0.0 <= r["hi"]:
            m += "‡"
        return m

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
    profiles = md.get("population_probit_icc") or {}
    general_profiles = {k: r for k, r in profiles.items() if r.get("method") == "general_population_probit"}
    fixed_cutoff_profiles = {k: r for k, r in profiles.items() if k not in general_profiles}
    profile_only = bool(profiles) and (cfg.get("attribution") or {}).get("enabled") is False
    if profile_only:
        computed = sum(r.get("complete", r.get("status") == "ok") for r in profiles.values())
        status.update(n_traits=len(profiles), n_refused=len(profiles)-computed,
                      lead=f"{computed} of {len(profiles)} population ICC profiles computed",
                      gates_line="Classical attribution was disabled. Population-model statuses are reported separately.",
                      flag=bool(computed < len(profiles)))
    mic_only = bool(cz_agents) and not n_all and not profiles
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
                      gates_line=f'{n_computed} computed, {n_full} full ranges, {n_missing} unavailable or incomplete. Each result retains its own target and data subset.',
                      flag=bool(n_missing or n_full))
    elif mic_only:
        status.update(n_traits=len(cz_agents), lead=f'{len(cz_agents)} antimicrobials with interval-censored MIC results',
                      gates_line='MIC analysis statuses are reported per antimicrobial.')
    display_groups: Any = groups_used if groups_used is not None else "an unrecorded number of"
    group_note = ""
    if profile_only:
        group_counts = [r["n_groups"] for r in profiles.values() if isinstance(r.get("n_groups"), int)]
        if group_counts:
            lower, upper = min(group_counts), max(group_counts)
            display_groups = str(lower) if lower == upper else f"{lower} to {upper}"
            group_note = " (agent-specific retained subsets)"

    S: List[Section] = []
    fig_no = [0]

    def fig_label() -> str:
        fig_no[0] += 1
        return "Figure %d" % fig_no[0]

    # ------------------------------------------------ 0 measurement summary
    s = Section("summary", "Measurement summary")
    if mic_only:
        s.p('**Quantity measured.** Interval-censored MIC variation associated with recorded lineage. '
            'The recorded panel and coarsening assumptions determine the interval interpretation.')
    elif profile_only:
        s.p("**Quantity measured.** Population liability ICC under a Gaussian "
            "random-intercept binomial-probit model. Classical attribution was disabled. "
            "This is a separate population target; finite-grid model calibration is "
            "not a universal coverage guarantee.")
    else:
        s.p("**Quantity measured.** The share of the variation in the recorded binary outcome "
            "across this collection that lineage membership accounts for, read "
            "from a lineage label and an interpreted result and scored on "
            "isolates the estimator did not see.")
    s.kv([
        ("Cohort", "%s isolates in %s lineages%s"
         % (record.get("n_isolates"), display_groups,
            (", typed by `%s`" % md["lineage_column"])
            + group_note if md.get("lineage_column") else "; no lineage column was given" + group_note)),
        ("Antimicrobials read", str(status['n_traits'])),
    ] + ([("Membership shares estimable", "%d of %d" % (n_scored - n_refused, n_all))]
         if n_scored else []))
    s.p("**Gates.** %s A withheld value failed a declared reporting condition "
        "of its estimator; the recorded reason distinguishes input limitations "
        "from numerical failure. Passing a gate does not verify the model."
        % status["gates_line"])
    if canonical:
        s.table('Table 0', 'Available analyses and their data subsets. A full range is a completed, uninformative result; an unavailable value is not zero.',
                [('Agent', 'lab'), ('Analysis', 'lab'), ('N', 'num'), ('Status', 'lab'), ('Data scope', 'lab')],
                [(r['agent'], r['analysis'], r['n'], r['status'], r['analysis_scope']) for r in canonical])
        attention = [r for r in canonical if r['status'] != 'computed' or r.get('reason')]
        if attention:
            s.table('Table 0b', 'Reasons and next checks for results requiring interpretation. These checks do not change the recorded status or justify changing reporting gates.',
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
        s.table("Table 1", "Every trait, ordered by share. The reading in the "
                "last column is fixed by two conditions the record holds: "
                "whether the interval excludes zero and whether the "
                "Benjamini-Yekutieli step-up on the permutation p-values selected "
                "the trait at a false-discovery level of %s, which holds whatever "
                "the dependence between traits."
                % num(selection.get("alpha"), 2),
                [("Trait", "lab"), ("Share", "num"), ("95 % interval", "num"),
                 ("Control", "num"), ("e-value", "num"), ("", "flag"),
                 ("Reading", "lab")],
                [(r["name"], num(r["pt"]), interval(r["lo"], r["hi"]),
                  num(r["null"], 2),
                  evalue((per_e.get(r["name"]) or {}).get("e_value")),
                  marks(r), reading(r)) for r in rows])
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
                 ("Reading", "lab")],
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
    s.table("Table 2", "Conditions the estimator requires before any result "
            "is reported.",
            [("Condition", "lab"), ("Observed", "num"), ("Required", "num"),
             ("Verdict", "lab")],
            [("Lineages with at least two isolates, fewest over the antimicrobials read",
              str(groups_used) if groups_used is not None else "not recorded",
              "≥ 2",
              "accepted" if groups_used is not None and groups_used >= 2
              else "refused" if groups_used is not None else "not evaluated"),
             ("Support, lowest over the antimicrobials read",
              pct(support) if support is not None else "not recorded",
              "reported",
              "singletons set aside" if support is not None else "not evaluated")])
    feasibility = (input_qc or {}).get("agent_feasibility") or {}
    if feasibility:
        from .qc import FEASIBILITY_HEADERS, FEASIBILITY_SCOPE, feasibility_rows
        s.p(FEASIBILITY_SCOPE)
        s.table("Table 2a", "Per-agent input feasibility.",
                [(name, "lab") for name in FEASIBILITY_HEADERS],
                feasibility_rows(feasibility))
    S.append(s)

    # ------------------------------------------------ 2 the share in detail
    s = Section("share", "The clonal share, trait by trait")
    if rows:
        n_show = min(len(rows), 18)
        shown = rows[:n_show]
        n_cross = sum(1 for r in shown if r["lo"] <= 0.0 <= r["hi"])
        s.p("The intervals have a nominal 95 %% level; their measured coverage "
            "depends on the cohort design and model assumptions. For "
            "**%d of the %d traits** shown the interval includes zero, so no "
            "lineage effect is distinguishable from none for that trait%s."
            % (n_cross, len(shown),
               ("; over the whole panel this holds for **%d of %d**"
                % (n_cross_all, len(rows))) if len(rows) > n_show else ""))
        s.p("The control column is the same estimator run on shuffled lineage "
            "labels; it should sit near zero, and a share is read against it "
            "rather than against zero.")
        with_species = any(finite(r.get("sp_lo")) and finite(r.get("sp_hi"))
                           for r in shown)
        s.figure(fig_label(), "intervals", shown,
                 "Clonal share by trait, point estimate with 95 %% interval, "
                 "%d largest of %d. Traits whose interval crosses zero are "
                 "drawn in grey and marked ‡.%s The intervals are drawn as "
                 "computed. The quantity lies between 0 and 1, so the part of "
                 "an interval below zero carries no information: reading each "
                 "interval as its overlap with that range leaves the coverage "
                 "unchanged, and a lower limit at or below zero means the same "
                 "thing either way."
                 % (n_show, len(rows),
                    " The thin line beneath each interval is the species "
                    "interval of Table 3, for a fresh draw of lineages; it is "
                    "floored at zero by construction." if with_species else ""))
        grid = _load_validation_grid()
        cov = ((grid or {}).get("coverage") or {}).get("clonal_share", {}).get(
            "binary") if grid else None
        if cov:
            s.p("On the release's validation grid, the interval this estimator "
                "prints for a binary trait contained the truth in %s of %s "
                "runs over %s simulated cohorts (by cohort, %s to %s). That "
                "figure belongs to the release, not to this run; it is what "
                "the phrase \"95 %% interval\" was measured to mean."
                % (num(cov.get("mean"), 3), "{:,}".format(cov.get("runs", 0)),
                   cov.get("cells"), num(cov.get("min_cell"), 3),
                   num(cov.get("max_cell"), 3)))
        if any(finite(r.get("sp_lo")) and finite(r.get("sp_hi")) for r in rows):
            n_lineages = groups_used
            s.p("Two intervals answer two questions. The interval above is for "
                "the share the lineages in this collection carry. The second "
                "interval below is for the share a fresh draw of lineages "
                "from the species would show, stated on the scale of the "
                "realised share of Table 3b below; it adds the sampling of "
                "the lineages themselves, on %s degrees of freedom, is "
                "widened to the envelope of the first taken on that scale "
                "(with a lower end no smaller than zero, since a "
                "species share is not negative), and is the one "
                "to quote when the figure is read as a property of the "
                "species rather than of this collection. With few lineages "
                "it is markedly wider; with many the two nearly coincide."
                % (str(int(n_lineages) - 1)
                   if n_lineages is not None and finite(n_lineages) else
                   "the number of lineages less one"))
            few = ((grid or {}).get("coverage") or {}).get("clonal_share", {}).get("few_lineages")
            if any(r.get("few") for r in rows):
                s.p("This collection holds fewer than ten lineages. %s"
                    "A resample of so few lineages has too few distinct "
                    "outcomes, and the species interval, wider in the "
                    "normal-effect cells, answers a different question and "
                    "does not repair uncertainty about this collection."
                    % ("On the release's validation grid the interval for lineage "
                       "membership contained the truth in only %s to %s of runs at "
                       "%d lineages. " % (num(few["min_cell"], 2), num(few["max_cell"], 2),
                                          int(few["lineages"])) if few else ""))
            dominated = [r for r in rows
                         if finite(r.get("dominant")) and r["dominant"] > 0.5]
            if dominated:
                s.p("For %s one lineage carries more than half of the "
                    "between-lineage variation (%s). The species interval "
                    "describes lineage effects drawn from one law, and a "
                    "collection in which one lineage carries the positive calls "
                    "is not that: on the validation grid a carrier law of "
                    "this kind took the species interval below its level "
                    "with few lineages while the interval for the lineages "
                    "in hand held. For %s read the first interval as the "
                    "statement about this collection and the species "
                    "interval with that reservation."
                    % ("one trait" if len(dominated) == 1 else
                       "%d traits" % len(dominated),
                       "; ".join("%s %s" % (r["name"], pct(r["dominant"]))
                                 for r in dominated[:6])
                       + ("; and others" if len(dominated) > 6 else ""),
                       "it" if len(dominated) == 1 else "them"))
            with_latent = any(finite(r.get("latent")) for r in rows)
            if with_latent:
                s.p("The last column is an exploratory, model-equivalent "
                    "restatement under a Gaussian probit threshold model. "
                    "It transforms the collection's observed-scale share "
                    "and interval at the estimated prevalence, treated as "
                    "fixed. It does not estimate the actual latent variance "
                    "of these particular lineages or a generic binomial "
                    "mixed-model intraclass correlation. Monotonicity "
                    "preserves coverage only for a matching target at a "
                    "known, fixed prevalence; it does not guarantee coverage "
                    "after estimating prevalence. Read the observed-scale "
                    "interval as primary. Values at or below zero are "
                    "displayed at the zero boundary.")
            cols = [("Trait", "lab"), ("Share", "num"),
                    ("95 % interval, these lineages", "num"),
                    ("Species share (95 % interval)", "num")]
            if with_latent:
                cols.append(("Latent restatement (exploratory interval)", "num"))
            s.table("Table 3", "The two intervals for every trait: for "
                    "lineage membership in this collection and for the "
                    "species%s."
                    % (", with the share on the latent scale" if with_latent
                       else ""),
                    cols,
                    [(r["name"], num(r["pt"]), interval(r["lo"], r["hi"]),
                      ("%s (%s)" % (num(r["sp_pt"]), interval(r["sp_lo"], r["sp_hi"]))
                       if finite(r.get("sp_pt")) else interval(r.get("sp_lo"), r.get("sp_hi"))))
                     + ((("%s (%s)" % (num(r["latent"]),
                                       interval(r["latent_lo"], r["latent_hi"]))
                          if finite(r.get("latent")) else "not a single call"),)
                        if with_latent else ())
                     for r in rows])
        realised = md.get("realised_share") or {}
        if realised:
            n_open = sum(1 for v in realised.values() if v.get("estimable"))
            s.p("**The realised share of the same call.** The variance-"
                "component ratio of the lineages in hand, with an interval "
                "that is exact under a Gaussian within-lineage law. Its gate "
                "reads the excess kurtosis of the within-lineage residuals "
                "and withholds the interval above 0.99. Rare binary calls "
                "often fail this check; prevalence alone does not determine "
                "the residual kurtosis. The interval is approximate on "
                "binary data even when the gate opens; "
                "the estimate is printed either way, the interval only where "
                "the gate opened, %s."
                % ("for %d of %d traits here" % (n_open, len(realised))
                   if n_open else "which it did for no trait here"))
            s.table("Table 3b", "The realised share per trait and the "
                    "verdict of its kurtosis gate.",
                    [("Trait", "lab"), ("Realised share", "num"),
                     ("Gaussian-model 95 % interval", "num"),
                     ("Residual excess kurtosis", "num"), ("Gate", "lab")],
                    [(name, num(v.get("kappa")),
                      interval(v.get("ci_low"), v.get("ci_high"))
                      if v.get("estimable") else "withheld",
                      num(v.get("residual_excess_kurtosis"), 2),
                      "open" if v.get("estimable") else "closed")
                     for name, v in sorted(realised.items(),
                                           key=lambda kv: -(kv[1].get("kappa")
                                                            if finite(kv[1].get("kappa"))
                                                            else -9))])
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
            "recomputed each time loses its guarantee. The e-value is evidence "
            "on a scale made for that: this run's e-value is a statement about "
            "this cohort, and a programme that adds an intake each year "
            "multiplies the e-value of each new intake, scored against the "
            "lineage rates learned from the earlier ones, into a running "
            "product whose guarantee holds at whatever intake it is read "
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
            s.p("This cohort also carries intakes, %s of them read in the "
                "order of %s, from %s to %s isolates each%s. Each intake is "
                "scored against the lineage rates learned from the intakes "
                "before it, and the running product is the sequential "
                "e-value: it may be read after any intake without spending "
                "the guarantee. %s of %s traits are selected on the product, "
                "against %s at this single look, which is the price of a "
                "guarantee that holds at whatever intake the programme is "
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
        s.p("Where a dilution was recorded, the share is also read from the "
            "dilution itself rather than from the call derived from it. A "
            "call and a dilution are both intervals on one concentration "
            "scale, so the dilution carries more of the reading; a reading at "
            "an end well is censored and enters as an interval. This table "
            "reads the share of that scale carried by lineage, beside the "
            "share of the binary call above. These are different quantities "
            "on independently retained readable and typed subsets; their denominators may differ.")
        crow = []
        for name, e in sorted(cz_agents.items()):
            sh = e.get("share") or {}
            cal = e.get("calibrated") or {}
            crow.append((name, num(sh.get("kappa")),
                         interval(sh.get("ci_low"), sh.get("ci_high")),
                         ("%s (%s)" % (num(cal.get("estimate")), interval(cal.get("low"), cal.get("high")))
                          if cal.get("estimable") else
                          ("not computed" if not cal else "not reportable")),
                         interval(sh.get("realised_low"),
                                  sh.get("realised_high")),
                         pct(sh.get("share_censored")),
                         _shape_reading(cal.get("shape_check")),
                         "yes" if sh.get("estimable")
                         else "no: %s" % sh.get("reason", "")))
        s.table("Table 5", "Share of the dilution scale carried by lineage, per "
                "agent, from the interval-censored likelihood. The share and its "
                "approximate F interval come from the moment iteration; that interval "
                "undercovered in simulation. The calibrated column gives the maximum-"
                "likelihood share with the interval validated by simulation for the "
                "Gaussian population model, subject to its sampling and coarsening "
                "assumptions. The realised interval is for the lineages this cohort holds; "
                "here it is an approximation built on the fitted mean-square ratio and "
                "not the exact interval of uncensored data.",
                [("Agent", "lab"), ("Share", "num"), ("Approximate F interval", "num"),
                 ("Calibrated share and 95 % interval", "num"),
                 ("Realised interval", "num"), ("Censored readings", "num"),
                 ("Gaussian model check", "lab"), ("Estimable", "lab")], crow)
        shape_rejected = sorted(name for name, e in cz_agents.items()
                                if ((e.get("calibrated") or {}).get("shape_check") or {}).get("rejected"))
        if shape_rejected:
            s.p("The check of the Gaussian model rejects it for %s. The calibrated "
                "interval was validated by simulation for Gaussian readings; for "
                "readings of another shape, such as a wild-type and a non-wild-type "
                "mode, its coverage has not been established and the interval is "
                "read as approximate. Where calls against an external cut-off are "
                "available, the lineage question for such an agent is read from the "
                "calls: the share of the binary call above and, if enabled, the "
                "population model." % ", ".join(shape_rejected))
        call_pt = {r["name"]: r["pt"] for r in rows}
        drows = []
        for name, e in sorted(cz_agents.items()):
            sh = e.get("share") or {}
            cal = e.get("calibrated") or {}
            # The calibrated interval is the validated reading of this quantity.
            # The approximate F interval undercovered in simulation, so it is
            # drawn only for an agent that has no calibrated interval, and it
            # is drawn dashed so that the two are never read as one series.
            validated = (bool(cal.get("estimable")) and finite(cal.get("estimate"))
                         and finite(cal.get("low")) and finite(cal.get("high")))
            point = cal.get("estimate") if validated else sh.get("kappa")
            if finite(point):
                drows.append({"name": name, "pt": point,
                              "lo": cal.get("low") if validated else sh.get("ci_low"),
                              "hi": cal.get("high") if validated else sh.get("ci_high"),
                              "calibrated": validated,
                              "r_lo": sh.get("realised_low"),
                              "r_hi": sh.get("realised_high"),
                              "censored": sh.get("share_censored"),
                              "estimable": bool(sh.get("estimable")),
                              "call": call_pt.get(name)})
        drows.sort(key=lambda r: -r["pt"])
        shown = drows[:18]
        n_fallback = sum(1 for r in shown if not r["calibrated"])
        if not shown or n_fallback == len(shown):
            method_note = (
                "No calibrated interval was computed on this run, so every "
                "interval drawn is the approximate F interval beside the moment "
                "estimate; that interval undercovered in simulation")
        elif n_fallback:
            method_note = (
                "A solid interval is the calibrated interval, validated by "
                "simulation for the Gaussian population model, drawn beside the "
                "maximum-likelihood estimate; a dashed interval is the "
                "approximate F interval beside the moment estimate, drawn for "
                "the %d %s with no calibrated interval, and it "
                "undercovered in simulation"
                % (n_fallback, "agent" if n_fallback == 1 else "agents"))
        else:
            method_note = (
                "Every interval drawn is the calibrated interval, validated by "
                "simulation for the Gaussian population model, beside the "
                "maximum-likelihood estimate")
        if drows:
            s.figure(fig_label(), "dilution", shown,
                     "Share of the dilution scale carried by lineage, per "
                     "agent, %d largest of %d: the point estimate with its 95 %% "
                     "interval (thin) and the realised "
                     "interval for these lineages (thick); an agent whose "
                     "reading is refused is drawn in grey. %s. The hollow diamond "
                     "is the share of the binary call from Table 1 for the "
                     "same agent, a different quantity whose retained isolates may differ, "
                     "placed here so that the two readings can be seen side "
                     "by side. The figure in parentheses at the right is the "
                     "share of readings at an end well, which enter as "
                     "censored."
                     % (min(len(drows), 18), len(drows), method_note))
        _calibration_detail(s, cz_agents, fig_label,
                            bool(((cfg.get("censored") or {}).get("calibrated_interval", True))))
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
            "seen only in the pooled run and in no stratum is carried by the "
            "differences between strata."
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
                cal = (mics.get(name) or {}).get("calibrated") or {}
                st_rows.append((str(lv), name, num(c.get("n"), 0) if finite(c.get("n")) else "",
                                num(c.get("kappa_adj")) if finite(c.get("kappa_adj")) else "not computed",
                                interval(c.get("ci_low"), c.get("ci_high"))
                                if finite(c.get("ci_low")) else "",
                                ("%s (%s)" % (num(cal.get("estimate")), interval(cal.get("low"), cal.get("high"))))
                                if cal.get("estimable") else (cal.get("status") or "not computed")))
        s.table("Table 5d", "Clonal share per stratum and agent: the share of the "
                "binary call with its interval, and the calibrated share of the "
                "dilution scale where a MIC table was supplied. The same rows are "
                "written to strata_results.csv.",
                [("Stratum", "lab"), ("Agent", "lab"), ("n", "num"),
                 ("Call share", "num"), ("95 % interval", "num"),
                 ("Calibrated MIC share (95 % interval)", "num")], st_rows)
        S.append(s)

    # ------------------------------------------ 5 composition or rate
    s = Section("decomposition", "A change of lineages or a change within them")
    sv = summary.get("surveillance") or {}
    dec = sv.get("decomposition")
    if dec:
        n_ref = int(dec.get("n_within_lineage_refused") or 0)
        n_feat = int(dec.get("n_features") or 0)
        s.p("Contrast %s: of %s traits, %s have a detected composition "
            "component (a change in lineage mix), and %s have a detected "
            "within-lineage component (a change in rate), after "
            "false-discovery control within each component family."
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

        s.table("Table 6", "The prevalence difference of each trait split into "
                "a change in lineage composition and a change in rate within "
                "lineages. The two components sum to the difference. The "
                "limits are bootstrap percentiles and are not adjusted for "
                "multiplicity; the last column names the components the "
                "false-discovery procedure selected.",
                [("Trait", "lab"), ("Difference", "num"), ("Composition", "num"),
                 ("95 % interval", "num"), ("Within lineage", "num"),
                 ("95 % interval", "num"), ("Selected", "lab")],
                [(name, num(block.get("difference")),
                  num(block.get("composition")), limits(block, "composition_ci95"),
                  num(block.get("within_lineage")),
                  limits(block, "within_lineage_ci95"),
                  selected(block) if block.get("status") == "ok"
                  else str(block.get("status")))
                 for name, block in sorted(per_feature.items())])
    else:
        s.p("A prevalence difference between two collections can be split "
            "into a change in lineage composition and a change in rate within "
            "lineages. That decomposition needs two collections and was not "
            "run here: this record holds one.")
    S.append(s)

    if general_profiles:
        s = Section("general_population", "Population liability ICC: general model method")
        s.p("For two or more repeated groups, the hull-coverage target is 95%. "
            "That branch uses a fitted-nuisance bootstrap with a conservative internal "
            "0.04 threshold. Neither fitted-nuisance nor continuous-parameter coverage "
            "is an exact guarantee. The primary result is the hull; all accepted "
            "components remain available as diagnostics.")
        accepted_general = all(r.get("validation_status") == "validation_accepted" and
                               ((r.get("provenance") or {}).get("protocol") or {}).get("accepted") is True
                               for r in general_profiles.values())
        s.p("The finite-panel simulation study validated the coverage target of this method."
            if accepted_general else "The coverage of this method is not validated for one or more result records.")
        s.p("With one repeated group, this method reports a one-sided upper limit. "
            "It has a separate 95% upper-limit target. That branch is a declared "
            "reporting policy. When every group is a "
            "singleton, the model cannot identify ICC and retains the full [0,1] range. "
            "The upper limit uses a separate finite-Monte-Carlo construction with "
            "coverage averaged over its random draws under the stated model.")
        s.p("The model assumes independent Gaussian group effects, conditional "
            "binomial counts, fixed or noninformative group sizes and ignorable "
            "selection. Results describe a population liability ICC, separate from "
            "the collection's lineage-membership share and its reporting gates.")
        kinds = {"model_confidence_set": "Model confidence set", "one_sided_upper_limit": "Upper confidence limit",
                 "structural_full_range": "Insufficient within-group information", "incomplete": "Calculation incomplete",
                 "unavailable": "Result unavailable"}
        def general_reading(r):
            if r.get("empty"):
                return "No parameter value retained"
            if not r.get("complete"):
                return "Calculation incomplete" if r.get("confidence_kind") == "incomplete" else "Result unavailable"
            if r.get("boundary_fill_applied"):
                return "Upper limit at zero (boundary convention)"
            if r.get("ci_low") == 0 and r.get("ci_high") == 1:
                return "Full range retained"
            return ("Computed; empirical validation accepted"
                    if r.get("validation_status") == "validation_accepted" and
                    ((r.get("provenance") or {}).get("protocol") or {}).get("accepted") is True
                    else "Computed; validation not accepted")
        s.table("Table G1", "Completed confidence objects and explicit computation status.",
                [("Agent", "lab"), ("Object", "lab"), ("Limits", "num"), ("Groups", "num"), ("Reading", "lab")],
                [(name, kinds.get(r.get("confidence_kind"), "Model result"),
                  "empty set" if r.get("empty") else interval(r.get("ci_low"), r.get("ci_high")),
                  str(r.get("n_groups", "unrecorded")), general_reading(r)) for name, r in general_profiles.items()])
        for name, r in general_profiles.items():
            if r.get("failure_reason"):
                s.p("%s: %s." % (name, r["failure_reason"]))
            if r.get("boundary_fill_applied"):
                s.p("%s: no raw upper-bound grid point was accepted. The reported [0,0] "
                    "set follows the declared boundary convention; it is distinct from "
                    "a genuinely empty bootstrap set." % name)
        first_general = next(iter(general_profiles.values()))
        protocol = (first_general.get("provenance") or {}).get("protocol") or {}
        s.kv([("Validation status", str(first_general.get("validation_status"))),
              ("Selected protocol", str(protocol.get("protocol_id"))),
              ("Protocol SHA-256", str(protocol.get("packaged_protocol_sha256")))])
        s.p("Numerically unresolved work has no completed confidence limits. The result "
            "record preserves provisional components, errors, exact retained counts, "
            "method streams and cache provenance. Memory settings are array planning "
            "allowances, not operating-system memory limits.")
        S.append(s)

    if fixed_cutoff_profiles:
        s = Section("population_probit", "Population liability ICC")
        first = next(iter(fixed_cutoff_profiles.values()))
        s.p("This optional result targets the Gaussian random-intercept population "
            "liability ICC, with marginal prevalence jointly fitted under a "
            "conditional binomial model. Independent Gaussian group effects and "
            "group sizes fixed independently of those effects are assumptions. "
            "Selection into the modeled population must be ignorable. The model "
            "does not identify transmission or intervention effects.")
        s.table("Table P1", "Population-model profile sets; finite-grid calibration "
                "is not a universal 95% coverage guarantee. These results are "
                "separate from classical lineage-membership estimates and admissions.",
                [("Agent", "lab"), ("Population liability ICC", "num"), ("Model profile interval", "num"),
                 ("Fitted prevalence", "num"), ("Gaussian law check", "lab"), ("Status", "lab")],
                [(name, num(r.get("rho_hat")), interval(r.get("ci_low"), r.get("ci_high")),
                  num(r.get("prevalence_hat")), _shape_reading(r.get("mixing_check")),
                  r.get("status", "not recorded"))
                 for name, r in fixed_cutoff_profiles.items()])
        flagged = [name for name, r in fixed_cutoff_profiles.items()
                   if (r.get("mixing_check") or {}).get("rejected")]
        if flagged:
            s.p("The check of the Gaussian law of lineage effects rejects it for %s. "
                "Lineage effects that fall into a few classes, such as a resistant "
                "clone among susceptible lineages, cause this; the interval's "
                "coverage has not been established for such data." % _join(flagged))
        s.kv([("Fixed LR cutoff", str(first.get("critical_value"))),
              ("Calibration", str(first.get("calibration_id"))),
              ("Calibration SHA-256", str(first.get("calibration_sha256"))),
              ("Independent validation", str(first.get("validation_id"))),
              ("Validation status", str(first.get("validation_status"))),
              ("Validation summary SHA-256", str(first.get("validation_summary_sha256")))])
        s.p("Constant outcomes or too few repeated groups retain an explicitly "
            "uninformative [0,1] set with no point estimate. Numerical failure "
            "has no estimate or interval. The repeated-group policy is separate "
            "from the clonal share, which sets singleton lineages aside. Counts, exclusions, numerical "
            "checks and the complete model-domain diagnostics are in the result record.")
        for name, r in fixed_cutoff_profiles.items():
            domain = r.get("domain_diagnostics") or {}
            warnings = domain.get("extrapolation_warnings") or []
            if warnings:
                s.p("%s: extrapolation warning — %s." % (name, "; ".join(warnings)))
            if r.get("failure_reason"):
                s.p("%s: %s." % (name, r["failure_reason"]))
        domain = first.get("domain_diagnostics") or {}
        s.p("Tested finite-grid group-count range: %s; maximum tested group size: %s; "
            "model prevalence range: %s. Range inclusion does not establish "
            "validation for every design within these bounds." % (
                domain.get("tested_group_count_range"), domain.get("tested_max_group_size"),
                domain.get("tested_prevalence_range")))
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
    interval_terms = []
    if not profile_only:
        interval_terms.append(
            "**Collection-bootstrap interval.** The lineage-membership interval "
            "targets the sampled collection. Separate Gaussian-model component "
            "and species intervals retain their own targets and assumptions; "
            "coverage measured on a finite simulation grid is not a universal guarantee.")
    if fixed_cutoff_profiles:
        interval_terms.append(
            "**Population-model interval.** The grouped-probit profile targets "
            "a Gaussian population liability ICC, using its recorded fixed cutoff "
            "and finite-grid calibration/validation status. It does not target "
            "the sampled collection's membership share.")
    if general_profiles:
        interval_terms.append(
            "**General population-model confidence object.** A 95% target using an "
            "approximate fitted-nuisance bootstrap with internal threshold 0.04, "
            "or a separately defined one-sided upper limit under the repeated-group "
            "policy. The method-specific validation status is recorded in the population section. A full set, an empty set and an "
            "incomplete calculation are different outcomes.")
    s.callout("Terms used in this report", [
        "**Share.** How much of the variation of the call between isolates "
        "is associated with their recorded lineage. Near 1: a strong "
        "association on the measured scale. Near 0: little association on "
        "that scale. This does not identify transmission or its mechanism.",
    ] + interval_terms + [
        "**Control.** The same calculation on shuffled lineage labels, which "
        "should return roughly zero; a share is read against it.",
        "**Support.** The share of isolates in lineages with at least two "
        "members. A singleton lineage cannot be predicted out of sample, so the "
        "share is scored on the other isolates and the singletons are set "
        "aside; the share then speaks for singleton lineages only if they are "
        "drawn from the same population of lineages as the repeated ones.",
        "**Gate.** A condition fixed before the run. If it does not hold, the "
        "software withholds the estimate and names the condition that failed.",
        "**e-value.** A fixed-look e-value and a sequential e-process have "
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
    if not profile_only:
        colophon.append("Classical collection-bootstrap endpoints are printed "
                        "without clipping; the separately labeled species interval "
                        "is floored at zero as part of its construction.")
    if fixed_cutoff_profiles:
        colophon.append("Population-model profile endpoints are computed within "
                        "[0,1]. A numerical failure has missing endpoints, and an "
                        "uninformative full set is labeled separately. Its population "
                        "target and finite-grid evidence are independent of the "
                        "classical collection-bootstrap interval.")
    if general_profiles:
        colophon.append("General population-model limits lie in [0,1]. The selected "
                        "95% target uses an approximate bootstrap with internal threshold "
                        "0.04, or a separately labeled one-sided upper limit. The recorded method-specific validation status applies. "
                        "Empty sets, full ranges and incomplete calculations "
                        "are reported separately; they do not alter classical admissions.")
    return Report(ident=ident, status=status, sections=S, symbols=symbols,
                  colophon=colophon)


#: lineages drawn per agent in the heat map, the largest first
HEATMAP_LINEAGES = 14


def _calibration_detail(s: "Section", cz_agents: dict, fig_label, switched_on: bool) -> None:
    """The panel each laboratory tested, the fixed-effect coefficients and the
    status of every calibrated interval, followed by the lineage by dilution
    counts as a heat map. ``switched_on`` says whether the run asked for the
    calibrated interval, so that an agent without one is read correctly."""
    geo_rows, coef_rows, status_rows = [], [], []
    for name, e in sorted(cz_agents.items()):
        g = e.get("panel") or {}
        per = g.get("per_panel") or ({"all readings": g} if g else {})
        cal = e.get("calibrated") or {}
        cuts = cal.get("panel_edges")
        for lab, geo in per.items():
            edges = (cuts.get(lab) if isinstance(cuts, dict) else cuts) or []
            geo_rows.append((name, str(lab), str(geo.get("n_wells", "")),
                             "%s to %s" % (num(geo.get("lowest"), 4), num(geo.get("highest"), 4)),
                             ("%s cut points" % len(edges)) if edges else "none",
                             "yes" if geo.get("doubling") else "no",
                             _join([str(m) for m in geo.get("admissible_modes") or ()]) or "none"))
        for block_ in cal.get("adjustment") or []:
            coefs = block_.get("coefficients") or {}
            for level, value in coefs.items():
                coef_rows.append((name, str(block_.get("column", "")), str(block_.get("reference_level")),
                                  str(level), num(value)))
        if cal:
            reason = cal.get("reason") or ("interval reported" if cal.get("estimable") else "")
        elif not switched_on:
            reason = "calibrated_interval is off in the configuration"
        else:
            reason = (e.get("share") or {}).get("reason") or "the agent did not reach the estimator"
        status_rows.append((name, str(cal.get("status") or "not computed"),
                            str(cal.get("method") or ""), str(cal.get("n_boot") or ""),
                            str(cal.get("bootstrap_failed", "")), reason))
    if geo_rows:
        s.p("The calibrated interval reads each laboratory on the wells that "
            "laboratory tested. The table below gives, per agent and laboratory, "
            "the panel the readings were taken to come from and the cut points "
            "the model uses; an end-well reading is censored at the panel edge "
            "of its own laboratory, not at the widest edge in the collection.")
        s.table("Table 5a", "Panel geometry per agent and testing laboratory. Cut "
                "points are the boundaries between adjacent wells on which the "
                "likelihood is written; the admissible modes say which readings "
                "of the end wells the panel supports.",
                [("Agent", "lab"), ("Laboratory", "lab"), ("Wells", "num"),
                 ("Range", "num"), ("Cut points", "num"), ("Doubling", "lab"),
                 ("Admissible modes", "lab")], geo_rows)
    if coef_rows:
        s.p("Where a covariate was declared, its levels enter the model as "
            "fixed effects on the log2 MIC. The first level in sorted order is "
            "the reference; each coefficient is the shift of a level from the "
            "reference, in log2 dilution steps, after the lineage effect. A "
            "coefficient near one doubling is a one-well shift between laboratories "
            "or countries that would otherwise be read as a lineage difference.")
        s.table("Table 5b", "Fixed-effect coefficients of the declared covariates, per "
                "agent, in log2 dilution steps relative to the reference level.",
                [("Agent", "lab"), ("Covariate", "lab"), ("Reference", "lab"),
                 ("Level", "lab"), ("Coefficient", "num")], coef_rows)
    if status_rows:
        s.table("Table 5c", "Status of the calibrated interval per agent. The "
                "status names the outcome of the interval inversion; the reason "
                "is given whenever no interval is reported.",
                [("Agent", "lab"), ("Status", "lab"), ("Method", "lab"),
                 ("Draws", "num"), ("Failed refits", "num"), ("Reason", "lab")],
                status_rows)
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
                   "lineages carry. Read the selection as the finding and the "
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
