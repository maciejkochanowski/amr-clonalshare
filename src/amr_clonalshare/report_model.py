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
recorded resolution adds, what the carriage pattern says for surveillance; then
provenance.

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
from .attribution import SUPPORT_THRESHOLD

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
    unread: List[tuple] = []
    for name, v in (md.get("clonal_share") or {}).items():
        if finite(v.get("kappa_adj")) and finite(v.get("ci_low")):
            rows.append({"name": name, "pt": v["kappa_adj"], "lo": v["ci_low"],
                         "hi": v["ci_high"], "null": v.get("null_mean"),
                         "sp_lo": v.get("superpopulation_low"),
                         "sp_hi": v.get("superpopulation_high"),
                         "p": v.get("p_value"), "support": v.get("support"),
                         "estimable": v.get("estimable"),
                         "latent": v.get("latent_share"),
                         "latent_lo": v.get("latent_low"),
                         "latent_hi": v.get("latent_high"),
                         "dominant": v.get("dominant_lineage_share"),
                         "few": bool(v.get("few_lineages"))})
        else:
            # a trait the estimator could not score at all: one lineage, a
            # constant call, or too few isolates; it is listed with the
            # reason rather than dropped from the report
            n_groups = v.get("n_groups")
            if isinstance(n_groups, int) and n_groups < 2:
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
    per_e = (md.get("lineage_evidence") or {}).get("per_feature") or {}
    seq_e = (md.get("lineage_evidence") or {}).get("sequential") or {}

    shares = [v for v in (md.get("clonal_share") or {}).values()
              if finite(v.get("support"))]
    support = min((v["support"] for v in shares), default=None)
    groups_used = max((v.get("n_groups") for v in shares
                       if isinstance(v.get("n_groups"), int)), default=None)

    cz = md.get("censored_share") or {}
    cz_agents = cz.get("per_agent") or {}
    cz_join = cz.get("join") or {}

    def reading(r) -> str:
        crosses = r["lo"] <= 0.0 <= r["hi"]
        if r.get("estimable") is False:
            return "refused"
        if crosses and r["name"] in rejected:
            return "lineage effect selected by e-BH; its size is not resolved"
        if crosses:
            return "no detectable lineage effect"
        if r["name"] in rejected:
            return "evidence of a lineage effect"
        return "interval excludes zero; not selected by e-BH"

    def marks(r) -> str:
        m = ""
        if r["lo"] <= 0.0 <= r["hi"]:
            m += "‡"
        if finite(r.get("support")) and r["support"] < gate:
            m += "†"
        return m

    n_refused = sum(1 for r in rows if r.get("estimable") is False)
    gate_recorded = ((input_qc or {}).get("lineage") or {}).get("support_threshold")
    gate = (float(gate_recorded) if isinstance(gate_recorded, (int, float))
            and finite(gate_recorded) else float(SUPPORT_THRESHOLD))
    ident = {
        "run": run_id, "issued": stamp,
        "software": "amr-clonalshare %s" % __version__,
        "record": "clonal_share_result.json",
        "isolates": str(record.get("n_isolates")),
        "lineage": str(md["lineage_column"]) if md.get("lineage_column")
        else "none",
        "seed": str(record.get("seed")),
        "configuration": "sha256:%s" % digest,
        "record_sha256": record_sha or "",
    }
    n_all = len(rows) + len(unread)
    status = {
        "n_traits": n_all, "n_refused": n_refused + len(unread),
        "lead": ("%d of %d antimicrobials estimable" % (len(rows) - n_refused,
                                                         n_all))
        if rows else ("no antimicrobial could be scored" if unread
                      else "no antimicrobial read"),
        "gates_line": (("The estimator withheld %d of %d; the condition that "
                        "failed is named per antimicrobial." % (n_refused, n_all)
                        if n_refused else "")
                       + ((" " if n_refused else "")
                          + "%d of %d could not be scored, for a reason named "
                          "below." % (len(unread), n_all) if unread else ""))
        if (n_refused or unread) else ("Every antimicrobial read cleared the "
                                       "support gate." if rows else
                                       "No antimicrobial carried a readable call."),
        "flag": bool(n_refused or unread),
    }

    S: List[Section] = []
    fig_no = [0]

    def fig_label() -> str:
        fig_no[0] += 1
        return "Figure %d" % fig_no[0]

    # ------------------------------------------------ 0 measurement summary
    s = Section("summary", "Measurement summary")
    s.p("**Quantity measured.** The share of the variation in resistance "
        "across this collection that lineage membership accounts for, read "
        "from a lineage label and a susceptibility result and scored on "
        "isolates the estimator did not see.")
    s.kv([
        ("Cohort", "%s isolates in %s lineages%s"
         % (record.get("n_isolates"), groups_used
            if groups_used is not None else "an unrecorded number of",
            (", typed by `%s`" % md["lineage_column"])
            if md.get("lineage_column") else "; no lineage column was given")),
        ("Antimicrobials read", str(n_all)),
    ] + ([("Estimable", "%d of %d" % (len(rows) - n_refused, n_all))]
         if rows else []))
    s.p("**Gates.** %s A withheld value is a completed reading: it says "
        "the collection cannot identify the quantity at this typing "
        "resolution. A value that clears the gate is not thereby correct."
        % status["gates_line"])
    if rows:
        s.table("Table 1", "Every trait, ordered by share. The reading in the "
                "last column is fixed by two conditions the record holds: "
                "whether the interval excludes zero and whether the e-BH "
                "procedure selected the trait at level %s."
                % num(evidence.get("alpha"), 2),
                [("Trait", "lab"), ("Share", "num"), ("95 % interval", "num"),
                 ("Control", "num"), ("e-value", "num"), ("", "flag"),
                 ("Reading", "lab")],
                [(r["name"], num(r["pt"]), interval(r["lo"], r["hi"]),
                  num(r["null"], 2),
                  evalue((per_e.get(r["name"]) or {}).get("e_value")),
                  marks(r), reading(r)) for r in rows])
    not_done = _not_evaluated(record, cz_agents)
    not_done = ["`%s`: not scored; %s." % (name, why) for name, why in unread] + not_done
    if not_done:
        s.callout("Not evaluated on this run, and why", not_done)
    inter = _interpretation(rows, rejected)
    if inter:
        s.callout("Interpretation (generated from Table 1 by fixed rules)", inter)
    S.append(s)

    # --------------------------------------------- 1 admissibility of input
    s = Section("admissibility", "Admissibility of the input")
    prov = []
    if md.get("lineage_column"):
        prov.append(("Lineage column", "`%s`" % md["lineage_column"]))
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
                "outcome: %d of %d. Below that count the estimate is still "
                "computed and its interval says how little it rests on."
                % (input_qc.get("min_minor_count"), weak, len(traits)))
        if lin:
            s.p("Lineage groups in the input: %s, of which %s hold a single "
                "isolate. Support %s against the %s the estimator needs: %s"
                % (lin.get("n_groups"), lin.get("n_singletons"),
                   pct(lin.get("support")), pct(lin.get("support_threshold")),
                   "accepted." if lin.get("estimable")
                   else "not accepted at this typing resolution."))
            sizes = lin.get("group_sizes") or {}
            if isinstance(sizes, dict) and sizes:
                ordered = sorted(((str(k), int(v)) for k, v in sizes.items()),
                                 key=lambda kv: (-kv[1], kv[0]))
                n_show = min(len(ordered), 40)
                s.figure(fig_label(), "lineages", {
                    "sizes": ordered[:n_show], "n_total": len(ordered),
                    "n_isolates": lin.get("n_typed", record.get("n_isolates")),
                    "support": lin.get("support"),
                    "threshold": lin.get("support_threshold"),
                    "min_group": input_qc.get("min_group_size", 2)},
                    "Isolates per lineage, largest first%s. A lineage of one "
                    "isolate, drawn in orange, cannot be predicted out of "
                    "sample and counts against support; support is the share "
                    "of isolates in the other lineages, %s here against the "
                    "%s the estimator requires. The largest lineage holds %s "
                    "of the isolates, which sets how much one lineage can "
                    "weigh in the share."
                    % ((", the %d largest of %d" % (n_show, len(ordered)))
                       if n_show < len(ordered) else "",
                       pct(lin.get("support")), pct(lin.get("support_threshold")),
                       pct(ordered[0][1] / float(lin.get("n_typed") or
                                                  record.get("n_isolates") or 1))))
    else:
        s.p("No input check was attached to this record.")
    s.table("Table 2", "Conditions the estimator requires before any result "
            "is reported.",
            [("Condition", "lab"), ("Observed", "num"), ("Required", "num"),
             ("Verdict", "lab")],
            [("Lineage support, lowest over the antimicrobials read",
              pct(support) if support is not None else "not recorded",
              "≥ " + pct(gate),
              "accepted" if support is not None and finite(support)
              and support >= gate
              else "refused" if support is not None else "not evaluated"),
             ("Lineage groups used",
              str(groups_used) if groups_used is not None else "not recorded",
              "reported",
              "accepted" if groups_used is not None else "not evaluated")])
    S.append(s)

    # ------------------------------------------------ 2 the share in detail
    s = Section("share", "The clonal share, trait by trait")
    if rows:
        n_show = min(len(rows), 18)
        shown = rows[:n_show]
        n_cross = sum(1 for r in shown if r["lo"] <= 0.0 <= r["hi"])
        s.p("Read the intervals as frequencies. If cohorts like this one were "
            "drawn again and again, the interval printed for a trait would "
            "cover that trait's true share on about 95 draws in 100. For "
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
            if any(r.get("few") for r in rows):
                s.p("This collection holds fewer than ten lineages. On the "
                    "release's validation grid the interval for lineage "
                    "membership contained the truth in only 0.83 to 0.96 of "
                    "runs at five lineages, because a resample of so few "
                    "lineages has too few distinct outcomes, while the "
                    "species interval held its level there. With this few "
                    "lineages read the species interval as the statement.")
            dominated = [r for r in rows
                         if finite(r.get("dominant")) and r["dominant"] > 0.5]
            if dominated:
                s.p("For %s one lineage carries more than half of the "
                    "between-lineage variation (%s). The species interval "
                    "describes lineage effects drawn from one law, and a "
                    "collection in which one lineage carries the resistance "
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
                s.p("The last column restates the share on the latent scale "
                    "of a threshold model, the scale on which a mixed model "
                    "with a binomial link reports an intraclass correlation. "
                    "The two scales are one-to-one at a given prevalence and "
                    "the map is monotone, so the interval keeps its coverage; "
                    "the latent figure is larger because a call discards the "
                    "part of the liability that does not cross the threshold. "
                    "A share at or below zero has no latent counterpart and "
                    "is printed as zero.")
            cols = [("Trait", "lab"), ("Share", "num"),
                    ("95 % interval, these lineages", "num"),
                    ("95 % interval, species", "num")]
            if with_latent:
                cols.append(("Latent scale (95 % interval)", "num"))
            s.table("Table 3", "The two intervals for every trait: for "
                    "lineage membership in this collection and for the "
                    "species%s."
                    % (", with the share on the latent scale" if with_latent
                       else ""),
                    cols,
                    [(r["name"], num(r["pt"]), interval(r["lo"], r["hi"]),
                      interval(r.get("sp_lo"), r.get("sp_hi")))
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
                "and withholds the interval above 0.99, which on a binary "
                "call closes wherever the prevalence is far from one half; "
                "the estimate is printed either way, the interval only where "
                "the gate opened, %s."
                % ("for %d of %d traits here" % (n_open, len(realised))
                   if n_open else "which it did for no trait here"))
            s.table("Table 3b", "The realised share per trait and the "
                    "verdict of its kurtosis gate.",
                    [("Trait", "lab"), ("Realised share", "num"),
                     ("Exact 95 % interval", "num"),
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
            % (evidence.get("n_rejected"), evidence.get("n_agents"),
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
                         key=lambda kv: -((kv[1].get("log_e") or [0.0])[-1]))])
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
                    "named at the right; the guarantee holds at whatever "
                    "intake the line is read. The frame is cut at −25: a line "
                    "that leaves it at the bottom has accumulated evidence "
                    "against a lineage effect and does not return in the "
                    "intakes shown."
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
            "share of the binary call above, which is a different quantity "
            "on the same isolates.")
        crow = []
        for name, e in sorted(cz_agents.items()):
            sh = e.get("share") or {}
            crow.append((name, num(sh.get("kappa")),
                         interval(sh.get("ci_low"), sh.get("ci_high")),
                         interval(sh.get("realised_low"),
                                  sh.get("realised_high")),
                         pct(sh.get("share_censored")),
                         "yes" if sh.get("estimable")
                         else "no: %s" % sh.get("reason", "")))
        s.table("Table 5", "Share of the dilution scale carried by lineage, per "
                "agent, from the interval-censored likelihood. The realised "
                "interval is for the lineages this cohort holds; the 95 % "
                "interval is for the species behind it.",
                [("Agent", "lab"), ("Share", "num"), ("95 % interval", "num"),
                 ("Realised interval", "num"), ("Censored readings", "num"),
                 ("Estimable", "lab")], crow)
        call_pt = {r["name"]: r["pt"] for r in rows}
        drows = []
        for name, e in sorted(cz_agents.items()):
            sh = e.get("share") or {}
            if finite(sh.get("kappa")):
                drows.append({"name": name, "pt": sh["kappa"],
                              "lo": sh.get("ci_low"), "hi": sh.get("ci_high"),
                              "r_lo": sh.get("realised_low"),
                              "r_hi": sh.get("realised_high"),
                              "censored": sh.get("share_censored"),
                              "estimable": bool(sh.get("estimable")),
                              "call": call_pt.get(name)})
        drows.sort(key=lambda r: -r["pt"])
        if drows:
            s.figure(fig_label(), "dilution", drows[:18],
                     "Share of the dilution scale carried by lineage, per "
                     "agent, %d largest of %d: point estimate with the 95 %% "
                     "interval for the species (thin) and the realised "
                     "interval for these lineages (thick); an agent whose "
                     "reading is refused is drawn in grey. The hollow diamond "
                     "is the share of the binary call from Table 1 for the "
                     "same agent, a different quantity on the same isolates, "
                     "placed here so that the two readings can be seen side "
                     "by side. The figure in parentheses at the right is the "
                     "share of readings at an end well, which enter as "
                     "censored."
                     % (min(len(drows), 18), len(drows)))
    else:
        s.p("No recorded dilutions were supplied on this run; the shares "
            "above are read from binary calls only.")
    S.append(s)

    # ------------------------------------------ 5 composition or rate
    s = Section("carriage", "How resistance is carried across lineages")
    sv = summary.get("surveillance") or {}
    dep = sv.get("departing_carriage_direction_counts") or {}
    n_read = sv.get("n_features_with_carriage_reading") or sum(
        (sv.get("carriage_direction_counts") or {}).values())
    n_dep = sv.get("n_features_departing_from_proportional_carriage") or 0
    if not n_read:
        s.p("No trait has a carriage reading on this run, so nothing is said "
            "here about whether resistance piles into a few lineages or is "
            "spread in proportion to their size.")
    else:
        parts = [t for t in (
            "%d in fewer lineages than chance gives" % dep["concentrated"]
            if dep.get("concentrated") else "",
            "%d in more" % dep["dispersed"] if dep.get("dispersed") else "",
            "%d in as many carrying lineages as chance gives but not in the "
            "same shares across them" % dep["proportional"]
            if dep.get("proportional") else "") if t]
        shape = _join(parts) or "the record does not break them down"
        s.p("Of the **%d** traits with a carriage reading, **%d** depart from "
            "carriage in proportion to lineage size far enough to matter for a "
            "prevalence reading: %s.%s"
            % (n_read, n_dep, shape,
               (" The remaining %d are consistent with proportional carriage."
                % (n_read - n_dep)) if n_read > n_dep else ""))
    if sv.get("lineage_prevalence_widest_gap_feature"):
        s.p("The widest gap between the per-isolate and the per-lineage "
            "prevalence is %s, on `%s`. Where the two differ, a small number "
            "of large lineages carry most of the resistance, and a change in "
            "prevalence may be a change in which lineages were sampled rather "
            "than a change in rate."
            % (num(sv.get("lineage_prevalence_widest_gap")),
               sv.get("lineage_prevalence_widest_gap_feature")))
    lrp = md.get("lineage_resolved_prevalence") or {}
    tc = md.get("trait_concentration") or {}
    c_rows = []
    for name in sorted(set(lrp) | set(tc)):
        a = lrp.get(name) or {}
        b = tc.get(name) or {}
        if a.get("status") != "ok" or b.get("status") != "ok":
            continue
        ci = a.get("prevalence_per_lineage_ci95") or [None, None]
        nci = b.get("null_effective_number_ci95") or [None, None]
        c_rows.append({"name": name,
                       "per_isolate": a.get("prevalence_per_isolate"),
                       "per_lineage": a.get("prevalence_per_lineage"),
                       "pl_lo": ci[0], "pl_hi": ci[1],
                       "n_carrying": b.get("n_lineages_carrying"),
                       "n_lineages": b.get("n_lineages"),
                       "eff": b.get("effective_number_of_lineages"),
                       "null_lo": nci[0], "null_hi": nci[1],
                       "direction": b.get("direction"),
                       "departs": finite(b.get("p_value"))
                       and b["p_value"] < 0.05,
                       "p": b.get("p_value"), "p_floor": b.get("p_value_floor")})
    if c_rows:
        c_rows.sort(key=lambda r: -(r["per_isolate"]
                                    if finite(r.get("per_isolate")) else -1))
        s.p("Two prevalences are reported for every trait, and they are two "
            "quantities rather than a biased and an unbiased one: per isolate "
            "is the clinical burden of the collection, per lineage is the "
            "share of its diversity that carries the trait. They separate "
            "whenever sampling across lineages is uneven. Beside them the "
            "effective number of carrying lineages, the reciprocal of a "
            "Herfindahl index over the carriers, is compared with the number "
            "that carriage in proportion to lineage size would give, by "
            "permutation. A trait departs from proportional carriage where "
            "that permutation p-value falls below 0.05; the direction says "
            "whether the carriers sit in fewer lineages than chance gives, in "
            "more, or in as many but in uneven shares across them.")
        s.table("Table 6", "Prevalence on two scales and the concentration of "
                "carriage, per trait. The per-lineage interval is a two-stage "
                "cluster bootstrap over lineages and then isolates; the null "
                "interval for the effective number is the permutation "
                "reference for carriage in proportion to lineage size.",
                [("Trait", "lab"), ("Per isolate", "num"),
                 ("Per lineage (95 % interval)", "num"),
                 ("Carrying lineages", "num"),
                 ("Effective number (null 95 %)", "num"),
                 ("Carriage", "lab")],
                [(r["name"], pct(r["per_isolate"]),
                  "%s (%s to %s)" % (pct(r["per_lineage"]), pct(r["pl_lo"]),
                                     pct(r["pl_hi"])),
                  "%s of %s" % (r["n_carrying"], r["n_lineages"]),
                  "%s (%s)" % (num(r["eff"], 1),
                               interval(r["null_lo"], r["null_hi"], 1)),
                  _carriage_word(r))
                 for r in c_rows])
        s.figure(fig_label(), "carriage", c_rows[:18],
                 "Left: prevalence per isolate (filled) and per lineage "
                 "(hollow, with its 95 %% interval) for each trait, %d "
                 "largest of %d by prevalence per isolate; a long connector "
                 "means that a few large lineages carry most of the "
                 "resistance, or that many small ones do. Right: the "
                 "effective number of carrying lineages (point) against the "
                 "permutation interval for carriage in proportion to lineage "
                 "size (grey band); a point below the band is carriage "
                 "concentrated in fewer lineages than chance gives, above it "
                 "carriage spread over more, and a point inside it with a "
                 "departure is carriage in as many lineages as chance gives "
                 "but in uneven shares. Traits that do not depart from "
                 "proportional carriage at level 0.05 are drawn in grey."
                 % (min(len(c_rows), 18), len(c_rows)))
    dec = sv.get("decomposition")
    if dec:
        n_ref = int(dec.get("n_within_lineage_refused") or 0)
        n_feat = int(dec.get("n_features") or 0)
        s.p("Contrast %s: of %s traits, %s differ because the two collections "
            "hold different lineages (a change in mix), and %s differ because "
            "resistance changed inside lineages (a change in rate), after "
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
    if not (md.get("prevalence_decomposition") or {}).get("per_feature"):
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
    s.callout("Terms used in this report", [
        "**Share.** How much of the difference in resistance between isolates "
        "is explained by which lineage they belong to. Near 1: resistance "
        "travels with lineages. Near 0: it moves between them.",
        "**Interval.** The range the true value is expected to lie in. It is "
        "an interval for the share the lineages in this collection carry, "
        "not for a fresh draw of lineages; the fraction of the time such "
        "intervals contain that truth was measured on the release's "
        "validation grid.",
        "**Control.** The same calculation on shuffled lineage labels, which "
        "should return roughly zero; a share is read against it.",
        "**Support.** The share of isolates sitting in lineages large enough "
        "to inform the estimate; below %s the software withholds the value."
        % pct(gate, 0),
        "**Gate.** A condition fixed before the run. If it does not hold, the "
        "software withholds the estimate and names the condition that failed.",
        "**e-value.** Evidence on a scale that stays honest when a panel is "
        "looked at again every year; larger is stronger, 1 is none.",
    ])
    S.append(s)

    symbols = [
        ("‡", "The 95 % interval includes zero. No lineage effect is "
              "distinguishable from none for this trait, and the point "
              "estimate must not be read on its own."),
        ("†", "Support below %s. Too few isolates sit in lineages large "
              "enough to inform the estimate." % pct(gate, 0)),
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
        "of the listed conditions fired. The interval for the lineages in "
        "hand is printed as computed, without clipping to the natural bounds "
        "of the quantity; clipping would leave its coverage unchanged and is "
        "left to the reader so that widths stay comparable across runs. The "
        "species interval is floored at zero, which is part of its "
        "construction.",
    ]
    return Report(ident=ident, status=status, sections=S, symbols=symbols,
                  colophon=colophon)


# ------------------------------------------------------------- the rules
def _carriage_word(r: dict) -> str:
    """The reading of one trait's carriage, from its permutation p-value and
    the direction the record gives."""
    if not r.get("departs"):
        return "proportional"
    d = str(r.get("direction") or "")
    if d == "concentrated":
        return "concentrated in fewer lineages"
    if d == "dispersed":
        return "spread over more lineages"
    return "as many lineages as chance, uneven shares"


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


def _interpretation(rows: List[dict], rejected: set) -> List[str]:
    """Three sentences fixed by the readings in Table 1. Empty sets drop out."""
    if not rows:
        return []
    admitted = [r for r in rows if r.get("estimable") is not False]
    with_evidence = [r["name"] for r in admitted
                     if not (r["lo"] <= 0.0 <= r["hi"]) and r["name"] in rejected
                     and r["pt"] >= 0.5]
    small_effect = [r["name"] for r in admitted
                    if not (r["lo"] <= 0.0 <= r["hi"]) and r["name"] in rejected
                    and r["pt"] < 0.5]
    unresolved = [r["name"] for r in admitted
                  if (r["lo"] <= 0.0 <= r["hi"]) and r["name"] in rejected]
    none_detectable = [r["name"] for r in admitted
                       if (r["lo"] <= 0.0 <= r["hi"]) and r["name"] not in rejected]
    interval_only = [r["name"] for r in admitted
                     if not (r["lo"] <= 0.0 <= r["hi"]) and r["name"] not in rejected]
    refused = [r["name"] for r in rows if r.get("estimable") is False]
    out: List[str] = []
    if with_evidence:
        out.append("For %s, resistance travels with lineage: the share is at "
                   "or above one half, so it is carried by particular clones "
                   "and will move when those clones move. Measures that act "
                   "on the clone, such as movement control and clone-directed "
                   "surveillance, act on these agents."
                   % _join(with_evidence))
    if small_effect:
        out.append("For %s, a lineage effect is detected but the share is "
                   "below one half: most of the variation in resistance sits "
                   "within lineages, so the resistance moves largely "
                   "independently of the clone, and selection pressure, dosing "
                   "and choice of agent are the levers that reach most of it."
                   % _join(small_effect))
    if unresolved:
        out.append("For %s, the e-value selects a lineage effect while the "
                   "interval for its size still reaches zero: there is "
                   "evidence that resistance is not spread evenly across the "
                   "lineages, and this collection is too small, or too "
                   "uneven across its lineages, to say how much of it the "
                   "lineages carry. Read the e-value as the finding and the "
                   "share as not yet resolved." % _join(unresolved))
    if interval_only:
        out.append("For %s, the interval for the share excludes zero while the "
                   "e-value does not select a lineage effect at the panel's "
                   "false-discovery level: the two readings disagree, and the "
                   "e-value, which carries the multiplicity correction, is the "
                   "one to read; treat the share as suggestive."
                   % _join(interval_only))
    if none_detectable:
        one = len(none_detectable) == 1
        out.append("For %s, no lineage effect is distinguishable from none in "
                   "this collection. %s resistance is not tied to any "
                   "lineage the data can see; this run gives no basis for "
                   "clone-directed action on %s, and a change in %s "
                   "prevalence is more plausibly a change in selection than "
                   "in which clones are present."
                   % (_join(none_detectable), "Its" if one else "Their",
                      "it" if one else "them", "its" if one else "their"))
    if refused:
        out.append("For %s, the estimator refused: the collection cannot "
                   "identify the quantity. No reading follows." % _join(refused))
    if out:
        out.append("These sentences follow from the table by fixed rules and "
                   "carry no judgement beyond it.")
    return out
