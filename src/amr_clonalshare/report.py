"""The run report in prose: the same structure as the page, rendered as text.

``clonal_share_result.json`` is the record. :mod:`.report_model` reads it once into
sections of typed blocks; this module writes those blocks as Markdown, and
:mod:`.report_html` writes the same blocks as a page with the figures drawn.
The two forms therefore cannot disagree about what a run found, and this one
diffs cleanly between runs.

Figures are not drawn here. Where the page carries a figure, this form carries
its caption, so a reader of the text knows what the page would have shown.
"""
from __future__ import annotations

from typing import List, Optional

from .report_model import Block, build_report, num, pct

__all__ = ["render_report"]


def _f(x, nd=3) -> str:
    """Kept for callers that format numbers the way the report does."""
    if isinstance(x, str):
        return x
    return num(x, nd)


def _pct(x) -> str:
    return pct(x)


def _cell(text) -> str:
    return str(text).replace("|", "\\|").replace("\n", " ")


def _block(b: Block) -> List[str]:
    if b.kind == "p":
        return [b.text, ""]
    if b.kind == "kv":
        return ["- **%s:** %s" % (k, v) for k, v in b.rows] + [""]
    if b.kind == "table":
        out = ["**%s.** %s" % (b.label, b.caption), ""]
        out.append("| " + " | ".join(_cell(h) for h, _ in b.headers) + " |")
        out.append("|" + "|".join("---:" if a == "num" else "---"
                                  for _, a in b.headers) + "|")
        for row in b.rows:
            out.append("| " + " | ".join(_cell(c) for c in row) + " |")
        return out + [""]
    if b.kind == "figure":
        return ["*%s (drawn on the page).* %s" % (b.label, b.caption), ""]
    if b.kind == "callout":
        return ["> **%s**" % b.title, ">"] + ["> %s" % line + "\n>" for line in b.lines] + [""]
    return []


def render_report(result: dict, summary: dict, *, input_qc: Optional[dict] = None,
                  record_bytes: Optional[bytes] = None) -> str:
    """Markdown report from the run record, in the fixed section order."""
    rep = build_report(result, summary, input_qc=input_qc,
                       record_bytes=record_bytes)
    lines = ["# amr-clonalshare run report", ""]
    lines.append("Every number below is read from `clonal_share_result.json`, "
                 "written by the same run; none is recomputed here.")
    lines.append("")
    ident = rep.ident
    lines += ["- **Run:** %s" % ident["run"],
              "- **Issued:** %s" % ident["issued"],
              "- **Software:** %s" % ident["software"],
              "- **Isolates:** %s" % ident["isolates"],
              "- **Lineage:** %s" % ident["lineage"],
              "- **Seed:** %s" % ident["seed"],
              "- **Configuration:** %s" % ident["configuration"]]
    if ident.get("record_sha256"):
        lines.append("- **Record digest:** sha256:%s" % ident["record_sha256"])
    st = rep.status
    lines += ["", "**%s.** %s" % (st["lead"][0].upper() + st["lead"][1:],
                                   st["gates_line"]), ""]
    n = 0
    for sec in rep.sections:
        n += 1
        lines += ["## %d. %s" % (n, sec.title), ""]
        for b in sec.blocks:
            lines += _block(b)
    lines += ["## Notes on reading this report", ""]
    for sym, text in rep.symbols:
        lines.append("- %s %s" % (sym, text))
    lines.append("")
    lines += rep.colophon
    return "\n".join(lines).rstrip() + "\n"
