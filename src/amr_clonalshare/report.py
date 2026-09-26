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

import html
from typing import List, Optional

from .report_model import Block, build_report

__all__ = ["render_report"]


def _cell(text) -> str:
    return html.escape(str(text), quote=False).replace("|", "\\|").replace("\n", " ")


def _block(b: Block) -> List[str]:
    if b.kind == "p":
        return [html.escape(b.text, quote=False), ""]
    if b.kind == "kv":
        return ["- **%s:** %s" % (_cell(k), _cell(v)) for k, v in b.rows] + [""]
    if b.kind == "table":
        out = ["**%s.** %s" % (_cell(b.label), html.escape(b.caption, quote=False)), ""]
        out.append("| " + " | ".join(_cell(h) for h, _ in b.headers) + " |")
        out.append("|" + "|".join("---:" if a == "num" else "---"
                                  for _, a in b.headers) + "|")
        for row in b.rows:
            out.append("| " + " | ".join(_cell(c) for c in row) + " |")
        return out + [""]
    if b.kind == "figure":
        return ["*%s (drawn on the page).* %s" % (_cell(b.label), html.escape(b.caption, quote=False)), ""]
    if b.kind == "callout":
        return ["> **%s**" % _cell(b.title), ">"] + ["> %s" % _cell(line) + "\n>" for line in b.lines] + [""]
    return []


def render_report(result: dict, summary: dict, *, input_qc: Optional[dict] = None,
                  record_bytes: Optional[bytes] = None) -> str:
    """Markdown report from the run record, in the fixed section order."""
    rep = build_report(result, summary, input_qc=input_qc,
                       record_bytes=record_bytes)
    lines = ["# amr-clonalshare run report", ""]
    lines.append("Every number below is read from `clonal_share_result.json`, "
                 "written by the same run, or from the release's validation grid "
                 "where the text says so; none is recomputed here.")
    lines.append("")
    ident = {key: _cell(value) for key, value in rep.ident.items()}
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
        lines += ["## %d. %s" % (n, _cell(sec.title)), ""]
        for b in sec.blocks:
            lines += _block(b)
    lines += ["## Notes on reading this report", ""]
    for sym, text in rep.symbols:
        lines.append("- %s %s" % (_cell(sym), html.escape(text, quote=False)))
    lines.append("")
    lines += [html.escape(line, quote=False) for line in rep.colophon]
    return "\n".join(lines).rstrip() + "\n"
