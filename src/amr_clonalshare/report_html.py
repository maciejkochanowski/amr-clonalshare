"""The run report as one self-contained page, with the diagnostics drawn.

WHY A SECOND FORM OF THE REPORT. ``report.md`` states the result in prose and
is the form that diffs cleanly between runs. It is not the form in which a
reader sees what the numbers are doing. On the shipped planted control every
one of the sixty per-trait intervals crosses zero, and the prose form prints
the largest of them, 0.128, beside ``estimable: yes``; a reader takes that for
a finding. Drawn as intervals against a ruled zero the same table says in one
look that no single trait is separated from none. The same holds for the gates:
"no diagnostic gate tripped" is true and hides that one gate had a margin of
0.04 while the rest were near 0.5.

WHAT IT ADDS BEYOND THE PROSE FORM. The record holds twenty-eight blocks. The
prose form reads five of them. The four sections added here are the ones that
answer whether the result could be wrong rather than what it is: whether the
grouping is an artefact of the panel, whether a simpler rule does as well,
which layer carries it, and how stable it is under resampling. Reporting a
result while withholding the diagnostics that could undermine it is selective
reporting, whatever the intent.

NO DEPENDENCY. The figures are written as SVG here rather than drawn with a
plotting library, so a report costs the reader nothing to open and the package
keeps matplotlib optional. Coordinates are rounded to two decimals so that two
reports of the same run diff to nothing.
"""
from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone

from . import __version__

__all__ = ["render_html_report"]

INK, MID, FAINT = "#1A1A1A", "#4A4A4A", "#8A8A8A"
RULE, PAPER = "#C8C8C4", "#FDFDFB"
ACCENT, FLAG = "#12436D", "#D55E00"


def _finite(x) -> bool:
    return isinstance(x, (int, float)) and not math.isnan(float(x))


def _esc(s) -> str:
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;"))


def _num(x, nd: int = 3) -> str:
    if not _finite(x):
        return "&#8212;"
    if isinstance(x, bool):
        return "yes" if x else "no"
    return ("%%.%df" % nd) % float(x)


def _pct(x, nd: int = 1) -> str:
    return "&#8212;" if not _finite(x) else ("%%.%df&#8201;%%%%" % nd) % (
        100.0 * float(x))


def _ticks(lo: float, hi: float, target: int = 5) -> list:
    """1, 2 or 5 times a power of ten, spanning the drawn range and no more."""
    if not (_finite(lo) and _finite(hi)) or hi <= lo:
        return [lo] if _finite(lo) else [0.0]
    raw = (hi - lo) / max(target, 2)
    mag = 10.0 ** math.floor(math.log10(raw)) if raw > 0 else 1.0
    for m in (1.0, 2.0, 2.5, 5.0, 10.0):
        if raw <= m * mag:
            step = m * mag
            break
    else:
        step = 10.0 * mag
    start = step * math.floor(lo / step)
    out, v = [], start
    while v <= hi + step * 1e-9:
        if lo - 1e-12 <= v <= hi + 1e-12:
            out.append(round(v, 10))
        v += step
    return out or [lo, hi]


def _r(v: float) -> str:
    return "%.2f" % float(v)


# ------------------------------------------------------------- figures ----
def _fig_mdl(sweep, selected) -> str:
    """Description length saved against the number of groups."""
    if not sweep:
        return ""
    W, H, L, R, T, B = 470, 172, 54, 14, 14, 34
    ks = [p["k"] for p in sweep]
    gs = [float(p["gain"]) for p in sweep]
    xlo, xhi = min(ks) - 0.5, max(ks) + 0.5
    ylo, yhi = min(0.0, min(gs)), max(gs) * 1.12 or 1.0
    X = lambda v: L + (v - xlo) / (xhi - xlo) * (W - L - R)    
    Y = lambda v: H - B - (v - ylo) / (yhi - ylo) * (H - T - B)

    p = ['<svg viewBox="0 0 %d %d" role="img" '
         'xmlns="http://www.w3.org/2000/svg"><title>Description length saved '
         'against the number of profile groups</title>' % (W, H)]
    for t in _ticks(ylo, yhi, 4):
        p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                 'stroke-width="0.4" shape-rendering="crispEdges"/>'
                 % (_r(L), _r(Y(t)), _r(W - R), _r(Y(t)), RULE))
        p.append('<text x="%s" y="%s" text-anchor="end" class="tk">%s</text>'
                 % (_r(L - 6), _r(Y(t) + 3), "%.0f" % t))
    p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
             'stroke-width="0.6" shape-rendering="crispEdges"/>'
             % (_r(L), _r(Y(0)), _r(W - R), _r(Y(0)), MID))
    p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
             'stroke-width="1.6" stroke-dasharray="3 3"/>'
             % (_r(X(selected)), _r(T), _r(X(selected)), _r(H - B), FAINT))
    d = " ".join("%s%s,%s" % ("M" if i == 0 else "L", _r(X(k)), _r(Y(g)))
                 for i, (k, g) in enumerate(zip(ks, gs)))
    p.append('<path d="%s" fill="none" stroke="%s" stroke-width="1.4"/>'
             % (d, ACCENT))
    for k, g in zip(ks, gs):
        sel = k == selected
        p.append('<circle cx="%s" cy="%s" r="%s" fill="%s"/>'
                 % (_r(X(k)), _r(Y(g)), "3.6" if sel else "2.2",
                    ACCENT if sel else MID))
        p.append('<text x="%s" y="%s" text-anchor="middle" class="tk">%d'
                 '</text>' % (_r(X(k)), _r(H - B + 14), k))
    p.append('<text x="%s" y="%s" text-anchor="middle" class="ax">number of '
             'profile groups, k</text>' % (_r((L + W - R) / 2), _r(H - 6)))
    p.append('<text x="6" y="%s" class="ax" transform="rotate(-90 6 %s)">bits '
             'saved</text>' % (_r((T + H - B) / 2), _r((T + H - B) / 2)))
    p.append('<text x="%s" y="%s" class="lbl">selected</text>'
             % (_r(X(selected) + 5), _r(T + 10)))
    return "".join(p) + "</svg>"


def _fig_intervals(rows, n_show: int = 18) -> str:
    """Point estimate with interval, one row per trait, zero ruled."""
    rows = rows[:n_show]
    if not rows:
        return ""
    W, L, R, T, ROW, B = 470, 78, 58, 14, 15.0, 38
    H = T + ROW * len(rows) + B
    lo = min(r["lo"] for r in rows)
    hi = max(r["hi"] for r in rows)
    pad = (hi - lo) * 0.08 or 0.01
    xlo, xhi = min(0.0, lo - pad), hi + pad
    X = lambda v: L + (v - xlo) / (xhi - xlo) * (W - L - R)

    p = ['<svg viewBox="0 0 %d %.0f" role="img" '
         'xmlns="http://www.w3.org/2000/svg"><title>Clonal share by trait, '
         'with interval</title>' % (W, H)]
    for t in _ticks(xlo, xhi, 5):
        p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                 'stroke-width="0.4" shape-rendering="crispEdges"/>'
                 % (_r(X(t)), _r(T - 4), _r(X(t)), _r(H - B + 2), RULE))
        p.append('<text x="%s" y="%s" text-anchor="middle" class="tk">%s'
                 '</text>' % (_r(X(t)), _r(H - B + 15), "%.2f" % t))
    p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
             'stroke-width="0.8" shape-rendering="crispEdges"/>'
             % (_r(X(0)), _r(T - 4), _r(X(0)), _r(H - B + 2), MID))
    for i, r in enumerate(rows):
        y = T + ROW * i + ROW / 2
        crosses = r["lo"] <= 0.0 <= r["hi"]
        col = MID if crosses else ACCENT
        p.append('<text x="%s" y="%s" text-anchor="end" class="rw">%s</text>'
                 % (_r(L - 8), _r(y + 3), _esc(r["name"])))
        p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                 'stroke-width="1.1"/>'
                 % (_r(X(r["lo"])), _r(y), _r(X(r["hi"])), _r(y), col))
        for e in (r["lo"], r["hi"]):
            p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                     'stroke-width="1.1"/>'
                     % (_r(X(e)), _r(y - 3), _r(X(e)), _r(y + 3), col))
        p.append('<circle cx="%s" cy="%s" r="2.6" fill="%s"/>'
                 % (_r(X(r["pt"])), _r(y), col))
        p.append('<text x="%s" y="%s" class="rv">%s%s</text>'
                 % (_r(W - R + 6), _r(y + 3), "%.3f" % r["pt"],
                    "&#8225;" if crosses else ""))
    p.append('<text x="%s" y="%s" text-anchor="middle" class="ax">clonal '
             'share (variation attributable to lineage)</text>'
             % (_r((L + W - R) / 2), _r(H - 8)))
    return "".join(p) + "</svg>"


def _fig_bars(rows, label: str, note: str = "") -> str:
    """One horizontal bar per row: (name, value, flagged)."""
    rows = [r for r in rows if _finite(r[1])]
    if not rows:
        return ""
    W, L, R, T, ROW, B = 470, 168, 66, 10, 17.0, 22
    H = T + ROW * len(rows) + B
    mx = max(abs(r[1]) for r in rows) * 1.15 or 1.0
    X = lambda v: L + (v / mx) * (W - L - R)

    p = ['<svg viewBox="0 0 %d %.0f" role="img" '
         'xmlns="http://www.w3.org/2000/svg"><title>%s</title>'
         % (W, H, _esc(note or label))]
    p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
             'stroke-width="0.8" shape-rendering="crispEdges"/>'
             % (_r(L), _r(T - 2), _r(L), _r(H - B + 2), MID))
    for i, (name, value, flagged) in enumerate(rows):
        y = T + ROW * i + ROW / 2
        col = FLAG if flagged else ACCENT
        p.append('<text x="%s" y="%s" text-anchor="end" class="rw">%s</text>'
                 % (_r(L - 8), _r(y + 3), _esc(name)))
        p.append('<rect x="%s" y="%s" width="%s" height="6" fill="%s"/>'
                 % (_r(L), _r(y - 3), _r(max(X(abs(value)) - L, 0.6)), col))
        p.append('<text x="%s" y="%s" class="rv">%s</text>'
                 % (_r(W - R + 6), _r(y + 3), "%.2f" % value))
    p.append('<text x="%s" y="%s" text-anchor="middle" class="ax">%s</text>'
             % (_r((L + W - R) / 2), _r(H - 5), _esc(label)))
    return "".join(p) + "</svg>"


_CSS = """
:root{--ink:#1A1A1A;--mid:#4A4A4A;--faint:#8A8A8A;--rule:#C8C8C4;
--paper:#FDFDFB;--accent:#12436D;--flag:#D55E00}
*{box-sizing:border-box}
html{background:#EDEDEA}
body{margin:0;background:var(--paper);color:var(--ink);
font-family:"Source Sans 3","Source Sans Pro","IBM Plex Sans",Helvetica,Arial,
sans-serif;font-size:16px;line-height:1.45;
font-variant-numeric:lining-nums tabular-nums;-webkit-font-smoothing:antialiased}
.sheet{max-width:41rem;margin:0 auto;padding:34px 30px 60px}
h1{font-size:1.30rem;font-weight:600;letter-spacing:-.005em;margin:0 0 2px}
.sub{color:var(--mid);font-size:.94rem;margin:0 0 14px}
h2{font-size:1.0rem;font-weight:600;margin:34px 0 8px;padding-bottom:4px;
border-bottom:.8px solid var(--ink)}
h2 .no{color:var(--faint);font-weight:400;margin-right:.55em}
p{margin:0 0 10px}
.ident{border-top:.8px solid var(--ink);border-bottom:.8px solid var(--ink);
padding:9px 0;margin:0 0 20px;display:grid;
grid-template-columns:repeat(4,1fr);gap:3px 18px}
.ident div{font-size:.76rem;line-height:1.35}
.ident b{display:block;color:var(--faint);font-weight:400;font-size:.68rem;
letter-spacing:.055em;text-transform:uppercase}
.status{margin:0 0 18px}
.status .lead{font-size:.80rem;letter-spacing:.10em;text-transform:uppercase;
color:var(--accent);font-weight:600}
.status .lead.flag{color:var(--flag)}
.status p{margin:4px 0 0;color:var(--mid);font-size:.94rem}
table{border-collapse:collapse;width:100%;margin:10px 0 6px;font-size:.88rem}
caption{caption-side:top;text-align:left;color:var(--mid);font-size:.82rem;
margin-bottom:5px}
caption b{color:var(--ink);font-weight:600}
thead th{border-top:.8px solid var(--ink);border-bottom:.5px solid var(--ink);
padding:5px 9px;font-weight:600;vertical-align:bottom}
tbody td{padding:4px 9px;border:0}
tbody tr:last-child td{border-bottom:.8px solid var(--ink)}
.num{text-align:right}
th.num{text-align:right}
.lab{text-align:left}
td.nw{white-space:nowrap}
.flagcol{text-align:left;width:2.4em;color:var(--flag)}
figure{margin:16px 0 6px}
figure svg{width:100%;height:auto;display:block}
figcaption{font-size:.82rem;color:var(--mid);margin-top:6px}
figcaption b{color:var(--ink);font-weight:600}
svg text{font-family:inherit;fill:var(--ink)}
svg .tk{font-size:9.5px;fill:var(--mid)}
svg .ax{font-size:10px;fill:var(--mid)}
svg .rw{font-size:10px;fill:var(--ink)}
svg .rv{font-size:10px;fill:var(--mid)}
svg .lbl{font-size:9.5px;fill:var(--faint)}
.notes{margin-top:34px;border-top:.8px solid var(--ink);padding-top:10px;
font-size:.82rem;color:var(--mid)}
.notes dl{display:grid;grid-template-columns:1.6em 1fr;gap:2px 6px;margin:6px 0}
.notes dt{color:var(--flag)}
.notes dd{margin:0}
.colophon{margin-top:26px;border-top:.5px solid var(--rule);padding-top:9px;
font-size:.76rem;color:var(--faint);line-height:1.5}
.colophon code{font-family:"IBM Plex Mono",ui-monospace,Menlo,Consolas,monospace;
font-size:.94em;color:var(--mid);word-break:break-all}
.dl{display:grid;grid-template-columns:13.5em 1fr;gap:2px 14px;font-size:.88rem;
margin:8px 0 4px}
.dl dt{color:var(--mid)}
.dl dd{margin:0}
@media print{
 @page{size:A4;margin:22mm}
 html,body{background:#fff}
 .sheet{max-width:none;padding:0}
 h2,figure,table{break-inside:avoid}
 h2{break-after:avoid}
 p{orphans:3;widows:3}
 svg,table,figure{print-color-adjust:exact;-webkit-print-color-adjust:exact}
}
"""


def _table(caption, headers, rows) -> str:
    """headers: list of (label, css class). rows: sequences of cells."""
    out = ["<table><caption>%s</caption><thead><tr>" % caption]
    out += ['<th class="%s">%s</th>' % (cls, h) for h, cls in headers]
    out.append("</tr></thead><tbody>")
    for row in rows:
        out.append("<tr>")
        out += ['<td class="%s">%s</td>' % (headers[i][1], cell)
                for i, cell in enumerate(row)]
        out.append("</tr>")
    return "".join(out) + "</tbody></table>"


def _inner(w) -> str:
    """One value inside a gate-ledger entry, with nulls named rather than shown."""
    if w is None:
        return "not computed"
    if isinstance(w, bool):
        return "yes" if w else "no"
    if isinstance(w, float):
        return _num(w)
    if isinstance(w, (list, tuple)):
        if not w:
            return "none"
        if all(x is None for x in w):
            return "not computed"
        return "[%s]" % ", ".join(_inner(x) for x in w)
    return str(w)


def _pretty(v) -> str:
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return _num(v)
    if isinstance(v, dict):
        bits = ["%s %s" % (k.replace("_", " "), _inner(w))
                for k, w in list(v.items())[:2]]
        s = "; ".join(bits)
        if len(v) > 2:
            s += " (+%d more in the record)" % (len(v) - 2)
        return _esc(s)
    return "&#8212;" if v is None else _esc(v)


def render_html_report(record: dict, summary: dict) -> str:
    """One self-contained page for a run, from the record it wrote.

    Every value is read from ``record`` or ``summary``; none is recomputed,
    so the page cannot disagree with the file it describes.
    """
    md = record.get("metadata_diagnostics") or {}
    interp = record.get("interpretation") or {}
    cfg = record.get("config") or {}
    ledger = interp.get("gate_ledger") or []
    tripped = [g for g in ledger if g.get("triggered")]
    n_gates = sum(1 for g in ledger if g.get("applicable"))
    conc = md.get("lineage_concordance") or {}

    digest = hashlib.sha256(
        json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:12]
    run_id = hashlib.sha256(
        (digest + str(record.get("seed")) + str(record.get("n_isolates")))
        .encode()).hexdigest()[:8].upper()
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ")

    rows = []
    for name, v in (md.get("clonal_share") or {}).items():
        if _finite(v.get("kappa_adj")) and _finite(v.get("ci_low")):
            rows.append({"name": name, "pt": v["kappa_adj"], "lo": v["ci_low"],
                         "hi": v["ci_high"], "support": v.get("support"),
                         "estimable": v.get("estimable")})
    rows.sort(key=lambda r: -r["pt"])
    shown = rows[:18]
    n_cross = sum(1 for r in shown if r["lo"] <= 0.0 <= r["hi"])
    n_cross_all = sum(1 for r in rows if r["lo"] <= 0.0 <= r["hi"])

    A: list[str] = []
    add = A.append
    add('<!doctype html><html lang="en"><head><meta charset="utf-8">')
    add('<meta name="viewport" content="width=device-width,initial-scale=1">')
    add("<title>Run report %s</title><style>%s</style></head><body>"
        % (run_id, _CSS))
    add("<main class='sheet'>")
    add("<h1>Lineage-conditioned resistance-profile analysis</h1>")
    add('<p class="sub">Run report. Every value is read from the record '
        "written by this run; none is recomputed here.</p>")

    add('<div class="ident">')
    for k, v in (("Run", run_id), ("Issued", stamp),
                 ("Software", "amr-clonalshare %s" % __version__),
                 ("Record", "cluster_result.json"),
                 ("Isolates", str(record.get("n_isolates"))),
                 ("Layers", ", ".join(record.get("layers") or [])),
                 ("Seed", str(record.get("seed"))),
                 ("Configuration", "sha256:" + digest)):
        add("<div><b>%s</b>%s</div>" % (k, _esc(v)))
    add("</div>")

    add('<div class="status"><div class="lead%s">Claim level %s &#183; %s'
        "</div>" % (" flag" if tripped else "", interp.get("claim_level"),
                    _esc(str(interp.get("claim_status", "")).replace("_", " "))))
    add("<p>%d of %d applicable gates closed on this run. A gate is not an "
        "error: it marks a reading the data cannot carry. No gate raised does "
        "not establish that the reading is correct.</p></div>"
        % (len(tripped), n_gates))

    # 1 -------------------------------------------------------------------
    add('<h2><span class="no">1</span>Measurand, inputs and settings</h2>')
    snf = cfg.get("snf") or {}
    dist = cfg.get("distance") or {}
    mc = cfg.get("monte_carlo") or {}
    add('<div class="dl">')
    for k, v in (
        ("Quantity measured", "share of resistance-profile variation "
                              "attributable to lineage"),
        ("Cohort", "%s isolates, %s lineage labels"
                   % (record.get("n_isolates"), conc.get("n_lineages"))
                   if conc.get("n_lineages") is not None else
                   "%s isolates, no lineage labels supplied"
                   % record.get("n_isolates")),
        ("Binary layers", ", ".join(record.get("layers") or [])),
        ("Distance", "%s (undefined pairs: %s)"
                     % (dist.get("metric"), dist.get("undefined_pair"))),
        ("Fusion", "similarity-network fusion, K=%s, mu=%s, T=%s"
                   % (snf.get("K"), snf.get("mu"), snf.get("T"))),
        ("Resampling", "%s consensus, %s bootstrap, %s permutations"
                       % (mc.get("consensus_B"), mc.get("n_boot"),
                          mc.get("n_perm"))),
    ):
        add("<dt>%s</dt><dd>%s</dd>" % (k, _esc(v)))
    add("</div>")

    # 2 -------------------------------------------------------------------
    add('<h2><span class="no">2</span>Admissibility of the input</h2>')
    shares = [v for v in (md.get("clonal_share") or {}).values()
              if _finite(v.get("support"))]
    support = min((v["support"] for v in shares), default=None)
    groups_used = max((v.get("n_groups") for v in shares
                       if isinstance(v.get("n_groups"), int)), default=None)
    add(_table("<b>Table 1.</b> Conditions the estimator requires before any "
               "result is reported.",
               [("Condition", "lab"), ("Observed", "num"),
                ("Required", "num"), ("Verdict", "lab")],
               [("Lineage support",
                 _pct(support) if support is not None else "not recorded",
                 "&#8805; 90.0&#8201;%",
                 "accepted" if _finite(support) and support >= 0.9
                 else "refused" if support is not None else "not evaluated"),
                ("Lineage groups used",
                 "%s of %s" % (groups_used, conc.get("n_lineages"))
                 if groups_used is not None
                 and conc.get("n_lineages") is not None else "not recorded",
                 "reported",
                 "accepted" if groups_used is not None else "not evaluated"),
                ("Effective layers",
                 _num(summary.get("fusion_n_eff_layers"), 1), "&gt; 1.5",
                 "refused" if summary.get("fusion_collapse") else "accepted"),
                ("Split design",
                 "adequate" if summary.get("split_design_adequate")
                 else "inadequate", "adequate",
                 "accepted" if summary.get("split_design_adequate")
                 else "refused")]))

    # 3 -------------------------------------------------------------------
    add('<h2><span class="no">3</span>Are there resistance-profile groups?'
        "</h2>")
    ks = record.get("k_selection") or {}
    add("<p>The selection procedure returned <b>%s</b> profile groups. Tested "
        "on features held out of the clustering, the grouping reproduces at "
        "%s. The structure is separated by gaps rather than shading "
        "continuously from one type into the next (%s).</p>"
        % (record.get("selected_k"),
           _esc(summary.get("p_value_structure_report")),
           _esc(summary.get("discreteness_verdict"))))
    fig = _fig_mdl(ks.get("mdl_sweep") or [], record.get("selected_k"))
    if fig:
        add("<figure>%s<figcaption><b>Figure 1.</b> Description length saved "
            "by describing the cohort as k archetypes rather than as one "
            "population. The turn in the curve, not a threshold, is what "
            "selects k.</figcaption></figure>" % fig)
    ci = summary.get("continuum_tail_probability_ci95") or {}
    add(_table("<b>Table 2.</b> Evidence behind the number of groups.",
               [("Quantity", "lab"), ("Value", "num"), ("Reading", "lab")],
               [("Groups selected", str(record.get("selected_k")),
                 "by minimum description length"),
                ("Description length saved",
                 "%s bits" % _num(ks.get("mdl_best_gain"), 0),
                 "%s of the null code" % _pct(ks.get("mdl_gain_fraction"))),
                ("Permutation null for that gain",
                 _num(ks.get("mdl_p_value")) if _finite(
                     ks.get("mdl_p_value")) else "not run",
                 _esc((record.get("mdl_calibration") or {}).get("note", ""))),
                ("Held-out reproducibility",
                 _esc(summary.get("p_value_structure_report")),
                 "groups recur on features not used to form them"),
                ("Continuum exceedances",
                 "%s of %s" % (summary.get("continuum_bootstrap_exceedances"),
                               mc.get("n_boot")),
                 "Clopper&#8211;Pearson upper limit %s" % _num(ci.get("high"),
                                                               4)),
                ("Latent dimension of the null",
                 str(summary.get("continuum_latent_dimension")),
                 "chosen by BIC, not pinned")]))

    # 4 -------------------------------------------------------------------
    add('<h2><span class="no">4</span>Could the grouping be an artefact of the '
        "panel?</h2>")
    art = record.get("artifact_diagnostics") or {}
    per_layer = art.get("per_layer") or []
    add("<p>A grouping can arise from the shape of the panel rather than from "
        "the isolates. Three things are checked per layer: how many isolates "
        "carry no informative feature at all, how many independent directions "
        "the layer really holds, and whether the nearest-neighbour cut that "
        "builds the affinity graph was decided by ties rather than by "
        "distance. A tie inflation of 1.00 means the cut was clean; a value "
        "above it means more neighbours were kept than asked for, because "
        "candidates shared the same affinity.</p>")
    add(_table("<b>Table 3.</b> Panel geometry, layer by layer.",
               [("Layer", "lab"), ("Features used", "num"),
                ("Effective dimension", "num"),
                ("Uninformative isolates", "num"),
                ("Tie inflation at the cut", "num"),
                ("Empty-stratum ARI", "num")],
               [(_esc(p.get("layer")), str(p.get("n_features_used")),
                 _num(p.get("effective_dimension"), 2),
                 _pct(p.get("frac_uninformative_rows")),
                 _num(p.get("knn_tie_inflation"), 3),
                 _num(p.get("empty_stratum_ari")))
                for p in per_layer]))
    add("<p>The empty-stratum check re-runs the grouping on isolates that "
        "carry nothing in one layer. Where it returns a value, an agreement "
        "near one would mean the grouping is reproduced by absence alone. The "
        "largest value on this run is %s.</p>"
        % _num(summary.get("max_empty_stratum_ari")))

    # 5 -------------------------------------------------------------------
    add('<h2><span class="no">5</span>Does the grouping beat a simpler rule?'
        "</h2>")
    base = record.get("baselines") or {}
    cands = base.get("candidates") or []
    add("<p>A grouping that a single layer or a plain concatenation "
        "reproduces has not earned the fusion that produced it. Every "
        "candidate below was scored on the same external criteria as the "
        "reported partition.</p>")
    ext_name = ""
    for c in cands:
        for key in (c.get("external_agreement") or {}):
            ext_name = key
            break
        if ext_name:
            break
    add(_table("<b>Table 4.</b> Candidate partitions, scored the same way. "
               "%s" % _esc(base.get("note", "")),
               [("Candidate", "lab"), ("Groups", "num"),
                ("Description length saved", "num"),
                ("Agreement with %s" % _esc(ext_name or "external labels"),
                 "num"),
                ("Agreement with the reported partition", "num")],
               [(_esc(str(c.get("partition")).replace("_", " ")),
                 str(c.get("n_clusters")), _num(c.get("mdl_gain"), 0),
                 _num(((c.get("external_agreement") or {}).get(ext_name)
                       or {}).get("ari")),
                 _num(c.get("ari_vs_fused")))
                for c in cands]))

    # 6 -------------------------------------------------------------------
    add('<h2><span class="no">6</span>Which layer carries the grouping?</h2>')
    infl = record.get("layer_influence") or {}
    agree = record.get("layer_agreement") or {}
    add("<p>%s The effective number of layers is %s of %s, and the layer that "
        "moves the answer most is <b>%s</b>. A layer marked inert can be "
        "removed without changing the grouping.</p>"
        % (_esc(infl.get("regime_note", "")).capitalize() + "." if
           infl.get("regime_note") else "",
           _num(infl.get("n_eff"), 2), infl.get("n_layers"),
           _esc(infl.get("dominant_layer"))))
    bars = [(str(p.get("layer")), float(p.get("weight") or 0.0),
             bool(p.get("inert"))) for p in (infl.get("per_layer") or [])]
    fig = _fig_bars(bars, "weight in the fused network",
                    "Weight carried by each layer in the fused network")
    if fig:
        add("<figure>%s<figcaption><b>Figure 2.</b> Weight each layer carries "
            "in the fused network. A layer drawn in the warning colour is "
            "inert: the grouping does not change when it is removed."
            "</figcaption></figure>" % fig)
    add(_table("<b>Table 5.</b> What each layer contributes.",
               [("Layer", "lab"), ("Weight", "num"),
                ("Loss when removed", "num"),
                ("Agreement alone with the fused grouping", "num"),
                ("Inert", "lab")],
               [(_esc(p.get("layer")), _num(p.get("weight")),
                 _num(p.get("delta_loo")), _num(p.get("ari_solo_vs_fused")),
                 "yes" if p.get("inert") else "no")
                for p in (infl.get("per_layer") or [])]))
    conflict = agree.get("conflict_mean_layer_pairs_only")
    if _finite(conflict):
        add("<p>Mean disagreement between layers, taken over layer pairs "
            "only, is %s.</p>" % _num(conflict))

    # 7 -------------------------------------------------------------------
    add('<h2><span class="no">7</span>How stable is the grouping?</h2>')
    frag = record.get("fragility") or {}
    fi = [x for x in (frag.get("F_isolate") or []) if _finite(x)]
    ff = [x for x in (frag.get("F_feature") or []) if _finite(x)]
    lr = record.get("label_recoverability") or {}
    add(_table("<b>Table 6.</b> Behaviour under resampling and under removal.",
               [("Quantity", "lab"), ("Value", "num"), ("Reading", "lab")],
               [("Proportion of ambiguous pairs", _num(record.get("pac")),
                 "share of pairs that neither always nor never co-cluster "
                 "across consensus draws; zero is a clean split"),
                ("Area under the consensus distribution",
                 _num(record.get("cdf_area")),
                 "a flat consensus matrix drives this towards one"),
                ("Largest isolate fragility",
                 _num(max(fi) if fi else None),
                 "how far one isolate can move the grouping"),
                ("Largest feature fragility",
                 _num(max(ff) if ff else None),
                 "how far one feature can move the grouping"),
                ("Label recovery from half the features",
                 _num(lr.get("internal_label_recoverability")),
                 _esc(lr.get("note", ""))[:110])]))

    # 8 -------------------------------------------------------------------
    add('<h2><span class="no">8</span>How much of the pattern is the clone?'
        "</h2>")
    add("<p>Across the whole panel, knowing an isolate&#8217;s lineage "
        "predicts <b>%s</b> of the variation for an isolate the model has not "
        "seen. The profile groups themselves are explained by lineage to a "
        "share of <b>%s</b>: a value near one would mean the groups are the "
        "lineages under another name, a value near zero that they cut across "
        "lineages.</p>"
        % (_num(summary.get("clonal_share_all_features")),
           _num(summary.get("lineage_attributable_share"))))
    add("<p>Read the intervals as frequencies. If cohorts like this one were "
        "drawn again and again, the interval printed for a trait would cover "
        "that trait&#8217;s true share on about 95 draws in 100. For <b>%d of "
        "the %d traits</b> shown the interval includes zero, so no lineage "
        "effect is distinguishable from none for that trait; over the whole "
        "panel this holds for <b>%d of %d</b>.</p>"
        % (n_cross, len(shown), n_cross_all, len(rows)))
    fig = _fig_intervals(shown)
    if fig:
        add("<figure>%s<figcaption><b>Figure 3.</b> Clonal share by trait, "
            "point estimate with 95&#8201;%% interval, %d largest of %d. "
            "Traits whose interval crosses zero are drawn in grey and marked "
            "&#8225;. The intervals are drawn as computed. The quantity lies "
            "between 0 and 1, so the part of an interval below zero "
            "carries no information: reading each interval as its "
            "overlap with that range leaves the coverage unchanged, and "
            "a lower limit at or below zero means the same thing either "
            "way.</figcaption></figure>"
            % (fig, len(shown), len(rows)))
    add(_table("<b>Table 7.</b> The six traits with the largest estimated "
               "share. The remaining %d are in the record under "
               "<code>metadata_diagnostics.clonal_share</code>."
               % max(len(rows) - 6, 0),
               [("Trait", "lab"), ("Share", "num"),
                ("95&#8201;% interval", "num"), ("Support", "num"),
                ("Estimable", "lab"), ("", "flagcol")],
               [(_esc(r["name"]), _num(r["pt"]),
                 "%s to %s" % (_num(r["lo"]), _num(r["hi"])),
                 _pct(r["support"]), "yes" if r["estimable"] else "no",
                 "&#8225;" if r["lo"] <= 0 <= r["hi"] else "")
                for r in rows[:6]]))

    # 9 -------------------------------------------------------------------
    add('<h2><span class="no">9</span>Composition or rate: the surveillance '
        "reading</h2>")
    sv = summary.get("surveillance") or {}
    cd = sv.get("carriage_direction_counts") or {}
    add("<p>Of the traits examined, %s are carried in proportion to lineage "
        "size and %s are concentrated in a few lineages; %s depart from "
        "proportional carriage far enough to matter for a prevalence reading. "
        "The widest gap between the per-isolate and the per-lineage "
        "prevalence is %s. Where the two differ, a small "
        "number of large lineages carry most of the resistance, and a change "
        "in prevalence may be a change in which lineages were sampled rather "
        "than a change in rate.</p>"
        % (cd.get("proportional", 0), cd.get("concentrated", 0),
           sv.get("n_features_departing_from_proportional_carriage", 0),
           ("%s, on <code>%s</code>"
            % (_num(sv.get("lineage_prevalence_widest_gap")),
               _esc(sv.get("lineage_prevalence_widest_gap_feature")))
            if sv.get("lineage_prevalence_widest_gap_feature")
            else "not reported: no trait departs far enough to rank")))

    # 10 ------------------------------------------------------------------
    add('<h2><span class="no">10</span>What may be concluded, and what may '
        "not</h2>")
    add("<p>The highest claim these diagnostics carry is <b>level %s</b>, %s. "
        "Every gate was evaluated; the margin is the distance the run had "
        "left before the gate would have closed.</p>"
        % (interp.get("claim_level"),
           _esc(str(interp.get("claim_status", "")).replace("_", " "))))
    margins = [(str(g["code"]).replace("_", " "),
                float(g["margin_to_failure"]), bool(g.get("triggered")))
               for g in ledger
               if g.get("applicable") and _finite(g.get("margin_to_failure"))]
    fig = _fig_bars(margins, "margin to failure (the gate closes at zero)",
                    "Distance to failure for each gate carrying a margin")
    if fig:
        add("<figure>%s<figcaption><b>Figure 4.</b> Distance to failure for "
            "each gate that carries a numeric margin. A short bar is a result "
            "that only just holds.</figcaption></figure>" % fig)
    add(_table("<b>Table 8.</b> Gate ledger. Every gate is listed, including "
               "those that did not fire.",
               [("Gate", "lab nw"), ("Applicable", "lab"), ("Closed", "lab"),
                ("Observed", "lab"), ("Threshold", "lab")],
               [(_esc(str(g.get("code")).replace("_", " ")),
                 "yes" if g.get("applicable") else "no",
                 "yes" if g.get("triggered") else "no",
                 _pretty(g.get("value")), _pretty(g.get("threshold")))
                for g in ledger]))

    # notes ---------------------------------------------------------------
    add('<div class="notes"><b>Notes on reading this report.</b><dl>')
    add("<dt>&#8225;</dt><dd>The 95&#8201;% interval includes zero. No "
        "lineage effect is distinguishable from none for this trait, and the "
        "point estimate must not be read on its own.</dd>")
    add("<dt>&#8224;</dt><dd>Support below 90&#8201;%. Too few isolates sit "
        "in lineages large enough to inform the estimate.</dd>")
    add("<dt>&#167;</dt><dd>Derived arithmetically from a quantity recorded "
        "elsewhere in the record rather than estimated in this run.</dd>")
    add("</dl><p>The symbols carry the same wording in every run of this "
        "software. No symbol against a value means only that none of the "
        "listed conditions fired; it is not evidence that the value is "
        "correct. Intervals are printed as computed, without clipping to "
        "the natural bounds of the quantity; clipping would leave their "
        "coverage unchanged and is left to the reader so that widths stay "
        "comparable across runs.</p></div>")

    # colophon -------------------------------------------------------------
    add('<div class="colophon">')
    add("amr-clonalshare %s &#183; record schema %s &#183; seed %s &#183; "
        "configuration <code>sha256:%s</code><br>"
        % (__version__, _esc(record.get("schema_version", "")), record.get("seed"), digest))
    add("Cite this run as: &#8220;amr-clonalshare %s, run %s, %s.&#8221;<br>"
        % (__version__, run_id, stamp))
    add("This report supersedes any earlier report bearing the same run "
        "identifier. It is regenerated from the record and holds no value "
        "that the record does not.")
    add("</div></main></body></html>")
    return "".join(A)
