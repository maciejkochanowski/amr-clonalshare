"""The run report as one self-contained page, with the diagnostics drawn.

Why a second form of the report
-------------------------------
``report.md`` states the result in prose and is the form that diffs cleanly
between runs. It is not the form in which a reader sees what the numbers are
doing: a panel of intervals that all cross zero is seen in one look when drawn
against a ruled zero, and counted only with effort when printed. The page is
the same report with its figures drawn, built from the same structure as the
prose, so the two cannot disagree.

What is drawn
-------------
Six figures, one per question the sections ask: the lineage sizes behind the
support gate, the share of every trait with both its intervals, the evidence
per trait against the two thresholds, the running product over intakes, the
share read from the recorded dilution beside the share of the call, and the two
prevalences with the concentration of carriage. A figure is drawn only where
the record holds what it needs.

No dependency
-------------
The figures are written as SVG here rather than drawn with a
plotting library, so a report costs the reader nothing to open and the package
keeps matplotlib optional. Coordinates are rounded to two decimals so that two
reports of the same run diff to nothing.
"""
from __future__ import annotations

import math
import re
from typing import Optional

from .report_model import Block, build_report

__all__ = ["render_html_report"]

INK, MID, FAINT = "#1A1A1A", "#4A4A4A", "#8A8A8A"
RULE, PAPER = "#C8C8C4", "#FDFDFB"
ACCENT, FLAG = "#12436D", "#D55E00"


def _finite(x) -> bool:
    return (isinstance(x, (int, float)) and not isinstance(x, bool)
            and math.isfinite(float(x)))


def _esc(s) -> str:
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;"))


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


def _left(names, base: int = 78, px: float = 5.3, cap: int = 170) -> int:
    """Left margin wide enough for the longest row label at the row font."""
    longest = max((len(_short(n)) for n in names), default=0)
    return int(min(cap, max(base, 10 + px * longest)))


def _short(name, n: int = 30) -> str:
    name = str(name)
    return name if len(name) <= n else name[:n - 1] + "\u2026"


# ------------------------------------------------------------- figures ----
def _fig_intervals(rows, n_show: int = 18) -> str:
    """Point estimate with interval, one row per trait, zero ruled."""
    rows = rows[:n_show]
    if not rows:
        return ""
    W, R, T, ROW, B = 470, 58, 14, 15.0, 38
    L = _left([r["name"] for r in rows])
    H = T + ROW * len(rows) + B
    lo = min(min(r["lo"], r["sp_lo"] if _finite(r.get("sp_lo")) else r["lo"])
             for r in rows)
    hi = max(max(r["hi"], r["sp_hi"] if _finite(r.get("sp_hi")) else r["hi"])
             for r in rows)
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
                 % (_r(L - 8), _r(y + 3), _esc(_short(r["name"]))))
        if _finite(r.get("sp_lo")) and _finite(r.get("sp_hi")):
            p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                     'stroke-width="0.6" opacity="0.7"/>'
                     % (_r(X(r["sp_lo"])), _r(y + 4.5), _r(X(r["sp_hi"])),
                        _r(y + 4.5), col))
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


def _frame(p, W, H, L, R, T, B, xlo, xhi, X, axis, zero=True, ticks=5):
    """Vertical rules at nice ticks, the tick labels and the axis title."""
    for t in _ticks(xlo, xhi, ticks):
        p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                 'stroke-width="0.4" shape-rendering="crispEdges"/>'
                 % (_r(X(t)), _r(T - 4), _r(X(t)), _r(H - B + 2), RULE))
        p.append('<text x="%s" y="%s" text-anchor="middle" class="tk">%s'
                 '</text>' % (_r(X(t)), _r(H - B + 15),
                              ("%g" % t) if abs(t) >= 1 or t == 0 else "%.2f" % t))
    if zero and xlo <= 0.0 <= xhi:
        p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                 'stroke-width="0.8" shape-rendering="crispEdges"/>'
                 % (_r(X(0)), _r(T - 4), _r(X(0)), _r(H - B + 2), MID))
    p.append('<text x="%s" y="%s" text-anchor="middle" class="ax">%s</text>'
             % (_r((L + W - R) / 2), _r(H - 8), _esc(axis)))


def _open(W, H, title) -> list:
    return ['<svg viewBox="0 0 %d %.0f" role="img" '
            'xmlns="http://www.w3.org/2000/svg"><title>%s</title>'
            % (W, H, _esc(title))]


def _fig_lineages(d) -> str:
    """Isolates per lineage, largest first; singletons flagged."""
    sizes = d.get("sizes") or []
    if not sizes:
        return ""
    W, R, T, ROW, B = 470, 62, 14, 11.0, 38
    L = _left([n for n, _ in sizes], px=4.6)
    H = T + ROW * len(sizes) + B
    xhi = max(n for _, n in sizes) * 1.06
    X = lambda v: L + v / xhi * (W - L - R)
    p = _open(W, H, "Isolates per lineage")
    _frame(p, W, H, L, R, T, B, 0.0, xhi, X, "isolates in the lineage",
           zero=False)
    small = int(d.get("min_group") or 2)
    for i, (name, n) in enumerate(sizes):
        y = T + ROW * i
        col = FLAG if n < small else ACCENT
        p.append('<rect x="%s" y="%s" width="%s" height="%s" fill="%s" '
                 'opacity="%s"/>' % (_r(X(0)), _r(y + 1.5), _r(X(n) - X(0)),
                                     _r(ROW - 3), col,
                                     "0.9" if n < small else "0.75"))
        if ROW >= 10 or i % 2 == 0:
            p.append('<text x="%s" y="%s" text-anchor="end" class="rw" '
                     'style="font-size:8.5px">%s</text>'
                     % (_r(L - 6), _r(y + ROW - 2.5), _esc(_short(name))))
        p.append('<text x="%s" y="%s" class="rv" style="font-size:8.5px">%d'
                 '</text>' % (_r(X(n) + 4), _r(y + ROW - 2.5), n))
    if _finite(d.get("support")):
        p.append('<text x="%s" y="%s" text-anchor="end" class="lbl">support '
                 '%s of %s isolates%s</text>'
                 % (_r(W - R + 56), _r(H - B - 4),
                    "%.1f&#8201;%%" % (100 * d["support"]),
                    d.get("n_isolates"),
                    (", %d lineages not drawn" % (d["n_total"] - len(sizes)))
                    if d.get("n_total", len(sizes)) > len(sizes) else ""))
    return "".join(p) + "</svg>"


def _fig_evalues(d) -> str:
    """Natural log of the e-value per trait, with the two thresholds ruled."""
    rows = d.get("rows") or []
    if not rows:
        return ""
    W, R, T, ROW, B = 470, 58, 14, 15.0, 38
    L = _left([r["name"] for r in rows])
    H = T + ROW * len(rows) + B
    vals = [r["log_e"] for r in rows]
    refs = [math.log(1.0 / d["alpha"]) if _finite(d.get("alpha"))
            and d["alpha"] > 0 else None,
            math.log(d["threshold"]) if _finite(d.get("threshold"))
            and d["threshold"] > 0 else None]
    lo = min([0.0] + vals + [v for v in refs if v is not None])
    hi = max([0.0] + vals + [v for v in refs if v is not None])
    pad = (hi - lo) * 0.08 or 1.0
    xlo, xhi = lo - pad, hi + pad
    X = lambda v: L + (v - xlo) / (xhi - xlo) * (W - L - R)
    p = _open(W, H, "Evidence per trait")
    _frame(p, W, H, L, R, T, B, xlo, xhi, X,
           "natural log of the e-value (1 = no evidence)")
    for ref, dash in zip(refs, ("3 2", "")):
        if ref is None:
            continue
        p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                 'stroke-width="0.9"%s/>'
                 % (_r(X(ref)), _r(T - 4), _r(X(ref)), _r(H - B + 2), FLAG,
                    (' stroke-dasharray="%s"' % dash) if dash else ""))
    for i, r in enumerate(rows):
        y = T + ROW * i + ROW / 2
        col = ACCENT if r.get("selected") else MID
        p.append('<text x="%s" y="%s" text-anchor="end" class="rw">%s</text>'
                 % (_r(L - 8), _r(y + 3), _esc(_short(r["name"]))))
        p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                 'stroke-width="1.1"/>'
                 % (_r(X(0)), _r(y), _r(X(r["log_e"])), _r(y), col))
        p.append('<circle cx="%s" cy="%s" r="2.8" fill="%s"/>'
                 % (_r(X(r["log_e"])), _r(y), col))
        p.append('<text x="%s" y="%s" class="rv">%.1f</text>'
                 % (_r(W - R + 6), _r(y + 3), r["log_e"]))
    return "".join(p) + "</svg>"


def _fig_sequential(d) -> str:
    """The running product per trait over the intakes, on the log scale."""
    rows = d.get("rows") or []
    batches = d.get("batches") or []
    if not rows or not batches:
        return ""
    W, H, L, R, T, B = 470, 240, 52, 96, 14, 40
    n = len(batches)
    vals = [v for r in rows for v in r["log_e"] if _finite(v)]
    refs = [math.log(1.0 / d["alpha"]) if _finite(d.get("alpha"))
            and d["alpha"] > 0 else None,
            math.log(d["threshold"]) if _finite(d.get("threshold"))
            and d["threshold"] > 0 else None]
    if not vals:
        return ""
    floor = -25.0
    ylo = max(min([0.0] + vals), floor)
    yhi = max([0.0] + vals + [v for v in refs if v])
    pad = (yhi - ylo) * 0.06 or 1.0
    ylo, yhi = ylo - pad, yhi + pad
    X = lambda i: L + (i / max(n - 1, 1)) * (W - L - R)
    Y = lambda v: T + (yhi - v) / (yhi - ylo) * (H - T - B)
    p = _open(W, H, "Sequential e-value over intakes")
    p.append('<defs><clipPath id="seqclip"><rect x="%s" y="%s" width="%s" '
             'height="%s"/></clipPath></defs>'
             % (_r(L - 1), _r(T - 4), _r(W - R - L + 2), _r(H - B - T + 6)))
    for t in _ticks(ylo, yhi, 5):
        p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                 'stroke-width="0.4" shape-rendering="crispEdges"/>'
                 % (_r(L), _r(Y(t)), _r(W - R), _r(Y(t)), RULE))
        p.append('<text x="%s" y="%s" text-anchor="end" class="tk">%g</text>'
                 % (_r(L - 5), _r(Y(t) + 3), t))
    if ylo <= 0.0 <= yhi:
        p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                 'stroke-width="0.8" shape-rendering="crispEdges"/>'
                 % (_r(L), _r(Y(0)), _r(W - R), _r(Y(0)), MID))
    for ref, dash in zip(refs, ("3 2", "")):
        if ref is None or not (ylo <= ref <= yhi):
            continue
        p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                 'stroke-width="0.9"%s/>'
                 % (_r(L), _r(Y(ref)), _r(W - R), _r(Y(ref)), FLAG,
                    (' stroke-dasharray="%s"' % dash) if dash else ""))
    step = max(1, int(math.ceil(n / 8.0)))
    for i, b in enumerate(batches):
        if i % step == 0 or i == n - 1:
            p.append('<text x="%s" y="%s" text-anchor="middle" class="tk">%s'
                     '</text>' % (_r(X(i)), _r(H - B + 14), _esc(b)))
    p.append('<text x="%s" y="%s" text-anchor="middle" class="ax">intake, in '
             'order of arrival</text>' % (_r((L + W - R) / 2), _r(H - 6)))
    p.append('<text transform="translate(11 %s) rotate(-90)" '
             'text-anchor="middle" class="ax">natural log of the running '
             'product</text>' % _r((T + H - B) / 2))
    ordered = sorted(rows, key=lambda r: -(r["log_e"][-1]
                                            if _finite(r["log_e"][-1]) else -9e9))
    labels_y: list = []
    for r in ordered:
        pts = [(X(i), Y(v)) for i, v in enumerate(r["log_e"][:n]) if _finite(v)]
        if not pts:
            continue
        col = ACCENT if r.get("selected") else MID
        p.append('<polyline fill="none" clip-path="url(#seqclip)" stroke="%s" '
                 'stroke-width="%s" opacity="%s" points="%s"/>'
                 % (col, "1.4" if r.get("selected") else "0.9",
                    "1" if r.get("selected") else "0.55",
                    " ".join("%s,%s" % (_r(x), _r(y)) for x, y in pts)))
        yl = min(pts[-1][1], H - B - 2)
        while any(abs(yl - o) < 9.5 for o in labels_y):
            yl += 9.5
        labels_y.append(yl)
        if r.get("selected") or len(rows) <= 6:
            p.append('<text x="%s" y="%s" class="rw" style="fill:%s">%s</text>'
                     % (_r(W - R + 5), _r(yl + 3), col, _esc(_short(r["name"]))))
    return "".join(p) + "</svg>"


def _fig_dilution(rows) -> str:
    """Share read from the recorded dilution, two intervals, with the call
    share of the same agent beside it."""
    if not rows:
        return ""
    W, R, T, ROW, B = 470, 84, 14, 16.0, 38
    L = _left([r["name"] for r in rows])
    H = T + ROW * len(rows) + B
    vals = [v for r in rows for v in (r.get("lo"), r.get("hi"), r.get("r_lo"),
                                       r.get("r_hi"), r.get("call"), r["pt"])
            if _finite(v)]
    lo, hi = min([0.0] + vals), max(vals)
    pad = (hi - lo) * 0.08 or 0.01
    xlo, xhi = lo - pad, hi + pad
    X = lambda v: L + (v - xlo) / (xhi - xlo) * (W - L - R)
    p = _open(W, H, "Share read from the recorded dilution")
    _frame(p, W, H, L, R, T, B, xlo, xhi, X,
           "share of the dilution scale carried by lineage")
    for i, r in enumerate(rows):
        y = T + ROW * i + ROW / 2
        col = ACCENT if r.get("estimable") else MID
        p.append('<text x="%s" y="%s" text-anchor="end" class="rw">%s</text>'
                 % (_r(L - 8), _r(y + 3), _esc(_short(r["name"]))))
        if _finite(r.get("lo")) and _finite(r.get("hi")):
            p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                     'stroke-width="0.7" opacity="0.8"/>'
                     % (_r(X(r["lo"])), _r(y), _r(X(r["hi"])), _r(y), col))
        if _finite(r.get("r_lo")) and _finite(r.get("r_hi")):
            p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                     'stroke-width="2.4"/>'
                     % (_r(X(r["r_lo"])), _r(y), _r(X(r["r_hi"])), _r(y), col))
        p.append('<circle cx="%s" cy="%s" r="2.6" fill="%s"/>'
                 % (_r(X(r["pt"])), _r(y), col))
        if _finite(r.get("call")):
            cx, cy = X(r["call"]), y
            p.append('<polygon points="%s,%s %s,%s %s,%s %s,%s" fill="%s" '
                     'stroke="%s" stroke-width="0.9"/>'
                     % (_r(cx), _r(cy - 3.6), _r(cx + 3.6), _r(cy),
                        _r(cx), _r(cy + 3.6), _r(cx - 3.6), _r(cy), PAPER, col))
        p.append('<text x="%s" y="%s" class="rv">%s%s</text>'
                 % (_r(W - R + 6), _r(y + 3), "%.3f" % r["pt"],
                    (" (%.0f&#8201;%%)" % (100 * r["censored"]))
                    if _finite(r.get("censored")) else ""))
    return "".join(p) + "</svg>"


def _fig_carriage(rows) -> str:
    """Left: prevalence per isolate and per lineage. Right: the effective
    number of carrying lineages against its permutation interval."""
    if not rows:
        return ""
    W, T, ROW, B, GAP = 470, 14, 16.0, 38, 22
    L = _left([r["name"] for r in rows])
    H = T + ROW * len(rows) + B
    L2 = 300 + (L - 78) // 2
    R2 = 30
    pv = [v for r in rows for v in (r.get("per_isolate"), r.get("per_lineage"),
                                     r.get("pl_lo"), r.get("pl_hi")) if _finite(v)]
    ev = [v for r in rows for v in (r.get("eff"), r.get("null_lo"),
                                     r.get("null_hi")) if _finite(v)]
    if not pv:
        return ""
    xhi1 = max(pv) * 1.08 or 0.01
    X1 = lambda v: L + v / xhi1 * (L2 - GAP - L)
    p = _open(W, H, "Prevalence on two scales and the concentration of carriage")
    for t in _ticks(0.0, xhi1, 4):
        p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                 'stroke-width="0.4" shape-rendering="crispEdges"/>'
                 % (_r(X1(t)), _r(T - 4), _r(X1(t)), _r(H - B + 2), RULE))
        p.append('<text x="%s" y="%s" text-anchor="middle" class="tk">%.0f&#8201;%%'
                 '</text>' % (_r(X1(t)), _r(H - B + 15), 100 * t))
    p.append('<text x="%s" y="%s" text-anchor="middle" class="ax">prevalence: '
             'per isolate &#9679;, per lineage &#9675;</text>'
             % (_r((L + L2 - GAP) / 2), _r(H - 8)))
    if ev:
        elo, ehi = 0.0, max(ev) * 1.1
        X2 = lambda v: L2 + (v - elo) / (ehi - elo) * (W - R2 - L2)
        for t in _ticks(elo, ehi, 3):
            p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                     'stroke-width="0.4" shape-rendering="crispEdges"/>'
                     % (_r(X2(t)), _r(T - 4), _r(X2(t)), _r(H - B + 2), RULE))
            p.append('<text x="%s" y="%s" text-anchor="middle" class="tk">%g'
                     '</text>' % (_r(X2(t)), _r(H - B + 15), t))
        p.append('<text x="%s" y="%s" text-anchor="middle" class="ax">effective '
                 'carrying lineages</text>' % (_r((L2 + W - R2) / 2), _r(H - 8)))
    for i, r in enumerate(rows):
        y = T + ROW * i + ROW / 2
        col = ACCENT if r.get("departs") else MID
        p.append('<text x="%s" y="%s" text-anchor="end" class="rw">%s</text>'
                 % (_r(L - 8), _r(y + 3), _esc(_short(r["name"]))))
        a, b = r.get("per_isolate"), r.get("per_lineage")
        if _finite(r.get("pl_lo")) and _finite(r.get("pl_hi")):
            p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                     'stroke-width="0.6" opacity="0.7"/>'
                     % (_r(X1(r["pl_lo"])), _r(y), _r(X1(r["pl_hi"])), _r(y), col))
        if _finite(a) and _finite(b):
            p.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                     'stroke-width="1.6"/>'
                     % (_r(X1(a)), _r(y), _r(X1(b)), _r(y), col))
        if _finite(a):
            p.append('<circle cx="%s" cy="%s" r="2.8" fill="%s"/>'
                     % (_r(X1(a)), _r(y), col))
        if _finite(b):
            p.append('<circle cx="%s" cy="%s" r="2.8" fill="%s" stroke="%s" '
                     'stroke-width="1.1"/>' % (_r(X1(b)), _r(y), PAPER, col))
        if ev and _finite(r.get("null_lo")) and _finite(r.get("null_hi")):
            p.append('<rect x="%s" y="%s" width="%s" height="%s" fill="%s" '
                     'opacity="0.55"/>'
                     % (_r(X2(r["null_lo"])), _r(y - 4.5),
                        _r(X2(r["null_hi"]) - X2(r["null_lo"])), _r(9), RULE))
        if ev and _finite(r.get("eff")):
            p.append('<circle cx="%s" cy="%s" r="2.8" fill="%s"/>'
                     % (_r(X2(r["eff"])), _r(y), col))
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
table{border-collapse:collapse;width:100%;margin:10px 0 6px;font-size:.84rem}
caption{caption-side:top;text-align:left;color:var(--mid);font-size:.82rem;
margin-bottom:5px}
caption b{color:var(--ink);font-weight:600}
thead th{border-top:.8px solid var(--ink);border-bottom:.5px solid var(--ink);
padding:5px 7px;font-weight:600;vertical-align:bottom}
tbody td{padding:4px 7px;border:0}
tbody tr:last-child td{border-bottom:.8px solid var(--ink)}
.num{text-align:right}
td.num{white-space:nowrap}
td.num br+*{font-size:.9em;color:var(--mid)}
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
.ident .wide{grid-column:1 / -1}
.ident code{font-family:"IBM Plex Mono",ui-monospace,Menlo,Consolas,monospace;
font-size:.72rem;color:var(--mid);word-break:break-all}
.callout{border-left:2.5px solid var(--accent);padding:6px 0 2px 12px;
margin:14px 0 12px}
.callout .ct{font-size:.80rem;letter-spacing:.06em;text-transform:uppercase;
color:var(--accent);font-weight:600;margin-bottom:4px}
.callout p{font-size:.92rem;margin:0 0 6px}
ul{margin:4px 0 10px;padding-left:1.2em}
@media print{
 @page{size:A4;margin:22mm}
 html,body{background:#fff}
 .sheet{max-width:none;padding:0}
 h2,figure,tr,.callout,.notes{break-inside:avoid}
 thead{display:table-header-group}
 h2,caption{break-after:avoid}
 p{orphans:3;widows:3}
 svg,table,figure{print-color-adjust:exact;-webkit-print-color-adjust:exact}
}
"""


def _table(caption, headers, rows) -> str:
    """headers: list of (label, css class). rows: sequences of cells."""
    out = ["<table><caption>%s</caption><thead><tr>" % caption]
    out += ['<th class="%s">%s</th>' % (cls, h.replace("-", "&#8209;"))
            for h, cls in headers]
    out.append("</tr></thead><tbody>")
    for row in rows:
        out.append("<tr>")
        out += ['<td class="%s">%s</td>' % (headers[i][1], cell)
                for i, cell in enumerate(row)]
        out.append("</tr>")
    return "".join(out) + "</tbody></table>"



def _inline(text: str) -> str:
    """Escape, then the two inline marks the model allows, and a line break."""
    t = _esc(text)
    t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t)
    t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
    return t.replace("\n", "<br>")


_ALIGN = {"lab": "lab", "num": "num", "flag": "flagcol"}


def _render_block(b: Block, figs: dict) -> str:
    if b.kind == "p":
        return "<p>%s</p>" % _inline(b.text)
    if b.kind == "kv":
        out = ['<div class="dl">']
        for k, v in b.rows:
            out.append("<dt>%s</dt><dd>%s</dd>" % (_inline(str(k)), _inline(str(v))))
        return "".join(out) + "</div>"
    if b.kind == "table":
        headers = [(_inline(h), _ALIGN.get(a, "lab")) for h, a in b.headers]
        rows = [[_inline(str(c)) for c in row] for row in b.rows]
        return _table("<b>%s.</b> %s" % (_esc(b.label), _inline(b.caption)),
                      headers, rows)
    if b.kind == "figure":
        svg = figs.get(b.figure, lambda _data: "")(b.data)
        if not svg:
            return ""
        return ("<figure>%s<figcaption><b>%s.</b> %s</figcaption></figure>"
                % (svg, _esc(b.label), _inline(b.caption)))
    if b.kind == "callout":
        out = ['<div class="callout"><div class="ct">%s</div>' % _inline(b.title)]
        for line in b.lines:
            out.append("<p>%s</p>" % _inline(line))
        return "".join(out) + "</div>"
    return ""


def render_html_report(record: dict, summary: dict, *,
                       input_qc: Optional[dict] = None,
                       record_bytes: Optional[bytes] = None) -> str:
    """One self-contained page for a run, from the record it wrote.

    Every value is read from ``record`` or ``summary``; none is recomputed,
    so the page cannot disagree with the file it describes.
    """
    rep = build_report(record, summary, input_qc=input_qc,
                       record_bytes=record_bytes)
    figs = {
        "intervals": lambda rows: _fig_intervals(rows),
        "lineages": _fig_lineages,
        "evalues": _fig_evalues,
        "sequential": _fig_sequential,
        "dilution": _fig_dilution,
        "carriage": _fig_carriage,
    }
    A: list[str] = []
    add = A.append
    add('<!doctype html><html lang="en"><head><meta charset="utf-8">')
    add('<meta name="viewport" content="width=device-width,initial-scale=1">')
    add("<title>Run report %s</title><style>%s</style></head><body>"
        % (rep.ident["run"], _CSS))
    add('<main class="sheet">')
    add("<h1>The share of resistance the lineages carry</h1>")
    add('<p class="sub">Run report. Every value is read from the record '
        "written by this run; none is recomputed here.</p>")
    add('<div class="ident">')
    for key, label in (("run", "Run"), ("issued", "Issued"),
                       ("software", "Software"), ("record", "Record"),
                       ("isolates", "Isolates"), ("lineage", "Lineage"),
                       ("seed", "Seed"), ("configuration", "Configuration")):
        add("<div><b>%s</b>%s</div>" % (label, _esc(rep.ident[key])))
    if rep.ident.get("record_sha256"):
        add('<div class="wide"><b>Record digest</b><code>sha256:%s</code></div>'
            % _esc(rep.ident["record_sha256"]))
    add("</div>")
    st = rep.status
    add('<div class="status"><span class="lead%s">%s</span><p>%s</p></div>'
        % (" flag" if st["flag"] else "", _esc(st["lead"]),
           _esc(st["gates_line"])))
    n = 0
    for sec in rep.sections:
        n += 1
        add('<h2><span class="no">%d</span>%s</h2>' % (n, _inline(sec.title)))
        for b in sec.blocks:
            add(_render_block(b, figs))
    add('<div class="notes"><b>Notes on reading this report.</b><dl>')
    for sym, text in rep.symbols:
        add("<dt>%s</dt><dd>%s</dd>" % (_esc(sym), _inline(text)))
    add("</dl></div>")
    add('<div class="colophon">')
    add("<br>".join(_inline(line) for line in rep.colophon))
    add("</div></main></body></html>")
    return "".join(A)
