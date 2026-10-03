"""One conclusion per agent on what the collection showed.

The record holds, per agent, a lineage share with its interval and
permutation p-value and, for recorded MICs, the bounds on the lineage share
of the latent ordering with the lower confidence limit. A reader who is not
a statistician needs the three to be read together once: whether the
collection shows lineage structure in that reading, with the numbers that
say so in parentheses. The verdict is a statement about the collection
analysed; it says nothing about transmission or a resistance mechanism.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, TypeGuard

__all__ = ["agent_verdict", "format_p", "format_interval", "summary_rows", "short_label", "VERDICT_NOTE", "SHORT_LABELS"]

#: The conclusion of a row in a table column, with the sentence of the label
#: given in the caption or the note of the table.
SHORT_LABELS = {"lineage structure established": "established",
                "interval spans the whole range": "whole range"}


def short_label(label: str) -> str:
    return SHORT_LABELS.get(label, label)


VERDICT_NOTE = ("The conclusion says whether this collection shows lineage structure in the "
                "measurement, with the values behind it. It is not a statement about "
                "transmission or a resistance mechanism.")


def _finite(x: object) -> TypeGuard[float]:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(float(x))


def _num(x, nd: int = 3) -> str:
    return ("%%.%df" % nd) % float(x) if _finite(x) else "—"


def format_p(p, floor=None) -> str:
    """A permutation p-value; the smallest value the permutations can give is
    written as a bound, since fewer than one permutation reached the
    observed statistic."""
    if not _finite(p):
        return "—"
    if _finite(floor) and float(p) <= float(floor) * (1 + 1e-9):
        return "≤ " + _num(floor)
    return _num(p)


def format_interval(lo, hi) -> str:
    """A 95 % interval; one that spans the whole range says so."""
    if not (_finite(lo) and _finite(hi)):
        return "—"
    if float(lo) <= 0.0 and float(hi) >= 1.0:
        return "spans the whole range"
    return "%s to %s" % (_num(lo), _num(hi))


def agent_verdict(block: Optional[dict], *, alpha: float = 0.05, q: Optional[float] = None) -> Dict[str, str]:
    """``{"label", "text"}`` for one agent's share block (a ``clonal_share``
    entry, or the ``order`` entry of ``censored_share``).

    Established: the interval of the share lies above zero and the agent is
    selected by the Benjamini-Yekutieli procedure across the panel (``q`` at
    or below ``alpha``; without a ``q``, the permutation p-value), or the
    lower confidence limit of the latent share is above zero. Whole range:
    the interval spans 0 to 1. Not established: any other computed share.
    Not estimable: the collection could not give the share, with the reason
    the record gives.
    """
    block = block or {}
    share = block.get("kappa_adj")
    if block.get("estimable") is False or not _finite(share):
        reason = block.get("reason") or "the share could not be scored on these readings"
        return {"label": "not estimable", "text": "Not estimable: %s." % reason}
    lo, hi = block.get("observed_low"), block.get("observed_high")
    p, floor = block.get("p_value"), block.get("p_floor")
    adjusted = q if _finite(q) else p
    # a call bounds the latent ordering too, with a coarser reading
    lower = block.get("latent_order_lower")
    limit = block.get("latent_order_lower_limit")
    parts = ["share %s" % _num(share), "95 %% interval %s" % format_interval(lo, hi)]
    if _finite(p):
        parts.append("p %s" % format_p(p, floor))
    if _finite(q):
        parts.append("q %s" % _num(q))
    if _finite(lower):
        parts.append("lower bound %s" % _num(lower))
    if _finite(limit):
        parts.append("lower confidence limit %s" % _num(limit))
    detail = "; ".join(parts)
    above = _finite(lo) and float(lo) > 0.0 and _finite(adjusted) and float(adjusted) <= alpha
    if above or (_finite(limit) and float(limit) > 0.0):
        return {"label": "lineage structure established",
                "text": "Lineage structure established (%s)." % detail}
    if _finite(lo) and _finite(hi) and float(lo) <= 0.0 and float(hi) >= 1.0:
        return {"label": "interval spans the whole range",
                "text": "The interval spans the whole range, so the collection does not resolve the share (%s)." % detail}
    return {"label": "not established", "text": "Not established on this collection (%s)." % detail}


def summary_rows(record: dict) -> List[dict]:
    """One row per agent for the results table of the form and the console:
    the call share and the MIC-ordering share with their intervals, p-values,
    q-values, bounds and conclusions, as text."""
    md = record.get("metadata_diagnostics") or {}
    calls = md.get("clonal_share") or {}
    orders = (md.get("censored_share") or {}).get("per_agent") or {}
    call_selection = md.get("lineage_selection") or {}
    order_selection = (md.get("censored_share") or {}).get("order_selection") or {}
    call_q = call_selection.get("q_values") or {}
    order_q = order_selection.get("q_values") or {}
    call_alpha = float(call_selection.get("alpha") or 0.05)
    order_alpha = float(order_selection.get("alpha") or 0.05)
    rows = []
    for agent in sorted(set(calls) | set(orders)):
        row: Dict[str, Any] = {"agent": str(agent)}
        if agent in calls:
            c = calls[agent] or {}
            row.update(call_share=_num(c.get("kappa_adj")) if c.get("estimable", True) else "—",
                       call_interval=format_interval(c.get("observed_low"), c.get("observed_high")),
                       call_p=format_p(c.get("p_value"), c.get("p_floor")),
                       call_q=_num(call_q.get(agent)),
                       call_lower_bound=_num(c.get("latent_order_lower")),
                       call_lower_limit=_num(c.get("latent_order_lower_limit")),
                       call_conclusion=agent_verdict(c, alpha=call_alpha, q=call_q.get(agent)))
        if agent in orders:
            o = (orders[agent] or {}).get("order") or {}
            row.update(mic_share=_num(o.get("kappa_adj")) if o.get("estimable", True) else "—",
                       mic_interval=format_interval(o.get("observed_low"), o.get("observed_high")),
                       mic_p=format_p(o.get("p_value"), o.get("p_floor")),
                       mic_q=_num(order_q.get(agent)),
                       lower_bound=_num(o.get("latent_order_lower")),
                       lower_limit=_num(o.get("latent_order_lower_limit")),
                       mic_conclusion=agent_verdict(o, alpha=order_alpha, q=order_q.get(agent)))
        rows.append(row)
    return rows
