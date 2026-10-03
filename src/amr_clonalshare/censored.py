"""MIC readings as intervals on the log2 scale, and what the panel tested.

A reading at a tested dilution means the MIC lies above the previous tested
dilution and at or below this one; a reading on an end well is censored.
These intervals are what the ordering share (:mod:`.mic_order`) scores and
what the bounds on the latent ordering (:mod:`.latent_order`) read. No
distribution is assumed for the MICs.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any, Dict, Sequence, Tuple

import numpy as np

__all__ = ["PanelGeometry", "panel_geometry", "intervals_from_mic",
           "intervals_by_panel", "intervals_from_binary"]


@dataclass(frozen=True)
class PanelGeometry:
    """What the panel tested, before any biology.

    ``n_wells`` counts the wells the panel is taken to have tested: every
    doubling from the lowest to the highest recorded value when the recorded
    values sit on a doubling lattice (``lattice == "doubling"``), the distinct
    recorded values otherwise (``lattice == "recorded"``), or the configured
    wells (``lattice == "given"``). ``n_wells_recorded`` counts the distinct
    recorded values before rounding variants of one dilution, such as 0.06
    and 0.064, are folded onto the same well. ``share_on_lowest`` and
    ``share_on_highest`` are the shares of readings on the end wells, which
    do not order their isolates.
    """

    n_wells: int
    lowest: float
    highest: float
    share_on_lowest: float
    share_on_highest: float
    lattice_ratios: Tuple[float, ...]
    doubling: bool
    n_wells_recorded: int = 0
    lattice: str = "recorded"

    def as_dict(self) -> Dict[str, Any]:
        return asdict(self)


LATTICE_TOLERANCE = 0.25


def _lattice(v: np.ndarray, wells=None):
    """The tested wells on the log2 scale and the index of each reading.

    With ``wells`` given, every reading must match a tested concentration.
    Without them, readings that all sit within a quarter of a doubling of an
    integer power of two are taken to come from a doubling panel, whose
    tested wells are every doubling from the lowest to the highest recorded
    value: a dilution no isolate landed on was still tested, and 0.06 and
    0.064 are one well recorded two ways. Readings off any such lattice are
    taken as their own wells, which is the reading of a panel whose steps are
    not known.
    """
    lv = np.log2(v)
    if wells is not None:
        w = np.asarray(wells, dtype=float)
        if (w.ndim != 1 or w.size == 0 or not np.isfinite(w).all()
                or np.any(w <= 0) or np.any(np.diff(w) <= 0)):
            raise ValueError("wells must be strictly increasing positive finite concentrations")
        # A reading matches the well it was read on, written either way:
        # 0.12 and 0.125 are one well, as are 0.015 and 0.015625. A reading
        # farther than a quarter of a doubling from every well is on none.
        lattice = np.log2(w)
        distance = np.abs(lv[:, None] - lattice[None, :])
        if not (distance.min(axis=1) <= LATTICE_TOLERANCE).all():
            raise ValueError("MIC values are outside configured wells; no panel snapping is allowed")
        idx = distance.argmin(axis=1)
        return lattice, idx, "given"
    k = np.round(lv)
    if v.size and np.all(np.abs(lv - k) <= LATTICE_TOLERANCE):
        lo_k, hi_k = int(k.min()), int(k.max())
        lattice = np.arange(lo_k, hi_k + 1, dtype=float)
        return lattice, (k - lo_k).astype(int), "doubling"
    lattice = np.unique(lv)
    return lattice, np.searchsorted(lattice, lv), "recorded"


def panel_geometry(values: Sequence[float], wells=None) -> PanelGeometry:
    """Describe the dilution panel the readings were made on.

    The lattice ratios between adjacent recorded values are reported as they
    are: ratios near 1.03 to 1.05 are rounding variants of one dilution and a
    ratio of 8 is a pair of dilutions no isolate landed on. Neither is a
    property of the panel, and :func:`_lattice` decides what the tested wells
    are; this function reports both readings so that a reader can see what
    was recorded and what was inferred from it.
    """
    v = np.asarray([x for x in np.asarray(values, dtype=float)
                    if np.isfinite(x) and x > 0])
    if v.size == 0:
        return PanelGeometry(0, float("nan"), float("nan"), float("nan"),
                             float("nan"), (), False, 0, "recorded")
    recorded = np.unique(v)
    ratios = tuple(round(float(r), 3)
                   for r in np.unique(recorded[1:] / recorded[:-1]))
    lattice, idx, kind = _lattice(v, wells)
    doubling = kind == "doubling" or (
        kind == "given" and lattice.size > 1
        and bool(np.allclose(np.diff(lattice), 1.0)))
    return PanelGeometry(int(lattice.size), float(2.0 ** lattice[0]), float(2.0 ** lattice[-1]),
                         float((idx == 0).mean()), float((idx == lattice.size - 1).mean()),
                         ratios, doubling, int(recorded.size), kind)


def intervals_from_binary(y, cutoff_log2: float = 0.0):
    """A dichotomised call as the interval it actually is.

    Wild type means the latent value did not exceed the cut-off; non-wild type
    means it did. Nothing else is known. This is the coarsest member of the
    family and exists so the same estimator can consume a collection that reports
    only an interpretation.
    """
    y = np.asarray(y, dtype=float).ravel()
    if not np.isfinite(cutoff_log2):
        raise ValueError("cutoff_log2 must be finite")
    missing = np.isnan(y)
    if not np.isin(y[~missing], (0.0, 1.0)).all():
        raise ValueError("values must be binary numeric 0/1 or NaN")
    lo = np.where(y == 1.0, cutoff_log2, -np.inf)
    hi = np.where(y == 1.0, np.inf, cutoff_log2)
    lo[missing] = np.nan
    hi[missing] = np.nan
    return lo, hi


#: Other spellings of the two inclusive signs, as instrument exports write them.
_SIGNS = {"\u2264": "<=", "\u2265": ">="}


def _normalise_operators(operators):
    """Recorded signs or an empty cell; an unknown sign is not an assumption."""
    op = np.asarray(operators, dtype=object).ravel()
    cleaned = np.array(["" if v is None or str(v).strip().lower() in
                        ("", "nan", "none", "na", "<na>") else _SIGNS.get(str(v).strip(), str(v).strip())
                        for v in op], dtype=str)
    invalid = sorted(set(cleaned) - {"", "=", "<", "<=", ">", ">="})
    if invalid:
        raise ValueError(f"unknown MIC operator(s): {invalid}; use =, <, <=, >, >= or an empty cell")
    return cleaned


def intervals_from_mic(values, *, operators=None, treat_end_wells_as_censored=True,
                       wells=None):
    """Intervals on the log2 scale from recorded MIC values.

    A reading at a tested dilution means the true value lies above the
    previous tested dilution and at or below this one. A reading on the lowest
    well is left-censored and one on the highest is right-censored, unless
    ``treat_end_wells_as_censored`` is turned off, which is the sensitivity
    arm for the coarsening assumption.

    ``wells`` are the concentrations the panel tested, when the laboratory
    recorded them; every reading is then placed on the nearest tested well.
    Without them the tested wells are inferred by :func:`_lattice`: every
    doubling from the lowest to the highest recorded value when the readings
    sit on a doubling lattice, so that a dilution no isolate landed on still
    bounds its neighbours and a rounding variant of one dilution is not read
    as a well of its own, and the recorded values themselves otherwise.

    ``operators`` may carry the recorded censoring signs; where present they
    win over the end-well heuristic, because a recorded operator is an
    observation and the heuristic is an assumption. ``>x`` places the MIC
    above x; ``>=x`` above the tested well below x; ``<x`` and ``<=x`` at or
    below x. The inclusive signs may also be written as the Unicode
    less-than-or-equal and greater-than-or-equal signs. An export that writes
    a sign only on censored readings leaves it blank on an exact one, so a
    blank reading on the highest tested well is the reading of that well,
    (previous well, well], whenever some reading carries ``>`` or ``>=`` at
    that well: growth at the well and no growth at it are two readings, and
    the sign tells them apart. On the lowest well the two coincide (no growth
    at the lowest well is what ``<=`` records), so a blank reading there
    stays as the end-well heuristic reads it.
    """
    v = np.asarray(values, dtype=float).ravel()
    ok = np.isfinite(v) & (v > 0)
    lo = np.full(v.shape, -np.inf)
    hi = np.full(v.shape, np.inf)
    if ok.any():
        lattice, idx, _kind = _lattice(v[ok], wells)
        # The lowest well: left-censored under the assumption, or one
        # doubling wide, the width of every other well, in the sensitivity
        # arm that reads the end wells as exact; the highest well likewise.
        bottom_lo = -np.inf if treat_end_wells_as_censored else lattice[0] - 1.0
        prev = np.where(idx > 0, lattice[np.maximum(idx - 1, 0)], bottom_lo)
        lo[ok] = prev
        hi[ok] = lattice[idx]
        if treat_end_wells_as_censored:
            top = np.zeros(v.shape, dtype=bool)
            top[ok] = idx == lattice.size - 1
            hi[top] = np.inf
    if operators is not None:
        op = _normalise_operators(operators)
        if op.size != v.size:
            raise ValueError(
                f"operators has {op.size} entries and values has {v.size}; "
                f"they must match")
        gt = (op == ">") & ok
        ge = (op == ">=") & ok
        le = np.isin(op, ("<", "<="))
        exact = (op == "=") & ok
        # The tested well below a reading; one doubling below the lowest well.
        below = np.full(v.shape, np.nan)
        full_hi = np.full(v.shape, np.nan)
        if ok.any():
            below[ok] = np.where(idx > 0, lattice[np.maximum(idx - 1, 0)], lattice[0] - 1.0)
            full_hi[ok] = lattice[idx]
        # An exact MIC call is a finite dilution interval, including at an end
        # well. Only missing signs invoke the end-well heuristic.
        lo[exact], hi[exact] = below[exact], full_hi[exact]
        # A blank sign on the highest well is an exact reading where the same
        # call carries > or >= at that well: the export wrote the sign only
        # on the readings with growth at the well, so a blank one is a
        # reading without it, (previous well, well].
        blank = (op == "") & ok
        if ok.any():
            top = np.zeros(v.shape, dtype=bool)
            top[ok] = idx == lattice.size - 1
            if np.any(top & (gt | ge)):
                lo[top & blank], hi[top & blank] = below[top & blank], full_hi[top & blank]
        # >x: growth at x, so the MIC lies above it. >=x: the MIC is x or
        # more, so it lies above the well below x.
        lo[gt] = np.log2(v[gt])
        hi[gt] = np.inf
        lo[ge] = below[ge]
        hi[ge] = np.inf
        # <x and <=x: at or below x. Exports print the lowest well either
        # way, and the interval contains the strict reading of <x.
        lo[le] = -np.inf
        hi[le] = np.log2(np.where(v[le] > 0, v[le], 1.0)) if le.any() else hi[le]
    lo[~ok] = np.nan
    hi[~ok] = np.nan
    return lo, hi


def intervals_by_panel(values, panel, *, operators=None,
                       treat_end_wells_as_censored=True, wells=None):
    """:func:`intervals_from_mic` applied within each panel.

    ``panel`` names the panel every reading was made on. The tested wells,
    and with them the end wells, are inferred from the readings of one panel
    at a time, so that the lowest well of one laboratory is not taken for an
    interior dilution because another laboratory tested lower. Recorded
    ``wells`` apply to every panel alike.
    """
    v = np.asarray(values, dtype=float).ravel()
    labels = np.asarray([str(x) for x in np.asarray(panel, dtype=object).ravel()])
    if labels.size != v.size:
        raise ValueError("panel must name the panel of every reading")
    ops = None if operators is None else np.asarray(operators, dtype=object).ravel()
    lo, hi = np.empty(v.size), np.empty(v.size)
    for name in np.unique(labels):
        m = labels == name
        lo[m], hi[m] = intervals_from_mic(
            v[m], operators=None if ops is None else ops[m],
            treat_end_wells_as_censored=treat_end_wells_as_censored, wells=wells)
    return lo, hi
