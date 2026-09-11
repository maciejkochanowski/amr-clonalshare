#!/usr/bin/env python3
"""One screening rule for both atlases, and a closed set of reasons.

WHY THIS EXISTS. The cross-species atlas and the veterinary atlas ran the same
loop over agents, written twice. Both copies passed over an agent that failed
the tested-isolate threshold without writing a row, so a species could report
fewer agents than it held and nothing in the file said where the rest went: 366
of 907 agents in the veterinary run and 21 in the cross-species run left no
trace. A count that does not reconcile is read as selection, and the reader has
no way to tell it apart from selection. Screening in one place, with an outcome
for every agent, removes the possibility rather than the instance.

THE RULE. An agent is analysed when enough isolates were tested against it, the
trait has both classes present, and the analysed subset holds at least two
lineages. Otherwise it is refused, and the refusal names which of those three
conditions failed. There is no fourth outcome: every agent present in a cohort
leaves exactly one record.

WHY A MINIMUM LINEAGE COUNT. A share is a contrast between lineages. On a
subset holding one lineage there is nothing to contrast, and the estimator
answers with zero, an interval of zero width and p = 1 -- arithmetic that reads
as an absence of lineage effect when nothing was compared. A public release
that assigns a whole species to one cluster produces exactly this, so the case
is reachable on real data.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np

#: Minimum isolates tested against an agent before its share is estimated.
#: Below this the bootstrap interval is wider than the range it is drawn on.
MIN_DRUG_N = 50

#: An agent needs both classes present to have variance to attribute.
MIN_MINOR = 0.02

#: A share needs a contrast between at least two lineages.
MIN_LINEAGES = 2

#: The closed set of refusals. A reason outside this set is a bug, not a new
#: category: the verifier in scripts/check_evidence_accounting.py rejects one.
REFUSAL_TOO_FEW_TESTED = "fewer than %d isolates tested" % MIN_DRUG_N
REFUSAL_NO_VARIANCE = "no variance"
REFUSAL_TOO_FEW_LINEAGES = ("fewer than %d lineages in the analysed subset"
                            % MIN_LINEAGES)
REFUSALS = (REFUSAL_TOO_FEW_TESTED, REFUSAL_NO_VARIANCE,
            REFUSAL_TOO_FEW_LINEAGES)


@dataclass(frozen=True)
class Screen:
    """The outcome for one agent in one cohort.

    ``analyse`` is true for exactly the agents that go on to the estimator.
    When it is false, ``record`` is the row to write and ``reason`` names which
    condition failed. When it is true, ``y`` and ``lineage`` are the analysed
    subset and ``record`` carries the fields every row shares.
    """

    analyse: bool
    record: dict
    reason: str = ""
    y: np.ndarray = field(default_factory=lambda: np.empty(0))
    lineage: tuple = ()


def screen_agent(values, lineage: Sequence, *, min_drug_n: int = MIN_DRUG_N,
                 min_minor: float = MIN_MINOR,
                 min_lineages: int = MIN_LINEAGES) -> Screen:
    """Decide whether one agent can be estimated, and say why when it cannot.

    Parameters
    ----------
    values : per-isolate calls for one agent, NaN where the isolate was not
        tested against it.
    lineage : lineage label per isolate, in the same order.

    Returns
    -------
    Screen
    """
    y = np.asarray(values, dtype=float)
    lineage = list(lineage)
    if y.size != len(lineage):
        raise ValueError(f"{y.size} calls against {len(lineage)} lineage "
                         f"labels; they must describe the same isolates")

    tested = np.isfinite(y)
    n_tested = int(tested.sum())
    if n_tested < min_drug_n:
        return Screen(False, {"n": n_tested, "skipped": REFUSAL_TOO_FEW_TESTED},
                      REFUSAL_TOO_FEW_TESTED)

    subset = y[tested]
    labels = tuple(lineage[i] for i in np.flatnonzero(tested))
    prevalence = float(subset.mean())
    if prevalence < min_minor or prevalence > 1.0 - min_minor:
        return Screen(False, {"n": n_tested, "prevalence": prevalence,
                              "skipped": REFUSAL_NO_VARIANCE},
                      REFUSAL_NO_VARIANCE)

    n_lineages = len(set(labels))
    if n_lineages < min_lineages:
        return Screen(False, {"n": n_tested, "prevalence": prevalence,
                              "n_lineages_analysed": n_lineages,
                              "skipped": REFUSAL_TOO_FEW_LINEAGES},
                      REFUSAL_TOO_FEW_LINEAGES)

    return Screen(True, {"n": n_tested, "prevalence": prevalence},
                  y=subset, lineage=labels)


def reconcile(n_agents_present: int, records: dict) -> None:
    """Every agent present in a cohort must have left exactly one record.

    Raised rather than reported, because a run that cannot account for its own
    units has nothing to report.
    """
    if len(records) != n_agents_present:
        raise AssertionError(
            f"{n_agents_present} agents present in the cohort but "
            f"{len(records)} records written; every agent must leave one")
    unknown = sorted({r["skipped"] for r in records.values()
                      if isinstance(r, dict) and r.get("skipped")}
                     - set(REFUSALS))
    if unknown:
        raise AssertionError(f"refusal reasons outside the closed set: "
                             f"{unknown}")
