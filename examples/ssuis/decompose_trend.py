#!/usr/bin/env python3
"""Split a change in non-wild-type prevalence into mix and rate components.

Runs :func:`amr_clonalshare.clonality.decompose_panel` over the shipped
*S. suis* panel for a contrast the caller chooses, controls the false discovery
rate within each component family, and reports the effective number of
independent agents beside the count of discoveries.

    python examples/ssuis/decompose_trend.py                     # UK periods
    python examples/ssuis/decompose_trend.py --lineage mlst      # the same
    python examples/ssuis/decompose_trend.py --contrast country  #    contrast
    python examples/ssuis/decompose_trend.py --json out.json     #    two ways

**The lineage column is the analysis decision, not a detail.** Two are shipped.
``baps_cluster`` is the hierBAPS population cluster reported by the source
study for every isolate; ``mlst`` is the sequence type served by BV-BRC, which
is missing for 219 of 677 isolates and missing differentially. Running the
United Kingdom period contrast both ways is the shortest demonstration of why
that matters, and of what the estimability gate is for:

* on ``baps_cluster``, 99 % of both periods carry a label, 17 of 26 clusters
  are shared, and the components are estimable. Ceftiofur non-wild-type
  prevalence rises 9.8 points, of which 13.3 points is a change in which
  lineages were sampled and -3.5 points a change in rate within them.
* on ``mlst``, 40 % of the later period carries a label against 97 % of the
  earlier one, and the unlabelled isolates are the more resistant. The same
  agent then appears to *fall* by 11 points. The gate fires and names the
  reason, because a difference computed on the labelled isolates is not the
  difference in the collection they came from.

The second reading is not a second opinion. The collection's own ceftiofur
prevalence, which needs no lineage label at all, rises by 9.5 points.

``--agents published_cutoff`` restricts the panel to the three agents that have
a published EUCAST cut-off, and asks whether a conclusion survives without the
ten cohort-derived ones. The two contrasts answer differently, and the
difference is the honest limitation of this cohort:

* the **country** conclusion survives. All three anchored agents show a
  within-lineage rate difference at q below 0.001, so "the same lineages are
  more often non-wild-type in Canada" does not rest on derived cut-offs.
* the **period** conclusion cannot be tested this way at all. The four agents
  that carried it are ceftiofur, penicillin, tiamulin and spectinomycin, and
  EUCAST publishes a cut-off for none of them in this species. The three
  anchored agents were not discoveries in the full panel either, so restricting
  to them removes the finding rather than testing it.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from amr_clonalshare.clonality import decompose_panel

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
#: The three agents with a published EUCAST (T)ECOFF. The other ten carry a
#: cut-off derived from this cohort's own MIC distribution, which is defensible
#: for this collection and is not a species-level value. Restricting the panel
#: to these three is the sensitivity analysis that asks whether a conclusion
#: survives without the derived cut-offs.
PUBLISHED_CUTOFF_AGENTS = ["doxycycline", "tetracycline", "erythromycin"]

CONTRASTS = {
    "period": ("collection_year >= 2013", "collection_year < 2013",
               "UK 2013-2014", "UK 2009-2011", "isolation_country"),
    "country": ("isolation_country == 'Canada'",
                "isolation_country == 'United Kingdom'",
                "Canada", "United Kingdom", None),
}


def load() -> tuple[pd.DataFrame, list[str]]:
    meta = pd.read_csv(DATA / "metadata.csv")
    ribo = pd.read_csv(DATA / "ribo.csv")
    cell = pd.read_csv(DATA / "cell.csv")
    frame = meta.merge(ribo, on="genome_id").merge(cell, on="genome_id")
    agents = [c for c in list(ribo.columns) + list(cell.columns)
              if c != "genome_id"]
    return frame, agents


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--contrast", choices=sorted(CONTRASTS),
                        default="period")
    parser.add_argument("--lineage", default="baps_cluster",
                        choices=["baps_cluster", "mlst"])
    parser.add_argument("--agents", default="all",
                        choices=["all", "published_cutoff"],
                        help="restrict the panel to the three agents with a "
                             "published EUCAST cut-off")
    parser.add_argument("--n-boot", type=int, default=4000)
    parser.add_argument("--q", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=23)
    parser.add_argument("--json", type=Path, default=None)
    args = parser.parse_args()

    frame, agents = load()
    if args.agents == "published_cutoff":
        agents = [a for a in agents if a in PUBLISHED_CUTOFF_AGENTS]
    query_a, query_b, name_a, name_b, restrict = CONTRASTS[args.contrast]
    if restrict == "isolation_country":
        frame = frame[(frame["isolation_country"] == "United Kingdom")
                      & frame["collection_year"].notna()]
    side_a, side_b = frame.query(query_a), frame.query(query_b)

    result = decompose_panel(
        side_a[agents], side_a[args.lineage],
        side_b[agents], side_b[args.lineage],
        n_boot=args.n_boot, q=args.q, rng=np.random.default_rng(args.seed))
    family = result["family"]
    per_agent = result["per_agent"]

    print(f"{name_a} (n={len(side_a)})  versus  {name_b} (n={len(side_b)})"
          f"   lineage: {args.lineage}   agents: {args.agents}")
    print(f"{'agent':14s} {'A%':>6s} {'B%':>6s} {'diff':>7s} "
          f"{'mix':>8s} {'[95% CI]':>16s} {'q':>6s} "
          f"{'rate':>8s} {'[95% CI]':>16s} {'q':>6s} {'turn':>5s} {'est':>4s}")
    for agent in agents:
        record = per_agent[agent]
        if record["status"] != "ok":
            print(f"{agent:14s} {record['reason']}")
            continue
        print(f"{agent:14s} {record['prevalence_a'] * 100:6.1f} "
              f"{record['prevalence_b'] * 100:6.1f} "
              f"{record['difference'] * 100:+7.1f} "
              f"{record['composition'] * 100:+8.1f} "
              f"[{record['composition_ci95'][0] * 100:+6.1f},"
              f"{record['composition_ci95'][1] * 100:+6.1f}] "
              f"{record['composition_q']:6.3f} "
              f"{record['within_lineage'] * 100:+8.1f} "
              f"[{record['within_lineage_ci95'][0] * 100:+6.1f},"
              f"{record['within_lineage_ci95'][1] * 100:+6.1f}] "
              f"{record['within_lineage_q']:6.3f} "
              f"{record['turnover_share']:5.2f} "
              f"{str(record['within_lineage_estimable'])[0]:>4s}")

    first = next(r for r in per_agent.values() if r["status"] == "ok")
    availability = first["lineage_label_availability"]
    print(f"\nfamily: {family['n_agents']} agents, effective independent "
          f"{family['effective_independent_agents']:.2f}; "
          f"false discovery rate controlled at q = {family['q']}")
    print(f"  composition   nominal {family['n_composition_nominal']:2d}"
          f"  ->  {family['n_composition_discoveries']:2d} discoveries")
    print(f"  within-lineage nominal {family['n_within_lineage_nominal']:2d}"
          f"  ->  {family['n_within_lineage_discoveries']:2d} discoveries")
    print(f"\nlineages {first['n_lineages']}: {first['n_lineages_shared']} "
          f"shared, {first['n_lineages_only_a']} only in {name_a}, "
          f"{first['n_lineages_only_b']} only in {name_b}; "
          f"{first['shared_support_isolate_share'] * 100:.1f} % of isolates "
          f"sit in shared lineages")
    print(f"lineage labels present: {availability['label_coverage_a'] * 100:.1f} % "
          f"of {name_a}, {availability['label_coverage_b'] * 100:.1f} % of "
          f"{name_b} (difference p = "
          f"{availability['label_coverage_differs_p']:.2g})")
    if not availability["labels_representative"]:
        print("  GATE: labels are differentially missing and the missingness "
              "is associated with the trait; these components describe the "
              "labelled isolates and not the collections")
    print(f"largest identity residual: "
          f"{max(abs(r['identity_residual']) for r in per_agent.values() if r['status'] == 'ok'):.1e}")
    print("  the contrast is described, not explained: collections that "
          "differ by period also differ in\n  submission behaviour and panel "
          "practice, and neither component is a causal effect")

    if args.json:
        args.json.write_text(json.dumps(result, indent=2, default=float))
        print(f"\nwritten: {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
