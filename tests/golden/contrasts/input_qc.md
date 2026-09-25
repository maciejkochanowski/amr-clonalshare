# Input check

Isolates in the raw call/MIC cohort: 96; antimicrobials recorded in the call table: 1.
Susceptibility rows read: 96; 0 unrecognized call(s), 0 missing call(s), and 0 intermediate call(s) set aside. Isolates without any readable call: 0.

## Antimicrobials

Antimicrobials with at least 20 isolates of the rarer outcome: 1 of 1. This count is a warning threshold; each method reports whether its own input requirements are met.

## Metadata

96 isolates (100.0 %) have a metadata row.

## Lineage groups (`lineage`)

12 lineages among 96 typed isolates (0 untyped). 0 lineages hold a single isolate; a single isolate cannot show how much a trait varies inside its lineage, so those isolates do not contribute to that part of the estimate.
Share of isolates in lineages of at least 2 (support): 100.0 %; the estimator needs 80.0 %. The clonal share can be estimated on this cohort.
Effective lineage size for the between-lineage variance: 8.0.


## Per-agent feasibility

Each row uses only isolates with both a readable result for that agent and a recorded lineage. The last column is an algebraic support scenario: one genuinely new tested isolate in each of that many distinct singleton lineages. It does not guarantee interval precision, coverage, or overall estimability; it cannot repair a constant trait or a missing lineage contrast. Duplicating existing rows adds no evidence. The rarer-outcome count is a warning threshold, not an additional estimator gate.

| Agent | Retained | Support | Input failures | Rarer outcome | New isolates for support |
| --- | --- | --- | --- | --- | --- |
| demo_agent | 96 | 100.0 % | none at input level | 39 | 0 |