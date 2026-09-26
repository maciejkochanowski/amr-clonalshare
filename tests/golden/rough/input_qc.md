# Input check

Isolates in the raw call/MIC cohort: 25; antimicrobials recorded in the call table: 3.
Susceptibility rows read: 74; 24 unrecognized call(s), 1 missing call(s), and 0 intermediate call(s) set aside. Isolates without any readable call: 0.
Values that could not be read as a call: `nt` (24). The declared vocabulary is `clinical_sir`.

## Antimicrobials

Antimicrobials with at least 20 isolates of the rarer outcome: 0 of 3. This count is a warning threshold; each method reports whether its own input requirements are met.
Without any readable call: agent_c.
Below the threshold: agent_a.
With a single value (nothing to explain): agent_b.

## Metadata

24 isolates (96.0 %) have a metadata row. Identifiers without one include ['iso_99'].

## Lineage groups (`lineage`)

6 lineages among 22 typed isolates (3 untyped). 2 lineages hold a single isolate; a single isolate cannot show how much a trait varies inside its lineage, so those isolates do not contribute to that part of the estimate.
Share of isolates in lineages of at least 2 (support): 90.9 %; the clonal share is scored on these isolates and the singletons are set aside. The clonal share can be estimated on this cohort.
Effective lineage size for the between-lineage variance: 3.5.


## Per-agent feasibility

Each row uses only isolates with both a readable result for that agent and a recorded lineage. Support is the share of them in lineages of at least two isolates; the share is scored on those isolates and the singletons are set aside. The rarer-outcome count is a warning threshold, not an estimator gate.

| Agent | Retained | Support | Input failures | Rarer outcome |
| --- | --- | --- | --- | --- |
| agent_a | 22 | 90.9 % | none at input level | 9 (below 20) |
| agent_b | 22 | 90.9 % | constant trait | 0 (below 20) |
| agent_c | 0 | not defined | fewer than 2 tested and typed isolates; fewer than 2 lineages with at least 2 isolates | 0 (below 20) |