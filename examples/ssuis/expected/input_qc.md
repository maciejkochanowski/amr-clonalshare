# Input check

Isolates in the raw call/MIC cohort: 677; antimicrobials recorded in the call table: 13.
Susceptibility rows read: 8801; 0 unrecognized call(s), 0 missing call(s), and 0 intermediate call(s) set aside. Isolates without any readable call: 0.
Recorded dilutions: 10832 rows for 677 of 677 isolates, 16 antimicrobials.

## Antimicrobials

Antimicrobials with at least 20 isolates of the rarer outcome: 13 of 13. This count is a warning threshold; each method reports whether its own input requirements are met.

## Metadata

677 isolates (100.0 %) have a metadata row.

## Lineage groups (`baps_cluster`)

30 lineages among 677 typed isolates (0 untyped). 3 lineages hold a single isolate; a single isolate cannot show how much a trait varies inside its lineage, so those isolates do not contribute to that part of the estimate.
Share of isolates in lineages of at least 2 (support): 99.6 %; the estimator needs 90.0 %. The clonal share can be estimated on this cohort.
Effective lineage size for the between-lineage variance: 21.1.


## Per-agent feasibility

Each row uses only isolates with both a readable result for that agent and a recorded lineage. The last column is an algebraic support scenario: one genuinely new tested isolate in each of that many distinct singleton lineages. It does not guarantee interval precision, coverage, or overall estimability; it cannot repair a constant trait or a missing lineage contrast. Duplicating existing rows adds no evidence. The rarer-outcome count is a warning threshold, not an additional estimator gate.

| Agent | Retained | Support | Input failures | Rarer outcome | New isolates for support |
| --- | --- | --- | --- | --- | --- |
| doxycycline | 677 | 99.6 % | none at input level | 106 | 0 |
| tetracycline | 677 | 99.6 % | none at input level | 103 | 0 |
| erythromycin | 677 | 99.6 % | none at input level | 311 | 0 |
| tylosin | 677 | 99.6 % | none at input level | 311 | 0 |
| tilmicosin | 677 | 99.6 % | none at input level | 313 | 0 |
| lincomycin | 677 | 99.6 % | none at input level | 253 | 0 |
| tiamulin | 677 | 99.6 % | none at input level | 131 | 0 |
| spectinomycin | 677 | 99.6 % | none at input level | 78 | 0 |
| amoxicillin | 677 | 99.6 % | none at input level | 37 | 0 |
| cefquinome | 677 | 99.6 % | none at input level | 26 | 0 |
| ceftiofur | 677 | 99.6 % | none at input level | 163 | 0 |
| penicillin | 677 | 99.6 % | none at input level | 156 | 0 |
| trimethoprim | 677 | 99.6 % | none at input level | 190 | 0 |