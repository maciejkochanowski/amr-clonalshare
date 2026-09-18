# Input check

Isolates in the raw call/MIC cohort: 7049; antimicrobials recorded in the call table: 22.
Susceptibility rows read: 101235; 0 unrecognized call(s), 0 missing call(s), and 0 intermediate call(s) set aside. Isolates without any readable call: 0.

## Antimicrobials

Antimicrobials with at least 20 isolates of the rarer outcome: 17 of 22. This count is a warning threshold; each method reports whether its own input requirements are met.
Below the threshold: amikacin.
With a single value (nothing to explain): meropenem, imipenem, spectinomycin, colistin.

## Metadata

7049 isolates (100.0 %) have a metadata row.

## Lineage groups (`pds_cluster`)

534 lineages among 7049 typed isolates (0 untyped). 207 lineages hold a single isolate; a single isolate cannot show how much a trait varies inside its lineage, so those isolates do not contribute to that part of the estimate.
Share of isolates in lineages of at least 2 (support): 97.1 %; the estimator needs 90.0 %. The clonal share can be estimated on this cohort.
Effective lineage size for the between-lineage variance: 12.9.


## Per-agent feasibility

Each row uses only isolates with both a readable result for that agent and a recorded lineage. The last column is an algebraic support scenario: one genuinely new tested isolate in each of that many distinct singleton lineages. It does not guarantee interval precision, coverage, or overall estimability; it cannot repair a constant trait or a missing lineage contrast. Duplicating existing rows adds no evidence. The rarer-outcome count is a warning threshold, not an additional estimator gate.

| Agent | Retained | Support | Input failures | Rarer outcome | New isolates for support |
| --- | --- | --- | --- | --- | --- |
| amikacin | 2401 | 95.4 % | none at input level | 1 (below 20) | 0 |
| amoxicillin-clavulanic acid | 6917 | 97.3 % | none at input level | 1163 | 0 |
| ampicillin | 7049 | 97.1 % | none at input level | 1864 | 0 |
| cefoxitin | 6917 | 97.3 % | none at input level | 649 | 0 |
| ceftiofur | 3616 | 96.2 % | none at input level | 529 | 0 |
| ceftriaxone | 6917 | 97.3 % | none at input level | 897 | 0 |
| chloramphenicol | 7049 | 97.1 % | none at input level | 425 | 0 |
| ciprofloxacin | 6935 | 97.0 % | none at input level | 606 | 0 |
| gentamicin | 7047 | 97.1 % | none at input level | 1008 | 0 |
| kanamycin | 3206 | 96.1 % | none at input level | 427 | 0 |
| nalidixic acid | 6915 | 97.3 % | none at input level | 584 | 0 |
| streptomycin | 6520 | 96.9 % | none at input level | 2659 | 0 |
| sulfamethoxazole | 253 | 88.9 % | support below 0.90 | 54 | 3 |
| tetracycline | 7047 | 97.1 % | none at input level | 3338 | 0 |
| trimethoprim-sulfamethoxazole | 6915 | 97.3 % | none at input level | 208 | 0 |
| sulfisoxazole | 6662 | 97.3 % | none at input level | 1927 | 0 |
| azithromycin | 4779 | 96.5 % | none at input level | 24 | 0 |
| meropenem | 3300 | 96.5 % | constant trait | 0 (below 20) | 0 |
| ceftazidime | 132 | 81.8 % | support below 0.90 | 37 | 10 |
| imipenem | 132 | 81.8 % | constant trait; support below 0.90 | 0 (below 20) | 10 |
| spectinomycin | 1 | not defined | fewer than 2 tested and typed isolates; fewer than 2 lineages; constant trait | 0 (below 20) | 1 |
| colistin | 525 | 93.5 % | constant trait | 0 (below 20) | 0 |