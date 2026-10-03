# Input check

Isolates in the raw call/MIC collection: 7049; antimicrobials recorded in the call table: 22.
Susceptibility rows read: 101235; 0 unrecognized call(s), 0 missing call(s), and 0 intermediate call(s) set aside. Isolates without any readable call: 0.

## Antimicrobials

Antimicrobials with at least 20 isolates of the rarer outcome: 17 of 22. This count is a warning threshold; each method reports whether its own input requirements are met.
Below the threshold: amikacin.
With a single value (nothing to explain): colistin, imipenem, meropenem, spectinomycin.

## Metadata

7049 isolates (100.0 %) have a metadata row.

## Lineage groups (`serovar`)

85 lineages among 7049 typed isolates (0 untyped). 20 lineages hold a single isolate; a single isolate cannot show how much a trait varies inside its lineage, so those isolates do not contribute to that part of the estimate.
Share of isolates in lineages of at least 2 (support): 99.7 %; the lineage share is scored on these isolates and the singletons are set aside. The lineage share can be estimated on this collection.
Effective lineage size for the between-lineage variance: 76.8.


## Per-agent feasibility

Each row uses only isolates with both a readable result for that agent and a recorded lineage. Support is the share of them in lineages of at least two isolates; the share is scored on those isolates and the singletons are set aside. The rarer-outcome count is a warning threshold, not a reporting condition of the estimator.

| Agent | Retained | Support | Input failures | Rarer outcome |
| --- | --- | --- | --- | --- |
| amikacin | 2401 | 99.7 % | none at input level | 1 (below 20) |
| amoxicillin-clavulanic acid | 6917 | 99.7 % | none at input level | 1163 |
| ampicillin | 7049 | 99.7 % | none at input level | 1864 |
| azithromycin | 4779 | 99.7 % | none at input level | 24 |
| cefoxitin | 6917 | 99.7 % | none at input level | 649 |
| ceftazidime | 132 | 96.2 % | none at input level | 37 |
| ceftiofur | 3616 | 99.7 % | none at input level | 529 |
| ceftriaxone | 6917 | 99.7 % | none at input level | 897 |
| chloramphenicol | 7049 | 99.7 % | none at input level | 425 |
| ciprofloxacin | 6935 | 99.7 % | none at input level | 606 |
| colistin | 525 | 97.3 % | constant trait | 0 (below 20) |
| gentamicin | 7047 | 99.7 % | none at input level | 1008 |
| imipenem | 132 | 96.2 % | constant trait | 0 (below 20) |
| kanamycin | 3206 | 99.8 % | none at input level | 427 |
| meropenem | 3300 | 99.5 % | constant trait | 0 (below 20) |
| nalidixic acid | 6915 | 99.7 % | none at input level | 584 |
| spectinomycin | 1 | not defined | fewer than 2 tested and typed isolates; fewer than 2 lineages with at least 2 isolates; constant trait | 0 (below 20) |
| streptomycin | 6520 | 99.7 % | none at input level | 2659 |
| sulfamethoxazole | 253 | 98.0 % | none at input level | 54 |
| sulfisoxazole | 6662 | 99.7 % | none at input level | 1927 |
| tetracycline | 7047 | 99.7 % | none at input level | 3338 |
| trimethoprim-sulfamethoxazole | 6915 | 99.7 % | none at input level | 208 |