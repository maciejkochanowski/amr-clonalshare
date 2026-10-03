# Input check

Isolates in the raw call/MIC collection: 677; antimicrobials recorded in the call table: 13.
Susceptibility rows read: 8801; 0 unrecognized call(s), 0 missing call(s), and 0 intermediate call(s) set aside. Isolates without any readable call: 0.
Recorded dilutions: 10832 rows for 677 of 677 isolates, 16 antimicrobials.

## MIC ranges

| Agent | Recorded concentrations (mg/L) | On the lowest or highest | Tested range declared |
| --- | --- | --- | --- |
| amoxicillin | 0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1, 2, 4 | 7.5 % | no |
| cefquinome | 0.002, 0.004, 0.008, 0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1, 2 | 1.5 % | no |
| ceftiofur | 0.03, 0.06, 0.12, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16 | 1.3 % | no |
| doxycycline | 0.03, 0.031, 0.06, 0.063, 0.12, 0.25, 0.5, 4, 8, 16, 32, 64 | 0.6 % | no |
| enrofloxacin | 0.008, 0.03, 0.06, 0.12, 0.25, 0.5, 1, 2, 4, 8 | 0.4 % | no |
| erythromycin | 0.015, 0.016, 0.03, 0.031, 0.06, 0.12, 0.25, 1, 2, 4, 8, 16, 32 | 56.1 % | no |
| florfenicol | 0.25, 0.5, 1, 2, 4 | 3.2 % | no |
| lincomycin | 0.06, 0.12, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128 | 55.8 % | no |
| marbofloxacin | 0.015, 0.25, 0.5, 1, 2, 16 | 4.3 % | no |
| penicillin | 0.03, 0.031, 0.06, 0.12, 0.25, 0.5, 1, 2, 4, 8, 16 | 65.7 % | no |
| spectinomycin | 1, 2, 4, 8, 16, 32, 64, 128, 256, 512 | 7.8 % | no |
| tetracycline | 0.06, 0.12, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128 | 9.9 % | no |
| tiamulin | 0.03, 0.12, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64 | 8.0 % | no |
| tilmicosin | 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024 | 1.9 % | no |
| trimethoprim | 0.015, 0.03, 0.06, 0.12, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32 | 9.9 % | no |
| tylosin | 0.12, 0.125, 0.25, 0.5, 1, 2, 4, 16, 32, 64, 128, 256 | 51.7 % | no |

Without the tested range, the lowest and the highest recorded concentration of an agent are taken as the end wells of its panel, and a reading on them as censored. For erythromycin (56.1 %), lincomycin (55.8 %), penicillin (65.7 %), tylosin (51.7 %) that is a large share of the readings: if the panel tested further dilutions, declare its range (`mic_wells`, or a shipped panel preset) so that the bounds are read on the wells the laboratory tested.

## Antimicrobials

Antimicrobials with at least 20 isolates of the rarer outcome: 13 of 13. This count is a warning threshold; each method reports whether its own input requirements are met.

## Metadata

677 isolates (100.0 %) have a metadata row.

## Lineage groups (`baps_cluster`)

30 lineages among 677 typed isolates (0 untyped). 3 lineages hold a single isolate; a single isolate cannot show how much a trait varies inside its lineage, so those isolates do not contribute to that part of the estimate.
Share of isolates in lineages of at least 2 (support): 99.6 %; the lineage share is scored on these isolates and the singletons are set aside. The lineage share can be estimated on this collection.
Effective lineage size for the between-lineage variance: 21.1.


## Per-agent feasibility

Each row uses only isolates with both a readable result for that agent and a recorded lineage. Support is the share of them in lineages of at least two isolates; the share is scored on those isolates and the singletons are set aside. The rarer-outcome count is a warning threshold, not a reporting condition of the estimator.

| Agent | Retained | Support | Input failures | Rarer outcome |
| --- | --- | --- | --- | --- |
| amoxicillin | 677 | 99.6 % | none at input level | 37 |
| cefquinome | 677 | 99.6 % | none at input level | 26 |
| ceftiofur | 677 | 99.6 % | none at input level | 163 |
| doxycycline | 677 | 99.6 % | none at input level | 106 |
| erythromycin | 677 | 99.6 % | none at input level | 311 |
| lincomycin | 677 | 99.6 % | none at input level | 253 |
| penicillin | 677 | 99.6 % | none at input level | 156 |
| spectinomycin | 677 | 99.6 % | none at input level | 78 |
| tetracycline | 677 | 99.6 % | none at input level | 103 |
| tiamulin | 677 | 99.6 % | none at input level | 131 |
| tilmicosin | 677 | 99.6 % | none at input level | 313 |
| trimethoprim | 677 | 99.6 % | none at input level | 190 |
| tylosin | 677 | 99.6 % | none at input level | 311 |