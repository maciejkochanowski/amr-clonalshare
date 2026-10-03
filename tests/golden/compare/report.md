# Matched-record lineage comparison

Descriptive, not causal. The middle difference changes the labels on the same scored isolates, the records that sit in a repeated lineage under both definitions. The terms beside it are the isolates each definition sets aside as singletons, and the outer terms change record inclusion. The difference on all common records is the sum of the three middle terms. Common records need not represent either original collection or a population.

Status: computed

Records supplied: 96; labelled by the first definition: 96, by the second: 88; common observed records: 88, of which scored under both definitions: 88. Seed 42, 999 permutations, 999 bootstrap draws.

| Definition | Records | n | Lineage share | 95% interval | p |
| --- | --- | ---: | ---: | ---: | ---: |
| lineage_a | all records the first definition labels | 96 | -0.084 | 0.000 to 0.014 | 0.995 |
| lineage_a | records both definitions label | 88 | -0.096 | 0.000 to 0.015 | 0.993 |
| lineage_a | records scored under both definitions | 88 | -0.096 | 0.000 to 0.015 | 0.993 |
| lineage_b | records scored under both definitions | 88 | -0.047 | 0.000 to 0.035 | 0.857 |
| lineage_b | records both definitions label | 88 | -0.047 | 0.000 to 0.035 | 0.857 |
| lineage_b | all records the second definition labels | 88 | -0.047 | 0.000 to 0.035 | 0.857 |

| Contrast | Difference in the lineage share | 95% interval |
| --- | ---: | ---: |
| Record selection into the shared frame | -0.0118 | -0.0427 to +0.0086 |
| Isolates set aside as singletons, first definition | +0.0000 | +0.0000 to +0.0000 |
| Relabelling on the common scored isolates | +0.0491 | -0.0195 to +0.0863 |
| Isolates set aside as singletons, second definition | +0.0000 | +0.0000 to +0.0000 |
| Record selection out of the shared frame | +0.0000 | +0.0000 to +0.0000 |
| Change of definition on the shared records | +0.0491 | -0.0195 to +0.0863 |
| Total difference | +0.0373 | -0.0248 to +0.0722 |
| Residual of the identity | +0.0000 |  |

The first five contrasts sum to the total difference, and the three middle ones to the change of definition on the shared records; the residual records that identity. Interval for a contrast: 95% bootstrap-t intervals, records and labels fixed, outcomes drawn again from the smoothed distribution of the records sharing both labels.

Read `arms.csv` with the eligibility flags; `comparison.json` records the six analyses, split settings and descriptive differences.
