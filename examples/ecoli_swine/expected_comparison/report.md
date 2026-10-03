# Matched-record lineage comparison

Descriptive, not causal. The middle difference changes the labels on the same scored isolates, the records that sit in a repeated lineage under both definitions. The terms beside it are the isolates each definition sets aside as singletons, and the outer terms change record inclusion. The difference on all common records is the sum of the three middle terms. Common records need not represent either original collection or a population.

Status: computed

Records supplied: 481; labelled by the first definition: 472, by the second: 481; common observed records: 472, of which scored under both definitions: 445. Seed 42, 999 permutations, 999 bootstrap draws.

| Definition | Records | n | Lineage share | 95% interval | p |
| --- | --- | ---: | ---: | ---: | ---: |
| sequence_type | all records the first definition labels | 472 | 0.196 | 0.115 to 0.280 | 0.001 |
| sequence_type | records both definitions label | 472 | 0.196 | 0.115 to 0.280 | 0.001 |
| sequence_type | records scored under both definitions | 445 | 0.196 | 0.115 to 0.280 | 0.001 |
| phylogroup | records scored under both definitions | 445 | 0.012 | 0.000 to 0.047 | 0.029 |
| phylogroup | records both definitions label | 472 | 0.010 | 0.000 to 0.041 | 0.076 |
| phylogroup | all records the second definition labels | 481 | 0.013 | 0.000 to 0.045 | 0.050 |

| Contrast | Difference in the lineage share | 95% interval |
| --- | ---: | ---: |
| Record selection into the shared frame | +0.0000 | +0.0000 to +0.0000 |
| Isolates set aside as singletons, first definition | +0.0000 | +0.0000 to +0.0000 |
| Relabelling on the common scored isolates | -0.1842 | -0.2713 to -0.1033 |
| Isolates set aside as singletons, second definition | -0.0026 | -0.0115 to +0.0065 |
| Record selection out of the shared frame | +0.0036 | -0.0019 to +0.0078 |
| Change of definition on the shared records | -0.1868 | -0.2733 to -0.1067 |
| Total difference | -0.1832 | -0.2695 to -0.1034 |
| Residual of the identity | +0.0000 |  |

The first five contrasts sum to the total difference, and the three middle ones to the change of definition on the shared records; the residual records that identity. Interval for a contrast: 95% bootstrap-t intervals, records and labels fixed, outcomes drawn again from the smoothed distribution of the records sharing both labels.

Read `arms.csv` with the eligibility flags; `comparison.json` records the six analyses, split settings and descriptive differences.
