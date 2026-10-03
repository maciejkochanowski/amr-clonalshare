# Matched-record lineage comparison

Descriptive, not causal. The middle difference changes the labels on the same scored isolates, the records that sit in a repeated lineage under both definitions. The terms beside it are the isolates each definition sets aside as singletons, and the outer terms change record inclusion. The difference on all common records is the sum of the three middle terms. Common records need not represent either original collection or a population.

Status: computed

Records supplied: 481; labelled by the first definition: 472, by the second: 481; common observed records: 472, of which scored under both definitions: 445. Seed 42, 999 permutations, 999 bootstrap draws.

| Definition | Records | n | Lineage share | 95% interval | p |
| --- | --- | ---: | ---: | ---: | ---: |
| sequence_type | all records the first definition labels | 472 | 0.196 | 0.115 to 0.280 | 0.001 |
| sequence_type | records both definitions label | 472 | 0.196 | 0.115 to 0.280 | 0.001 |
| sequence_type | records scored under both definitions | 445 | 0.196 | 0.115 to 0.280 | 0.001 |
| phylogroup_clermont2013 | records scored under both definitions | 445 | 0.055 | 0.011 to 0.110 | 0.001 |
| phylogroup_clermont2013 | records both definitions label | 472 | 0.053 | 0.011 to 0.113 | 0.001 |
| phylogroup_clermont2013 | all records the second definition labels | 481 | 0.062 | 0.019 to 0.122 | 0.001 |

| Contrast | Difference in the lineage share | 95% interval |
| --- | ---: | ---: |
| Record selection into the shared frame | +0.0000 | +0.0000 to +0.0000 |
| Isolates set aside as singletons, first definition | +0.0000 | +0.0000 to +0.0000 |
| Relabelling on the common scored isolates | -0.1417 | -0.2263 to -0.0624 |
| Isolates set aside as singletons, second definition | -0.0016 | -0.0144 to +0.0103 |
| Record selection out of the shared frame | +0.0089 | +0.0034 to +0.0148 |
| Change of definition on the shared records | -0.1433 | -0.2255 to -0.0595 |
| Total difference | -0.1344 | -0.2171 to -0.0499 |
| Residual of the identity | +0.0000 |  |

The first five contrasts sum to the total difference, and the three middle ones to the change of definition on the shared records; the residual records that identity. Interval for a contrast: 95% bootstrap-t intervals, records and labels fixed, outcomes drawn again from the smoothed distribution of the records sharing both labels.

Read `arms.csv` with the eligibility flags; `comparison.json` records the six analyses, split settings and descriptive differences.
