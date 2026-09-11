# 7. Mix or rate: two collections

A prevalence of 30 % tetracycline resistance last year and 40 % this year can
mean two different things: the lineages in the collection changed (a lineage
that carries the trait became more common), or resistance rose inside the
lineages already there. The first is a change in **mix**, the second a change
in **rate**, and they call for different responses. The Kitagawa
decomposition splits the observed difference into the two.

## Configuration

```yaml
dataset:
  metadata: metadata.csv
  lineage_column: lineage
  contrast_column: period
  contrast_levels: [2019, 2023]
surveillance:
  enabled: true
  n_boot: 400
  min_shared_support: 0.8
```

The two levels of `contrast_column` are the two collections. Every
antimicrobial in the panel is decomposed.

## What is reported, per trait

| field | meaning |
|---|---|
| `difference` | prevalence in the second collection minus the first |
| `composition` | the part due to a change in lineage mix, with its interval |
| `composition_shared` | the part of `composition` carried by lineages present in both collections |
| `nonshared` | the part of `composition` carried by lineages present in one collection only: the whole contribution of a lineage the other collection never held, reported apart so that a mix finding can be read for how much of it rests on such lineages |
| `within_lineage` | the part due to a change in rate inside lineages, with its interval |
| `status` | `ok`, or a refusal with its reason |

`composition` equals `composition_shared + nonshared`, and `difference` equals
`composition + within_lineage`. A lineage seen in one collection only is a
statement about the mix, not about a rate, which is why its contribution sits
in `composition`; a lineage first seen in the later collection may have entered
it, or may have been present and unsampled before, and the number does not say
which.

The report counts discoveries after the step-up, never nominal intervals,
and a component the gate refused is counted as refused whatever its interval
says: the within-lineage component where shared support falls short
(`within_lineage_estimable`), both components where the lineage labels are
informatively missing (`composition_estimable`). It flags **offsetting** traits: both components
discoveries, opposite in sign, each larger than the difference they produce.
A prevalence table shows such a trait as unchanged although both components
moved.

## The gates

**Shared support.** The decomposition compares rates inside lineages present
in both collections. When fewer than `min_shared_support` of the isolates sit
in shared lineages, the within-lineage term is not identified and the trait is
refused instead of being reported on the lineages that happen to overlap.

**Informatively missing lineage labels.** If the untyped share differs between
the two collections and untyped isolates carry a different rate, the contrast
is refused for both components: the decomposition would then measure the
typing process. On the
shipped *S. suis* cohort the period contrast at sequence-type resolution
(`python examples/ssuis/decompose_trend.py --contrast period --lineage mlst`)
meets both gates at once: the type is present for 97 % of one period and
40 % of the other, and the isolates the two periods share in a typed lineage
fall below the support the within-lineage term needs, so the script prints
the refusal and the reasons beside the table. `config_contrast.yaml`, the
country contrast at the same resolution, is refused by the shared-support
gate alone.

**False-discovery control across the panel.** Benjamini-Yekutieli, valid
under arbitrary dependence between antimicrobials, is the reported control;
the count under the independence assumption is given beside it.

## Prevalence on two scales

Every trait is also reported as prevalence per isolate and prevalence per
lineage. When the two differ widely, a few large lineages carry most of the
resistance, and the per-isolate figure is a statement about sampling as much
as about the population.
