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

The two levels of `contrast_column` are the two collections. Every trait in
the layers is decomposed.

## What is reported, per trait

| field | meaning |
|---|---|
| `difference` | prevalence in the second collection minus the first |
| `composition` | the part due to a change in lineage mix, with its interval |
| `within_lineage` | the part due to a change in rate inside lineages, with its interval |
| `status` | `ok`, or a refusal with its reason |

The report flags **offsetting** traits: both components significant, opposite
in sign, each larger than the difference they produce. A prevalence table
calls such a trait stable. It is not.

## The gates

**Shared support.** The decomposition compares rates inside lineages present
in both collections. When fewer than `min_shared_support` of the isolates sit
in shared lineages, the within-lineage term is not identified and the trait is
refused rather than reported on the lineages that happen to overlap.

**Informatively missing lineage labels.** If the untyped share differs between
the two collections and untyped isolates carry a different rate, the contrast
is refused: the decomposition would then measure the typing process. The
shipped *S. suis* configuration `config_contrast.yaml` demonstrates this
refusal on MLST labels that are present for 97 % of one period and 40 % of the
other.

**False-discovery control across the panel.** Benjamini-Yekutieli, valid
under arbitrary dependence between antimicrobials, is the reported control;
the count under the independence assumption is given beside it.

## Prevalence on two scales

Every trait is also reported as prevalence per isolate and prevalence per
lineage. When the two differ widely, a few large lineages carry most of the
resistance, and the per-isolate figure is a statement about sampling as much
as about the population.
