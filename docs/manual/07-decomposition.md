# 7 Compare two collections

The Kitagawa identity describes how a prevalence difference partitions into lineage composition and within-lineage rates. It does not identify causes of the difference. Set a metadata contrast column and exactly two levels; the reported difference is the first level minus the second, so list the later period first.

```yaml
dataset:
  contrast_column: period
  contrast_levels: [later, earlier]
surveillance:
  enabled: true
  min_shared_support: 0.9
```

This fragment belongs in a complete configuration; see the contrast recipe. `composition` equals `composition_shared + nonshared`, and `difference` equals `composition + within_lineage`. A lineage observed in one collection only contributes to the nonshared term. It may have been unsampled in the other collection; its appearance does not establish introduction.

Shared support, the mean of the two collections' shares of isolates in lineages both hold, controls whether the comparison has adequate overlap: below `min_shared_support` (0.9 by default) both components are marked not estimable. Read `composition_estimable`, `within_lineage_estimable`, scope and reasons separately. Missing labels do not erase a supported arithmetic description of the typed subset. They do prevent a claim that this description generalizes to the full collection: `collection_generalization_supported` remains false whenever labels are missing. Association tests for typing status and phenotype are diagnostics, not certificates of representativeness, even if typing fractions are equal or the tests are nonsignificant.

The limits of both components are bootstrap percentiles. Isolates are resampled within each collection, or whole sampling units where `unit_column` names them, since isolates of one farm or host need not be independent.

Panel discoveries use Benjamini–Yekutieli adjustment within each stated component family. Multiplicity control does not repair selective observation or turn the components into causal effects. Report both component signs, intervals, support and subset denominators. Opposite-signed components can offset each other despite a small total difference.

Counts from other runs were generated under their own interpretation and configuration. Rerun an analysis before presenting its counts as new results. The *S. suis* and poultry examples keep their source provenance.
