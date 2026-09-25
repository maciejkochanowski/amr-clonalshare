# 6 Interpret a lineage share

The lineage label, retained isolates and phenotype definition are part of the quantity. `clonal_share` reports an observed-scale lineage-membership share based on held-out predictions, permutation correction and lineage resampling. A permuted-label control helps characterize this estimator; it does not establish a biological mechanism.

The design-corrected realised component ratio has a different conditional target and Gaussian-model interval. The optional population model has a separate latent-liability ICC target. Compare uncertainty only for like targets. A lineage-bootstrap interval and a population interval cannot be ranked solely by width.

The attribution support fraction is the share of retained isolates in lineages with at least two members. Below 0.80, the output is flagged as unsupported. Fewer than ten lineages and dominant-lineage concentration also require attention. A wider population interval is not a general remedy for a small or selective collection. Historical validation contains substantial failures under sparse, nonnormal and latent-target settings.

Read the method's status, retained denominator and reasons before the point estimate. Missing or withheld limits are not zero limits. Lineage association does not identify transmission, causality, resistance mechanisms or clinical utility. See [methods](../methodology.md) for distinct targets and preserved negative findings.
