# 3 Check the input

```bash
amr-clonalshare --config examples/workflows/calls.yaml --check-input
```

Review the phenotype definition and duplicate policy first. Then inspect the raw observation frame, metadata-only identifiers, outcome and label missingness, empty agents, lineage sizes and exclusions. A successful input check means the files can be interpreted under the configuration. It does not mean every method is supported or that the collection is representative.

The raw frame includes call-only and MIC-only identifiers. Denominators later depend on the method and agent. Blank calls, unrecognized vocabulary and records absent for an agent are different reasons for missing data. Identical duplicate rows collapse; conflicting keys require resolution or explicit exclusion. Keep the QC JSON and Markdown with the run.

Attribution reports support as the share of retained isolates in lineages with at least two members. Other modules have their own requirements. Small datasets can retain useful prevalence and missing-outcome descriptions while refusing an interval, contrast or population fit. Read the method-specific reason rather than substituting zero or declaring the whole dataset unusable.
