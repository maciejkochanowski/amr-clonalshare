# 2 Prepare the input

Use a metadata CSV with one isolate identifier and its lineage label, and a long-format calls CSV, a long-format MIC CSV, or both. Identifier spelling must match exactly. The recorded observation frame is the union of identifiers in the calls and MIC tables. Metadata-only identifiers are inventoried; they are not automatically sampled isolates. Each analytical module then retains the outcomes and labels it needs, so MIC-only isolates are not discarded because they lack calls.

```yaml
dataset:
  name: laboratory_collection
  data_dir: data
  metadata: metadata.csv
  strain_id_column: isolate_id
  lineage_column: lineage
  phenotype: calls.csv
  phenotype_id_column: isolate_id
  phenotype_antibiotic_column: agent
  phenotype_call_column: call
  phenotype_kind: clinical_sir
  phenotype_positive_definition: R only
  phenotype_source: laboratory export
  duplicate_policy: error
```

`name` labels the collection; it is written into the record and the reports and enters the run identifier, and is read nowhere else. `data_dir` is relative to the configuration file. The recipe directory supplies complete examples. Record the label definition, clustering pipeline/version or serotyping vocabulary with the data. The software cannot recover that provenance from a short label.

## Declare the phenotype

| `phenotype_kind` | Intended input | Positive outcome |
|---|---|---|
| `clinical_sir` | Clinical susceptibility S/I/R | R only by default; state any explicit intermediate policy |
| `wt_nwt` | Wild type/non-wild type | NWT; do not relabel this clinical resistance |
| `binary` | Recorded 0/1 | 1; describe the meaning in `phenotype_positive_definition` |
| `legacy` | Historical package vocabulary | Exact historical interpretation, labelled as legacy |

Use `phenotype_source`, `ast_standard` and `ast_version` when applicable and known. Do not invent a standard for cohort-derived thresholds. `phenotype_positive_definition` documents the meaning; it does not calculate clinical breakpoints. For clinical data, `phenotype_intermediate: non_susceptible` explicitly combines I and R; `susceptible` maps I to zero and `drop` excludes I. The definition and policy must agree. WT/NWT and binary configurations do not accept an intermediate policy.

Blank calls, unrecognized strings, absent agent records, untyped isolates and empty agents are counted separately in QC. None is automatically a negative outcome. Inspect these counts before interpreting a prevalence denominator.

## Resolve duplicates explicitly

Identical repeated records collapse. Conflicting records raise an error by default (`duplicate_policy: error`). `drop_conflicts` excludes conflicting keys and records the exclusions. `legacy` is for reproduction of historical handling. Do not select a permissive mode merely to remove an error: resolve the source discrepancy or state why exclusion is justified. These checks apply to the relevant metadata, call and MIC keys.

## Preserve MIC measurement information

```yaml
  mic: mic.csv
  mic_id_column: isolate_id
  mic_antibiotic_column: agent
  mic_value_column: measurement
  mic_operator_column: operator
  mic_panel_column: laboratory
  mic_covariate_columns: [country, year]
  stratify_by: country
  mic_wells:
    demo_agent: [0.25, 0.5, 1, 2, 4, 8]
  mic_units:
    demo_agent: mg/L
```

Place these keys inside `dataset`. The well list is agent-specific, strictly increasing, positive and finite. The operator column can record `<`, `<=`, `>` or `>=`; keep the numerical concentration in the measurement column. Explicit operators take precedence over end-well assumptions. Units are recorded, not automatically converted. Supply actual panel wells; wells inferred from observed concentrations can omit tested but unobserved concentrations. Check the reported geometry and sensitivity result.

`mic_panel_column` names the panel a reading was made on: the laboratory, the panel product, or any label under which every isolate was tested on the same wells. It may be a column of the MIC table or, per isolate, of the metadata. Without it, the wells are inferred from the pooled readings of an agent, and when two laboratories tested different ranges the lowest well of one passes as an interior dilution of the other, so its left-censored readings are read as exact. With it, the wells and the end wells are inferred within each label, the geometry is reported per label, and the calibrated interval simulates every reading on the panel of the record it stands for. Every reading needs a label; `mic_wells`, when given, apply to every label alike.

`mic_covariate_column` names one categorical covariate of the MIC model, in the MIC table or, per isolate, in the metadata: the testing laboratory, the country, or another factor whose levels may shift MICs. Its levels enter the calibrated interval as fixed effects, with the first level in sorted order as the reference. The coefficients are estimated jointly with the mean, the scale and ρ, re-estimated under every candidate ρ and in every simulated dataset, so that a shift between levels is removed from both variance components rather than absorbed into whichever of them it is aligned with; the fitted coefficients are recorded with the interval. A level that follows lineage exactly is not separable from the lineage effects, and the model then rests on the isolates that break the alignment. The moment estimator and the end-well sensitivity do not use the covariate, and exact readings with a covariate use the bootstrap rather than the exact pivot.

`mic_covariate_columns` lists several covariates, for instance the country and the year of isolation. Each enters as its own set of fixed effects with the first level in sorted order as its reference; the offsets add, so the model is additive across covariates on the log2 scale and does not fit interactions between them. `mic_covariate_column` names one covariate and may be combined with the list. The report gives every coefficient with its covariate and reference level.

`stratify_by` names a metadata column whose levels define strata. The whole run is repeated on the isolates of each level, with every analysis on its own lineages and, for a MIC table, its own panel; the pooled run is unchanged and the strata are added to the record under `strata`, to the report as a summary table, and to `strata_results.csv`. Isolates without a level are left out of the strata and counted. Stratification answers whether a share holds within each level; the fixed effect of `mic_covariate_columns` answers whether a shift between levels explains it. The two can be used together.

Metadata may also supply `contrast_column` with exactly two `contrast_levels`, or a `batch_column` that sorts into true arrival order. A sequential analysis requires that ordering and the method's null assumptions; retrospective reordering changes the analysis.
