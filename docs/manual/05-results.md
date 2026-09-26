# 5 Read the results

Start with phenotype meaning, retained cohort and method status. A numerical estimate, a supported estimate and an interval admitted by its diagnostic are distinct outputs. A withheld or incomplete result is not a zero effect.

The canonical file is `clonal_share_result.json`. Schema 2.0 adds explicit collection bounds and module results. CSV and Markdown/HTML reports are derived from this record; they are not independent analyses. QC files explain joins, duplicate handling, missing outcomes and labels, and per-agent exclusions. The completion manifest identifies a finished bundle and records file identities.

Both reports show recorded reasons and next checks for results requiring attention. The status vocabulary is unchanged:

| Status | Meaning and next check |
|---|---|
| `computed` | A reportable value exists. Read its target, interval kind, assumptions and retained subset. |
| `full_range` | A completed result spans [0, 1] and is uninformative. For collection bounds, inspect missing outcomes; for a model interval, inspect repeated-group support and outcome variation. |
| `unavailable` | The method has no reportable result under its recorded conditions. Read the reason and inspect relevant input QC, support or method domain. |
| `computation_incomplete` | The calculation did not finish. Inspect numerical/resource diagnostics before retrying in a new directory. |
| `empty_confidence_set` | A completed calculation retained no parameter value. Inspect numerical diagnostics and model compatibility; this is not a zero estimate. |

These checks explain the result; they do not authorize relaxing a support gate or changing the cohort to obtain a preferred answer.

| Result family | What to read |
|---|---|
| Attribution | Observed-scale lineage-membership share, bootstrap, permutation control and support |
| Realised component | Conditional variance-component target and model-specific interval status |
| MIC | Dilution-scale component with its calibrated 95% interval (`calibrated`), the moment estimate with its approximate F interval as a comparison (`share`), panel geometry and the dilution table |
| Strata | With `stratify_by`, the same families repeated per level of the column (`strata` in the record, Table 5d in the report, `strata_results.csv`); each stratum stands on its own cohort and is not adjusted for the others |
| Contrast | Typed-subset scope, component terms, shared support and missingness diagnostics |
| Population probit | Gaussian latent-liability ICC target, explicit interval method, computation status |
| Collection bounds | Recorded-frame denominator and exact range under missing binary outcomes |
| Evidence | Fixed-look or sequential construction and its null-model interpretation |

For a MIC table the report adds, after the share table, the panel each laboratory tested with its cut points (Table 5a), the fixed-effect coefficients of every declared covariate with their reference level (Table 5b), the status and reason of every calibrated interval with the simulation budget it was given (Table 5c), and a heat map of the readings per lineage and dilution interval. Read the heat map first: a lineage whose readings fill one or two adjacent cells carries the lineage share, and a row spread along the panel does not.

## Missing outcomes and labels

If there are `r` recorded positive outcomes, `m` missing outcomes and `N` isolates in the stated frame, finite-collection prevalence lies between `r/N` and `(r+m)/N`. These exact bounds allow every missing outcome to be zero or one. They contain no sampling confidence level, do not impute a value and say nothing about isolates outside the frame. The frame and counts must accompany the bounds. The prevalence point estimate describes observed outcomes only; it is not a point estimate for the full frame when outcomes are missing.

A supported decomposition with missing lineage labels describes the typed subset. `collection_generalization_supported` remains false when labels are missing. Missingness-association tests are descriptive diagnostics; equal typing fractions or a nonsignificant test cannot show that untyped isolates have the same lineage distribution.

The result reader accepts records of schema 1.0 and refuses any other schema. Keep original records and manifests with any downstream reader.

## Compare compatible quantities

Compare like targets, phenotype definitions, units, cohorts and lineage resolutions. The lineage-membership share, conditional component ratio, dilution-scale component and Gaussian liability ICC are not interchangeable. Narrower intervals for different targets do not demonstrate a better estimator. No result alone establishes transmission, a mechanism, an intervention effect or clinical validity.
