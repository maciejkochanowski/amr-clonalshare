# 5 Read the results

Start with phenotype meaning, retained collection and method status. A numerical estimate, an interval and a bound are distinct outputs, each with its own target. A result that was not computed is not a zero effect.

The canonical file is `clonal_share_result.json`, a record of schema 1.0 carrying the configuration as run, the software versions, the input check, the collection bounds and the method results. CSV and Markdown/HTML reports are derived from this record; they are not independent analyses. QC files explain joins, duplicate handling, missing outcomes and labels, and per-agent exclusions. The completion manifest identifies a finished bundle and records file identities.

Both reports show recorded reasons and next checks for results requiring attention. Every row of `results.csv` carries one of three statuses:

| Status | Meaning and next check |
|---|---|
| `computed` | A reportable value exists. Read its target, interval kind, assumptions and retained subset. |
| `full_range` | A completed interval or pair of bounds spans [0, 1] and is uninformative. For collection bounds, inspect missing outcomes; otherwise inspect repeated-lineage support and outcome variation. |
| `unavailable` | The method has no reportable result under its recorded conditions. Read the reason and inspect relevant input QC, support or method domain. |

These checks explain the result; they do not authorize relaxing a reporting condition or changing the collection to obtain a preferred answer.

| Result family | What to read |
|---|---|
| Lineage share of the call | `kappa_adj` with its interval for the represented lineages (`observed_low`, `observed_high`, `observed_se`), the permuted control (`null_mean`), the permutation p-value with its floor, support and the lineage structure; the selection across agents (`lineage_selection`) |
| Latent ordering | Beside every call and every MIC ordering, the bounds on the lineage share of the latent ordering (`latent_order_lower`, `latent_order_upper_bound`, with `latent_order_upper` the largest share attained and `latent_order_upper_exact` whether the two agree), the one-sided 95% lower confidence limit (`latent_order_lower_limit`), `latent_order_sharp`, `panel_resolution` and, with several strata, the lower bound within each |
| MIC | Lineage share of the MIC ordering with its interval for the represented lineages, permutation p-value, strata, distinct readings and share of end-well readings (`order`); the selection across agents (`order_selection`); panel geometry and the dilution table |
| Strata | With `stratify_by`, the same families repeated per level of the column (`strata` in the record, Table 5d in the report, `strata_results.csv`); each stratum stands on its own collection and is not adjusted for the others |
| Contrast | Typed-subset scope, component terms, shared support and missingness diagnostics |
| Collection bounds | Recorded-frame denominator and exact range under missing binary outcomes |
| Evidence | Fixed-look or sequential construction and its null-model interpretation |

In `results.csv` every analysis is a row per agent: `collection_membership` for the share of the call and `mic_order` for the share of the MIC ordering, with the interval kind `95% interval for the represented lineages (isolates drawn again, lineages fixed)`; `call_latent_bounds` and `mic_latent_bounds` for the bounds, with no point estimate and the interval kind `bounds the readings allow; not a confidence interval; lower end established by the readings`, followed by `upper end sharp` or `upper end a certified bound, not sharp`; `call_latent_lower_limit` and `mic_latent_lower_limit`, with the interval kind `one-sided 95% lower confidence limit for the lower bound, and so for the share`; and the rows of the collection bounds (`finite_collection_prevalence`), the e-values (`lineage_evidence`, `sequential_lineage_evidence`) and the decomposition (`decomposition_composition`, `decomposition_within_lineage`). The lower limit is the share of the latent ordering that the lineages of the collection are shown to account for. Every row also carries its permutation `p_value`, the Benjamini-Yekutieli `q_value` across agents where one is computed, and `selected` (yes or no) where the row enters a selection. The share itself can be negative near the permutation null; its interval is reported within 0 and 1.

The rows of the two shares also carry a `conclusion`: one sentence that reads the share, its interval, the p- and q-values and the bounds together. *Lineage structure established* means the interval lies above zero with a q-value (or, without a selection, a permutation p-value) at or below the level of the selection, or the lower confidence limit of the latent share is above zero. *Not established* is a computed share that meets neither condition; *interval spans the whole range* a share the collection cannot resolve; *not estimable* a share the collection could not give, with the reason. The same conclusion stands in the last column of the results table of the local form, of the console of a command-line run and of Tables 1 and 5 of the report. It is a statement about the collection analysed, not about transmission or a resistance mechanism. A permutation p-value equal to the smallest value the permutations can give, 1/(permutations + 1), is written with ≤ in the reports, since no permutation reached the observed statistic; `results.csv` keeps the number.

The report opens with a measurement summary (Tables 0 to 1b: the share of every call with its interval, p-value, Benjamini-Yekutieli q-value, selection and reading, and the lineage structure behind it), the input conditions (Tables 2 and 2a) and the share of every call with the bounds it places on the latent ordering (Table 3), followed by the evidence (Tables 4 and 4b). For a MIC table it gives the lineage share of the MIC ordering per agent with its interval, p-value, q-value, selection and bounds (Table 5, with the mean share of the permuted labels as the control), the lineage structure behind it (Table 5b), the panel each laboratory tested (Table 5a) and a heat map of the readings per lineage and dilution interval. Read the heat map first: a lineage whose readings fill one or two adjacent cells accounts for the lineage share, and a row spread along the panel does not. Table 5d gives the shares within each stratum, Table 6 the decomposition of a prevalence difference, later level minus earlier as the configuration orders them, and Table 6b the prevalence of each collection with the lineages they share.

## Missing outcomes and labels

If there are `r` recorded positive outcomes, `m` missing outcomes and `N` isolates in the stated frame, finite-collection prevalence lies between `r/N` and `(r+m)/N`. These exact bounds allow every missing outcome to be zero or one. They contain no sampling confidence level, do not impute a value and say nothing about isolates outside the frame. The frame and counts must accompany the bounds. The prevalence point estimate describes observed outcomes only; it is not a point estimate for the full frame when outcomes are missing.

A supported decomposition with missing lineage labels describes the typed subset. `collection_generalization_supported` remains false when labels are missing. Missingness-association tests are descriptive diagnostics; equal typing fractions or a nonsignificant test cannot show that untyped isolates have the same lineage distribution.

The result reader accepts records of schema 1.0 and refuses any other schema. Keep original records and manifests with any downstream reader.

## Compare compatible quantities

Compare like targets, phenotype definitions, units, collections and lineage resolutions. The share of a call, the share of the MIC ordering, the bounds on the share of the latent ordering and the components of a prevalence difference are not interchangeable: the first two are estimated for the represented lineages from what was recorded, the bounds are what the readings allow for the latent ordering, and the components describe two collections. Narrower intervals for different targets do not demonstrate a better estimator. No result alone establishes transmission, a mechanism, an intervention effect or clinical validity.
