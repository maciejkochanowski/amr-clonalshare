# 5. Read the results

Start with `report.md`. It has the same seven sections on every run, in the
same order, and the first page is meant to be read on its own: the measurand,
the gates that fired, one row per trait, what was not evaluated on this run
and why, and an interpretation generated from that table by fixed rules.
Everything after the first page defends it. Both forms of the report,
`report.md` and `report.html`, are rendered from one structure, read every
value from `clonal_share_result.json`, and recompute nothing; the page carries
the record's SHA-256 so a report can be matched to the record it came from.
Open the record only when a specific number has to be traced.

![The run report on the shipped cohort](../img/05_report.png)

## The sections of the report

**1. Measurement summary.** The quantity measured, the cohort, the gates
that fired, and one row per trait: the share, its interval, the shuffled-label
control, the e-value and a one-phrase reading. Below the table: what was not
evaluated on this run, with the reason, and the interpretation the fixed rules
draw from the table. This is the page a control programme reads.

**2. Admissibility of the input.** The lineage column and its resolution, the
panel and the dilutions joined, the resampling budgets, the input check, the
lineage sizes drawn largest first with the singletons that count against
support flagged, and the conditions the estimator requires before a value is
reported, each with its observed value, its threshold and a verdict.

**3. The clonal share, trait by trait.** The intervals drawn, each with the
species interval beneath it, with the control and the calibration line: how often intervals of this estimator contained the
truth on the validation grid, so the stated 95 % can be read as a measured
figure. Then the two intervals side by side: for the share the lineages in
this collection carry, and for the share a fresh draw of lineages from the
species would show, which is wider when the collection holds few lineages;
for a binary trait the same table restates the share on the latent scale a
binomial mixed model reports.

**4. Evidence that survives re-reading.** The e-value per trait for this run,
drawn against the single-look threshold 1/α and the e-BH threshold of the
panel, and, when intakes are declared, the running product of every trait
over the intakes, drawn against the same two rules, so that a programme sees
at which intake a trait's evidence crossed the line and whether it stayed
there.

**5. Reading at the recorded resolution.** Present when the panel carries
dilution values: the share read from the dilutions, with its realised and its
species interval, and the share of the binary call, drawn side by side per
agent with the share of readings that sat at an end well.

**6. How resistance is carried across lineages.** Per-lineage against
per-isolate prevalence for every trait, in a table and drawn, beside the
effective number of carrying lineages against its permutation interval for
carriage in proportion to lineage size, with the direction of any departure;
and the decomposition into lineage composition and within-lineage change
when two collections were contrasted.

**7. Provenance and terms.** Software, seed, configuration and record
digests, and the terms used, defined once.

## clonal_share_result.json

The full record: the configuration as run (`config`), the versions of the
package and its numerical dependencies (`versions`), the input check, and
under `metadata_diagnostics` every estimate with its inputs: the per-trait
clonal shares (`clonal_share`, each with `kappa_adj`, `ci_low`, `ci_high`, the
species interval, the latent restatement, `null_mean`, `p_value`, `p_floor`,
`support` and `estimable`), the realised share per trait with its exact
interval and the verdict of its kurtosis gate (`realised_share`), the e-values (`lineage_evidence`, with
`sequential` when an intake column was declared), the censored-panel reading
(`censored_share`), prevalence per isolate and per lineage
(`lineage_resolved_prevalence`), the concentration of carriage
(`trait_concentration`) and the decomposition (`prevalence_decomposition`)
when two collections were contrasted. It is strict JSON: no `NaN`, no
`Infinity`, so any parser reads it. The command line prints a digest of the
same record on standard output: one row per antimicrobial with the share, its
interval, the control, the p-value, the support and the verdict, and the
agents e-BH selected.
