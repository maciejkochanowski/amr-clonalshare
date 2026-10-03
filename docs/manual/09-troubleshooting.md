# 9 Exit codes and troubleshooting

| Code | Meaning | Next step |
|---|---|---|
| `0` | Completed run, including methods that are unavailable on this collection | Inspect the completion manifest, canonical method statuses and QC |
| `2` | Configuration, input or output error | Read the error; correct its cause and use a fresh destination |

`amr-clonalshare` and `amr-clonalshare-compare` use the same two codes. A configuration may be valid while a particular estimator is unsupported. A directory without a completion manifest is not a finished bundle, even if an intermediate file exists. The results folder of every run holds `run.log`, the progress lines and warnings of the run with the time of each; a run that failed leaves `run.log` with its error as the only file of the folder. When asking for help, send that file with the configuration. `amr-clonalshare doctor` prints the versions of Python and the numerical stack, the processor threads, whether workbooks can be read and where the examples are.

## Input and output errors

Unknown configuration keys are rejected, with the nearest known key named. Check spelling and the current input manual. A missing-column error means the configured column name does not occur in its source table; changing the name does not change the phenotype vocabulary. Declare `clinical_sir`, `wt_nwt` or `binary`, and inspect unrecognized strings. Empty agents remain in QC instead of silently becoming negative outcomes.

Conflicting records are rejected by default. Resolve the source discrepancy or explicitly select a documented conflict policy. Identical records collapse. Missing identifiers are input errors; the software does not repair names. Call-only and MIC-only identifiers belong to the raw union frame, so a lack of overlap with calls is not by itself a reason to remove an MIC record. A MIC reading without a value in the panel column or in a declared covariate column is refused.

An existing output destination is protected. Use a new directory or an intentional `--overwrite`, which replaces the earlier results files; copy the earlier manifest first when comparing runs.

## Numerical statuses

An unavailable share may reflect insufficient repeated-lineage support, few usable outcomes, a constant outcome or another method-specific condition. Read the recorded reason. A coarser label changes the target; it is not an automatic fix that preserves the original question.

A slightly negative share can arise near the permutation null; interpret it with its interval and support. An interval that is not reported with a computed share means that fewer than 90% of the bootstrap draws gave a share; the reason says how many did. An interval of [0, 0] lies wholly below zero before it is bounded. A p-value at `p_floor` is the smallest `n_perm` permutations can return; raise `n_perm` rather than read it as exact.

Bounds on the share of the latent ordering whose upper end is marked as a certified bound, not sharp, did not meet the largest share attained; the bound still holds. Bounds that are not sharp (`latent_order_sharp` false) come from a stratum whose readings overlap, as when two panels share a stratum; where the panels are known, `mic_panel_column` reads each within a stratum of its own.

An e-value of one and a nonzero share are not contradictory. Evidence and magnitude are different outputs with different constructions. Sequential validity requires actual batch ordering and the stated conditional null model.

## Report an issue

Include the software version, configuration, error text and a small non-sensitive reproducer. Share QC and manifest information when appropriate. The repository support contact is recorded in the code metadata.
