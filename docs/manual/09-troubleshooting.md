# 9 Exit codes and troubleshooting

| Code | Meaning | Next step |
|---|---|---|
| `0` | Completed run, including withheld methods | Inspect the completion manifest, canonical method statuses and QC |
| `2` | Configuration, input or output error | Read the error; correct its cause and use a fresh destination |
| `3` | `amr-clonalshare-mic` only: the computation completed but the interval is not reportable (unresolved inversion, full range) | Read the status in the output file; raise `--bootstrap` or inspect the tested values |

A configuration may be valid while a particular estimator is unsupported. A directory without a completion manifest is not a finished bundle, even if an intermediate file exists.

## Input and output errors

Unknown configuration keys are rejected. Check spelling and the current input manual. A missing-column error means the configured column name does not occur in its source table; changing the name does not change the phenotype vocabulary. Declare `clinical_sir`, `wt_nwt`, `binary` or explicit `legacy`, and inspect unrecognized strings. Empty agents remain in QC instead of silently becoming negative outcomes.

Conflicting records are rejected by default. Resolve the source discrepancy or explicitly select a documented conflict policy. Identical records collapse. Missing identifiers are input errors; the software does not repair names. Call-only and MIC-only identifiers belong to the raw union frame, so a lack of overlap with calls is not by itself a reason to remove an MIC record.

An existing output destination is protected. Use a new directory or intentional `--overwrite`, which preserves a backup. Retain the earlier manifest when comparing runs.

## Numerical statuses

An unsupported attribution may reflect insufficient repeated-lineage support, few usable outcomes or another method-specific condition. Read the recorded reason. A coarser label changes the target; it is not an automatic fix that preserves the original question.

A slightly negative debiased attribution can arise near the permutation null; interpret it with its uncertainty and support. A withheld realised interval means its reporting criterion was not met. Passing the kurtosis criterion does not prove Gaussian errors or ensure nominal coverage, and historical stress tests retain serious failures.

An e-value of one and a nonzero share are not contradictory. Evidence and magnitude are different outputs with different constructions. Sequential validity requires actual batch ordering and the stated conditional null model.

For the optional general population route, a full range, a one-sided upper object, no retained data and incomplete computation are separate states. Do not change missing limits to zero, or interpret a full interval as an estimated absence of association. See the general-method resource and status documentation.

## Report an issue

Include the software version, configuration, error text and a small non-sensitive reproducer. Share QC and manifest information when appropriate. The repository support contact is recorded in the code metadata.
