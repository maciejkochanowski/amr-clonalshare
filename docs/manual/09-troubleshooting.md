# 9. Exit codes and troubleshooting

| code | meaning | what to do |
|---|---|---|
| `0` | the run completed, including a run whose estimate was withheld | read `report.md`, whose first page names every condition that failed |
| `2` | configuration or input error | the message on standard error names the key, file, row or column; nothing was computed |

## Messages you may see at exit code 2

**`unknown key 'lineage_colum' in section 'dataset'; did you mean 'lineage_column'?`**
The configuration is a contract; a key the loader does not read is refused,
not ignored, because an ignored key is a silently disabled analysis.

**`phenotype table ... has no column 'call'`**
Name the columns in `dataset.phenotype_id_column`,
`phenotype_antibiotic_column`, `phenotype_call_column`. The call column must
hold `susceptible`, `non-susceptible`, `resistant` or `intermediate`; a row
with any other word is not read, and an agent no isolate carries a readable
call for is absent from the run.

**`metadata ... has no column 'ST'`**
The lineage column is named in `dataset.lineage_column` and must be a column
of the metadata table.

**`MIC table ... joined none of the 677 aligned strains`**
The identifiers of the dilution table do not match those of the phenotype
table. Compare them as text; a trailing space or a different case is a
different isolate. Nothing is repaired here, because a loader that repaired
identifiers could join the wrong isolates without saying so.

**`MIC table ... has no column 'measurement'`**
Name the columns in `dataset.mic_id_column`, `mic_antibiotic_column`,
`mic_value_column`.

## Readings that look wrong and are not

**Every clonal share is "not estimable".** The lineage definition is too fine
for this collection: too many isolates sit alone in their lineage. The input
check reports the support; a serovar or a coarser cluster level raises it.

**The clonal share of a trait is slightly negative.** The out-of-sample
estimate is debiased, and a trait the lineage carries no information about
scatters around zero. Read the interval.

**The realised interval is withheld and the kurtosis is printed.** The exact
interval rests on the within-lineage residuals having tails no heavier than
the Gaussian, and the gate protects by refusing where they do; the share and
its lineage bootstrap are still reported. A binary call at low prevalence is
the usual case.

**The e-value is 1 and the share is not zero.** The e-value asks whether the
lineage carries any information about the trait at this look, scored on
held-out folds; a share whose interval crosses zero can sit beside an e-value
of one, and both are saying the same thing.

## Getting help

Open an issue in the repository with the configuration, `input_qc.md`, and
the message on standard error. Do not attach isolate-level data unless it is
already public.
