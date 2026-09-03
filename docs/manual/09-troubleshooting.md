# 9. Exit codes and troubleshooting

| code | meaning | what to do |
|---|---|---|
| `0` | the run completed and no diagnostic gate tripped | read `report.md` |
| `1` | `--validate` ran the planted-truth recovery and failed its threshold | the environment does not reproduce the shipped control; check the library versions against `requirements-lock.txt` |
| `2` | configuration or input error | the message on standard error names the key, file, row or column; nothing was computed |
| `3` | the analysis ran but a diagnostic gate tripped | read section 6 of `report.md`; the partition is not interpretable as reported and the record says why |

## Messages you may see at exit code 2

**`unknown key 'lineage_colum' in section 'dataset'; did you mean 'lineage_column'?`**
The configuration is a contract; a key the loader does not read is refused
rather than ignored, because an ignored key is a silently disabled analysis.

**`binary layers hold cells with no value; a 0 written there would read as absence of the trait`**
Set `dataset.missing_policy` to `drop_rows` or `drop_columns`, or repair the
file. The message lists example rows and columns.

**`layer 'amr' holds values other than 0 and 1: {'2': 3}`**
A count, a MIC, or a text flag has reached a binary layer. Recode, or move the
column to the MIC table.

**`layer 'amr' has duplicate strain IDs`**
De-duplicate before the run; the tool does not guess which row is right.

**`strict_n policy: aligned 640 strains, expected 677`**
Some identifiers do not match across layers. Compare the identifier columns
as text; a trailing space or a different case is a different isolate.

**`MIC table ... has no column 'measurement'`**
Name the columns in `dataset.mic_id_column`, `mic_antibiotic_column`,
`mic_value_column`.

## Readings that look wrong and are not

**Every clonal share is "not estimable".** The lineage definition is too fine
for this collection: too many isolates sit alone in their lineage. The input
check reports the support; a serovar or a coarser cluster level raises it.

**`selected_k` is 1.** The profiles do not separate into types. This is the
null result and it is reported as one; the lineage question in section 4 of
the report is still answered.

**The clonal share of a trait is slightly negative.** The out-of-sample
estimate is debiased, and a trait the lineage carries no information about
scatters around zero. Read the interval.

**`inference_status` is `withheld_inadequate_split_design`.** The panel is
too small or too correlated to split into training and test halves, so whether
the groups reproduce could not be tested. The groups are descriptive only, and
the record says so rather than releasing a p-value.

## Getting help

Open an issue in the repository with the configuration, `input_qc.md`, and
the message on standard error. Do not attach isolate-level data unless it is
already public.
