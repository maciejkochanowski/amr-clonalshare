# 2. Prepare the input

Three kinds of file, all plain CSV, all declared in one YAML configuration.

## 2.1 Trait layers (required)

One file per layer. A layer is a group of traits measured the same way and
sharing a mechanism grain: for example every acquired resistance gene, or every
non-wild-type call from one susceptibility panel, or a virulence locus panel.

| rule | why |
|---|---|
| one row per isolate, one column per trait | the row is the unit the clustering and the shares are computed on |
| a strain-identifier column, named in the config (`strain_id_column`) | identifiers are joined as written; nothing is stripped or case-folded, so an identifier that differs by a space is a different isolate |
| every cell exactly `0` or `1` | `1` is the state the tool treats as presence (the gene, the non-wild-type call); any other value is refused |
| no empty cells, or a stated policy for them | a binary layer has no state for "not measured"; see below |
| identifiers unique within a file | a duplicated identifier is refused with the duplicates named |

```text
isolate,tetM,ermB,lnuB,aac6-aph2
S1,1,1,0,0
S2,1,0,0,1
S3,0,0,0,0
```

!!! warning "Empty cells are not zeros"
    If a trait was not measured for an isolate, leave the cell empty and set
    `dataset.missing_policy` (below). Writing `0` there tells the tool the
    trait is absent, and the estimate will treat it as such. The default policy
    is to refuse a file with empty cells and name the rows and columns.

A `one_hot` layer is a categorical variable such as capsule type spread over
columns, at most one `1` per row. Declare it as `kind: one_hot`; it is then
given a matching-coefficient distance, validated at load, and counted as one
hypothesis rather than one per level.

## 2.2 Metadata (needed for the lineage question)

One row per isolate, the same identifier column, and at least the lineage
column: a sequence type, a hierBAPS or PopPUNK cluster, an NCBI Pathogen
Detection SNP cluster, or a serovar. Read every column as text; the tool does.

```text
isolate,lineage,host,country,year
S1,ST1,pig,PL,2021
S2,ST1,pig,PL,2021
S3,ST28,pig,DE,2022
```

Which lineage definition to use is a scientific choice, not a setting, and the
tool measures its consequence: a finer definition (SNP cluster) leaves more
isolates in single-member lineages and lowers the support the estimator needs;
a coarser one (serovar) raises support at the price of a lower share. Run both
if you hold both; the report says which cells each can and cannot estimate.

An isolate with an empty lineage cell is kept and counted as untyped. It joins
its own level, the share of untyped isolates is reported, and if typing failed
non-randomly that fact is visible in the record rather than hidden by dropping
the rows.

## 2.3 Recorded MICs (optional)

A long table, one row per isolate and antimicrobial, with the measured
concentration and, where the laboratory recorded one, the censoring operator
(`<=`, `>`). The tool reads a recorded dilution as an interval on the
log2 scale and reports the panel geometry before any share, so a value on the
top well of the panel is treated as "at least" rather than "equal to".

```text
genome_id,antibiotic,measurement,operator
S1,tetracycline,128,>
S1,erythromycin,0.25,<=
```

## 2.4 The configuration

```yaml
dataset:
  name: my_cohort
  strain_id_column: isolate
  data_dir: data                # every path below is relative to this
  metadata: metadata.csv
  lineage_column: lineage
  missing_policy: refuse        # refuse | drop_rows | drop_columns
  strain_alignment_policy: intersect_core

files:
  amr: { path: amr.csv, kind: wide_binary }
  vir: { path: vir.csv, kind: wide_binary }

trait_cluster:
  layers: [amr, vir]
  k_range: [2, 3, 4, 5]
```

Every key the loader does not read is refused with the nearest known key
named. A misspelt `lineage_colum` would otherwise silently switch the lineage
diagnostics off, and the output would be indistinguishable from a cohort with
no lineage labels.

`strain_alignment_policy` decides which isolates enter when layers do not hold
the same identifiers: `intersect_core` keeps those present in every layer
(default); `union` keeps every isolate and leaves the cells of the layers that
lack it empty, so it only runs under `missing_policy: drop_rows` or
`drop_columns`; `strict_n` refuses unless exactly `expected_n` align.

The shipped `examples/ssuis/config.yaml` is a fully annotated real
configuration: 677 *Streptococcus suis* isolates, 13 antimicrobials in two
layers split by molecular target, a hierBAPS lineage, and a 16-agent MIC panel.
Each decision in it is stated together with the measurement on that cohort
that drove it.
