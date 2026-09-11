# 2. Prepare the input

Two tables, plain CSV, declared in one YAML configuration, and an optional
third.

## 2.1 Susceptibility calls (required)

One long table, one row per isolate and antimicrobial, which is the shape a
laboratory exports and the shape NCBI Pathogen Detection and BV-BRC serve.

| rule | why |
|---|---|
| one row per isolate and antimicrobial | the isolates of the run are the isolates this table holds; an isolate tested for no agent is not in the run |
| an identifier column, named in the config (`phenotype_id_column`) | identifiers are joined as written; nothing is stripped or case-folded, so an identifier that differs by a space is a different isolate |
| a call column holding `susceptible`, `non-susceptible`, `resistant` or `intermediate`, in any case | `1` is the state the tool treats as non-wild-type; `phenotype_intermediate` states whether an intermediate counts as non-susceptible (default), susceptible, or is dropped; any other word is not read |
| a repeated (isolate, antimicrobial) pair | the non-susceptible result wins, which is the conservative reading for surveillance |

```text
isolate,antibiotic,call
S1,tetracycline,non-susceptible
S1,erythromycin,susceptible
S2,tetracycline,non-susceptible
```

A pair that was not tested is simply absent, and stays distinguishable from
"tested and susceptible"; an agent no isolate carries a readable call for is
not read.

## 2.2 Metadata (required)

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

An isolate with an empty lineage cell is counted as untyped and set aside
before the share is estimated: an untyped isolate is not a member of any
lineage, and a level made of every isolate whose typing failed would measure
the typing process rather than lineage membership. The share of untyped
isolates (`missing_share`) and their number (`n_dropped_untyped`) are in the
record, so an estimate that rests on a typed subset says so; if typing failed
non-randomly, the decomposition's label-availability test names it.

## 2.3 Recorded MICs (optional)

A long table, one row per isolate and antimicrobial, with the measured
concentration and, where the laboratory recorded one, the censoring operator
(`<=`, `>`). The tool reads a recorded dilution as an interval on the
log2 scale and reports the panel geometry before any share, so a value on the
top well of the panel is treated as "at least" rather than "equal to". The
tested wells are taken as every doubling from the lowest to the highest
recorded value when the values sit on a doubling lattice, so a dilution no
isolate landed on still bounds its neighbours and 0.06 and 0.064 are read as
one well; a laboratory that knows its panel can pass the tested
concentrations to `intervals_from_mic(..., wells=...)`.

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
  batch_column: year            # optional: the intake, sorting into arrival order
  phenotype: calls.csv
  phenotype_id_column: isolate
  phenotype_antibiotic_column: antibiotic
  phenotype_call_column: call
  mic: mic.csv                  # optional: recorded dilutions
  mic_id_column: isolate
  mic_antibiotic_column: antibiotic
  mic_value_column: measurement
```

Four further sections are optional and carry the budgets and thresholds of
the estimators; a section left out takes the defaults shown, and a section
given with a value of the wrong type is refused:

```yaml
attribution:                    # the clonal share per antimicrobial
  folds: 5
  repeats: 20                   # fold draws averaged per estimate
  n_boot: 400                   # lineage bootstrap draws
  n_perm: 200                   # permutations for the control and the p-value
surveillance:                   # the lineage-resolved readings and the decomposition
  n_boot: 2000
  n_perm: 2000
  q_fdr: 0.05                   # target false discovery rate per component family
  min_shared_support: 0.8       # the shared-support gate
  label_alpha: 0.05             # the two label-availability tests
censored:                       # the reading of recorded dilutions
  n_boot: 200
  end_wells_censored: true      # a reading on an end well is an interval
  sensitivity: true             # also read with the end wells taken as points
evidence:                       # the e-value per antimicrobial
  folds: 5
  repeats: 20
  alpha: 0.05                   # the e-BH level
```

Every key the loader does not read is refused with the nearest known key
named. A misspelt `batch_colum` would otherwise silently switch the sequential
evidence off, and the output would be indistinguishable from a cohort with no
intakes.

The shipped `examples/ssuis/config.yaml` is a fully annotated real
configuration: 677 *Streptococcus suis* isolates, 13 antimicrobials called
against an epidemiological cut-off, a hierBAPS lineage, the yearly intakes and
a 16-agent MIC panel. `config_contrast.yaml` reads the same cohort at
sequence-type resolution with two countries declared as the collections to
contrast. `examples/salmonella_poultry/config.yaml` and `config_cluster.yaml`
read the poultry-meat *Salmonella* cell of one NCBI Pathogen Detection
release, 7,049 isolates with the calls the release records, at serovar and at
SNP-cluster resolution: the same two tables, a different lineage column.
