# 2 Prepare the input

Use a metadata CSV with one isolate identifier and its lineage label, and a calls CSV, a MIC CSV, or both; each may be long (one row per isolate and agent) or wide (one column per agent). Identifier spelling must match exactly. The recorded observation frame is the union of identifiers in the calls and MIC tables. Metadata-only identifiers are inventoried; they are not automatically sampled isolates. Each analytical module then retains the outcomes and labels it needs, so MIC-only isolates are not discarded because they lack calls.

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

`name` labels the collection; it is written into the record and the reports and enters the run identifier, and is read nowhere else. `data_dir` is relative to the configuration file.

Tables are read as a laboratory system or a spreadsheet exports them. The separator of a text table is the one its header line uses most, so a semicolon-separated export with decimal commas (`0,5`) and a tab-separated table copied from a spreadsheet read like a comma-separated one; a byte-order mark, spaces around column names and empty rows at the end are removed; `n.d.`, `NA`, `-` and a blank are read as a value that was not recorded. A file in an encoding other than UTF-8 is refused with the instruction to save it as UTF-8, since guessing an encoding would guess the lineage labels.

Agent names are unified across the tables: a WHONET or EUCAST code (`CIP`, `CIP_NM`, `SXT`) or an abbreviation is read as the agent it stands for, so that the call table and the MIC table meet on one name, and the names changed are listed in the input check and the record. A name the package does not know is kept as written. The recipe directory supplies complete examples. Record the label definition, clustering pipeline/version or serotyping vocabulary with the data. The software cannot recover that provenance from a short label.

## One table is enough

The two tables may be one file. A long export that carries the lineage beside every call, one row per isolate and antimicrobial, is named as both `metadata` and `phenotype`:

```yaml
dataset:
  name: laboratory_collection
  data_dir: data
  metadata: calls.csv
  phenotype: calls.csv
  strain_id_column: isolate_id
  lineage_column: lineage
  phenotype_id_column: isolate_id
  phenotype_antibiotic_column: agent
  phenotype_call_column: call
  phenotype_kind: clinical_sir
```

The isolate then repeats in the metadata table, once per antimicrobial. Repeated rows collapse when they agree in the columns this analysis reads: the lineage and any contrast, batch, stratifying, sampling-unit, panel or covariate column the configuration names. A repeat that disagrees in one of those is still refused, and the error names the column. A column the analysis never opens cannot make a repeat a conflict. Identifiers cannot fail to match across tables when there is only one, which removes the most common preparation error.

## A first configuration from the column names

`--init` reads each table it is given, prints a configuration with the mapping it guessed, and stops. It writes no file and runs no analysis:

```bash
amr-clonalshare --init data/metadata.csv data/mic.csv > analysis.yaml
```

The vocabulary is read from the recorded values rather than guessed from the column name, so a table of `S` and `R` is declared `clinical_sir` and not left undeclared. A table whose values are concentrations on the twofold dilution series, with or without a sign, is recognised as the MIC table, in the long shape (an agent column beside a concentration column) or the wide shape (one concentration column per agent), and the draft writes the `mic:` block for it with `mic_wells` filled by the concentrations observed per agent. Those are the concentrations the table holds, not the range the panel tested: the comment above them says to replace them with the plate's range, lowest to highest well, or to delete the block and let the range be taken from the readings. A table that carries the lineage column beside the readings is named as the metadata too, so one sheet is a complete input. Commented templates for strata, sampling units and the two-collection comparison follow. Anything the names do not settle is written as `REPLACE_ME`; nothing is invented. `data_dir` is written as an absolute path and the tables relative to it, so the draft runs wherever it is saved, and every value is quoted, so a column named `Isolate #` or `2019` reads back unchanged. Read the result before running it: the guesses are mechanical, and only the recorded provenance of the lineage labels can tell you whether the column named `cluster` is the one you mean.

## A sheet with one column per agent

A laboratory sheet usually holds one row per isolate and one column per antimicrobial. Name those columns and the table is read in that shape:

```yaml
dataset:
  phenotype: calls.csv
  phenotype_id_column: isolate_id
  phenotype_agent_columns: [ampicillin, tetracycline]
  phenotype_kind: clinical_sir
```

The table is melted to the long shape before anything reads it, so the vocabulary, the duplicate policy, the exclusion counts and every reported number are the same as for the long shape; a wide table and the long table it melts to give identical panels. `phenotype_agent_columns` replaces `phenotype_antibiotic_column` and `phenotype_call_column`, which name the two columns of the long shape instead. A named column the table does not carry is refused. `--init` proposes this form when no column names a call and the remaining columns hold values that fit one vocabulary.

The MIC table has the same two shapes. In the wide shape each cell holds the recorded concentration with its sign (`<=0.5`, `8`, `>16`):

```yaml
dataset:
  mic: mic.csv
  mic_id_column: isolate_id
  mic_agent_columns: [CIP, GEN, TET]
```

`mic_agent_columns` replaces `mic_antibiotic_column`, `mic_value_column` and `mic_operator_column`; the sign is read inside the cell, so an operator column has no place in this shape and is refused. A panel, unit or covariate column of the table is kept beside every reading. The wide and the long table of the same readings give the same record.

## Workbooks

Any table may be an `.xlsx` or `.xlsm` workbook instead of a CSV. The workbook must hold one sheet: a workbook with several sheets is refused and its sheets are named, so that a table is never read from a sheet the user did not mean; save the sheet to analyse as its own workbook or as CSV. A damaged or old-format file is refused with the same instruction. Reading a workbook needs `openpyxl`, which the core install leaves out to keep the dependency stack to the numerical one:

```bash
pip install "amr-clonalshare[excel]"
```

A run given a workbook without it is refused with that instruction rather than a traceback. A workbook and the CSV of the same sheet give the same numbers; the CSV remains the format to archive, because it is the one a reader can open in ten years without a library.



## Declare the phenotype

| `phenotype_kind` | Intended input | Positive outcome |
|---|---|---|
| `clinical_sir` | Clinical susceptibility S/I/R | R only by default; state any explicit intermediate policy |
| `wt_nwt` | Wild type/non-wild type | NWT; do not relabel this clinical resistance |
| `binary` | Recorded 0/1 | 1; describe the meaning in `phenotype_positive_definition` |
| `undeclared` | Resistant, intermediate, susceptible and non-susceptible, read when no kind is declared | R or I; the run warns |

Use `phenotype_source`, `ast_standard` and `ast_version` when applicable and known. Do not invent a standard for collection-derived thresholds. `phenotype_positive_definition` documents the meaning; it does not calculate clinical breakpoints. For clinical data, `phenotype_intermediate: non_susceptible` explicitly combines I and R; `susceptible` maps I to zero and `drop` excludes I. The definition and policy must agree. WT/NWT and binary configurations do not accept an intermediate policy.

Blank calls, unrecognized strings, absent agent records, untyped isolates and empty agents are counted separately in QC. None is automatically a negative outcome. Inspect these counts before interpreting a prevalence denominator.

## Resolve duplicates explicitly

Identical repeated records collapse. Conflicting records raise an error by default (`duplicate_policy: error`). `drop_conflicts` excludes conflicting keys and records the exclusions. `positive_wins` counts a conflicting call as positive if any of its records is, and keeps the first record of other tables; the run warns. Do not select a permissive mode merely to remove an error: resolve the source discrepancy or state why exclusion is justified. These checks apply to the relevant metadata, call and MIC keys.

## Preserve MIC measurement information

```yaml
  mic: mic.csv
  mic_id_column: isolate_id
  mic_antibiotic_column: agent
  mic_value_column: measurement
  mic_operator_column: operator
  mic_unit_column: unit
  mic_panel_column: laboratory
  mic_covariate_columns: [country, year]
  stratify_by: country
  mic_wells:
    demo_agent: [0.25, 0.5, 1, 2, 4, 8]
  mic_units:
    demo_agent: mg/L
```

Place these keys inside `dataset`. The well list is agent-specific, strictly increasing, positive and finite; a reading is matched to the well it was read on when it is written either way, so `0.12` and `0.125`, or `0.015` and `0.015625`, are one well, and a reading farther than a quarter of a doubling from every well is refused. For the panels the European Union prescribes for the monitoring of resistance in food-producing animals, the ranges are shipped as presets transcribed from Commission Implementing Decision (EU) 2020/1729 (Annex, Part A, Tables 2 to 5):

```yaml
  mic_panel_preset: eu-2020-1729-salmonella-ecoli
```

The presets are `eu-2020-1729-salmonella-ecoli`, `eu-2020-1729-salmonella-ecoli-second`, `eu-2020-1729-campylobacter` and `eu-2020-1729-enterococcus`. A preset fills `mic_wells` for every agent of the panel that `mic_wells` does not name itself; an agent written explicitly keeps its own range, and an agent of the preset that the table does not carry is ignored. The record and the report name the preset as the source of the wells and say that the plate sheet of the laboratory remains the authority, since a product or a lot may differ from the Decision. Without `mic_wells` or a preset, the lowest and the highest recorded concentration of an agent are taken as the end wells of its panel; the input check lists the recorded concentrations of every agent with the share of readings on those two, and asks for the tested range where that share is large. The operator column can record `<`, `<=`, `>` or `>=` (the Unicode signs for less than or equal and greater than or equal are also read); keep the numerical concentration in the measurement column. `>x` places the MIC above x and `>=x` above the tested well below x; `<x` and `<=x` both place it at or below x, since exports print the lowest well either way. A sign may also lead the value itself, and a combination written as `0.5/9.5` is read on its first component. A censored value that does not parse is refused rather than set aside. Explicit operators take precedence over end-well assumptions: where the export writes a sign only on censored readings, a blank reading on the highest tested well beside a `>` or `>=` reading at that well is read as the well itself, (previous well, well], since growth at the well and no growth at it are two readings. `mic_unit_column` names a column holding the unit of every reading: a reading in mg/L (also written µg/mL or ug/mL) is read, a row in any other unit, such as a zone diameter in mm that an export keeps beside the MICs, is set aside and counted per agent under `mic_join`, and a row without a unit is read and counted; the unit read is recorded for every agent that `mic_units` does not name. Units are recorded, not automatically converted. Supply actual panel wells; wells inferred from observed concentrations can omit tested but latent concentrations. Check the reported geometry.

`mic_panel_column` names the panel a reading was made on: the laboratory, the panel product, or any label under which every isolate was tested on the same wells. It may be a column of the MIC table or, per isolate, of the metadata. Without it, the wells are inferred from the pooled readings of an agent, and when two laboratories tested different ranges the lowest well of one passes as an interior dilution of the other, so its left-censored readings are read as exact. With it, the wells and the end wells are inferred within each label, the geometry is reported per label, and the share of the MIC ordering places every reading among the readings of its own label and permutes lineage labels only within it. Every reading needs a label; `mic_wells`, when given, apply to every label alike.

`mic_covariate_columns` lists categorical covariates of the MIC reading, in the MIC table or, per isolate, in the metadata: the country, the year of isolation, or another factor whose levels may shift MICs. Their levels, crossed with each other and with the panel label, define the strata of the MIC analysis: every reading is scored within its stratum, lineage labels are permuted only within strata, and the latent ordering the bounds describe is the ordering within strata, so that a shift between levels is not read as a lineage difference. `mic_covariate_column` names one covariate and may be combined with the list; the names must be distinct. Every reading needs a value in every covariate. A level that follows lineage exactly removes the lineage difference with it, since lineages read only in different strata are not compared.

`stratify_by` names a metadata column whose levels define strata of another kind. The whole run is repeated on the isolates of each level, with every analysis on its own lineages and, for a MIC table, its own panel; the pooled run is unchanged and the strata are added to the record under `strata`, to the report as a summary table, and to `strata_results.csv`. Isolates without a level are left out of the strata and counted. Stratification asks whether a share holds within each level on its own; a MIC covariate keeps one analysis and reads the ordering within its levels. The two can be used together.

Metadata may also supply `contrast_column` with exactly two `contrast_levels`, or a `batch_column` that sorts into true arrival order. A sequential analysis requires that ordering and the method's null assumptions; retrospective reordering changes the analysis.

## Declare strata for the calls

```yaml
  phenotype_covariate_column: laboratory
```

`phenotype_covariate_column`, inside `dataset`, names a metadata column whose levels are the strata of the analyses of a call: the testing laboratory, the country, or another factor whose levels may shift the calls. Each call is then centred on the prevalence of its level, the score of a reading with one cut point, so the share is the share of the variance within levels that lies between the lineages; the permutation test exchanges lineage labels only within a level, the interval draws each isolate's level with its call and recomputes the scores in every draw, and the bounds a call places on the latent ordering are the bounds within levels, so that a shift between levels is not read as a lineage difference; the record names the column and the number of levels under `call_strata`. Every isolate with a call needs a value. The e-values test one common probability and take no strata; the decomposition of a prevalence difference and the comparison of two lineage definitions read the calls as recorded. A level that follows lineage exactly removes the lineage difference with it, as it does for the MIC readings.

## Declare sampling units

```yaml
  unit_column: farm
```

`unit_column`, inside `dataset`, names a metadata column holding the sampling unit of every isolate, such as the farm or the host. The permutation tests of the lineage shares then exchange lineage labels only within units (and within the strata of a MIC analysis): a test of lineage differences within units, which stays exact when the isolates of one unit share more than their lineage. The decomposition of a prevalence difference then resamples whole units. The shares, their intervals and the bounds are unchanged and describe the lineages and the represented units. An isolate without a recorded unit is a unit of its own; the record counts them under `sampling_units`.
