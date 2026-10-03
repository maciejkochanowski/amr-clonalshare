# *Escherichia coli* from diseased pigs: provenance and reading of the panel

481 isolates x 6 antimicrobials, one broth microdilution panel, typed by MLST
sequence type. Diseased pigs, China, 2011-2017; 472 isolates carry a sequence
type (75 types, 27 seen once) and 9 do not.

## 1. Source

Li X, Hu H, Zhu Y, Wang S, Jiang Y, Yang D, Li M, Wang Z, Wang X, Liu B, et
al. Population structure and antibiotic resistance of swine extraintestinal
pathogenic *Escherichia coli* from China. **Nature Communications**
2024;15:5811. doi:10.1038/s41467-024-50268-2. Licence CC BY 4.0.

The tables are derived from the workbook of Supplementary Data 1-14 served
beside the article (`41467_2024_50268_MOESM4_ESM.xlsx`, SHA-256
`b723215ec2dd20d889a7284858620094bd39b960f4ff5cdca9bee3365a5d27c9`):
Supplementary Data 1 (499 isolates with phylogroup, sequence type, tissue,
year, province, serotype and BioSample accession) and Supplementary Data 13
(the MIC of 20 antimicrobials for the 485 isolates with a susceptibility
record, with the tested range of every agent). `build_tables.py` reads the
workbook, checks its digest and writes the four tables of `data/`; the
workbook itself is not shipped.

## 2. Agents and the reading of the panel

Six agents are kept, each on the range the source reports for it: cefotaxime
1-32, ceftazidime 1-16, meropenem 1-8, ciprofloxacin 0.5-2, gentamicin 2-8
and tetracycline 2-8 mg/L, in twofold steps. These ranges are declared as
`mic_wells` in the configuration. A value written "<=1" lies at or below the
lowest well, ">16" above the highest, and an unprefixed value is the twofold
dilution interval ending at that concentration. All readings are in mg/L on
one panel (BD Phoenix), so the readings form one stratum.

Four isolates are set aside before any analysis, so that every agent is
analysed on the same isolates: A24 has no gentamicin reading, A62 and A71 no
ciprofloxacin reading, and the meropenem reading of A154, 3 mg/L, lies on no
well of the declared range. No value is imputed, rounded or moved. The nine
isolates whose sequence type is "-" in Supplementary Data 1 keep a blank
lineage and are counted by the run as untyped; sequence types seen once are
set aside by the lineage share and counted (27 isolates), which leaves 445
isolates in 48 sequence types.

## 3. The tables and configurations

| File | Content |
| --- | --- |
| `data/metadata.csv` | One row per isolate of Supplementary Data 1 (499): sequence type (`st`, blank where the source has none), phylogroup, year, tissue, province, serotype, BioSample |
| `data/mic_long.csv` | 2,886 readings: isolate, agent, operator, measurement, unit, panel |
| `data/mic_two_categories.csv` | The same readings merged into two categories at one well per agent: ciprofloxacin 1, meropenem 2, the other four agents 4 mg/L (at or below the well, or above it) |
| `data/mic_long_without_st410.csv` | The readings of `mic_long.csv` without the 74 isolates of ST410 (73 of them with readings) |
| `data/build_receipt.json` | Counts and settings written by `build_tables.py` |
| `data/metadata_periods.csv` | `data/metadata.csv` with a column `period`: 2011-2013 or 2014-2017 from the year of isolation, blank for the 11 isolates without a year (none of them has readings) |
| `data/calls_two_categories.csv` | The readings of `mic_two_categories.csv` written as 0/1 calls: 1 when the MIC lies above the cut of its agent |
| `input.csv` | One row per isolate with readings (481): the ciprofloxacin call, the sequence type, the phylogroup as published and the phylogroup with ST23, ST88 and ST410 in phylogroup C (`phylogroup_clermont2013`) |
| `data/derived_receipt.json` | Counts and the SHA-256 of every file read and written by `build_derived_tables.py` |

`config.yaml` reads the full panel, `config_two_categories.yaml` the merged
readings on a one-well range per agent, and `config_without_st410.yaml` the
full panel without ST410. `build_derived_tables.py` writes the last three
tables from the shipped ones, without the workbook. `config_periods.yaml`
splits the change in the share of calls above the cut from 2011-2013 (197
isolates) to 2014-2017 (284) into a change in the mix of phylogroups and a
change within phylogroups; the phylogroup is read because every
isolate carries one and both periods share all seven, whereas the sequence
types seen in both periods hold too few of the isolates (shared support
0.82, below the 0.9 the decomposition requires). `input.csv` compares the
sequence type with the phylogroup on the ciprofloxacin call, with the phylogroup as published and with clonal complex 23 in phylogroup C (`README.md`). The cuts of the two-category reading are
technical, chosen inside every range before the results were examined; they
are not clinical breakpoints or epidemiological cut-offs.
