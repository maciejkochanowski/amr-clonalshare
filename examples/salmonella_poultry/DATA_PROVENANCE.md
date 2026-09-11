# Poultry-meat *Salmonella* cell: provenance

7,049 isolates x 22 antimicrobials, one NCBI Pathogen Detection release,
poultry meat, United States (6,902) and Mexico (132) with 15 unplaced,
collected 2002 to 2021; 85 serovars, 534 SNP clusters. This is the cell
Section 3.1 of the article reads, and the two configurations beside this
file read it at the two typing resolutions the release carries.

## 1. Source

NCBI Pathogen Detection, organism group *Salmonella*, release
PDG000000002.4210, the accession pinned in
`benchmarks/pathogen_detection_releases.tsv` and retrieved by
`benchmarks/fetch_pathogen_detection.sh`. Two tables of that release are
read: the isolate metadata with its `AST_phenotypes` field, which carries the
susceptibility calls the submitting laboratories deposited, one
`agent=category` pair per antimicrobial, and the cluster assignment file,
which gives every isolate its SNP cluster (`PDS` accession) at that release.
The data are public and carry no licence restriction; the isolates were
deposited by the United States National Antimicrobial Resistance Monitoring
System and other submitters, and the release is the version the atlases of
supplement sections S11, S13 and S15 were run on.

NCBI recomputes the SNP clusters at every release, so a later release is a
different lineage variable rather than a longer version of this one. A reader
who wants this cell as shipped fetches the pinned release; a reader who
fetches `latest_snps` gets a new retrieval, which is a new cohort at cluster
resolution and the same one, up to added isolates, at serovar resolution.

## 2. How the cell is cut

`build_cell.py` keeps an isolate when all four conditions hold: the release
places it in a SNP cluster; the host and sampling-matrix rules of
`benchmarks/vet_source_taxonomy.py`, applied to the `host` and
`isolation_source` fields in the fixed order those rules document, assign it
to poultry (chicken, turkey or duck) and to the food matrix, which is retail
or processing-plant meat rather than a live animal, a farm or a clinical
sample; the organism name carries a serovar; and the collection date carries
a year. No isolate is removed by any other decision, and no call is imputed:
an isolate not tested against an agent has no row for it.

The nalidixic-acid reading of the article is the 6,915 of these isolates
with a nalidixic-acid call, of which 584 are non-susceptible, in 84 serovars
and 496 clusters; those counts are what `--check-input` and the run report
for that agent.

## 3. Columns

`data/metadata.csv`, one row per isolate:

| column | source |
|---|---|
| `genome_id` | the isolate accession (`target_acc`) |
| `serovar` | the text after "serovar" in the release's organism name, as written; antigenic formulae such as `4,[5],12:i:-` are kept as they are, and `-:r:1,5` is not merged into Infantis here (supplement S20 measures what that merge is worth) |
| `pds_cluster` | the SNP cluster accession of the release |
| `collection_year` | the first four characters of the collection date |
| `period` | `before_2016` or `from_2016`, the cut the article fixed in advance on the 2014 detection of the emergent Infantis clone in United States retail meat |
| `country` | the part of `geo_loc_name` before the first colon |
| `isolation_source` | the free-text field the matrix rule read |

`data/calls_long.csv`, one row per isolate and antimicrobial: `genome_id`,
`antibiotic` as the release spells it, and `call`, `susceptible` for S and
`non-susceptible` for I or R, the reading NCBI applies to its own categories.
`data/cell_receipt.json` records the counts at the time the tables were
written.

## 4. What the shipped runs show

`config.yaml` reads the serovar as the lineage and `config_cluster.yaml` the
SNP cluster, with the budgets the article's reading used (five folds, ten
repeats, 300 bootstrap draws, 150 permutations). The expected outputs of both
are in `expected_serovar/` and `expected_cluster/`. For nalidixic acid the
serovar run returns 0.765 (0.017 to 0.851) and the cluster run 0.937 (0.233
to 0.975), against 0.766 (0.016 to 0.853) and 0.935 (0.247 to 0.977) in
Section 3.1 of the article, with the permuted control near zero at both; the
exact realised interval is withheld on nalidixic acid and on most of the
panel, because a call at a few per cent prevalence has within-lineage
residuals far heavier-tailed than the Gaussian law that interval assumes, and
opens only on the agents whose prevalence sits nearer one half; the report
says which. The article's
figures were computed by `benchmarks/period_split.py` on the same tables with
the same settings and a fixed seed; the shipped run draws its random stream
from the master seed of the command line, so the two agree to the second
decimal place and not to the last digit.

Adding `batch_column: collection_year` to either configuration turns the
cell into a programme of yearly intakes and adds the sequential e-process
over them; adding `contrast_column: period` with the two levels asks for the
Kitagawa split between the periods. Neither is part of the article's reading
of this cell, so neither is switched on here.
