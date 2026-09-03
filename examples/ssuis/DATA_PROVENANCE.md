# *Streptococcus suis* cohort: provenance, cut-offs and repairs

677 isolates x 13 antimicrobials, no missing cells. Pig 650, human 27. United
Kingdom 423, Canada 205, Vietnam 49. Collected 1983-2016; 108 sequence types
among the 458 isolates that carry one.

## 1. Source

Hadjirin NF, Miller EL, Murray GGR, Yen PLK, Phuc HD, Wileman TM,
Hernandez-Garcia J, Williamson SM, Parkhill J, Maskell DJ, Zhou R, Fittipaldi N,
Gottschalk M, Tucker AWD, Hoa NT, Welch JJ, Weinert LA. Large-scale genomic
analysis of antimicrobial resistance in the zoonotic pathogen *Streptococcus
suis*. **BMC Biology** 2021;19:191. doi:10.1186/s12915-021-01094-1
(PMID 34493269). Licence CC BY 4.0.

Measurements are replicated broth microdilution performed to CLSI M100-S25
(2015) and VET01S 3rd edition (2015). The isolates and their MIC records are
served by the Bacterial and Viral Bioinformatics Resource Center (BV-BRC,
taxon 1307) and were retrieved from there; every MIC row in `data/mic_long.csv`
carries `pmid = 34493269`, so the laboratory measurements are those of the
publication above and BV-BRC is the distribution channel rather than a second
source.

BV-BRC serves 677 of the 678 published isolates with a laboratory-method MIC
record. The cohort is the 677 that BV-BRC holds; no isolate was removed by a
decision made here.

### Retrieval

The BV-BRC field `evidence` is set to `Laboratory Method` on only a fraction of
the rows that carry a measured MIC, and filtering on it discards most of the
real laboratory data. The correct filter is `laboratory_typing_method` in
(`Broth dilution`, `Agar dilution`, `MIC`) with a populated
`measurement_value`, then dropping `evidence = Computational Method`, which is
predicted from the genome and would be circular as a phenotype anchor.

## 2. From MIC to non-wild-type

The phenotype is microbiological throughout: an MIC compared with an
epidemiological cut-off, giving wild type or non-wild-type. No S/I/R category
and no composite resistance index is produced anywhere in this example.

**Coding is `0` = wild type, `1` = non-wild-type.** This is the opposite of the
reading "1 means the isolate satisfies the named condition" and is stated here
because it has been misread before.

**Three of the thirteen agents have a published EUCAST cut-off and ten do
not.** For the ten the cut-off is derived from the observed MIC distribution by
the method EUCAST uses to set its own: iterative log-normal fitting of the
wild-type subpopulation at 99 % coverage, rounded up to the next two-fold
dilution (Turnidge, Kahlmeter & Kronvall 2006, *Clin Microbiol Infect*
12:418-425). `derive_ecoffs.py` performs both steps and writes
`data/ecoff_derived.json`.

| agent | cut-off (mg/L) | source | non-wild-type |
|---|---|---|---|
| doxycycline | 0.5 | published, tentative | 571 (84.3 %) |
| tetracycline | 4 | published | 574 (84.8 %) |
| erythromycin | 0.25 | published, tentative | 366 (54.1 %) |
| lincomycin | 2 | derived | 424 (62.6 %) |
| tylosin | 2 | derived | 366 (54.1 %) |
| tilmicosin | 8 | derived | 364 (53.8 %) |
| trimethoprim | 0.25 | derived | 190 (28.1 %) |
| ceftiofur | 0.25 | derived | 163 (24.1 %) |
| penicillin | 0.0625 | derived | 156 (23.0 %) |
| tiamulin | 8 | derived | 131 (19.4 %) |
| spectinomycin | 32 | derived | 78 (11.5 %) |
| amoxicillin | 0.0625 | derived | 37 (5.5 %) |
| cefquinome | 0.125 | derived | 26 (3.8 %) |

A derived cut-off is specific to the collection it was derived from. It is a
defensible cut-off for this cohort and is not proposed as a species-level
value. Any analysis that turns on the ten derived agents inherits that, which
is why the article reports the three published agents separately as a
sensitivity analysis.

### What validates the derivation

Where a published value exists it is used, and the derived value is compared
with it rather than discarded. Across the four agents in the full 16-agent
panel that have one:

| agent | derived | published EUCAST | agreement |
|---|---|---|---|
| florfenicol | 4 | 4 | exact |
| tetracycline | 8 | 4 | +1 dilution |
| doxycycline | 0.25 | (0.5) | -1 dilution |
| erythromycin | 0.0625 | (0.25) | -2 dilutions |

Three of four agree within one two-fold dilution, which is the reproducibility
of the MIC method itself (ISO 20776-1). Parenthesised values are tentative
EUCAST cut-offs.

Two implementation details changed the answer and are recorded because they
would otherwise be invisible. Submitters report the same dilution at different
precision (`0.03` and `0.031`, `0.12` and `0.125`, `0.015` and `0.016`); left
as separate bins they widen every fitted distribution and push the 99 % cut-off
one to two dilutions too high, and florfenicol came out at 16 instead of 4. The
endpoint search is anchored on the wild-type peak; without that anchor it walks
past the resistant population and fits everything as one wild type, which for
doxycycline, where 84 % of isolates are non-wild-type, returned 64 mg/L against
a published 0.5.

### Boundary rule

An isolate is non-wild-type when its MIC is **strictly greater** than the
cut-off. Every one of the 10,832 MIC records in this cohort is an exact
measurement; none is censored, so no interval rule is needed.

### Reproducing the matrix

`tests/test_ssuis_provenance.py` rebuilds all 8,801 cells of `ribo.csv` and
`cell.csv` from `data/mic_long.csv` and the cut-offs in
`data/ecoff_derived.json` and requires an exact match.

Read the identifier column as a string. `1307.7510` parsed as a float loses its
trailing zero and silently fails to join, which is the same trap that once
disabled the lineage diagnostic on this cohort: 458 sequence-type labels were
present and none of them matched.

## 3. Agents excluded

Three of the sixteen measured agents carry a minority class below 20 isolates
and are not in the matrix: florfenicol (0 non-wild-type), marbofloxacin (2),
enrofloxacin (16). A binary column with a class that small supports no
comparison and would enter every multiplicity correction as if it did.

## 4. Repairs applied to the metadata

`repair_metadata.py` performs these and is idempotent; run it with `--check` to
confirm the shipped file is in the repaired state.

**Serotype 1/2 arrives as the string `1-Feb`, in 39 of 677 rows.** *S. suis*
serotype 1/2 is a recognised serotype, distinct from 1 and from 2. Upstream of
the public record the string `1/2` passed through a spreadsheet, which read it
as a date and wrote it back as `1-Feb`. The corruption is present in the BV-BRC
dump as received and was not introduced here. It is the error Ziemann, Eren and
El-Osta documented for gene symbols (2016, *Genome Biology* 17:177), reaching a
serotype field. `serovar` is an external variable in this example and never
enters the clustering, so the repair changes reported cluster composition and
external agreement and changes no partition.

**Twelve rows carry free-text serotyping notes** rather than a serotype:
`1 + 6`, `8 (+2 partial)`, `(4 partial)`, and one reading
`(26 partia/31 partiall)`. These are left exactly as received. A note about a
partial agglutination is information, and guessing which serotype it meant
would not be.

## 5. Columns derived here

Neither is a repair. Both exist so that a two-collection contrast can be named
in a YAML configuration rather than computed inside a script.

* `collection_period`: `early` before 2013, `late` from 2013. The boundary is
  the gap in the United Kingdom series, which runs 2009-2011 and then 2013-2015
  with nothing in 2012.
* `country_period`: `isolation_country` joined to `collection_period`.

Among United Kingdom isolates that carry a sequence type, `late` is exactly
2013-2014 (n = 93) and `early` is exactly 2009-2011 (n = 177); the 14 United
Kingdom isolates from 2015 carry no sequence type and leave any
lineage-resolved analysis at that step.

## 6. Two lineage definitions, and why both are shipped

**`mlst` is missing for 219 of 677 isolates (32.3 %), and it is not missing at
random.** Among United Kingdom isolates collected from 2013 only 40.3 % carry a
sequence type, against 97.3 % of those collected before 2013
(Fisher exact p = 1.8e-39), and the untyped isolates are the more resistant:
across the thirteen agents they carry a mean of 6.7 non-wild-type results
against 4.3 for the typed ones. A contrast between those two periods computed
on sequence-typed isolates is therefore a contrast between two differently
selected subsets. For ceftiofur the sequence-typed subset falls by 11.0 points
while the collection it is drawn from rises by 9.5.

**`baps_cluster` is present for all 677.** It is the hierBAPS population
cluster reported by the source study, inferred from the core-genome alignment
(Cheng, Connor, Siren, Aanensen & Corander 2013, *Mol Biol Evol* 30:1224-1228),
and it is attached by `link_source_lineages.py`. It gives 29 clusters against
108 sequence types, which is what a core-genome population cluster should be
relative to a seven-locus type, and it raises the shared support of the United
Kingdom period contrast from 76.8 % of isolates to 94.5 %.

Both are shipped because the pair is the shortest demonstration of what the
estimability gate is for. The analyses in the article use `baps_cluster`.

### How the join is made and what validates it

BV-BRC carries the source study's own strain name inside `genome_name`
("Streptococcus suis SS1038"), and the longest source identifier contained in
that string identifies the isolate. This is an identifier join, not a
statistical match: all 677 isolates match, each to a distinct strain, with no
collisions.

Two variables validate it and neither takes part in it.

| check | result |
|---|---|
| all 10,832 MIC cells, BV-BRC against the published table | agree, 1.000 |
| collection year, where both tables hold it (651 isolates) | agree, 1.000 |
| adjusted mutual information, hierBAPS cluster against MLST | 0.720, z = 59.8 against its permutation null |

The first two are enforced as preconditions: `link_source_lineages.py` writes
nothing unless both are exact. The third is a consistency check rather than a
gate, and it recovers known population structure: ST1, the dominant clone,
falls entirely inside one cluster.

The MIC check is also a provenance result in its own right. The two tables were
prepared by different routes from the same laboratory work, and every one of
10,832 measurements agrees.

## 7. The determinant layer

The source study reports presence or absence of 43 resistance determinants,
called from the genome assemblies of the same isolates. These are used in one
place only, `validate_with_determinants.py`, and never in the clustering or in
the phenotype analysis. They are an independent measurement of the same
isolates by a different instrument, which is what makes them a check on the
decomposition rather than an illustration of it.

## 8. Files

| file | contents |
|---|---|
| `data/mic_long.csv` | 10,832 raw MIC records, 677 isolates x 16 agents, complete |
| `data/ecoff_derived.json` | per-agent fit: mu, sigma, R-squared, endpoint, derived cut-off, published value |
| `data/ribo.csv` | 677 x 8 non-wild-type calls, ribosomal target site |
| `data/cell.csv` | 677 x 5 non-wild-type calls, cell wall and folate targets |
| `data/metadata.csv` | strain name, host, country, year, serotype, MLST, hierBAPS cluster, derived period keys |
| `data/source_table_s1.csv` | Additional file 1 of the source publication, flattened |
| `data/lineage_link_receipt.json` | the join and its validation |
| `fetch_source_table.py` | retrieves Additional file 1 and flattens it |
| `derive_ecoffs.py` | cut-off retrieval, derivation and validation; rewrites `ecoff_derived.json` |
| `repair_metadata.py` | the repairs and derived columns in sections 4 and 5 |
| `link_source_lineages.py` | the hierBAPS join and its two validations |
| `decompose_trend.py` | worked prevalence-difference decomposition over the panel |
| `validate_with_determinants.py` | the same decomposition on the genotype layer |

## 9. Indication of changes, CC BY 4.0 section 3(a)(1)(B)

The licence is <https://creativecommons.org/licenses/by/4.0/>. The material was
modified as follows and in no other way: MIC values were dichotomised against
the cut-offs in section 2; three agents were dropped as in section 3; the
serotype field was repaired and two period keys were derived as in sections 4
and 5. No MIC value was altered.
