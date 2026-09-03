# Provenance and licence of the *Klebsiella pneumoniae* case-study data

## What these files are

`data/` contains Kleborate genotype calls for **1500 *Klebsiella pneumoniae***
isolates, together with isolate metadata:

| file | shape | content |
|---|---|---|
| `amr.csv` | 1500 × 17 | **acquired** antimicrobial-resistance **class** calls (12 classes + 5 β-lactamase subclasses) |
| `vir.csv` | 1500 × 38 | virulence **gene** calls across the *ybt*, *clb*, *iuc*, *iro* and *rmp* loci |
| `kloc.csv` | 1500 × 78 | capsule **K locus**, one-hot |
| `otype.csv` | 1500 × 12 | **O antigen**, one-hot |
| `metadata.csv` | 1500 × 10 | ST, country, region, year, source, K locus, O type, Kleborate virulence score (0–5), resistance score (0–3), convergence-event label |

`amr.csv`, `vir.csv`, `kloc.csv` and `otype.csv` are derived from the same
underlying Kleborate output; `kloc.csv` and `otype.csv` were produced by
splitting the single capsule matrix of earlier releases into its two
categorical variables, because K locus and O antigen vary quasi-independently
and a one-hot block of both at once has no usable set-overlap structure (see
`docs/methodology.md`, "One-hot layers").

## What a zero in `amr.csv` means — read this before interpreting any cluster

**A zero is "no acquired determinant of this class in Kleborate's database",
not "susceptible".** The panel is acquired-gene-only. `Col` captures *mcr*
genes; `Flq` captures *qnr* and *aac(6')-Ib-cr*. Chromosomal and mutational
mechanisms are invisible to it: *mgrB* / *pmrB* for colistin, *ramR* / *rpsJ* /
*tet(A)* and efflux for tigecycline, *gyrA* / *parC* for fluoroquinolones,
OmpK35/36 loss for carbapenems.

These files contain the counter-example directly. Of the 38 isolates whose
Kleborate `resistance_score` is 3 — carbapenemase **plus** colistin resistance —
**33 have `Col = 0`**, because Kleborate's score uses a broader colistin
definition than the acquired-gene column does. In *K. pneumoniae*, *mgrB*
inactivation rather than *mcr* carriage is the predominant colistin-resistance
mechanism (Hu et al. 2023, *Int J Antimicrob Agents* 62:106873,
doi:[10.1016/j.ijantimicag.2023.106873](https://doi.org/10.1016/j.ijantimicag.2023.106873)),
and tigecycline resistance is predominantly efflux-mediated through *ramR* and
*tet(A)* (Li et al. 2025, *Front Cell Infect Microbiol* 15:1540967,
doi:[10.3389/fcimb.2025.1540967](https://doi.org/10.3389/fcimb.2025.1540967)).
`Tgc` is 0/1500 here, and that says nothing whatever about tigecycline
resistance in this cohort.

Any statement that a cluster is "resistance-poor" is therefore a statement
about acquired determinants only. Kleborate's `Col_mutations`,
`Flq_mutations`, `Omp_mutations` and `Bla_chr` columns are **not** shipped, so
the case study cannot address mutational resistance at all.

## The two Kleborate scores are not external validation

`virulence_score` and `resistance_score` are deterministic functions of the
same columns being clustered: reconstructing them from the shipped layers using
Kleborate's published rules reproduces `virulence_score` for 1495/1500 isolates
(99.7 %) and `resistance_score` for 1467/1500 (97.8 %, every discrepancy being
the colistin definition above). Scoring a partition against them measures
internal consistency, not external agreement, and "the virulence layer alone
recovers the published virulence score at ARI 0.875" is close to a tautology.

That gap is now partly closed. `fetch_phenotypes.py` retrieves
laboratory-measured antimicrobial susceptibility from BV-BRC (Olson et al. 2023,
*Nucleic Acids Research* 51:D678-D689,
doi:[10.1093/nar/gkac1003](https://doi.org/10.1093/nar/gkac1003)) for the
isolates whose `Strain_ID` is an assembly accession, writing
`data/ast_phenotypes.csv`.

| | |
|---|---|
| isolates with an assembly accession | 1074 / 1500 |
| with >= 1 laboratory AST record in BV-BRC | 333 |
| **with >= 1 usable S/I/R call** | **207 (13.8 % of the cohort)** |
| AST records retrieved | 2205 |
| of those, carrying an MIC but no S/I/R string | **508 (23 %), discarded** |
| usable S/I/R calls | 1697 |
| antibiotics | 37 (30 scorable, 19 in the FDR family, 12 significant) |

`resolve_genome_ids` keeps **one** BV-BRC genome per assembly accession, and
about 50 accessions map to several (up to seven). Which one is kept depends on
the order BV-BRC returns them, so the retrieval is not bit-reproducible; keeping
all of them adds roughly two isolates and forty calls. That is a known defect,
recorded here rather than fixed.

Two filters matter. Only `evidence = "Laboratory Method"` records are kept: the
same BV-BRC table carries 1.8 M AdaBoost *predictions* for *K. pneumoniae*, and
a predicted phenotype is a function of the genome, so scoring against it would
reproduce exactly the circularity above. And BV-BRC keys on GenBank (`GCA_`)
accessions while this cohort mostly carries RefSeq (`GCF_`) ones; paired
assemblies share the numeric part, so `GCF_003830175.1` is looked up as
`GCA_003830175.1`.

**Coverage is 13.8 %, and it is not a random 13.8 %.** Isolates reach BV-BRC's
AST table because someone measured and deposited susceptibility, which is more
likely for clinical than environmental isolates and more likely for resistant
than susceptible ones. Every phenotype statistic in the case study is therefore
conditional on that ascertainment, and none of them is a prevalence estimate.

What the anchor says is in §5.4 of the manuscript. In short: the partition does
predict measured non-susceptibility (12 of the 19 antibiotics in the FDR family
are significant), and it is **significantly beaten by one ordinal metadata
column** — Kleborate's `resistance_score`, +0.090 balanced accuracy, bootstrap
95 % CI [+0.048, +0.125]. It is *not* significantly beaten by the one-bit
"carries any acquired determinant" (+0.010, CI [-0.019, +0.041]); an earlier
draft claimed it was, on a point estimate with no uncertainty attached.
Balanced accuracy here has a permutation floor of about 0.55, not 0.50, because
it is maximised over both label orientations. Colistin, whose resistance in this
species is chromosomal, is predicted by nothing in the panel (0.54 against a
floor of 0.55, *p* = 1.0) - the "a zero is not susceptibility" point measured
rather than argued.

**The 207 are not a random 13.8 %.** 67 of them come from one Caribbean
collection, against a whole-cohort composition dominated by the USA, the UK,
Australia and China, and 191 of 207 are human. The retained calls also pool CLSI
and EUCAST across breakpoint years, and BV-BRC's `antibiotic` field is free text
(trimethoprim/sulfamethoxazole appears under two submitter strings with zero
isolate overlap). None of the phenotype numbers is a prevalence estimate and all
of them are conditional on which collections happened to deposit AST.

## Source

These are **derived data, not original observations**. They are Kleborate
genotype calls for public *K. pneumoniae* genomes, of the kind published in:

> Lam MMC, Wick RR, Watts SC, Cerdeira LT, Wyres KL, Holt KE (2021).
> *A genomic surveillance framework and genotyping tool for Klebsiella
> pneumoniae and its related species complex.*
> **Nature Communications** 12:4188. doi:[10.1038/s41467-021-24448-3](https://doi.org/10.1038/s41467-021-24448-3)

The `virulence_score` (0–5, ordered *ybt* < *clb* < *iuc*) and
`resistance_score` (0–3, ESBL < carbapenemase < carbapenemase + colistin)
columns are Kleborate's own summary scores as defined in that paper, and the
`Convergence_event` labels correspond to the AMR–virulence convergence events
enumerated there.

Nature Communications articles are published under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), which permits
redistribution of the associated data with attribution. **The MIT licence in
this repository covers the source code, tests and documentation only. The files
in `examples/klebsiella/data/` are redistributed under CC BY 4.0 and must be
attributed to Lam et al. (2021) in any downstream use.**

## Reproducibility gaps — read this before citing the case study

The following are **not** currently recorded, and cannot be reconstructed from
the files themselves. They must be supplied before this case study is used as
anything other than an illustration of the software:

1. **Accessions are partial, not absent.** 1261 of the 1500 `Strain_ID` values
   (84.1 %) are resolvable public accessions - 1074 RefSeq/GenBank assembly
   accessions (`GCF_003830175.1`) and 187 ENA/SRA run accessions
   (`ERR560516`). The remaining **239** are internal identifiers
   (`16703_5_23`, `Klebsiella_pneumoniae_INF145`) that cannot be resolved to a
   public record from these files. Those 239 isolates cannot be re-derived from
   primary data; the other 1261 can, and the accession list should be shipped
   as a TSV alongside these files.
2. **Kleborate and database versions.** The Kleborate release and the versions
   of its AMR/virulence/capsule databases determine the column set and the
   calls. They are not recorded.
3. **Subsetting procedure.** An earlier release of this repository shipped
   metadata for 8254 isolates alongside layers for 1500, with no script
   connecting them. `subset.py` in this directory reproduces the 1500-isolate
   selection *from a supplied full cohort*, but the selection actually used to
   produce the committed files was not recorded and `subset.py` is therefore a
   specification of intent, not a replay of history.

What *is* verifiable from the shipped files, and is reported here so a reader
can judge representativeness:

* 1500 isolates in **768 distinct sequence types**; the largest ST contributes
  8 isolates (0.53 %). The subset is heavily de-replicated relative to a raw
  surveillance collection, which suppresses — but, as the lineage diagnostics
  in the case study show, does not remove — clonal confounding.
* 78 K loci and 12 O antigens are represented.
* Prevalence of the acquired AMR classes ranges from 0 (tigecycline,
  ESBL+inhibitor-resistant β-lactamase) to 0.49 (aminoglycoside).
* 39 isolates carry a named AMR–virulence convergence label.

## How the case study uses these data

`config.yaml` treats AMR and virulence as `wide_binary` layers and K locus as a
`one_hot` layer, collapses each virulence locus (*ybt*, *clb*, *iuc*, *iro*,
*rmp*) to a single presence call — the biological unit is the mobile element,
not the gene, and the uncollapsed encoding weights each locus by its gene count
— and keeps every AMR class regardless of prevalence, because the rare
determinant is the surveillance signal.

`metadata.csv` is **not** used for clustering. It supplies the lineage column
for the population-structure diagnostics and the external columns
(`virulence_score`, `resistance_score`, `K_locus`) that the discovered partition
is scored against, so that the reader can see how much of the result was
already available from one column of the input.

## Regenerating

```bash
# from a full Kleborate output table (not shipped):
python examples/klebsiella/subset.py --kleborate FULL.tsv --n 1500 --seed 42 \
    --out examples/klebsiella/data/
```
