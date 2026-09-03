# Novelty evidence

This file supports the priority words used in the article and in the source
docstrings. Each claim names the search behind it. A search is a measurement
with a date on it, so the queries are in `docs/novelty_search.py`, the returned
records are in `docs/novelty_search_results.json`, and re-running the script
will not reproduce the counts because the corpus grows. That is the point: a
reader who re-runs it sees what has appeared since 2026-09-02.

## Claim 1

> No published work provides a packaged estimator that splits an observed
> difference in antimicrobial non-wild-type prevalence between two bacterial
> collections into a lineage-composition component and a within-lineage rate
> component, with a bootstrap interval on each component, a stated convention
> for lineages seen in only one collection, and estimability gates whose
> thresholds are measured.

**Withdrawn, and what replaced it.** An earlier form of this claim said that
no published work makes the split at all. That was false, and it is corrected
here rather than quietly narrowed. Plumb and colleagues (Vaccine
2020;38(27):4273-4280, doi:10.1016/j.vaccine.2020.04.048) split a rise in
pneumococcal non-susceptibility in Alaska between 2008 and 2015 into a
serotype-redistribution part and a within-serotype part, by comparing the
observed data with modelled data removing either factor, and report that the
overall three-point rise came from within-serotype increases that outweighed a
protective redistribution which alone would have produced a three-point fall.
That is direct standardisation, it is this decomposition, and it is published.

**Why the search could not have found it.** The queries behind the earlier
claim screened for the estimator under the names it carries in demography and
economics: Kitagawa, Oaxaca, Blinder, Fairlie, shift-share, composition effect.
The epidemiological literature that performs the same arithmetic calls it
standardisation, serotype replacement, serotype redistribution or clonal
replacement, and none of those words was in the screen. A search that can
return only the vocabulary it was given is not a search, and the absence it
produced was a property of the query rather than of the corpus.

**The search as it now stands.** A fourth query family was added on
2026-09-02 using the epidemiological vocabulary and is recorded under
`europepmc_epidemiological` in `docs/novelty_search_results.json`:
5 queries, 250 titles screened, 50 carrying both a
standardisation or replacement term and a microbiology term. Reading those
50 shows almost all of them to be surveillance reports that describe
serotype distribution and resistance side by side without splitting either.
One further fact belongs here: the screen still does not return Plumb and
colleagues, whose title names neither the method nor the split. A reader found
it. A title screen cannot find a paper that does not announce itself in its
title, and stating that limit is worth more than a count that implies the
search was exhaustive.

**What survives, and the evidence for it.** Thirteen queries across Europe PMC (which indexes PubMed, PubMed
Central, Agricola and the bioRxiv and medRxiv preprint servers) and Crossref,
run on 2026-09-02. Every returned title was screened for carrying both a
decomposition term (Kitagawa, Oaxaca, Blinder, Fairlie, shift-share,
"components of a difference", "composition effect") and a microbiology term
(bacterial, antimicrobial, lineage, sequence type, serotype, clone, isolate).

| system | queries | titles screened | carrying both |
|---|---|---|---|
| Europe PMC | 10 | 464 | 0 |
| Crossref | 3 | 150 | 0 |

A Crossref hit count is not evidence and is not quoted as such: the API scores
an implicit OR of every term, and these three queries report 1.2, 6.3 and 1.9
million hits while none of the returned records is on the subject. The screen
is the measurement.

The surviving claim is about the packaged estimator rather than about the
arithmetic. What the located precedent does not carry, and what the package
adds, is: an interval on each of the two components rather than on a trend
within one of them; a stated convention for lineages present in only one
collection, which direct standardisation leaves to the analyst; a shared-support
gate that refuses the within-lineage component when the two collections hold too
few lineages in common, with the threshold read off a coverage curve; and a
label-completeness gate that refuses the comparison when lineage typing is
missing differentially between the collections. Sections S2b and S3 of the
article supplement report all four.

A separate direct query on PubMed on the same date,

```
(Kitagawa decomposition OR Oaxaca-Blinder OR "composition effect")
AND (antimicrobial resistance OR antibiotic resistance)
AND (lineage OR clone OR serotype OR "sequence type")
```

returned **0 records**.

## Claim 2

> The Kitagawa decomposition and its descendants have not been carried into
> bacterial population data.

**Evidence.** Query Q7, `"Kitagawa decomposition" OR "Kitagawa-Blinder-Oaxaca"`,
returns 75 records in Europe PMC. The 50 returned are demographic and
health-services work: survival after dementia diagnosis, caesarean delivery
rates, hospitalisation risk under population ageing, fertility, disparities in
medication management. None concerns a bacterial population, a pathogen, or
antimicrobial susceptibility.

## Claim 3

> A pass-or-fail gate built on a same-lineage pair co-assignment permutation
> test cannot fail on a cohort large enough to be worth analysing, because the
> statistic is computed over O(n^2) pairs while the structure it measures is
> not a function of n.

**Evidence, and the part of the claim that was withdrawn.** The sample
dependence of pairwise co-assignment statistics is documented, and an earlier
draft of the article that called it unnamed was wrong. Pinto, Melo-Cristino and
Ramirez (PLoS ONE 2008;3:e3696) state the sample dependence of the Wallace
coefficient and propose a confidence interval and an independence baseline;
Severiano, Pinto, Ramirez and Carrico (J Clin Microbiol 2011;49:3997) add the
adjusted coefficient with intervals; Kelly and colleagues (Bioinformatics
2015;31:2461) build a power and sample-size framework for permutation tests on
pairwise distances. Harmon and Glor (Evolution 2010;64:2173) and Guillot and
Rousset (Methods Ecol Evol 2013;4:336) make the general point for
distance-matrix permutation tests.

What no located paper does is apply the critique to such a test used as a
**gate** on a candidate trait partition, or state the O(n^2) pair-count framing
in these terms. The article claims that and nothing more, and quotes the
subsampling series that measures it on the shipped cohort: z = 1.5, 7.2, 10.7,
19.0, 25.5 at n = 67, 135, 270, 473, 677, against an attributable share of
0.139, 0.209, 0.299, 0.298, 0.296 on the same partition and the same labels.

## Claim 4

> The two-player Shapley decomposition of the variance-explained game against a
> lineage label has not been applied in bacterial population genomics.

**Evidence, and the part of the claim that was withdrawn.** Dividing bacterial
phenotype variance between lineage and other effects is established and an
earlier draft that implied otherwise was wrong. Lees and colleagues (eLife
2017;6:e26255) partition pneumococcal carriage duration into genomic variation,
serotype, drug resistance and other locus effects with numbers attached, and
MacFadden and colleagues (J Clin Microbiol 2019;57:e01780-18) compare
sequence-type-only against resistance-locus prediction. What was not located is
the Shapley form of that split, which is what the article claims, together with
its application to a routine susceptibility panel and its use as a magnitude
where a significance test stood.

A terminology collision is pre-empted on first use rather than left to a
reader: the Shapley value already names a diversity index on phylogenies
equivalent to Fair Proportion (Haake, Kashiwada and Su, J Math Biol
2008;56:479). The game decomposed here is the variance explained, not the tree.
These are also not SHAP attributions, which are Shapley values of a prediction
with respect to input features.

## Claim 5

> No published work reports, for a bacterial collection, both the share of a
> trait's variance that a fresh draw of lineages would show and the share the
> lineages in hand do, from one fit and with an interval on each.

**Evidence.** 7 queries, 5 on Europe PMC and
2 on Crossref, run on 2026-09-02 and recorded in
`docs/novelty_search_results.json` under `europepmc_estimand` and
`crossref_estimand`. They are kept apart from the queries behind Claims 1 and 2
so that adding them could not move the counts already reported there. Titles
were screened for carrying both a conditional-share term (noncentral F,
noncentrality, noncentral chi-square, finite population, narrow or broad
inference, superpopulation, conditional inference, intraclass correlation,
variance component, proportion of variance, variance explained, heritability)
and a microbiology term.

| system | queries | titles screened | carrying both |
|---|---|---|---|
| Europe PMC | 5 | 250 | 1 |
| Crossref | 2 | 100 | 0 |

The one Europe PMC title carrying both is a terminology collision and is
recorded rather than ignored. Hafez and Abbas (Insects 2026;17(5):480,
doi:10.3390/insects17050480) estimate the *realized heritability* of
diflubenzuron resistance in *Musca domestica* over forty-two generations of
selection. Realized heritability in that sense is the response to selection
divided by the selection differential, a quantity of an experiment that
selected, and it is neither an intraclass correlation nor a statement about the
groups a collection holds. The word is the same and the estimand is not.

The two searches that would have refuted the claim outright returned nothing:
E2 asks for a noncentral F or noncentrality interval used on bacterial data at
all, and E3 asks for the two estimands named together anywhere in the life
sciences. Neither returned a title coupling the method with the subject.

**What this claim is not.** The arithmetic is old and is cited in the article
and in `docs/methodology.md` section 13. Bounding a noncentrality parameter by
inverting the noncentral distribution in it is Venables (J R Stat Soc B
1975;37:406, doi:10.1111/j.2517-6161.1975.tb01554.x). Carrying that inversion
to an exact interval for a fixed-effects analysis-of-variance effect size is
Steiger (Psychol Methods 2004;9:164, doi:10.1037/1082-989X.9.2.164). Treating
the group effects a design holds as a finite population with a variance of its
own, rather than as a draw from an infinite one, is Cornfield and Tukey (Ann
Math Stat 1956;27:907, doi:10.1214/aoms/1177728067). The contribution is the
pairing and the gate: the two shares reported side by side from one fit of a
bacterial cohort, with the exact one refused where the residual kurtosis
exceeds the limit its own coverage curve set.

## What the nearest existing work does

The claim is that the estimator is absent, not that the question is
unrecognised. It is recognised, and answered by other means.

**Phylodynamic fitness estimation.** The closest published approach fits a
hierarchical Bayesian phylodynamic model to a dated phylogeny and estimates the
growth-rate effect of each resistance determinant as antimicrobial use changes,
demonstrated on a twenty-year *Neisseria gonorrhoeae* collection with national
treatment data (Helekal and colleagues, Nature Microbiology 2026;11:375,
doi:10.1038/s41564-025-02235-w). It answers a deeper question than
the decomposition does and needs a dated whole-genome phylogeny, a determinant
catalogue and a usage series to do it. The decomposition needs a lineage label
and a susceptibility result, which is what a diagnostic laboratory has.

**Qualitative reading off a phylogeny.** Genomic surveillance reports routinely
distinguish clonal expansion of a resistant lineage from independent
acquisition, and do so by inspecting a tree. That reading is a statement about
mechanism made without an estimate or an interval, and it is not available at
all to a laboratory that types by MLST rather than by sequencing.

**Regression adjustment for lineage.** A logistic model with lineage terms and
a collection indicator is the obvious alternative and is the one a reader will
propose. It is implemented and measured against the decomposition in
`benchmarks/decomposition_vs_regression.py` rather than dismissed.

## What is not claimed

The identity is Kitagawa's, published in 1955, and no part of the arithmetic is
new. The same holds for the second estimand: the noncentral inversion and the
finite-population reading of a variance component are both old, and are cited
to their sources under Claim 5. The Herfindahl-Hirschman index is the competition-authority concentration
measure and the departure statistic is the Theil form. The contribution is the
transfer: these estimators applied to a bacterial lineage variable, with the
estimability conditions that a bacterial cohort violates made explicit, gated
and measured.

Three priority words were withdrawn on the evidence recorded under Claims 1, 3
and 4, and the surviving forms are narrower than the ones first drafted. A
claim that a literature search cuts down is a claim the search was worth
running. Claim 1 was cut down not by its own search but by a reader, and the
search was then repaired so that it could have found what the reader did.

Priority words are used only where a row above supports them. The article does
not claim novelty for similarity network fusion, consensus clustering,
post-clustering selective inference, latent-trait modelling, or any component
of the diagnostic pipeline, all of which are cited to their sources.
