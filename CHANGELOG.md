# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/);
versioning: [SemVer](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-09-10

First public release. Everything below describes what the package contains;
where a paragraph mentions an earlier state, that state was never released.

### Added

**The run.** A configuration names a susceptibility table and a metadata
table with the lineage column, and nothing else. The run reads the cohort for
the share its lineages carry and stops: the per-antimicrobial
lineage-membership share against its permuted control, the e-value and the
sequential e-process over the declared intakes, the reading of a recorded
dilution, and the decomposition of a prevalence difference, each behind its own
gate. The record is `clonal_share_result.json`, and both reports are rendered
from it. `examples/ssuis/config.yaml` reads the shipped cohort this way.

**Input check.** Every run starts by reading the phenotype table, joining the
metadata and the dilutions, and writing `input_qc.json` with a plain-language
`input_qc.md`: per antimicrobial the count of the rarer outcome, the join,
and per lineage the group sizes, the share of isolates in lineages of at
least two, and whether the clonal-share estimator will accept the cohort.
`--check-input` runs this step alone. A refused input leaves with exit code 2
and a sentence naming the position, not a traceback.

**Run report, written and drawn.** `report.md` and `report.html` beside
`clonal_share_result.json`, both rendered from one structure in a fixed order
of seven sections. The first page is the measurement summary: the measurand,
the gates that fired, one row per trait with share, interval, control,
e-value and a reading fixed by two conditions the record holds, what was not
evaluated on this run and why, and an interpretation generated from that
table by fixed rules. The pages after it defend the first: the admissibility
conditions with observed value, threshold and verdict, the lineage sizes
behind the support gate, the intervals drawn with a calibration line from the
release's validation grid and the species interval beneath each, the
e-values against their two thresholds and the running product over intakes,
the reading at recorded resolution beside the share of the call, and the two
prevalences with the concentration of carriage, six figures in all. Every value
on either page is read from the record and none is recomputed; both carry the
SHA-256 of the record; the drawn page is one file with no script and nothing
to fetch.

**Three shares, named apart.** One collection admits three quantities and the
package reports them as three. The lineage-membership share
`eta_c = B_w / (B_w + sigma^2)`, with `B_w` the isolate-weighted variance of
the lineage means, is what a rule scored out of sample reaches and is the
target of `clonal_share`. The realised share `rho_c` of equation (1) is the
same numerator on the scale of a variance component,
`S_a^2 = B_w / (1 - sum w_g^2)`, and is the target of `realised_share`, whose
noncentral-F inversion is exact for it. The superpopulation share `rho` is
what a fresh draw of lineages would show. The first two differ by the design
factor alone, `(G - 1) / G` at equal lineage sizes, so they part company where
lineages are few; the estimator grid carries all three truths per cell and
scores each estimator against its own.

**Two intervals for the clonal share.** `clonal_share` reports its lineage
bootstrap as the interval for the lineage-membership share
(`interval_target`), and a second interval, `superpopulation_low` and
`superpopulation_high`, for the share a fresh draw of lineages from the
species would show: a within-lineage bootstrap under a chi-square layer on
the number of lineages less one. The chi-square law is written for the
design-corrected component ratio, so the estimate and its within-lineage
draws are carried onto that scale by the design factor `1 - sum w_g^2` on the
odds scale before the law is inverted, the species interval is stated on
that scale, and it is widened to the envelope of the lineage bootstrap taken
on the same scale. On the estimator grid the first covers its
target at 0.92 to 1.00 from ten lineages to three hundred and the second
covers the fresh-draw target at 0.94 to 1.00 with ten lineages and 0.89 to
1.00 with five, binary and continuous traits alike; the lineage bootstrap
alone covers the fresh-draw target at 0.71 to 0.96 with ten lineages, 0.82
on average. Below ten
lineages the record flags the collection (`few_lineages`) and the report
names the species interval as the one to read; for a single trait the record
carries the largest single-lineage share of the between-lineage variation
(`dominant_lineage_share`), the diagnostic of the carrier structure under
which the species interval, like the variance component and the mixed model,
loses its level with few lineages while the interval for the lineages in
hand does not. The report prints both in its Table 3, and
every number computed before the second interval is unchanged by it. The
grid runs from five lineages to a thousand and from three isolates a lineage
to twenty, holds cells whose lineage effects follow a two-point carrier law
rather than the normal law the chi-square layer assumes, and a binary arm at
eight per cent prevalence, so that each interval's level is measured at both
ends of the range and against the violation of its own assumption.

**The share on the latent scale.** For a single binary trait `clonal_share`
restates its estimate and both intervals on the latent scale of a threshold
model (`latent_share`, `latent_low`, `latent_high`, `latent_species_low`,
`latent_species_high`; `latent.py`), the scale on which a mixed model with a
binomial link reports an intraclass correlation. The map is the exact
liability-scale transformation, monotone at a given prevalence, so the
intervals keep their measured coverage; the report prints the latent figure
as the last column of its Table 3. Every interval remains measured on the
observed scale, which is the scale a prevalence table uses.

**The floor of the permutation p-value.** `clonal_share` carries
`p_floor = 1 / (n_perm + 1)`, and the report prints a p-value at the floor as
the floor of that many permutations rather than as a number.

**The S. suis cohort against its determinants and at three lineage
definitions.** `benchmarks/ssuis_mechanism.py` runs the per-agent share, the
share of each mechanism block, and the same share on carriage of the
determinants the source study calls for each agent, a layer the estimate
never sees. `benchmarks/ssuis_resolution.py` runs the per-agent share with
the lineage defined as the hierBAPS cluster, the serotype and the sequence
type, and recomputes the ordering at each. Each agent's determinants are
placed as the source study's determinant table places them, chromosomal,
on a mobile element, or on both with the chromosome carrying most, and the
ordering is reported as pair counts between the classes and between the
seven determinant sets the panel holds, with no test attached, since agents
scored on one determinant set are one observation.

**Anytime-valid evidence over a programme of looks.**
`evalues.sequential_e_process` takes the intakes of a programme as separate
batches in arrival order, fits the lineage probabilities on every earlier
batch, scores the new batch against the null maximised on it, and multiplies
the factors. Each factor has conditional expectation at most one given the
past, so the running product is a test supermartingale: it may be inspected
after any intake and stopped at any time, and e-BH on the running products of
a panel holds its level at whatever intake a programme stops. `e_process`
remains the e-value of one look, and its documentation says what
recomputing it on a growing cohort does and does not guarantee: a valid
e-value at each look, no guarantee over the sequence. A run given
`dataset.batch_column`, the metadata column naming the intake an isolate
arrived in, carries the programme in the record: the batches in ascending
order of that column, the running product per trait, and a second e-BH
decision taken on the log scale by `evalues.e_bh_log`, so that a product over
many intakes does not overflow before the decision is read. The repeated-looks
benchmark scores the sequential rule on the same simulated panels as the
others (`ebh_sequential`; `--evalue-only` reruns that arm alone on the same
seed stream).

**The dilution panel read on its lattice.** `intervals_from_mic` takes the
tested wells to be every doubling from the lowest to the highest recorded
value when the readings sit on a doubling lattice, so that a dilution no
isolate landed on still bounds its neighbours and a rounding variant of one
dilution, 0.06 beside 0.064, is not read as a well of its own; readings off
any such lattice are their own wells, and a laboratory that recorded its panel
passes it as `wells`. `PanelGeometry` reports the reading made (`lattice`),
the wells the panel is taken to have (`n_wells`) and the distinct values
recorded (`n_wells_recorded`). On the shipped cohort this folds rounding
variants on seven agents and fills unoccupied wells on six; the penicillin
share read from the dilution moves from 0.655 to 0.679 with 73 rather than
66 per cent of its readings censored, and the agreement of the dilution with
the call rises from Pearson 0.93 to 0.94 over the thirteen paired agents. The
calibration, the design grid and the cut-off sensitivity envelope were run
again on the corrected reading.

**The conditional moments of a truncated normal in the upper tail.** The
moments of a normal restricted to an interval are computed by reflection when
the interval lies above the mean, because `log_ndtr` saturates at zero for
large positive arguments. Before the reflection an interval more than about
seven residual standard deviations above the lineage mean returned the mean
itself as its conditional expectation, so an isolate far beyond the top of the
panel was read as sitting at its lineage mean. On the shipped cohort the
amoxicillin share read from the dilution moves from 0.758 to 0.713, the one
agent whose readings reach that far above their lineage means; the
calibration of the interval, run again on the corrected reading, is unchanged
within Monte-Carlo error. Both tails are now checked against an independent
implementation.

**Untyped isolates are set aside.** No estimator places isolates without a
lineage label in a level of their own: an untyped isolate is not a member of
any lineage, and a level made of every isolate whose typing failed would
measure the typing process. `clonal_share`, `realised_share`,
`censored_clonal_share`, the e-value and the sequential e-process all set
such isolates aside, and each records `n_dropped_untyped`; `clonal_share`
also records `missing_share`. No cohort in the article carried an untyped
isolate at the resolution it was read, so no number in the article moves.

**The nonshared contribution of a prevalence difference.**
`decompose_prevalence_difference` reports `nonshared`, the part of the
difference carried by lineages seen in one collection only, and
`composition_shared`, the rest of the composition term, so that a composition
finding can be read for how much of it rests on lineages the other collection
never held. The identity, the intervals and the gate are unchanged.

**The estimator grid's binary intercept.** The benchmark's binary cohorts
place the liability intercept at the value that makes the prevalence marginal
over the lineage effects equal to the cell's nominal prevalence; the first
release placed the quantile of the nominal prevalence on the liability
directly, which gave a marginal prevalence that rose with the share, 0.32
rather than 0.25 at a share of 0.5. The seventy binary cells were run again;
the continuous cells are bit-identical. On the corrected cells `clonal_share`
covers 0.965 at 25 per cent prevalence where it covered 0.970, the realised
interval is estimable in 0.725 of runs where it was in 0.819, and at 8 per
cent prevalence it is refused in every run; `validation_grid.json` carries
the new figures.

**Lineage attribution.** A per-antimicrobial clonal share estimated out of
sample so that a lineage variable is not paid for its levels, with a
cross-validation debias, a cluster bootstrap interval that draws lineages
whole, and a support gate whose threshold is read off a coverage curve. A
cohort with a single lineage, like a cohort with no variance in the trait,
returns a non-estimable result with the reason named, not a share of zero.

**An exact permutation p-value beside the clonal share.** The reported share is
a mean over `repeats` fold draws, so the permutation test compares the mean of
the first `null_repeats` observed draws with the mean of `null_repeats` draws
under each permutation: the same function of the labels on both sides, which
is what makes the test exact under exchangeability (`null_repeats=5` by
default). The extra draws come from a second random stream taken after every
recorded quantity, so the estimate, its interval and the penalty do not depend
on the setting. `benchmarks/null_uniformity.py` draws cohorts with no lineage
effect and checks that the p-value is uniform on the permutation grid, then
draws cohorts with a graded lineage effect and records the power curve; the
check runs in the test suite on a small design and its full record is shipped
under `benchmarks/results_null_uniformity_2026-09-04/`.

**Interval-censored panel reading.** The same share from a dichotomised
non-wild-type call, a recorded minimum inhibitory concentration or a censored
reading, which enter one interval likelihood at three widths. Panel geometry is
reported before any share, an end-well reading is treated as censored under a
stated coarsening assumption with the alternative computed beside it, a
recorded censoring operator overrides that assumption, and a cohort with more
than half its isolates in lineages wholly beyond the panel is refused. Three
widths are given and only one is an interval: the variance-ratio inversion in
`ci_low` to `ci_high`, with a cluster bootstrap and a profile likelihood from
the marginal likelihood beside it as checks that answer different questions;
`profile_interval` cuts the profile at the level `alpha` names.

**Anytime-valid evidence.** An e-value per antimicrobial, in the betting sense,
built as the split likelihood ratio of universal inference on the same folds
the share already uses, with e-BH across the panel controlling the false
discovery rate under arbitrary dependence. A surveillance programme may inspect
the running value at every intake without spending an error budget.

**Prevalence decomposition.** A Kitagawa split of a prevalence difference into
a lineage-composition component and a within-lineage rate component, with an
exact residual, Benjamini-Yekutieli control within each component family, the
effective number of independent agents, and two estimability gates: shared
lineage support, and the pair of conditions under which missing lineage labels
make the two collections incomparable.

**Surveillance readings that need no contrast.** Prevalence per isolate and per
lineage, the effective number of lineages carrying a feature, and the direction
of a departure from proportional carriage.

**Contract and reproducibility.** A strict YAML configuration that refuses an
unknown key rather than silently disabling the analysis it names, and refuses
a value of the wrong type, a permutation budget of zero, or a contrast or
intake column the metadata does not hold; one record under public output
schema 1.0 carrying the configuration as run (`config`) and the versions of
the package and its numerical stack (`versions`), with the realised share and
the verdict of its kurtosis gate per trait (`realised_share`) beside the
clonal share, and two reports rendered from it; exit codes 0
and 2 separating a completed run, in which a withheld estimate is a result,
from a configuration or input error; a master seed spawning every stochastic
stage; and two verifiers, `scripts/check_evidence_accounting.py` for the
atlas evidence trees, which checks that every agent of every cohort leaves
exactly one row with a refusal from the closed set, and
`scripts/check_cited_hashes.py`, which checks that every digest a document
quotes still matches the file it names; both run on the evidence deposit
rather than in continuous integration, since the deposit is not in the
repository.

**Property tests.** `tests/test_properties.py` states the identities the
arithmetic must satisfy for every input and lets Hypothesis search for a
counterexample: the two components of the prevalence decomposition sum to the
difference exactly, checked against rational arithmetic; the exact binomial
interval brackets the proportion and moves with the count; the permutation
p-value is the Phipson-Smyth ratio; Benjamini-Hochberg matches its step-up
definition; and a lineage label is a label, so renaming the lineages changes
nothing.

**Evidence accounting.** `benchmarks/agent_screen.py` is the one screening rule
for both atlases, with a closed set of refusal reasons; `reconcile()` raises
when a cohort cannot account for every agent it holds. The cross-species and
veterinary atlas records account for every agent present in every cohort, name
the reason for every refusal from the closed set, and carry the estimator's own
reason for a refused realised interval. Two verifiers,
`scripts/check_evidence_accounting.py` and `scripts/check_cited_hashes.py`,
hold the evidence deposit to that accounting and to the digests the
documents quote; they run on the deposit rather than in continuous
integration.

**Engineering gates.** `ruff` and `mypy` run in CI and both pass over the
package; the coverage gate is per module (`scripts/check_coverage_floors.py`);
the documentation builds strictly in CI and carries an API reference generated
from the docstrings; releases are published from GitHub releases through PyPI
Trusted Publishing with PEP 740 attestations and no stored token; actions are
pinned to commit digests; mutation testing is configured for the estimator
modules and its score is recorded in `CONTRIBUTING.md`.

**Supported stack.** Python 3.11 or later, NumPy 2.0, pandas 2.2 and SciPy
1.13 or later, following the scientific-Python support window; CI runs 3.11
to 3.13. Every reported number was produced inside these floors.

**Examples.** Two real cohorts. 677 *Streptococcus suis* isolates with
thirteen antimicrobials called against an epidemiological cut-off, sixteen of
recorded concentrations, a hierBAPS population cluster for every isolate and
the yearly intakes, read at population-cluster resolution and, with two
countries declared as the collections to contrast, at sequence-type
resolution. The poultry-meat *Salmonella* cell of one NCBI Pathogen Detection
release, 7,049 isolates with the calls the release records for 22 agents, cut
from the pinned release by `examples/salmonella_poultry/build_cell.py` under
the host and matrix rules of `benchmarks/vet_source_taxonomy.py`, read at
serovar and at SNP-cluster resolution from the same two tables, with the
expected outputs of both runs. And a planted control the continuous
integration builds and runs.

**Inputs that cannot be read as they stand.** `censored_clonal_share` and
`profile_interval` set aside isolates with no lineage label and count them
(`n_dropped_untyped`) instead of coding the missing label as a lineage of its
own, which is the policy `clonal_share` and `realised_share` already applied,
and the input check reads the support on the typed isolates alone; a dilution
value written with its censoring sign (`<=0.25`, `>64`) is read as the
operator it carries, a value that cannot be read as a number is counted
(`n_unparseable_values`) and named in the input check rather than dropped in
silence, a repeated isolate identifier in the metadata is counted
(`metadata_join`) and warned about with the first row kept, and a table that
cannot be decoded, a join that matches no isolate or a contrast level with no
isolates is refused with exit code 2 and a sentence; a within-lineage
component whose shared-support gate refused it is not counted as a discovery
whatever its interval says, nor is either component where the lineage labels
are missing informatively between the two collections
(`composition_estimable`), and the decomposition family reports the smallest
attainable q of its bootstrap beside the target;
an interval whose lower bound exceeds its upper is refused as an input error
rather than fitted; and a set of readings that are all one value is refused
with the reason that there is no variance to split. `evalues.e_bh` ranks
an infinite e-value first rather than treating it as zero, and delegates to
`e_bh_log`, which takes the whole decision on the log scale so that evidence
strong enough to overflow a float is still evidence.

**The veterinary reading of the same source.** `benchmarks/vet_atlas.py` cuts
the Pathogen Detection isolates of one organism by host and sampling matrix
(`benchmarks/vet_source_taxonomy.py`) and runs every cell at two definitions
of a lineage, the SNP cluster and the serovar, with the permuted-label control
at both (`benchmarks/resolution_atlas.py`); `benchmarks/period_split.py`
splits one poultry cell at a year fixed in advance and reports the share on
each side of the cut with its control, and `benchmarks/verify_vet_claims.py`
checks the counts the split rests on before it is run.

**A snapshot the cross-species atlas can be rerun on.**
`benchmarks/pathogen_detection_releases.tsv` records the Pathogen Detection
release accession per organism, and `fetch_pathogen_detection.sh` fetches that
release rather than the current one, so the retrieval reads the lineage
variable the published atlas was written against. `atlas_cross_species.py`
carries the support gate, the support and the refusal reason into every row of
`atlas_table.json` and summarises two layers: the cells the gate admits, which
is what a claim about a species rests on, and every cell that reached the
estimator, which is what the permuted-label control is measured on.

### Known limits

- The clonal share and the censored share are different estimands, one
  predictive skill and one population variance, and are reported as such.
- A single cut point does not identify the latent scale. The censored share is
  then reported on the liability scale by convention, and the calibration in
  `benchmarks/censored_calibration.py` gives its bias as a function of the
  prevalence at the cut.
- The shipped examples are two cohorts; the cross-species reading of the
  article is in its evidence deposit, not in the repository.
