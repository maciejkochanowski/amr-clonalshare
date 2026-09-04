# amr-clonalshare

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![Python](https://img.shields.io/badge/python-3.11%2B-blue.svg)
![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22306354.svg)](https://doi.org/10.5281/zenodo.22306354)

**Release status:** public product version `1.0.0`.

**How much of the antimicrobial resistance in a collection travels with the
clone?** Give it a lineage label per isolate (sequence type, core-genome or
SNP cluster, or serovar) and a susceptibility call or a recorded dilution, and
it returns the share of the resistance variance the lineages carry, per
antimicrobial, as two estimands with their own intervals: a realised share for
the collection in hand and a superpopulation share for the species it was
drawn from. Each sits behind a gate that refuses where the collection cannot
support an estimate, each carries an anytime-valid e-value for panels that
surveillance re-reads every year, and a prevalence difference between two
collections is split into lineage mix and within-lineage rate. The section
[How much of this is the clone?](#how-much-of-this-is-the-clone) is where to
start.

The package also fuses one or more wide binary matrices per isolate (acquired
AMR determinants, virulence loci, capsule type, plasmid replicons) into a
similarity network, clusters it, and profiles the result. So will a dozen other
tools. The reason to use this one is that it also answers the questions that
decide whether the partition means anything:

- **Is the largest cluster just the isolates with no detected features?**
  Jaccard is undefined for two all-zero rows and the usual convention calls them
  identical. On real accessory-genome data that block is often the majority of
  the cohort. Reported as `empty_stratum_ari`, with a gate.
- **Would a different row order in the CSV give a different answer?** With
  binary data the *K*-nearest-neighbour graph is cut in the middle of a huge tie
  block unless neighbour selection is made order-invariant, which it is here,
  and `tie_inflation` tells you how close to degenerate the graph is.
- **Did the fusion actually fuse, or did one layer win?** Leave-one-layer-out
  influence, an effective number of contributing layers, and a `regime` label
  (`complementary` / `collapsed` / `redundant` / `unstructured` / `uninformative`),
  because `n_eff` alone gives the same answer to three concordant layers and to
  three noise layers.
- **Are these clusters just clones?** Same-lineage pair concordance against a
  permutation null. Adjusted mutual information will tell you there is no clonal
  confounding when there is; the pair statistic will not.
- **Is the structure real, and is it *discrete*?** Two different questions,
  answered by two different tests, because a valid post-clustering p-value
  rejects on a smooth gradient just as readily as on genuine archetypes. The
  discreteness null is a *q*-factor latent trait with *q* chosen by BIC; pinned
  at one factor it calls any two-dimensional gradient discrete 100 % of the time.
  The default fits dimensions 1--3, conditionally extends to 4 when the
  provisional optimum is on the boundary, and returns a third, non-estimable
  state if the optimum remains on the final boundary.
- **Does the partition beat one column of the input?** If you can supply measured
  phenotype, that is the only criterion here that is not a function of what was
  clustered. On the shipped *Klebsiella* cohort the answer is no, and the CLI
  exits non-zero saying so.

If any of these trips, the CLI says so on stderr and exits non-zero.

## What the controls and real cohorts show

Two independently seeded planted controls provide positive examples that do
not end in refusal. In each, the pipeline recovers all three planted groups
(ARI 1.000), detects held-out structure, rejects the fitted two-factor
continuum (0/99 exceedances; corrected `p = 0.01`; exact tail-probability 95%
CI 0-0.0366), validates 60 defining features, finds no lineage confounding,
reaches claim level 4 (`archetype_candidate`), and exits 0.

On the shipped *Klebsiella* cohort a discrete signal is detected beyond the
BIC-selected three-factor continuum (`p = 0.01`); with the BIC sweep run to
four factors, `q = 3` is an interior optimum. This statistical result
does not establish biological archetypes: the two-group partition remains
strongly lineage-associated, correlates with the all-zero virulence stratum,
and is inferior to a one-dimensional resistance score for measured
susceptibility (mean balanced accuracy 0.790 versus 0.695; bootstrap difference
+0.090, 95 % CI [+0.048, +0.125]). The CLI therefore exits 3 while preserving
the positive discreteness result. The *S. suis* example likewise detects a
discrete signal beyond a three-factor continuum (`p = 0.01`) but exits 3 because
the five-group partition is lineage-confounded. These are independent biological
interpretation gates, not refusals by the continuum test.

## Install

From PyPI, or from an unpacked source tree:

```bash
python -m pip install amr-clonalshare       # released version from PyPI
python -m pip install .                     # core, from source
python -m pip install ".[phylo]"            # + the phylogenetic convergence test
python -m pip install -e ".[dev]"           # editable source + test suite
```

NumPy, pandas, SciPy, scikit-learn, PyYAML and pandera. `dendropy` is optional
and only needed for `amr_clonalshare.archephy`.

## Quickstart

Run the controls first. They are the reason the real-cohort diagnostics are
interpretable.

```bash
python examples/synthetic/make_data.py

# positive control: 3 planted archetypes
amr-clonalshare --config examples/synthetic/planted.yaml --results-dir out/planted

# independent positive confirmation: a different random seed, same design
amr-clonalshare --config examples/synthetic/planted_confirmation.yaml --results-dir out/planted-confirmation

# negative control: no clusters at all -> selected_k = 1
amr-clonalshare --config examples/synthetic/null.yaml    --results-dir out/null

# adversarial control: all structure is clonal -> flagged as lineage-confounded
amr-clonalshare --config examples/synthetic/clonal.yaml  --results-dir out/clonal

# the case study: 1500 Klebsiella pneumoniae genomes
amr-clonalshare --config examples/klebsiella/config.yaml --results-dir out/kp
```

Each run writes `summary.json` (headline numbers, claim level and active gates),
`cluster_result.json` (everything), `archetype_profiles.tsv` (per-cluster effect
sizes) and `assignment.tsv` (one row per isolate, joined to every metadata
column, which is the table you actually want to open). It also writes the run
report in two forms, `report.md` and `report.html`: the same sections and the
same numbers, the second with the diagnostics drawn, so that a panel of
intervals that all cross zero is seen rather than counted. Both read every
value from `cluster_result.json` and recompute nothing, and `report.html` is a
single file with no script and nothing to fetch, so it opens offline.

## Reading the output

`summary.json` is designed to be read top to bottom before anything else:

```json
{
  "selected_k": 2,
  "p_value_structure_report": "p < 1e-100 across all checked split orderings",
  "diagnostic_p_value_structure_exact": 1.3978130041463598e-193,
  "structure_detected": true,
  "discrete_beyond_a_gradient": true,
  "discreteness_status": "ok",
  "discreteness_verdict": "discrete",
  "discreteness_p_value": 0.01,
  "continuum_bootstrap_exceedances": 0,
  "continuum_tail_probability_ci95": {"low": 0.0, "high": 0.03657574498347894},
  "continuum_decision_resolved_at_alpha_0_05": true,
  "continuum_latent_dimension": 3,
  "continuum_null_under_dimensioned": false,
  "n_defining_validated": 17,
  "fusion_n_eff_layers": 1.7893160074349297,
  "fusion_collapse": false,
  "max_empty_stratum_ari": 0.6328112682300324,
  "lineage_concordance_z": 14.0029718182725,
  "lineage_confounded": true,
  "claim_level": 3,
  "claim_status": "statistically_discrete_partition",
  "active_gate_codes": ["empty_stratum", "lineage_confounding", "phenotype_superiority"],
  "mdl_gain_fraction_of_null_code": 0.14891050144447857
}
```

The thresholded `p_value_structure_report` is the publication-facing result;
the exact floating-point value is retained only as a diagnostic. Reporting the
decision over checked split orderings avoids presenting a seed/order-dependent
tail value as a stable measurement.

The claim ladder is deliberately monotone: 0 `no_structure`, 1
`descriptive_partition`, 2 `validated_structure`, 3
`statistically_discrete_partition`, and 4 `archetype_candidate`. The complete
result stores every applicable gate with its trigger, threshold, value and
signed margin. Level 4 means the applicable diagnostics passed; it is not
causal, clinical or cross-cohort biological confirmation.

`structure_detected` alone is **not** a claim that archetypes exist. The
feature-split test rejects whenever the held-out features depend on the
discovery labels, which a single continuous gradient also produces.
`discrete_beyond_a_gradient` is the one that distinguishes them, via a
parametric bootstrap from a latent-trait model whose dimension is chosen by BIC.
Here BIC selects `q = 3`; the adaptive `q = 4` fit turns upward, so the optimum
is interior. Zero of 99 null replicates exceed the observed statistic, and the
exact 95% interval for the underlying tail probability is 0-0.0366: enough to
resolve alpha 0.05, but not a claim that the floor `p = 0.01` is precisely
estimated. If a future
dataset instead selects the largest permitted dimension, the software returns
`discreteness_status = "withheld_under_dimensioned"`, a null verdict and no
bootstrap p-value. It never converts an under-dimensioned null into either a
positive or a negative decision.
`mdl_gain_fraction_of_null_code` puts the description-length gain on an
interpretable scale: 0.63 on the planted control, 0.15 on the *Klebsiella*
cohort, -0.00001 on the null control. (The block above is the *Klebsiella*
`summary.json` verbatim: it detects reproducible, statistically discrete
structure and independently trips artifact, lineage and phenotype gates. That
separation is the case the diagnostics exist for.)

Exit codes: `0` clean, `1` `--validate` failed its threshold, `2` config or
input error (a refused file, an empty cell under the default policy, a value
that is not 0 or 1), `3` the analysis ran but a diagnostic gate tripped.

## How much of this is the clone?

Ask whether a partition of resistance profiles is really a lineage relabelling
and the honest answer in bacteria is "partly, always". The usual instrument is
a permutation test on same-lineage pair co-assignment, computed over the O(n²)
pairs of the cohort, so its null narrows as the cohort grows while the
structure it measures does not. On the shipped *S. suis* example, with the
partition and the lineage variable held fixed and the cohort subsampled, the
statistic reads

| isolates | attributable share `λ` | concordance `z` |
|---|---|---|
| 67 | 0.139 | 1.5 |
| 270 | 0.299 | 10.7 |
| 677 | 0.296 | 25.5 |

so a threshold on `z` is a threshold on cohort size wearing a biological name.
The package keeps the `z` as a descriptor and gates on the magnitude instead.

```python
from amr_clonalshare.attribution import clonal_share, attribute_partition

r = clonal_share(non_wild_type_calls, sequence_types)   # one agent
r.kappa_adj, (r.ci_low, r.ci_high), r.estimable

a = attribute_partition(trait_matrix, cluster_labels, sequence_types)
a.lam, a.shapley_lineage, a.shapley_partition
```

`clonal_share` returns the share of one agent's non-wild-type variance that the
lineage label explains, estimated out of sample so that a typing scheme is not
paid for its levels, debiased on its own permutation null, with a two-stage
bootstrap interval that draws lineages rather than isolates. Read it as the
same fork the mix-versus-rate section below draws, applied to the level rather
than to a difference: **near 1** the resistance travels with the clone and the
lever is biosecurity, movement and mixing; **near 0** it travels independently
of the clone and the lever is selection pressure, dosing and choice of agent.

On the 677 *S. suis* isolates that separates the panel: tetracycline 0.036
(−0.051 to 0.072) and doxycycline 0.029, against penicillin 0.515 (0.247 to
0.706) and ceftiofur 0.504. The determinant calls in the source study, which
take no part in the estimate, agree at Pearson 0.76 and put the *erm* and *tet*
families near zero and the penicillin-binding-protein haplotypes at 0.756.

`attribute_partition` does the same for a whole partition, splitting its
explanatory power against the lineage label by commonality analysis and by the
two-player Shapley value of the variance-explained game. `lam` is the
lineage-attributable share and carries the gate; the default threshold is 0.5
and is set in `attribution.lambda_gate`. A run whose point estimate is below
the threshold but whose interval crosses it is reported as unresolved rather
than cleared.

A singleton lineage cannot be predicted out of sample, so `support`, the share
of isolates in a lineage with at least two members, is reported and
`estimable` is withheld below 0.90, a threshold read off the coverage curve in
`benchmarks/attribution_calibration.py` rather than chosen. The *K. pneumoniae*
sequence types sit at 0.687 and the whole-cohort estimate is withheld.

None of the statistics is new: it is an out-of-sample intraclass correlation
and a standard *R*² decomposition. What is new is the two-player Shapley split
of that decomposition against a lineage label, its application to a routine
susceptibility panel, and a magnitude in the place where a gate that cannot
fail used to stand.

### The same share from the recorded concentration

A non-wild-type call keeps which side of a cut-off an isolate fell and throws
away where in the panel it fell. Both are intervals on the log2 concentration
scale, so one likelihood serves both, and the reading changes how much each
isolate carries rather than what is estimated.

```python
from amr_clonalshare.censored import (censored_clonal_share,
                                         intervals_from_mic, panel_geometry)

panel_geometry(mic_values)                 # what the panel can support
lo, hi = intervals_from_mic(mic_values)    # end wells treated as censored
r = censored_clonal_share(lo, hi, sequence_types)
r.kappa, (r.ci_low, r.ci_high), r.estimable      # the interval to report
r.boot_low, r.boot_high, r.profile_low, r.profile_high   # two checks beside it
```

Three widths are returned and only one is the interval. `ci_low` to `ci_high`
inverts the variance ratio the two mean squares form, which over a grid of
24,000 simulated cohorts covers 0.92 to 0.99 of the time. The cluster bootstrap
is reported beside it and is **not** the interval: it was measured on the same
grid and covers a median 0.61, because resampling thirty lineages does not
reproduce the sampling distribution of a component built from thirty lineages.
The likelihood-ratio width is not an interval either; it measures how much
information a reading carries.

That last one is what improves when a dilution is read rather than a call. On
the shipped *S. suis* panel it falls from 0.250 to 0.200, the reported interval
narrows for ten of the eleven comparable agents, three agents with no
epidemiological cut-off become analysable at all, and two more whose call sits
at 3.8 % and 5.5 % non-wild-type are refused as calls and estimable as
dilutions.

Reading a value on the lowest or highest tested well as censored is an
assumption about the coarsening, not an observation. `panel_geometry` reports
the mass on each end well, a recorded operator column overrides the assumption,
`sensitivity_endpoints` computes the share both ways, and a cohort with more
than half its isolates in lineages wholly beyond the panel is refused.

### Evidence you may re-read every year

Surveillance re-analyses the same panel whenever a year of isolates arrives,
and a false-discovery procedure recomputed at an unplanned number of looks
controls nothing.

```python
from amr_clonalshare.evalues import e_process, e_bh, combine_independent

e = e_process(non_wild_type_calls, sequence_types)
e.e_value, e.reject_05                       # 1/alpha is the threshold
e_bh([x.e_value for x in panel], alpha=0.05) # arbitrary dependence
combine_independent([last_year, this_year])  # a test martingale
```

The e-value here is the betting sense of Vovk and Wang, not the BLAST
expectation value and not the sensitivity measure of the same name. It is the
split likelihood ratio of universal inference on the folds the share already
uses, so its expectation under the null is at most one however often it is
inspected. On the shipped panel e-BH rejects eleven of thirteen agents at 0.05
and leaves tetracycline and doxycycline, the two whose shares are
indistinguishable from zero.

## Surveillance reading: mix versus rate

Every statistic above asks about a *partition*. A surveillance laboratory asks
about an *agent*, and the reported prevalence cannot answer it. When
non-wild-type prevalence differs between two collections, whether two periods,
two countries or two hosts, the difference has two mechanisms with opposite
responses:

* the **lineage mix** changed, because a lineage that already carried the trait
  became a larger share of what was sampled. The response is transmission
  control;
* the **within-lineage rate** changed, because the same lineages became more
  often non-wild-type. The response is selection-pressure control, that is,
  stewardship.

Declare the two collections and every agent is split by the Kitagawa identity,
with a bootstrap interval on each component and Benjamini-Yekutieli control
within each component family, valid whatever the dependence between agents:

```yaml
dataset:
  lineage_column: baps_cluster
  contrast_column: country_period
  contrast_levels: ["United Kingdom late", "United Kingdom early"]
```

On the shipped *S. suis* panel the same thirteen agents give opposite answers
to two contrasts, which is what makes the split worth trusting:

| contrast | mix discoveries | rate discoveries | reading |
|---|---|---|---|
| United Kingdom, 2013-2014 against 2009-2011 | 4 | 0 | the lineages changed |
| Canada against the United Kingdom | 0 | 12 | the same lineages differ in rate |

Ceftiofur non-wild-type prevalence rose 9.5 points between the two United
Kingdom periods, of which 13.1 points is a change in which lineages were
sampled. Erythromycin is 28.5 points higher in Canada than in the United
Kingdom, of which 27.0 points is a higher rate inside the same lineages.

```bash
python examples/ssuis/decompose_trend.py                     # UK periods
python examples/ssuis/decompose_trend.py --contrast country
python examples/ssuis/decompose_trend.py --lineage mlst      # refused, see below
```

### Two gates, and why they exist

**Shared support.** The within-lineage component is identified on lineages
present in both collections. Below `min_shared_support`, 0.8 by default, it is
returned with `within_lineage_estimable` false and a signed margin. The
threshold is measured rather than chosen: across 720,000 simulated
decompositions in `benchmarks/decomposition_calibration.py`, coverage of the
nominal 95 % interval runs 0.82 below a support of 0.4 and reaches 0.94 above
0.8.

**Label availability.** A lineage-resolved statistic describes the labelled
isolates, and describes the collections only if labelling was unrelated to the
trait. Two things must go wrong together: coverage must differ between the
collections, and the trait must differ between labelled and unlabelled
isolates. Either alone is harmless. Both together mean the difference is
computed between two differently selected subsets, and it can carry the wrong
sign.

That is not hypothetical. Repeat the period contrast on `mlst` instead of
`baps_cluster` and ceftiofur appears to *fall* by 11.0 points, because only
40 % of the later period carries a sequence type and the untyped isolates are
the more resistant ones. The collection it is drawn from rises by 9.5 points.
The software refuses the run and names the reason.

Two further readings need no configuration and are always emitted:

* **prevalence per isolate and per lineage.** Two estimands, not a biased and
  an unbiased one: clinical burden against diversity. They separate whenever
  sampling across lineages is uneven. On the *S. suis* cohort ceftiofur
  non-wild-type is 11.1 % per isolate and 32.2 % per lineage, because the
  dominant sequence type is nearly free of it. The interval is a two-stage
  cluster bootstrap over lineages and then isolates within them.
* **effective number of carrying lineages**, the reciprocal of a
  Herfindahl-Hirschman index, with a permutation null on the departure from
  proportional carriage. The departure is a magnitude and fires on dispersion
  as readily as on clonality, so `direction` is reported with it and must be
  read with it, together with the permutation floor and an exact interval for
  the tail.

None of this needs a phylogeny. A sequence type is enough, and
`benchmarks/decomposition_vs_regression.py` measures the estimator against the
lineage-adjusted logistic regression that is the obvious alternative.

## Input format

Each layer is a wide CSV: a strain-ID column plus one 0/1 column per feature.
The config declares the layers, what kind each is, and every choice that changes
the answer:

```yaml
dataset:
  name: my_cohort
  strain_id_column: isolate
  data_dir: data
  metadata: metadata.csv        # optional; not used for clustering
  lineage_column: ST            # enables the population-structure diagnostics
  external_columns: [published_score]

files:
  amr: { path: amr.csv, kind: wide_binary, protected: [Bla_Carb, Col] }
  vir:
    path: vir.csv
    kind: wide_binary
    groups:                     # co-inherited loci, collapsed before distances
      ybt: [ybtS, ybtX, irp1, irp2, fyuA]
      iuc: [iucA, iucB, iucC, iucD, iutA]
  kloc: { path: kloc.csv, kind: one_hot }   # categorical, validated at load

distance:
  metric: jaccard
  undefined_pair: identical     # what two all-zero rows mean. State it.

trait_cluster:
  layers: [amr, vir]
  k_range: [2, 3, 4, 5]         # k = 1 is always added
  collapse_feature_groups: true
```

`kind: one_hot` is not cosmetic: such a layer gets a matching coefficient rather
than Jaccard, is validated to have at most one positive per row, is excluded
from count aggregation, and counts as one hypothesis rather than one per level.

### Empty cells, and what the run does with them

A 0/1 layer has no state for "not measured". An empty field in the CSV, or an
isolate that one layer holds and another does not under
`strain_alignment_policy: union`, is therefore never written as 0, because a 0
reads as absence of the trait. `dataset.missing_policy` decides what happens
instead: `refuse` (default) stops before any estimate and names the rows and
columns; `drop_rows` removes the affected isolates from every layer;
`drop_columns` removes the affected features. Any value other than 0 or 1 is
refused. The choice and its consequences are written to `input_qc.json`.

### Checking the input before the run

```bash
amr-clonalshare --config config.yaml --results-dir out --check-input
```

loads the data, applies the policy above, and writes `input_qc.json` and a
plain-language `input_qc.md`, then stops. The same two files are written at
the start of every full run. They record, per layer, the count of isolates of
the rarer outcome for each trait against the project's reporting threshold of
20, and, when `lineage_column` is set, the size of every lineage group, the
number of single-isolate lineages, the resulting support, and whether the
clonal-share estimator will accept the cohort at that typing resolution. There
is no single minimum sample size: what binds is the number of lineages with at
least two isolates and the count of the rarer outcome, and the check reports
both so the verdict is visible before the estimate is made.

## What is in the box

| module | what it is for |
|---|---|
| `qc` | the input check: empty cells and the policy applied, per-trait counts, per-lineage group sizes and the estimator's verdict |
| `core` | the pipeline: gating, distances, fusion, k selection, consensus, profiles |
| `fusion` | Wang et al. (2014) affinity and cross-diffusion, order-invariant; both published update rules |
| `inference` | feature splitting, p-value merging, count splitting, the continuum null |
| `influence` | leave-one-layer-out influence, effective number of layers |
| `lineage` | same-lineage pair concordance, cluster composition, de-replication |
| `baselines` | concatenation, single layers, Bernoulli mixture, external agreement |
| `archephy` | exact phylogenetic convergence test for small cohorts (needs dendropy) |
| `stats` | binary distances with an explicit empty-union policy, BH, permutation p |
| `small_n` | gap statistic with a matched observed/reference clusterer |

## Methods and calibration

`docs/methodology.md` states every statistic, threshold and assumption with its
citation. `benchmarks/calibration_study.py` regenerates every type-I and power
number the paper quotes:

```bash
python benchmarks/calibration_study.py --quick        # ~2 min smoke run
python benchmarks/calibration_study.py                # the full study
```

A clonal share answers one of two questions, and they do not share an interval:

```python
from amr_clonalshare.realised import realised_share
result = realised_share(y, lineage)      # this collection's own lineages
result.ci_low, result.ci_high            # exact, by inverting the noncentral F
result.superpopulation_low, result.superpopulation_high   # a new draw
```

The same estimator is also run on public data from thirteen species in NCBI
Pathogen Detection, one run as recorded and one with lineage labels permuted:

```bash
bash benchmarks/fetch_pathogen_detection.sh <raw_dir>
python benchmarks/atlas_cross_species.py --raw <raw_dir> --out <out_dir>
```

The manuscript is `paper/manuscript.md`.

## Tests

```bash
pytest -m "not slow"      # a few minutes
pytest                    # + the Monte-Carlo calibration tests
```

## Citing

See `CITATION.cff`. The tagged 1.0.0 release is archived at
https://doi.org/10.5281/zenodo.22306354; the derived tables from which every
number in the article and its supplement was read are archived separately at
https://doi.org/10.5281/zenodo.22307388 (CC BY 4.0). The *Klebsiella* data in `examples/klebsiella/data/` is
**not** covered by this repository's MIT licence. It derives from Lam et al.
(2021) and is redistributed under CC BY 4.0. See
`examples/klebsiella/DATA_PROVENANCE.md`, which also lists the reproducibility
gaps we could not close.

## License

MIT for the code; see [`LICENSE`](LICENSE).
