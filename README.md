# amr-clonalshare

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![Python](https://img.shields.io/badge/python-3.11%2B-blue.svg)
![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22306353.svg)](https://doi.org/10.5281/zenodo.22306353)

**Release status:** public product version `1.0.0`.

**How much of the antimicrobial resistance in a collection travels with the
clone?** Give it a lineage label per isolate (sequence type, core-genome or
SNP cluster, or serovar) and a susceptibility call or a recorded dilution. It
returns the share of the resistance variance the lineages carry, per
antimicrobial, as three estimands with their own intervals. Two describe the
collection in hand, the design-corrected component ratio `realised_share`
estimates and the lineage-membership share `clonal_share` estimates; beside
them stands a superpopulation share for the species the collection was drawn
from. The exact interval of `realised_share` sits behind a kurtosis gate that
on a binary call opens only for prevalences between about 0.17 and 0.83, so on
a routine panel it is a minority of agents: across the three shipped records it
opens for 12 of 57 agent readings, 5 of 13 on the *S. suis* panel and 2 and 5
of 22 on the two *Salmonella* runs. Each estimate sits behind a gate that refuses where the collection
cannot support it, each carries an e-value, with a sequential e-value for
panels that surveillance re-reads every year, and a prevalence difference
between two collections is split into lineage mix and within-lineage rate. It
reads two tables and never a genome.

## Install

From PyPI, or from an unpacked source tree:

```bash
python -m pip install amr-clonalshare       # released version from PyPI
python -m pip install .                     # from source
python -m pip install -e ".[dev]"           # editable source + test suite
```

The dependencies are NumPy, pandas, SciPy and PyYAML.

## Quickstart

The shipped cohort is 677 *Streptococcus suis* isolates, thirteen agents
called non-wild-type against an epidemiological cut-off, typed by the source
study's population cluster, with the yearly intakes declared so that the
sequential e-process is exercised as well:

```bash
# does the estimator accept the cohort at this typing resolution?
amr-clonalshare --config examples/ssuis/config.yaml --results-dir out/ssuis --check-input

# the run
amr-clonalshare --config examples/ssuis/config.yaml --results-dir out/ssuis
```

A run writes the record, `clonal_share_result.json`, and two reports rendered
from it, `report.md` and `report.html`. The two have the same seven sections
in the same order and the same numbers; the second also has six figures
drawn, so that a panel of intervals that all cross zero is seen at a glance
instead of counted. The figures show the lineage sizes behind the support
gate, every share with both its intervals, the evidence per agent against the
two thresholds, the running product over intakes, the share read from the
dilution beside the share of the call, and the two prevalences with the
concentration of carriage. The first page
is the measurement summary: the quantity measured, the gates, one row per
antimicrobial with its share, interval, permuted-label control, e-value and a
reading fixed by two conditions the record holds, what was not evaluated on
this run and why, and an interpretation drawn from that table by fixed rules;
the pages after it defend the first. Both forms read every value from the
record and recompute nothing, both carry the record's SHA-256, and
`report.html` is a single file with no script and nothing to fetch, so it
opens offline. `input_qc.json` and `input_qc.md`, written first, say whether
the cohort can support the estimate at all. The record the shipped run wrote
is in `examples/ssuis/expected/`.

The second shipped cohort is the poultry-meat *Salmonella* cell of the
article, 7,049 isolates of one NCBI Pathogen Detection release with the
susceptibility calls the release records for 22 agents, typed by serovar and
by SNP cluster, so the same isolates can be read at both resolutions:

```bash
amr-clonalshare --config examples/salmonella_poultry/config.yaml --results-dir out/poultry
amr-clonalshare --config examples/salmonella_poultry/config_cluster.yaml --results-dir out/poultry_cluster
```

`examples/salmonella_poultry/DATA_PROVENANCE.md` says how the cell was cut
from the pinned release and what the two runs show.

Exit codes: `0` the run completed, including a run whose estimate was
withheld; `2` config or input error.

## How much of this is the clone?

In a bacterial collection some of the resistance is always carried by its
clones and some is not. What the package returns is how much, as a share of variance with an interval, against the same estimate
on a permuted labelling, and a verdict on whether the collection can identify
the quantity at all.

```python
from amr_clonalshare.attribution import clonal_share

r = clonal_share(non_wild_type_calls, sequence_types)   # one agent
r.kappa_adj, (r.ci_low, r.ci_high), r.estimable
```

`clonal_share` returns the share of one agent's non-wild-type variance that the
lineage label explains, estimated out of sample so that a typing scheme is not
paid for its levels, debiased on its own permutation null, with two intervals:
a lineage bootstrap for the lineage-membership share of the collection, the
part of its variance that knowing the label carries, and a
superpopulation interval, from a within-lineage bootstrap under a chi-square
layer on the number of lineages, for the share a fresh draw of lineages from
the species would show. For a binary trait the record also restates the share
and both intervals on the latent scale a binomial mixed model reports
(`latent_share`), and it carries the floor the permutation p-value cannot fall
below (`p_floor`). Read the share as the
same fork the mix-versus-rate section below draws, applied to the level rather
than to a difference: **near 1** the resistance travels with the clone and the
lever is biosecurity, movement and mixing; **near 0** it travels independently
of the clone and the lever is selection pressure, dosing and choice of agent.

On the 677 *S. suis* isolates the estimated share separates the panel: tetracycline 0.028
(−0.049 to 0.074) and doxycycline 0.022, against penicillin 0.517 (0.237 to
0.698) and ceftiofur 0.507. The determinant calls in the source study, which
take no part in the estimate, put the *erm* and *tet* families near zero and
the penicillin-binding-protein haplotypes high, and the six agents whose
determinants sit on the chromosome, entirely or for most carriers, sit above
the seven mobile-borne ones in 41 of the 42 pairs at this seed; the article's
run of the same cohort, `benchmarks/ssuis_mechanism.py` at its own seed, puts
amoxicillin a hair below lincomycin as well and reads 40, which is the width
of the Monte Carlo noise on a boundary pair.

A singleton lineage cannot be predicted out of sample, so `support`, the share
of isolates in a lineage with at least two members, is reported and
`estimable` is withheld below 0.90, a threshold read off the coverage curve in
`benchmarks/attribution_calibration.py`. Read the same
cohort at sequence-type resolution, where a third of the isolates carry no
type, and the gate says so (`examples/ssuis/config_contrast.yaml`).

None of the statistics is new: it is an out-of-sample intraclass correlation
and a standard *R*² decomposition. What is new is its application to a
routine susceptibility panel at a stated typing resolution, with the
conditions a bacterial collection violates made explicit, gated and measured.

### The same share from the recorded concentration

A non-wild-type call keeps which side of a cut-off an isolate fell and throws
away where in the panel it fell. Both are intervals on the log2 concentration
scale, so one likelihood serves both, and the reading changes only how much
each isolate carries; what is estimated stays the same.

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
inverts the variance ratio the two mean squares form, which for a recorded
value or a dilution interval covers 0.93 to 0.98 of the time by true share and
0.85 to 1.00 cell by cell over 40 cells of 200 simulated cohorts, the lowest
at a hundred lineages of six isolates; read from a single cut point instead,
the same interval keeps its level only where the cut sits near the middle of
the distribution, which is why a call is read by `clonal_share` and the
dilution by this function. The cluster bootstrap
is reported beside it and is **not** the interval: on the first run of the
same grid, with 200 draws a cohort, it covered a median 0.61, because
resampling thirty lineages does not reproduce the sampling distribution of a
component built from thirty lineages.
The likelihood-ratio width is not an interval either; it measures how much
information a reading carries.

That last one is what improves when a dilution is read instead of a call. On
the shipped *S. suis* panel it falls from 0.249 to 0.175 and the reported
interval narrows for ten of the thirteen agents read both ways
(`benchmarks/censored_real_cohort.py`, whose result ships in
`benchmarks/results_censored/`); three agents with no epidemiological cut-off
become analysable at all, and two more whose call sits at 3.8 % and 5.5 %
non-wild-type are refused as calls and estimable as dilutions.

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
from amr_clonalshare.evalues import e_process, e_bh, sequential_e_process

e = e_process(non_wild_type_calls, sequence_types)   # one look
e.e_value, e.reject_05                               # 1/alpha is the threshold
e_bh([x.e_value for x in panel], alpha=0.05)         # arbitrary dependence
seq = sequential_e_process([calls_2019, calls_2020, calls_2021],
                           [types_2019, types_2020, types_2021])
seq.e_value[-1], seq.reject_05                       # a test supermartingale
```

The e-value here is the betting sense of Vovk and Wang, not the BLAST
expectation value and not the sensitivity measure of the same name. It is the
split likelihood ratio of universal inference on the folds the share already
uses, so its expectation under the null is at most one. One call is one look;
for a programme of intakes `sequential_e_process` fits on every earlier intake
and scores the new one, and its running product may be inspected after any
intake and stopped at any time. The command line reaches the same construction
through `dataset.batch_column`, the metadata column naming the intake an
isolate arrived in:

```yaml
dataset:
  lineage_column: baps_cluster
  batch_column: collection_year      # must sort into arrival order
```

The record then carries `lineage_evidence.sequential`: the batches, the running
product per trait, and an e-BH decision taken on the log scale, where a product
over many intakes does not overflow. On the shipped panel e-BH rejects eleven of
thirteen agents at 0.05 at a single look and leaves tetracycline and
doxycycline, the two whose shares are indistinguishable from zero; the
sequential product over the twenty-five yearly intakes, many of them a
handful of isolates, selects penicillin and ceftiofur only. The sequential
guarantee holds at whatever intake the programme is read at, and the smaller
selection is what that guarantee costs.

## Surveillance reading: mix versus rate

Every statistic above asks about a level. A surveillance laboratory also asks
about a change, and the reported prevalence cannot answer it. When
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
to two contrasts:

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
amr-clonalshare --config examples/ssuis/config_contrast.yaml --results-dir out/ssuis-contrast
```

### Two gates, and why they exist

**Shared support.** The within-lineage component is identified on lineages
present in both collections. Below `min_shared_support`, 0.8 by default, it is
returned with `within_lineage_estimable` false and a signed margin. The
threshold is read off the coverage curve in
`benchmarks/decomposition_calibration.py`: across 720,000 simulated
decompositions, coverage of the
nominal 95 % interval runs 0.82 below a support of 0.4 and reaches 0.94 above
0.8.

**Label availability.** A lineage-resolved statistic describes the labelled
isolates, and describes the collections only if labelling was unrelated to the
trait. Two things must go wrong together: coverage must differ between the
collections, and the trait must differ between labelled and unlabelled
isolates. Either alone is harmless. Both together mean the difference is
computed between two differently selected subsets, and it can carry the wrong
sign.

Repeat the period contrast on `mlst` instead of
`baps_cluster` and ceftiofur appears to *fall* by 11.0 points, because only
40 % of the later period carries a sequence type and the untyped isolates are
the more resistant ones. The collection it is drawn from rises by 9.5 points.
The software refuses the run and names the reason.

Two further readings need no configuration and are always emitted:

* **prevalence per isolate and per lineage.** Two estimands, not a biased and
  an unbiased one: clinical burden against diversity. They separate whenever
  sampling across lineages is uneven. On the *S. suis* cohort ceftiofur
  non-wild-type is 24.1 % per isolate and 41.9 % per lineage, because the
  largest population cluster is nearly free of it. The interval is a two-stage
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

Two tables, plain CSV, declared in one YAML configuration. The phenotype
table is long, one row per isolate and antimicrobial, which is the shape a
laboratory exports and NCBI Pathogen Detection and BV-BRC serve; the call
column holds `susceptible`, `non-susceptible`, `resistant` or `intermediate`
in any case, and `phenotype_intermediate` states what an intermediate counts
as. The metadata table holds one row per isolate with the lineage column and,
where declared, the intake and the contrast column. A recorded dilution
table, long as well, is optional.

```yaml
dataset:
  name: my_cohort
  strain_id_column: isolate
  data_dir: data
  metadata: metadata.csv
  lineage_column: ST             # sequence type, cluster or serovar
  batch_column: year             # optional: the intake, in arrival order
  phenotype: calls.csv
  phenotype_id_column: isolate
  phenotype_antibiotic_column: antibiotic
  phenotype_call_column: call
  mic: mic.csv                   # optional: recorded dilutions
  mic_id_column: isolate
```

Every key the loader does not read is refused with the nearest known key
named: a misspelt `lineage_colum` would otherwise silently switch the
estimate off. `examples/ssuis/config.yaml` is a fully annotated real
configuration.

### Checking the input before the run

```bash
amr-clonalshare --config config.yaml --results-dir out --check-input
```

loads the tables, joins them, and writes `input_qc.json` and a plain-language
`input_qc.md`, then stops. The same two files are written at the start of
every full run. They record, per antimicrobial, the count of isolates of the
rarer outcome against the project's reporting threshold of 20, and, for the
lineage column, the size of every lineage group, the number of single-isolate
lineages, the resulting support, and whether the estimator will accept the
cohort at that typing resolution. There is no single minimum sample size:
what binds is the number of lineages with at least two isolates and the count
of the rarer outcome, and the check reports both so the verdict is visible
before the estimate is made.

## What is in the box

| module | what it is for |
|---|---|
| `attribution` | the out-of-sample lineage-membership share against its permuted null, with its two intervals |
| `realised` | the design-corrected component ratio with its exact interval, and the superpopulation interval beside it |
| `latent` | the share restated on the liability scale of a threshold model |
| `censored` | the same share read from a dilution panel as an interval likelihood |
| `evalues` | the e-value per agent, the sequential e-process over intakes, e-BH |
| `clonality` | the Kitagawa decomposition, prevalence on two scales, concentration of carriage |
| `phenotype` | the reader of a long susceptibility table |
| `stats` | the permutation p-value, the Benjamini-Hochberg and Benjamini-Yekutieli step-ups, Fisher's exact test and the effective number of independent agents |
| `qc` | the input check: per-antimicrobial counts, per-lineage group sizes and the estimator's verdict |
| `core`, `cli` | the run, the record and the command line |
| `report`, `report_model`, `report_html` | the two reports, rendered from one structure |

## The examples the article reads

Three public collections, none of which shares the others' choices. Every
number below is read from a run of the shipped scripts, and the runs
themselves are in the evidence deposit of the article
(https://doi.org/10.5281/zenodo.22307388), which names each file by path and
digest.

**A surveillance signal read from a serovar label.** The *Salmonella*
isolates of one NCBI Pathogen Detection release, cut by host and sampling
point, with the lineage defined as the SNP cluster and again as the serovar
(`benchmarks/vet_atlas.py`, `benchmarks/resolution_atlas.py`). The package
reports the difference between the two resolutions as refusals: at
SNP-cluster resolution the median support is 0.735, short of the 0.90 the
gate requires, and 28 of 131 cells are estimable against 117 at serovar
resolution; where both are, the serovar figure is lower in every cell, by a
median of 0.178, and the permuted control excludes zero in 4 of 131 cells
against 6.6 expected. Among 6,915 serotyped poultry-meat isolates with a
nalidixic-acid result the share is 0.935 at SNP-cluster and 0.766 at serovar
resolution (the shipped example, run at its own seed, gives 0.937 and 0.765);
split at 2016, a cut written into the analysis script before the run
(`benchmarks/period_split.py`), the serovar share is 0.012 before the cut and
0.846 after, while prevalence rose 29-fold, which places the rise in one
serovar of the collection rather than across the serovars it already held.
That cell ships in `examples/salmonella_poultry/`, cut from the pinned
release by its own script, with a configuration for each resolution.

**A swine pathogen typed by population cluster.** The 677 *Streptococcus
suis* isolates shipped in `examples/ssuis/`, read above: the six agents
whose determinants the source study places on the chromosome, entirely or
for most carriers, average 0.356 and the seven whose determinants sit on
mobile elements 0.128, the chromosomal agent being higher in 40 of the 42
pairs, and the determinant calls of the source study, which the estimate
never sees, order the agents the same way (`benchmarks/ssuis_mechanism.py`);
the same cohort read at serotype resolution attenuates the chromosomal
shares (`benchmarks/ssuis_resolution.py`), and the Kitagawa split reverses
between the period contrast and the country contrast.

**The same estimator across species.** Fifty organisms of one Pathogen
Detection release, each at the release accession the retrieval recorded
(`benchmarks/fetch_pathogen_detection.sh`, `benchmarks/atlas_cross_species.py`;
the same fetch feeds the *Salmonella* scripts above): every agent of every cohort that reaches the estimator leaves one row,
estimated or refused with the condition that failed; 192 cells are estimated
in twelve species and the support gate admits 44 in six, each with a
permuted-label control.

The article and its supplement are not part of this repository. The
supplement, with every table the article reads from, is in the data record
cited below.

## Methods and calibration

`docs/methodology.md` states every statistic, threshold and assumption with its
citation. Every estimator was scored against the truth it estimates on a grid
of 135 simulated cells, five lineages to a thousand:

```bash
python benchmarks/attribution_calibration.py --quick   # the share and its support gate
python benchmarks/ssuis_mechanism.py                   # the S. suis share against its determinants
python benchmarks/estimator_benchmark.py --plan        # the estimator grid, cell by cell
python benchmarks/decomposition_calibration.py         # the shared-support gate
python benchmarks/censored_calibration.py              # the dilution reading
```

A clonal share answers one of two questions, and they do not share an interval:

```python
from amr_clonalshare.realised import realised_share
result = realised_share(y, lineage)      # this collection's own lineages
result.ci_low, result.ci_high            # exact under the Gaussian one-way model
result.superpopulation_low, result.superpopulation_high   # a new draw
```

## Tests

```bash
pytest -m "not slow"      # a few minutes
pytest                    # + the Monte-Carlo calibration tests
```

## Citing

See `CITATION.cff`. The tagged 1.0.0 release is archived at
https://doi.org/10.5281/zenodo.22306353, the Zenodo concept DOI, which always
resolves to the version record deposited from this tag; the derived tables from which every
number in the article and its supplement was read are archived separately at
https://doi.org/10.5281/zenodo.22307388 (CC BY 4.0). The *S. suis* panel in
`examples/ssuis/data/` derives from Hadjirin et al. (2021) and is redistributed
under CC BY 4.0; see `examples/ssuis/DATA_PROVENANCE.md`. The *Salmonella*
cell in `examples/salmonella_poultry/data/` is cut from a public NCBI
Pathogen Detection release; see `examples/salmonella_poultry/DATA_PROVENANCE.md`.

## License

MIT for the code; see [`LICENSE`](LICENSE).
