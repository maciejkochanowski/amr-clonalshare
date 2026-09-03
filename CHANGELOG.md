# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/);
versioning: [SemVer](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-09-02

First public release. Everything below is what the package contains, not a
record of how it was built: the development history is in the repository log
and the pre-submission audit is in `audit_2026-08-31/FINDINGS.md`.

### Added

**Input check.** Every run starts by loading the layers, applying the stated
policy to cells that hold no value (refuse, drop the isolate, or drop the
feature; never fill in a 0), refusing any value other than 0 and 1, and writing
`input_qc.json` with a plain-language `input_qc.md`: per trait the count of the
rarer outcome, and per lineage the group sizes, the share of isolates in
lineages of at least two, and whether the clonal-share estimator will accept
the cohort. `--check-input` runs this step alone. A refused input leaves with
exit code 2 and a sentence naming the position, not a traceback.

**Run report.** `report.md` beside `cluster_result.json`, in a fixed section
order, reading every number from the record and stating in plain language what
was computed, what the value means, and which gate limits the reading.

**Partition diagnostics.** Multi-view consensus clustering over binary trait
layers with similarity network fusion, order-invariant neighbour selection and
a description-length criterion whose sweep includes `k = 1`, so a cohort with
no structure can be reported as having none. An ordered sequence of artefact,
layer, held-out-feature and adaptively dimensioned continuum diagnostics runs
before any partition is released, and each is tied to a status code rather than
to prose.

**Lineage attribution.** A per-antimicrobial clonal share and a per-partition
attributable share, both estimated out of sample so that a lineage variable is
not paid for its levels, with a cross-validation debias, a cluster bootstrap
interval that draws lineages whole, a commonality and Shapley split of the
joint term, and a support gate whose threshold is read off a coverage curve.
The gate on a partition is a magnitude, not a significance test.

**Interval-censored panel reading.** The same share from a dichotomised
non-wild-type call, a recorded minimum inhibitory concentration or a censored
reading, which enter one interval likelihood at three widths. Panel geometry is
reported before any share, an end-well reading is treated as censored under a
stated coarsening assumption with the alternative computed beside it, a
recorded censoring operator overrides that assumption, and a cohort with more
than half its isolates in lineages wholly beyond the panel is refused. Two
intervals are given: a cluster bootstrap and a profile likelihood from the
marginal likelihood, which answer different questions.

**Anytime-valid evidence.** An e-value per antimicrobial, in the betting sense,
built as the split likelihood ratio of universal inference on the same folds
the share already uses, with e-BH across the panel controlling the false
discovery rate under arbitrary dependence. A surveillance programme may inspect
the running value at every intake without spending an error budget.

**Prevalence decomposition.** A Kitagawa split of a prevalence difference into
a lineage-composition component and a within-lineage rate component, with an
exact residual, Benjamini-Hochberg control within each component family, the
effective number of independent agents, and two estimability gates: shared
lineage support, and the pair of conditions under which missing lineage labels
make the two collections incomparable.

**Surveillance readings that need no contrast.** Prevalence per isolate and per
lineage, the effective number of lineages carrying a feature, and the direction
of a departure from proportional carriage.

**Contract and reproducibility.** A strict YAML configuration that refuses an
unknown key rather than silently disabling the analysis it names; four output
artifacts under public output schema 1.0; exit codes 0, 1, 2 and 3 separating
success, a failed planted-truth validation, a configuration error and an
informative refusal; a master seed spawning every stochastic stage; and a
release pipeline that records input hashes, an environment manifest and a
claim register.

**Examples.** Two real cohorts, 677 *Streptococcus suis* isolates with 16
antimicrobials of recorded concentrations and a hierBAPS population cluster for
every isolate, and a *Klebsiella pneumoniae* subset; three planted synthetic
controls; and a holdout arm.

### Known limits

- The clonal share and the censored share are different estimands, one
  predictive skill and one population variance, and are reported as such.
- A single cut point does not identify the latent scale. The censored share is
  then reported on the liability scale by convention, and the calibration in
  `benchmarks/censored_calibration.py` gives its bias as a function of the
  prevalence at the cut.
- The examples are not a cross-species validation study, and no superiority
  over all multi-view methods is claimed.
