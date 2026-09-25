# Changelog

## [1.0.0]

First release.

- Share of a binary phenotype carried by lineage: a permutation-adjusted, cross-validated Brier skill score of lineage-specific rates, with lineage-bootstrap limits, a species interval on the scale of a population intraclass correlation, a permuted-label control, the support of repeated lineages and a flag below 80% support.
- Held-out e-values per antimicrobial with e-BH selection, and a running product over intakes that may be read after any of them.
- MIC reading from dilution intervals: a Gaussian random-intercept model on the log2 scale with the panel read within each testing laboratory, categorical covariates (a laboratory, a country, a year) as additive fixed effects, and a calibrated 95% interval for the lineage share of latent variance: Wald's generalised F pivot for exact readings, and for intervals a likelihood-ratio test inverted over candidate values, each calibrated by 199 datasets simulated from the fit restricted to that value. A moment estimate with an approximate F interval is reported beside it as a labelled comparison. `amr-clonalshare-mic` computes the interval for a single table.
- Kitagawa decomposition of a prevalence difference between two collections into composition and within-lineage terms, with bootstrap limits, a shared-support gate and exact recorded-frame bounds for missing outcomes.
- Matched-record comparison of two lineage definitions (`amr-clonalshare-compare`), separating record selection from relabelling.
- Runs repeated within the levels of a metadata column (`stratify_by`).
- Optional population liability intraclass correlation under a grouped probit model (`population_probit`), with a fixed-cut-off profile-likelihood interval and a general bootstrap method.
- Explicit phenotype declarations (clinical S/I/R, WT/NWT, binary), duplicate policies, exact identifier joins with hints for identifiers that differ only in spaces, leading zeros or case, CSV and `.xlsx` input, long or wide call tables, and a first configuration drafted from the column names (`--init`).
- A canonical JSON record, CSV summaries, input diagnostics, Markdown and self-contained HTML reports rendered from the record, and a manifest of file digests; existing output is protected and `--overwrite` replaces only a results directory.
- A result is a function of the data and the seed: it does not depend on the order of rows or antimicrobials, on the number of antimicrobials beside an agent or on the number of worker processes.
- Tests against independent references and forty-digit recomputations, invariances, data mutation, whole-artefact comparisons, refusal wording, cost laws and a recorded mutation-testing baseline (`CONTRIBUTING.md`).
- A validation campaign run from one commit on one pinned stack, with its commands in `benchmarks/campaign/CAMPAIGN.md` and a receipt beside every result.
