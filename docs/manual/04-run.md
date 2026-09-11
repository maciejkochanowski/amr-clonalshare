# 4. Run

```bash
amr-clonalshare --config config.yaml --results-dir out --seed 42 --threads 4
```

| option | meaning |
|---|---|
| `--config` | the YAML configuration (required) |
| `--results-dir` | where every output file goes; created if absent |
| `--seed` | the master seed; every stochastic stage is spawned from it, so the same seed on the same stack gives the same record |
| `--threads` | BLAS and OpenMP thread limit |
| `--check-input` | stop after the input check |
| `--quiet` | suppress the digest on standard output |
| `--no-check-files` | validate the configuration without requiring the named tables to exist |

![A full run on the shipped cohort](../img/04_run.png)

## What happens, in order

1. **Input check.** As in the previous page; a refusal stops here.
2. **The share, per antimicrobial.** The lineage-membership share estimated
   out of sample against its permuted control, with the lineage bootstrap
   interval and the species interval beside it, restated on the latent scale
   for a binary call, and the support gate that withholds it where too few
   isolates sit in lineages of two or more.
3. **Evidence.** An e-value per antimicrobial for this look and, when
   `batch_column` is declared, the running product over the intakes, with
   e-BH control across the panel.
4. **The recorded dilutions.** Where `mic` names a table, the same share read
   from the dilution as an interval likelihood, with the panel geometry and
   the end-well sensitivity.
5. **Carriage.** Prevalence per isolate and per lineage, the effective number
   of carrying lineages, and, when a contrast column is configured, the
   mix-versus-rate decomposition.
6. **Outputs.** `clonal_share_result.json` (the full record), `report.md`
   and `report.html`, `input_qc.json` and `input_qc.md`.

## Run time

The shipped cohort, 677 isolates and thirteen antimicrobials with the
dilution reading and the sequential e-process, runs in about a minute and a
half on one thread of a 2.1 GHz Xeon core; the poultry cell, 7,049 isolates
and 22 antimicrobials at serovar resolution, in under a minute. The cost is dominated by the bootstrap and permutation counts in
the `attribution`, `surveillance` and `evidence` sections of the config, which
scale linearly; halve them for a first look and restore them for the run you
report.

## Reproducing a run

The record carries the seed, the configuration, and the library versions.
Two runs on the same stack with the same seed give the same record. Two runs
on different BLAS builds agree on every reported decision and every share to
the precision reported, and can differ in the last digits of a bootstrap
interval; the container removes that difference.
