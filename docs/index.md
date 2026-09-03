# amr-clonalshare

A command-line tool that answers one question about a collection of bacterial
isolates: **how much of the resistance you see is carried by the lineages the
collection holds, and how much moves between them?**

The answer is a share between 0 and 1 with an interval, per antimicrobial, and
a decision about whether your collection can support an estimate at all. Around
it sit the diagnostics that decide whether a resistance-profile clustering means
anything, a decomposition of a prevalence difference between two collections
into a change in lineage mix and a change in within-lineage rate, an
interval-censored reading of recorded MICs, and an e-value per agent for panels
that surveillance re-reads every year.

## Who it is for

Veterinary and public-health laboratories that hold a typed collection with a
susceptibility panel, and want a number, not a tree, to say whether resistance
is clone-locked or mobile. No statistics beyond reading an interval is needed
to use it; the run writes a report in plain language.

## Where to start

1. [Install](manual/01-install.md) with `pip` or run the container.
2. [Prepare the input](manual/02-input.md): one CSV per layer, one row per
   isolate, one column per trait, 0 or 1 in every cell, plus a metadata table
   with the lineage.
3. [Check the input](manual/03-check-input.md) before spending compute.
4. [Run](manual/04-run.md) and [read the results](manual/05-results.md).

## What the shipped controls show

The package ships a positive control (three planted archetypes), an independent
confirmation with a second seed, a negative control with no structure, and an
adversarial control where every difference is clonal. The first recovers the
planted truth, the second recovers it again, the third returns one group, the
fourth is flagged as lineage-confounded. These four runs are the first thing to
reproduce on a new machine; the [container page](manual/08-container.md) does
exactly that.
