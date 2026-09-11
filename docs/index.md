# amr-clonalshare

A command-line tool that answers one question about a collection of bacterial
isolates: **how much of the resistance you see is carried by the lineages the
collection holds, and how much moves between them?**

The answer is a share between 0 and 1 with an interval, per antimicrobial, and
a decision about whether your collection can support an estimate at all. The
exact interval for the lineages the collection holds sits behind a kurtosis
gate that on a binary call opens only for prevalences between about 0.17 and
0.83, so on a routine panel it is a minority of agents: across the three
shipped records it opens for 12 of 57 agent readings, 5 of 13 on the *S. suis*
panel and 2 and 5 of 22 on the two *Salmonella* runs. Around
it sit a decomposition of a prevalence difference between two collections
into a change in lineage mix and a change in within-lineage rate, an
interval-censored reading of recorded MICs, and an e-value per agent, with a
sequential e-value for panels that surveillance re-reads every year. It takes
a susceptibility table and a lineage column as input; it does not read
sequence data.

## Who it is for

Veterinary and public-health laboratories that hold a typed collection with a
susceptibility panel, and want a number, not a tree, to say whether resistance
is clone-locked or mobile. No statistics beyond reading an interval is needed
to use it; the run writes a report in plain language.

## Where to start

1. [Install](manual/01-install.md) with `pip` or run the container.
2. [Prepare the input](manual/02-input.md): a susceptibility table, one row
   per isolate and antimicrobial, and a metadata table with the lineage.
3. [Check the input](manual/03-check-input.md) before spending compute.
4. [Run](manual/04-run.md) and [read the results](manual/05-results.md).

## What the shipped cohort shows

The package ships 677 *Streptococcus suis* isolates typed by population
cluster, with thirteen antimicrobials called against an epidemiological
cut-off and the recorded dilutions beside them. The run reports which agents
travel with the clone and which do not, the same reading from the dilutions,
the e-value per agent and its running product over the yearly intakes, and,
at sequence-type resolution, what the support gate refuses and why. This run
is the first thing to reproduce on a new machine; the
[container page](manual/08-container.md) does exactly that.
