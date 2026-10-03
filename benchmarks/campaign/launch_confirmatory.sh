#!/bin/bash
# Submits the confirmatory campaign (CONFIRMATORY_PROTOCOL.md) from a prepared
# campaign directory, which holds source/ (a checkout of the commit that
# carries the protocol), env/ (or AMR_CAMPAIGN_ENV) and rlib/:
#   AMR_CAMPAIGN_ROOT=/path/to/campaign benchmarks/campaign/launch_confirmatory.sh \
#       ACCOUNT PARTITION "module load GCC/13.2.0 R/4.3.3"
# The queues are written; Studies 1 to 3 and Study 4 run at once; the merges,
# the summaries and the verdicts run when both have ended well. Job ids go to
# jobs.txt.
set -euo pipefail
R=${AMR_CAMPAIGN_ROOT:?export AMR_CAMPAIGN_ROOT, the campaign directory}
cd "$R/source"
PY=${AMR_CAMPAIGN_ENV:-$R/env}/bin/python
S="sbatch --parsable -A ${1:?account} -p ${2:?partition} --export=ALL,AMR_CAMPAIGN_ROOT=$R"
C=benchmarks/campaign
Q=$C/commands
L=$R/logs
J=$R/jobs.txt
mkdir -p "$L" "$R/outputs"
: > "$J"
export PYTHONPATH="$R/source/src:$R/source"
"$PY" $C/make_confirmatory_queues.py "${3:?the command that sets up R}"
last() { echo $(( $(cut -f1 "$1" | sort -n | tail -1) )); }
sub() { local name=$1; shift; local id; id=$($S "$@"); echo "$name=$id" >> "$J"; echo "$id"; }
MAIN=$(sub main -o "$L/main_%A_%a.out" --array=0-$(last $Q/confirmatory_main.txt) \
    $C/queue.sbatch $Q/confirmatory_main.txt)
COMP=$(sub competitors -o "$L/competitors_%A_%a.out" --array=0-$(last $Q/confirmatory_competitors.txt) \
    $C/queue.sbatch $Q/confirmatory_competitors.txt)
sub summaries -o "$L/summaries_%j.out" --dependency=afterok:$MAIN:$COMP --nodes=1 --ntasks=1 \
    --cpus-per-task=8 --mem=32G --time=04:00:00 --job-name=amr-summaries \
    $C/run_summaries.sh $Q/confirmatory_summaries.txt >/dev/null
cat "$J"
