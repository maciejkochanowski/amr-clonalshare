#!/bin/bash
# Submits the runs of the re-issue of 1.0.0 from a prepared campaign directory,
# which holds source/ (a checkout of the final commit) and env/ (or
# AMR_CAMPAIGN_ENV, a virtual environment built from requirements-lock.txt):
#   AMR_CAMPAIGN_ROOT=/path/to/campaign benchmarks/campaign/launch_reissue.sh \
#       ACCOUNT PARTITION [--goldens | --mutation]
# With --goldens it submits one short job only: the golden artefacts
# regenerated on the cluster's stack (outputs/goldens/goldens.patch). They are
# committed first, so that the runs proper start from a commit whose own test
# suite passes; mutmut refuses to start on a failing suite.
# Four jobs, each one node of the partition: the empirical queue (the eight
# configured analyses and the comparison, the four empirical scripts, the seed-stability check
# and the profile), the auxiliary queue (the decomposition calibration and its
# gate, the regression comparison, the null uniformity and the 24 cells of the
# repeated looks), the checks (the test suite and the high-precision oracles
# on the cluster's stack) and, with --mutation, the mutation run (about half an
# hour on 192 cores; needs mutmut in the environment). Job ids are appended to
# jobs.txt under the commit they run from. Nothing is copied back: when the
# jobs have ended, benchmarks/campaign/collect_reissue.py places the outputs
# in the results folders and checks every receipt against the source.
set -euo pipefail
R=${AMR_CAMPAIGN_ROOT:?export AMR_CAMPAIGN_ROOT, the campaign directory}
cd "$R/source"
PY=${AMR_CAMPAIGN_ENV:-$R/env}/bin/python
ACCOUNT=${1:?account}
PARTITION=${2:?partition}
MUTATION=${3:-}
S="sbatch --parsable -A $ACCOUNT -p $PARTITION --export=ALL,AMR_CAMPAIGN_ROOT=$R"
C=benchmarks/campaign
L=$R/logs
J=$R/jobs.txt
mkdir -p "$L" "$R/outputs"
COMMIT=$(git rev-parse --short HEAD)
[ -z "$(git status --porcelain)" ] || { echo "source/ has uncommitted changes; the receipts would name a commit the tree is not" >&2; exit 2; }
"$PY" -c "import numpy, scipy, pandas; print('stack', numpy.__version__, scipy.__version__, pandas.__version__)"
if [ "$MUTATION" = "--mutation" ]; then
    "$PY" -c "import mutmut" 2>/dev/null || { echo "mutmut is not installed in the environment: $PY -m pip install mutmut" >&2; exit 2; }
fi
echo "## reissue" >> "$J"
sub() { local name=$1; shift; local id; id=$($S "$@"); echo "commit=$COMMIT $name=$id" >> "$J"; echo "$name=$id"; LAST=$id; }
if [ "$MUTATION" = "--goldens" ]; then
    sub goldens -o "$L/goldens_%A_%a.out" --array=0-0 --cpus-per-task=4 --mem=16G --time=00:30:00 $C/queue.sbatch $C/goldens_queue.txt
    tail -2 "$J"
    exit 0
fi
sub empirical -o "$L/empirical_%A_%a.out" --array=0-0 --time=04:00:00 $C/queue.sbatch $C/empirical_queue.txt
sub auxiliary -o "$L/auxiliary_%A_%a.out" --array=0-0 --time=04:00:00 $C/queue.sbatch $C/auxiliary_queue.txt
sub checks -o "$L/checks_%A_%a.out" --array=0-0 --cpus-per-task=8 --mem=32G --time=02:00:00 $C/queue.sbatch $C/checks_queue.txt
if [ "$MUTATION" = "--mutation" ]; then
    # after the checks: they rewrite tests/golden for a moment in the same checkout
    sub mutation -o "$L/mutation_%A_%a.out" --array=0-0 --time=04:00:00 --dependency=afterany:$LAST $C/queue.sbatch $C/mutation_queue.txt
fi
tail -5 "$J"
