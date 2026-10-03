#!/bin/bash
# Usage: run_summaries.sh COMMANDS_FILE  - runs every line in order, stops on the first failure.
source "${AMR_CAMPAIGN_ROOT:?}/source/benchmarks/campaign/env.sh"
while IFS= read -r LINE; do
  [ -n "$LINE" ] || continue
  echo "== $LINE"
  eval "$LINE" || { echo "FAILED: $LINE"; exit 1; }
done < "${1:?commands file}"
