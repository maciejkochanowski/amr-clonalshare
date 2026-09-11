#!/usr/bin/env bash
# Run the S. suis decomposition arms the article and supplement read, into one
# evidence tree with a receipt: the Kitagawa split of the period and the
# country contrast at population-cluster and at sequence-type resolution, the
# two published-cut-off arms, and the determinant check. Every command, its
# exit status and the sha256 of every file it wrote are recorded.
#
#   bash scripts/run_decomposition_evidence.sh <out-dir>
set -euo pipefail
OUT=${1:?out dir}
mkdir -p "$OUT/results" "$OUT/logs"
PY=${PYTHON:-python}
RECORDS="$OUT/command_records.json"
echo "{" > "$RECORDS"
first=1
run() {
  name=$1; shift
  start=$(date -u +%Y-%m-%dT%H:%M:%SZ)
  set +e
  "$PY" "$@" > "$OUT/logs/$name.stdout.log" 2> "$OUT/logs/$name.stderr.log"
  rc=$?
  set -e
  end=$(date -u +%Y-%m-%dT%H:%M:%SZ)
  [ $first -eq 1 ] || echo "," >> "$RECORDS"
  first=0
  "$PY" - "$name" "$rc" "$start" "$end" "$PY" "$@" >> "$RECORDS" <<'EOF'
import json, sys
name, rc, start, end = sys.argv[1:5]
argv = sys.argv[5:]
print(json.dumps({name: {"argv": argv, "returncode": int(rc), "started_utc": start,
                         "finished_utc": end}})[1:-1])
EOF
  echo "$name rc=$rc"
}
run metadata_repair examples/ssuis/repair_metadata.py --check
run provenance examples/ssuis/link_source_lineages.py --check
run decompose_period examples/ssuis/decompose_trend.py --contrast period --lineage baps_cluster --n-boot 4000 --seed 23 --json "$OUT/results/decompose_period.json"
run decompose_country examples/ssuis/decompose_trend.py --contrast country --lineage baps_cluster --n-boot 4000 --seed 23 --json "$OUT/results/decompose_country.json"
run decompose_period_mlst examples/ssuis/decompose_trend.py --contrast period --lineage mlst --n-boot 4000 --seed 23 --json "$OUT/results/decompose_period_mlst.json"
run decompose_country_mlst examples/ssuis/decompose_trend.py --contrast country --lineage mlst --n-boot 4000 --seed 23 --json "$OUT/results/decompose_country_mlst.json"
run published_cutoff_period examples/ssuis/decompose_trend.py --contrast period --agents published_cutoff --n-boot 4000 --seed 23 --json "$OUT/results/published_cutoff_period.json"
run published_cutoff_country examples/ssuis/decompose_trend.py --contrast country --agents published_cutoff --n-boot 4000 --seed 23 --json "$OUT/results/published_cutoff_country.json"
run genotype_check examples/ssuis/validate_with_determinants.py --contrast country --n-boot 4000 --seed 7 --json "$OUT/results/genotype_check.json"
echo "}" >> "$RECORDS"
"$PY" - "$OUT" <<'EOF'
import hashlib, json, platform, sys
from pathlib import Path
import numpy, pandas, scipy
import amr_clonalshare
out = Path(sys.argv[1])
records = json.load(open(out / "command_records.json"))
files = sorted(p for p in out.rglob("*") if p.is_file() and p.name not in ("SHA256SUMS", "RUN_RECEIPT.json"))
with open(out / "SHA256SUMS", "w") as fh:
    for p in files:
        fh.write(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}\n")
json.dump({
    "schema": "amr-clonalshare-evidence-receipt-1.0",
    "product_version": amr_clonalshare.__version__,
    "python": platform.python_version(), "numpy": numpy.__version__,
    "pandas": pandas.__version__, "scipy": scipy.__version__,
    "commands": records,
    "status": "PASS" if all(r["returncode"] == 0 for r in records.values()) else "FAIL",
    "n_files": len(files),
}, open(out / "RUN_RECEIPT.json", "w"), indent=1)
print("receipt written,", len(files), "files")
EOF
