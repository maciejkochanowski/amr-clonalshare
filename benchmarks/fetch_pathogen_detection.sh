#!/bin/bash
# Retrieve the NCBI Pathogen Detection tables the cross-species atlas reads.
#
#     bash benchmarks/fetch_pathogen_detection.sh <raw_dir> [organism ...]
#
# For each organism the script takes the current latest_snps release, keeps
# the metadata rows that carry an AST_phenotypes value, and keeps the cluster
# assignment file whole. The release accession is written to the log, because
# the release is what a rerun would have to match: NCBI recomputes the SNP
# clusters at every release, so a later release is a different lineage
# variable, not a longer version of the same one.
#
# The column names are read from the live header rather than assumed by
# position. Downloads run eight at a time; NCBI throttles per connection, so
# eight streams finish in about the time one takes.
set -u

RAW=${1:?usage: fetch_pathogen_detection.sh <raw_dir> [organism ...]}
shift || true
NCBI=https://ftp.ncbi.nlm.nih.gov/pathogen/Results
LOGS="$RAW/../logs"
mkdir -p "$RAW" "$LOGS/org"

DEFAULT_ORGS="Salmonella Campylobacter Escherichia_coli_Shigella Klebsiella
Acinetobacter Enterococcus_faecalis Enterococcus_faecium
Staphylococcus_aureus Staphylococcus_pseudintermedius Listeria
Pseudomonas_aeruginosa Neisseria_gonorrhoeae Streptococcus_pneumoniae
Clostridioides_difficile Enterobacter_hormaechei Citrobacter_freundii
Vibrio_cholerae Vibrio_parahaemolyticus"
ORGS=${*:-$DEFAULT_ORGS}

one() {
  local ORG="$1"
  if [ -s "$RAW/${ORG}.ast.tsv" ] && [ -s "$RAW/${ORG}.clusters.tsv" ]; then
    echo "$ORG already present"; return 0
  fi
  local MU="$NCBI/$ORG/latest_snps/Metadata/"
  local MF
  MF=$(timeout 120 curl -sSL "$MU" | grep -oE 'PDG[0-9.]+\.metadata\.tsv' | head -1)
  [ -z "$MF" ] && { echo "$ORG: no metadata release, skipping"; return 1; }
  local CU="$NCBI/$ORG/latest_snps/Clusters/"
  local CF
  CF=$(timeout 120 curl -sSL "$CU" | grep -oE 'PDG[0-9.]+\.reference_target\.all_isolates\.tsv' | head -1)
  [ -z "$CF" ] && { echo "$ORG: no cluster release, skipping"; return 1; }
  echo "$ORG start $(date -Is) metadata=$MF cluster=$CF"

  timeout 10800 curl -sSL "$MU$MF" \
    | awk -F'\t' -v OFS='\t' '
        NR==1 { for(i=1;i<=NF;i++) h[$i]=i
                print "target_acc","AST_phenotypes","scientific_name",
                      "collection_date","geo_loc_name","host",
                      "isolation_source","epi_type"
                next }
        h["AST_phenotypes"] && $h["AST_phenotypes"] != "" && $h["AST_phenotypes"] != "NULL" {
                print $h["target_acc"], $h["AST_phenotypes"], $h["scientific_name"],
                      $h["collection_date"], $h["geo_loc_name"], $h["host"],
                      $h["isolation_source"], $h["epi_type"] }
      ' > "$RAW/${ORG}.ast.tsv.part" && mv "$RAW/${ORG}.ast.tsv.part" "$RAW/${ORG}.ast.tsv"

  timeout 7200 curl -sSL "$CU$CF" > "$RAW/${ORG}.clusters.tsv.part" \
    && mv "$RAW/${ORG}.clusters.tsv.part" "$RAW/${ORG}.clusters.tsv"

  echo "$ORG done $(date -Is) AST rows=$(($(wc -l < "$RAW/${ORG}.ast.tsv")-1))" \
       "cluster rows=$(($(wc -l < "$RAW/${ORG}.clusters.tsv")-1))"
}

for ORG in $ORGS; do
  one "$ORG" > "$LOGS/org/$ORG.log" 2>&1 &
  while [ "$(jobs -rp | wc -l)" -ge 8 ]; do sleep 10; done
done
wait
echo "all organisms done $(date -Is)"
cat "$LOGS"/org/*.log
