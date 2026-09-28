#!/usr/bin/env bash
#===============================================================================
# verify_farci_houst.sh
#
# Closes the taxonomy question for the clonal cluster with two decisive checks:
#   (1) Are M. houstonense and M. farcinogenes the same genomospecies?
#       -> all-vs-all fastANI between their type strains (+ the cluster).
#   (2) Why did OUR farcinogenes reference (GCF_000723385.1) behave anomalously?
#       -> ANI of GCF_000723385.1 vs the other DSM 43637 assembly (GCA_025821245)
#          and vs the whole fortuitum-complex type panel (what does it match?),
#          plus optional CheckM2 to test for contamination.
#
# Type-strain accessions are taken from the TYGS Table 4 (authoritative DSMZ set).
#
# Usage:
#   ./verify_farci_houst.sh -q CLUSTER_REP.fna -o OUT_DIR [-c] [-t N]
#
# Required:
#   -q FILE  one genome from the clonal cluster (any of the 7; they are clonal)
#   -o DIR   output directory
# Options:
#   -t INT   threads (default 8)
#   -c       also run CheckM2 (completeness/contamination) on every genome  (off)
#   -h       help
#
# Dependencies (bioconda): ncbi-datasets-cli fastani seqkit [checkm2 if -c]
#===============================================================================
set -euo pipefail
THREADS=8; RUN_CHECKM2=0; QREP=""; OUT=""
usage(){ sed -n '2,33p' "$0"; exit "${1:-0}"; }
log(){ echo -e "[$(date '+%F %T')] $*"; }
while getopts ":q:o:t:ch" opt; do case $opt in
  q) QREP=$OPTARG;; o) OUT=$OPTARG;; t) THREADS=$OPTARG;; c) RUN_CHECKM2=1;;
  h) usage 0;; \?) echo "bad opt -$OPTARG">&2; usage 1;; :) echo "-$OPTARG needs arg">&2; usage 1;;
esac; done
[[ -z "$QREP" || -z "$OUT" ]] && { echo "ERROR: -q and -o required">&2; usage 1; }
[[ -s "$QREP" ]] || { echo "ERROR: cluster rep not found: $QREP">&2; exit 1; }
for t in datasets fastANI seqkit; do command -v $t >/dev/null || { echo "MISSING: $t">&2; exit 1; }; done
[[ $RUN_CHECKM2 -eq 1 ]] && { command -v checkm2 >/dev/null || { echo "MISSING: checkm2 (needed for -c)">&2; exit 1; }; }
mkdir -p "$OUT"/{dl,work,logs}; OUT=$(realpath "$OUT")

# --- type-strain / reference set (accession -> label), from TYGS Table 4 -------
declare -A ACC=(
  ["GCA_900078665.2"]="houstonense_TYPE_ATCC49403"
  ["GCA_025821245.1"]="farcinogenes_TYPE_DSM43637"
  ["GCF_000723385.1"]="farcinogenes_OURS_454_SUSPECT"
  ["GCA_025822895.1"]="senegalense_TYPE"
  ["GCA_025823105.1"]="porcinum_TYPE"
  ["GCA_001245615.1"]="neworleansense_TYPE"
  ["GCA_010731295.1"]="boenickei_TYPE"
  ["GCA_001052995.1"]="conceptionense_TYPE_D16"
  ["GCA_000455325.1"]="septicum_TYPE_DSM44393"
  ["GCA_000805385.1"]="setense_TYPE"
  ["GCF_002102345.1"]="peregrinum_TYPE"
  ["GCA_000295855.2"]="fortuitum_TYPE_DSM46621"
)
LABELS="$OUT/work/labels.tsv"; : > "$LABELS"
addlabel(){ printf "%s\t%s\n" "$(realpath "$1")" "$2" >> "$LABELS"; }

for acc in "${!ACC[@]}"; do
  d="$OUT/dl/$acc"
  if [[ ! -d "$d" ]]; then
    log "download $acc (${ACC[$acc]})"
    if datasets download genome accession "$acc" --include genome --filename "$OUT/dl/$acc.zip" >"$OUT/dl/$acc.log" 2>&1; then
      unzip -oq "$OUT/dl/$acc.zip" -d "$d"
    else log "WARN: download failed for $acc (try adding/removing the version suffix)"; continue; fi
  fi
  fna=$(find "$d" -name "*.fna" | head -n1)
  [[ -n "$fna" ]] && addlabel "$fna" "${ACC[$acc]}" || log "WARN: no .fna for $acc"
done
# add the cluster representative
cp "$QREP" "$OUT/work/CLUSTER_rep.fna"; addlabel "$OUT/work/CLUSTER_rep.fna" "CLUSTER_rep"

cut -f1 "$LABELS" > "$OUT/work/genomes.txt"
n=$(wc -l < "$OUT/work/genomes.txt"); log "genomes in comparison: $n"
[[ "$n" -lt 3 ]] && { echo "ERROR: too few genomes downloaded (network?)">&2; exit 1; }

# --- sizes / GC (context for the GC-anomaly of the suspect) --------------------
echo -e "label\tsize\tGC" > "$OUT/size_gc.tsv"
while IFS=$'\t' read -r path label; do
  read -r size gc < <(seqkit stats -a -T "$path" 2>/dev/null | awk -F'\t' '
    NR==1{for(i=1;i<=NF;i++){if($i=="sum_len")s=i; if($i ~ /^GC/)g=i}} NR==2{print $s"\t"$g}')
  echo -e "${label}\t${size}\t${gc}" >> "$OUT/size_gc.tsv"
done < "$LABELS"

# --- all-vs-all fastANI --------------------------------------------------------
log "fastANI all-vs-all"
fastANI --ql "$OUT/work/genomes.txt" --rl "$OUT/work/genomes.txt" \
        -o "$OUT/work/ani_raw.tsv" -t "$THREADS" >"$OUT/logs/fastani.log" 2>&1 || true
# labeled long table (path->label), averaged both directions
awk -F'\t' 'FNR==NR{lab[$1]=$2; next}{ if($1 in lab && $2 in lab) print lab[$1]"\t"lab[$2]"\t"$3 }' \
  "$LABELS" "$OUT/work/ani_raw.tsv" > "$OUT/ani_labeled.tsv"

ani(){ # mean ANI between label $1 and $2 (both directions)
  awk -F'\t' -v a="$1" -v b="$2" '($1==a&&$2==b)||($1==b&&$2==a){s+=$3;n++} END{if(n) printf "%.2f",s/n; else printf "NA"}' "$OUT/ani_labeled.tsv"
}
besthit(){ # best non-self hit for label $1 -> "label ANI"
  awk -F'\t' -v a="$1" '$1==a && $2!=a{if($3>b){b=$3;r=$2}} END{printf "%s\t%.2f",r,b}' "$OUT/ani_labeled.tsv"
}

# --- optional CheckM2 ----------------------------------------------------------
CK=""
if [[ $RUN_CHECKM2 -eq 1 ]]; then
  log "CheckM2"
  checkm2 predict --threads "$THREADS" --force --input $(cut -f1 "$LABELS" | tr '\n' ' ') \
    --output-directory "$OUT/checkm2" >"$OUT/logs/checkm2.log" 2>&1 || log "WARN: checkm2 failed"
  CK="$OUT/checkm2/quality_report.tsv"
fi

# --- report --------------------------------------------------------------------
R="$OUT/VERDICT.txt"
{
echo "=================== TAXONOMY VERIFICATION ==================="
echo
echo "Q1. Are M. houstonense and M. farcinogenes the same genomospecies?"
echo "    (species boundary ~95-96% ANI ; 70% dDDH)"
printf "    houstonense_TYPE  vs  farcinogenes_TYPE (GCA_025821245)  = %s%% ANI\n" "$(ani houstonense_TYPE_ATCC49403 farcinogenes_TYPE_DSM43637)"
printf "    CLUSTER_rep       vs  houstonense_TYPE                    = %s%% ANI\n" "$(ani CLUSTER_rep houstonense_TYPE_ATCC49403)"
printf "    CLUSTER_rep       vs  farcinogenes_TYPE (GCA_025821245)   = %s%% ANI\n" "$(ani CLUSTER_rep farcinogenes_TYPE_DSM43637)"
echo
echo "Q2. Is OUR farcinogenes reference (GCF_000723385.1) sound?"
printf "    farcinogenes_OURS vs farcinogenes_TYPE (same strain DSM 43637) = %s%% ANI\n" "$(ani farcinogenes_OURS_454_SUSPECT farcinogenes_TYPE_DSM43637)"
printf "    CLUSTER_rep       vs farcinogenes_OURS (GCF_000723385.1)       = %s%% ANI\n" "$(ani CLUSTER_rep farcinogenes_OURS_454_SUSPECT)"
printf "    farcinogenes_OURS best hit in the whole panel                 = %s\n" "$(besthit farcinogenes_OURS_454_SUSPECT)"
echo "    (two clean assemblies of the same type strain should be >99% ANI;"
echo "     a much lower value indicates the 454 draft is degraded/contaminated)"
echo
echo "Size / GC (a >1% GC gap between two assemblies of one strain is a red flag):"
if command -v column >/dev/null; then column -t -s $'\t' "$OUT/size_gc.tsv"; else cat "$OUT/size_gc.tsv"; fi
if [[ -n "$CK" && -s "$CK" ]]; then
  echo; echo "CheckM2 (completeness / contamination):"
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="Name")n=i;if($i=="Completeness")c=i;if($i=="Contamination")x=i};next}
              {printf "    %-34s comp=%s  contam=%s\n",$n,$c,$x}' "$CK"
fi
echo
echo "Full labeled ANI matrix: $OUT/ani_labeled.tsv"
echo "============================================================"
} | tee "$R"
log "DONE. Verdict written to $R"
