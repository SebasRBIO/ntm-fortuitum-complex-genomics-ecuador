#!/usr/bin/env bash
#===============================================================================
# checkm2_ani.sh
#
# Two things in one pass over a set of confirmed genomes:
#   1) CheckM2  -> completeness / contamination (MIMAG HQ check)
#   2) all-vs-all fastANI -> pairwise ANI + single-linkage clustering at a
#      "clonal" threshold, to test whether the houstonense group (or any group)
#      is a clonal cluster (point-source) vs. diverse isolates.
#
# Output (in OUT/):
#   checkm2/quality_report.tsv     raw CheckM2 output
#   ani/pairs.tsv                  every genome pair, ANI, sorted high->low
#   ani/ani.matrix                 fastANI matrix (phylip-like)
#   clusters.tsv                   sample -> clonal_cluster_id
#   summary.tsv                    sample | completeness | contamination | MIMAG | cluster | cluster_size
#   cluster_ani_stats.tsv          per multi-member cluster: min/mean/max intra-cluster ANI
#
# Usage:
#   ./checkm2_ani.sh -g GENOME_DIR -o OUT_DIR [options]
#
# Required:
#   -g DIR   directory with confirmed genomes (*.fna|*.fa|*.fasta), one per isolate
#   -o DIR   output directory
#
# Options:
#   -t INT   threads                                   (default: 8)
#   -n FLOAT ANI threshold to call a pair "clonal"     (default: 99.9)
#   -D PATH  CheckM2 diamond DB (.dmnd); sets CHECKM2DB (else uses configured one)
#   -w       download the CheckM2 DB first (one-time, large)             (off)
#   -h       help
#
# Notes:
#   * fastANI saturates near ~99.9%, so it flags likely-clonal groups but is not a
#     SNP count. To quantify transmission (SNP distances) use a mapping-based SNP
#     pipeline (e.g. snippy) on each cluster afterwards.
#   * Run this on ALL usable genomes at once (e.g. 5 fortuitum + 7 houstonense) so
#     cross-group pairs stay low-ANI and only true clonal groups collapse together.
#
# Dependencies (bioconda): checkm2 fastani seqkit
#===============================================================================
set -euo pipefail

THREADS=8
CLONAL_ANI=99.9
GDIR=""; OUT=""; CHECKM2_DB=""; DOWNLOAD_DB=0

usage() { sed -n '2,45p' "$0"; exit "${1:-0}"; }
log() { echo -e "[$(date '+%F %T')] $*"; }

while getopts ":g:o:t:n:D:wh" opt; do
  case $opt in
    g) GDIR=$OPTARG ;;
    o) OUT=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    n) CLONAL_ANI=$OPTARG ;;
    D) CHECKM2_DB=$OPTARG ;;
    w) DOWNLOAD_DB=1 ;;
    h) usage 0 ;;
    \?) echo "Unknown option: -$OPTARG" >&2; usage 1 ;;
    :) echo "Option -$OPTARG needs an argument" >&2; usage 1 ;;
  esac
done
[[ -z "$GDIR" || -z "$OUT" ]] && { echo "ERROR: -g and -o are required." >&2; usage 1; }
[[ -d "$GDIR" ]] || { echo "ERROR: genome dir not found: $GDIR" >&2; exit 1; }

for t in checkm2 fastANI; do
  command -v "$t" >/dev/null 2>&1 || { echo "MISSING tool: $t (bioconda)" >&2; exit 1; }
done
[[ -n "$CHECKM2_DB" ]] && export CHECKM2DB="$CHECKM2_DB"

mkdir -p "$OUT"/{checkm2,ani,logs}

# gather genomes + a stable sample list
shopt -s nullglob
GENOMES=()
for ext in fna fa fasta; do for f in "$GDIR"/*.$ext; do [[ -e "$f" ]] && GENOMES+=("$f"); done; done
shopt -u nullglob
[[ ${#GENOMES[@]} -eq 0 ]] && { echo "No genomes (*.fna|*.fa|*.fasta) in $GDIR" >&2; exit 1; }
sample_of() { local b; b=$(basename "$1"); echo "${b%.*}"; }
: > "$OUT/ani/genomes.txt"; : > "$OUT/ani/samples.txt"
for g in "${GENOMES[@]}"; do echo "$g" >> "$OUT/ani/genomes.txt"; sample_of "$g" >> "$OUT/ani/samples.txt"; done
log "Genomes: ${#GENOMES[@]}"

#===============================================================================
# 1) CheckM2
#===============================================================================
if [[ $DOWNLOAD_DB -eq 1 ]]; then
  log "CheckM2: downloading database (one-time)"
  checkm2 database --download >"$OUT/logs/checkm2_db.log" 2>&1 || { echo "ERROR: DB download failed (see log)"; exit 1; }
fi
QR="$OUT/checkm2/quality_report.tsv"
if [[ ! -s "$QR" ]]; then
  log "CheckM2 predict"
  checkm2 predict --threads "$THREADS" --force \
    --input "${GENOMES[@]}" --output-directory "$OUT/checkm2" \
    >"$OUT/logs/checkm2.log" 2>&1 || { echo "ERROR: checkm2 predict failed — see $OUT/logs/checkm2.log"; exit 1; }
fi
[[ -s "$QR" ]] || { echo "ERROR: no CheckM2 quality_report.tsv produced."; exit 1; }

#===============================================================================
# 2) all-vs-all fastANI
#===============================================================================
log "fastANI all-vs-all"
fastANI --ql "$OUT/ani/genomes.txt" --rl "$OUT/ani/genomes.txt" \
        -o "$OUT/ani/ani_raw.tsv" --matrix -t "$THREADS" >"$OUT/logs/fastani.log" 2>&1 || true
[[ -f "$OUT/ani/ani_raw.tsv.matrix" ]] && mv -f "$OUT/ani/ani_raw.tsv.matrix" "$OUT/ani/ani.matrix"

# named pairs (basename, exclude self), sorted by ANI desc
awk 'BEGIN{OFS="\t"}
  { a=$1; b=$2; sub(/.*\//,"",a); sub(/\.[^.]*$/,"",a);
             sub(/.*\//,"",b); sub(/\.[^.]*$/,"",b);
    if(a!=b) print a,b,$3 }' "$OUT/ani/ani_raw.tsv" \
  | sort -k3,3nr > "$OUT/ani/pairs.tsv"

#===============================================================================
# 3) single-linkage clustering at the clonal threshold (union-find)
#===============================================================================
awk -v thr="$CLONAL_ANI" '
  function find(x,   r,nx){ r=x; while(parent[r]!=r) r=parent[r];
    while(parent[x]!=r){ nx=parent[x]; parent[x]=r; x=nx } return r }
  FNR==NR{ s=$1; if(!(s in parent)) parent[s]=s; nodes[s]=1; next }
  { a=$1; b=$2; ani=$3+0;
    if(!(a in parent)) parent[a]=a; if(!(b in parent)) parent[b]=b;
    nodes[a]=1; nodes[b]=1;
    if(ani>=thr){ ra=find(a); rb=find(b); if(ra!=rb) parent[rb]=ra } }
  END{ cid=0;
    for(n in nodes){ r=find(n); if(!(r in cl)){ cid++; cl[r]=cid } }
    for(n in nodes){ root=find(n); print n"\t"cl[root] } }
' "$OUT/ani/samples.txt" "$OUT/ani/pairs.tsv" | sort -k2,2n -k1,1 > "$OUT/clusters.tsv"

# cluster sizes
awk -F'\t' '{c[$2]++} END{for(k in c) print k"\t"c[k]}' "$OUT/clusters.tsv" | sort -k1,1n > "$OUT/ani/cluster_sizes.tsv"

# per multi-member cluster: intra-cluster ANI min/mean/max
awk -F'\t' '
  FNR==NR{ cl[$1]=$2; next }
  { a=$1; b=$2; ani=$3+0;
    if((a in cl)&&(b in cl)&&cl[a]==cl[b]){
      c=cl[a]; n[c]++; sum[c]+=ani;
      if(!(c in mn)||ani<mn[c]) mn[c]=ani;
      if(!(c in mx)||ani>mx[c]) mx[c]=ani } }
  END{ print "cluster\tn_pairs\tmin_ANI\tmean_ANI\tmax_ANI";
    for(c in n) printf "%s\t%d\t%.3f\t%.3f\t%.3f\n", c, n[c], mn[c], sum[c]/n[c], mx[c] }
' "$OUT/clusters.tsv" "$OUT/ani/pairs.tsv" | sort -k1,1n > "$OUT/cluster_ani_stats.tsv"

#===============================================================================
# 4) combined summary (join CheckM2 + cluster)
#===============================================================================
# CheckM2 quality_report columns: Name Completeness Contamination ... (tab-sep, header)
SUM="$OUT/summary.tsv"
awk -F'\t' -v clf="$OUT/clusters.tsv" -v szf="$OUT/ani/cluster_sizes.tsv" '
  BEGIN{
    while((getline l < clf)>0){ split(l,a,"\t"); CL[a[1]]=a[2] }
    while((getline l < szf)>0){ split(l,a,"\t"); SZ[a[1]]=a[2] }
    print "sample\tcompleteness\tcontamination\tMIMAG_HQ\tcluster\tcluster_size"
  }
  NR==1{ for(i=1;i<=NF;i++){ if($i=="Name")cN=i; if($i=="Completeness")cC=i; if($i=="Contamination")cX=i } next }
  {
    name=$cN; comp=$cC+0; cont=$cX+0;
    hq=(comp>90 && cont<5)?"pass":"FLAG";
    c=(name in CL)?CL[name]:"NA"; sz=(c in SZ)?SZ[c]:1;
    printf "%s\t%.2f\t%.2f\t%s\t%s\t%s\n", name, comp, cont, hq, c, sz
  }
' "$QR" | sort -k5,5n -k1,1 > "$SUM"

#===============================================================================
log "DONE."
echo; echo "===================== QUALITY + CLONALITY ($SUM) ====================="
if command -v column >/dev/null; then column -t -s $'\t' "$SUM"; else cat "$SUM"; fi
echo
echo "Intra-cluster ANI (clonal threshold = ${CLONAL_ANI}%):"
if command -v column >/dev/null; then column -t -s $'\t' "$OUT/cluster_ani_stats.tsv"; else cat "$OUT/cluster_ani_stats.tsv"; fi
echo
echo "Reading it:  a multi-member cluster with min_ANI >= ${CLONAL_ANI} = candidate clonal"
echo "point-source group. MIMAG_HQ=FLAG means completeness<=90 or contamination>=5."
echo "Top pairs: $OUT/ani/pairs.tsv   |   full matrix: $OUT/ani/ani.matrix"
