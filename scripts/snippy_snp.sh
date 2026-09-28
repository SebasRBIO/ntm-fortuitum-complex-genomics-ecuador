#!/usr/bin/env bash
#===============================================================================
# snippy_snp.sh
#
# SNP analysis of a (candidate) clonal cluster:
#   reads -> snippy mapping/variant calling -> snippy-core (core alignment)
#         -> snp-dists (pairwise SNP distance matrix)
#         -> cleaned SNP alignment -> IQ-TREE (ML tree, ultrafast bootstrap)
#
# Built for the 7-genome houstonense cluster, but works for any read set.
#
# Usage:
#   ./snippy_snp.sh -i READS_DIR -o OUT_DIR -r REFERENCE [options]
#
# Required:
#   -i DIR   directory with paired reads (*.fastq.gz) for the cluster members ONLY
#   -o DIR   output directory
#   -r FILE  reference genome (FASTA or GenBank). See "Choosing a reference".
#
# Options:
#   -t INT   threads                                   (default: 8)
#   -1 STR   R1 filename tag                            (default: _R1)
#   -2 STR   R2 filename tag                            (default: _R2)
#   -c INT   snippy --mincov (min depth to call)        (default: 10)
#   -s INT   SNP distance to flag a pair as "clonal"    (default: 20)
#   -g       remove recombination with Gubbins before the tree (off)
#   -h       help
#
# Choosing a reference (-r):
#   * Best for accurate intra-cluster SNP distances: a CLOSE, contiguous genome.
#     Options, in order of preference:
#       (a) a complete M. houstonense genome, if available;
#       (b) the most contiguous assembly from within the cluster (e.g. 2409131-001,
#           N50 ~239 kb). If you use an internal assembly as the reference, do NOT
#           also put that isolate's reads in -i (avoids a redundant "Reference" vs.
#           same-isolate tip); its column will appear as "Reference".
#   * A GenBank (.gbk) reference lets snippy annotate variant effects; FASTA works too.
#   * A very divergent reference (e.g. a distant type strain at ~98.5% ANI) shrinks the
#     core and can inflate apparent SNPs — prefer a within-cluster / same-species ref.
#
# Notes:
#   * "Reference" is included as a pseudo-sample in the core alignment; it is kept in
#     the full matrix but EXCLUDED from the intra-cluster clonality statistics.
#   * fastANI flagged clonality; this step QUANTIFIES it in SNPs. Interpret the SNP
#     threshold (-s) with care: outbreak cutoffs are context-specific and, for slow
#     divergence, a handful of SNPs across isolates is consistent with a point source.
#
# Dependencies (bioconda): snippy snp-dists snp-sites iqtree  [optional: gubbins]
#===============================================================================
set -euo pipefail

THREADS=8; R1_TAG="_R1"; R2_TAG="_R2"; MINCOV=10; SNP_FLAG=20; RUN_GUBBINS=0
IN=""; OUT=""; REF=""

usage() { sed -n '2,45p' "$0"; exit "${1:-0}"; }
log() { echo -e "[$(date '+%F %T')] $*"; }

while getopts ":i:o:r:t:1:2:c:s:gh" opt; do
  case $opt in
    i) IN=$OPTARG ;;
    o) OUT=$OPTARG ;;
    r) REF=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    1) R1_TAG=$OPTARG ;;
    2) R2_TAG=$OPTARG ;;
    c) MINCOV=$OPTARG ;;
    s) SNP_FLAG=$OPTARG ;;
    g) RUN_GUBBINS=1 ;;
    h) usage 0 ;;
    \?) echo "Unknown option: -$OPTARG" >&2; usage 1 ;;
    :) echo "Option -$OPTARG needs an argument" >&2; usage 1 ;;
  esac
done
[[ -z "$IN" || -z "$OUT" || -z "$REF" ]] && { echo "ERROR: -i, -o and -r are required." >&2; usage 1; }
[[ -d "$IN" ]] || { echo "ERROR: reads dir not found: $IN" >&2; exit 1; }
[[ -s "$REF" ]] || { echo "ERROR: reference not found: $REF" >&2; exit 1; }

check_tools() {
  local miss=0
  for t in snippy snippy-multi snippy-core snippy-clean_full_aln snp-dists snp-sites; do
    command -v "$t" >/dev/null 2>&1 || { echo "MISSING tool: $t" >&2; miss=1; }
  done
  IQTREE=""
  for c in iqtree2 iqtree; do command -v "$c" >/dev/null 2>&1 && { IQTREE=$c; break; }; done
  [[ -z "$IQTREE" ]] && { echo "MISSING tool: iqtree/iqtree2" >&2; miss=1; }
  [[ $RUN_GUBBINS -eq 1 ]] && { command -v run_gubbins.py >/dev/null 2>&1 || { echo "MISSING (needed for -g): run_gubbins.py" >&2; miss=1; }; }
  [[ $miss -eq 1 ]] && { echo "Install missing tools (bioconda) and re-run." >&2; exit 1; }
}
check_tools

mkdir -p "$OUT"/{snippy,tree,logs}
REF=$(realpath "$REF")

#===============================================================================
# 1) build the snippy-multi input table (ID <TAB> R1 <TAB> R2)
#===============================================================================
TAB="$OUT/snippy/input.tab"; : > "$TAB"
shopt -s nullglob
n=0
for R1 in "$IN"/*"$R1_TAG"*.fastq.gz "$IN"/*"$R1_TAG"*.fq.gz; do
  [[ -e "$R1" ]] || continue
  base=$(basename "$R1"); SAMPLE=${base%%${R1_TAG}*}
  R2=${R1/${R1_TAG}/${R2_TAG}}
  [[ -e "$R2" ]] || { log "WARN: no R2 for $SAMPLE — skipping"; continue; }
  printf "%s\t%s\t%s\n" "$SAMPLE" "$(realpath "$R1")" "$(realpath "$R2")" >> "$TAB"
  n=$((n+1))
done
shopt -u nullglob
[[ $n -lt 2 ]] && { echo "ERROR: need >=2 paired samples in $IN (found $n)." >&2; exit 1; }
log "cluster members: $n"
cut -f1 "$TAB" | sort > "$OUT/snippy/members.txt"

#===============================================================================
# 2) snippy per-sample + snippy-core  (run from the snippy dir)
#===============================================================================
CORE="$OUT/snippy/core.aln"
if [[ ! -s "$CORE" ]]; then
  log "snippy-multi + snippy-core (mincov=$MINCOV)"
  ( cd "$OUT/snippy"
    snippy-multi "$TAB" --ref "$REF" --mincov "$MINCOV" --cpus "$THREADS" > runme.sh 2>"$OUT/logs/snippy-multi.log"
    bash runme.sh >"$OUT/logs/snippy-run.log" 2>&1
  )
fi
[[ -s "$CORE" ]] || { echo "ERROR: snippy-core did not produce core.aln (see $OUT/logs/)."; exit 1; }

#===============================================================================
# 3) pairwise SNP distance matrix (+ molten long form) with snp-dists
#===============================================================================
log "snp-dists"
snp-dists    "$CORE" > "$OUT/snp_distances.matrix.tsv" 2>"$OUT/logs/snp-dists.log"
snp-dists -m "$CORE" > "$OUT/snp_distances.long.tsv"   2>>"$OUT/logs/snp-dists.log"

#===============================================================================
# 4) cleaned SNP alignment for the tree (constant-site-free -> safe for +ASC)
#===============================================================================
log "clean full alignment -> SNP sites"
snippy-clean_full_aln "$OUT/snippy/core.full.aln" > "$OUT/tree/clean.full.aln" 2>"$OUT/logs/clean.log"
ALN_FOR_TREE="$OUT/tree/core.snps.aln"
if [[ $RUN_GUBBINS -eq 1 ]]; then
  log "Gubbins (recombination removal)"
  ( cd "$OUT/tree"
    run_gubbins.py --threads "$THREADS" --prefix gubbins clean.full.aln >"$OUT/logs/gubbins.log" 2>&1
  ) || { echo "WARN: Gubbins failed (see log); falling back to non-recombination-filtered SNPs."; RUN_GUBBINS=0; }
  if [[ $RUN_GUBBINS -eq 1 && -s "$OUT/tree/gubbins.filtered_polymorphic_sites.fasta" ]]; then
    snp-sites -c -o "$ALN_FOR_TREE" "$OUT/tree/gubbins.filtered_polymorphic_sites.fasta" 2>>"$OUT/logs/snp-sites.log" || cp "$OUT/tree/gubbins.filtered_polymorphic_sites.fasta" "$ALN_FOR_TREE"
  fi
fi
if [[ ! -s "$ALN_FOR_TREE" ]]; then
  snp-sites -c -o "$ALN_FOR_TREE" "$OUT/tree/clean.full.aln" 2>>"$OUT/logs/snp-sites.log"
fi
[[ -s "$ALN_FOR_TREE" ]] || { echo "ERROR: empty SNP alignment for tree."; exit 1; }

#===============================================================================
# 5) IQ-TREE ML tree (GTR+G, ascertainment-bias correction, UFBoot)
#===============================================================================
log "IQ-TREE ($IQTREE)"
( cd "$OUT/tree"
  "$IQTREE" -s "$(basename "$ALN_FOR_TREE")" -m GTR+G+ASC -B 1000 -T AUTO --prefix cluster \
    >"$OUT/logs/iqtree.log" 2>&1
) || log "WARN: IQ-TREE failed — see $OUT/logs/iqtree.log (if it reports invariant sites, the SNP alignment may be too small for +ASC)."

#===============================================================================
# 6) clonality summary from the SNP distances (exclude the 'Reference' pseudo-sample)
#===============================================================================
SUM="$OUT/clonality_summary.txt"
awk -F'\t' -v thr="$SNP_FLAG" '
  $1!=$2 && $1!="Reference" && $2!="Reference" {
    d=$3+0; n++; sum+=d;
    if(!seen++){mn=d;mx=d} else {if(d<mn)mn=d; if(d>mx)mx=d}
    if(d>thr && $1<$2){ over[++o]=$1" <-> "$2" : "d } }
  END{
    if(n==0){ print "No non-Reference pairs found."; exit }
    # long form lists each unordered pair twice; report unique pair count
    printf "Intra-cluster SNP distances (excluding Reference):\n";
    printf "  pairs (unordered): %d\n", n/2;
    printf "  min  = %d\n  mean = %.1f\n  max  = %d\n", mn, sum/n, mx;
    if(o>0){ printf "\nPairs above threshold (%d SNPs):\n", thr; for(i=1;i<=o;i++) print "  "over[i] }
    else   { printf "\nAll pairs <= %d SNPs.\n", thr }
  }' "$OUT/snp_distances.long.tsv" > "$SUM"

log "DONE."
echo; echo "===================== SNP CLONALITY SUMMARY ====================="
cat "$SUM"
echo "================================================================"
echo "Matrix:      $OUT/snp_distances.matrix.tsv"
echo "Long form:   $OUT/snp_distances.long.tsv"
echo "Tree:        $OUT/tree/cluster.treefile (+ .contree, UFBoot support)"
echo "SNP align:   $ALN_FOR_TREE"
[[ $RUN_GUBBINS -eq 1 ]] && echo "Recombination-filtered via Gubbins: $OUT/tree/gubbins.*"
echo "Reminder: 'Reference' distances are in the matrix but excluded from the stats above."
