#!/usr/bin/env bash
#===============================================================================
# context_tree.sh
#
# Places your focal genomes in the global context of their species:
#   download global set (NCBI, by taxon) -> QC filter -> [optional dereplicate]
#   -> annotate (Prokka) -> core genome (Panaroo | Roary) -> IQ-TREE (ML + UFBoot)
#
# The "global set" is fetched automatically by taxon — no manual curation needed.
# Defaults to M. fortuitum (taxid 1766) + M. houstonense (taxid 146021).
#
# Usage:
#   ./context_tree.sh -g FOCAL_DIR -o OUT_DIR [options]
#
# Required:
#   -g DIR   directory with YOUR focal genomes (*.fna|*.fa|*.fasta), one per isolate
#   -o DIR   output directory
#
# Options:
#   -t INT   threads                                   (default: 8)
#   -T LIST  comma-separated taxa to download           (default: "Mycolicibacterium fortuitum,Mycolicibacterium houstonense")
#   -s SRC   assembly source: RefSeq | GenBank | all    (default: RefSeq)
#   -O ACC   outgroup accession (GCF_/GCA_…) for rooting (e.g. GCF_000015005.1 M. smegmatis) (default: none)
#   -A LIST  extra reference type strains to include, as acc:label pairs, comma-separated
#            and added by ACCESSION (bypass QC). Default anchors the houstonense/
#            farcinogenes conspecificity: M. farcinogenes DSM 43637 + M. senegalense DSM 43656.
#            IMPORTANT: do NOT add farcinogenes via -T (by taxon) — the public
#            'M. farcinogenes' set contains a mislabeled assembly (GCF_000723385.1,
#            which is actually M. senegalense); pulling by exact accession avoids it.
#   -D FLOAT dereplicate context genomes at this Mash distance (e.g. 0.001 ≈ 99.9% ANI); keep 1/cluster (default: off)
#   -j INT   parallel Prokka annotation jobs (default: THREADS/2); each job uses THREADS/j cpus.
#            This is usually the biggest speed-up. Lower it if RAM-limited (~1-2 GB per job).
#   -R       use Roary instead of Panaroo               (default: Panaroo)
#   -h       help
#
# QC filter (edit in CONFIG): size window, min N50, max contigs. Drops junk public assemblies.
#
# Notes:
#   * Annotation + core-genome on many genomes is the heavy step. Prefer -s RefSeq
#     and -D to keep the set tractable; run on a multi-core machine.
#   * Panaroo is preferred (robust to annotation error / fragmented assemblies).
#   * The core alignment has variable AND conserved sites, so IQ-TREE runs WITHOUT
#     ascertainment-bias correction (unlike a SNP-only alignment).
#
# Dependencies (bioconda): ncbi-datasets-cli seqkit prokka panaroo iqtree  [+ mash if -D] [+ roary if -R]
#===============================================================================
set -euo pipefail

#------------------------------- CONFIG ----------------------------------------
THREADS=8
TAXA_CSV="Mycolicibacterium fortuitum,Mycolicibacterium houstonense"
SRC="RefSeq"
OUTGROUP=""
REF_ACC_CSV="GCA_025821245.1:M_farcinogenes_TYPE_DSM43637,GCA_025822895.1:M_senegalense_TYPE_DSM43656"
DEREP=""
USE_ROARY=0
PAR_JOBS=""     # parallel Prokka jobs (default: THREADS/2); cpus/job = THREADS/PAR_JOBS
GDIR=""; OUT=""
# QC thresholds for public assemblies
MINSIZE=5000000; MAXSIZE=8000000; MIN_N50=10000; MAX_CONTIGS=1000
#-------------------------------------------------------------------------------

usage() { sed -n '2,45p' "$0"; exit "${1:-0}"; }
log() { echo -e "[$(date '+%F %T')] $*"; }

while getopts ":g:o:t:T:s:O:A:D:j:Rh" opt; do
  case $opt in
    g) GDIR=$OPTARG ;;
    o) OUT=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    T) TAXA_CSV=$OPTARG ;;
    s) SRC=$OPTARG ;;
    O) OUTGROUP=$OPTARG ;;
    A) REF_ACC_CSV=$OPTARG ;;
    D) DEREP=$OPTARG ;;
    j) PAR_JOBS=$OPTARG ;;
    R) USE_ROARY=1 ;;
    h) usage 0 ;;
    \?) echo "Unknown option: -$OPTARG" >&2; usage 1 ;;
    :) echo "Option -$OPTARG needs an argument" >&2; usage 1 ;;
  esac
done
[[ -z "$GDIR" || -z "$OUT" ]] && { echo "ERROR: -g and -o are required." >&2; usage 1; }
[[ -d "$GDIR" ]] || { echo "ERROR: focal dir not found: $GDIR" >&2; exit 1; }

check_tools() {
  local miss=0
  for t in datasets seqkit prokka iqtree2; do
    command -v "$t" >/dev/null 2>&1 || { [[ "$t" == iqtree2 ]] && command -v iqtree >/dev/null 2>&1 || { echo "MISSING tool: $t" >&2; miss=1; }; }
  done
  IQTREE=$(command -v iqtree2 || command -v iqtree)
  if [[ $USE_ROARY -eq 1 ]]; then command -v roary >/dev/null 2>&1 || { echo "MISSING (needed for -R): roary" >&2; miss=1; }
  else command -v panaroo >/dev/null 2>&1 || { echo "MISSING: panaroo" >&2; miss=1; }; fi
  [[ -n "$DEREP" ]] && { command -v mash >/dev/null 2>&1 || { echo "MISSING (needed for -D): mash" >&2; miss=1; }; }
  [[ $miss -eq 1 ]] && { echo "Install missing tools (bioconda) and re-run." >&2; exit 1; }
}
check_tools

mkdir -p "$OUT"/{download,context_raw,work,annot,logs}
OUT=$(realpath "$OUT")
META="$OUT/genome_metadata.tsv"; echo -e "genome\tcategory\tsource" > "$META"
STATS="$OUT/context_qc.tsv"; echo -e "genome\tsize\tN50\tcontigs\tqc" > "$STATS"

sanitize() { echo "$1" | tr -c 'A-Za-z0-9_.-' '_' | sed 's/\.[^.]*$//; s/[._-]\+$//'; }

#===============================================================================
# 1) download global set by taxon
#===============================================================================
IFS=',' read -ra TAXA <<< "$TAXA_CSV"
srcflag=""; [[ "$SRC" != "all" ]] && srcflag="--assembly-source $SRC"
for TAX in "${TAXA[@]}"; do
  TAX="$(echo "$TAX" | sed 's/^ *//; s/ *$//')"; N=${TAX// /_}
  if [[ ! -d "$OUT/download/$N" ]]; then
    log "download: $TAX ($SRC)"
    datasets download genome taxon "$TAX" $srcflag --include genome \
      --filename "$OUT/download/$N.zip" >"$OUT/logs/dl_$N.log" 2>&1 || { log "WARN: download failed for $TAX"; continue; }
    unzip -oq "$OUT/download/$N.zip" -d "$OUT/download/$N"
  fi
  find "$OUT/download/$N" -name "*.fna" -exec cp {} "$OUT/context_raw/" \;
done
shopt -s nullglob
CTX=("$OUT"/context_raw/*.fna)
shopt -u nullglob
[[ ${#CTX[@]} -eq 0 ]] && { echo "ERROR: no context genomes downloaded (network/datasets?)." >&2; exit 1; }
log "context genomes downloaded: ${#CTX[@]}"

#===============================================================================
# 2) QC filter public assemblies (size / N50 / contigs)
#===============================================================================
PASS_DIR="$OUT/context_pass"; mkdir -p "$PASS_DIR"
for f in "${CTX[@]}"; do
  read -r size n50 ctg < <(seqkit stats -a -T "$f" 2>/dev/null | awk -F'\t' '
    NR==1{for(i=1;i<=NF;i++){if($i=="sum_len")cs=i; if($i=="N50")cn=i; if($i=="num_seqs")cc=i}}
    NR==2{print $cs"\t"$cn"\t"$cc}')
  size=${size:-0}; n50=${n50:-0}; ctg=${ctg:-999999}
  gname=$(basename "$f" .fna)
  if awk -v s="$size" -v n="$n50" -v c="$ctg" -v smin="$MINSIZE" -v smax="$MAXSIZE" \
         -v nmin="$MIN_N50" -v cmax="$MAX_CONTIGS" \
         'BEGIN{exit !(s>=smin && s<=smax && n>=nmin && c<=cmax)}'; then
    cp "$f" "$PASS_DIR/"; echo -e "${gname}\t${size}\t${n50}\t${ctg}\tPASS" >> "$STATS"
  else
    echo -e "${gname}\t${size}\t${n50}\t${ctg}\tDROP" >> "$STATS"
  fi
done
shopt -s nullglob
PASS=("$PASS_DIR"/*.fna); shopt -u nullglob
log "context genomes passing QC: ${#PASS[@]} (dropped $(( ${#CTX[@]} - ${#PASS[@]} )))"
[[ ${#PASS[@]} -eq 0 ]] && { echo "ERROR: no context genome passed QC (relax thresholds?)." >&2; exit 1; }

#===============================================================================
# 3) optional Mash dereplication of the context set (keep best N50 per cluster)
#===============================================================================
REP_DIR="$PASS_DIR"
if [[ -n "$DEREP" ]]; then
  log "dereplicating context at Mash distance $DEREP"
  REP_DIR="$OUT/context_derep"; mkdir -p "$REP_DIR"
  mash sketch -p "$THREADS" -o "$OUT/work/ctx" "${PASS[@]}" >"$OUT/logs/mash.log" 2>&1
  mash dist -p "$THREADS" "$OUT/work/ctx.msh" "$OUT/work/ctx.msh" > "$OUT/work/mashdist.tsv" 2>>"$OUT/logs/mash.log"
  # basename pairs (exclude self), single-linkage clustering at DEREP
  awk 'BEGIN{OFS="\t"}{a=$1;b=$2; sub(/.*\//,"",a);sub(/\.[^.]*$/,"",a); sub(/.*\//,"",b);sub(/\.[^.]*$/,"",b); if(a!=b) print a,b,$3}' \
    "$OUT/work/mashdist.tsv" > "$OUT/work/mash_pairs.tsv"
  ls "$PASS_DIR" | sed 's/\.fna$//' > "$OUT/work/ctx_names.txt"
  awk -v thr="$DEREP" '
    function find(x, r,nx){r=x;while(parent[r]!=r)r=parent[r];while(parent[x]!=r){nx=parent[x];parent[x]=r;x=nx}return r}
    FNR==NR{if(!($1 in parent))parent[$1]=$1; nodes[$1]=1; next}
    {a=$1;b=$2;d=$3+0; if(!(a in parent))parent[a]=a; if(!(b in parent))parent[b]=b; nodes[a]=1;nodes[b]=1;
     if(d<=thr){ra=find(a);rb=find(b); if(ra!=rb)parent[rb]=ra}}
    END{for(n in nodes) print n"\t"find(n)}' \
    "$OUT/work/ctx_names.txt" "$OUT/work/mash_pairs.tsv" > "$OUT/work/ctx_clusters.tsv"
  # pick best N50 representative per cluster (N50 from STATS)
  awk -F'\t' 'FNR==NR{n50[$1]=$3; next} {print $2"\t"$1"\t"(($1 in n50)?n50[$1]:0)}' \
    "$STATS" "$OUT/work/ctx_clusters.tsv" | sort -k1,1 -k3,3nr \
    | awk -F'\t' '!seen[$1]++{print $2}' > "$OUT/work/ctx_reps.txt"
  while read -r rep; do cp "$PASS_DIR/${rep}.fna" "$REP_DIR/" 2>/dev/null || true; done < "$OUT/work/ctx_reps.txt"
  shopt -s nullglob; REPS=("$REP_DIR"/*.fna); shopt -u nullglob
  log "context representatives after dereplication: ${#REPS[@]}"
fi

#===============================================================================
# 4) assemble working set: focal + context(+reps) + outgroup, sanitized names
#===============================================================================
add_genome() {  # $1 = fasta ; $2 = category
  local f=$1 cat=$2 name; name=$(sanitize "$(basename "$f")")
  cp "$f" "$OUT/work/${name}.fna"
  echo -e "${name}\t${cat}\t$(basename "$f")" >> "$META"
}
shopt -s nullglob
for ext in fna fa fasta; do for f in "$GDIR"/*.$ext; do [[ -e "$f" ]] && add_genome "$f" focal; done; done
for f in "$REP_DIR"/*.fna; do add_genome "$f" context; done
shopt -u nullglob

# reference type strains added by ACCESSION (bypass QC) to anchor the taxonomy.
# Default: M. farcinogenes DSM 43637 + M. senegalense DSM 43656, so the tree SHOWS
# the houstonense/farcinogenes conspecificity rather than only asserting it.
if [[ -n "$REF_ACC_CSV" ]]; then
  IFS=',' read -ra REFS <<< "$REF_ACC_CSV"
  for entry in "${REFS[@]}"; do
    acc="${entry%%:*}"; lab="${entry##*:}"; [[ "$lab" == "$acc" ]] && lab="ref_$acc"
    d="$OUT/download/ref_$acc"
    if [[ ! -d "$d" ]]; then
      log "download reference $acc ($lab)"
      if datasets download genome accession "$acc" --include genome --filename "$OUT/download/ref_$acc.zip" >"$OUT/logs/dl_ref_$acc.log" 2>&1; then
        unzip -oq "$OUT/download/ref_$acc.zip" -d "$d"
      else log "WARN: reference download failed for $acc"; continue; fi
    fi
    rfna=$(find "$d" -name "*.fna" | head -n1)
    if [[ -n "$rfna" ]]; then
      tmp="$OUT/download/${lab}.fna"; cp "$rfna" "$tmp"; add_genome "$tmp" reference
    else log "WARN: no .fna for reference $acc"; fi
  done
fi

OUTNAME=""
if [[ -n "$OUTGROUP" ]]; then
  log "downloading outgroup $OUTGROUP"
  datasets download genome accession "$OUTGROUP" --include genome --filename "$OUT/download/outgroup.zip" >"$OUT/logs/dl_outgroup.log" 2>&1 \
    && unzip -oq "$OUT/download/outgroup.zip" -d "$OUT/download/outgroup" \
    && { ofna=$(find "$OUT/download/outgroup" -name "*.fna" | head -n1); [[ -n "$ofna" ]] && { add_genome "$ofna" outgroup; OUTNAME=$(sanitize "$(basename "$ofna")"); }; } \
    || log "WARN: outgroup download failed."
fi
shopt -s nullglob; WORK=("$OUT"/work/*.fna); shopt -u nullglob
log "total genomes in tree: ${#WORK[@]}"

#===============================================================================
# 5) annotate all with Prokka -> GFF (Panaroo/Roary-compatible), IN PARALLEL
#===============================================================================
# Prokka over many genomes is embarrassingly parallel. Instead of one Prokka at a
# time using all cores, run PAR_JOBS Prokka jobs concurrently, each with
# CPUS_PER = THREADS/PAR_JOBS cores. Default PAR_JOBS = THREADS/2 (>=1).
# Lower PAR_JOBS if RAM-limited (~1-2 GB per concurrent Prokka).
[[ -z "$PAR_JOBS" ]] && PAR_JOBS=$(( THREADS/2 )); [[ "$PAR_JOBS" -lt 1 ]] && PAR_JOBS=1
CPUS_PER=$(( THREADS/PAR_JOBS )); [[ "$CPUS_PER" -lt 1 ]] && CPUS_PER=1
log "Prokka annotation: $PAR_JOBS parallel jobs x $CPUS_PER cpus each"

# worklist: only genomes not yet annotated (idempotent / resumable)
: > "$OUT/work/annot_todo.txt"
for f in "${WORK[@]}"; do
  s=$(basename "$f" .fna)
  [[ -s "$OUT/annot/$s/$s.gff" ]] || printf '%s\n' "$f" >> "$OUT/work/annot_todo.txt"
done

run_prokka() {  # $1 = genome fasta
  local f="$1" s; s=$(basename "$f" .fna)
  prokka --outdir "$OUT/annot/$s" --prefix "$s" --genus Mycolicibacterium \
         --cpus "$CPUS_PER" --force "$f" >"$OUT/logs/prokka_$s.log" 2>&1 \
    || echo "[WARN] prokka failed for $s (see logs/prokka_$s.log)"
}
export -f run_prokka; export OUT CPUS_PER

n_todo=$(wc -l < "$OUT/work/annot_todo.txt")
if [[ "$n_todo" -gt 0 ]]; then
  log "annotating $n_todo genomes..."
  xargs -a "$OUT/work/annot_todo.txt" -P "$PAR_JOBS" -I{} bash -c 'run_prokka "$@"' _ {} \
    2>>"$OUT/logs/prokka_parallel.log" || log "WARN: some Prokka jobs reported errors (see logs)"
else
  log "all genomes already annotated; skipping."
fi
shopt -s nullglob; GFFS=("$OUT"/annot/*/*.gff); shopt -u nullglob
[[ ${#GFFS[@]} -lt 4 ]] && { echo "ERROR: <4 GFFs produced; cannot build a tree." >&2; exit 1; }
log "GFFs ready: ${#GFFS[@]}"

#===============================================================================
# 6) core genome: Panaroo (default) or Roary
#===============================================================================
if [[ $USE_ROARY -eq 1 ]]; then
  log "Roary"
  roary -e --mafft -p "$THREADS" -f "$OUT/roary" "${GFFS[@]}" >"$OUT/logs/roary.log" 2>&1 || true
  CORE_ALN=$(find "$OUT/roary" -name "core_gene_alignment.aln" | head -n1)
else
  log "Panaroo"
  panaroo -i "${GFFS[@]}" -o "$OUT/panaroo" --clean-mode strict --aligner mafft -a core -t "$THREADS" \
    >"$OUT/logs/panaroo.log" 2>&1 || true
  CORE_ALN=$(find "$OUT/panaroo" -name "core_gene_alignment_filtered.aln" -o -name "core_gene_alignment.aln" | head -n1)
fi
[[ -s "$CORE_ALN" ]] || { echo "ERROR: no core gene alignment produced (see logs)." >&2; exit 1; }
log "core alignment: $CORE_ALN"

#===============================================================================
# 7) IQ-TREE (ModelFinder + ultrafast bootstrap)
#===============================================================================
log "IQ-TREE"
oflag=""; [[ -n "$OUTNAME" ]] && oflag="-o $OUTNAME"
( cd "$OUT" && "$IQTREE" -s "$CORE_ALN" -m MFP -B 1000 -T AUTO $oflag --prefix context_tree \
    >"$OUT/logs/iqtree.log" 2>&1 ) || log "WARN: IQ-TREE failed (see logs/iqtree.log)."

log "DONE."
echo; echo "===================== CONTEXT TREE ====================="
echo "Genomes in tree : ${#WORK[@]}  (focal + context$([[ -n "$OUTNAME" ]] && echo " + outgroup"))"
echo "Metadata        : $META   (category = focal/context/outgroup — use to colour tips)"
echo "QC table        : $STATS"
echo "Core alignment  : $CORE_ALN"
echo "Tree            : $OUT/context_tree.treefile   (UFBoot support; .contree consensus)"
echo "Open the .treefile in iTOL/FigTree and colour by the 'category' column of the metadata."
