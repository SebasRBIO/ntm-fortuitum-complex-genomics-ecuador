#!/usr/bin/env bash
#===============================================================================
# ntm_diagnose.sh
#
# Diagnostic triage for paired-end WGS of putative NTM isolates that gave a
# discordant Kraken/DRAGEN re-identification (e.g. MALDI = M. fortuitum vs.
# reads classified as Neisseria mucosa).
#
# Per sample it runs:  QC/trim -> assembly -> size+GC -> fastANI -> verdict
# and, only for MIXED samples, optional taxon-directed read recovery + reassembly.
# It writes one summary TSV with a per-sample verdict:
#   MYCOBACTERIUM  | NEISSERIA | MIXED | RECOVERED_MYCOBACTERIUM | REVIEW
#
# Dependencies (one conda/mamba env):
#   fastp shovill fastani seqkit ncbi-datasets-cli   [optional: kraken2 krakentools checkm-genome quast]
#
# Usage:
#   ./ntm_diagnose.sh -i RAW_DIR -o OUT_DIR [options]
#
# Required:
#   -i DIR   directory with paired fastq.gz
#   -o DIR   output directory (created if absent)
#
# Options:
#   -t INT   threads                              (default: 8)
#   -m INT   RAM in GB for shovill                (default: 32)
#   -1 STR   R1 filename tag                      (default: _R1)
#   -2 STR   R2 filename tag                      (default: _R2)
#   -k DIR   Kraken2 DB  -> enables MIXED recovery (default: none = recovery off)
#   -x INT   taxid to extract in recovery         (default: 1866885 = Mycolicibacterium;
#                                                   use 1763 if your DB keeps the old
#                                                   "Mycobacterium" taxonomy)
#   -c       run CheckM (completeness/contamination)  (default: off)
#   -R DIR   use an external type-strain panel built by build_type_panel.sh
#            (expects DIR/labels.tsv); skips the built-in taxon download.
#            RECOMMENDED so the whole workflow shares one anchored panel.
#   -h       help
#
# Notes:
#   * Default references are pulled with `datasets ... --reference` (ONE
#     representative per taxon, not necessarily the type strain). For manuscript-
#     grade ANI, build an anchored type-strain panel with build_type_panel.sh and
#     pass it via -R (the whole complex, incl. M. farcinogenes, is anchored there).
#   * Thresholds are variables at the top of the CONFIG block — tune to taste.
#===============================================================================

set -euo pipefail

#------------------------------- CONFIG ----------------------------------------
THREADS=8
RAM=32
R1_TAG="_R1"
R2_TAG="_R2"
KRAKEN_DB=""
MYCO_TAXID=1866885
RUN_CHECKM=0
RAW=""
OUT=""
EXT_PANEL=""

# Reference taxa for the ANI panel (fortuitum complex + the confounder)
REF_TAXA=(
  "Mycolicibacterium fortuitum"
  "Mycolicibacterium peregrinum"
  "Mycolicibacterium porcinum"
  "Mycolicibacterium senegalense"
  "Neisseria mucosa"
)

# Decision thresholds
ANI_SPECIES=95.0        # >= this ANI to call a species
MYCO_GC_LO=63.0         # mycobacterial GC window
MYCO_GC_HI=68.0
MYCO_SIZE_LO=5200000    # mycobacterial size window (bp)
MYCO_SIZE_HI=7200000
NEIS_GC_LO=48.0
NEIS_GC_HI=54.0
NEIS_SIZE_LO=1900000
NEIS_SIZE_HI=3200000
HIGHGC_CUT=59.0         # per-contig GC cut separating the two organisms
MIX_FRAC_LO=0.15        # if high-GC fraction is between LO and HI -> bimodal -> MIXED
MIX_FRAC_HI=0.85
MAPFRAC_MIN=0.60        # min fastANI mapped-fragment fraction for a confident call
#-------------------------------------------------------------------------------

usage() { sed -n '2,40p' "$0"; exit "${1:-0}"; }

while getopts ":i:o:t:m:1:2:k:x:cR:h" opt; do
  case $opt in
    i) RAW=$OPTARG ;;
    o) OUT=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    m) RAM=$OPTARG ;;
    1) R1_TAG=$OPTARG ;;
    2) R2_TAG=$OPTARG ;;
    k) KRAKEN_DB=$OPTARG ;;
    x) MYCO_TAXID=$OPTARG ;;
    c) RUN_CHECKM=1 ;;
    R) EXT_PANEL=$OPTARG ;;
    h) usage 0 ;;
    \?) echo "Unknown option: -$OPTARG" >&2; usage 1 ;;
    :) echo "Option -$OPTARG needs an argument" >&2; usage 1 ;;
  esac
done

[[ -z "$RAW" || -z "$OUT" ]] && { echo "ERROR: -i and -o are required." >&2; usage 1; }
[[ -d "$RAW" ]] || { echo "ERROR: input dir not found: $RAW" >&2; exit 1; }

log() { echo -e "[$(date '+%F %T')] $*"; }

check_tools() {
  local missing=0
  for t in fastp shovill fastANI seqkit datasets; do
    command -v "$t" >/dev/null 2>&1 || { echo "MISSING tool: $t" >&2; missing=1; }
  done
  if [[ -n "$KRAKEN_DB" ]]; then
    for t in kraken2 extract_kraken_reads.py; do
      command -v "$t" >/dev/null 2>&1 || { echo "MISSING (needed for -k): $t" >&2; missing=1; }
    done
  fi
  [[ $RUN_CHECKM -eq 1 ]] && { command -v checkm >/dev/null 2>&1 || { echo "MISSING (needed for -c): checkm" >&2; missing=1; }; }
  [[ $missing -eq 1 ]] && { echo "Install missing tools and re-run." >&2; exit 1; }
}

mkdir -p "$OUT"/{trimmed,assembly,stats,refs,ani,myco_reads,reassembly,checkm,logs}
check_tools

#--- helper: size, overall GC, and high-GC fraction of one assembly ------------
# prints: total_len<TAB>gc_pct<TAB>high_gc_frac
asm_stats() {
  local fa=$1
  seqkit fx2tab -nlg "$fa" | awk -v cut="$HIGHGC_CUT" '
    { len+=$2; gc+=$2*$3/100; if($3>=cut) hi+=$2 }
    END { if(len==0){print "0\t0\t0"; exit}
          printf "%d\t%.2f\t%.4f\n", len, (gc/len)*100, hi/len }'
}

#--- helper: best fastANI hit for a query path ---------------------------------
# args: query_fasta_path  ani_results.tsv  ref_labels.tsv
# prints: label<TAB>ani<TAB>mapped_frac   (or "NA 0 0")
best_ani() {
  local q=$1 res=$2 labels=$3
  awk -v q="$q" -v lab="$labels" '
    BEGIN{ while((getline l < lab)>0){ split(l,a,"\t"); L[a[1]]=a[2] } }
    $1==q { af=$4/$5; if($3>best){best=$3; bref=$2; bfrac=af} }
    END{ if(best==""){print "NA\t0\t0"} else {
           name=(bref in L)?L[bref]:bref
           printf "%s\t%.2f\t%.4f\n", name, best, bfrac } }' "$res"
}

#===============================================================================
# STEP 1-2: trim + assemble each sample
#===============================================================================
shopt -s nullglob
SAMPLES=()
for R1 in "$RAW"/*"$R1_TAG"*.fastq.gz "$RAW"/*"$R1_TAG"*.fq.gz; do
  [[ -e "$R1" ]] || continue
  base=$(basename "$R1")
  SAMPLE=${base%%${R1_TAG}*}
  R2=${R1/${R1_TAG}/${R2_TAG}}
  [[ -e "$R2" ]] || { log "WARN: no R2 for $SAMPLE ($R2) — skipping"; continue; }
  SAMPLES+=("$SAMPLE")

  tR1="$OUT/trimmed/${SAMPLE}_R1.trim.fastq.gz"
  tR2="$OUT/trimmed/${SAMPLE}_R2.trim.fastq.gz"
  if [[ ! -s "$tR1" ]]; then
    log "[$SAMPLE] fastp"
    fastp -i "$R1" -I "$R2" -o "$tR1" -O "$tR2" -w "$THREADS" \
      -j "$OUT/stats/${SAMPLE}.fastp.json" -h "$OUT/stats/${SAMPLE}.fastp.html" \
      >"$OUT/logs/${SAMPLE}.fastp.log" 2>&1
  fi

  asm="$OUT/assembly/${SAMPLE}/contigs.fa"
  if [[ ! -s "$asm" ]]; then
    log "[$SAMPLE] shovill"
    shovill --R1 "$tR1" --R2 "$tR2" --outdir "$OUT/assembly/${SAMPLE}" \
      --cpus "$THREADS" --ram "$RAM" --force \
      >"$OUT/logs/${SAMPLE}.shovill.log" 2>&1
  fi
done
shopt -u nullglob
[[ ${#SAMPLES[@]} -eq 0 ]] && { echo "No paired samples found with tags $R1_TAG/$R2_TAG in $RAW" >&2; exit 1; }
log "Samples: ${SAMPLES[*]}"

#===============================================================================
# STEP 4a: build reference panel (cached)
#===============================================================================
LABELS="$OUT/refs/ref_labels.tsv"
if [[ -n "$EXT_PANEL" ]]; then
  # Use the shared, well-anchored panel from build_type_panel.sh
  [[ -s "$EXT_PANEL/labels.tsv" ]] || { echo "ERROR: $EXT_PANEL/labels.tsv not found. Run build_type_panel.sh first." >&2; exit 1; }
  log "using external type-strain panel: $EXT_PANEL"
  cp "$EXT_PANEL/labels.tsv" "$LABELS"
elif [[ ! -s "$LABELS" ]]; then
  : > "$LABELS"
  for TAX in "${REF_TAXA[@]}"; do
    NAME=${TAX// /_}
    zip="$OUT/refs/${NAME}.zip"
    if [[ ! -d "$OUT/refs/${NAME}" ]]; then
      log "download ref: $TAX"
      datasets download genome taxon "$TAX" --reference --include genome \
        --filename "$zip" >"$OUT/logs/datasets_${NAME}.log" 2>&1 || { log "WARN: download failed for $TAX"; continue; }
      unzip -oq "$zip" -d "$OUT/refs/${NAME}"
    fi
    while IFS= read -r fna; do
      printf "%s\t%s\n" "$fna" "$NAME" >> "$LABELS"
    done < <(find "$OUT/refs/${NAME}" -name "*.fna")
  done
fi
cut -f1 "$LABELS" > "$OUT/ani/ref_list.txt"
[[ -s "$OUT/ani/ref_list.txt" ]] || { echo "ERROR: no reference genomes downloaded (check network/datasets)." >&2; exit 1; }

#===============================================================================
# STEP 3+4: stats + fastANI on primary assemblies -> verdict
#===============================================================================
: > "$OUT/ani/query_list.txt"
for S in "${SAMPLES[@]}"; do echo "$OUT/assembly/${S}/contigs.fa" >> "$OUT/ani/query_list.txt"; done

log "fastANI (primary assemblies)"
fastANI --ql "$OUT/ani/query_list.txt" --rl "$OUT/ani/ref_list.txt" \
        -o "$OUT/ani/ani_results.tsv" -t "$THREADS" \
        >"$OUT/logs/fastani.log" 2>&1 || true

SUMMARY="$OUT/summary.tsv"
echo -e "sample\tstage\ttotal_len_bp\tGC_pct\thigh_gc_frac\tbest_hit\tANI\tmapped_frac\tverdict" > "$SUMMARY"

# verdict function shared by primary and recovered assemblies
decide() {
  local len=$1 gc=$2 hifrac=$3 hit=$4 ani=$5 mfrac=$6
  awk -v len="$len" -v gc="$gc" -v hf="$hifrac" -v hit="$hit" -v ani="$ani" -v mf="$mfrac" \
      -v aS="$ANI_SPECIES" -v gLo="$MYCO_GC_LO" -v gHi="$MYCO_GC_HI" \
      -v sLo="$MYCO_SIZE_LO" -v sHi="$MYCO_SIZE_HI" \
      -v nLo="$NEIS_GC_LO" -v nHi="$NEIS_GC_HI" -v nsLo="$NEIS_SIZE_LO" -v nsHi="$NEIS_SIZE_HI" \
      -v mixLo="$MIX_FRAC_LO" -v mixHi="$MIX_FRAC_HI" -v mfMin="$MAPFRAC_MIN" 'BEGIN{
    myco = (index(hit,"Mycolicibacterium")>0 && ani>=aS)
    neis = (index(hit,"Neisseria")>0 && ani>=aS)
    bimodal = (hf>mixLo && hf<mixHi)
    inflated = (len>sHi)
    if (bimodal || inflated || (myco && mf<mfMin)) { print "MIXED"; exit }
    if (myco && gc>=gLo && gc<=gHi && len>=sLo && len<=sHi) { print "MYCOBACTERIUM"; exit }
    if (neis && gc>=nLo && gc<=nHi && len>=nsLo && len<=nsHi) { print "NEISSERIA"; exit }
    print "REVIEW"
  }'
}

declare -a MIXED_SAMPLES=()
for S in "${SAMPLES[@]}"; do
  fa="$OUT/assembly/${S}/contigs.fa"
  read -r LEN GC HF < <(asm_stats "$fa")
  read -r HIT ANI MF < <(best_ani "$fa" "$OUT/ani/ani_results.tsv" "$LABELS")
  V=$(decide "$LEN" "$GC" "$HF" "$HIT" "$ANI" "$MF")
  echo -e "${S}\tprimary\t${LEN}\t${GC}\t${HF}\t${HIT}\t${ANI}\t${MF}\t${V}" >> "$SUMMARY"
  [[ "$V" == "MIXED" ]] && MIXED_SAMPLES+=("$S")
done

#===============================================================================
# STEP 5: recover mycobacterial fraction for MIXED samples (only if -k given)
#===============================================================================
if [[ ${#MIXED_SAMPLES[@]} -gt 0 ]]; then
  if [[ -z "$KRAKEN_DB" ]]; then
    log "MIXED samples found (${MIXED_SAMPLES[*]}) but no -k DB given: skipping recovery."
  else
    for S in "${MIXED_SAMPLES[@]}"; do
      log "[$S] recovery: kraken2 + extract taxid $MYCO_TAXID"
      tR1="$OUT/trimmed/${S}_R1.trim.fastq.gz"; tR2="$OUT/trimmed/${S}_R2.trim.fastq.gz"
      kraken2 --db "$KRAKEN_DB" --paired "$tR1" "$tR2" --threads "$THREADS" \
        --report "$OUT/myco_reads/${S}.kreport" --output "$OUT/myco_reads/${S}.kout" \
        >"$OUT/logs/${S}.kraken.log" 2>&1
      extract_kraken_reads.py -k "$OUT/myco_reads/${S}.kout" -r "$OUT/myco_reads/${S}.kreport" \
        -s1 "$tR1" -s2 "$tR2" \
        -o "$OUT/myco_reads/${S}_myco_R1.fq" -o2 "$OUT/myco_reads/${S}_myco_R2.fq" \
        --taxid "$MYCO_TAXID" --include-children --fastq-output \
        >"$OUT/logs/${S}.extract.log" 2>&1

      nread=$( ( [[ -s "$OUT/myco_reads/${S}_myco_R1.fq" ]] && wc -l < "$OUT/myco_reads/${S}_myco_R1.fq" ) || echo 0 )
      if [[ "$nread" -lt 4 ]]; then
        log "[$S] no mycobacterial reads extracted — leaving as MIXED"
        echo -e "${S}\trecovery\t0\t0\t0\tNA\t0\t0\tMIXED_NO_RECOVERY" >> "$SUMMARY"
        continue
      fi

      shovill --R1 "$OUT/myco_reads/${S}_myco_R1.fq" --R2 "$OUT/myco_reads/${S}_myco_R2.fq" \
        --outdir "$OUT/reassembly/${S}" --cpus "$THREADS" --ram "$RAM" --force \
        >"$OUT/logs/${S}.reassembly.log" 2>&1

      rfa="$OUT/reassembly/${S}/contigs.fa"
      echo "$rfa" > "$OUT/ani/${S}.requery.txt"
      fastANI --ql "$OUT/ani/${S}.requery.txt" --rl "$OUT/ani/ref_list.txt" \
        -o "$OUT/ani/${S}.reani.tsv" -t "$THREADS" >>"$OUT/logs/fastani.log" 2>&1 || true
      read -r LEN GC HF < <(asm_stats "$rfa")
      read -r HIT ANI MF < <(best_ani "$rfa" "$OUT/ani/${S}.reani.tsv" "$LABELS")
      V=$(decide "$LEN" "$GC" "$HF" "$HIT" "$ANI" "$MF")
      [[ "$V" == "MYCOBACTERIUM" ]] && V="RECOVERED_MYCOBACTERIUM"
      echo -e "${S}\trecovery\t${LEN}\t${GC}\t${HF}\t${HIT}\t${ANI}\t${MF}\t${V}" >> "$SUMMARY"
    done
  fi
fi

#===============================================================================
# Optional: CheckM on all primary assemblies
#===============================================================================
if [[ $RUN_CHECKM -eq 1 ]]; then
  log "CheckM (completeness/contamination)"
  bins="$OUT/checkm/bins"; mkdir -p "$bins"
  for S in "${SAMPLES[@]}"; do cp "$OUT/assembly/${S}/contigs.fa" "$bins/${S}.fna"; done
  checkm lineage_wf -x fna -t "$THREADS" "$bins" "$OUT/checkm/out" \
    --tab_table -f "$OUT/checkm/checkm_results.tsv" >"$OUT/logs/checkm.log" 2>&1 || log "WARN: CheckM failed (see log)"
fi

#===============================================================================
log "DONE."
echo
echo "======================= SUMMARY ($SUMMARY) ======================="
column -t -s $'\t' "$SUMMARY"
echo "=================================================================="
echo "Legend: MYCOBACTERIUM=clean isolate (supports 'contaminated reference DB');"
echo "        NEISSERIA=real contamination; MIXED=bimodal/inflated assembly;"
echo "        RECOVERED_MYCOBACTERIUM=rescued after read extraction; REVIEW=inspect manually."
