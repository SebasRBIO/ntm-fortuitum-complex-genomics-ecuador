#!/usr/bin/env bash
#===============================================================================
# ntm_type_resistome.sh
#
# Step 1 (typing):  type-strain ANI panel -> fastANI -> best-hit assignment
#                   + prep folder for dDDH (TYGS/GGDC, web tools)
# Step 2 (annot+resistome):
#                   Bakta annotation -> AMRFinderPlus + Abricate (acquired genes)
#                   -> TARGETED screen for the M. fortuitum determinants
#                      (blaF, erm(39), aac(2'), sul, dfr) using alleles EXTRACTED
#                      from the annotated M. fortuitum type strain (no hardcoded
#                      sequences) -> rrs/rrl/gyrA extraction for point-mutation review
#
# Input: a directory with your CONFIRMED genomes, one FASTA per isolate
#        (e.g. 240216-MAG-020.fna, 240216-MAG-023.fna, 2408141-027.fna, 2409131-008.fna)
#
# Usage:
#   ./ntm_type_resistome.sh -g GENOME_DIR -o OUT_DIR -d BAKTA_DB [options]
#
# Required:
#   -g DIR   directory with confirmed isolate genomes (*.fna|*.fa|*.fasta)
#   -o DIR   output directory
#   -d DIR   Bakta database directory (the folder that contains db/ ; run once:
#              bakta_db download --output <parent> --type full)
#
# Options:
#   -t INT   threads                              (default: 8)
#   -a FILE  custom curated allele FASTA for the targeted resistome screen.
#            If given, overrides the extract-from-type-strain default for ALL.
#   -H REF   houstonense allele source (accession GCF_/GCA_… or a FASTA path).
#            Isolates assigned to M. houstonense by ANI are screened against
#            determinant alleles extracted from THIS genome instead of the
#            M. fortuitum type strain, removing the fortuitum-centric identity
#            bias for that group. Ignored if -a is given.
#   -r DIR   reuse an external panel from build_type_panel.sh (expects
#            DIR/labels.tsv); skips the built-in download. RECOMMENDED so this
#            step shares the exact same anchored panel as ntm_diagnose.sh.
#   -s       skip type-strain panel download (reuse OUT/refs_type)          (off)
#   -h       help
#
# Type-strain ANI panel (verified accessions + species-level fallbacks).
# CONFIRM the fallback assemblies resolve to the type strain in the table you were
# given before citing them in Methods. Edit PANEL_* below to taste.
#
# Dependencies (bioconda): ncbi-datasets-cli fastani bakta ncbi-amrfinderplus
#                          abricate seqkit barrnap blast
#===============================================================================

set -euo pipefail

#------------------------------- CONFIG ----------------------------------------
THREADS=8
GDIR=""; OUT=""; BAKTA_DB=""; ALLELES=""; SKIP_PANEL=0; EXT_PANEL=""; HOUST_REF=""

# Verified type-strain assemblies (download by accession). Full complex, incl.
# M. farcinogenes. Labels: _type=confirmed type material; _confirm=verify strain.
declare -A PANEL_ACC=(
  ["GCF_022179545.1"]="M_fortuitum__ATCC6841_type"
  ["GCF_000723385.1"]="M_farcinogenes__NCTC10955_type"
  ["GCF_000455325.1"]="M_septicum__DSM44393_type"
  ["GCF_000805385.1"]="M_setense__CIP109395_type"
  ["GCF_001052995.1"]="M_conceptionense__D16_type"
  ["GCF_019645875.1"]="M_senegalense__ATCC35796_confirm"
  ["GCF_900289185.1"]="M_porcinum__confirm_typematerial"
)
# Species pulled as RefSeq reference (CONFIRM strain = type strain after download)
PANEL_TAXA=(
  "Mycolicibacterium peregrinum"
  "Mycolicibacterium houstonense"
  "Mycolicibacterium neworleansense"
  "Mycolicibacterium boenickei"
  "Mycolicibacterium brisbanense"
)
# The assembly used as the resistome allele SOURCE (M. fortuitum type strain)
TYPE_STRAIN_ACC="GCF_022179545.1"

# Targeted determinant keywords (matched against Bakta gene/product, case-insensitive)
# key = short tag used in the summary ; value = extended regex over gene|product
#
# sul / dfr target ACQUIRED sulfonamide/trimethoprim resistance genes only. The
# earlier broad patterns ("dihydropteroate synthase", "dihydrofolate reductase")
# also matched the essential housekeeping folate-pathway genes folP / folA — the
# drug TARGETS present in every mycobacterium, NOT resistance determinants — which
# produced universal false-positive "yes" calls. The patterns below require the
# acquired-gene symbols (sul1-4, dfrA/B/G/K...) or an explicit "…-resistant"
# product. Expected outcome for these isolates: no hit ("-"), i.e. no acquired
# sulfonamide/trimethoprim resistance. CONFIRM any hit against AMRFinderPlus/Abricate
# (OUT/amr/), which are the appropriate curated sources for acquired/mobile genes.
declare -A DETERMINANTS=(
  ["blaF"]="blaF|class A beta-lactamase|beta-lactamase"
  ["erm39"]="erm\\(39\\)|erm39|ermF|rRNA .*methyltransferase|23S rRNA .*adenine.*methyltransferase"
  ["aac"]="aac\\(2|aminoglycoside.*acetyltransferase|acetyltransferase.*aminoglycoside"
  ["sul"]="sul[0-9]|[Ss]ulfonamide[ -]resistant"
  ["dfr"]="dfr[a-zA-Z]|[Tt]rimethoprim[ -]resistant"
)
ID_MIN=85       # blastn %identity threshold for calling a determinant present
COV_MIN=70      # blastn query coverage threshold (%)
#-------------------------------------------------------------------------------

usage() { sed -n '2,52p' "$0"; exit "${1:-0}"; }
log() { echo -e "[$(date '+%F %T')] $*"; }

while getopts ":g:o:d:t:a:r:H:sh" opt; do
  case $opt in
    g) GDIR=$OPTARG ;;
    o) OUT=$OPTARG ;;
    d) BAKTA_DB=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    a) ALLELES=$OPTARG ;;
    r) EXT_PANEL=$OPTARG ;;
    H) HOUST_REF=$OPTARG ;;
    s) SKIP_PANEL=1 ;;
    h) usage 0 ;;
    \?) echo "Unknown option: -$OPTARG" >&2; usage 1 ;;
    :) echo "Option -$OPTARG needs an argument" >&2; usage 1 ;;
  esac
done

[[ -z "$GDIR" || -z "$OUT" || -z "$BAKTA_DB" ]] && { echo "ERROR: -g, -o and -d are required." >&2; usage 1; }
[[ -d "$GDIR" ]] || { echo "ERROR: genome dir not found: $GDIR" >&2; exit 1; }

check_tools() {
  local miss=0
  for t in datasets fastANI bakta amrfinder abricate seqkit barrnap blastn makeblastdb; do
    command -v "$t" >/dev/null 2>&1 || { echo "MISSING tool: $t" >&2; miss=1; }
  done
  [[ $miss -eq 1 ]] && { echo "Install missing tools (bioconda) and re-run." >&2; exit 1; }
}
check_tools

mkdir -p "$OUT"/{refs_type,ani_type,dddh_upload,annot,amr,resistome,logs}

# Collect isolate genomes
shopt -s nullglob
GENOMES=()
for ext in fna fa fasta; do for f in "$GDIR"/*.$ext; do [[ -e "$f" ]] && GENOMES+=("$f"); done; done
shopt -u nullglob
[[ ${#GENOMES[@]} -eq 0 ]] && { echo "No genomes (*.fna|*.fa|*.fasta) in $GDIR" >&2; exit 1; }
log "Isolate genomes: ${#GENOMES[@]}"
sample_of() { local b; b=$(basename "$1"); echo "${b%.*}"; }

#===============================================================================
# STEP 1a: build type-strain panel
#===============================================================================
LABELS="$OUT/refs_type/labels.tsv"
if [[ -n "$EXT_PANEL" ]]; then
  [[ -s "$EXT_PANEL/labels.tsv" ]] || { echo "ERROR: $EXT_PANEL/labels.tsv not found. Run build_type_panel.sh first." >&2; exit 1; }
  log "using external type-strain panel: $EXT_PANEL"
  cp "$EXT_PANEL/labels.tsv" "$LABELS"
elif [[ $SKIP_PANEL -eq 0 || ! -s "$LABELS" ]]; then
  : > "$LABELS"
  cd "$OUT/refs_type"
  # verified accessions
  ACC_LIST=$(printf "%s " "${!PANEL_ACC[@]}")
  log "download type-strain accessions: $ACC_LIST"
  datasets download genome accession $ACC_LIST --include genome \
    --filename acc.zip >"$OUT/logs/datasets_acc.log" 2>&1 || { echo "ERROR: accession download failed (network?)"; exit 1; }
  unzip -oq acc.zip -d acc
  for ACC in "${!PANEL_ACC[@]}"; do
    fna=$(find acc -path "*${ACC}*" -name "*.fna" | head -n1)
    [[ -n "$fna" ]] && printf "%s\t%s\n" "$fna" "${PANEL_ACC[$ACC]}" >> "$LABELS"
  done
  # species-level references (confirm strain later)
  for TAX in "${PANEL_TAXA[@]}"; do
    N=${TAX// /_}
    datasets download genome taxon "$TAX" --reference --include genome,seq-report \
      --filename "$N.zip" >"$OUT/logs/datasets_${N}.log" 2>&1 || { log "WARN: no reference for $TAX"; continue; }
    unzip -oq "$N.zip" -d "$N"
    fna=$(find "$N" -name "*.fna" | head -n1)
    [[ -n "$fna" ]] && printf "%s\t%s\n" "$fna" "${N}__REF_confirm_strain" >> "$LABELS"
  done
  cd - >/dev/null
fi
cut -f1 "$LABELS" > "$OUT/ani_type/ref_list.txt"
[[ -s "$OUT/ani_type/ref_list.txt" ]] || { echo "ERROR: empty type-strain panel."; exit 1; }

#===============================================================================
# STEP 1b: fastANI + best-hit assignment (best vs second-best margin)
#===============================================================================
printf "%s\n" "${GENOMES[@]}" > "$OUT/ani_type/query_list.txt"
log "fastANI vs type-strain panel"
fastANI --ql "$OUT/ani_type/query_list.txt" --rl "$OUT/ani_type/ref_list.txt" \
        -o "$OUT/ani_type/ani.tsv" --matrix -t "$THREADS" >"$OUT/logs/fastani.log" 2>&1 || true

# best_hit(query): prints best_label<TAB>best_ani<TAB>second_label<TAB>second_ani
best_two() {
  local q=$1
  awk -v q="$q" -v lab="$LABELS" '
    BEGIN{ while((getline l<lab)>0){split(l,a,"\t"); L[a[1]]=a[2]} }
    $1==q { if($3>=b1){b2=b1;b2r=b1r;b1=$3;b1r=$2} else if($3>b2){b2=$3;b2r=$2} }
    END{ n1=(b1r in L)?L[b1r]:"NA"; n2=(b2r in L)?L[b2r]:"NA";
         printf "%s\t%.2f\t%s\t%.2f\n", (b1==""?"NA":n1), (b1==""?0:b1),
                                        (b2==""?"NA":n2), (b2==""?0:b2) }' "$OUT/ani_type/ani.tsv"
}

#===============================================================================
# STEP 1c: prepare genomes for dDDH (TYGS / GGDC are web tools -> can't automate)
#===============================================================================
for g in "${GENOMES[@]}"; do cp "$g" "$OUT/dddh_upload/$(sample_of "$g").fasta"; done
cat > "$OUT/dddh_upload/README.txt" <<'TXT'
dDDH is computed on the web:
  TYGS  https://tygs.dsmz.de   -> upload the FASTA files here; it compares against
        ALL type strains automatically and returns dDDH (GBDP) + species + tree.
  GGDC  https://ggdc.dsmz.de   -> pairwise; use formula 2. Same-species cutoff: dDDH >= 70%.
TXT
log "dDDH: genomes staged in $OUT/dddh_upload (upload to TYGS/GGDC)"

#===============================================================================
# STEP 2a: Bakta annotation (isolates + the type strain used for alleles)
#===============================================================================
annotate() {  # $1 = genome fasta ; $2 = sample name ; $3 = species
  local g=$1 s=$2 sp=$3
  [[ -s "$OUT/annot/$s/$s.tsv" ]] && return 0
  log "bakta: $s"
  bakta --db "$BAKTA_DB" --genus Mycolicibacterium --species "$sp" \
        --prefix "$s" --output "$OUT/annot/$s" --threads "$THREADS" --force \
        "$g" >"$OUT/logs/bakta_$s.log" 2>&1
}
for g in "${GENOMES[@]}"; do annotate "$g" "$(sample_of "$g")" "fortuitum"; done

# annotate the type strain (allele source), unless a custom allele FASTA was given
TS_FNA=$(grep -P "\t${PANEL_ACC[$TYPE_STRAIN_ACC]}$" "$LABELS" | cut -f1 | head -n1)
if [[ -z "$ALLELES" ]]; then
  [[ -n "$TS_FNA" ]] || { echo "ERROR: type-strain FASTA not found for allele extraction."; exit 1; }
  annotate "$TS_FNA" "TYPE_fortuitum" "fortuitum"
fi

#===============================================================================
# STEP 2b: acquired-resistance screen (AMRFinderPlus + Abricate)
#===============================================================================
for g in "${GENOMES[@]}"; do
  s=$(sample_of "$g")
  log "AMRFinderPlus + Abricate: $s"
  amrfinder -n "$g" --plus -o "$OUT/amr/${s}.amrfinder.tsv" >"$OUT/logs/amrfinder_$s.log" 2>&1 || log "WARN amrfinder $s"
  for DB in card resfinder ncbi megares; do
    abricate --db "$DB" "$g" > "$OUT/amr/${s}.${DB}.tsv" 2>>"$OUT/logs/abricate_$s.log" || true
  done
done

#===============================================================================
# STEP 2c: TARGETED determinant screen (alleles from type strain or custom FASTA)
#===============================================================================
ALLELE_DB="$OUT/resistome/alleles.ffn"

# --- helper: extract determinant alleles from an annotated reference ---------
# $1 = annotation dir name (under $OUT/annot) ; $2 = output prefix
# writes ${2}.ffn (allele sequences) and ${2}.loci.tsv (locus <TAB> tag)
build_allele_db() {
  local adir=$1 pref=$2 tag
  local tsv="$OUT/annot/$adir/$adir.tsv" ffn="$OUT/annot/$adir/$adir.ffn"
  : > "${pref}.loci.tsv"
  for tag in "${!DETERMINANTS[@]}"; do
    awk -F'\t' -v re="${DETERMINANTS[$tag]}" -v t="$tag" 'BEGIN{IGNORECASE=1}
      $0 ~ re { print $6"\t"t }' "$tsv" >> "${pref}.loci.tsv" 2>/dev/null || true
  done
  cut -f1 "${pref}.loci.tsv" | sort -u > "${pref}.loci_ids.txt"
  [[ -s "${pref}.loci_ids.txt" ]] && \
    seqkit grep -f "${pref}.loci_ids.txt" "$ffn" > "${pref}.ffn" 2>>"$OUT/logs/seqkit.log" || true
}

# --- fortuitum allele source (default) or a global custom FASTA via -a --------
FORT_PREF="$OUT/resistome/alleles_fortuitum"
if [[ -n "$ALLELES" ]]; then
  cp "$ALLELES" "${FORT_PREF}.ffn"
  # name custom records erm39/blaF/aac/sul/dfr so they map to the summary columns
  seqkit fx2tab -n -i "${FORT_PREF}.ffn" 2>/dev/null | awk '{print $1"\t"$1}' > "${FORT_PREF}.loci.tsv" || : > "${FORT_PREF}.loci.tsv"
else
  log "extracting reference alleles from M. fortuitum type strain"
  build_allele_db "TYPE_fortuitum" "$FORT_PREF"
  [[ -s "${FORT_PREF}.ffn" ]] || log "WARN: no fortuitum alleles extracted (check keywords/Bakta naming)."
fi

# --- houstonense allele source (optional, -H): removes fortuitum bias ---------
HOUST_PREF=""
if [[ -n "$HOUST_REF" && -z "$ALLELES" ]]; then
  log "preparing houstonense allele source: $HOUST_REF"
  HREF_FNA=""
  if [[ "$HOUST_REF" =~ ^GC[AF]_ ]]; then
    mkdir -p "$OUT/refs_houst"
    if datasets download genome accession "$HOUST_REF" --include genome \
         --filename "$OUT/refs_houst/h.zip" >"$OUT/logs/datasets_houst.log" 2>&1; then
      unzip -oq "$OUT/refs_houst/h.zip" -d "$OUT/refs_houst"
      HREF_FNA=$(find "$OUT/refs_houst" -name "*.fna" | head -n1)
    else
      log "WARN: could not download $HOUST_REF (see log)."
    fi
  else
    HREF_FNA="$HOUST_REF"
  fi
  if [[ -s "$HREF_FNA" ]]; then
    annotate "$HREF_FNA" "TYPE_houstonense" "houstonense"
    HOUST_PREF="$OUT/resistome/alleles_houstonense"
    build_allele_db "TYPE_houstonense" "$HOUST_PREF"
    [[ -s "${HOUST_PREF}.ffn" ]] || { log "WARN: no houstonense alleles extracted; using fortuitum alleles for that group."; HOUST_PREF=""; }
  else
    log "WARN: houstonense reference not found; using fortuitum alleles for all groups."
  fi
fi

# --- blastn each isolate vs the allele source matching its ANI species group --
declare -A HIT   # HIT[sample|tag]=identity
: > "$OUT/resistome/allele_source.tsv"
for g in "${GENOMES[@]}"; do
  s=$(sample_of "$g")
  read -r GB1 _ _ _ < <(best_two "$g")
  if [[ -n "$HOUST_PREF" && "$GB1" == *houstonense* ]]; then
    APREF="$HOUST_PREF"; SRC="houstonense"
  else
    APREF="$FORT_PREF"; SRC="fortuitum"
  fi
  printf "%s\t%s\n" "$s" "$SRC" >> "$OUT/resistome/allele_source.tsv"
  [[ -s "${APREF}.ffn" ]] || { log "WARN: no allele DB for $s ($SRC-source)"; continue; }
  makeblastdb -in "$g" -dbtype nucl -out "$OUT/resistome/${s}_db" >/dev/null 2>&1
  blastn -query "${APREF}.ffn" -db "$OUT/resistome/${s}_db" \
    -outfmt '6 qseqid pident length qlen qcovs sstart send' \
    -perc_identity "$ID_MIN" > "$OUT/resistome/${s}.blast.tsv" 2>>"$OUT/logs/blast_$s.log" || true
  while IFS=$'\t' read -r qid pid len qlen qcov ss se; do
    tag=$(awk -v q="$qid" '$1==q{print $2; exit}' "${APREF}.loci.tsv")
    [[ -z "$tag" ]] && continue
    awk -v c="$qcov" -v m="$COV_MIN" 'BEGIN{exit !(c>=m)}' && HIT["$s|$tag"]="$pid"
  done < "$OUT/resistome/${s}.blast.tsv"
done

#===============================================================================
# STEP 2d: rRNA (rrs/rrl) + gyrA extraction for point-mutation review
#===============================================================================
for g in "${GENOMES[@]}"; do
  s=$(sample_of "$g")
  barrnap --kingdom bac --threads "$THREADS" --outseq "$OUT/resistome/${s}.rRNA.fa" \
    "$g" >"$OUT/logs/barrnap_$s.log" 2>&1 || true
  # gyrA/gyrB from Bakta CDS
  awk -F'\t' 'BEGIN{IGNORECASE=1} $0 ~ /gyrA|gyrB|DNA gyrase/ {print $6}' \
    "$OUT/annot/$s/$s.tsv" | sort -u > "$OUT/resistome/${s}.gyr_loci.txt" || true
  [[ -s "$OUT/resistome/${s}.gyr_loci.txt" ]] && \
    seqkit grep -f "$OUT/resistome/${s}.gyr_loci.txt" "$OUT/annot/$s/$s.ffn" \
      > "$OUT/resistome/${s}.gyr.fa" 2>>"$OUT/logs/seqkit.log" || true
done
cat > "$OUT/resistome/POINT_MUTATIONS_README.txt" <<'TXT'
Point-mutation review (align isolate genes to the type-strain alleles and inspect):
  rrl (23S): macrolides -> positions 2058/2059 (E. coli numbering)
  rrs (16S): amikacin   -> position 1408 (E. coli numbering)
  gyrA (QRDR): fluoroquinolones -> codons ~83/87 (E. coli numbering)
Per isolate: <sample>.rRNA.fa (16S/23S from barrnap) and <sample>.gyr.fa (gyrA/gyrB CDS).
Suggested: mafft-align each gene across isolates + type strain, then map to E. coli
coordinates to call the exact residue. Automated coordinate calling is intentionally
left manual to avoid mis-numbering.
TXT

#===============================================================================
# Final summary
#===============================================================================
SUM="$OUT/typing_resistome_summary.tsv"
echo -e "sample\tbest_species\tbest_ANI\t2nd_species\t2nd_ANI\tANI_margin\term39\tblaF\taac\tsul\tdfr\tallele_src" > "$SUM"
for g in "${GENOMES[@]}"; do
  s=$(sample_of "$g")
  read -r B1 A1 B2 A2 < <(best_two "$g")
  MARGIN=$(awk -v a="$A1" -v b="$A2" 'BEGIN{printf "%.2f", a-b}')
  SRC=$(awk -v s="$s" '$1==s{print $2; exit}' "$OUT/resistome/allele_source.tsv" 2>/dev/null); SRC=${SRC:-fortuitum}
  cell() { [[ -n "${HIT[$s|$1]:-}" ]] && echo "yes(${HIT[$s|$1]}%)" || echo "-"; }
  echo -e "${s}\t${B1}\t${A1}\t${B2}\t${A2}\t${MARGIN}\t$(cell erm39)\t$(cell blaF)\t$(cell aac)\t$(cell sul)\t$(cell dfr)\t${SRC}" >> "$SUM"
done

log "DONE."
echo; echo "===================== SUMMARY ($SUM) ====================="
column -t -s $'\t' "$SUM"
echo "========================================================================"
echo "ANI: best_species should be M. fortuitum with ANI>=95-96 and a clear margin."
echo "dDDH: upload $OUT/dddh_upload/*.fasta to TYGS/GGDC (>=70% = same species)."
echo "Acquired genes: $OUT/amr/  |  Targeted determinants: $OUT/resistome/"
echo "Point mutations: see $OUT/resistome/POINT_MUTATIONS_README.txt"
