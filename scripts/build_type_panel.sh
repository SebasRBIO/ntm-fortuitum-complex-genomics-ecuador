#!/usr/bin/env bash
#===============================================================================
# build_type_panel.sh
#
# Builds ONE well-anchored M. fortuitum complex type-strain ANI panel that both
# ntm_diagnose.sh (-R) and ntm_type_resistome.sh (-r) can consume, so the whole
# workflow shares a single, citable reference set.
#
# Downloads the complex type strains by ACCESSION (verified type material where
# possible) + the confounder Neisseria mucosa, and emits:
#     <out>/ref_list.txt   (absolute paths, one per line)
#     <out>/labels.tsv     (path <TAB> label)
#
# Usage:  ./build_type_panel.sh -o PANEL_DIR
#
# Confidence of each accession is encoded in its label:
#   *_type    = confirmed assembly from type material (cite as-is)
#   *_confirm = correct species, but VERIFY it is the type strain before citing
#   *_REF...  = pulled as species RefSeq reference (NOT guaranteed type strain)
#
# Dependency: ncbi-datasets-cli (datasets), unzip
#===============================================================================
set -euo pipefail
OUT=""
while getopts ":o:h" opt; do case $opt in
  o) OUT=$OPTARG ;;
  h) sed -n '2,25p' "$0"; exit 0 ;;
  \?) echo "Unknown option: -$OPTARG" >&2; exit 1 ;;
  :) echo "Option -$OPTARG needs an argument" >&2; exit 1 ;;
esac; done
[[ -z "$OUT" ]] && { echo "ERROR: -o PANEL_DIR is required." >&2; exit 1; }
command -v datasets >/dev/null || { echo "ERROR: 'datasets' (ncbi-datasets-cli) not found." >&2; exit 1; }

mkdir -p "$OUT/dl"
LABELS="$OUT/labels.tsv"; : > "$LABELS"
log() { echo -e "[$(date '+%F %T')] $*"; }

# --- Type strains by accession (label carries strain + confidence) -------------
declare -A ACC=(
  ["GCF_022179545.1"]="M_fortuitum__ATCC6841_type"        # anchor, complete genome
  ["GCF_000723385.1"]="M_farcinogenes__NCTC10955_type"
  ["GCF_000455325.1"]="M_septicum__DSM44393_type"
  ["GCF_000805385.1"]="M_setense__CIP109395_type"
  ["GCF_001052995.1"]="M_conceptionense__D16_type"
  ["GCF_019645875.1"]="M_senegalense__ATCC35796_confirm"
  ["GCF_900289185.1"]="M_porcinum__confirm_typematerial"
)
# --- Complex members pulled as species RefSeq reference (confirm strain) --------
declare -A TAXA=(
  ["Mycolicibacterium peregrinum"]="M_peregrinum__REF_confirm_ATCC14467"
  ["Mycolicibacterium houstonense"]="M_houstonense__REF_confirm_ATCC49403"
  ["Mycolicibacterium neworleansense"]="M_neworleansense__REF_confirm_ATCC49404"
  ["Mycolicibacterium boenickei"]="M_boenickei__REF_confirm_ATCC49935"
  ["Mycolicibacterium brisbanense"]="M_brisbanense__REF_confirm_ATCC49938"
)
# --- Confounder ----------------------------------------------------------------
declare -A CONTAM=(
  ["Neisseria mucosa"]="Neisseria_mucosa__contaminant_ref"
)

add_label() { printf "%s\t%s\n" "$(realpath "$1")" "$2" >> "$LABELS"; }

# download each accession independently so one failure doesn't abort the batch
for acc in "${!ACC[@]}"; do
  d="$OUT/dl/$acc"
  if [[ ! -d "$d" ]]; then
    log "accession: $acc (${ACC[$acc]})"
    if datasets download genome accession "$acc" --include genome \
         --filename "$OUT/dl/$acc.zip" >"$OUT/dl/$acc.log" 2>&1; then
      unzip -oq "$OUT/dl/$acc.zip" -d "$d"
    else
      log "WARN: download failed for $acc — check the accession/version and re-run."
      continue
    fi
  fi
  fna=$(find "$d" -name "*.fna" | head -n1)
  [[ -n "$fna" ]] && add_label "$fna" "${ACC[$acc]}" || log "WARN: no .fna for $acc"
done

# species-level references (+ confounder)
download_taxon() {  # $1 = taxon ; $2 = label
  local tax=$1 lab=$2 name d
  name=${1// /_}
  d="$OUT/dl/$name"
  if [[ ! -d "$d" ]]; then
    log "taxon: $tax ($lab)"
    if datasets download genome taxon "$tax" --reference --include genome,seq-report \
         --filename "$OUT/dl/$name.zip" >"$OUT/dl/$name.log" 2>&1; then
      unzip -oq "$OUT/dl/$name.zip" -d "$d"
    else
      log "WARN: no reference for $tax"; return
    fi
  fi
  local fna; fna=$(find "$d" -name "*.fna" | head -n1)
  [[ -n "$fna" ]] && add_label "$fna" "$lab" || log "WARN: no .fna for $tax"
}
for tax in "${!TAXA[@]}"; do download_taxon "$tax" "${TAXA[$tax]}"; done
for tax in "${!CONTAM[@]}"; do download_taxon "$tax" "${CONTAM[$tax]}"; done

cut -f1 "$LABELS" > "$OUT/ref_list.txt"
n=$(wc -l < "$OUT/ref_list.txt")
[[ "$n" -gt 0 ]] || { echo "ERROR: panel is empty (network/datasets?)." >&2; exit 1; }

log "Panel built: $n references"
echo; echo "== $LABELS =="
if command -v column >/dev/null; then column -t -s $'\t' "$LABELS"; else cat "$LABELS"; fi
echo
echo "Consume it with:"
echo "  ./ntm_diagnose.sh       -i RAW -o OUT -R $OUT   [other opts]"
echo "  ./ntm_type_resistome.sh -g GENOMES -o OUT -d BAKTA_DB -r $OUT/ref_list.txt"
echo
echo "Labels ending in _confirm / _REF_confirm: verify the assembly is the type"
echo "strain (see the strain hint in the label) before citing it in Methods."
