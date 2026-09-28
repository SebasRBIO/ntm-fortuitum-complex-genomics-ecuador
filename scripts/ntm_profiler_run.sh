#!/usr/bin/env bash
#===============================================================================
# ntm_profiler_run.sh
#
# Orthogonal NTM species identification with NTM-Profiler (companion to
# ntm_type_resistome.sh). Runs on paired reads and collates per-sample calls.
#
# ROLE / SCOPE (read this):
#   * Used here for SPECIES prediction only — an independent cross-check to
#     fastANI and dDDH (TYGS). NTM-Profiler species prediction is k-mer based,
#     with a mash-vs-GTDB fallback, so it covers the M. fortuitum complex.
#   * NTM-Profiler's built-in RESISTANCE database currently covers only
#     M. leprae and the M. abscessus complex (abscessus/bolletii/massiliense).
#     It therefore returns NO resistance calls for M. fortuitum / M. houstonense.
#     The resistome for these species comes from ntm_type_resistome.sh
#     (Bakta + AMRFinderPlus/Abricate + targeted blaF/erm(39)/aac screen).
#   * NTM-Profiler is in alpha; treat its output as a complementary line of
#     evidence, not a sole authority.
#   * To obtain NTM-Profiler resistance calls for the fortuitum complex you would
#     need a custom resistance DB (ntm-profiler create_resistance_db): a reference
#     genome, gff, variables.json, and a curated Gene/Mutation/Drug CSV. Out of
#     scope here; flagged for later.
#
# Usage:
#   ./ntm_profiler_run.sh -i READS_DIR -o OUT_DIR [options]
#
# Required:
#   -i DIR   directory with paired reads (*.fastq.gz) for the isolates
#   -o DIR   output directory
#
# Options:
#   -t INT   threads                          (default: 8)
#   -1 STR   R1 filename tag                   (default: _R1)
#   -2 STR   R2 filename tag                   (default: _R2)
#   -u       run 'ntm-profiler update_db' first (one-time DB download)   (off)
#   -h       help
#
# Dependency (bioconda): ntm-profiler   (then: ntm-profiler update_db)
#===============================================================================
set -euo pipefail

THREADS=8; R1_TAG="_R1"; R2_TAG="_R2"; UPDATE_DB=0
IN=""; OUT=""

usage() { sed -n '2,45p' "$0"; exit "${1:-0}"; }
log() { echo -e "[$(date '+%F %T')] $*"; }

while getopts ":i:o:t:1:2:uh" opt; do
  case $opt in
    i) IN=$OPTARG ;;
    o) OUT=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    1) R1_TAG=$OPTARG ;;
    2) R2_TAG=$OPTARG ;;
    u) UPDATE_DB=1 ;;
    h) usage 0 ;;
    \?) echo "Unknown option: -$OPTARG" >&2; usage 1 ;;
    :) echo "Option -$OPTARG needs an argument" >&2; usage 1 ;;
  esac
done
[[ -z "$IN" || -z "$OUT" ]] && { echo "ERROR: -i and -o are required." >&2; usage 1; }
[[ -d "$IN" ]] || { echo "ERROR: reads dir not found: $IN" >&2; exit 1; }
command -v ntm-profiler >/dev/null 2>&1 || { echo "MISSING tool: ntm-profiler (conda install bioconda::ntm-profiler)" >&2; exit 1; }

mkdir -p "$OUT"/logs
OUT=$(realpath "$OUT"); IN=$(realpath "$IN")

if [[ $UPDATE_DB -eq 1 ]]; then
  log "ntm-profiler update_db (one-time)"
  ntm-profiler update_db >"$OUT/logs/update_db.log" 2>&1 || { echo "ERROR: update_db failed (see log)."; exit 1; }
fi

#===============================================================================
# per-sample profiling (species; resistance only for supported species)
#===============================================================================
shopt -s nullglob
n=0
for R1 in "$IN"/*"$R1_TAG"*.fastq.gz "$IN"/*"$R1_TAG"*.fq.gz; do
  [[ -e "$R1" ]] || continue
  base=$(basename "$R1"); SAMPLE=${base%%${R1_TAG}*}
  R2=${R1/${R1_TAG}/${R2_TAG}}
  [[ -e "$R2" ]] || { log "WARN: no R2 for $SAMPLE — skipping"; continue; }
  n=$((n+1))
  if [[ -s "$OUT/results/${SAMPLE}.results.json" || -s "$OUT/results/${SAMPLE}.results.txt" ]]; then
    log "[$SAMPLE] already profiled — skipping"; continue
  fi
  log "[$SAMPLE] ntm-profiler profile"
  ( cd "$OUT"
    ntm-profiler profile -1 "$R1" -2 "$R2" -p "$SAMPLE" --threads "$THREADS" --txt \
      >"$OUT/logs/${SAMPLE}.ntmp.log" 2>&1
  ) || log "WARN: profiling failed for $SAMPLE (see $OUT/logs/${SAMPLE}.ntmp.log)"
done
shopt -u nullglob
[[ $n -eq 0 ]] && { echo "No paired reads found in $IN with tags $R1_TAG/$R2_TAG." >&2; exit 1; }
log "profiled samples: $n"

#===============================================================================
# collate all runs into one table
#===============================================================================
log "ntm-profiler collate"
( cd "$OUT" && ntm-profiler collate --prefix ntmprofiler >"$OUT/logs/collate.log" 2>&1 ) \
  || log "WARN: collate failed (see $OUT/logs/collate.log)"

log "DONE."
echo; echo "===================== NTM-PROFILER RESULTS ====================="
COLL=""
for f in "$OUT/ntmprofiler.txt" "$OUT/ntmprofiler.csv"; do [[ -s "$f" ]] && { COLL=$f; break; }; done
if [[ -n "$COLL" ]]; then
  echo "Collated report: $COLL"
  echo "----------------------------------------------------------------"
  if command -v column >/dev/null && [[ "$COLL" == *.txt ]]; then column -t -s $'\t' "$COLL"; else cat "$COLL"; fi
else
  echo "No collated report found; per-sample reports are in $OUT/results/"
fi
echo "================================================================"
echo "Use the species column as an independent cross-check against fastANI/dDDH."
echo "Resistance is expected to be EMPTY for M. fortuitum / M. houstonense"
echo "(NTM-Profiler resistance DB = M. leprae + M. abscessus complex only)."
echo "Resistome for these species: see ntm_type_resistome.sh output."
