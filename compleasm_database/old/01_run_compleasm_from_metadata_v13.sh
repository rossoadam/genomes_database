#!/usr/bin/env bash
###############################################################################
# 01_run_compleasm_from_metadata_v13.sh
#
# PURPOSE
#   Run Compleasm end-to-end (including internal miniprot) on genomes listed
#   in a metadata CSV, extract CDS sequences, and maintain a compleasm/metadata.csv
#   mapping species to outputs.
#
# KEY CHANGE vs older versions:
#   - NO manual miniprot calls.
#   - Uses: compleasm run ... (Compleasm runs miniprot internally)
#
# OUTPUT METADATA (metadata.csv)
#   Headers:
#     genus_species,accession,organism_name,lineage,full_table,cds_fasta
#
# USAGE
#   bash 01_run_compleasm_from_metadata_v13.sh <query|allgenomes> [options]
#
# ARGUMENTS
#   <query>       Accession OR organism name (case-insensitive, underscores ok)
#   allgenomes    Run Compleasm on all genomes in the metadata CSV
#
# OPTIONS
#   -t INT        Threads (default: 8)
#   -l STR        Lineage folder name (default: sauropsida_odb12)
#   -f            Force re-run (archives existing lineage dir first)
#   -h            Help
#
# ENV OVERRIDES (optional)
#   METADATA_CSV=...
#   LINEAGES_DIR=...
#   OUT_ROOT=...
#   METADATA_OUT=...
#   GET_CDS_PY=...
#   HMMSEARCH_BIN=...
#   EXPECTED_ENV=compleasm_env
#
###############################################################################

set -euo pipefail

# ---- Defaults (match your current layout) ----
METADATA_CSV="${METADATA_CSV:-/Users/rossoaa/projects/genomes/records/genomes_metadata.csv}"
LINEAGES_DIR="${LINEAGES_DIR:-/Users/rossoaa/projects/genomes/records/compleasm/busco_lineages/mb_downloads}"
OUT_ROOT="${OUT_ROOT:-/Users/rossoaa/projects/genomes/records/compleasm}"

LINEAGE_ID="${LINEAGE_ID:-sauropsida_odb12}"
THREADS="${THREADS:-8}"
FORCE_RERUN="${FORCE_RERUN:-0}"

METADATA_OUT="${METADATA_OUT:-${OUT_ROOT}/metadata.csv}"

# Your CDS extraction script (same idea as before)
GET_CDS_PY="${GET_CDS_PY:-$(dirname "$0")/02a_get_cds_from_compleasm_v6.py}"

# Optional: point to hmmsearch explicitly; otherwise rely on PATH inside env
HMMSEARCH_BIN="${HMMSEARCH_BIN:-}"

# Optional: purely a warning (no legacy env logic anymore)
EXPECTED_ENV="${EXPECTED_ENV:-compleasm_env}"

QUERY="${1:-}"
shift || true

while getopts ":t:l:fh" opt; do
  case "$opt" in
    t) THREADS="$OPTARG" ;;
    l) LINEAGE_ID="$OPTARG" ;;
    f) FORCE_RERUN=1 ;;
    h)
      cat <<EOF
Usage: $(basename "$0") <query|allgenomes> [-t threads] [-l lineage] [-f]

Defaults:
  -l sauropsida_odb12
  -t 8

Notes:
  - This script uses "compleasm run" (internal miniprot). No manual miniprot.
  - It will SKIP genomes that already have:
      <outdir>/<lineage>/miniprot_output.gff
      <outdir>/<lineage>/full_table.tsv (or .csv)
      <outdir>/<lineage>/<genus_species>_cds_compleasm.fasta
    unless you use -f.

Env overrides:
  METADATA_CSV=..., LINEAGES_DIR=..., OUT_ROOT=..., METADATA_OUT=...
  GET_CDS_PY=..., HMMSEARCH_BIN=..., EXPECTED_ENV=...

EOF
      exit 0
      ;;
    \?) echo "ERROR: invalid option -$OPTARG" >&2; exit 2 ;;
    :)  echo "ERROR: option -$OPTARG requires an argument" >&2; exit 2 ;;
  esac
done

if [[ -z "$QUERY" ]]; then
  echo "ERROR: missing argument (<query> or allgenomes)" >&2
  exit 2
fi

mkdir -p "$OUT_ROOT" "$LINEAGES_DIR"

# ---- Sanity checks ----
if [[ ! -d "${LINEAGES_DIR}/${LINEAGE_ID}" ]]; then
  echo "ERROR: lineage directory not found: ${LINEAGES_DIR}/${LINEAGE_ID}" >&2
  exit 1
fi

if [[ ! -f "$GET_CDS_PY" ]]; then
  echo "ERROR: Could not find CDS extraction script: $GET_CDS_PY" >&2
  echo "       Set GET_CDS_PY=/full/path/to/02a_get_cds_from_compleasm_v6.py" >&2
  exit 1
fi

if ! command -v compleasm >/dev/null 2>&1 && ! command -v compleasm.py >/dev/null 2>&1; then
  echo "ERROR: compleasm not found in PATH. Activate your latest compleasm env." >&2
  exit 1
fi

# Prefer the installed entrypoint
if command -v compleasm >/dev/null 2>&1; then
  COMPLEASM="compleasm"
else
  COMPLEASM="compleasm.py"
fi

# Optional warning about env
if [[ -n "${CONDA_DEFAULT_ENV:-}" && -n "$EXPECTED_ENV" && "${CONDA_DEFAULT_ENV}" != "$EXPECTED_ENV" ]]; then
  echo "[WARN] Active conda env: ${CONDA_DEFAULT_ENV} (expected: ${EXPECTED_ENV})" >&2
  echo "[WARN] Not fatal, but may explain missing tools/odd behavior." >&2
fi

# Initialize metadata output (fresh header if missing)
if [[ ! -f "$METADATA_OUT" ]]; then
  echo "genus_species,accession,organism_name,lineage,full_table,cds_fasta" > "$METADATA_OUT"
else
  echo "[Info] Appending to existing: $METADATA_OUT"
fi

append_metadata_row() {
  local genus_species="$1"
  local accession="$2"
  local organism="$3"
  local lineage="$4"
  local full_table="$5"
  local cds_fasta="$6"

  # Prevent exact duplicate rows
  if grep -Fq "\"${genus_species}\",\"${accession}\",\"${organism}\",\"${lineage}\",\"${full_table}\",\"${cds_fasta}\"" "$METADATA_OUT" 2>/dev/null; then
    return
  fi

  printf '"%s","%s","%s","%s","%s","%s"\n' \
    "$genus_species" "$accession" "$organism" "$lineage" "$full_table" "$cds_fasta" >> "$METADATA_OUT"
}

archive_lineage_dir() {
  local outdir="$1"
  local lineage="$2"

  local ts archive_root
  ts="$(date +%Y-%m-%d_%H-%M-%S)"
  archive_root="${outdir}/archive"
  mkdir -p "$archive_root"

  local d="${outdir}/${lineage}"
  if [[ -d "$d" ]]; then
    echo "[ARCHIVE] Moving ${d} -> ${archive_root}/${lineage}__archived_${ts}"
    mv "$d" "${archive_root}/${lineage}__archived_${ts}"
  fi
}

# Also archive any legacy odb10/odb12 lineage dirs if you force re-run
archive_any_lineage_dirs() {
  local outdir="$1"
  local ts archive_root
  ts="$(date +%Y-%m-%d_%H-%M-%S)"
  archive_root="${outdir}/archive"
  mkdir -p "$archive_root"

  local existing
  existing="$(find "$outdir" -maxdepth 1 -type d \( -name "*_odb10" -o -name "*_odb12" \) 2>/dev/null || true)"
  if [[ -n "$existing" ]]; then
    while IFS= read -r d; do
      [[ -z "$d" ]] && continue
      local base
      base="$(basename "$d")"
      echo "[ARCHIVE] Moving existing ${base} -> ${archive_root}/${base}__archived_${ts}"
      mv "$d" "${archive_root}/${base}__archived_${ts}"
    done <<< "$existing"
  fi
}

find_full_table() {
  local lin_dir="$1"
  if [[ -f "${lin_dir}/full_table.tsv" ]]; then
    echo "${lin_dir}/full_table.tsv"
    return
  fi
  if [[ -f "${lin_dir}/full_table.csv" ]]; then
    echo "${lin_dir}/full_table.csv"
    return
  fi
  local any
  any="$(ls -1 "${lin_dir}"/full_table.* 2>/dev/null | head -n 1 || true)"
  echo "$any"
}

run_one() {
  local accession="$1"
  local organism="$2"
  local genome="$3"

  if [[ -z "$genome" || ! -f "$genome" ]]; then
    echo "[WARN] Genome missing, skipping: $genome"
    return
  fi

  local safe_org outdir lin_dir
  safe_org="$(echo "$organism" | tr ' /' '__' | tr -cd '[:alnum:]_-.')"
  outdir="${OUT_ROOT}/${accession}__${safe_org}"
  lin_dir="${outdir}/${LINEAGE_ID}"

  # Derive genus_species from organism name (first two words)
  local genus_species
  genus_species="$(echo "$organism" | awk '{print $1"_"$2}' | tr '[:upper:]' '[:lower:]' | tr -cd '[:alnum:]_')"
  if [[ -z "$genus_species" ]]; then
    echo "ERROR: Failed to derive genus_species from organism: $organism" >&2
    exit 1
  fi

  mkdir -p "$outdir"

  local gff="${lin_dir}/miniprot_output.gff"
  local full_table
  full_table="$(find_full_table "$lin_dir")"
  local cds_fasta="${lin_dir}/${genus_species}_cds_compleasm.fasta"

  # ---- Skip logic ----
  if [[ "$FORCE_RERUN" -ne 1 ]]; then
    if [[ -s "$gff" && -n "$full_table" && -s "$full_table" && -s "$cds_fasta" ]]; then
      echo "[SKIP] Found existing outputs for ${accession} (${LINEAGE_ID}); not re-running."
      append_metadata_row "$genus_species" "$accession" "$organism" "$LINEAGE_ID" "$full_table" "$cds_fasta"
      return
    fi

    # If lineage dir exists but incomplete, archive it and rebuild cleanly
    if [[ -d "$lin_dir" ]]; then
      echo "[INFO] Incomplete outputs detected for ${accession} (${LINEAGE_ID}); archiving and regenerating."
      archive_lineage_dir "$outdir" "$LINEAGE_ID"
    fi
  else
    # Force: archive any existing lineage dirs (odb10/odb12) to keep history
    echo "[INFO] Force rerun requested (-f). Archiving existing lineage dirs under: $outdir"
    archive_any_lineage_dirs "$outdir"
  fi

  echo "[RUN ] $accession | $organism | lineage=${LINEAGE_ID}"

  # ---- Compleasm run (internal miniprot) ----
  # NOTE: `compleasm run` expects:
  #   -a / --assembly for genome fasta
  #   -l lineage name
  #   -o output dir
  #   -t threads
  #   -L lineage library path (your local mb_downloads folder)
  #
  # If hmmsearch isn't on PATH in this env, set HMMSEARCH_BIN=/full/path/to/hmmsearch
  set -x
  if [[ -n "$HMMSEARCH_BIN" ]]; then
    "$COMPLEASM" run -a "$genome" -l "$LINEAGE_ID" -o "$outdir" -t "$THREADS" -L "$LINEAGES_DIR" --hmmsearch_execute_path "$HMMSEARCH_BIN"
  else
    "$COMPLEASM" run -a "$genome" -l "$LINEAGE_ID" -o "$outdir" -t "$THREADS" -L "$LINEAGES_DIR"
  fi
  set +x

  # ---- Validate expected outputs ----
  if [[ ! -s "$gff" ]]; then
    echo "ERROR: Expected GFF not created or empty: $gff" >&2
    echo "       This suggests compleasm's internal miniprot step failed or produced no stdout." >&2
    exit 1
  fi

  full_table="$(find_full_table "$lin_dir")"
  if [[ -z "$full_table" || ! -s "$full_table" ]]; then
    echo "ERROR: Could not find non-empty full_table in: $lin_dir" >&2
    exit 1
  fi

  # ---- Extract CDS (your existing Python) ----
  echo "[CDS ] python3 $(basename "$GET_CDS_PY") <genome> <full_table> ${genus_species}"
  python3 "$GET_CDS_PY" "$genome" "$full_table" "$genus_species"

  if [[ ! -s "$cds_fasta" ]]; then
    echo "ERROR: CDS FASTA not created or empty: $cds_fasta" >&2
    exit 1
  fi

  append_metadata_row "$genus_species" "$accession" "$organism" "$LINEAGE_ID" "$full_table" "$cds_fasta"
}

# -----------------------
# Batch mode
# -----------------------
if [[ "$QUERY" == "allgenomes" ]]; then
  echo "[INFO] Running compleasm on ALL genomes with lineage=${LINEAGE_ID}"
  awk -F',' '
    function trim(s){ gsub(/^[ \t"]+|[ \t"]+$/, "", s); return s }
    NR==1{ for(i=1;i<=NF;i++) h[trim($i)]=i; next }
    {
      acc=trim($(h["accession"]))
      org=trim($(h["organism_name"]))
      fna=trim($(h["path_to_fna"]))
      if(acc!="" && fna!="") print acc "\t" org "\t" fna
    }
  ' "$METADATA_CSV" |
  while IFS=$'\t' read -r acc org fna; do
    run_one "$acc" "$org" "$fna"
  done
  echo "[INFO] allgenomes run complete"
  exit 0
fi

# -----------------------
# Single-genome mode (robust matching)
# -----------------------
match="$(
awk -F',' -v q="$QUERY" '
  function trim(s){ gsub(/^[ \t"]+|[ \t"]+$/, "", s); return s }
  function norm(s){ s=tolower(trim(s)); gsub(/[ \t_-]+/, "_", s); gsub(/[^a-z0-9_.]/, "", s); return s }
  BEGIN{ qn=norm(q) }
  NR==1{ for(i=1;i<=NF;i++) h[trim($i)]=i; next }
  {
    acc=trim($(h["accession"]))
    org=trim($(h["organism_name"]))
    fna=trim($(h["path_to_fna"]))
    if(acc==q || fna==q || norm(org)==qn)
      print acc "\t" org "\t" fna
  }
' "$METADATA_CSV"
)"

n=$(echo "$match" | awk 'NF{c++} END{print c+0}')
if [[ "$n" -ne 1 ]]; then
  echo "ERROR: Query matched $n rows; use accession for uniqueness." >&2
  echo "$match" >&2
  exit 1
fi

IFS=$'\t' read -r acc org fna <<< "$match"
run_one "$acc" "$org" "$fna"

echo "[INFO] Done"
