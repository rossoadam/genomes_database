#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") <query|allgenomes> [-t threads] [-l lineage] [-f] [-h]

Positional:
  query        Accession OR organism name (case-insensitive; spaces/underscores treated the same)
  allgenomes   Run compleasm for every row in METADATA_CSV

Options:
  -t INT       Threads (default: ${THREADS:-8})
  -l STR       Lineage ID (default: ${LINEAGE_ID:-rickettsiales_odb12})
  -f           Force rerun (archive lineage dir even if outputs exist)
  -h           Help

Env overrides:
  METADATA_CSV, OUT_ROOT, LIBDIR, METADATA_OUT, GET_CDS_PY, THREADS, LINEAGE_ID
EOF
}

# If user asks for help, show it and exit before any other logic
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

QUERY="${1:-}"
if [[ -z "$QUERY" ]]; then
  echo "ERROR: missing required argument <query|allgenomes>" >&2
  usage >&2
  exit 2
fi
shift || true

# Defaults (can still be overridden by env or flags)
####### THREADS #######
THREADS="${THREADS:-8}"
####### LINEAGE #######
LINEAGE_ID="${LINEAGE_ID:-sauropsida_odb12}"
#######
FORCE_RERUN="${FORCE_RERUN:-0}"

# Parse flags
while getopts ":t:l:fh" opt; do
  case "$opt" in
    t) THREADS="$OPTARG" ;;
    l) LINEAGE_ID="$OPTARG" ;;
    f) FORCE_RERUN=1 ;;
    h) usage; exit 0 ;;
    \?) echo "ERROR: invalid option -$OPTARG" >&2; usage >&2; exit 2 ;;
    :)  echo "ERROR: option -$OPTARG requires an argument" >&2; usage >&2; exit 2 ;;
  esac
done
shift $((OPTIND - 1))

# -----------------
########## genomes metadata ##########
METADATA_CSV="${METADATA_CSV:-/Users/rossoaa/projects/genomes/records/genomes_metadata.csv}"

########## singular genome path ##########
# GENOME=""

########## out directory ##############
OUT_ROOT="/Users/rossoaa/projects/genomes/records/compleasm" 

########## directory w/ protein dbs ###########
LIBDIR="/Users/rossoaa/projects/genomes/records/compleasm/mb_downloads"

########## compleasm out metadata ##########
METADATA_OUT="${METADATA_OUT:-${OUT_ROOT}/metadata.csv}"

########## CDS extraction script (update path if needed)
GET_CDS_PY="${GET_CDS_PY:-$(dirname "$0")/02a_get_cds_from_compleasm_v6.py}"

# -----------------
# Metadata output (step 4)
# -----------------
# v13 metadata headers:
# genus_species,accession,organism_name,lineage,full_table,cds_fasta

# Initialize metadata output (fresh header if missing). Note plain csv no quoting fields
if [[ ! -f "$METADATA_OUT" ]]; then
  mkdir -p "$(dirname "$METADATA_OUT")"
  echo "genus_species,accession,organism_name,lineage,full_table,cds_fasta" > "$METADATA_OUT"
else
  echo "[INFO] Appending to existing: $METADATA_OUT" >&2
fi

append_metadata_row() {
  local genus_species="$1"
  local accession="$2"
  local organism_name="$3"
  local lineage="$4"
  local full_table="$5"
  local cds_fasta="$6"

# prevent exact duplicate rows
  local line="${genus_species},${accession},${organism_name},${lineage},${full_table},${cds_fasta}"
  if grep -Fqx "$line" "$METADATA_OUT" 2>/dev/null; then
    return 0
  fi

  # Remove any existing row with same accession + lineage (keeps metadata “one row per run target”)
  # (We keep header line intact)
  tmp="${METADATA_OUT}.tmp.$$"
  awk -F',' -v acc="$accession" -v lin="$lineage" '
    NR==1 {print; next}
    { if($2==acc && $4==lin) next; print }
  ' "$METADATA_OUT" > "$tmp" && mv "$tmp" "$METADATA_OUT"

  printf '%s,%s,%s,%s,%s,%s\n' "$genus_species" "$accession" "$organism_name" "$lineage" "$full_table" "$cds_fasta" >> "$METADATA_OUT"

}

find_full_table() {
  local lin_dir="$1"

  if [[ -f "${lin_dir}/full_table.csv" ]]; then
    echo "${lin_dir}/full_table.csv"
    return 0
  fi
  if [[ -f "${lin_dir}/full_table.tsv" ]]; then
    echo "${lin_dir}/full_table.tsv"
    return 0
  fi

  # Fallback: first match the full_table.*
  local any
  any="$(ls -1 "${lin_dir}"/full_table.* 2>/dev/null | head -n 1 || true)"
  echo "$any"
}

archive_lineage_dir() {
  local outdir="$1"
  local lineage="$2"

  local ts archive_root d
  ts="$(date +%Y-%m-%d_%H-%M-%S)"
  archive_root="${outdir}/archive"
  mkdir -p "$archive_root"

  d="${outdir}/${lineage}"
  if [[ -d "$d" ]]; then
    echo "[ARCHIVE] Moving ${d} -> ${archive_root}/${lineage}__archived_${ts}" >&2
    mv "$d" "${archive_root}/${lineage}__archived_${ts}"
  fi
}

run_one() {
  local accession="$1"
  local organism_name="$2"
  local genome="$3"

  if [[ -z "$genome" || ! -f "$genome" ]]; then
    echo "[WARN] Genome missing; skipping: $genome" >&2
    return 0
  fi

  # genus_species (use underscore name directly if present; else first two words)
  local genus_species
  if [[ "$organism_name" == *_* ]]; then
    genus_species="$(echo "$organism_name" | tr '[:upper:]' '[:lower:]' | tr -cd '[:alnum:]_')"
  else
    genus_species="$(echo "$organism_name" | awk '{print $1"_"$2}' | tr '[:upper:]' '[:lower:]' | tr -cd '[:alnum:]_')"
  fi
  genus_species="${genus_species%_}"
  if [[ -z "$genus_species" ]]; then
    echo "ERROR: Could not derive genus_species from: $organism_name" >&2
    return 1
  fi

  local safe_org
  safe_org="$(echo "$organism_name" | tr ' /' '__' | tr -cd '[:alnum:]_-.')"

  local outdir lin_dir gff full_table cds_fasta
  outdir="${OUT_ROOT}/${accession}__${safe_org}"
  lin_dir="${outdir}/${LINEAGE_ID}"
  
  gff="${lin_dir}/miniprot_output.gff"
  full_table="$(find_full_table "$lin_dir")"
  cds_fasta="${lin_dir}/${genus_species}_cds_compleasm.fasta"

  echo "[INFO] accession=$accession" >&2
  echo "[INFO] genus_species=$genus_species" >&2
  echo "[INFO] OUTDIR=$outdir" >&2
  echo "[INFO] LIN_DIR=$lin_dir" >&2

  mkdir -p "$outdir"

  # Skip logic unless forcing
  if [[ "$FORCE_RERUN" -ne 1 ]]; then
    if [[ -s "$gff" && -n "$full_table" && -s "$full_table" && -s "$cds_fasta" ]]; then
      echo "[SKIP] Found existing outputs for ${accession} (${LINEAGE_ID})" >&2
      append_metadata_row "$genus_species" "$accession" "$organism_name" "$LINEAGE_ID" "$full_table" "$cds_fasta"
      return 0
    fi

    # If lineage dir exists but incomplete, archive it
    if [[ -d "$lin_dir" ]]; then
      echo "[INFO] Incomplete lineage outputs detected; archiving and re-running $lin_dir" >&2
      archive_lineage_dir "$outdir" "$LINEAGE_ID"
    fi
  else
    # Force: archive lineage dir if present
    if [[ -d "$lin_dir" ]]; then
      echo "[INFO] Force rerun (-f); archiving existing lineage dir: $lin_dir" >&2
      archive_lineage_dir "$outdir" "$LINEAGE_ID"
    fi
  fi

  # Run compleasm
  compleasm run --assembly_path "$genome" --output_dir "$outdir" --threads "$THREADS" --lineage "$LINEAGE_ID" --library "$LIBDIR"

  # Validate + locate full_table
  full_table="$(find_full_table "$lin_dir")"
  echo "[INFO] FULL_TABLE=$full_table" >&2

  if [[ ! -s "$gff" ]]; then
    echo "ERROR: Expected GFF missing/empty: $gff" >&2
    return 1
  fi
  if [[ -z "$full_table" || ! -s "$full_table" ]]; then
    echo "ERROR: Expected full_table missing/empty under: $lin_dir" >&2
    return 1
  fi

  # CDS extraction
  echo "[CDS ] python3 $(basename "$GET_CDS_PY") <genome> <full_table> <genus_species>" >&2
  python3 "$GET_CDS_PY" "$genome" "$full_table" "$genus_species"

  if [[ ! -s "$cds_fasta" ]]; then
    echo "ERROR: CDS FASTA missing/empty: $cds_fasta" >&2
    return 1
  fi

  append_metadata_row "$genus_species" "$accession" "$organism_name" "$LINEAGE_ID" "$full_table" "$cds_fasta"
  echo "[INFO] Updated metadata: $METADATA_OUT" >&2
}

# -----------------
# Sanity checks (step 2)
# -----------------

# Prefer installed entrypoint check
if ! command -v compleasm >/dev/null 2>&1; then
 echo "ERROR: compleasm not found in PATH. Activate the correct env." >&2
 exit 1
fi

# Library base dir must exist
if [[ ! -d "$LIBDIR" ]]; then
  echo "ERROR: LIBDIR not found: $LIBDIR" >&2
  exit 1
fi

# Lineage folder must exist inside library dir (my current v15 expects unpacked lineage folders here)
if [[ ! -d "${LIBDIR}/${LINEAGE_ID}" ]]; then
  echo "ERROR: lineage not found in LIBDIR: ${LIBDIR}/${LINEAGE_ID}" >&2
  echo "       Check LINEAGE_ID (-l) or LIBDIR." >&2
  exit 1
fi

# Output root should be wrtiable/creatable
mkdir -p "$OUT_ROOT"
if [[ ! -d "$OUT_ROOT" ]]; then
  echo "ERROR: could not create OUT_ROOT: $OUT_ROOT" >&2
  exit 1
fi

if [[ ! -f "$GET_CDS_PY" ]]; then
  echo "ERROR: CDS  extraction script not found: $GET_CDS_PY" >&2
  echo "       SET GET_CDS_PY=Set GET_CDS_PY=/full/path/to/02a_get_cds_from_compleasm_v6.py" >&2
  exit 1
fi

# -----------------
# Main dispatch (v20)
# -----------------

if [[ "$QUERY" == "allgenomes" ]]; then
  echo "[INFO] Running compleasm on ALL genomes with lineage=${LINEAGE_ID}" >&2
  awk -F',' '
    function trim(s){ gsub(/^[ \t"]+|[ \t"]+$/, "", s); return s}
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

  echo "[INFO] allgenomes run complete" >&2
  exit 0
fi

# single-genome mode: robust matching
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
      print( acc "\t" org "\t" fna )
   }

' "$METADATA_CSV"
)"

n=$(echo "$match" | awk 'NF{c++} END{print c+0}')
if [[ "$n" -ne 1 ]]; then
  echo "ERROR: query matched $n rows; use accession for uniqueness." >&2
  echo "$match" >&2
  exit 1
fi

IFS=$'\t' read -r acc org fna <<< "$match"
run_one "$acc" "$org" "$fna"
