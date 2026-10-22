#!/bin/bash

set -euo pipefail

#
# Run nhPhyml on every .phylip file in a directory, optionally in parallel.
#
# Usage:
#   bash 13_nhphyml_parallel.sh \
#     /path/to/phylip_input_dir \
#     /path/to/rooted_tree.tre \
#     40 \
#     [/path/to/nhPhyml_binary] \
#     [/path/to/output_dir]
#
# Arguments:
#   1. phylip_input_dir   Required. Directory containing .phylip files.
#   2. tree_file          Required. Rooted tree file to pass to nhPhyml.
#   3. threads            Required. Number of parallel nhPhyml jobs to run.
#   4. nhphyml_binary     Optional. Path to nhPhyml binary.
#   5. output_dir         Optional. Directory where all nhPhyml outputs go.
#
# Notes:
# - This script assumes one nhPhyml process uses one core effectively.
# - Each job runs in its own temporary directory to avoid file collisions.
# - All outputs are then moved into one flat output directory.
#

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <phylip_input_dir> <tree_file> <threads> [nhphyml_binary] [output_dir]"
    exit 1
fi

INPUT_DIR="$1"
TREE_FILE="$2"
THREADS="$3"
NHPHYML_BIN="${4:-/home/lepidodactylus/bin/nhPhyml/nhPhyml}"

if [ -n "${5:-}" ]; then
    OUTPUT_DIR="$5"
else
    PARENT_DIR=$(dirname "${INPUT_DIR%/}")
    OUTPUT_DIR="$PARENT_DIR/09_nhphyml_output"
fi

# Validate inputs
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory not found: $INPUT_DIR"
    exit 1
fi

if [ ! -f "$TREE_FILE" ]; then
    echo "Error: Tree file not found: $TREE_FILE"
    exit 1
fi

if [ ! -x "$NHPHYML_BIN" ]; then
    echo "Error: nhPhyml binary not executable: $NHPHYML_BIN"
    exit 1
fi

if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: threads must be a positive integer. Got: $THREADS"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "Input directory : $INPUT_DIR"
echo "Tree file       : $TREE_FILE"
echo "Threads         : $THREADS"
echo "nhPhyml binary  : $NHPHYML_BIN"
echo "Output directory: $OUTPUT_DIR"
echo "-----------------------------------"

# Confirm there is work to do before launching xargs.
PHYLIP_COUNT=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "*.phylip" | wc -l)
PHYLIP_COUNT=$(echo "$PHYLIP_COUNT" | tr -d '[:space:]')

if [ "$PHYLIP_COUNT" -eq 0 ]; then
    echo "No .phylip files found in: $INPUT_DIR"
    exit 1
fi

echo "Found $PHYLIP_COUNT phylip files."
echo "Launching parallel nhPhyml jobs..."
echo "-----------------------------------"

run_one() {
    local phylip_file="$1"
    local base_name
    local job_tmp
    local status_file
    local tree_copy
    local seq_copy

    base_name=$(basename "$phylip_file" .phylip)
    job_tmp=$(mktemp -d "${TMPDIR:-/tmp}/nhphyml_${base_name}_XXXXXX")
    status_file="$OUTPUT_DIR/${base_name}.status"

    cleanup() {
        rm -rf "$job_tmp"
    }
    trap cleanup RETURN

    seq_copy="$job_tmp/${base_name}.phylip"
    tree_copy="$job_tmp/$(basename "$TREE_FILE")"

    cp "$phylip_file" "$seq_copy"
    cp "$TREE_FILE" "$tree_copy"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] START $base_name"

    if (
        cd "$job_tmp"
        "$NHPHYML_BIN" \
            -sequences="$seq_copy" \
            -tree="$tree_copy" \
            -format=i \
            -positions=3 \
            -tstv=e \
            -rates=1 \
            -alpha=e \
            -topology=k \
            -outseqs=n \
            -eqfreq=lim \
            -numeqfreq=5 \
            -treefile="${base_name}.treefile"
    ); then
        # Move every output into one shared output directory.
        # Keep the main treefile name clean; prefix ancillary files to avoid collisions.
        shopt -s nullglob
        for f in "$job_tmp"/*; do
            [ -f "$f" ] || continue
            local fname
            fname=$(basename "$f")

            case "$fname" in
                "${base_name}.phylip"|"$(basename "$TREE_FILE")")
                    continue
                    ;;
                "${base_name}.treefile")
                    mv "$f" "$OUTPUT_DIR/${base_name}.treefile"
                    ;;
                *)
                    mv "$f" "$OUTPUT_DIR/${base_name}__${fname}"
                    ;;
            esac
        done
        shopt -u nullglob

        echo "SUCCESS" > "$status_file"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE  $base_name"
    else
        echo "FAILED" > "$status_file"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] FAIL  $base_name" >&2
        return 1
    fi
}

export INPUT_DIR TREE_FILE OUTPUT_DIR NHPHYML_BIN
export -f run_one

find "$INPUT_DIR" -maxdepth 1 -type f -name "*.phylip" -print0 |
    xargs -0 -n 1 -P "$THREADS" bash -c 'run_one "$1"' _

SUCCESS_COUNT=$(grep -l '^SUCCESS$' "$OUTPUT_DIR"/*.status 2>/dev/null | wc -l | tr -d '[:space:]')
FAILED_COUNT=$(grep -l '^FAILED$' "$OUTPUT_DIR"/*.status 2>/dev/null | wc -l | tr -d '[:space:]')

echo "-----------------------------------"
echo "Completed nhPhyml runs."
echo "Successful: $SUCCESS_COUNT"
echo "Failed    : $FAILED_COUNT"
echo "Outputs   : $OUTPUT_DIR"

if [ "$FAILED_COUNT" -gt 0 ]; then
    exit 1
fi
