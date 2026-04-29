#!/bin/bash

set -e

#
# bash 13_nhphyml.sh \
#  /path/to/phylip_dir \
#  /path/to/rooted_tree.tre
#

# Check required arguments
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <phylip_input_dir> <tree_file> [nhphyml_binary] [output_dir]"
    exit 1
fi

INPUT_DIR="$1"
TREE_FILE="$2"

# Optional args
#NHPHYML_BIN="/Users/rossoaa/src/nhPhyml/nhPhyml"
NHPHYML_BIN="/home/lepidodactylus/bin/nhPhyml/nhPhyml"
# Default output directory (sister directory)
if [ -z "$4" ]; then
    PARENT_DIR=$(dirname "${INPUT_DIR%/}")
    OUTPUT_DIR="$PARENT_DIR/09_nhphyml_output"
else
    OUTPUT_DIR="$4"
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

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

echo "Input directory: $INPUT_DIR"
echo "Tree file: $TREE_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "nhPhyml binary: $NHPHYML_BIN"
echo "-----------------------------------"

# Loop through each .phylip file
for phylip_file in "$INPUT_DIR"/*.phylip; do

    [ -e "$phylip_file" ] || continue

    base_name=$(basename "$phylip_file" .phylip)

    output_tree="$OUTPUT_DIR/${base_name}.treefile"

    echo "Running nhPhyml on $base_name..."

    "$NHPHYML_BIN" \
        -sequences="$phylip_file" \
        -tree="$TREE_FILE" \
        -format=i \
        -positions=3 \
        -tstv=e \
        -rates=1 \
        -alpha=e \
        -topology=k \
        -outseqs=n \
        -eqfreq=lim \
        -numeqfreq=5 \
        -treefile="$output_tree"

    if [ $? -eq 0 ]; then
        echo "✔ Success: $base_name"
    else
        echo "✖ Failed: $base_name"
    fi

done

echo "All runs complete."
