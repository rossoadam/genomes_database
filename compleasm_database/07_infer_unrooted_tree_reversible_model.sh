#!/bin/bash

set -e

###### USAGE ######
# bash GENOMES_DIR PHY_FILE
# bash ~/projects/genomes/ /Users/rossoaa/projects/genomes/records/compleasm/alignments/04_concatenated_clean_alignments_t95_e1_o2/concatenated.phy

# Check inputs
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <genomes_dir> <phylip_file>"
    exit 1
fi

GENOMES_DIR="$1"
PHY_FILE="$2"

# Define output directory
OUT_DIR="$GENOMES_DIR/records/compleasm/alignments/05_unrooted_tree"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Move into output directory so IQ-TREE writes files here
cd "$OUT_DIR"

# Run IQ-TREE
iqtree3 -s "$PHY_FILE" \
        -p "${PHY_FILE}.nex" \
        -B 1000 \
        -T AUTO \
        --prefix rev_dna

echo "Tree inference complete:"
echo "Output directory: $OUT_DIR"