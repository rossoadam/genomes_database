#!/bin/bash

set -e

###### USAGE ######
# bash GENOMES_DIR PHY_FILE
# bash ~/projects/genomes/ /Users/rossoaa/projects/genomes/records/compleasm/alignments/05_unrooted_tree_t100_e1_o7/rev_dna.treefile

# Check inputs
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <genomes_dir> <phylip_file> <outgroup_taxa>"
    exit 1
fi

GENOMES_DIR="$1"
PHY_FILE="$2"
OUT_TAXA="$3"

# Define output directory
OUT_DIR="$GENOMES_DIR/records/compleasm/alignments/06_rooted_treefile"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Move into output directory so IQ-TREE writes files here
cd "$OUT_DIR"

# Run IQ-TREE
iqtree3 -s "$PHY_FILE" \
        -o "$OUT_TAXA"

echo "OUTGROUP ESTABILISHED with $OUT_TAXA"
echo "Output directory: $OUT_DIR"