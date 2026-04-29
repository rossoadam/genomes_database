#!/bin/bash

set -e

######USAGE######
# bash genomes_dir/ /cleaned_alignment_dir/
# bash /Users/rosssoaa/projects/genomes/ /Users/rossoaa/projects/genomes/records/compleasm/alignments/03_clean_alignments_t95_e1_o2
#################

# Check input
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <genomes_dir> <03_align_dir>"
    exit 1
fi

GENOMES_DIR="$1"

# Define input and output paths
ALIGN_DIR="$2"
OUT_DIR="$GENOMES_DIR/records/compleasm/alignments/04_concatenated"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Output file
OUT_PHY="$OUT_DIR/concatenated.phy"

# Run IQ-TREE concatenation
iqtree3 -p "$ALIGN_DIR" --out-aln "$OUT_PHY"

echo "Concatenation complete:"
echo "Output: $OUT_PHY"
