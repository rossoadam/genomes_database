#!/bin/bash

set -e

###### USAGE ######
# bash SCRIPT.sh GENOMES_DIR PHY_FILE [OUTGROUP_TXT]
# bash ~/projects/genomes/ /Users/rossoaa/projects/genomes/records/compleasm/alignments/04_concatenated_clean_alignments_t95_e1_o2/concatenated.phy
# bash ~/projects/genomes/ /path/to/concatenated.phy /path/to/outgroup.txt

# Check inputs
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <genomes_dir> <phylip_file> [outgroup_txt]"
    exit 1
fi

GENOMES_DIR="$1"
PHY_FILE="$2"
OUTGROUP_TXT="${3:-}"

# Define output directory
OUT_DIR="$GENOMES_DIR/records/compleasm/alignments/05_unrooted_tree"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Move into output directory so IQ-TREE writes files here
cd "$OUT_DIR"

# Optional outgroup handling:
# outgroup.txt should contain one taxon name per line.
IQTREE_OUTGROUP_ARGS=()
if [ -n "$OUTGROUP_TXT" ]; then
    if [ ! -f "$OUTGROUP_TXT" ]; then
        echo "Error: outgroup file not found: $OUTGROUP_TXT"
        exit 1
    fi

    OUTGROUPS=$(paste -sd, "$OUTGROUP_TXT")

    if [ -z "$OUTGROUPS" ]; then
        echo "Error: outgroup file is empty: $OUTGROUP_TXT"
        exit 1
    fi

    IQTREE_OUTGROUP_ARGS=(-o "$OUTGROUPS")
    echo "Using outgroup taxa: $OUTGROUPS"
fi

# Run IQ-TREE
iqtree3 -s "$PHY_FILE" \
        -p "${PHY_FILE}.nex" \
        -B 1000 \
        -T AUTO \
        --prefix rev_dna \
        "${IQTREE_OUTGROUP_ARGS[@]}"

echo "Tree inference complete:"
echo "Output directory: $OUT_DIR"
