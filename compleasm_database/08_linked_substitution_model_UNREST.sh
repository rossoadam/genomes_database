#!/bin/bash

#
# usage:
# bash 08_linked_substitution_model_UNREST.sh GENOMES_DIR CONCAT_DIR SCHEME_DIR
# bash 08_linked_substitution_model_UNREST.sh /Users/rossoaa/projects/genomes /Users/rossoaa/projects/genomes/records/compleasm/alignments/04_concatenated_clean_alignments_t95_e1_o2 /Users/rossoaa/projects/genomes/records/compleasm/alignments/05_unrooted_tree/
#
set -e

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
	echo "Usage: $0 <genomes_dir> <concat_dir> <scheme_dir>"
	exit 1
fi

GENOMES_DIR="$1"
CONCAT_DIR="$2"
SCHEME_DIR="$3"

PHY_FILE="$CONCAT_DIR/concatenated.phy"
SCHEME_FILE="$SCHEME_DIR/rev_dna.best_scheme.nex"

OUT_DIR="$GENOMES_DIR/records/compleasm/alignments/06_unrooted_tree_nonreversible_unrestricted"
mkdir -p "$OUT_DIR"

if [ ! -f "$PHY_FILE" ]; then
    echo "Error: PHYLIP file not found: $PHY_FILE"
    exit 1
fi

if [ ! -f "$SCHEME_FILE" ]; then
    echo "Error: scheme file not found: $SCHEME_FILE"
    exit 1
fi

cd "$OUT_DIR"

iqtree3 -s "$PHY_FILE" \
        -p "$SCHEME_FILE" \
        --model-joint UNREST \
        -B 1000 \
        -T AUTO \
        --prefix nonrev_dna

echo "Rooted nonreversible tree inference complete."
echo "Output directory: $OUT_DIR"