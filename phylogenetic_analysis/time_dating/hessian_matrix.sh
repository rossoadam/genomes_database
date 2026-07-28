#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Help menu block
show_help() {
    echo "Usage: bash run_iqtree_hessian.sh <alignment.phy> <model> <rooted_tree.nwk>"
    echo ""
    echo "Arguments:"
    echo "  1. alignment.phy  : Path to the multiple sequence alignment file"
    echo "  2. model          : Substitution model (e.g., GTR+G4 or 'MIX{GTR,HKY}+G4')"
    echo "  3. rooted_tree.nwk: Path to the rooted Newick tree topology file"
    echo ""
    echo "Example:"
    echo "  bash run_iqtree_hessian.sh my_data.phy GTR+G4 my_rooted_tree.nwk"
    exit 1
}

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Error: Invalid number of arguments."
    show_help
fi

# Assign terminal positional parameters to readable variables
ALIGNMENT=$1
MODEL=$2
TREE=$3

# Check if input files actually exist before running the compute
if [ ! -f "$ALIGNMENT" ]; then
    echo "Error: Alignment file '$ALIGNMENT' not found."
    exit 1
fi

if [ ! -f "$TREE" ]; then
    echo "Error: Tree file '$TREE' not found."
    exit 1
fi

echo "========================================================"
echo "Starting IQ-TREE 3 Hessian Matrix Generation"
echo "Alignment:  $ALIGNMENT"
echo "Model:      $MODEL"
echo "Tree:       $TREE"
echo "========================================================"

# Run IQ-TREE 3 with the terminal arguments
iqtree3 -s "$ALIGNMENT" \
        -m "$MODEL" \
        -te "$TREE" \
        --dating mcmctree \
        --mcmc-clock IND \
        --mcmc-iter 20000,200,50000 \
        --mcmc-birth-death-sampling 1,1,0.5

echo "========================================================"
echo "IQ-TREE 3 execution finished successfully."
echo "Generated outputs: mcmctree.hessian, rooted.nwk, mcmctree.ctl"
echo "========================================================"
