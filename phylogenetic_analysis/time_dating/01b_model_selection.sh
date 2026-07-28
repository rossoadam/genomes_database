#!/usr/bin/env bash

# Run IQ-TREE 3 ModelFinder on a fixed rooted topology using
# edge-unlinked codon-position partitions.

set -euo pipefail

show_help() {
    cat <<'EOF'
Usage:
  bash run_iqtree_codon_model_finder.sh \
      <alignment.phy> \
      <codon_partitions> \
      <rooted_tree.nwk> \
      [output_prefix] \
      [threads]

Arguments:
  1. alignment.phy
       Concatenated nucleotide alignment in PHYLIP format.

  2. codon_partitions
       Partition file defining codon positions 1, 2, and 3.

  3. rooted_tree.nwk
       Fixed rooted Newick topology. Taxon names must match the alignment.

  4. output_prefix
       Optional IQ-TREE output prefix.
       Default: codon_position_models

  5. threads
       Optional number of CPU threads, or AUTO.
       Default: AUTO

Example:
  bash run_iqtree_codon_model_finder.sh \
      concatenated.phy \
      concatenated.codon_positions.partitions \
      rooted_tree.nwk \
      codon_position_models \
      AUTO
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    show_help
    exit 0
fi

if [[ "$#" -lt 3 || "$#" -gt 5 ]]; then
    echo "Error: expected 3 to 5 arguments; received $#." >&2
    echo >&2
    show_help >&2
    exit 1
fi

ALIGNMENT=$1
PARTITIONS=$2
TREE=$3
PREFIX=${4:-codon_position_models}
THREADS=${5:-AUTO}

if [[ ! -f "$ALIGNMENT" ]]; then
    echo "Error: alignment file not found: $ALIGNMENT" >&2
    exit 1
fi

if [[ ! -f "$PARTITIONS" ]]; then
    echo "Error: partition file not found: $PARTITIONS" >&2
    exit 1
fi

if [[ ! -f "$TREE" ]]; then
    echo "Error: rooted tree file not found: $TREE" >&2
    exit 1
fi

# Locate IQ-TREE automatically unless the user specifies IQTREE_BIN

if [[ -z "${IQTREE_BIN:-}" ]]; then
    if command -v iqtree3 >/dev/null 2>&1; then
        IQTREE_BIN=$(command -v iqtree3)
    elif command -v iqtree3_intel >/dev/null 2>&1; then
        IQTREE_BIN=$(command -v iqtree3_intel)
    elif command -v iqtree >/dev/null 2>&1; then
        IQTREE_BIN=$(command -v iqtree)
    else
        echo "Error: Could not find IQ-TREE in PATH." >&2
        exit 1
    fi
fi

echo "Using IQ-TREE executable:"
echo "  $IQTREE_BIN"


PREFIX_DIR=$(dirname "$PREFIX")
mkdir -p "$PREFIX_DIR"

echo "========================================================"
echo "Starting IQ-TREE 3 codon-position ModelFinder analysis"
echo "Alignment:       $ALIGNMENT"
echo "Partitions:      $PARTITIONS"
echo "Fixed topology:  $TREE"
echo "Model search:    MFP+MERGE"
echo "Partition mode:  edge-unlinked (-Q)"
echo "Threads:         $THREADS"
echo "Output prefix:   $PREFIX"
echo "IQ-TREE binary:  $IQTREE_BIN"
echo "========================================================"

"$IQTREE_BIN" \
    -s "$ALIGNMENT" \
    -Q "$PARTITIONS" \
    -m MFP+MERGE \
    -te "$TREE" \
    -T "$THREADS" \
    --prefix "$PREFIX"

echo "========================================================"
echo "IQ-TREE 3 ModelFinder analysis finished successfully."
echo
echo "Inspect these outputs:"
echo "  ${PREFIX}.iqtree"
echo "  ${PREFIX}.best_model.nex"
echo "  ${PREFIX}.best_scheme"
echo "  ${PREFIX}.best_scheme.nex"
echo "  ${PREFIX}.log"
echo
echo "Before using the selected scheme for MCMCtree Hessian"
echo "generation, confirm that the final scheme still contains"
echo "the intended codon-position partitions."
echo "========================================================"
