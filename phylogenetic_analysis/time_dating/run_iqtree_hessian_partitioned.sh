#!/usr/bin/env bash

# Generate codon-position partitions from a concatenated codon alignment
# and calculate partition-specific IQ-TREE 3 Hessians for MCMCtree.
#
# Assumption:
#   Every locus begins at codon position 1 and contributes a number of
#   alignment columns divisible by 3. This must remain true after trimming
#   and concatenation.

set -euo pipefail

show_help() {
    cat <<'EOF'
Usage:
  bash run_iqtree_hessian_partitioned.sh \
      <alignment.phy> <model> <rooted_tree.nwk> [output_prefix]

Arguments:
  1. alignment.phy
       Concatenated nucleotide alignment in PHYLIP format.

  2. model
       IQ-TREE DNA substitution model applied separately to each
       codon-position partition, for example:
         GTR+G4
         GTR+F+G4
         MIX{GTR,HKY}+G4

       Quote models containing shell-special characters:
         'MIX{GTR,HKY}+G4'

  3. rooted_tree.nwk
       Fixed rooted Newick topology. Taxon names must match the alignment.

  4. output_prefix
       Optional IQ-TREE output prefix.
       Default: <alignment_without_extension>_codon_partitioned

Example:
  bash run_iqtree_hessian_partitioned.sh \
      concatenated.phy \
      'GTR+G4' \
      rooted_calibrated_tree.nwk \
      results/squamate_codon_partitioned
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    show_help
    exit 0
fi

if [[ "$#" -lt 3 || "$#" -gt 4 ]]; then
    echo "Error: expected 3 or 4 arguments; received $#." >&2
    echo >&2
    show_help >&2
    exit 1
fi

ALIGNMENT=$1
MODEL=$2
TREE=$3
PREFIX=${4:-"${ALIGNMENT%.*}_codon_partitioned"}

if [[ ! -f "$ALIGNMENT" ]]; then
    echo "Error: alignment file not found: $ALIGNMENT" >&2
    exit 1
fi

if [[ ! -f "$TREE" ]]; then
    echo "Error: tree file not found: $TREE" >&2
    exit 1
fi

if ! command -v iqtree3 >/dev/null 2>&1; then
    echo "Error: iqtree3 was not found in PATH." >&2
    exit 1
fi

# Read the first non-empty PHYLIP header line: <number_taxa> <number_sites>.
HEADER=$(awk 'NF >= 2 {print $1, $2; exit}' "$ALIGNMENT")
read -r NTAX NCHAR <<< "$HEADER"

# Remove a possible Windows carriage return.
NTAX=${NTAX//$'\r'/}
NCHAR=${NCHAR//$'\r'/}

if [[ ! "$NTAX" =~ ^[0-9]+$ || ! "$NCHAR" =~ ^[0-9]+$ ]]; then
    echo "Error: could not read numeric taxon/site counts from the PHYLIP header." >&2
    echo "Header interpreted as: '$HEADER'" >&2
    exit 1
fi

if (( NCHAR == 0 )); then
    echo "Error: alignment contains zero sites." >&2
    exit 1
fi

if (( NCHAR % 3 != 0 )); then
    echo "Error: alignment length ($NCHAR) is not divisible by 3." >&2
    echo "A global codon-position partition would therefore be invalid." >&2
    echo "Check trimming, concatenation, and locus boundaries." >&2
    exit 1
fi

PREFIX_DIR=$(dirname "$PREFIX")
mkdir -p "$PREFIX_DIR"

PARTITION_FILE="${PREFIX}.codon_positions.nex"

# Use the complete alignment length rather than relying on an open-ended
# partition range. IQ-TREE uses backslash notation for every third site.
printf '#nexus\n' > "$PARTITION_FILE"
printf 'begin sets;\n' >> "$PARTITION_FILE"
printf '    charset codon_pos1 = 1-%s\\3;\n' "$NCHAR" >> "$PARTITION_FILE"
printf '    charset codon_pos2 = 2-%s\\3;\n' "$NCHAR" >> "$PARTITION_FILE"
printf '    charset codon_pos3 = 3-%s\\3;\n' "$NCHAR" >> "$PARTITION_FILE"
printf 'end;\n' >> "$PARTITION_FILE"

echo "========================================================"
echo "Starting partitioned IQ-TREE 3 Hessian generation"
echo "Alignment:       $ALIGNMENT"
echo "Taxa:            $NTAX"
echo "Sites:           $NCHAR"
echo "Model:           $MODEL"
echo "Rooted tree:     $TREE"
echo "Partition file:  $PARTITION_FILE"
echo "Output prefix:   $PREFIX"
echo "Partition mode:  edge-unlinked (-Q)"
echo "========================================================"

cat "$PARTITION_FILE"

# -Q:
#   Treats codon positions as edge-unlinked partitions and writes separate
#   gradient/Hessian components for the three partitions.
#
# -te:
#   Fixes the supplied rooted topology while optimizing branch lengths and
#   model parameters.
#
# --mcmc-iter:
#   burnin, sample frequency, number of retained samples written into the
#   generated MCMCtree control file.
iqtree3 \
    -s "$ALIGNMENT" \
    -Q "$PARTITION_FILE" \
    -m "$MODEL" \
    -te "$TREE" \
    --dating mcmctree \
    --mcmc-clock IND \
    --mcmc-iter 20000,200,50000 \
    --mcmc-bds 1,1,0.5 \
    --prefix "$PREFIX"

echo "========================================================"
echo "IQ-TREE 3 finished successfully."
echo
echo "Important outputs:"
echo "  ${PREFIX}.mcmctree.hessian"
echo "  ${PREFIX}.rooted.nwk"
echo "  ${PREFIX}.mcmctree.ctl"
echo "  ${PREFIX}.dummy.aln"
echo
echo "Use the IQ-TREE-generated rooted tree with its matching Hessian."
echo "Check that the generated control file reports ndata = 3."
echo "========================================================"
