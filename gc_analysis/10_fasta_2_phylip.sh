#!/bin/bash

# use conda environment samtools for this to access seqmagick cli library
# usage
#
# bash 11_fasta_2_phylip.sh /Users/rossoaa/projects/genomes/records/compleasm/alignments/07_renamed_fastas
#

set -e

# Check input argument
if [ -z "$1" ]; then
    echo "Usage: $0 <input_dir> [output_dir]"
    exit 1
fi

INPUT_DIR="$1"

# Default output directory if not provided
if [ -z "$2" ]; then
    PARENT_DIR=$(dirname "${INPUT_DIR%/}")
	OUTPUT_DIR="$PARENT_DIR/08_fasta_converted_2_phylip"
else
    OUTPUT_DIR="$2"
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"

# Loop through fasta files
for input_fasta in "$INPUT_DIR"/*.fasta; do

    # Skip if no files match
    [ -e "$input_fasta" ] || continue

    # Output filename
    base_name=$(basename "$input_fasta" .fasta)
    output_phylip="$OUTPUT_DIR/${base_name}.phylip"

    # Convert
    seqmagick convert --output-format phylip-relaxed "$input_fasta" "$output_phylip"

    if [ $? -eq 0 ]; then
        echo "Converted: $base_name.fasta → $output_phylip"
    else
        echo "Failed: $base_name.fasta"
    fi

done

echo "Done."