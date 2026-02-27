#!/bin/bash

# Directory paths
input_dir="/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/genomes_2/busco_run_05_results/02_consolidate_copy/02_with_names/"
output_dir="/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/genomes_2/busco_run_05_results/02_consolidate_copy/02_with_names/"

# Loop through each fasta file in the input directory
for input_fasta in "$input_dir"*.fasta; do

  # Define output filename by replacing .fasta with .phylip
  output_phylip="${output_dir}$(basename "$input_fasta" .fasta).phylip"

  # Convert fasta to relaxed phylip format using seqmagick
  seqmagick convert --output-format phylip-relaxed "$input_fasta" "$output_phylip"

  # Check if conversion was successful
  if [ $? -eq 0 ]; then
    echo "Conversion to relaxed PHYLIP successful for $(basename "$input_fasta")! Output saved to $output_phylip"
  else
    echo "Conversion failed for $(basename "$input_fasta")."
  fi

done
