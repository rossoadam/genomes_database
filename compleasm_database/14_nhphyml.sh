#!/bin/bash

# Directory containing the phylip files
input_dir="/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/genomes_2/busco_run_05_results/02_consolidate_copy/02_with_names/"

# Loop through each .phylip file in the input directory
for phylip_file in "$input_dir"*.phylip; do

  # Get the base name of the file (without directory or extension)
  base_name=$(basename "$phylip_file" .phylip)

  # Run nhPhyml command with the current phylip file and appropriate output treefile
  /home/lepidodactylus/bin/nhPhyml/nhPhyml \
  -sequences="$phylip_file" \
  -tree=nonrev_tree_rooted_with_geckos_copy.tre \
  -format=i \
  -positions=3 \
  -tstv=e \
  -rates=1 \
  -alpha=e \
  -topology=k \
  -outseqs=n \
  -eqfreq=lim \
  -numeqfreq=5 \
  -treefile="${base_name}_.treefile"

  # Check if the command was successful
  if [ $? -eq 0 ]; then
    echo "nhPhyml run successful for $base_name!"
  else
    echo "nhPhyml run failed for $base_name."
  fi

done

