#!/bin/bash
cd ./01_pre_alignments
for d in */; do
 cd /media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/genomes_2/busco_run_05_results/01_pre_alignments/"$d"
 mv codon_align.fasta renamed_"${d[@]::-1}".fasta
done
