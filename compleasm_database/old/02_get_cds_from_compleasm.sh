#!/usr/bin/env bash

source /Users/rossoaa/miniconda3/etc/profile.d/conda.sh
conda activate liftoff_2

OUTDIR_COMPLEASM="/Users/rossoaa/projects/genomes/records/compleasm/"
GENOME="/Users/rossoaa/projects/genomes/GCF_004329235.1/ncbi_dataset/data/GCF_004329235.1/GCF_004329235.1_PodMur_1.0_genomic.fna"
#FULL_TSV="/Users/rossoaa/projects/genomes/records/compleasm/GCF_004329235.1__Podarcis_muralis/sauropsida_odb12/full_table.tsv"
FULL_TSV="/Users/rossoaa/projects/genomes/records/compleasm/GCF_004329235.1__Podarcis_muralis/sauropsida_odb10/full_table.tsv"

python3 get_cds_from_compleasm_v6.py \
 "$GENOME" \
 "$FULL_TSV" \
 podarcis_muralis
