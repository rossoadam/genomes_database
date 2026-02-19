#!/usr/bin/env bash

set -euo pipefail

# paths
GENOME="/Users/rossoaa/Desktop/dummy/genomes/GCA_966213895.1/ncbi_dataset/data/GCA_966213895.1/GCA_966213895.1_Wolbachia_endosymbiont_of_Gonepteryx_cleopatra_genomic.fna"
OUTDIR="/Users/rossoaa/Desktop/dummy/genomes/records/compleasm"
THREADS=7
# use a short name - not path, for lineage; I think this needs to be exact match to what is found in..
# compleasm list --remote
#SLINEAGE="eudicots_odb12"
LINEAGE="rickettsiales_odb12"
# directory containing unpacked folder version
LIBDIR="/Users/rossoaa/Desktop/dummy/genomes/records/compleasm/mb_downloads"
#MINIPROT="/Users/rossoaa/miniprot/miniprot" now part of the conda environment
#hmmsearch is currently part of the conda environment

compleasm run --assembly_path "$GENOME" --output_dir "$OUTDIR" --threads "$THREADS" --lineage "$LINEAGE" --library "$LIBDIR"
