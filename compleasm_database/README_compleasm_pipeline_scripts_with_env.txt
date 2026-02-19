README: Compleasm Pipeline Scripts (Updated with Conda Environment for ODB12)
===============================================================================

This README describes three core scripts used in the Compleasm-based genome
completeness and CDS extraction pipeline, and documents the Conda environment
used for running Compleasm with BUSCO odb12 lineages.

---------------------------------------------------------------------------
PIPELINE OVERVIEW
---------------------------------------------------------------------------

This pipeline performs the following high-level steps:

1. Select genomes from a metadata CSV file
2. Run Compleasm (miniprot + analyze) using BUSCO lineages (odb10 or odb12)
3. Generate completeness summaries
4. Extract BUSCO-complete CDS sequences for downstream analyses

The pipeline is designed to be reproducible, lineage-flexible, and compatible
with both legacy odb10 runs and newer odb12 runs.

---------------------------------------------------------------------------
SCRIPTS
---------------------------------------------------------------------------

1) run_compleasm_from_metadata_v6.sh
-----------------------------------

Purpose:
- Reads a metadata CSV containing genome accessions, organism names, and paths
  to genome FASTA files.
- Matches a query species name (or 'allgenomes') against metadata.
- Runs Compleasm in two explicit stages:
    a) miniprot alignment against the BUSCO reference proteome
    b) compleasm analyze on the resulting GFF
- Archives previous lineage runs to avoid overwriting results.

Key features:
- Case-insensitive species matching (e.g., podarcis_muralis, Podarcis muralis)
- Supports odb10 and odb12 lineages
- Explicit paths for miniprot and hmmsearch binaries
- Thread control via -t
- Safety checks for missing genomes or empty GFF output

Primary outputs:
- miniprot_output.gff
- compleasm summary statistics
- BUSCO gene classification tables

2) get_cds_from_compleasm_v6.py
-------------------------------

Purpose:
- Extracts CDS sequences corresponding to BUSCO-complete genes from Compleasm
  output.
- Uses the miniprot-generated GFF and the original genome FASTA.
- Outputs CDS FASTA files suitable for orthology, alignment, or GC analyses.

Key features:
- Uses pyfaidx for efficient FASTA indexing
- Handles strand orientation and multi-exon CDS
- Produces one CDS per BUSCO gene per genome

Inputs:
- Genome FASTA
- Compleasm full_table.tsv or equivalent summary
- miniprot_output.gff

Outputs:
- CDS FASTA files
- Optional per-gene logging

3) get_cds_from_compleasm.sh
----------------------------

Purpose:
- Lightweight wrapper script to run get_cds_from_compleasm_v6.py.
- Simplifies batch execution across genomes or lineages.

---------------------------------------------------------------------------
CONDA ENVIRONMENT (ODB12)
---------------------------------------------------------------------------

For odb12 analyses, the pipeline uses a dedicated Conda environment to ensure
tool compatibility and reproducibility.

Environment name:
- compleasm_env (or similar, as defined in compleasm_env.yml)

Key installed tools:
- compleasm == 0.2.6
- miniprot == 0.18
- hmmer == 3.4
- python == 3.10
- pandas, numpy, pyfaidx, matplotlib (dependencies)

Verified binary locations (within environment):
- compleasm:  <env>/bin/compleasm
- miniprot:   <env>/bin/miniprot
- hmmsearch:  <env>/bin/hmmsearch

This environment is compatible with:
- BUSCO odb12 lineage libraries
- Local lineage directories specified via -L
- Apple Silicon (osx-arm64) systems

To create the environment:
--------------------------
conda env create -f compleasm_env.yml
conda activate compleasm_env

---------------------------------------------------------------------------
LINEAGE LIBRARIES
---------------------------------------------------------------------------

BUSCO lineage directories (odb10 or odb12) must be unpacked locally and pointed
to explicitly using the -L flag in Compleasm.

Example:
- sauropsida_odb12/
- sauropsida_odb10/

Each directory should contain:
- refseq_db.faa.gz
- hmms/
- dataset.cfg
- scores_cutoff, lengths_cutoff

---------------------------------------------------------------------------
NOTES
---------------------------------------------------------------------------

- Direct calls to `compleasm analyze` were used to validate script behavior.
- Empty GFF issues were resolved by ensuring correct miniprot invocation and
  passing the correct GFF path via -g.
- This pipeline preserves backward compatibility with odb10 while enabling
  forward use of odb12 for updated analyses.

---------------------------------------------------------------------------
AUTHOR / MAINTAINER
---------------------------------------------------------------------------

Pipeline maintained by:
Adam Rosso
Quantitative Biology PhD
University of Texas at Arlington
