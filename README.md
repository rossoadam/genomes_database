2026SEP26
	Creating genomes_file_recovery_v3.py resulted in unintended behavior.
	The reference flag was dropped while querying "SQUAMATA".
	This resulted in the total number of genomes ballooning to 223.
	I had chat write a clean_genomes.py script for removing non-reference genomes with a set of rules.
	I am going to run the clean_genomes.py script now. 
	
	Cleanup plan
	------------
	Taxon: squamata
	Candidate accession directories: 52
	Estimated disk space recoverable: 51.417 GB
	
	Protected:
	  NCBI references: 159
	  Project manifests: 87
	  Local preferred assemblies: 101
	  Non-NCBI assemblies: 1
	
	52 directories will be MOVED TO TRASH.
	No files will be permanently deleted.
	
	Cleanup execution complete
	--------------------------
	Moved to trash: 52
	Move failures:  0
	Trash archive:  /Users/rossoaa/projects/genomes/trash/clean_genomes_20260926_131934
	Metadata backup: /Users/rossoaa/projects/genomes/trash/clean_genomes_20260926_131934/metadata_backup
	Audit report:   /Users/rossoaa/projects/genomes/records/cleanup/clean_genomes_20260926_131934.tsv
	
	Manifest reconciliation complete.
	  Accessions checked: 85
	  FASTAs downloaded: 0
	  FASTAs decompressed: 0
	  GFFs downloaded: 0
	  GFFs decompressed: 0
	  Repository GFF unavailable/not recovered: 50
	  Failed accessions: 1
    - GWHFILM00000000.1: 
    The gwh failure was expected 

2025SEP02_genomes_hpc.py will take any argument that ncbi_datasets CLI can understand and download genomes accordingly.
    to get this working you will need to create the respective conda environment genomes_management.yml

our shared directory
    /work/10950/axc5473/sharedirectory


Actually, I need to align the compleasm output to the busco output to make sure that these are correct
I want to next update the genomes database repository to include a compleasm argument
I imagine this working by importing script to the genomes_hpc.py script and then including an argument that runs compleasm