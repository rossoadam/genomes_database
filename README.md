2025SEP02_genomes_hpc.py will take any argument that ncbi_datasets CLI can understand and download genomes accordingly.
    to get this working you will need to create the respective conda environment genomes_management.yml

our shared directory
    /work/10950/axc5473/sharedirectory


Actually, I need to align the compleasm output to the busco output to make sure that these are correct
I want to next update the genomes database repository to include a compleasm argument
I imagine this working by importing script to the genomes_hpc.py script and then including an argument that runs compleasm


What do you think is the best way to integrate these two 