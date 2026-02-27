#!/usr/bin/env python 
import csv

# Load metadata from CSV file
def load_metadata(csv_file):
    accession_to_organism = {}
    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            accession = row['accession']
            organism_name = row['organismName']
            accession_to_organism[accession] = organism_name
    return accession_to_organism

# Replace accession numbers with organism names in the PHYLIP file
def replace_taxa_in_phylip(phylip_file, accession_to_organism, output_file):
    with open(phylip_file, 'r') as file:
        lines = file.readlines()

    with open(output_file, 'w') as out_file:
        # Write the first line unchanged (header of the phylip file)
        out_file.write(lines[0])

        # Replace taxa names (accession numbers) with organism names in the sequences
        for line in lines[1:]:
            parts = line.split()
            accession = parts[0]
            sequence = parts[1]
            organism_name = accession_to_organism.get(accession, accession)  # Default to accession if no match
            out_file.write(f"{organism_name:<10} {sequence}\n")

if __name__ == "__main__":
    csv_file = '/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/genomes_2/metadata_genomes_2.csv'
    phylip_file = 'concatenated_names.phylip'
    output_file = 'updated_concatenated_names.phylip'

    # Load metadata from CSV
    accession_to_organism = load_metadata(csv_file)

    # Replace taxa names in the PHYLIP file
    replace_taxa_in_phylip(phylip_file, accession_to_organism, output_file)

    print(f"Updated PHYLIP file saved as {output_file}")
