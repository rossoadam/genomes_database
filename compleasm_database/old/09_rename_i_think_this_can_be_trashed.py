#!/usr/bin/env python3

import csv
import argparse
from pathlib import Path


# Load metadata from CSV file
def load_metadata(csv_file):
    accession_to_organism = {}
    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            accession = row['accession']
            organism_name = row['organism_name']
            accession_to_organism[accession] = organism_name
    return accession_to_organism


# Replace accession numbers with organism names in the PHYLIP file
def replace_taxa_in_phylip(phylip_dir, accession_to_organism, output_file):
    phylip_file = phylip_dir / "concatenated.phy"
    with open(phylip_file, 'r') as file:
        lines = file.readlines()

    with open(output_file, 'w') as out_file:
        # Write header unchanged
        out_file.write(lines[0])

        for line in lines[1:]:
            parts = line.split()
            accession = parts[0]
            sequence = parts[1]

            organism_name = accession_to_organism.get(accession, accession)

            # Replace spaces with underscores for PHYLIP safety
            organism_name = organism_name.replace(" ", "_")

            out_file.write(f"{organism_name:<15} {sequence}\n")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Rename PHYLIP taxa from accessions to organism names "
            "(genus_species) using metadata.csv"
        )
    )

    parser.add_argument(
        "genomes_dir",
        help="Path to genomes directory (must end with 'genomes')"
    )

    parser.add_argument(
        "--phylip",
        required=True,
        help="Path to input PHYLIP file"
    )

    parser.add_argument(
        "--outdir",
        help="Optional output directory (default: genomes/records/.../07_renamed)"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    genomes_dir = Path(args.genomes_dir).resolve()

    # Automatically locate metadata
    metadata_csv = genomes_dir / "records/compleasm/records/metadata.csv"

    if not metadata_csv.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_csv}")

    # Input PHYLIP
    phylip_dir = Path(args.phylip).resolve()

    if not phylip_dir.exists():
        raise FileNotFoundError(f"PHYLIP dir not found: {phylip_dir}")

    # Output directory
    if args.outdir:
        outdir = Path(args.outdir).resolve()
    else:
        outdir = genomes_dir / "records/compleasm/alignments/07_renamed"

    outdir.mkdir(parents=True, exist_ok=True)

    # Output file name
    output_file = outdir / f"concatenated_renamed.phy"

    print(f"Metadata: {metadata_csv}")
    print(f"Input PHYLIP directory: {phylip_dir}")
    print(f"Output: {output_file}")

    # Load metadata
    accession_to_organism = load_metadata(metadata_csv)

    # Replace taxa
    replace_taxa_in_phylip(phylip_dir, accession_to_organism, output_file)

    print(f"Updated PHYLIP file saved as {output_file}")