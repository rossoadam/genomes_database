#!/usr/bin/env python3

import csv
import re
import argparse
from pathlib import Path
"""
python3 12_rename_treefile.py \
  /Users/rossoaa/projects/genomes \
  --tree /Users/rossoaa/projects/genomes/records/compleasm/alignments/06_rooted_tree_nonreversible/nonrev_dna.treefile
"""

def load_metadata(metadata_csv):
    accession_to_species = {}
    with open(metadata_csv, "r", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            accession = row["accession"].strip()
            genus_species = row["genus_species"].strip().replace(" ", "_")
            accession_to_species[accession] = genus_species
    return accession_to_species


def build_accession_pattern(accession_to_species):
    accessions = sorted(accession_to_species.keys(), key=len, reverse=True)
    escaped = [re.escape(acc) for acc in accessions]
    pattern = r'(?<=[(,])(' + "|".join(escaped) + r')(?=[:),;\s])'
    return re.compile(pattern)


def rename_tree_labels(tree_text, accession_to_species, pattern):
    def replacer(match):
        accession = match.group(1)
        return accession_to_species.get(accession, accession)
    return pattern.sub(replacer, tree_text)


def remove_branch_lengths(tree_text):
    """
    Remove branch lengths like :0.123, :1e-05, :-0.02
    """
    return re.sub(r':-?\d+(\.\d+)?([eE][-+]?\d+)?', '', tree_text)


def remove_bootstrap_values(tree_text):
    """
    Remove internal node support labels that appear immediately after ')'
    Examples:
      )95
      )100
      )99.5
      )1e2
    """
    return re.sub(r'\)(\d+(\.\d+)?([eE][-+]?\d+)?)', ')', tree_text)


def clean_tree_for_nhphyml(tree_text):
    tree_text = remove_branch_lengths(tree_text)
    tree_text = remove_bootstrap_values(tree_text)
    return tree_text


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Rename accession labels in a Newick tree using metadata.csv, "
            "and remove branch lengths and bootstrap/support values for NHPhyML."
        )
    )
    parser.add_argument("genomes_dir", help="Path to genomes directory")
    parser.add_argument("--tree", required=True, help="Path to input tree file")
    parser.add_argument(
        "--outdir",
        help="Optional output directory "
             "(default: genomes/records/compleasm/alignments/07_renamed_treefile)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    genomes_dir = Path(args.genomes_dir).resolve()
    tree_file = Path(args.tree).resolve()

    metadata_csv = genomes_dir / "records/compleasm/records/metadata.csv"
    if not metadata_csv.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_csv}")

    if not tree_file.exists():
        raise FileNotFoundError(f"Tree file not found: {tree_file}")

    if args.outdir:
        outdir = Path(args.outdir).resolve()
    else:
        outdir = genomes_dir / "records/compleasm/alignments/08_renamed_treefile"

    outdir.mkdir(parents=True, exist_ok=True)
    output_tree = outdir / tree_file.name

    accession_to_species = load_metadata(metadata_csv)
    pattern = build_accession_pattern(accession_to_species)

    with open(tree_file, "r") as handle:
        tree_text = handle.read().strip()

    renamed_tree = rename_tree_labels(tree_text, accession_to_species, pattern)
    cleaned_tree = clean_tree_for_nhphyml(renamed_tree)

    with open(output_tree, "w") as handle:
        handle.write(cleaned_tree + "\n")

    print(f"Metadata: {metadata_csv}")
    print(f"Input tree: {tree_file}")
    print(f"Output tree: {output_tree}")
    print("Done. Tree written in topology-only form for NHPhyML.")


if __name__ == "__main__":
    main()