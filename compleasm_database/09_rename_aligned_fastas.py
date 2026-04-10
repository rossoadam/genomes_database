#!/usr/bin/env python3

import csv
import argparse
from pathlib import Path

"""
python3 10_rename_aligned_fasta.py \
  /Users/rossoaa/projects/genomes \
  --fastadir /Users/rossoaa/projects/genomes/records/compleasm/alignments/04_concatenated_clean_alignments_t100_e1_o7
"""


class RenameAlignedFastas:
    def __init__(self, genomes_dir, fasta_dir, outdir):
        self.genomes_dir = Path(genomes_dir).resolve()
        self.fasta_dir = Path(fasta_dir).resolve()
        self.outdir = Path(outdir).resolve()

        self.metadata_csv = self.genomes_dir / "records/compleasm/records/metadata.csv"
        self.accession_to_organism = {}
        self.fasta_list = []

    def load_metadata(self):
        """Load accession -> genus_species mapping from metadata.csv"""
        with open(self.metadata_csv, "r", newline="") as file:
            reader = csv.DictReader(file)
            for row in reader:
                accession = row["accession"].strip()
                organism_name = row["genus_species"].strip()
                self.accession_to_organism[accession] = organism_name

    def load_fasta_list(self):
        """Collect fasta files from the input directory"""
        self.fasta_list = sorted(
            [
                path for path in self.fasta_dir.iterdir()
                if path.is_file() and path.suffix in {".fasta", ".fa", ".fas", ".fna"}
            ]
        )

    def rename_headers_in_fasta(self, fasta_path):
        """
        Read one FASTA and return a list of (header, sequence) tuples
        with accession headers replaced by genus_species when possible.
        """
        renamed_records = []
        current_header = None
        current_seq = []

        with open(fasta_path, "r") as handle:
            for line in handle:
                line = line.rstrip()

                if not line:
                    continue

                if line.startswith(">"):
                    if current_header is not None:
                        renamed_records.append((current_header, "".join(current_seq)))

                    raw_header = line[1:].strip()
                    accession = raw_header.split()[0]

                    new_header = self.accession_to_organism.get(accession, accession)
                    current_header = new_header
                    current_seq = []
                else:
                    current_seq.append(line)

        if current_header is not None:
            renamed_records.append((current_header, "".join(current_seq)))

        return renamed_records

    def write_new_fasta(self, fasta_path):
        """Write renamed FASTA to output directory"""
        renamed_records = self.rename_headers_in_fasta(fasta_path)

        output_path = self.outdir / fasta_path.name
        with open(output_path, "w") as out_handle:
            for header, seq in renamed_records:
                out_handle.write(f">{header}\n")
                out_handle.write(f"{seq}\n")

        print(f"Wrote renamed fasta: {output_path}")

    def loop(self):
        """Rename all FASTA files in the input directory"""
        for fasta_path in self.fasta_list:
            print(f"Renaming {fasta_path.name}")
            self.write_new_fasta(fasta_path)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Rename aligned FASTA taxa from accession to genus_species "
            "using genomes/records/compleasm/records/metadata.csv"
        )
    )

    parser.add_argument(
        "genomes_dir",
        help="Path to genomes directory"
    )

    parser.add_argument(
        "--fastadir",
        required=True,
        help="Path to input FASTA directory"
    )

    parser.add_argument(
        "--outdir",
        help="Optional output directory "
             "(default: genomes/records/compleasm/alignments/07_renamed_fastas)"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    genomes_dir = Path(args.genomes_dir).resolve()
    metadata_csv = genomes_dir / "records/compleasm/records/metadata.csv"

    if not metadata_csv.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_csv}")

    fasta_dir = Path(args.fastadir).resolve()
    if not fasta_dir.exists():
        raise FileNotFoundError(f"FASTA directory not found: {fasta_dir}")
    if not fasta_dir.is_dir():
        raise NotADirectoryError(f"Expected FASTA directory, got: {fasta_dir}")

    if args.outdir:
        outdir = Path(args.outdir).resolve()
    else:
        outdir = genomes_dir / "records/compleasm/alignments/07_renamed_fastas"

    outdir.mkdir(parents=True, exist_ok=True)

    renamer = RenameAlignedFastas(
        genomes_dir=genomes_dir,
        fasta_dir=fasta_dir,
        outdir=outdir,
    )

    renamer.load_metadata()
    renamer.load_fasta_list()
    renamer.loop()