#!/usr/bin/env python3

"""
validate_intron_slicing.py

Validate intron slices produced by extract_introns_for_one_ortholog.py.

This script checks that each intron sequence in the intron TSV matches a direct
pyfaidx slice from the genome FASTA.

It also checks basic coordinate math:
    intron_length == intron_end_1_based - intron_start_1_based + 1

For minus-strand introns, the genomic slice is reverse-complemented before
comparison, because extract_introns_for_one_ortholog.py writes introns in
transcript orientation.

Example usage:
    python3 validate_intron_slicing.py \
        --genome /path/to/genome.fna \
        --introns_tsv /path/to/sphenodon_punctatus_423306at8457_introns.tsv

Validate only one intron:
    python3 validate_intron_slicing.py \
        --genome /path/to/genome.fna \
        --introns_tsv /path/to/sphenodon_punctatus_423306at8457_introns.tsv \
        --intron_id i_423306at8457_1

Print boundary windows around each intron:
    python3 validate_intron_slicing.py \
        --genome /path/to/genome.fna \
        --introns_tsv /path/to/sphenodon_punctatus_423306at8457_introns.tsv \
        --boundary_window 10
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List

import pyfaidx


def load_genome(input_fasta: Path) -> pyfaidx.Fasta:
    # 1. Index the FASTA file with pyfaidx for fast access
    try:
        genome = pyfaidx.Fasta(
            str(input_fasta),
            as_raw=True,        # slices return plain strings
            build_index=True     # makes .fai if missing
        )
        print(f"Status: Genome FASTA indexed successfully: {input_fasta}")
        return genome

    except pyfaidx.FaidxException as e:
        print(
            f"CRITICAL ERROR: Failed to load genome FASTA with pyfaidx. "
            f"Ensure the .fai file exists for {input_fasta}. Error: {e}"
        )
        sys.exit(1)

    except Exception as e:
        print(f"CRITICAL ERROR: Failed to load genome FASTA. Error: {e}")
        sys.exit(1)


def reverse_complement(seq):
    complement = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
        'D': 'H', 'H': 'D', 'N': 'N',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
        'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
        'k': 'm', 'm': 'k', 'b': 'v', 'v': 'b',
        'd': 'h', 'h': 'd', 'n': 'n'
    }

    return ''.join(complement.get(base, base) for base in reversed(seq))


def read_introns_tsv(introns_tsv: Path) -> List[Dict[str, str]]:
    """
    Read intron TSV records.
    """
    with open(introns_tsv, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    required_columns = {
        "intron_id",
        "busco_id",
        "chromosome",
        "strand",
        "intron_start_1_based",
        "intron_end_1_based",
        "intron_length",
        "sequence",
    }

    missing_columns = required_columns.difference(reader.fieldnames or [])

    if missing_columns:
        raise ValueError(
            "Intron TSV is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    return rows


def get_boundary_window(
    genome: pyfaidx.Fasta,
    chromosome: str,
    start_1_based: int,
    end_1_based: int,
    boundary_window: int,
) -> Dict[str, str]:
    """
    Return small genomic windows around the intron start and end.

    Coordinates are clipped at 1 for the left side.
    """
    left_start_1 = max(1, start_1_based - boundary_window)
    left_end_1 = start_1_based + boundary_window - 1

    right_start_1 = max(1, end_1_based - boundary_window + 1)
    right_end_1 = end_1_based + boundary_window

    left_seq = genome[chromosome][left_start_1 - 1:left_end_1]
    right_seq = genome[chromosome][right_start_1 - 1:right_end_1]

    return {
        "left_boundary_start_1_based": str(left_start_1),
        "left_boundary_end_1_based": str(left_end_1),
        "left_boundary_sequence": left_seq,
        "right_boundary_start_1_based": str(right_start_1),
        "right_boundary_end_1_based": str(right_end_1),
        "right_boundary_sequence": right_seq,
    }


def validate_intron(
    genome: pyfaidx.Fasta,
    row: Dict[str, str],
    boundary_window: int = 0,
) -> Dict[str, object]:
    """
    Validate one intron record against the genome FASTA.
    """
    intron_id = row["intron_id"]
    chromosome = row["chromosome"]
    strand = row["strand"]

    intron_start_1_based = int(row["intron_start_1_based"])
    intron_end_1_based = int(row["intron_end_1_based"])
    reported_length = int(row["intron_length"])

    reported_sequence = str(row["sequence"]).upper()

    expected_length = intron_end_1_based - intron_start_1_based + 1
    length_math_ok = expected_length == reported_length

    start_0_based = intron_start_1_based - 1
    end_0_based = intron_end_1_based

    genomic_slice = genome[chromosome][start_0_based:end_0_based].upper()

    if strand == "-":
        expected_sequence = reverse_complement(genomic_slice).upper()
    else:
        expected_sequence = genomic_slice

    sequence_match = expected_sequence == reported_sequence
    sequence_length_match = len(reported_sequence) == expected_length

    result = {
        "intron_id": intron_id,
        "busco_id": row["busco_id"],
        "chromosome": chromosome,
        "strand": strand,
        "intron_start_1_based": intron_start_1_based,
        "intron_end_1_based": intron_end_1_based,
        "reported_length": reported_length,
        "expected_length": expected_length,
        "reported_sequence_length": len(reported_sequence),
        "length_math_ok": length_math_ok,
        "sequence_length_match": sequence_length_match,
        "sequence_match": sequence_match,
        "reported_first_20": reported_sequence[:20],
        "expected_first_20": expected_sequence[:20],
        "reported_last_20": reported_sequence[-20:],
        "expected_last_20": expected_sequence[-20:],
    }

    if boundary_window > 0:
        boundary = get_boundary_window(
            genome=genome,
            chromosome=chromosome,
            start_1_based=intron_start_1_based,
            end_1_based=intron_end_1_based,
            boundary_window=boundary_window,
        )
        result.update(boundary)

    return result


def write_tsv(rows: List[Dict[str, object]], output_path: Path) -> None:
    """
    Write validation results to TSV.
    """
    if not rows:
        raise ValueError("No validation rows to write.")

    columns = list(rows[0].keys())

    with open(output_path, "w") as out_handle:
        out_handle.write("\t".join(columns) + "\n")

        for row in rows:
            out_handle.write(
                "\t".join(str(row[column]) for column in columns) + "\n"
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate intron slicing against the genome FASTA."
    )

    parser.add_argument(
        "--genome",
        required=True,
        help="Genome FASTA file used to extract introns.",
    )
    parser.add_argument(
        "--introns_tsv",
        required=True,
        help="Intron TSV produced by extract_introns_for_one_ortholog.py.",
    )
    parser.add_argument(
        "--intron_id",
        required=False,
        help="Optional intron_id to validate only one intron.",
    )
    parser.add_argument(
        "--boundary_window",
        type=int,
        default=0,
        help=(
            "Optional number of bases to print around each intron boundary. "
            "Default: 0, meaning no boundary windows are written."
        ),
    )
    parser.add_argument(
        "--out",
        required=False,
        help="Optional output TSV path. Default: <introns_tsv stem>_slicing_validation.tsv",
    )

    args = parser.parse_args()

    genome_path = Path(args.genome).resolve()
    introns_tsv = Path(args.introns_tsv).resolve()

    if not genome_path.exists():
        raise FileNotFoundError(f"Genome FASTA not found: {genome_path}")

    if not introns_tsv.exists():
        raise FileNotFoundError(f"Introns TSV not found: {introns_tsv}")

    genome = load_genome(genome_path)

    rows = read_introns_tsv(introns_tsv)

    if args.intron_id:
        rows = [
            row for row in rows
            if row["intron_id"] == args.intron_id
        ]

        if not rows:
            raise ValueError(
                f"No rows matched intron_id={args.intron_id}"
            )

    validation_rows = [
        validate_intron(
            genome=genome,
            row=row,
            boundary_window=args.boundary_window,
        )
        for row in rows
    ]

    if args.out:
        output_path = Path(args.out).resolve()
    else:
        output_path = introns_tsv.with_name(
            f"{introns_tsv.stem}_slicing_validation.tsv"
        )

    write_tsv(validation_rows, output_path)

    total = len(validation_rows)
    sequence_matches = sum(
        1 for row in validation_rows
        if row["sequence_match"]
    )
    length_matches = sum(
        1 for row in validation_rows
        if row["sequence_length_match"]
    )
    length_math_ok = sum(
        1 for row in validation_rows
        if row["length_math_ok"]
    )

    print("\nFinished validating intron slicing.")
    print(f"Introns validated: {total}")
    print(f"Length math OK: {length_math_ok}/{total}")
    print(f"Sequence length match: {length_matches}/{total}")
    print(f"Sequence match: {sequence_matches}/{total}")
    print(f"Validation TSV written to: {output_path}")


if __name__ == "__main__":
    main()
