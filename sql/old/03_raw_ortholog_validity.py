#!/usr/bin/env python3

"""
03_raw_ortholog_validity.py

Evaluate raw ortholog CDS sequences created by 01a_get_cds_from_compleasm_v6.py.

Each input sequence is expected to be one ortholog sequence composed of multiple
CDS fragments extracted from Compleasm output.

This script checks each sequence for:
    - length divisibility by 3
    - internal stop codons
    - terminal stop codon
    - ambiguous bases
    - valid nucleotide content
    - translated protein length using translate_dna()

The script writes one TSV summary file next to the input FASTA.

Example usage:
    python3 03_raw_ortholog_validity.py \
        /Users/rossoaa/projects/genomes/records/compleasm/GCA_003113815.1__Sphenodon_punctatus/sauropsida_odb12/sphenodon_punctatus_cds_compleasm.fasta
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import SeqIO


STOP_CODONS = {"TAA", "TAG", "TGA"}
VALID_BASES = {"A", "T", "G", "C"}
AMBIGUOUS_BASES = {"N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"}


def translate_dna(dna: str) -> str:
    """
    Translate DNA into amino acids using the standard genetic code.

    Unknown or ambiguous codons are translated as X.
    Stop codons are translated as *.
    """
    table = {
        "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
        "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
        "TCT":"S","TCC":"S","TCA":"S","TCG":"S","AGT":"S","AGC":"S",
        "CCT":"P","CCC":"P","CCA":"P","CCG":"P","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
        "GCT":"A","GCC":"A","GCA":"A","GCG":"A","TAT":"Y","TAC":"Y",
        "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","AAT":"N","AAC":"N","AAA":"K","AAG":"K",
        "GAT":"D","GAC":"D","GAA":"E","GAG":"E","TGT":"C","TGC":"C","TGG":"W",
        "CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
        "TAA":"*","TAG":"*","TGA":"*"
    }
    up = dna.upper()
    stop = len(up) - (len(up) % 3)
    return "".join(table.get(up[i:i+3], "X") for i in range(0, stop, 3))


def count_bases(sequence: str) -> Dict[str, int]:
    """
    Count canonical, ambiguous, and invalid characters in a DNA sequence.
    """
    sequence = str(sequence).upper().replace("-", "")

    a_count = sequence.count("A")
    t_count = sequence.count("T")
    g_count = sequence.count("G")
    c_count = sequence.count("C")

    n_count = sequence.count("N")

    ambiguous_count = sum(
        1 for base in sequence
        if base in AMBIGUOUS_BASES
    )

    invalid_count = sum(
        1 for base in sequence
        if base not in VALID_BASES and base not in AMBIGUOUS_BASES
    )

    total_valid_bases = a_count + t_count + g_count + c_count

    return {
        "a_count": a_count,
        "t_count": t_count,
        "g_count": g_count,
        "c_count": c_count,
        "n_count": n_count,
        "ambiguous_count": ambiguous_count,
        "invalid_count": invalid_count,
        "total_valid_bases": total_valid_bases,
    }


def check_stop_codons(sequence: str) -> Tuple[int, bool]:
    """
    Count internal stop codons and determine whether a terminal stop is present.

    The sequence is expected to already be in frame. If length is not divisible
    by 3, stop-codon results should be interpreted cautiously.
    """
    sequence = str(sequence).upper().replace("-", "")

    codons = [
        sequence[i:i + 3]
        for i in range(0, len(sequence) - 2, 3)
    ]

    if not codons:
        return 0, False

    internal_stops = sum(
        1 for codon in codons[:-1]
        if codon in STOP_CODONS
    )

    terminal_stop = codons[-1] in STOP_CODONS

    return internal_stops, terminal_stop


def translate_sequence(sequence: str, terminal_stop: bool) -> str:
    """
    Translate a DNA sequence.

    If a terminal stop codon is present, remove it before translation length
    summaries so protein_length reflects amino acids rather than the stop.
    """
    sequence = str(sequence).upper().replace("-", "")

    if terminal_stop:
        sequence = sequence[:-3]

    trimmed_length = len(sequence) - (len(sequence) % 3)
    sequence = sequence[:trimmed_length]

    if not sequence:
        return ""

    return translate_dna(sequence)


def evaluate_sequence(record_id: str, sequence: str) -> Dict[str, object]:
    """
    Evaluate one raw ortholog CDS sequence.
    """
    sequence = str(sequence).upper().replace("-", "")

    seq_length = len(sequence)
    length_mod_3 = seq_length % 3

    base_counts = count_bases(sequence)

    ambiguous_fraction = 0.0
    invalid_fraction = 0.0

    if seq_length > 0:
        ambiguous_fraction = base_counts["ambiguous_count"] / seq_length
        invalid_fraction = base_counts["invalid_count"] / seq_length

    internal_stop_count, terminal_stop = check_stop_codons(sequence)

    protein = translate_sequence(sequence, terminal_stop=terminal_stop)
    protein_length = len(protein)

    passes_length_mod_3 = length_mod_3 == 0
    passes_internal_stop_check = internal_stop_count == 0
    passes_invalid_base_check = base_counts["invalid_count"] == 0

    passes_raw_cds_qc = (
        passes_length_mod_3
        and passes_internal_stop_check
        and passes_invalid_base_check
    )

    failure_reasons: List[str] = []

    if not passes_length_mod_3:
        failure_reasons.append("length_not_divisible_by_3")

    if not passes_internal_stop_check:
        failure_reasons.append("internal_stop_codons")

    if not passes_invalid_base_check:
        failure_reasons.append("invalid_characters")

    if not failure_reasons:
        failure_reasons.append("pass")

    return {
        "sequence_id": record_id,
        "sequence_length": seq_length,
        "length_mod_3": length_mod_3,
        "internal_stop_count": internal_stop_count,
        "terminal_stop": terminal_stop,
        "protein_length": protein_length,
        "a_count": base_counts["a_count"],
        "t_count": base_counts["t_count"],
        "g_count": base_counts["g_count"],
        "c_count": base_counts["c_count"],
        "n_count": base_counts["n_count"],
        "ambiguous_count": base_counts["ambiguous_count"],
        "ambiguous_fraction": ambiguous_fraction,
        "invalid_count": base_counts["invalid_count"],
        "invalid_fraction": invalid_fraction,
        "total_valid_bases": base_counts["total_valid_bases"],
        "passes_length_mod_3": passes_length_mod_3,
        "passes_internal_stop_check": passes_internal_stop_check,
        "passes_invalid_base_check": passes_invalid_base_check,
        "passes_raw_cds_qc": passes_raw_cds_qc,
        "failure_reasons": ";".join(failure_reasons),
    }


def write_tsv(rows: List[Dict[str, object]], output_path: Path) -> None:
    """
    Write QC results to a TSV file.
    """
    if not rows:
        raise ValueError("No rows to write.")

    columns = list(rows[0].keys())

    with open(output_path, "w") as out_handle:
        out_handle.write("\t".join(columns) + "\n")

        for row in rows:
            values = [
                str(row[column])
                for column in columns
            ]
            out_handle.write("\t".join(values) + "\n")


def summarize(rows: List[Dict[str, object]]) -> None:
    """
    Print a concise terminal summary.
    """
    total_sequences = len(rows)

    passed = sum(
        1 for row in rows
        if row["passes_raw_cds_qc"]
    )

    failed = total_sequences - passed

    mod3_failed = sum(
        1 for row in rows
        if not row["passes_length_mod_3"]
    )

    internal_stop_failed = sum(
        1 for row in rows
        if not row["passes_internal_stop_check"]
    )

    invalid_base_failed = sum(
        1 for row in rows
        if not row["passes_invalid_base_check"]
    )

    terminal_stop_count = sum(
        1 for row in rows
        if row["terminal_stop"]
    )

    print("\nFinished raw ortholog CDS QC.")
    print(f"Total sequences: {total_sequences}")
    print(f"Passed raw CDS QC: {passed}")
    print(f"Failed raw CDS QC: {failed}")
    print(f"Failed length % 3 check: {mod3_failed}")
    print(f"Failed internal stop check: {internal_stop_failed}")
    print(f"Failed invalid base check: {invalid_base_failed}")
    print(f"Sequences with terminal stop codon: {terminal_stop_count}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate raw ortholog CDS sequences produced by "
            "01a_get_cds_from_compleasm_v6.py."
        )
    )
    parser.add_argument(
        "fasta",
        help="Input genus_species_cds_compleasm.fasta file.",
    )

    args = parser.parse_args()

    fasta_path = Path(args.fasta).resolve()

    if not fasta_path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {fasta_path}")

    rows = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        row = evaluate_sequence(
            record_id=record.id,
            sequence=str(record.seq),
        )
        rows.append(row)

    if not rows:
        raise ValueError(f"No FASTA records found in: {fasta_path}")

    output_path = fasta_path.with_name(
        f"{fasta_path.stem}_raw_ortholog_validity.tsv"
    )

    write_tsv(rows, output_path)
    summarize(rows)

    print(f"QC table written to: {output_path}")


if __name__ == "__main__":
    main()
