#!/usr/bin/env python3

"""
extract_introns_for_one_ortholog.py

Extract intron sequences for one BUSCO/ortholog from a Compleasm full_table.tsv
and a genome FASTA.

This script is intentionally simple and mirrors the language/style used in
01a_get_cds_from_compleasm_v6.py.

The input ortholog is expected to be represented in Compleasm full_table.tsv,
where column 13 (index 12) contains CDS/exon coordinate tokens like:

    start_end_strand|start_end_strand|...

For a BUSCO with multiple CDS fragments, introns are inferred as the genomic
intervals between adjacent CDS fragments after sorting fragments by genomic
start coordinate.

Important:
    - Coordinates in full_table.tsv are treated as 1-based inclusive.
    - pyfaidx slices are 0-based, end-exclusive.
    - Introns are output in transcript orientation.
    - For minus-strand genes, genomic intron intervals are still calculated
      between adjacent genomic exons, but the extracted intron sequence is
      reverse-complemented so the FASTA sequence is in coding/transcript
      orientation.

Example usage:
    python3 extract_introns_for_one_ortholog.py \
        --genome /Users/rossoaa/projects/genomes/GCA_003113815.1/ncbi_dataset/data/GCA_003113815.1/GCA_003113815.1_ASM311381v1_genomic.fna \
        --full_table /Users/rossoaa/projects/genomes/records/compleasm/GCA_003113815.1__Sphenodon_punctatus/sauropsida_odb12/full_table.tsv \
        --busco_id 423306at8457 \
        --species sphenodon_punctatus
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

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


def complement(sequence):
    """Return the complement of a DNA sequence, including IUPAC ambiguity codes."""
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
        'H': 'D', 'V': 'B', 'N': 'N',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
        'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
        'k': 'm', 'm': 'k', 'b': 'v', 'd': 'h',
        'h': 'd', 'v': 'b', 'n': 'n'
    }
    return ''.join(complement.get(base, base) for base in sequence)


def calculate_gc(sequence: str) -> Tuple[float, int, int, int, int, int, int]:
    """
    Calculate GC content and nucleotide counts for a DNA sequence.

    GC content is calculated as:
        (G + C) / (A + T + G + C)

    N bases are counted separately but excluded from the denominator.
    """
    sequence = str(sequence).upper()

    g_count = sequence.count("G")
    c_count = sequence.count("C")
    a_count = sequence.count("A")
    t_count = sequence.count("T")
    n_count = sequence.count("N")

    total_valid_bases = a_count + t_count + g_count + c_count

    gc_content_float = 0.0
    if total_valid_bases > 0:
        gc_content_float = (g_count + c_count) / total_valid_bases

    return (
        gc_content_float,
        g_count,
        c_count,
        a_count,
        t_count,
        n_count,
        total_valid_bases,
    )


def parse_coordinate_token(coordinates: str, busco_id: str) -> Dict[str, object]:
    """
    Parse one Compleasm coordinate token.

    Expected token:
        start_end_strand
    """
    parts = coordinates.split("_")

    if len(parts) != 3:
        raise ValueError(
            f"WARNING: {busco_id}: unexpected exon token '{coordinates}'"
        )

    start_1_based = int(parts[0])
    end_1_based = int(parts[1])
    exon_strand = parts[2]

    exon = {
        "busco_id": busco_id,
        "start": start_1_based,
        "end": end_1_based,
        "strand": exon_strand,
    }

    return exon


def infer_introns(sorted_exons: List[Dict[str, object]], strand: str) -> List[Dict[str, object]]:
    """
    Infer intron coordinates from sorted exon/CDS fragments.

    Exons should be sorted low-to-high by genomic start coordinate.

    For both plus and minus strand genes, the genomic intron interval is:
        previous_exon_end + 1  through  next_exon_start - 1

    For minus strand genes, transcript order is opposite genomic order, so the
    final intron list is reversed for transcript-oriented numbering.
    """
    introns = []

    if len(sorted_exons) <= 1:
        return introns

    genomic_introns = []

    for i in range(len(sorted_exons) - 1):
        left_exon = sorted_exons[i]
        right_exon = sorted_exons[i + 1]

        intron_start = int(left_exon["end"]) + 1
        intron_end = int(right_exon["start"]) - 1

        # Adjacent or overlapping CDS fragments do not produce a positive-length intron.
        if intron_start > intron_end:
            continue

        genomic_introns.append({
            "busco_id": sorted_exons[0]["busco_id"],
            "genomic_left_exon": i + 1,
            "genomic_right_exon": i + 2,
            "intron_start": intron_start,
            "intron_end": intron_end,
            "intron_length": intron_end - intron_start + 1,
            "strand": strand,
        })

    if strand == "-":
        genomic_introns = list(reversed(genomic_introns))

    for ij, intron in enumerate(genomic_introns, start=1):
        intron["intron_number"] = ij
        intron["intron_id"] = f"i_{intron['busco_id']}_{ij}"

    return genomic_introns


def find_busco_row(full_table: Path, busco_id_to_find: str) -> List[str]:
    """
    Find one BUSCO row in a Compleasm full_table.tsv file.
    """
    with open(full_table, "r") as compleasmhandle:
        for line in compleasmhandle:
            line = line.rstrip()

            if not line:
                continue

            if line.startswith("#"):
                continue

            if line.startswith("Gene"):
                continue

            data = line.split("\t")

            if len(data) < 13:
                continue

            busco_id = data[0]

            if busco_id == busco_id_to_find:
                return data

    raise ValueError(
        f"BUSCO ID {busco_id_to_find} was not found in {full_table}"
    )


def extract_introns_for_busco(
    genome: pyfaidx.Fasta,
    full_table: Path,
    busco_id_to_find: str,
    species: str,
) -> List[Dict[str, object]]:
    """
    Extract intron sequences for one BUSCO ID.
    """
    data = find_busco_row(
        full_table=full_table,
        busco_id_to_find=busco_id_to_find,
    )

    status = data[1]
    busco_id = data[0]

    if status != "Single":
        raise ValueError(
            f"{busco_id} has status {status}; this first version only extracts introns for Single BUSCOs."
        )

    chromosome = data[2]
    strand = data[5]
    location = data[2]
    coordinate_data = data[12].split("|")

    sorted_exons = []

    for coordinates in coordinate_data:
        exon = parse_coordinate_token(
            coordinates=coordinates,
            busco_id=busco_id,
        )

        exon_strand = exon["strand"]

        if exon_strand != strand:
            raise ValueError(
                f"Skipping {busco_id} because exon strand != gene strand"
            )

        sorted_exons.append(exon)

    sorted_exons.sort(key=lambda x: int(x["start"]))

    introns = infer_introns(
        sorted_exons=sorted_exons,
        strand=strand,
    )

    output_rows = []

    for intron in introns:
        intron_start_1_based = int(intron["intron_start"])
        intron_end_1_based = int(intron["intron_end"])

        start_0_based = intron_start_1_based - 1
        end_0_based = intron_end_1_based

        try:
            intron_sequence = genome[chromosome][start_0_based:end_0_based]
        except KeyError:
            raise KeyError(
                f"WARNING: {busco_id}: contig '{chromosome}' not found in FASTA"
            )

        if strand == "-":
            intron_sequence = reverse_complement(intron_sequence)

        (
            gc_content_float,
            g_count,
            c_count,
            a_count,
            t_count,
            n_count,
            total_valid_bases,
        ) = calculate_gc(intron_sequence)

        output_row = {
            "intron_id": intron["intron_id"],
            "busco_id": busco_id,
            "species": species,
            "location": location,
            "chromosome": chromosome,
            "strand": strand,
            "intron_number": intron["intron_number"],
            "intron_start_1_based": intron_start_1_based,
            "intron_end_1_based": intron_end_1_based,
            "intron_length": intron["intron_length"],
            "gc": gc_content_float,
            "g_count": g_count,
            "c_count": c_count,
            "a_count": a_count,
            "t_count": t_count,
            "n_count": n_count,
            "total_valid_bases": total_valid_bases,
            "sequence": intron_sequence,
        }

        output_rows.append(output_row)

    return output_rows


def write_introns_tsv(rows: List[Dict[str, object]], output_tsv: Path) -> None:
    """
    Write intron records to TSV.
    """
    if not rows:
        columns = [
            "intron_id",
            "busco_id",
            "species",
            "location",
            "chromosome",
            "strand",
            "intron_number",
            "intron_start_1_based",
            "intron_end_1_based",
            "intron_length",
            "gc",
            "g_count",
            "c_count",
            "a_count",
            "t_count",
            "n_count",
            "total_valid_bases",
            "sequence",
        ]
    else:
        columns = list(rows[0].keys())

    with open(output_tsv, "w") as outhandle:
        outhandle.write("\t".join(columns) + "\n")

        for row in rows:
            values = [
                str(row[column])
                for column in columns
            ]
            outhandle.write("\t".join(values) + "\n")


def write_introns_fasta(rows: List[Dict[str, object]], output_fasta: Path) -> None:
    """
    Write intron records to FASTA.
    """
    with open(output_fasta, "w") as fastahandle:
        for row in rows:
            header = (
                f">{row['intron_id']} "
                f"busco_id={row['busco_id']} "
                f"species={row['species']} "
                f"chromosome={row['chromosome']} "
                f"start={row['intron_start_1_based']} "
                f"end={row['intron_end_1_based']} "
                f"strand={row['strand']} "
                f"gc={row['gc']}"
            )

            fastahandle.write(header + "\n")
            fastahandle.write(str(row["sequence"]) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Extract introns for one BUSCO/ortholog from a Compleasm "
            "full_table.tsv and genome FASTA."
        )
    )

    parser.add_argument(
        "--genome",
        required=True,
        help="Genome FASTA file.",
    )
    parser.add_argument(
        "--full_table",
        required=True,
        help="Compleasm full_table.tsv file.",
    )
    parser.add_argument(
        "--busco_id",
        required=True,
        help="BUSCO/ortholog ID to extract introns for.",
    )
    parser.add_argument(
        "--species",
        required=True,
        help="Species label to include in output.",
    )
    parser.add_argument(
        "--outdir",
        required=False,
        help="Output directory. Default: parent directory of full_table.tsv.",
    )

    args = parser.parse_args()

    genome_path = Path(args.genome).resolve()
    full_table = Path(args.full_table).resolve()

    if not genome_path.exists():
        raise FileNotFoundError(f"Genome FASTA not found: {genome_path}")

    if not full_table.exists():
        raise FileNotFoundError(f"Compleasm full_table.tsv not found: {full_table}")

    if args.outdir:
        outdir = Path(args.outdir).resolve()
    else:
        outdir = full_table.parent

    outdir.mkdir(parents=True, exist_ok=True)

    genome = load_genome(genome_path)

    rows = extract_introns_for_busco(
        genome=genome,
        full_table=full_table,
        busco_id_to_find=args.busco_id,
        species=args.species,
    )

    output_tsv = outdir / f"{args.species}_{args.busco_id}_introns.tsv"
    output_fasta = outdir / f"{args.species}_{args.busco_id}_introns.fasta"

    write_introns_tsv(rows, output_tsv)
    write_introns_fasta(rows, output_fasta)

    print("\nFinished extracting introns.")
    print(f"BUSCO ID: {args.busco_id}")
    print(f"Species: {args.species}")
    print(f"Introns extracted: {len(rows)}")
    print(f"TSV written to: {output_tsv}")
    print(f"FASTA written to: {output_fasta}")


if __name__ == "__main__":
    main()
