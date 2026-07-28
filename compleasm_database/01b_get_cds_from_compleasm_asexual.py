#!/usr/bin/env python3

"""
Extract one coding sequence per BUSCO ID from a genome using a Compleasm
full_table.tsv file.

Both Single and Duplicated BUSCO records are considered. When multiple valid
records have the same BUSCO ID, only the first valid sequence encountered in
the Compleasm table is written.

Usage:
    python3 get_cds_from_compleasm_v6.py \
        genome.fasta \
        full_table.tsv \
        species_name
"""

import sys
from datetime import datetime
from pathlib import Path

import pyfaidx


def load_genome(input_fasta: Path) -> pyfaidx.Fasta:
    """Index and load the genome FASTA with pyfaidx."""
    try:
        genome = pyfaidx.Fasta(
            str(input_fasta),
            as_raw=True,
            build_index=True,
        )

        print(f"Status: Genome FASTA indexed successfully: {input_fasta}")
        return genome

    except pyfaidx.FaidxException as error:
        print(
            "CRITICAL ERROR: Failed to load genome FASTA with pyfaidx. "
            f"Ensure the FASTA is readable and can be indexed: {input_fasta}. "
            f"Error: {error}",
            file=sys.stderr,
        )
        sys.exit(1)

    except Exception as error:
        print(
            f"CRITICAL ERROR: Failed to load genome FASTA. Error: {error}",
            file=sys.stderr,
        )
        sys.exit(1)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement_table = str.maketrans(
        {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
            "R": "Y",
            "Y": "R",
            "S": "S",
            "W": "W",
            "K": "M",
            "M": "K",
            "B": "V",
            "V": "B",
            "D": "H",
            "H": "D",
            "N": "N",
            "a": "t",
            "t": "a",
            "g": "c",
            "c": "g",
            "r": "y",
            "y": "r",
            "s": "s",
            "w": "w",
            "k": "m",
            "m": "k",
            "b": "v",
            "v": "b",
            "d": "h",
            "h": "d",
            "n": "n",
        }
    )

    return sequence.translate(complement_table)[::-1]


def parse_coordinate_data(
    coordinate_string: str,
    busco_id: str,
    record_handle,
):
    """
    Parse Compleasm exon coordinates.

    Expected token format:
        start_end_strand

    Multiple exon tokens are separated by:
        |
    """
    parsed_coordinates = []

    for coordinate_token in coordinate_string.split("|"):
        coordinate_token = coordinate_token.strip()

        if not coordinate_token:
            print(
                f"WARNING: {busco_id}: encountered an empty exon token.",
                file=record_handle,
            )
            return None

        parts = coordinate_token.split("_")

        if len(parts) != 3:
            print(
                f"WARNING: {busco_id}: unexpected exon token "
                f"'{coordinate_token}'. Expected start_end_strand.",
                file=record_handle,
            )
            return None

        try:
            start_1_based = int(parts[0])
            end_1_based = int(parts[1])
        except ValueError:
            print(
                f"WARNING: {busco_id}: non-integer coordinates in exon token "
                f"'{coordinate_token}'.",
                file=record_handle,
            )
            return None

        exon_strand = parts[2]

        if start_1_based < 1:
            print(
                f"WARNING: {busco_id}: exon start must be at least 1 in "
                f"'{coordinate_token}'.",
                file=record_handle,
            )
            return None

        if end_1_based < start_1_based:
            print(
                f"WARNING: {busco_id}: exon end precedes exon start in "
                f"'{coordinate_token}'.",
                file=record_handle,
            )
            return None

        parsed_coordinates.append(
            {
                "start_1_based": start_1_based,
                "end_1_based": end_1_based,
                "strand": exon_strand,
                "token": coordinate_token,
            }
        )

    # Build the genomic sequence from low to high coordinates.
    # Negative-strand genes are reverse-complemented after concatenation.
    parsed_coordinates.sort(key=lambda exon: exon["start_1_based"])

    return parsed_coordinates


def main() -> None:
    if len(sys.argv) != 4:
        print(
            "Usage:\n"
            "    python3 get_cds_from_compleasm_v6.py "
            "<genome.fasta> <full_table.tsv> <species_name>",
            file=sys.stderr,
        )
        sys.exit(1)

    genome_path = Path(sys.argv[1]).expanduser().resolve()
    compleasm_path = Path(sys.argv[2]).expanduser().resolve()
    species = sys.argv[3]

    if not genome_path.is_file():
        print(
            f"CRITICAL ERROR: Genome FASTA not found: {genome_path}",
            file=sys.stderr,
        )
        sys.exit(1)

    if not compleasm_path.is_file():
        print(
            f"CRITICAL ERROR: Compleasm table not found: {compleasm_path}",
            file=sys.stderr,
        )
        sys.exit(1)

    genome = load_genome(genome_path)

    outdir = compleasm_path.parent
    cds_filename = outdir / f"{species}_cds_compleasm.tsv"
    cds_fasta = outdir / f"{species}_cds_compleasm.fasta"
    record_filename = outdir / f"{species}_cds_compleasm_records.txt"

    # This set represents the BUSCO headers already written to the output FASTA.
    # Checking a set is much faster than rescanning the FASTA for every row.
    written_busco_ids = set()

    total_rows = 0
    eligible_rows = 0
    written_sequences = 0
    skipped_duplicate_ids = 0
    skipped_invalid_records = 0

    try:
        with (
            open(cds_filename, "w") as output_handle,
            open(cds_fasta, "w") as fasta_handle,
            open(compleasm_path, "r") as compleasm_handle,
            open(record_filename, "w") as record_handle,
        ):
            output_handle.write(
                "Busco ID\tspecies\tlocation\tstrand\tsequence\n"
            )

            start_time = datetime.now()

            print(f"Run started: {start_time}", file=record_handle)
            print(f"Genome FASTA: {genome_path}", file=record_handle)
            print(f"Compleasm table: {compleasm_path}", file=record_handle)
            print(f"Species label: {species}", file=record_handle)
            print("", file=record_handle)

            print("Running BUSCO CDS extraction...")

            for line_number, line in enumerate(compleasm_handle, start=1):
                line = line.rstrip("\n\r")

                if not line:
                    continue

                # Skip common header and comment lines.
                if line.startswith("Gene") or line.startswith("#"):
                    continue

                total_rows += 1
                data = line.split("\t")

                if len(data) < 13:
                    skipped_invalid_records += 1
                    print(
                        f"WARNING: Line {line_number}: skipping malformed line "
                        f"with {len(data)} columns: {line}",
                        file=record_handle,
                    )
                    continue

                busco_id = data[0].strip()
                status = data[1].strip()

                if status not in {"Single", "Duplicated"}:
                    print(
                        f"{busco_id} is {status}; not extracted.",
                        file=record_handle,
                    )
                    continue

                eligible_rows += 1

                # Do not extract another sequence if this BUSCO ID has already
                # been successfully written.
                if busco_id in written_busco_ids:
                    skipped_duplicate_ids += 1
                    print(
                        f"Skipping additional {status} copy of {busco_id}: "
                        "a valid sequence for this BUSCO ID has already been "
                        "written.",
                        file=record_handle,
                    )
                    continue

                chromosome = data[2].strip()
                strand = data[5].strip()
                coordinate_string = data[12].strip()
                location = chromosome

                if strand not in {"+", "-"}:
                    skipped_invalid_records += 1
                    print(
                        f"WARNING: {busco_id}: unexpected gene strand "
                        f"'{strand}'; skipping.",
                        file=record_handle,
                    )
                    continue

                parsed_coordinates = parse_coordinate_data(
                    coordinate_string=coordinate_string,
                    busco_id=busco_id,
                    record_handle=record_handle,
                )

                if parsed_coordinates is None:
                    skipped_invalid_records += 1
                    continue

                sequence_parts = []
                good_busco = True

                for exon_data in parsed_coordinates:
                    exon_strand = exon_data["strand"]

                    if exon_strand != strand:
                        print(
                            f"WARNING: Skipping {busco_id}: exon strand "
                            f"'{exon_strand}' does not match gene strand "
                            f"'{strand}' in token '{exon_data['token']}'.",
                            file=record_handle,
                        )
                        good_busco = False
                        break

                    start_0_based = exon_data["start_1_based"] - 1
                    end_0_based = exon_data["end_1_based"]

                    try:
                        exon_sequence = genome[chromosome][
                            start_0_based:end_0_based
                        ]
                    except KeyError:
                        print(
                            f"WARNING: {busco_id}: contig '{chromosome}' was "
                            "not found in the genome FASTA; skipping this copy.",
                            file=record_handle,
                        )
                        good_busco = False
                        break
                    except Exception as error:
                        print(
                            f"WARNING: {busco_id}: failed to extract "
                            f"{chromosome}:{start_0_based + 1}-{end_0_based}. "
                            f"Error: {error}",
                            file=record_handle,
                        )
                        good_busco = False
                        break

                    sequence_parts.append(str(exon_sequence))

                if not good_busco:
                    skipped_invalid_records += 1
                    continue

                sequence = "".join(sequence_parts)

                if not sequence:
                    skipped_invalid_records += 1
                    print(
                        f"WARNING: {busco_id}: extracted sequence was empty; "
                        "skipping this copy.",
                        file=record_handle,
                    )
                    continue

                if strand == "-":
                    sequence = reverse_complement(sequence)

                output_handle.write(
                    f"{busco_id}\t{species}\t{location}\t{strand}\t"
                    f"{sequence}\n"
                )
                fasta_handle.write(f">{busco_id}\n{sequence}\n")

                # Add the ID only after successful sequence extraction and
                # writing. A later copy can therefore be tried if this copy
                # failed.
                written_busco_ids.add(busco_id)
                written_sequences += 1

                print(
                    f"Wrote {busco_id} ({status}) from "
                    f"{chromosome}, strand {strand}, length {len(sequence)} bp.",
                    file=record_handle,
                )

            end_time = datetime.now()

            summary_lines = [
                "",
                "Extraction summary",
                "------------------",
                f"Run started: {start_time}",
                f"Run finished: {end_time}",
                f"Compleasm data rows examined: {total_rows}",
                f"Single or Duplicated rows considered: {eligible_rows}",
                f"Unique BUSCO sequences written: {written_sequences}",
                (
                    "Additional copies skipped after BUSCO was written: "
                    f"{skipped_duplicate_ids}"
                ),
                f"Invalid records or failed extractions: {skipped_invalid_records}",
                f"TSV output: {cds_filename}",
                f"FASTA output: {cds_fasta}",
            ]

            for summary_line in summary_lines:
                print(summary_line, file=record_handle)

        print("Extraction complete.")
        print(f"Unique BUSCO sequences written: {written_sequences}")
        print(
            "Additional duplicated BUSCO copies skipped: "
            f"{skipped_duplicate_ids}"
        )
        print(f"Failed or invalid records: {skipped_invalid_records}")
        print(f"FASTA: {cds_fasta}")
        print(f"TSV:   {cds_filename}")
        print(f"Log:   {record_filename}")

    finally:
        genome.close()


if __name__ == "__main__":
    main()
