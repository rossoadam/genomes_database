#!/usr/bin/env python3
"""
00c_build_tsvs_from_gff.py

Build an intergenic-sequence (IGS) TSV from project GFF and FASTA files.

This script is a focused refactor of the legacy
20251111_numos-genbank_introns.py workflow. It performs only the steps needed
to identify intergenic regions between genes, calculate their GC content, map
them to sequence_pk values, and write records directly to igs.tsv.

Workflow
--------
1. Read the project manifest.
2. Left join manifest accessions to genomes_metadata.csv on accession.
3. Read gene coordinates from each available GFF file.
4. Group genes by GFF sequence ID (chromosome/contig ID).
5. Sort and merge overlapping gene intervals.
6. Define IGS intervals between consecutive merged gene intervals.
7. Match each GFF sequence ID to sequences.tsv using:
       GFF chrom_id == sequences.tsv sequence_id
8. Extract each IGS sequence with pyfaidx and calculate GC content as:
       (G + C) / (A + T + G + C)
   N and other non-ATGC characters are excluded from the denominator.
9. Write:
       <genomes_dir>/records/sql_tsvs/igs.tsv

Important behavior
------------------
- Empty or missing path_to_gff values are skipped with a warning; they do not
  terminate the run.
- Empty or missing path_to_fna values do not terminate the run. IGS coordinates
  are still written when sequence_pk can be resolved, but gc is left blank.
- GFF sequence IDs that cannot be mapped uniquely to sequence_pk are reported
  and skipped rather than writing an invalid foreign key.
- Only IGSs BETWEEN consecutive merged gene intervals are emitted. Regions
  before the first gene and after the last gene are not included.
- --debug has been removed and debug_rows has been replaced with mapping_rows

Example
-------
python 00c_build_tsvs_from_gff.py /path/to/genomes --debug
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


IGS_COLUMNS = [
    "igs_pk",
    "sequence_pk",
    "igs_seq_start",
    "igs_seq_end",
    "length",
    "gc",
]

def warn(message: str) -> None:
    """Write a warning message to stderr."""
    print(f"WARNING: {message}", file=sys.stderr)


def calculate_gc(sequence: str) -> Tuple[float, int, int, int, int, int, int]:
    """
    Calculate GC content and nucleotide counts for a DNA sequence.

    This function incorporates the logic from calculate_gc.py.

    GC content is:
        (G + C) / (A + T + G + C)

    N bases are counted separately and excluded from the denominator.
    Other non-ATGC characters are also excluded because the denominator is
    calculated explicitly from A, T, G, and C counts.

    Returns
    -------
    tuple
        gc_content_float,
        g_count,
        c_count,
        a_count,
        t_count,
        n_count,
        total_valid_bases
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


def open_text_auto(path: Path):
    """Open a plain-text or gzip-compressed text file."""
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def require_columns(
    fieldnames: Optional[Sequence[str]],
    required: Iterable[str],
    path: Path,
) -> None:
    """Raise a clear error when a required input column is missing."""
    available = set(fieldnames or [])
    missing = [column for column in required if column not in available]
    if missing:
        raise ValueError(
            f"{path} is missing required column(s): {', '.join(missing)}"
        )


def read_delimited_rows(
    path: Path,
    delimiter: str,
    required_columns: Iterable[str],
) -> List[Dict[str, str]]:
    """Read a CSV/TSV into a list of string-valued dictionaries."""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        require_columns(reader.fieldnames, required_columns, path)
        return [
            {key: (value if value is not None else "") for key, value in row.items()}
            for row in reader
        ]


def clean(value: object) -> str:
    """Return a stripped string, treating None as empty."""
    if value is None:
        return ""
    return str(value).strip()


def resolve_project_path(raw_path: object, genomes_dir: Path) -> Optional[Path]:
    """
    Resolve a path stored in metadata.

    Absolute paths are used directly. Relative paths are interpreted relative
    to genomes_dir. Blank values return None.
    """
    text = clean(raw_path)
    if not text:
        return None

    path = Path(text).expanduser()
    if path.is_absolute():
        return path
    return genomes_dir / path


def left_join_manifest_metadata(
    manifest_rows: Sequence[Mapping[str, str]],
    metadata_rows: Sequence[Mapping[str, str]],
) -> List[Dict[str, str]]:
    """
    Left join manifest rows to genomes metadata on accession.

    Metadata values are added to the manifest row. When the same column occurs
    in both files, a non-empty metadata value takes precedence.
    """
    metadata_by_accession: Dict[str, Dict[str, str]] = {}

    for row in metadata_rows:
        accession = clean(row.get("accession"))
        if not accession:
            continue

        if accession in metadata_by_accession:
            warn(
                f"Duplicate accession {accession!r} found in genomes metadata; "
                "using the first row."
            )
            continue

        metadata_by_accession[accession] = dict(row)

    joined: List[Dict[str, str]] = []
    seen_accessions = set()

    for manifest_row in manifest_rows:
        accession = clean(manifest_row.get("accession"))
        if not accession:
            warn("Manifest row with empty accession was skipped.")
            continue

        # The output is accession-driven, so process a repeated accession once.
        if accession in seen_accessions:
            warn(
                f"Accession {accession!r} occurs more than once in the manifest; "
                "processing it once."
            )
            continue
        seen_accessions.add(accession)

        merged = dict(manifest_row)
        metadata_row = metadata_by_accession.get(accession)

        if metadata_row is None:
            warn(
                f"Accession {accession!r} was not found in genomes_metadata.csv."
            )
        else:
            for key, value in metadata_row.items():
                if clean(value):
                    merged[key] = value

        joined.append(merged)

    return joined


def load_sequence_index(
    sequences_path: Path,
) -> Dict[str, List[Dict[str, str]]]:
    """
    Index sequences.tsv by sequence_id.

    A list is retained for each sequence_id so ambiguous duplicate IDs can be
    detected rather than silently mapped to the wrong sequence_pk.
    """
    rows = read_delimited_rows(
        sequences_path,
        delimiter="\t",
        required_columns=["sequence_pk", "sequence_id"],
    )

    index: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in rows:
        sequence_id = clean(row.get("sequence_id"))
        sequence_pk = clean(row.get("sequence_pk"))

        if not sequence_id or not sequence_pk:
            continue

        index[sequence_id].append(row)

    return index


def resolve_sequence_record(
    chrom_id: str,
    sequence_index: Mapping[str, Sequence[Mapping[str, str]]],
    metadata_genome_pk: str = "",
) -> Tuple[Optional[Mapping[str, str]], str]:
    """
    Resolve a GFF chrom_id to one sequences.tsv row.

    Primary rule:
        chrom_id == sequence_id

    If the same sequence_id occurs multiple times and genomes metadata contains
    genome_pk, genome_pk is used only to disambiguate the exact sequence-ID
    matches.

    Returns
    -------
    (record, status)
        status is one of: matched, unmatched, ambiguous.
    """
    candidates = list(sequence_index.get(chrom_id, []))

    if not candidates:
        return None, "unmatched"

    if len(candidates) == 1:
        return candidates[0], "matched"

    if metadata_genome_pk:
        narrowed = [
            row
            for row in candidates
            if clean(row.get("genome_pk")) == metadata_genome_pk
        ]
        if len(narrowed) == 1:
            return narrowed[0], "matched"

    return None, "ambiguous"


def parse_gene_intervals(gff_path: Path) -> Dict[str, List[Tuple[int, int]]]:
    """
    Parse gene features from a GFF/GFF3 file and group intervals by seqid.

    GFF coordinates are retained as 1-based inclusive coordinates.
    """
    genes_by_chrom: Dict[str, List[Tuple[int, int]]] = defaultdict(list)

    with open_text_auto(gff_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                warn(
                    f"{gff_path}: line {line_number} has fewer than 9 GFF fields; "
                    "skipping."
                )
                continue

            chrom_id = clean(fields[0])
            feature_type = clean(fields[2]).lower()

            if feature_type != "gene":
                continue

            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                warn(
                    f"{gff_path}: line {line_number} has non-integer gene "
                    "coordinates; skipping."
                )
                continue

            if not chrom_id:
                warn(
                    f"{gff_path}: line {line_number} has an empty sequence ID; "
                    "skipping."
                )
                continue

            if start < 1 or end < 1:
                warn(
                    f"{gff_path}: line {line_number} has a coordinate < 1; "
                    "skipping."
                )
                continue

            # GFF should already have start <= end, but normalizing here avoids
            # a malformed record breaking interval sorting.
            if start > end:
                start, end = end, start

            genes_by_chrom[chrom_id].append((start, end))

    return genes_by_chrom


def fix_overlap(intervals: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Sort gene intervals and merge overlapping bounds.

    This preserves the intent of the legacy fix_overlap function. Adjacent
    intervals are not merged unless they actually overlap.
    """
    if not intervals:
        return []

    sorted_intervals = sorted(intervals, key=lambda item: (item[0], item[1]))
    merged: List[Tuple[int, int]] = [sorted_intervals[0]]

    for current_start, current_end in sorted_intervals[1:]:
        last_start, last_end = merged[-1]

        if last_end >= current_start:
            merged[-1] = (last_start, max(last_end, current_end))
        else:
            merged.append((current_start, current_end))

    return merged


def calculate_igs_intervals(
    merged_gene_intervals: Sequence[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    """
    Return IGS coordinates between consecutive merged gene intervals.

    Coordinates are 1-based and inclusive.
    """
    igs_intervals: List[Tuple[int, int]] = []

    for left, right in zip(
        merged_gene_intervals,
        merged_gene_intervals[1:],
    ):
        igs_start = left[1] + 1
        igs_end = right[0] - 1

        # Adjacent genes produce no intergenic bases.
        if igs_start <= igs_end:
            igs_intervals.append((igs_start, igs_end))

    return igs_intervals


def open_fasta(fna_path: Path):
    """
    Open a FASTA with pyfaidx.

    Import is delayed so a project with missing FASTA paths can still complete
    coordinate generation without failing at program startup.
    """
    try:
        from pyfaidx import Fasta
    except ImportError as exc:
        raise RuntimeError(
            "pyfaidx is required to calculate GC content. Install pyfaidx in "
            "the environment used to run this script."
        ) from exc

    return Fasta(str(fna_path))


def igs_gc(
    fasta,
    chrom_id: str,
    start_1based: int,
    end_1based: int,
) -> float:
    """
    Extract an IGS substring with pyfaidx and return its GC ratio.

    GFF/IGS coordinates:
        1-based, inclusive

    pyfaidx slicing:
        0-based start, end-exclusive
    """
    sequence = str(fasta[chrom_id][start_1based - 1 : end_1based])
    gc_value, *_ = calculate_gc(sequence)
    return gc_value


def build_igs_rows(
    joined_genomes: Sequence[Mapping[str, str]],
    sequence_index: Mapping[str, Sequence[Mapping[str, str]]],
    genomes_dir: Path,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    """
    Build IGS data rows and sequence-mapping debug rows in memory.
    """
    igs_rows: List[Dict[str, object]] = []
    mapping_rows: List[Dict[str, object]] = []

    for genome_row in joined_genomes:
        accession = clean(genome_row.get("accession"))
        metadata_genome_pk = clean(genome_row.get("genome_pk"))

        gff_path = resolve_project_path(genome_row.get("path_to_gff"), genomes_dir)
        fna_path = resolve_project_path(genome_row.get("path_to_fna"), genomes_dir)

        if gff_path is None:
            warn(
                f"{accession}: path_to_gff is empty; accession skipped."
            )
            mapping_rows.append(
                {
                    "accession": accession,
                    "chrom_id": "",
                    "sequence_id": "",
                    "matched_sequence_pk": "",
                    "igs_count": 0,
                }
            )
            continue

        if not gff_path.is_file():
            warn(
                f"{accession}: GFF file not found at {gff_path}; accession skipped."
            )
            mapping_rows.append(
                {
                    "accession": accession,
                    "chrom_id": "",
                    "sequence_id": "",
                    "matched_sequence_pk": "",
                    "igs_count": 0,
                }
            )
            continue

        genes_by_chrom = parse_gene_intervals(gff_path)

        if not genes_by_chrom:
            warn(f"{accession}: no gene features were found in {gff_path}.")

        fasta = None
        fasta_available = False

        if fna_path is None:
            warn(
                f"{accession}: path_to_fna is empty; GC values will be blank."
            )
        elif not fna_path.is_file():
            warn(
                f"{accession}: FASTA file not found at {fna_path}; "
                "GC values will be blank."
            )
        else:
            try:
                fasta = open_fasta(fna_path)
                fasta_available = True
            except Exception as exc:
                warn(
                    f"{accession}: could not open FASTA {fna_path}: {exc}. "
                    "GC values will be blank."
                )

        try:
            for chrom_id in sorted(genes_by_chrom):
                merged_genes = fix_overlap(genes_by_chrom[chrom_id])
                intervals = calculate_igs_intervals(merged_genes)

                sequence_record, match_status = resolve_sequence_record(
                    chrom_id,
                    sequence_index,
                    metadata_genome_pk=metadata_genome_pk,
                )

                if sequence_record is None:
                    if match_status == "ambiguous":
                        warn(
                            f"{accession}: chrom_id {chrom_id!r} matches more than "
                            "one sequences.tsv row and could not be disambiguated; "
                            f"{len(intervals)} IGS row(s) skipped."
                        )
                    else:
                        warn(
                            f"{accession}: chrom_id {chrom_id!r} was not found in "
                            f"sequences.tsv; {len(intervals)} IGS row(s) skipped."
                        )

                    mapping_rows.append(
                        {
                            "accession": accession,
                            "chrom_id": chrom_id,
                            "sequence_id": "",
                            "matched_sequence_pk": "",
                            "igs_count": len(intervals),
                        }
                    )
                    continue

                sequence_pk = clean(sequence_record.get("sequence_pk"))
                sequence_id = clean(sequence_record.get("sequence_id"))

                mapping_rows.append(
                    {
                        "accession": accession,
                        "chrom_id": chrom_id,
                        "sequence_id": sequence_id,
                        "matched_sequence_pk": sequence_pk,
                        "igs_count": len(intervals),
                    }
                )

                chrom_in_fasta = False
                if fasta_available and fasta is not None:
                    try:
                        chrom_in_fasta = chrom_id in fasta
                    except Exception:
                        chrom_in_fasta = False

                    if not chrom_in_fasta and intervals:
                        warn(
                            f"{accession}: chrom_id {chrom_id!r} is not present "
                            "in the FASTA; GC values for this sequence will be blank."
                        )

                for igs_start, igs_end in intervals:
                    gc_value: object = ""

                    if fasta_available and fasta is not None and chrom_in_fasta:
                        try:
                            gc_value = igs_gc(
                                fasta,
                                chrom_id,
                                igs_start,
                                igs_end,
                            )
                        except Exception as exc:
                            warn(
                                f"{accession}: failed to calculate GC for "
                                f"{chrom_id}:{igs_start}-{igs_end}: {exc}. "
                                "GC left blank."
                            )
                            gc_value = ""

                    igs_rows.append(
                        {
                            "sequence_pk": int(sequence_pk),
                            "igs_seq_start": igs_start,
                            "igs_seq_end": igs_end,
                            "length": igs_end - igs_start + 1,
                            "gc": gc_value,
                        }
                    )
        finally:
            if fasta is not None and hasattr(fasta, "close"):
                fasta.close()

    # Assign PKs only after all skipped/unmatched rows have been removed.
    for igs_pk, row in enumerate(igs_rows, start=1):
        row["igs_pk"] = igs_pk

    return igs_rows, mapping_rows

def write_igs_tsv(
    output_path: Path,
    rows: Sequence[Mapping[str, object]],
) -> None:
    """Write the final igs.tsv file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=IGS_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()

        for row in rows:
            writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build records/sql_tsvs/igs.tsv from project GFF and FASTA files. "
            "Gene intervals are merged before intergenic coordinates are "
            "calculated, and sequence_pk is mapped from sequences.tsv."
        )
    )

    parser.add_argument(
        "genomes_dir",
        type=Path,
        help=(
            "Root genomes directory. Default input/output paths are resolved "
            "under this directory."
        ),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help=(
            "Project manifest CSV. Default: "
            "<genomes_dir>/records/project_manifests/"
            "body-size-dna-dynamics_manifest.csv"
        ),
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=None,
        help=(
            "Genomes metadata CSV. Default: "
            "<genomes_dir>/records/genomes_metadata.csv"
        ),
    )
    parser.add_argument(
        "--sequences",
        type=Path,
        default=None,
        help=(
            "Sequences reference TSV. Default: "
            "<genomes_dir>/records/sql_tsvs/sequences.tsv"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Output IGS TSV. Default: "
            "<genomes_dir>/records/sql_tsvs/igs.tsv"
        ),
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    genomes_dir = args.genomes_dir.expanduser().resolve()

    manifest_path = (
        args.manifest.expanduser()
        if args.manifest is not None
        else genomes_dir
        / "records"
        / "project_manifests"
        / "body-size-dna-dynamics_manifest.csv"
    )
    metadata_path = (
        args.metadata.expanduser()
        if args.metadata is not None
        else genomes_dir / "records" / "genomes_metadata.csv"
    )
    sequences_path = (
        args.sequences.expanduser()
        if args.sequences is not None
        else genomes_dir / "records" / "sql_tsvs" / "sequences.tsv"
    )
    output_path = (
        args.output.expanduser()
        if args.output is not None
        else genomes_dir / "records" / "sql_tsvs" / "igs.tsv"
    )

    for required_path, label in [
        (manifest_path, "manifest"),
        (metadata_path, "genomes metadata"),
        (sequences_path, "sequences TSV"),
    ]:
        if not required_path.is_file():
            print(
                f"ERROR: Required {label} file not found: {required_path}",
                file=sys.stderr,
            )
            return 1

    try:
        manifest_rows = read_delimited_rows(
            manifest_path,
            delimiter=",",
            required_columns=["accession"],
        )
        metadata_rows = read_delimited_rows(
            metadata_path,
            delimiter=",",
            required_columns=["accession", "path_to_gff", "path_to_fna"],
        )
        joined_genomes = left_join_manifest_metadata(
            manifest_rows,
            metadata_rows,
        )
        sequence_index = load_sequence_index(sequences_path)

        igs_rows, mapping_rows = build_igs_rows(
            joined_genomes,
            sequence_index,
            genomes_dir,
        )

        write_igs_tsv(output_path, igs_rows)

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    matched_mapping_rows = sum(
        1 for row in mapping_rows if clean(row.get("matched_sequence_pk"))
    )
    unmatched_mapping_rows = len(mapping_rows) - matched_mapping_rows

    print(f"Wrote {len(igs_rows):,} IGS rows to: {output_path}")
    print(
        "Sequence mapping summary: "
        f"{matched_mapping_rows:,} matched; "
        f"{unmatched_mapping_rows:,} unmatched/skipped."
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
