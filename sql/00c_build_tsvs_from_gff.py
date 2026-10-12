#!/usr/bin/env python3
"""
00c_build_tsvs_from_gff.py

Build gene and intergenic-sequence (IGS) relational TSVs from project GFF and
FASTA files, plus genome-level summary tables.

Outputs (default: <genomes_dir>/records/sql_tsvs/)
--------------------------------------------------
- genes.tsv
- igs.tsv
- gene_summary.tsv
- igs_summary.tsv

The script:
1. Reads the project manifest.
2. Left joins genomes_metadata.csv on accession to resolve path_to_gff and
   path_to_fna.
3. Streams gene features from each GFF using pure Python.
4. Maps GFF sequence IDs to sequence_pk values from sequences.tsv.
5. Writes original gene features to genes.tsv.
6. Merges overlapping gene bounds per sequence and calculates IGS intervals
   only BETWEEN consecutive merged gene intervals.
7. Uses pyfaidx to extract sequences and calculate GC as
       (G + C) / (A + T + G + C)
   with N/other characters excluded from the GC denominator.
8. Aggregates genome-level gene and IGS summaries.

Important behavior
------------------
- Missing/empty GFF paths do not terminate the project-wide run.
- Missing/empty FASTA paths do not prevent coordinate rows from being written;
  GC is left blank and summary recovery_status becomes "partial".
- GFF sequence IDs that cannot be mapped uniquely to sequence_pk are skipped
  from relational detail tables and recorded in recovery_notes.
- Gene coordinates are 1-based inclusive, as supplied by the GFF.
- IGS coordinates are 1-based inclusive.
- Adjacent genes are not merged. Because they have no bases between them, they
  simply do not create an IGS interval.
- Terminal regions before the first gene and after the last gene are not IGSs.

Example production run
----------------------
python 00c_build_tsvs_from_gff.py /path/to/genomes

Example test run with explicit outputs
--------------------------------------
python 00c_build_tsvs_from_gff.py /path/to/genomes \
    --manifest /path/to/test_manifest.csv \
    --metadata /path/to/test_genomes_metadata.csv \
    --sequences /path/to/test_sequences.tsv \
    --output /path/to/test_igs.tsv \
    --output-genes /path/to/test_genes.tsv \
    --output-gene-summary /path/to/test_gene_summary.tsv \
    --output-igs-summary /path/to/test_igs_summary.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import statistics
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
from urllib.parse import unquote


IGS_COLUMNS = [
    "igs_pk",
    "sequence_pk",
    "igs_seq_start",
    "igs_seq_end",
    "length",
    "gc",
]

GENE_COLUMNS = [
    "gene_pk",
    "sequence_pk",
    "gene_id",
    "g_strand",
    "g_seq_start",
    "g_seq_end",
    "length",
    "gc",
]

GENE_SUMMARY_COLUMNS = [
    "gene_summary_pk",
    "genome_pk",
    "species_pk",
    "n_genes",
    "callable_bp_total",
    "mean_gc",
    "weighted_mean_gc",
    "sd_gc",
    "var_gc",
    "median_gc",
    "mad_gc",
    "iqr_gc",
    "q05_gc",
    "q25_gc",
    "q75_gc",
    "q95_gc",
    "mean_gene_length",
    "median_gene_length",
    "recovery_status",
    "recovery_notes",
    "created_at",
]

IGS_SUMMARY_COLUMNS = [
    "igs_summary_pk",
    "genome_pk",
    "species_pk",
    "n_igs",
    "callable_bp_total",
    "mean_gc",
    "weighted_mean_gc",
    "sd_gc",
    "var_gc",
    "median_gc",
    "mad_gc",
    "iqr_gc",
    "q05_gc",
    "q25_gc",
    "q75_gc",
    "q95_gc",
    "mean_igs_length",
    "median_igs_length",
    "recovery_status",
    "recovery_notes",
    "created_at",
]


def warn(message: str) -> None:
    """Write a warning to stderr."""
    print(f"WARNING: {message}", file=sys.stderr)


def clean(value: object) -> str:
    """Return a stripped string, treating None as empty."""
    if value is None:
        return ""
    return str(value).strip()


def now_utc_iso() -> str:
    """Current UTC timestamp in ISO 8601 format."""
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace(
        "+00:00", "Z"
    )


def calculate_gc(sequence: str) -> Tuple[float, int, int, int, int, int, int]:
    """
    Calculate GC content using the supplied calculate_gc.py logic.

    GC = (G + C) / (A + T + G + C)
    N bases are counted separately and excluded from the denominator.
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
    """Open plain-text or gzip-compressed GFF text."""
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def require_columns(
    fieldnames: Optional[Sequence[str]],
    required: Iterable[str],
    path: Path,
) -> None:
    """Raise a clear error when required input columns are absent."""
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
    """Read a CSV/TSV file into dictionaries."""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        require_columns(reader.fieldnames, required_columns, path)
        return [
            {key: (value if value is not None else "") for key, value in row.items()}
            for row in reader
        ]


def resolve_project_path(raw_path: object, genomes_dir: Path) -> Optional[Path]:
    """Resolve metadata paths; blank values return None."""
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
    """Left join manifest to metadata on accession."""
    metadata_by_accession: Dict[str, Dict[str, str]] = {}

    for row in metadata_rows:
        accession = clean(row.get("accession"))
        if not accession:
            continue
        if accession in metadata_by_accession:
            warn(
                f"Duplicate accession {accession!r} in genomes metadata; "
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
            warn(f"Accession {accession!r} was not found in genomes_metadata.csv.")
        else:
            for key, value in metadata_row.items():
                if clean(value):
                    merged[key] = value
        joined.append(merged)

    return joined


def load_sequence_index(
    sequences_path: Path,
) -> Tuple[
    Dict[str, List[Dict[str, str]]],
    Dict[str, List[Dict[str, str]]],
]:
    """Index sequences.tsv by sequence_id and by genome_pk."""
    rows = read_delimited_rows(
        sequences_path,
        delimiter="\t",
        required_columns=["sequence_pk", "genome_pk", "sequence_id"],
    )

    by_sequence_id: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    by_genome_pk: Dict[str, List[Dict[str, str]]] = defaultdict(list)

    for row in rows:
        sequence_id = clean(row.get("sequence_id"))
        sequence_pk = clean(row.get("sequence_pk"))
        genome_pk = clean(row.get("genome_pk"))
        if not sequence_id or not sequence_pk:
            continue
        by_sequence_id[sequence_id].append(row)
        if genome_pk:
            by_genome_pk[genome_pk].append(row)

    return by_sequence_id, by_genome_pk


def load_optional_genomes_index(
    genomes_path: Path,
) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Dict[str, str]]]:
    """
    Load genomes.tsv when available.

    The DBML uses accession_id, while some starter files may use accession.
    Both are accepted for lookup.
    """
    if not genomes_path.is_file():
        return {}, {}

    rows = read_delimited_rows(
        genomes_path,
        delimiter="\t",
        required_columns=["genome_pk", "species_pk"],
    )

    by_accession: Dict[str, Dict[str, str]] = {}
    by_genome_pk: Dict[str, Dict[str, str]] = {}
    for row in rows:
        genome_pk = clean(row.get("genome_pk"))
        accession = clean(row.get("accession_id")) or clean(row.get("accession"))
        if genome_pk and genome_pk not in by_genome_pk:
            by_genome_pk[genome_pk] = dict(row)
        if accession and accession not in by_accession:
            by_accession[accession] = dict(row)

    return by_accession, by_genome_pk


def resolve_summary_ids(
    joined_row: Mapping[str, str],
    genomes_by_accession: Mapping[str, Mapping[str, str]],
) -> Tuple[str, str]:
    """Resolve genome_pk/species_pk from joined inputs or genomes.tsv."""
    accession = clean(joined_row.get("accession"))
    genome_pk = clean(joined_row.get("genome_pk"))
    species_pk = clean(joined_row.get("species_pk"))

    genome_row = genomes_by_accession.get(accession)
    if genome_row:
        genome_pk = clean(genome_row.get("genome_pk")) or genome_pk
        species_pk = clean(genome_row.get("species_pk")) or species_pk

    return genome_pk, species_pk


def resolve_sequence_record(
    chrom_id: str,
    sequence_index: Mapping[str, Sequence[Mapping[str, str]]],
    expected_genome_pk: str = "",
) -> Tuple[Optional[Mapping[str, str]], str]:
    """
    Resolve GFF chrom_id to one row of sequences.tsv.

    Exact sequence_id matching is always required. genome_pk is used only to
    disambiguate duplicate sequence IDs.
    """
    candidates = list(sequence_index.get(chrom_id, []))
    if not candidates:
        return None, "unmatched"

    if expected_genome_pk:
        narrowed = [
            row
            for row in candidates
            if clean(row.get("genome_pk")) == expected_genome_pk
        ]
        if len(narrowed) == 1:
            return narrowed[0], "matched"
        if len(narrowed) > 1:
            return None, "ambiguous"
        # If the supplied genome_pk does not match, do not silently map the
        # same sequence_id to another genome.
        return None, "genome_mismatch"

    if len(candidates) == 1:
        return candidates[0], "matched"

    return None, "ambiguous"


def parse_gff_attributes(attribute_text: str) -> Dict[str, str]:
    """Parse semicolon-delimited GFF3 column-9 attributes."""
    attributes: Dict[str, str] = {}
    for item in attribute_text.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key.strip()] = unquote(value.strip())
    return attributes


def extract_gene_id(attribute_text: str) -> str:
    """
    Extract a gene identifier from GFF column 9.

    Preferred field is ID=; fallbacks support common gene/gene_id/Name fields.
    The returned value excludes the attribute key prefix itself.
    """
    attrs = parse_gff_attributes(attribute_text)
    for key in ("ID", "gene", "gene_id", "Name"):
        value = clean(attrs.get(key))
        if value:
            return value
    return ""


def parse_gene_features(gff_path: Path) -> Dict[str, List[Dict[str, object]]]:
    """
    Stream gene features from GFF/GFF3 and group them by sequence ID.
    """
    genes_by_chrom: Dict[str, List[Dict[str, object]]] = defaultdict(list)

    with open_text_auto(gff_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            fields = line.rstrip("\r\n").split("\t")
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
                    f"{gff_path}: line {line_number} has coordinate < 1; skipping."
                )
                continue
            if start > end:
                warn(
                    f"{gff_path}: line {line_number} has start > end; "
                    "coordinates were normalized."
                )
                start, end = end, start

            gene_id = extract_gene_id(fields[8])
            if not gene_id:
                warn(
                    f"{gff_path}: line {line_number} gene has no usable ID; skipping."
                )
                continue

            strand = clean(fields[6])
            if strand not in {"+", "-"}:
                warn(
                    f"{gff_path}: line {line_number} gene {gene_id!r} has strand "
                    f"{strand!r}; retaining the value as supplied."
                )

            genes_by_chrom[chrom_id].append(
                {
                    "gene_id": gene_id,
                    "strand": strand,
                    "start": start,
                    "end": end,
                }
            )

    return genes_by_chrom


def fix_overlap(intervals: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Sort and merge overlapping gene intervals, preserving adjacency."""
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
    """Calculate 1-based inclusive IGSs between consecutive merged genes."""
    intervals: List[Tuple[int, int]] = []
    for left, right in zip(merged_gene_intervals, merged_gene_intervals[1:]):
        start = left[1] + 1
        end = right[0] - 1
        if start <= end:
            intervals.append((start, end))
    return intervals


def open_fasta(fna_path: Path):
    """Open a genome FASTA with pyfaidx."""
    try:
        from pyfaidx import Fasta
    except ImportError as exc:
        raise RuntimeError(
            "pyfaidx is required to calculate GC content. Install pyfaidx in "
            "the environment used to run this script."
        ) from exc
    return Fasta(str(fna_path))


def interval_gc(
    fasta,
    chrom_id: str,
    start_1based: int,
    end_1based: int,
) -> float:
    """Extract one 1-based inclusive interval and calculate GC."""
    sequence_length = len(fasta[chrom_id])
    if start_1based < 1 or end_1based > sequence_length or start_1based > end_1based:
        raise ValueError(
            f"interval {start_1based}-{end_1based} is outside FASTA sequence "
            f"length {sequence_length}"
        )

    sequence = str(fasta[chrom_id][start_1based - 1 : end_1based])
    gc_value, *_ = calculate_gc(sequence)
    return gc_value


def quantile_linear(values: Sequence[float], q: float) -> float:
    """NumPy-style linear quantile interpolation using standard Python."""
    if not values:
        raise ValueError("quantile requires at least one value")
    ordered = sorted(float(v) for v in values)
    if len(ordered) == 1:
        return ordered[0]

    position = (len(ordered) - 1) * q
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] + (ordered[upper] - ordered[lower]) * fraction


def summarize_intervals(
    rows: Sequence[Mapping[str, object]],
    length_field_name: str,
) -> Dict[str, object]:
    """
    Calculate summary statistics for interval rows.

    n_intervals counts all relational rows. GC statistics and
    callable_bp_total use rows with a successfully calculated numeric GC.
    callable_bp_total is therefore the denominator used for weighted_mean_gc,
    as requested.
    """
    lengths_all = [int(row["length"]) for row in rows]

    callable_rows: List[Tuple[int, float]] = []
    for row in rows:
        raw_gc = row.get("gc", "")
        if raw_gc in (None, ""):
            continue
        try:
            gc_value = float(raw_gc)
            length = int(row["length"])
        except (TypeError, ValueError):
            continue
        if math.isfinite(gc_value) and length > 0:
            callable_rows.append((length, gc_value))

    gc_values = [gc for _, gc in callable_rows]
    callable_bp_total = sum(length for length, _ in callable_rows)

    result: Dict[str, object] = {
        "n_intervals": len(rows),
        "callable_bp_total": callable_bp_total,
        "mean_gc": "",
        "weighted_mean_gc": "",
        "sd_gc": "",
        "var_gc": "",
        "median_gc": "",
        "mad_gc": "",
        "iqr_gc": "",
        "q05_gc": "",
        "q25_gc": "",
        "q75_gc": "",
        "q95_gc": "",
        length_field_name: statistics.mean(lengths_all) if lengths_all else "",
        f"median_{length_field_name[5:]}": (
            statistics.median(lengths_all) if lengths_all else ""
        ),
    }

    if gc_values:
        mean_gc = statistics.mean(gc_values)
        median_gc = statistics.median(gc_values)
        q05 = quantile_linear(gc_values, 0.05)
        q25 = quantile_linear(gc_values, 0.25)
        q75 = quantile_linear(gc_values, 0.75)
        q95 = quantile_linear(gc_values, 0.95)

        result.update(
            {
                "mean_gc": mean_gc,
                "weighted_mean_gc": (
                    sum(length * gc for length, gc in callable_rows)
                    / callable_bp_total
                    if callable_bp_total > 0
                    else ""
                ),
                "sd_gc": statistics.stdev(gc_values) if len(gc_values) > 1 else 0.0,
                "var_gc": statistics.variance(gc_values) if len(gc_values) > 1 else 0.0,
                "median_gc": median_gc,
                "mad_gc": statistics.median(
                    [abs(value - median_gc) for value in gc_values]
                ),
                "iqr_gc": q75 - q25,
                "q05_gc": q05,
                "q25_gc": q25,
                "q75_gc": q75,
                "q95_gc": q95,
            }
        )

    return result


def add_note(state: Dict[str, object], note: str) -> None:
    """Add a de-duplicated recovery note and mark the genome partial."""
    state["partial"] = True
    notes = state.setdefault("notes", [])
    if note not in notes:
        notes.append(note)


def build_all_rows(
    joined_genomes: Sequence[Mapping[str, str]],
    sequence_index: Mapping[str, Sequence[Mapping[str, str]]],
    genomes_by_accession: Mapping[str, Mapping[str, str]],
    genomes_by_genome_pk: Mapping[str, Mapping[str, str]],
    genomes_dir: Path,
) -> Tuple[
    List[Dict[str, object]],
    List[Dict[str, object]],
    Dict[str, Dict[str, object]],
]:
    """Build genes, IGSs, and per-accession recovery state in one GFF pass."""
    gene_rows: List[Dict[str, object]] = []
    igs_rows: List[Dict[str, object]] = []
    states: Dict[str, Dict[str, object]] = {}

    seen_gene_keys = set()
    seen_igs_keys = set()

    for joined_row in joined_genomes:
        accession = clean(joined_row.get("accession"))
        genome_pk, species_pk = resolve_summary_ids(
            joined_row, genomes_by_accession
        )

        state: Dict[str, object] = {
            "accession": accession,
            "genome_pk": genome_pk,
            "species_pk": species_pk,
            "partial": False,
            "notes": [],
        }
        states[accession] = state

        gff_path = resolve_project_path(joined_row.get("path_to_gff"), genomes_dir)
        fna_path = resolve_project_path(joined_row.get("path_to_fna"), genomes_dir)

        if gff_path is None:
            warn(f"{accession}: path_to_gff is empty; accession skipped.")
            add_note(state, "path_to_gff is empty")
            continue
        if not gff_path.is_file():
            warn(f"{accession}: GFF file not found at {gff_path}; accession skipped.")
            add_note(state, "GFF file not found")
            continue

        genes_by_chrom = parse_gene_features(gff_path)
        if not genes_by_chrom:
            warn(f"{accession}: no gene features were found in {gff_path}.")
            add_note(state, "no gene features found in GFF")

        fasta = None
        fasta_available = False
        if fna_path is None:
            warn(f"{accession}: path_to_fna is empty; GC values will be blank.")
            add_note(state, "path_to_fna is empty; GC unavailable")
        elif not fna_path.is_file():
            warn(
                f"{accession}: FASTA file not found at {fna_path}; "
                "GC values will be blank."
            )
            add_note(state, "FASTA file not found; GC unavailable")
        else:
            try:
                fasta = open_fasta(fna_path)
                fasta_available = True
            except Exception as exc:
                warn(
                    f"{accession}: could not open FASTA {fna_path}: {exc}. "
                    "GC values will be blank."
                )
                add_note(state, f"FASTA could not be opened: {exc}")

        try:
            for chrom_id in sorted(genes_by_chrom):
                sequence_record, match_status = resolve_sequence_record(
                    chrom_id,
                    sequence_index,
                    expected_genome_pk=clean(state.get("genome_pk")),
                )

                if sequence_record is None:
                    if match_status == "ambiguous":
                        message = (
                            f"chrom_id {chrom_id!r} has ambiguous sequence_id mapping"
                        )
                    elif match_status == "genome_mismatch":
                        message = (
                            f"chrom_id {chrom_id!r} matched sequence_id but not "
                            "the expected genome_pk"
                        )
                    else:
                        message = f"chrom_id {chrom_id!r} not found in sequences.tsv"
                    warn(f"{accession}: {message}; sequence skipped.")
                    add_note(state, message)
                    continue

                sequence_pk = clean(sequence_record.get("sequence_pk"))
                sequence_genome_pk = clean(sequence_record.get("genome_pk"))

                if not clean(state.get("genome_pk")) and sequence_genome_pk:
                    state["genome_pk"] = sequence_genome_pk
                    genome_row = genomes_by_genome_pk.get(sequence_genome_pk)
                    if genome_row and not clean(state.get("species_pk")):
                        state["species_pk"] = clean(genome_row.get("species_pk"))

                chrom_in_fasta = False
                if fasta_available and fasta is not None:
                    try:
                        chrom_in_fasta = chrom_id in fasta
                    except Exception:
                        chrom_in_fasta = False
                    if not chrom_in_fasta:
                        warn(
                            f"{accession}: chrom_id {chrom_id!r} is not present in "
                            "the FASTA; GC values will be blank for this sequence."
                        )
                        add_note(
                            state,
                            f"FASTA sequence {chrom_id!r} could not be resolved",
                        )

                chrom_genes = genes_by_chrom[chrom_id]

                # Original gene features are written without overlap merging.
                for gene in chrom_genes:
                    gene_id = clean(gene.get("gene_id"))
                    start = int(gene["start"])
                    end = int(gene["end"])
                    gene_key = (sequence_pk, gene_id)
                    if gene_key in seen_gene_keys:
                        warn(
                            f"{accession}: duplicate (sequence_pk, gene_id) "
                            f"{gene_key!r}; duplicate gene row skipped."
                        )
                        add_note(state, f"duplicate gene_id {gene_id!r} skipped")
                        continue
                    seen_gene_keys.add(gene_key)

                    gc_value: object = ""
                    if fasta_available and fasta is not None and chrom_in_fasta:
                        try:
                            gc_value = interval_gc(fasta, chrom_id, start, end)
                        except Exception as exc:
                            warn(
                                f"{accession}: failed GC for gene {gene_id!r} "
                                f"({chrom_id}:{start}-{end}): {exc}; GC left blank."
                            )
                            add_note(
                                state,
                                f"GC unavailable for one or more genes on {chrom_id!r}",
                            )

                    gene_rows.append(
                        {
                            "sequence_pk": int(sequence_pk),
                            "gene_id": gene_id,
                            "g_strand": clean(gene.get("strand")),
                            "g_seq_start": start,
                            "g_seq_end": end,
                            "length": end - start + 1,
                            "gc": gc_value,
                            "_accession": accession,
                            "_genome_pk": clean(state.get("genome_pk")),
                        }
                    )

                merged_genes = fix_overlap(
                    [(int(g["start"]), int(g["end"])) for g in chrom_genes]
                )
                for igs_start, igs_end in calculate_igs_intervals(merged_genes):
                    igs_key = (sequence_pk, igs_start, igs_end)
                    if igs_key in seen_igs_keys:
                        warn(
                            f"{accession}: duplicate IGS key {igs_key!r}; row skipped."
                        )
                        add_note(state, "duplicate IGS interval skipped")
                        continue
                    seen_igs_keys.add(igs_key)

                    gc_value = ""
                    if fasta_available and fasta is not None and chrom_in_fasta:
                        try:
                            gc_value = interval_gc(
                                fasta, chrom_id, igs_start, igs_end
                            )
                        except Exception as exc:
                            warn(
                                f"{accession}: failed GC for IGS "
                                f"{chrom_id}:{igs_start}-{igs_end}: {exc}; "
                                "GC left blank."
                            )
                            add_note(
                                state,
                                f"GC unavailable for one or more IGSs on {chrom_id!r}",
                            )

                    igs_rows.append(
                        {
                            "sequence_pk": int(sequence_pk),
                            "igs_seq_start": igs_start,
                            "igs_seq_end": igs_end,
                            "length": igs_end - igs_start + 1,
                            "gc": gc_value,
                            "_accession": accession,
                            "_genome_pk": clean(state.get("genome_pk")),
                        }
                    )
        finally:
            if fasta is not None and hasattr(fasta, "close"):
                fasta.close()

    # Some rows may have been created before genome_pk was inferred from a
    # sequence record. Fill internal grouping keys from final accession state.
    for row in gene_rows:
        if not clean(row.get("_genome_pk")):
            row["_genome_pk"] = clean(states[clean(row.get("_accession"))].get("genome_pk"))
    for row in igs_rows:
        if not clean(row.get("_genome_pk")):
            row["_genome_pk"] = clean(states[clean(row.get("_accession"))].get("genome_pk"))

    for gene_pk, row in enumerate(gene_rows, start=1):
        row["gene_pk"] = gene_pk
    for igs_pk, row in enumerate(igs_rows, start=1):
        row["igs_pk"] = igs_pk

    return gene_rows, igs_rows, states


def build_summary_rows(
    states: Mapping[str, Mapping[str, object]],
    detail_rows: Sequence[Mapping[str, object]],
    summary_type: str,
    created_at: str,
) -> List[Dict[str, object]]:
    """Build one genome-level summary row per resolvable accession."""
    if summary_type not in {"gene", "igs"}:
        raise ValueError("summary_type must be 'gene' or 'igs'")

    grouped: Dict[str, List[Mapping[str, object]]] = defaultdict(list)
    for row in detail_rows:
        accession = clean(row.get("_accession"))
        grouped[accession].append(row)

    output: List[Dict[str, object]] = []
    pk_name = "gene_summary_pk" if summary_type == "gene" else "igs_summary_pk"
    n_name = "n_genes" if summary_type == "gene" else "n_igs"
    mean_length_name = (
        "mean_gene_length" if summary_type == "gene" else "mean_igs_length"
    )
    median_length_name = (
        "median_gene_length" if summary_type == "gene" else "median_igs_length"
    )

    next_pk = 1
    for accession, state in states.items():
        genome_pk = clean(state.get("genome_pk"))
        species_pk = clean(state.get("species_pk"))

        if not genome_pk:
            warn(
                f"{accession}: genome_pk could not be resolved; "
                f"{summary_type}_summary row skipped."
            )
            continue

        notes = list(state.get("notes", []))
        partial = bool(state.get("partial", False))
        if not species_pk:
            partial = True
            note = "species_pk could not be resolved"
            if note not in notes:
                notes.append(note)

        rows = grouped.get(accession, [])
        stats = summarize_intervals(rows, mean_length_name)

        # Missing GC for any emitted interval is explicitly partial.
        n_callable = sum(1 for row in rows if row.get("gc", "") not in (None, ""))
        if n_callable < len(rows):
            partial = True
            note = f"GC unavailable for {len(rows) - n_callable} of {len(rows)} {summary_type} intervals"
            if note not in notes:
                notes.append(note)

        summary_row: Dict[str, object] = {
            pk_name: next_pk,
            "genome_pk": int(genome_pk) if genome_pk else "",
            "species_pk": int(species_pk) if species_pk else "",
            n_name: stats["n_intervals"],
            "callable_bp_total": stats["callable_bp_total"],
            "mean_gc": stats["mean_gc"],
            "weighted_mean_gc": stats["weighted_mean_gc"],
            "sd_gc": stats["sd_gc"],
            "var_gc": stats["var_gc"],
            "median_gc": stats["median_gc"],
            "mad_gc": stats["mad_gc"],
            "iqr_gc": stats["iqr_gc"],
            "q05_gc": stats["q05_gc"],
            "q25_gc": stats["q25_gc"],
            "q75_gc": stats["q75_gc"],
            "q95_gc": stats["q95_gc"],
            mean_length_name: stats[mean_length_name],
            median_length_name: stats[median_length_name],
            "recovery_status": "partial" if partial else "success",
            "recovery_notes": "; ".join(notes),
            "created_at": created_at,
        }
        output.append(summary_row)
        next_pk += 1

    return output


def write_tsv(
    output_path: Path,
    rows: Sequence[Mapping[str, object]],
    columns: Sequence[str],
) -> None:
    """Write a relational TSV with an exact ordered schema."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(columns),
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
            "Build genes.tsv, igs.tsv, gene_summary.tsv, and igs_summary.tsv "
            "from project GFF/FASTA files."
        )
    )

    parser.add_argument(
        "genomes_dir",
        type=Path,
        help="Root genomes directory used to resolve default project paths.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help=(
            "Project manifest CSV. Default: <genomes_dir>/records/"
            "project_manifests/body-size-dna-dynamics_manifest.csv"
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
        "--genomes",
        type=Path,
        default=None,
        help=(
            "Optional genomes.tsv used to resolve genome_pk/species_pk. "
            "Default: <genomes_dir>/records/sql_tsvs/genomes.tsv when present."
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
    parser.add_argument(
        "--output-genes",
        type=Path,
        default=None,
        help=(
            "Output genes TSV. Default: "
            "<genomes_dir>/records/sql_tsvs/genes.tsv"
        ),
    )
    parser.add_argument(
        "--output-gene-summary",
        type=Path,
        default=None,
        help=(
            "Output gene summary TSV. Default: "
            "<genomes_dir>/records/sql_tsvs/gene_summary.tsv"
        ),
    )
    parser.add_argument(
        "--output-igs-summary",
        type=Path,
        default=None,
        help=(
            "Output IGS summary TSV. Default: "
            "<genomes_dir>/records/sql_tsvs/igs_summary.tsv"
        ),
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()
    genomes_dir = args.genomes_dir.expanduser().resolve()
    sql_dir = genomes_dir / "records" / "sql_tsvs"

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
        else sql_dir / "sequences.tsv"
    )
    genomes_path = (
        args.genomes.expanduser()
        if args.genomes is not None
        else sql_dir / "genomes.tsv"
    )

    output_igs = args.output.expanduser() if args.output else sql_dir / "igs.tsv"
    output_genes = (
        args.output_genes.expanduser() if args.output_genes else sql_dir / "genes.tsv"
    )
    output_gene_summary = (
        args.output_gene_summary.expanduser()
        if args.output_gene_summary
        else sql_dir / "gene_summary.tsv"
    )
    output_igs_summary = (
        args.output_igs_summary.expanduser()
        if args.output_igs_summary
        else sql_dir / "igs_summary.tsv"
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
            manifest_path, delimiter=",", required_columns=["accession"]
        )
        metadata_rows = read_delimited_rows(
            metadata_path,
            delimiter=",",
            required_columns=["accession", "path_to_gff", "path_to_fna"],
        )
        joined_genomes = left_join_manifest_metadata(manifest_rows, metadata_rows)
        sequence_index, _ = load_sequence_index(sequences_path)
        genomes_by_accession, genomes_by_genome_pk = load_optional_genomes_index(
            genomes_path
        )

        gene_rows, igs_rows, states = build_all_rows(
            joined_genomes=joined_genomes,
            sequence_index=sequence_index,
            genomes_by_accession=genomes_by_accession,
            genomes_by_genome_pk=genomes_by_genome_pk,
            genomes_dir=genomes_dir,
        )

        created_at = now_utc_iso()
        gene_summary_rows = build_summary_rows(
            states, gene_rows, "gene", created_at
        )
        igs_summary_rows = build_summary_rows(
            states, igs_rows, "igs", created_at
        )

        write_tsv(output_genes, gene_rows, GENE_COLUMNS)
        write_tsv(output_igs, igs_rows, IGS_COLUMNS)
        write_tsv(output_gene_summary, gene_summary_rows, GENE_SUMMARY_COLUMNS)
        write_tsv(output_igs_summary, igs_summary_rows, IGS_SUMMARY_COLUMNS)

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote {len(gene_rows):,} gene rows to: {output_genes}")
    print(f"Wrote {len(igs_rows):,} IGS rows to: {output_igs}")
    print(
        f"Wrote {len(gene_summary_rows):,} gene summary rows to: "
        f"{output_gene_summary}"
    )
    print(
        f"Wrote {len(igs_summary_rows):,} IGS summary rows to: "
        f"{output_igs_summary}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
