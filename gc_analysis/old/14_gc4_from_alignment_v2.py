#!/usr/bin/env python3
"""14_gc4_from_alignment_v2.py

Calculate GC4 from codon-aware FASTA alignments using cumulative masks across
all sequences in each alignment. Output SQL-friendly detailed and per-species
TSVs with metadata identifying one nhPhyML analysis run.

Version 2 preserves the v1 biological masking and pooled-count calculations.
The only GC4 value change is that a zero denominator is written as \\N
(SQL NULL), rather than being reported as a GC4 fraction of zero.

Example:
  python 14_gc4_from_alignment_v2.py \
    --input-dir /path/to/alignments/80_2026oct14/01_clean_ortholog_alignments \
    --manifest /path/to/manifest.csv \
    --nhphyml-runs /path/to/records/sql_tsvs/nhphyml_runs.tsv \
    --output-detailed /path/to/alignments/80_2026oct14/11_gc4/detailed.tsv \
    --output-species-summary /path/to/alignments/80_2026oct14/11_gc4/species_summary.tsv
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np

FASTA_SUFFIXES = {".fa", ".faa", ".fasta", ".fas", ".fna"}
SQL_NULL = r"\N"
GC4_PREFIX_TO_FAMILY = {
    "GT": "Val", "CC": "Pro", "AC": "Thr", "GC": "Ala",
    "GG": "Gly", "CT": "Leu", "CG": "Arg", "TC": "Ser",
}

DETAIL_FIELDS = [
    "nhphyml_run_pk", "run_name", "alignment_name", "odb12_id", "species_key_short",
    "gc4", "g_gc4", "c_gc4", "a_gc4", "t_gc4", "n_gc4", "valid_gc4_bases",
    "retained_gc4_sites", "retained_gc4_codons", "n_sequences_in_alignment",
    "alignment_length_nt", "total_sites_nt", "complete_sites_nt",
    "dropped_sites_incomplete_nt", "total_codons", "complete_codons",
    "dropped_incomplete_codons", "candidate_gc4_codons_all_sequences",
    "dropped_not_gc4_all_sequences", "dropped_mixed_gc4_family", "status",
    "invalid_char_gap", "invalid_char_frameshift", "invalid_char_N",
]
SUMMARY_FIELDS = [
    "nhphyml_run_pk", "run_name", "alignment_name", "species_key_short",
    "gc4", "g_gc4", "c_gc4", "a_gc4", "t_gc4", "n_gc4", "valid_gc4_bases",
    "retained_gc4_sites", "retained_gc4_codons", "n_alignments",
    "alignments_with_retained_gc4", "alignment_length_nt", "total_sites_nt",
    "complete_sites_nt", "dropped_sites_incomplete_nt", "total_codons",
    "complete_codons", "dropped_incomplete_codons",
    "candidate_gc4_codons_all_sequences", "dropped_not_gc4_all_sequences",
    "dropped_mixed_gc4_family", "invalid_char_gap", "invalid_char_frameshift",
    "invalid_char_N",
]
SUM_FIELDS = [
    "g_gc4", "c_gc4", "a_gc4", "t_gc4", "n_gc4", "valid_gc4_bases",
    "retained_gc4_sites", "retained_gc4_codons", "alignment_length_nt",
    "total_sites_nt", "complete_sites_nt", "dropped_sites_incomplete_nt",
    "total_codons", "complete_codons", "dropped_incomplete_codons",
    "candidate_gc4_codons_all_sequences", "dropped_not_gc4_all_sequences",
    "dropped_mixed_gc4_family", "invalid_char_gap", "invalid_char_frameshift",
    "invalid_char_N",
]

# -----------------------------------------------------------------------------
# FASTA parsing, cumulative masks and counts: unchanged from v1
# -----------------------------------------------------------------------------

def parse_fasta(fasta_path: Path) -> Dict[str, str]:
    """Parse an aligned FASTA file into {record_id: aligned_sequence}."""
    records: Dict[str, List[str]] = {}
    current_id: str | None = None

    with fasta_path.open("r") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                current_id = line[1:].strip().split()[0]

                if not current_id:
                    raise ValueError(
                        f"{fasta_path}: empty FASTA header on line {line_number}"
                    )

                if current_id in records:
                    raise ValueError(
                        f"{fasta_path}: duplicate FASTA ID: {current_id}"
                    )

                records[current_id] = []
                continue

            if current_id is None:
                raise ValueError(
                    f"{fasta_path}: sequence data before first FASTA header "
                    f"on line {line_number}"
                )

            records[current_id].append(line)

    if not records:
        raise ValueError(f"{fasta_path}: no FASTA records found")

    return {
        record_id: "".join(seq_lines).upper()
        for record_id, seq_lines in records.items()
    }


def validate_alignment(records: Dict[str, str], alignment_name: str = "alignment") -> int:
    """Confirm all sequences have same length and length is divisible by 3."""
    lengths = {record_id: len(seq) for record_id, seq in records.items()}
    unique_lengths = set(lengths.values())

    if len(unique_lengths) != 1:
        examples = ", ".join(
            f"{record_id}={length}"
            for record_id, length in list(lengths.items())[:10]
        )
        raise ValueError(
            f"{alignment_name}: aligned sequences do not all have the same "
            f"length. First lengths: {examples}"
        )

    alignment_length = next(iter(unique_lengths))

    if alignment_length % 3 != 0:
        raise ValueError(
            f"{alignment_name}: alignment length ({alignment_length}) is not "
            "divisible by 3. This script expects codon-aware alignments."
        )

    return alignment_length


# -----------------------------------------------------------------------------
# Mask creation
# -----------------------------------------------------------------------------


def make_complete_site_mask(
    records: Dict[str, str],
    invalid_chars: str = "-!Nn",
    alignment_name: str = "alignment",
) -> np.ndarray:
    """
    Build a cumulative nucleotide/site mask.

    Starts as all True. For each sequence, positions containing invalid_chars
    are set to False. Once False, a position stays False.
    """
    alignment_length = validate_alignment(records, alignment_name=alignment_name)
    complete_site_mask = np.ones(alignment_length, dtype=bool)

    invalid_byte_values = np.frombuffer(
        invalid_chars.encode("ascii"),
        dtype=np.uint8,
    )

    for sequence in records.values():
        seq_array = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
        invalid_positions = np.isin(seq_array, invalid_byte_values)
        complete_site_mask[invalid_positions] = False

    return complete_site_mask


def make_complete_codon_mask(site_mask: np.ndarray) -> np.ndarray:
    """
    Convert a nucleotide/site mask into a codon mask.

    If any of the three nucleotide positions in a codon is False, all three
    positions for that codon are set to False.
    """
    if len(site_mask) % 3 != 0:
        raise ValueError(f"Mask length ({len(site_mask)}) is not divisible by 3")

    codon_triplets = site_mask.reshape(-1, 3)
    complete_codons = np.all(codon_triplets, axis=1)
    return np.repeat(complete_codons, 3)


def codon_is_gc4(codon: str) -> bool:
    """Return True if codon belongs to one of the 8 GC4 codon families."""
    return codon[:2].upper() in GC4_PREFIX_TO_FAMILY


def gc4_family(codon: str) -> str:
    """Return the GC4 family name for a codon, or empty string if not GC4."""
    return GC4_PREFIX_TO_FAMILY.get(codon[:2].upper(), "")


def make_cumulative_gc4_mask(
    records: Dict[str, str],
    complete_codon_mask: np.ndarray,
    require_same_gc4_family: bool = False,
    alignment_name: str = "alignment",
) -> Tuple[np.ndarray, Dict[str, int]]:
    """
    Build a cumulative GC4 mask across all sequences.

    A codon is retained only when:
        1. The codon is complete across all species according to complete_codon_mask.
        2. Every sequence has a GC4-family codon at that codon position.

    If require_same_gc4_family=True, all sequences must also belong to the same
    GC4 codon family at that codon position.
    """
    alignment_length = validate_alignment(records, alignment_name=alignment_name)
    n_codons = alignment_length // 3

    complete_codons = complete_codon_mask.reshape(-1, 3).all(axis=1)
    final_gc4_site_mask = np.zeros(alignment_length, dtype=bool)

    complete_codon_count = 0
    dropped_incomplete_codon_count = 0
    candidate_gc4_all_sequences_count = 0
    dropped_not_gc4_all_sequences_count = 0
    dropped_mixed_gc4_family_count = 0
    retained_gc4_codon_count = 0

    sequences = list(records.values())

    for codon_index in range(n_codons):
        start = codon_index * 3
        end = start + 3

        if not complete_codons[codon_index]:
            dropped_incomplete_codon_count += 1
            continue

        complete_codon_count += 1
        codons = [sequence[start:end] for sequence in sequences]

        if not all(codon_is_gc4(codon) for codon in codons):
            dropped_not_gc4_all_sequences_count += 1
            continue

        candidate_gc4_all_sequences_count += 1

        if require_same_gc4_family:
            families = {gc4_family(codon) for codon in codons}
            if len(families) != 1:
                dropped_mixed_gc4_family_count += 1
                continue

        final_gc4_site_mask[start + 2] = True
        retained_gc4_codon_count += 1

    reason_counts = {
        "total_codons": n_codons,
        "complete_codons": complete_codon_count,
        "dropped_incomplete_codons": dropped_incomplete_codon_count,
        "candidate_gc4_codons_all_sequences": candidate_gc4_all_sequences_count,
        "dropped_not_gc4_all_sequences": dropped_not_gc4_all_sequences_count,
        "dropped_mixed_gc4_family": dropped_mixed_gc4_family_count,
        "retained_gc4_codons": retained_gc4_codon_count,
    }

    return final_gc4_site_mask, reason_counts


# -----------------------------------------------------------------------------
# GC calculation and summaries
# -----------------------------------------------------------------------------


def calculate_gc(sequence: str) -> Tuple[float, int, int, int, int, int, int]:
    """Calculate GC content as (G + C) / (A + T + G + C)."""
    sequence = str(sequence).upper()

    g_count = sequence.count("G")
    c_count = sequence.count("C")
    a_count = sequence.count("A")
    t_count = sequence.count("T")
    n_count = sequence.count("N")
    total_valid_bases = a_count + t_count + g_count + c_count
    gc_content = (g_count + c_count) / total_valid_bases if total_valid_bases else 0.0

    return gc_content, g_count, c_count, a_count, t_count, n_count, total_valid_bases


def extract_masked_sequence(sequence: str, mask: np.ndarray) -> str:
    """Extract characters from sequence where mask is True."""
    if len(sequence) != len(mask):
        raise ValueError(
            f"Sequence length ({len(sequence)}) does not match mask length ({len(mask)})"
        )

    seq_array = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
    return seq_array[mask].tobytes().decode("ascii")


def summarize_invalid_sites_by_character(records: Dict[str, str], invalid_chars: str) -> Dict[str, int]:
    """
    Count sequence-level invalid characters.

    These are character counts across all sequences, not unique alignment-column
    counts. Because parse_fasta() uppercases sequences, N/n are summarized as N.
    """
    labels: Dict[str, str] = {}
    for char in invalid_chars:
        normalized = char.upper()
        if normalized not in labels:
            if normalized == "-":
                labels[normalized] = "invalid_char_gap"
            elif normalized == "!":
                labels[normalized] = "invalid_char_frameshift"
            else:
                labels[normalized] = f"invalid_char_{normalized}"

    counts = {label: 0 for label in labels.values()}

    for sequence in records.values():
        for normalized, label in labels.items():
            counts[label] += sequence.count(normalized)

    return counts


# -----------------------------------------------------------------------------
# Manifest, run metadata and paths
# -----------------------------------------------------------------------------


def short_species_key(full_key: str, manifest_path: Path, row_number: int) -> str:
    """Derive an exact Genus_species key from the manifest's species_key."""
    components = full_key.split("_")
    if (len(components) < 2 or not components[0] or not components[1]
            or any(not item for item in components)):
        raise ValueError(
            f"{manifest_path}: row {row_number}: invalid species_key {full_key!r}; "
            "expected Genus_species or Genus_species_subspecies."
        )
    if any(any(character.isspace() for character in part) for part in components):
        raise ValueError(
            f"{manifest_path}: row {row_number}: species_key {full_key!r} "
            "contains whitespace."
        )
    return "_".join(components[:2])


def load_manifest(manifest_path: Path) -> dict[str, list[tuple[int, str]]]:
    """Index manifest rows by shortened species key; keep duplicates for QC."""
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Project manifest not found: {manifest_path}")

    lookup: dict[str, list[tuple[int, str]]] = defaultdict(list)
    with manifest_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "species_key" not in reader.fieldnames:
            raise ValueError(f"{manifest_path}: missing required column 'species_key'.")
        for row_number, row in enumerate(reader, start=2):
            value = (row.get("species_key") or "").strip()
            if not value:
                raise ValueError(f"{manifest_path}: row {row_number}: empty species_key.")
            shorter = short_species_key(value, manifest_path, row_number)
            lookup[shorter].append((row_number, value))
    if not lookup:
        raise ValueError(f"{manifest_path}: manifest contains no species records.")
    return dict(lookup)


def validate_fasta_species(
    records: Dict[str, str], manifest_lookup: dict[str, list[tuple[int, str]]],
    manifest_path: Path, fasta_path: Path,
) -> None:
    """Require each FASTA identifier to correspond to exactly one manifest row."""
    missing = sorted(species for species in records if species not in manifest_lookup)
    ambiguous = {
        species: manifest_lookup[species]
        for species in sorted(records)
        if len(manifest_lookup.get(species, [])) > 1
    }
    issues: list[str] = []
    if missing:
        issues.append("FASTA species_key_short absent from manifest: " + ", ".join(missing))
    if ambiguous:
        for short_key, entries in ambiguous.items():
            conflicts = ", ".join(f"row {row}: {full}" for row, full in entries)
            issues.append(f"Ambiguous species_key_short {short_key!r}: {conflicts}")
    if issues:
        raise ValueError(
            f"{fasta_path}: manifest validation failed ({manifest_path}):\n    "
            + "\n    ".join(issues)
        )


def path_alignment_group(path: Path) -> str | None:
    """Find 'records/compleasm/alignments/<analysis-group>' in a path."""
    parts = path.resolve().parts
    for idx in range(len(parts) - 3):
        if tuple(parts[idx:idx + 3]) == ("records", "compleasm", "alignments"):
            return parts[idx + 3]
    return None


def select_alignment_name(
    manual_name: str | None, input_path: Path, detailed_path: Path,
    summary_path: Path,
) -> tuple[str, str]:
    """Determine one run group from known paths unless manually overridden."""
    if manual_name is not None:
        manual_name = manual_name.strip()
        if not manual_name:
            raise ValueError("--alignment-name must be nonempty.")
        return manual_name, "Manually specified (--alignment-name)"

    detected: dict[str, str] = {}
    for label, path in (
        ("input", input_path),
        ("detailed output", detailed_path),
        ("species summary output", summary_path),
    ):
        name = path_alignment_group(path)
        if name is not None:
            detected[label] = name
    groups = set(detected.values())
    if not groups:
        raise ValueError(
            "Unable to detect alignment_name from input/output directories. "
            "Expected records/compleasm/alignments/<alignment_name>/... . "
            "Supply --alignment-name manually."
        )
    if len(groups) != 1:
        raise ValueError(
            "Conflicting alignment groups detected from paths: "
            + "; ".join(f"{label}={value}" for label, value in detected.items())
            + ". Supply --alignment-name to override automatic detection."
        )
    return next(iter(groups)), "Automatically detected"


def infer_nhphyml_runs(*paths: Path) -> Path:
    """Infer genomes/records/sql_tsvs/nhphyml_runs.tsv from a standard path."""
    candidates: set[Path] = set()
    for path in paths:
        parts = path.resolve().parts
        for idx in range(len(parts) - 3):
            if tuple(parts[idx:idx + 3]) == ("records", "compleasm", "alignments"):
                records_path = Path(*parts[:idx + 1])
                candidates.add(records_path / "sql_tsvs" / "nhphyml_runs.tsv")
                break
    if len(candidates) == 1:
        candidate = next(iter(candidates))
        if candidate.is_file():
            return candidate
        raise FileNotFoundError(
            f"Inferred nhphyml_runs.tsv is missing: {candidate}. "
            "Supply --nhphyml-runs explicitly."
        )
    if candidates:
        raise ValueError(
            "Input/output paths imply multiple genomes roots. Supply --nhphyml-runs "
            "explicitly: " + ", ".join(str(path) for path in sorted(candidates))
        )
    raise ValueError(
        "Cannot infer nhphyml_runs.tsv without a standard directory path. "
        "Supply --nhphyml-runs explicitly."
    )


def load_nhphyml_run(tsv_path: Path, alignment_name: str) -> dict[str, object]:
    """Retrieve exactly one matching run, including optional expected gene count."""
    if not tsv_path.is_file():
        raise FileNotFoundError(f"nhPhyML runs TSV not found: {tsv_path}")
    matches: list[tuple[int, dict[str, str]]] = []
    with tsv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"nhphyml_run_pk", "run_name", "alignment_name"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{tsv_path}: required columns missing: {', '.join(sorted(missing))}.")
        for row_number, row in enumerate(reader, start=2):
            # Preserve all other rows; only the selected run is relevant here.
            if (row.get("alignment_name") or "").strip() == alignment_name:
                matches.append((row_number, row))
    if not matches:
        raise ValueError(
            f"No row with alignment_name={alignment_name!r} in {tsv_path}."
        )
    if len(matches) > 1:
        detail = ", ".join(
            f"row {number}: nhphyml_run_pk={row.get('nhphyml_run_pk')!r}"
            for number, row in matches
        )
        raise ValueError(
            f"Duplicate alignment_name={alignment_name!r} in {tsv_path}: {detail}. "
            "Exactly one matching run is required."
        )
    line, row = matches[0]
    pk_text = (row.get("nhphyml_run_pk") or "").strip()
    try:
        pk = int(pk_text)
    except ValueError as exc:
        raise ValueError(
            f"{tsv_path}: row {line}: invalid nhphyml_run_pk={pk_text!r}."
        ) from exc
    if pk < 1:
        raise ValueError(f"{tsv_path}: row {line}: nhphyml_run_pk must be positive.")
    run_name = (row.get("run_name") or "").strip()
    if not run_name:
        raise ValueError(f"{tsv_path}: row {line}: run_name cannot be empty.")
    n_genes: int | None = None
    raw_genes = (row.get("n_genes") or "").strip()
    if raw_genes:
        try:
            n_genes = int(raw_genes)
            if n_genes < 0:
                raise ValueError
        except ValueError:
            print(
                f"[WARNING] Invalid optional n_genes={raw_genes!r} at "
                f"{tsv_path} row {line}; skipping gene-count comparison.",
                file=sys.stderr,
            )
            n_genes = None
    return {
        "nhphyml_run_pk": pk,
        "run_name": run_name,
        "alignment_name": alignment_name,
        "n_genes": n_genes,
    }


# -----------------------------------------------------------------------------
# Calculation and summaries: same cumulative masks and pooled GC4 as v1
# -----------------------------------------------------------------------------


def calculate_gc4_for_alignment(
    fasta_path: Path,
    odb12_id: str,
    manifest_lookup: dict[str, list[tuple[int, str]]],
    manifest_path: Path,
    run_metadata: dict[str, object],
    invalid_chars: str = "-!Nn",
    require_same_gc4_family: bool = False,
) -> List[Dict[str, object]]:
    """Calculate per-species GC4 for one ortholog using v1 cumulative masks."""
    records = parse_fasta(fasta_path)
    validate_fasta_species(records, manifest_lookup, manifest_path, fasta_path)
    alignment_file_name = fasta_path.name
    alignment_length = validate_alignment(records, alignment_name=alignment_file_name)
    n_sequences = len(records)

    site_mask = make_complete_site_mask(
        records=records,
        invalid_chars=invalid_chars,
        alignment_name=alignment_file_name,
    )
    codon_mask = make_complete_codon_mask(site_mask)
    gc4_site_mask, codon_reason_counts = make_cumulative_gc4_mask(
        records=records,
        complete_codon_mask=codon_mask,
        require_same_gc4_family=require_same_gc4_family,
        alignment_name=alignment_file_name,
    )

    total_sites = alignment_length
    complete_sites = int(site_mask.sum())
    dropped_sites_incomplete = total_sites - complete_sites

    total_codons = alignment_length // 3
    complete_codons = int(codon_mask.sum() // 3)
    dropped_incomplete_codons = total_codons - complete_codons

    retained_gc4_sites = int(gc4_site_mask.sum())
    retained_gc4_codons = retained_gc4_sites

    # As in v1, these counts cover *all* sequences, then repeat in detailed
    # rows. They are not per-sequence counts or distinct alignment columns.
    invalid_counts = summarize_invalid_sites_by_character(records, invalid_chars)
    # The SQL schema fixes these three columns, even for a custom invalid mask.
    fixed_invalid_counts = {
        col: invalid_counts.get(col, 0)
        for col in ("invalid_char_gap", "invalid_char_frameshift", "invalid_char_N")
    }

    rows: List[Dict[str, object]] = []
    for species_key_short, sequence in records.items():
        gc4_sequence = extract_masked_sequence(sequence, gc4_site_mask)
        (
            gc4_content, g_count, c_count, a_count, t_count, n_count,
            total_valid_bases,
        ) = calculate_gc(gc4_sequence)
        rows.append({
            "nhphyml_run_pk": run_metadata["nhphyml_run_pk"],
            "run_name": run_metadata["run_name"],
            "alignment_name": run_metadata["alignment_name"],
            "odb12_id": odb12_id,
            "species_key_short": species_key_short,
            "gc4": gc4_content if total_valid_bases else SQL_NULL,
            "g_gc4": g_count,
            "c_gc4": c_count,
            "a_gc4": a_count,
            "t_gc4": t_count,
            "n_gc4": n_count,
            "valid_gc4_bases": total_valid_bases,
            "retained_gc4_sites": retained_gc4_sites,
            "retained_gc4_codons": retained_gc4_codons,
            "n_sequences_in_alignment": n_sequences,
            "alignment_length_nt": alignment_length,
            "total_sites_nt": total_sites,
            "complete_sites_nt": complete_sites,
            "dropped_sites_incomplete_nt": dropped_sites_incomplete,
            "total_codons": total_codons,
            "complete_codons": complete_codons,
            "dropped_incomplete_codons": dropped_incomplete_codons,
            "candidate_gc4_codons_all_sequences": codon_reason_counts[
                "candidate_gc4_codons_all_sequences"
            ],
            "dropped_not_gc4_all_sequences": codon_reason_counts[
                "dropped_not_gc4_all_sequences"
            ],
            "dropped_mixed_gc4_family": codon_reason_counts[
                "dropped_mixed_gc4_family"
            ],
            "status": "ok",
            **fixed_invalid_counts,
        })
    return rows


def summarize_by_species(detail_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    """Group by (run PK, species); compute pooled (G+C)/valid base count."""
    grouped: dict[tuple[int, str], dict[str, object]] = {}
    for row in detail_rows:
        key = (int(row["nhphyml_run_pk"]), str(row["species_key_short"]))
        if key not in grouped:
            grouped[key] = {
                "nhphyml_run_pk": row["nhphyml_run_pk"],
                "run_name": row["run_name"],
                "alignment_name": row["alignment_name"],
                "species_key_short": row["species_key_short"],
                "n_alignments": 0,
                "alignments_with_retained_gc4": 0,
                **{field: 0 for field in SUM_FIELDS},
            }
        out = grouped[key]
        out["n_alignments"] = int(out["n_alignments"]) + 1
        if int(row["valid_gc4_bases"]) > 0:
            out["alignments_with_retained_gc4"] = (
                int(out["alignments_with_retained_gc4"]) + 1
            )
        for field in SUM_FIELDS:
            out[field] = int(out[field]) + int(row[field])

    summary_rows: List[Dict[str, object]] = []
    for key in sorted(grouped):
        out = grouped[key]
        valid = int(out["valid_gc4_bases"])
        # Identical to v1 when valid > 0. Never average per-ortholog GC4.
        out["gc4"] = (
            (int(out["g_gc4"]) + int(out["c_gc4"])) / valid
            if valid else SQL_NULL
        )
        summary_rows.append(out)
    return summary_rows


# -----------------------------------------------------------------------------
# Files and paired output writing
# -----------------------------------------------------------------------------


def collect_input_fastas(input_path: Path | None, input_dir: Path | None) -> list[Path]:
    if input_path is not None:
        if not input_path.is_file():
            raise FileNotFoundError(f"Input FASTA does not exist: {input_path}")
        if input_path.suffix.lower() not in FASTA_SUFFIXES:
            raise ValueError(
                f"Unsupported FASTA extension for {input_path}; expected "
                + ", ".join(sorted(FASTA_SUFFIXES))
            )
        return [input_path]
    if input_dir is None:
        raise ValueError("Supply --input or --input-dir.")
    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input directory not found: {input_dir}")
    fasta_files = [
        path for path in sorted(input_dir.iterdir())
        if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES
    ]
    if not fasta_files:
        raise ValueError(f"No supported FASTA alignment files in {input_dir}")
    return fasta_files


def ortholog_ids(fasta_files: list[Path]) -> dict[Path, str]:
    result: dict[Path, str] = {}
    previous: dict[str, Path] = {}
    for path in fasta_files:
        identifier = path.stem
        if not identifier.strip():
            raise ValueError(f"Empty odb12_id extracted from {path}")
        if identifier in previous:
            raise ValueError(
                f"Duplicate odb12_id {identifier!r} from {previous[identifier]} "
                f"and {path}. The input directory must contain one file per ortholog."
            )
        previous[identifier] = path
        result[path] = identifier
    return result


def write_tsv(rows: list[dict[str, object]], path: Path, columns: list[str]) -> None:
    """Write one staged TSV with an explicit, stable SQL-oriented schema."""
    if not rows:
        raise ValueError(f"No rows available for {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", extrasaction="raise",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
        handle.flush()
        os.fsync(handle.fileno())


def reserve_temp(parent: Path, label: str) -> Path:
    descriptor, name = tempfile.mkstemp(prefix=f".gc4_{label}_", suffix=".tmp", dir=parent)
    os.close(descriptor)
    return Path(name)


def write_paired_tsvs(
    detail_rows: list[dict[str, object]], summary_rows: list[dict[str, object]],
    detail_path: Path, summary_path: Path,
) -> None:
    """Stage both files, snapshot originals, and roll back if a swap fails.

    This is a best-effort paired file transaction, not a filesystem-level atomic
    transaction. Staging and backup files live alongside their final outputs.
    """
    targets = [(detail_path, detail_rows, DETAIL_FIELDS),
               (summary_path, summary_rows, SUMMARY_FIELDS)]
    for target, _, _ in targets:
        target.parent.mkdir(parents=True, exist_ok=True)
    staged: dict[Path, Path] = {}
    backups: dict[Path, Path] = {}
    replaced: list[Path] = []
    preserve_backups = False
    try:
        for target, rows, fields in targets:
            tmp = reserve_temp(target.parent, "stage")
            staged[target] = tmp
            write_tsv(rows, tmp, fields)
        # No existing result changes before both files are staged successfully.
        for target, _, _ in targets:
            if target.exists():
                if not target.is_file():
                    raise ValueError(f"Output target exists but is not a file: {target}")
                backup = reserve_temp(target.parent, "backup")
                backups[target] = backup
                shutil.copy2(target, backup)
        try:
            for target, _, _ in targets:
                os.replace(staged[target], target)
                replaced.append(target)
        except Exception as original_error:
            rollback_errors = []
            for target in reversed(replaced):
                try:
                    if target in backups:
                        os.replace(backups[target], target)
                    else:
                        target.unlink(missing_ok=True)
                except OSError as rollback_error:
                    rollback_errors.append(f"{target}: {rollback_error}")
            if rollback_errors:
                # Never delete a backup that might be needed for manual recovery.
                preserve_backups = True
                backup_locations = ", ".join(str(path) for path in backups.values())
                raise RuntimeError(
                    f"Failed installing paired output: {original_error}. "
                    "Rollback also failed: " + "; ".join(rollback_errors)
                    + f". Original-file backups retained at: {backup_locations}"
                ) from original_error
            raise
    finally:
        for temp_file in staged.values():
            temp_file.unlink(missing_ok=True)
        if not preserve_backups:
            for temp_file in backups.values():
                temp_file.unlink(missing_ok=True)


def check_output_paths(
    detail_path: Path, summary_path: Path, protected_files: Iterable[Path],
) -> None:
    detailed = detail_path.resolve()
    summary = summary_path.resolve()
    if detailed == summary:
        raise ValueError("Detailed and species-summary output paths must be different.")
    protected = {item.resolve() for item in protected_files}
    for target in (detailed, summary):
        if target in protected:
            raise ValueError(
                f"Output path {target} overlaps an input FASTA or metadata file. "
                "Input files must not be overwritten."
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="GC4 from cumulative codon-aware masks, with validated species "
                    "and nhPhyML run identifiers for SQL TSV output."
    )
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--input", type=Path, help="Single codon-aware aligned FASTA.")
    inputs.add_argument("--input-dir", type=Path, help="Directory with aligned FASTAs.")
    parser.add_argument("--manifest", type=Path, required=True,
                        help="CSV manifest with species_key column.")
    parser.add_argument("--nhphyml-runs", type=Path,
                        help="TSV with nhphyml_run_pk, run_name and alignment_name. "
                             "Inferred from genomes/records/sql_tsvs when possible.")
    parser.add_argument("--alignment-name", type=str,
                        help="Override group name detection, e.g. 80_2026oct14.")
    parser.add_argument("--output-detailed", type=Path, required=True,
                        help="Destination for one row per ortholog and species.")
    parser.add_argument("--output-species-summary", type=Path, required=True,
                        help="Destination for per-species pooled GC4 summary.")
    parser.add_argument("--invalid-chars", default="-!Nn",
                        help="Alignment characters invalidating a site (default: -!Nn).")
    parser.add_argument("--require-same-gc4-family", action="store_true",
                        help="Require same qualifying fourfold-degenerate codon family "
                             "in all sequences at each retained codon.")
    args = parser.parse_args()

    try:
        fasta_files = collect_input_fastas(args.input, args.input_dir)
        ids = ortholog_ids(fasta_files)
        source = args.input if args.input is not None else args.input_dir
        assert source is not None
        alignment_name, detection_source = select_alignment_name(
            args.alignment_name, source, args.output_detailed, args.output_species_summary,
        )
        nhphyml_path = args.nhphyml_runs or infer_nhphyml_runs(
            source, args.output_detailed, args.output_species_summary,
        )
        check_output_paths(
            args.output_detailed, args.output_species_summary,
            [*fasta_files, args.manifest, nhphyml_path],
        )
        manifest_lookup = load_manifest(args.manifest)
        run_metadata = load_nhphyml_run(nhphyml_path, alignment_name)

        custom_chars = {c.upper() for c in args.invalid_chars} - {"-", "!", "N"}
        if custom_chars:
            print(
                "[WARNING] Additional --invalid-chars are used for masking, but "
                "the fixed SQL output schema reports only gap (-), frameshift (!) "
                "and N counts: " + ", ".join(sorted(custom_chars)), file=sys.stderr,
            )

        detailed_rows: list[dict[str, object]] = []
        failed_alignments: list[tuple[Path, str]] = []
        alignments_with_gc4 = 0
        for fasta_path in fasta_files:
            try:
                rows = calculate_gc4_for_alignment(
                    fasta_path=fasta_path,
                    odb12_id=ids[fasta_path],
                    manifest_lookup=manifest_lookup,
                    manifest_path=args.manifest,
                    run_metadata=run_metadata,
                    invalid_chars=args.invalid_chars,
                    require_same_gc4_family=args.require_same_gc4_family,
                )
                detailed_rows.extend(rows)
                if rows and int(rows[0]["retained_gc4_sites"]) > 0:
                    alignments_with_gc4 += 1
            except Exception as exc:
                failed_alignments.append((fasta_path, str(exc)))

        if failed_alignments:
            failures = "\n".join(
                f"  {path}: {problem}" for path, problem in failed_alignments
            )
            raise RuntimeError(
                f"{len(failed_alignments)} alignment(s) failed; no output was written:\n"
                + failures
            )

        detailed_rows.sort(key=lambda row: (str(row["odb12_id"]), str(row["species_key_short"])))
        summary_rows = summarize_by_species(detailed_rows)
        zero_details = sum(int(row["valid_gc4_bases"]) == 0 for row in detailed_rows)
        zero_species = sum(int(row["valid_gc4_bases"]) == 0 for row in summary_rows)
        actual_genes = len(ids)
        expected_genes = run_metadata.get("n_genes")
        if expected_genes is not None and actual_genes != expected_genes:
            print(
                f"[WARNING] {actual_genes} ortholog alignments processed, but "
                f"nhphyml_runs.tsv lists n_genes={expected_genes} for "
                f"{alignment_name}. This may be intentional for a subset.",
                file=sys.stderr,
            )

        write_paired_tsvs(
            detailed_rows, summary_rows, args.output_detailed, args.output_species_summary
        )

        print("14_gc4_from_alignment_v2.py")
        print("=" * 60)
        print(f"Manifest:                  {args.manifest}")
        print(f"nhPhyML runs TSV:          {nhphyml_path}")
        print(f"Alignment name:            {alignment_name}")
        print(f"Alignment name source:     {detection_source}")
        print(f"nhphyml_run_pk:            {run_metadata['nhphyml_run_pk']}")
        print(f"run_name:                  {run_metadata['run_name']}")
        print(f"Input alignment files:     {len(fasta_files)}")
        print(f"Distinct orthologs:        {actual_genes}")
        print(f"Distinct species:          {len(summary_rows)}")
        print(f"Alignments with GC4:       {alignments_with_gc4}")
        print(f"Detailed zero-GC4 bases:   {zero_details}")
        print(f"Species undefined GC4:     {zero_species}")
        print(f"Detailed output rows:      {len(detailed_rows)}")
        print(f"Species summary rows:      {len(summary_rows)}")
        print(f"Detailed TSV:              {args.output_detailed}")
        print(f"Species summary TSV:       {args.output_species_summary}")
    except (OSError, ValueError, RuntimeError, UnicodeError) as exc:
        parser.exit(1, f"[ERROR] {exc}\n")


if __name__ == "__main__":
    main()
