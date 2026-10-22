#!/usr/bin/env python3
"""
gc4_alignment_cumulative_v2.py

Calculate per-sequence GC4 content from one codon-aware FASTA alignment or
from every FASTA alignment in a directory using cumulative masks across all
sequences in each alignment.

For each alignment:
    1. Parse FASTA.
    2. Make a cumulative valid-site mask:
       keep a nucleotide column only if every sequence has a valid nucleotide.
    3. Make a complete-codon mask:
       keep a codon only if all three nucleotide columns are valid across all
       sequences.
    4. Make a cumulative GC4 mask:
       keep a codon only if every sequence has a GC4-family codon at that
       codon position. Only the third nucleotide position is used for GC4.
    5. Calculate GC4 for each sequence using the final retained GC4 third sites.
    6. Write:
       a) detailed per-alignment/per-sequence TSV
       b) species-level summary TSV aggregated across all alignments

Important:
    This script preserves alignment coordinates. It does not remove gaps from
    each sequence independently before codon parsing.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np


FASTA_SUFFIXES = {".fa", ".faa", ".fasta", ".fas", ".fna"}

GC4_PREFIX_TO_FAMILY = {
    "GT": "Val",
    "CC": "Pro",
    "AC": "Thr",
    "GC": "Ala",
    "GG": "Gly",
    "CT": "Leu",
    "CG": "Arg",
    "TC": "Ser",
}


# -----------------------------------------------------------------------------
# FASTA parsing and alignment validation
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


def calculate_gc4_for_alignment(
    fasta_path: Path,
    invalid_chars: str = "-!Nn",
    require_same_gc4_family: bool = False,
) -> List[Dict[str, object]]:
    """Calculate per-sequence GC4 using cumulative site/codon/GC4 masks."""
    records = parse_fasta(fasta_path)
    alignment_name = fasta_path.name
    alignment_length = validate_alignment(records, alignment_name=alignment_name)
    n_sequences = len(records)

    site_mask = make_complete_site_mask(
        records=records,
        invalid_chars=invalid_chars,
        alignment_name=alignment_name,
    )
    codon_mask = make_complete_codon_mask(site_mask)
    gc4_site_mask, codon_reason_counts = make_cumulative_gc4_mask(
        records=records,
        complete_codon_mask=codon_mask,
        require_same_gc4_family=require_same_gc4_family,
        alignment_name=alignment_name,
    )

    total_sites = alignment_length
    complete_sites = int(site_mask.sum())
    dropped_sites_incomplete = total_sites - complete_sites

    total_codons = alignment_length // 3
    complete_codons = int(codon_mask.sum() // 3)
    dropped_incomplete_codons = total_codons - complete_codons

    retained_gc4_sites = int(gc4_site_mask.sum())
    retained_gc4_codons = retained_gc4_sites

    invalid_char_counts = summarize_invalid_sites_by_character(records, invalid_chars)

    rows: List[Dict[str, object]] = []

    for sequence_id, sequence in records.items():
        gc4_sequence = extract_masked_sequence(sequence, gc4_site_mask)
        (
            gc4_content,
            g_count,
            c_count,
            a_count,
            t_count,
            n_count,
            total_valid_bases,
        ) = calculate_gc(gc4_sequence)

        rows.append(
            {
                "alignment": alignment_name,
                "sequence_id": sequence_id,
                "gc4": gc4_content,
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
                **invalid_char_counts,
            }
        )

    return rows


def summarize_by_species(detail_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    """
    Aggregate per-alignment rows into one species-level GC4 summary.

    GC4 is recalculated from summed counts:
        sum(G + C) / sum(A + T + G + C)
    not averaged across alignments.
    """
    grouped: Dict[str, Dict[str, object]] = {}

    numeric_sum_fields = [
        "g_gc4",
        "c_gc4",
        "a_gc4",
        "t_gc4",
        "n_gc4",
        "valid_gc4_bases",
        "retained_gc4_sites",
        "retained_gc4_codons",
        "alignment_length_nt",
        "total_sites_nt",
        "complete_sites_nt",
        "dropped_sites_incomplete_nt",
        "total_codons",
        "complete_codons",
        "dropped_incomplete_codons",
        "candidate_gc4_codons_all_sequences",
        "dropped_not_gc4_all_sequences",
        "dropped_mixed_gc4_family",
        "invalid_char_gap",
        "invalid_char_frameshift",
        "invalid_char_N",
    ]

    for row in detail_rows:
        sequence_id = str(row["sequence_id"])
        if sequence_id not in grouped:
            grouped[sequence_id] = {
                "sequence_id": sequence_id,
                "n_alignments": 0,
                "alignments_with_retained_gc4": 0,
            }
            for field in numeric_sum_fields:
                grouped[sequence_id][field] = 0

        out = grouped[sequence_id]
        out["n_alignments"] = int(out["n_alignments"]) + 1

        if int(row.get("valid_gc4_bases", 0)) > 0:
            out["alignments_with_retained_gc4"] = int(out["alignments_with_retained_gc4"]) + 1

        for field in numeric_sum_fields:
            out[field] = int(out.get(field, 0)) + int(row.get(field, 0) or 0)

    summary_rows: List[Dict[str, object]] = []
    for sequence_id in sorted(grouped):
        out = grouped[sequence_id]
        g_count = int(out["g_gc4"])
        c_count = int(out["c_gc4"])
        valid = int(out["valid_gc4_bases"])
        out["gc4"] = (g_count + c_count) / valid if valid else 0.0
        summary_rows.append(out)

    # Put GC4 near the front.
    preferred_order = [
        "sequence_id",
        "gc4",
        "g_gc4",
        "c_gc4",
        "a_gc4",
        "t_gc4",
        "n_gc4",
        "valid_gc4_bases",
        "retained_gc4_sites",
        "retained_gc4_codons",
        "n_alignments",
        "alignments_with_retained_gc4",
    ]

    reordered_rows: List[Dict[str, object]] = []
    for row in summary_rows:
        reordered = {key: row[key] for key in preferred_order if key in row}
        for key, value in row.items():
            if key not in reordered:
                reordered[key] = value
        reordered_rows.append(reordered)

    return reordered_rows


# -----------------------------------------------------------------------------
# Input/output
# -----------------------------------------------------------------------------


def iter_fasta_files(input_dir: Path) -> Iterable[Path]:
    """Yield FASTA files from a directory."""
    for path in sorted(input_dir.iterdir()):
        if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES:
            yield path


def collect_input_fastas(input_path: Path | None, input_dir: Path | None) -> List[Path]:
    """Return the list of FASTA files to process."""
    if input_path is not None:
        return [input_path]

    if input_dir is None:
        raise ValueError("Either --input or --input-dir is required")

    fasta_files = list(iter_fasta_files(input_dir))
    if not fasta_files:
        raise ValueError(f"No FASTA files found in {input_dir}")

    return fasta_files


def write_tsv(rows: List[Dict[str, object]], output_path: Path) -> None:
    """Write rows to TSV."""
    if not rows:
        raise ValueError(f"No rows to write for {output_path}")

    fieldnames: List[str] = []
    seen = set()
    for row in rows:
        for key in row.keys():
            if key not in seen:
                fieldnames.append(key)
                seen.add(key)

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate per-sequence GC4 from codon-aware alignments using "
            "cumulative masks across all sequences."
        )
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--input",
        type=Path,
        help="One codon-aware aligned FASTA file.",
    )
    input_group.add_argument(
        "--input-dir",
        type=Path,
        help="Directory containing codon-aware aligned FASTA files.",
    )

    parser.add_argument(
        "--output-detailed",
        type=Path,
        required=True,
        help="Detailed output TSV: one row per sequence per alignment.",
    )
    parser.add_argument(
        "--output-species-summary",
        type=Path,
        required=True,
        help="Species summary TSV: one row per sequence/species across all alignments.",
    )
    parser.add_argument(
        "--invalid-chars",
        default="-!Nn",
        help="Characters that make an alignment site invalid. Default: '-!Nn'.",
    )
    parser.add_argument(
        "--require-same-gc4-family",
        action="store_true",
        help=(
            "Require all sequences to have the same GC4 codon family at a "
            "retained codon. Default only requires each sequence to have some "
            "GC4-family codon."
        ),
    )

    args = parser.parse_args()

    fasta_files = collect_input_fastas(args.input, args.input_dir)

    detailed_rows: List[Dict[str, object]] = []
    failed_alignments: List[Tuple[str, str]] = []

    for fasta_path in fasta_files:
        try:
            detailed_rows.extend(
                calculate_gc4_for_alignment(
                    fasta_path=fasta_path,
                    invalid_chars=args.invalid_chars,
                    require_same_gc4_family=args.require_same_gc4_family,
                )
            )
        except Exception as exc:
            failed_alignments.append((str(fasta_path), str(exc)))

    if failed_alignments:
        failed_text = "\n".join(f"  {path}: {error}" for path, error in failed_alignments)
        raise RuntimeError(
            "One or more alignments failed. No output was written.\n" + failed_text
        )

    species_summary_rows = summarize_by_species(detailed_rows)

    write_tsv(detailed_rows, args.output_detailed)
    write_tsv(species_summary_rows, args.output_species_summary)

    print(f"Processed {len(fasta_files)} alignment file(s)")
    print(f"Wrote {len(detailed_rows)} detailed rows to {args.output_detailed}")
    print(
        f"Wrote {len(species_summary_rows)} species summary rows to "
        f"{args.output_species_summary}"
    )


if __name__ == "__main__":
    main()
