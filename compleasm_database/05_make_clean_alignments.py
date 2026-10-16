#!/usr/bin/env python3
"""
05_make_clean_alignments.py

Re-align cohort-specific shared orthologs for the Compleasm -> MACSE -> nhPhyML
pipeline.

For each requested cohort (80, 70, 60, and/or 50 genomes), this script treats:

    <cohort>/included_accessions.txt
    <cohort>/shared_orthologs.txt

as authoritative. It retrieves the original unaligned CDS sequences from
01_pre_alignments, creates a cohort-specific unaligned FASTA for every requested
ortholog, re-runs MACSE independently for that taxon set, validates the new
alignment, and only publishes the cohort when every requested accession x
ortholog combination passes final MACSE QC.

Important design choices
------------------------
* Cohorts are re-aligned independently. Existing 85-genome MACSE alignments are
  NOT subsetted.
* Accession matching is exact. There is no accession-root/version fallback.
* MACSE frameshifts are counted at the aligned-codon level: a complete codon
  containing one or more "!" characters is one frameshift event. A trailing
  partial aligned codon containing "!" is also one event.
* Internal stops are unambiguous frame-0 TAA/TAG/TGA codons before the terminal
  biological complete codon. A terminal stop is recorded but does not fail QC.
* Equal aligned lengths are required within each ortholog.
* Deliberate difference from 04b quality control: aligned length modulo 3 != 0
  is reported as a warning and does NOT by itself fail final_alignment_qc.
* Production publication is transactional: partial or failed work never
  replaces a previously completed cohort.
* Test mode uses a separate output namespace and cannot replace production
  outputs.

Production outputs per cohort
-----------------------------
    <cohort>/01_clean_ortholog_alignments/<odb12_id>.fasta
    <cohort>/final_macse_validation.tsv

Test outputs per cohort
-----------------------
    <cohort>/05_make_clean_alignments_test/01_clean_ortholog_alignments/
    <cohort>/05_make_clean_alignments_test/final_macse_validation.tsv

Examples
--------
Full production run:
    python 05_make_clean_alignments.py \
        /Users/rossoaa/projects/genomes \
        -j 8

Selected groups:
    python 05_make_clean_alignments.py \
        /Users/rossoaa/projects/genomes \
        --groups 80 70 \
        -j 8

Explicit source override:
    python 05_make_clean_alignments.py \
        /Users/rossoaa/projects/genomes \
        --source-prealignment-dir \
        /Users/rossoaa/projects/genomes/records/compleasm/alignments/01_pre_alignments \
        -j 8

Deterministic test mode:
    python 05_make_clean_alignments.py \
        /Users/rossoaa/projects/genomes \
        --groups 80 \
        --test \
        -j 4
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import traceback
import uuid
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple


# =============================================================================
# Fixed pipeline settings
# =============================================================================

GROUP_DIRS: Dict[int, str] = {
    80: "80_2026oct14",
    70: "70_2026oct14",
    60: "60_2026oct14",
    50: "50_2026oct14",
}

DEFAULT_GROUPS: Tuple[int, ...] = (80, 70, 60, 50)
TEST_ORTHOLOG_COUNT = 5

QC_VERSION = "final_macse_validation_v1"
SQL_NULL = r"\N"

CANONICAL_DNA: Set[str] = set("ACGT")
IUPAC_AMBIGUITY: Set[str] = set("RYSWKMBDHVN")
ACCEPTED_MACSE_ALPHABET: Set[str] = (
    CANONICAL_DNA | IUPAC_AMBIGUITY | {"-", "!"}
)
STOP_CODONS: Set[str] = {"TAA", "TAG", "TGA"}  # NCBI genetic code 1

OUTPUT_COLUMNS = [
    "accession",
    "odb12_id",
    "final_alignment_qc",
    "macse_qc",
    "macse_frameshift",
    "macse_internal_stop",
    "macse_terminal_stop",
    "macse_failure_reasons",
    "qc_version",
]


class PipelineError(RuntimeError):
    """Raised when an integrity condition makes cohort publication unsafe."""


# =============================================================================
# Generic helpers
# =============================================================================

def positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Expected an integer, got: {value!r}") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("Value must be >= 1.")
    return parsed


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def short_list(values: Iterable[str], limit: int = 10) -> str:
    vals = list(values)
    if not vals:
        return "(none)"
    shown = ", ".join(vals[:limit])
    if len(vals) > limit:
        shown += f", ... (+{len(vals) - limit} more)"
    return shown


def ensure_safe_token(value: str, label: str, path: Path, line_number: int) -> None:
    if not value:
        raise PipelineError(f"Blank {label} in {path} at line {line_number}.")
    if any(ch.isspace() for ch in value):
        raise PipelineError(
            f"{label} must be one whitespace-free token in {path} at line "
            f"{line_number}: {value!r}"
        )
    if "\x00" in value:
        raise PipelineError(
            f"NUL character is not allowed in {label} in {path} at line {line_number}."
        )


def ensure_safe_filename_stem(
    value: str, label: str, path: Path, line_number: int
) -> None:
    ensure_safe_token(value, label, path, line_number)
    if value in {".", ".."} or "/" in value or "\\" in value:
        raise PipelineError(
            f"Unsafe filename-forming {label} in {path} at line {line_number}: "
            f"{value!r}"
        )


def read_authoritative_list(
    path: Path,
    label: str,
    *,
    filename_stem: bool = False,
) -> List[str]:
    """
    Read a strict one-token-per-line authoritative list.

    Interior blank lines are errors; duplicates are errors. A final newline is
    normal and is not considered a blank record by splitlines().
    """
    if not path.is_file():
        raise FileNotFoundError(f"Required {label} file not found: {path}")

    text = path.read_text(encoding="utf-8-sig")
    lines = text.splitlines()
    if not lines:
        raise PipelineError(f"Required {label} file is empty: {path}")

    items: List[str] = []
    seen: Set[str] = set()

    for line_number, raw_line in enumerate(lines, start=1):
        value = raw_line.strip()
        if filename_stem:
            ensure_safe_filename_stem(value, label, path, line_number)
        else:
            ensure_safe_token(value, label, path, line_number)

        if value in seen:
            raise PipelineError(
                f"Duplicate {label} {value!r} in {path} at line {line_number}."
            )
        seen.add(value)
        items.append(value)

    return items


# =============================================================================
# Strict FASTA parsing and writing
# =============================================================================

def parse_fasta_strict(path: Path, *, uppercase: bool) -> Dict[str, str]:
    """
    Parse FASTA using the first whitespace-delimited header token as the ID.

    Enforces:
      * nonempty file content/records;
      * nonblank headers;
      * unique first-token identifiers;
      * no sequence before the first header;
      * nonempty sequences.

    Sequence-line whitespace is formatting and is removed. When uppercase=False,
    nucleotide letter case is otherwise preserved so raw CDS sequence text is
    not case-normalized before MACSE. MACSE output is parsed with uppercase=True
    to reproduce the upstream QC convention.
    """
    if not path.is_file():
        raise FileNotFoundError(f"FASTA not found: {path}")
    if path.stat().st_size == 0:
        raise PipelineError(f"FASTA is empty: {path}")

    records: Dict[str, str] = {}
    current_id: str | None = None
    sequence_parts: List[str] = []

    def finalize_record() -> None:
        nonlocal current_id, sequence_parts
        if current_id is None:
            return
        sequence = "".join(sequence_parts)
        if uppercase:
            sequence = sequence.upper()
        if not sequence:
            raise PipelineError(
                f"Empty FASTA sequence for identifier {current_id!r} in {path}"
            )
        if current_id in records:
            raise PipelineError(
                f"Duplicate FASTA identifier {current_id!r} in {path}"
            )
        records[current_id] = sequence

    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue

            if stripped.startswith(">"):
                finalize_record()
                title = stripped[1:].strip()
                if not title:
                    raise PipelineError(
                        f"Empty FASTA header in {path} at line {line_number}"
                    )
                identifier = title.split()[0]
                if not identifier:
                    raise PipelineError(
                        f"Blank first-token FASTA identifier in {path} "
                        f"at line {line_number}"
                    )
                # Catch a repeated identifier immediately, including the case in
                # which finalize_record() just inserted the prior occurrence.
                if identifier in records:
                    raise PipelineError(
                        f"Duplicate FASTA identifier {identifier!r} in {path}"
                    )
                current_id = identifier
                sequence_parts = []
            else:
                if current_id is None:
                    raise PipelineError(
                        f"Sequence content before first FASTA header in {path} "
                        f"at line {line_number}"
                    )
                sequence_parts.append(re.sub(r"\s+", "", stripped))

    finalize_record()

    if not records:
        raise PipelineError(f"FASTA contains no records: {path}")
    return records


def write_fasta(path: Path, records: Dict[str, str]) -> None:
    """Write exact first-token headers in deterministic alphabetical ID order."""
    if not records:
        raise PipelineError(f"Refusing to write empty FASTA: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for accession in sorted(records):
            sequence = records[accession]
            if not sequence:
                raise PipelineError(
                    f"Refusing to write empty sequence for {accession!r}: {path}"
                )
            handle.write(f">{accession}\n{sequence}\n")


def validate_expected_accessions(
    records: Dict[str, str],
    expected_accessions: Sequence[str],
    path: Path,
) -> None:
    expected = set(expected_accessions)
    observed = set(records)

    missing = sorted(expected - observed)
    extra = sorted(observed - expected)

    if missing or extra or len(records) != len(expected_accessions):
        raise PipelineError(
            f"FASTA accession set/count mismatch in {path}. "
            f"Expected={len(expected_accessions)}; observed={len(records)}; "
            f"missing={short_list(missing)}; extra={short_list(extra)}"
        )


def validate_macse_alphabet(sequence: str, path: Path, accession: str) -> None:
    invalid = sorted(set(sequence.upper()) - ACCEPTED_MACSE_ALPHABET)
    if invalid:
        printable = ", ".join(repr(x) for x in invalid)
        raise PipelineError(
            f"Unrecognized MACSE NT character(s) for accession {accession!r} "
            f"in {path}: {printable}. Accepted alphabet is ACGT, IUPAC "
            "RYSWKMBDHVN, '-' and '!'."
        )


def validate_alignment_structure(
    records: Dict[str, str],
    expected_accessions: Sequence[str],
    path: Path,
) -> int:
    """Validate exact occupancy, accepted alphabet, and one common aligned length."""
    validate_expected_accessions(records, expected_accessions, path)

    lengths: Dict[int, int] = {}
    for accession, sequence in records.items():
        if not sequence:
            raise PipelineError(
                f"Empty aligned sequence for accession {accession!r} in {path}"
            )
        validate_macse_alphabet(sequence, path, accession)
        lengths[len(sequence)] = lengths.get(len(sequence), 0) + 1

    if len(lengths) != 1:
        formatted = ", ".join(
            f"{length}:{count}" for length, count in sorted(lengths.items())
        )
        raise PipelineError(
            f"MACSE FASTA does not have one common aligned nucleotide length: "
            f"{path}. length:count = {formatted}"
        )

    return next(iter(lengths))


# =============================================================================
# MACSE biological QC
# =============================================================================

def evaluate_macse_sequence(sequence: str, path: Path, accession: str) -> Dict[str, object]:
    """
    Evaluate one MACSE nucleotide sequence under NCBI genetic code 1.

    Only complete frame-0 codons are evaluated for stops. If the aligned length
    leaves 1-2 trailing positions, those positions do not cause failure merely
    by existing. A trailing partial codon containing one or more "!" symbols is
    nevertheless counted as one frameshift event.
    """
    sequence = sequence.upper()
    validate_macse_alphabet(sequence, path, accession)

    complete_length = len(sequence) - (len(sequence) % 3)
    complete_codons = [
        sequence[i : i + 3] for i in range(0, complete_length, 3)
    ]
    trailing_partial = sequence[complete_length:]

    frameshift_count = sum(1 for codon in complete_codons if "!" in codon)
    if trailing_partial and "!" in trailing_partial:
        frameshift_count += 1

    non_gap_indices = [
        i for i, codon in enumerate(complete_codons) if codon != "---"
    ]

    terminal_stop = False
    internal_stop_count = 0

    if non_gap_indices:
        last_biological_index = non_gap_indices[-1]
        terminal_stop = complete_codons[last_biological_index] in STOP_CODONS

        for i, codon in enumerate(complete_codons):
            # Ambiguous/gapped/frameshift codons are never guessed to be stops.
            if codon in STOP_CODONS:
                if terminal_stop and i == last_biological_index:
                    continue
                internal_stop_count += 1

    reasons: List[str] = []
    if frameshift_count > 0:
        reasons.append("frameshift")
    if internal_stop_count > 0:
        reasons.append("internal_stop")

    macse_qc = "fail" if reasons else "pass"

    return {
        "macse_qc": macse_qc,
        "macse_frameshift": frameshift_count,
        "macse_internal_stop": internal_stop_count,
        "macse_terminal_stop": bool_text(terminal_stop),
        "macse_failure_reasons": ";".join(reasons) if reasons else SQL_NULL,
    }


def build_validation_row(
    accession: str,
    odb12_id: str,
    metrics: Dict[str, object],
) -> Dict[str, object]:
    macse_qc = str(metrics["macse_qc"])
    return {
        "accession": accession,
        "odb12_id": odb12_id,
        "final_alignment_qc": macse_qc,
        "macse_qc": macse_qc,
        "macse_frameshift": metrics["macse_frameshift"],
        "macse_internal_stop": metrics["macse_internal_stop"],
        "macse_terminal_stop": metrics["macse_terminal_stop"],
        "macse_failure_reasons": metrics["macse_failure_reasons"],
        "qc_version": QC_VERSION,
    }


# =============================================================================
# TSV helpers
# =============================================================================

def write_validation_tsv(rows: Sequence[Dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = path.with_name(path.name + ".tmp")

    with temp_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            if list(row.keys()) != OUTPUT_COLUMNS:
                missing = [c for c in OUTPUT_COLUMNS if c not in row]
                extra = [c for c in row if c not in OUTPUT_COLUMNS]
                raise PipelineError(
                    f"Internal TSV column mismatch. missing={missing}; extra={extra}"
                )
            writer.writerow(row)

    os.replace(temp_path, path)


def read_validation_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise PipelineError(f"Validation TSV missing or empty: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != OUTPUT_COLUMNS:
            raise PipelineError(
                f"Unexpected columns in {path}. Expected {OUTPUT_COLUMNS}; "
                f"found {reader.fieldnames}"
            )
        return list(reader)


# =============================================================================
# One-ortholog worker
# =============================================================================

def process_ortholog(
    odb12_id: str,
    included_accessions: Sequence[str],
    source_prealignment_dir: str,
    temp_prealignment_dir: str,
    temp_macse_root: str,
    staged_clean_dir: str,
    macse_executable: str,
) -> Dict[str, object]:
    """
    Build, align, validate, and stage one ortholog.

    The worker never appends to a shared TSV. It returns structured QC rows to
    the parent process.
    """
    try:
        source_dir = Path(source_prealignment_dir)
        pre_dir = Path(temp_prealignment_dir)
        macse_root = Path(temp_macse_root)
        clean_dir = Path(staged_clean_dir)

        source_path = source_dir / f"{odb12_id}_unaligned.fasta"
        if not source_path.is_file():
            raise FileNotFoundError(
                f"Source unaligned FASTA not found for odb12_id={odb12_id}: "
                f"{source_path}"
            )
        if source_path.stat().st_size == 0:
            raise PipelineError(
                f"Source unaligned FASTA is empty for odb12_id={odb12_id}: "
                f"{source_path}"
            )

        source_records = parse_fasta_strict(source_path, uppercase=False)

        missing = [
            accession
            for accession in included_accessions
            if accession not in source_records
        ]
        if missing:
            raise PipelineError(
                f"odb12_id={odb12_id}: exact accession(s) missing from source "
                f"unaligned FASTA {source_path}: {short_list(sorted(missing))}"
            )

        selected_records = {
            accession: source_records[accession]
            for accession in sorted(included_accessions)
        }

        cohort_pre_path = pre_dir / f"{odb12_id}_unaligned.fasta"
        write_fasta(cohort_pre_path, selected_records)

        # Re-parse and validate the cohort-specific raw input before MACSE.
        reparsed_pre = parse_fasta_strict(cohort_pre_path, uppercase=False)
        validate_expected_accessions(
            reparsed_pre, included_accessions, cohort_pre_path
        )

        # Ensure the script did not alter raw sequence characters while staging.
        for accession in included_accessions:
            if reparsed_pre[accession] != source_records[accession]:
                raise PipelineError(
                    f"odb12_id={odb12_id}: raw sequence changed while creating "
                    f"cohort pre-alignment for accession={accession}"
                )

        ortholog_macse_dir = macse_root / odb12_id
        ortholog_macse_dir.mkdir(parents=True, exist_ok=False)
        nt_path = ortholog_macse_dir / f"{odb12_id}_NT.fasta"
        aa_path = ortholog_macse_dir / f"{odb12_id}_AA.fasta"

        command = [
            macse_executable,
            "-prog",
            "alignSequences",
            "-seq",
            str(cohort_pre_path),
            "-out_NT",
            str(nt_path),
            "-out_AA",
            str(aa_path),
        ]

        completed = subprocess.run(
            command,
            cwd=ortholog_macse_dir,
            check=False,
            capture_output=True,
            text=True,
        )

        if completed.returncode != 0:
            raise PipelineError(
                f"MACSE exited nonzero for odb12_id={odb12_id}; "
                f"returncode={completed.returncode}; "
                f"stderr={completed.stderr.strip()!r}; "
                f"stdout={completed.stdout.strip()!r}"
            )

        if not nt_path.is_file() or nt_path.stat().st_size == 0:
            raise PipelineError(
                f"MACSE reported success but did not produce a nonempty NT FASTA "
                f"for odb12_id={odb12_id}: {nt_path}"
            )

        macse_records = parse_fasta_strict(nt_path, uppercase=True)
        aligned_length = validate_alignment_structure(
            macse_records, included_accessions, nt_path
        )

        rows: List[Dict[str, object]] = []
        failed_accessions: List[str] = []

        for accession in sorted(included_accessions):
            metrics = evaluate_macse_sequence(
                macse_records[accession], nt_path, accession
            )
            row = build_validation_row(accession, odb12_id, metrics)
            rows.append(row)
            if row["macse_qc"] != "pass":
                failed_accessions.append(accession)

        # Stage a deterministic final FASTA. It remains unpublished unless the
        # entire cohort passes every later integrity check.
        staged_path = clean_dir / f"{odb12_id}.fasta"
        deterministic_records = {
            accession: macse_records[accession]
            for accession in sorted(included_accessions)
        }
        write_fasta(staged_path, deterministic_records)

        # Validate what was actually staged, not merely the MACSE source file.
        staged_records = parse_fasta_strict(staged_path, uppercase=True)
        staged_length = validate_alignment_structure(
            staged_records, included_accessions, staged_path
        )
        if staged_length != aligned_length:
            raise PipelineError(
                f"odb12_id={odb12_id}: staged aligned length changed unexpectedly "
                f"({aligned_length} -> {staged_length})"
            )

        return {
            "status": "ok",
            "odb12_id": odb12_id,
            "rows": rows,
            "aligned_length": aligned_length,
            "alignment_length_mod_3": aligned_length % 3,
            "failed_accessions": failed_accessions,
            "error": None,
        }

    except Exception as exc:
        return {
            "status": "error",
            "odb12_id": odb12_id,
            "rows": [],
            "aligned_length": None,
            "alignment_length_mod_3": None,
            "failed_accessions": [],
            "error": (
                f"{type(exc).__name__}: {exc}\n"
                f"{traceback.format_exc()}"
            ),
        }


# =============================================================================
# Cohort-level validation and publication
# =============================================================================

def validate_row_integrity(
    rows: Sequence[Dict[str, object]],
    ortholog_ids: Sequence[str],
    included_accessions: Sequence[str],
) -> None:
    expected_rows = len(ortholog_ids) * len(included_accessions)
    if len(rows) != expected_rows:
        raise PipelineError(
            f"Validation row-count mismatch. Expected {expected_rows}; "
            f"observed {len(rows)}."
        )

    expected_pairs = {
        (accession, odb12_id)
        for odb12_id in ortholog_ids
        for accession in included_accessions
    }
    observed_pairs: Set[Tuple[str, str]] = set()

    for row in rows:
        pair = (str(row["accession"]), str(row["odb12_id"]))
        if pair in observed_pairs:
            raise PipelineError(
                f"Duplicate (accession, odb12_id) validation pair: {pair}"
            )
        observed_pairs.add(pair)

    missing_pairs = expected_pairs - observed_pairs
    extra_pairs = observed_pairs - expected_pairs
    if missing_pairs or extra_pairs:
        raise PipelineError(
            "Validation pair-set mismatch. "
            f"missing={short_list(sorted(f'{a}|{o}' for a, o in missing_pairs))}; "
            f"extra={short_list(sorted(f'{a}|{o}' for a, o in extra_pairs))}"
        )


def validate_complete_output(
    alignment_dir: Path,
    validation_tsv: Path,
    ortholog_ids: Sequence[str],
    included_accessions: Sequence[str],
    *,
    require_all_pass: bool,
) -> Dict[str, int]:
    """Re-validate the exact completed output set from disk."""
    if not alignment_dir.is_dir():
        raise PipelineError(f"Final alignment directory not found: {alignment_dir}")

    expected_orthologs = set(ortholog_ids)
    fasta_paths = sorted(alignment_dir.glob("*.fasta"))
    observed_orthologs = {path.stem for path in fasta_paths}

    if len(fasta_paths) != len(expected_orthologs):
        raise PipelineError(
            f"Final FASTA count mismatch in {alignment_dir}. "
            f"Expected={len(expected_orthologs)}; observed={len(fasta_paths)}."
        )

    missing = sorted(expected_orthologs - observed_orthologs)
    stale = sorted(observed_orthologs - expected_orthologs)
    if missing or stale:
        raise PipelineError(
            f"Final FASTA filename set mismatch in {alignment_dir}. "
            f"missing={short_list(missing)}; stale={short_list(stale)}"
        )

    warning_mod3 = 0
    recomputed_failures = 0
    recomputed_rows: Dict[Tuple[str, str], Dict[str, object]] = {}

    for fasta_path in fasta_paths:
        records = parse_fasta_strict(fasta_path, uppercase=True)
        aligned_length = validate_alignment_structure(
            records, included_accessions, fasta_path
        )
        if aligned_length % 3 != 0:
            warning_mod3 += 1

        odb12_id = fasta_path.stem
        for accession in included_accessions:
            metrics = evaluate_macse_sequence(
                records[accession], fasta_path, accession
            )
            expected_row = build_validation_row(accession, odb12_id, metrics)
            recomputed_rows[(accession, odb12_id)] = expected_row
            if metrics["macse_qc"] != "pass":
                recomputed_failures += 1

    tsv_rows = read_validation_tsv(validation_tsv)
    expected_rows = len(ortholog_ids) * len(included_accessions)
    if len(tsv_rows) != expected_rows:
        raise PipelineError(
            f"Final TSV row-count mismatch in {validation_tsv}. "
            f"Expected={expected_rows}; observed={len(tsv_rows)}."
        )

    seen_pairs: Set[Tuple[str, str]] = set()
    tsv_failures = 0
    for row in tsv_rows:
        pair = (row["accession"], row["odb12_id"])
        if pair in seen_pairs:
            raise PipelineError(
                f"Duplicate pair in final TSV {validation_tsv}: {pair}"
            )
        seen_pairs.add(pair)

        if row["qc_version"] != QC_VERSION:
            raise PipelineError(
                f"Unexpected qc_version for pair {pair} in {validation_tsv}: "
                f"{row['qc_version']!r}"
            )

        expected_row = recomputed_rows.get(pair)
        if expected_row is None:
            raise PipelineError(
                f"TSV pair has no corresponding validated FASTA record: {pair}"
            )
        mismatched_fields = [
            column
            for column in OUTPUT_COLUMNS
            if str(row[column]) != str(expected_row[column])
        ]
        if mismatched_fields:
            details = ", ".join(
                f"{column}:tsv={row[column]!r}/recomputed={str(expected_row[column])!r}"
                for column in mismatched_fields
            )
            raise PipelineError(
                f"Final TSV does not match recomputed MACSE QC for pair {pair}: "
                f"{details}"
            )

        if row["macse_qc"] != "pass" or row["final_alignment_qc"] != "pass":
            tsv_failures += 1

        if row["macse_qc"] == "pass":
            if row["macse_frameshift"] != "0":
                raise PipelineError(
                    f"macse_qc=pass but frameshift != 0 for pair {pair}"
                )
            if row["macse_internal_stop"] != "0":
                raise PipelineError(
                    f"macse_qc=pass but internal_stop != 0 for pair {pair}"
                )

    expected_pairs = {
        (accession, odb12_id)
        for odb12_id in ortholog_ids
        for accession in included_accessions
    }
    if seen_pairs != expected_pairs:
        missing_pairs = expected_pairs - seen_pairs
        extra_pairs = seen_pairs - expected_pairs
        raise PipelineError(
            f"Final TSV pair-set mismatch. "
            f"missing={short_list(sorted(f'{a}|{o}' for a, o in missing_pairs))}; "
            f"extra={short_list(sorted(f'{a}|{o}' for a, o in extra_pairs))}"
        )

    if require_all_pass and (recomputed_failures or tsv_failures):
        raise PipelineError(
            f"Completed output contains biological QC failures. "
            f"Recomputed failures={recomputed_failures}; "
            f"TSV failures={tsv_failures}."
        )

    return {
        "fasta_count": len(fasta_paths),
        "validation_rows": len(tsv_rows),
        "warning_mod3": warning_mod3,
        "recomputed_failures": recomputed_failures,
        "tsv_failures": tsv_failures,
    }


def report_existing_output_state(
    final_dir: Path,
    expected_ortholog_ids: Sequence[str],
) -> None:
    """Detect stale/missing production/test FASTA names before a replacement run."""
    if not final_dir.exists():
        return
    if not final_dir.is_dir():
        print(
            f"[WARNING] Existing final alignment path is not a directory: "
            f"{final_dir}"
        )
        return

    expected = set(expected_ortholog_ids)
    observed = {path.stem for path in final_dir.glob("*.fasta")}
    stale = sorted(observed - expected)
    missing = sorted(expected - observed)

    if stale or missing or len(observed) != len(expected):
        print(
            "[WARNING] Existing completed-output namespace does not exactly "
            "match the current shared_orthologs definition."
        )
        print(f"          Existing FASTAs: {len(observed)}")
        print(f"          Expected FASTAs: {len(expected)}")
        if stale:
            print(f"          Stale FASTAs:    {short_list(stale)}")
        if missing:
            print(f"          Missing FASTAs:  {short_list(missing)}")
        print(
            "          It will NOT be touched unless the new cohort run passes "
            "all validation."
        )


def publish_transactionally(
    staged_dir: Path,
    staged_tsv: Path,
    final_dir: Path,
    final_tsv: Path,
    ortholog_ids: Sequence[str],
    included_accessions: Sequence[str],
    run_id: str,
) -> Dict[str, int]:
    """
    Replace the two completed outputs only after staging is fully valid.

    Existing completed outputs are temporarily backed up. If publication or
    post-publication validation fails, the previous completed outputs are
    restored.
    """
    final_dir.parent.mkdir(parents=True, exist_ok=True)
    final_tsv.parent.mkdir(parents=True, exist_ok=True)

    backup_dir = final_dir.parent / f".{final_dir.name}.backup.{run_id}"
    backup_tsv = final_tsv.parent / f".{final_tsv.name}.backup.{run_id}"

    if backup_dir.exists() or backup_tsv.exists():
        raise PipelineError(
            f"Refusing publication because backup path already exists: "
            f"{backup_dir} or {backup_tsv}"
        )

    if final_dir.exists() and not final_dir.is_dir():
        raise PipelineError(
            f"Existing final alignment path is not a directory: {final_dir}"
        )
    if final_tsv.exists() and not final_tsv.is_file():
        raise PipelineError(
            f"Existing final validation TSV path is not a regular file: {final_tsv}"
        )

    had_dir = final_dir.exists()
    had_tsv = final_tsv.exists()

    try:
        if had_dir:
            os.replace(final_dir, backup_dir)
        if had_tsv:
            os.replace(final_tsv, backup_tsv)

        os.replace(staged_dir, final_dir)
        os.replace(staged_tsv, final_tsv)

        stats = validate_complete_output(
            final_dir,
            final_tsv,
            ortholog_ids,
            included_accessions,
            require_all_pass=True,
        )

    except BaseException:
        # Remove any newly published outputs, then restore the previous valid
        # namespace when it existed.
        if final_dir.exists():
            if final_dir.is_dir():
                shutil.rmtree(final_dir)
            else:
                final_dir.unlink()
        if final_tsv.exists():
            if final_tsv.is_dir():
                shutil.rmtree(final_tsv)
            else:
                final_tsv.unlink()

        if backup_dir.exists():
            os.replace(backup_dir, final_dir)
        if backup_tsv.exists():
            os.replace(backup_tsv, final_tsv)
        raise

    else:
        if backup_dir.exists():
            shutil.rmtree(backup_dir)
        if backup_tsv.exists():
            backup_tsv.unlink()
        return stats


def archive_failed_run(
    work_root: Path,
    group_dir: Path,
    run_id: str,
    reason_lines: Sequence[str],
    *,
    test_mode: bool,
) -> Path:
    """Preserve failed work under an explicitly non-production namespace."""
    if work_root.exists():
        failure_note = work_root / "FAILED_RUN.txt"
        with failure_note.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("05_make_clean_alignments.py failed cohort run\n")
            handle.write(f"run_id: {run_id}\n")
            handle.write(f"test_mode: {test_mode}\n\n")
            for line in reason_lines:
                handle.write(str(line).rstrip() + "\n")

    failed_parent = group_dir / (
        ".05_make_clean_alignments_test_failed"
        if test_mode
        else ".05_make_clean_alignments_failed"
    )
    failed_parent.mkdir(parents=True, exist_ok=True)
    destination = failed_parent / run_id

    suffix = 1
    while destination.exists():
        destination = failed_parent / f"{run_id}.{suffix}"
        suffix += 1

    if work_root.exists():
        os.replace(work_root, destination)
    else:
        destination.mkdir(parents=True, exist_ok=False)
        (destination / "FAILED_RUN.txt").write_text(
            "\n".join(reason_lines) + "\n", encoding="utf-8"
        )

    return destination


# =============================================================================
# Cohort processing
# =============================================================================

def run_group(
    group_size: int,
    alignments_root: Path,
    source_prealignment_dir: Path,
    threads: int,
    macse_executable: str,
    *,
    test_mode: bool,
) -> Dict[str, object]:
    cohort_dir = alignments_root / GROUP_DIRS[group_size]
    included_path = cohort_dir / "included_accessions.txt"
    shared_path = cohort_dir / "shared_orthologs.txt"

    if not cohort_dir.is_dir():
        raise NotADirectoryError(f"Cohort directory not found: {cohort_dir}")

    included_accessions = read_authoritative_list(
        included_path, "accession", filename_stem=False
    )
    if len(included_accessions) != group_size:
        raise PipelineError(
            f"{included_path} must contain exactly {group_size} unique accessions; "
            f"found {len(included_accessions)}."
        )

    all_ortholog_ids = read_authoritative_list(
        shared_path, "ortholog ID", filename_stem=True
    )
    if not all_ortholog_ids:
        raise PipelineError(f"No ortholog IDs found in {shared_path}")

    requested_orthologs = sorted(all_ortholog_ids)
    if test_mode:
        requested_orthologs = requested_orthologs[:TEST_ORTHOLOG_COUNT]
        if not requested_orthologs:
            raise PipelineError(
                f"Test mode selected zero orthologs from {shared_path}"
            )

    output_root = (
        cohort_dir / "05_make_clean_alignments_test"
        if test_mode
        else cohort_dir
    )
    final_dir = output_root / "01_clean_ortholog_alignments"
    final_tsv = output_root / "final_macse_validation.tsv"

    run_id = (
        f"{'test' if test_mode else 'prod'}."
        f"{os.getpid()}.{uuid.uuid4().hex[:10]}"
    )
    work_root = cohort_dir / f".05_make_clean_alignments_tmp.{run_id}"
    pre_dir = work_root / "pre_alignments"
    macse_root = work_root / "macse"
    stage_root = work_root / "publish"
    staged_clean_dir = stage_root / "01_clean_ortholog_alignments"
    staged_tsv = stage_root / "final_macse_validation.tsv"

    pre_dir.mkdir(parents=True, exist_ok=False)
    macse_root.mkdir(parents=True, exist_ok=False)
    staged_clean_dir.mkdir(parents=True, exist_ok=False)

    print()
    print("=" * 78)
    print(f"Group: {group_size}")
    print("=" * 78)
    print(f"Cohort directory:                  {cohort_dir}")
    print(f"included_accessions.txt:           {included_path}")
    print(f"Included accessions:               {len(included_accessions)}")
    print(f"shared_orthologs.txt:              {shared_path}")
    print(f"Orthologs in authoritative file:   {len(all_ortholog_ids)}")
    if test_mode:
        print(f"Requested orthologs (TEST subset): {len(requested_orthologs)}")
    else:
        print(f"Requested orthologs:               {len(requested_orthologs)}")
    print(f"Temporary pre-alignment directory: {pre_dir}")
    print(f"Final alignment directory:         {final_dir}")
    print(f"Final validation TSV:              {final_tsv}")

    report_existing_output_state(final_dir, requested_orthologs)

    # Fail before launching MACSE if any required unaligned source file is absent.
    missing_sources = [
        source_prealignment_dir / f"{odb12_id}_unaligned.fasta"
        for odb12_id in requested_orthologs
        if not (source_prealignment_dir / f"{odb12_id}_unaligned.fasta").is_file()
    ]
    if missing_sources:
        reason = [
            f"Missing {len(missing_sources)} required source unaligned FASTA(s):",
            *[str(path) for path in missing_sources],
        ]
        failed_dir = archive_failed_run(
            work_root,
            cohort_dir,
            run_id,
            reason,
            test_mode=test_mode,
        )
        raise PipelineError(
            f"Group {group_size} failed source preflight. "
            f"Failed run archived at: {failed_dir}"
        )

    results: List[Dict[str, object]] = []
    total = len(requested_orthologs)
    completed_count = 0
    progress_every = max(1, min(50, total // 10 if total >= 10 else 1))

    with ProcessPoolExecutor(max_workers=min(threads, total)) as executor:
        future_map = {
            executor.submit(
                process_ortholog,
                odb12_id,
                tuple(included_accessions),
                str(source_prealignment_dir),
                str(pre_dir),
                str(macse_root),
                str(staged_clean_dir),
                macse_executable,
            ): odb12_id
            for odb12_id in requested_orthologs
        }

        for future in as_completed(future_map):
            odb12_id = future_map[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {
                    "status": "error",
                    "odb12_id": odb12_id,
                    "rows": [],
                    "aligned_length": None,
                    "alignment_length_mod_3": None,
                    "failed_accessions": [],
                    "error": (
                        f"Worker process exception: {type(exc).__name__}: {exc}\n"
                        f"{traceback.format_exc()}"
                    ),
                }

            results.append(result)
            completed_count += 1

            if (
                completed_count == total
                or completed_count % progress_every == 0
            ):
                print(
                    f"  Progress: {completed_count}/{total} ortholog workers "
                    "completed"
                )

    results.sort(key=lambda item: str(item["odb12_id"]))
    worker_errors = [
        result for result in results if result["status"] != "ok"
    ]

    successful_rows: List[Dict[str, object]] = []
    for result in results:
        successful_rows.extend(result.get("rows", []))

    successful_rows.sort(
        key=lambda row: (str(row["odb12_id"]), str(row["accession"]))
    )

    if worker_errors:
        # Preserve all rows that were safely produced, but label them partial so
        # they cannot be mistaken for a complete cohort.
        if successful_rows:
            write_validation_tsv(
                successful_rows,
                stage_root / "final_macse_validation.partial.tsv",
            )

        error_lines = [
            f"Group {group_size} worker failure(s): {len(worker_errors)}"
        ]
        for result in worker_errors:
            error_lines.append(
                f"\n--- odb12_id={result['odb12_id']} ---\n{result['error']}"
            )

        failed_dir = archive_failed_run(
            work_root,
            cohort_dir,
            run_id,
            error_lines,
            test_mode=test_mode,
        )
        raise PipelineError(
            f"Group {group_size} failed because {len(worker_errors)} "
            f"ortholog worker(s) failed. Failed run archived at: {failed_dir}"
        )

    # Every ortholog worker succeeded structurally; now the parent alone builds
    # the cohort-level validation table.
    validate_row_integrity(
        successful_rows, requested_orthologs, included_accessions
    )
    write_validation_tsv(successful_rows, staged_tsv)

    failed_rows = [
        row
        for row in successful_rows
        if row["macse_qc"] != "pass"
        or row["final_alignment_qc"] != "pass"
    ]

    warning_mod3_results = [
        result
        for result in results
        if int(result["alignment_length_mod_3"]) != 0
    ]

    if warning_mod3_results:
        print(
            f"[WARNING] Group {group_size}: "
            f"{len(warning_mod3_results)} alignment(s) have a common aligned "
            "length not divisible by 3. This is reported but does NOT fail QC "
            "by itself."
        )
        for result in warning_mod3_results[:10]:
            print(
                f"          {result['odb12_id']}: "
                f"length={result['aligned_length']} "
                f"(mod 3 = {result['alignment_length_mod_3']})"
            )
        if len(warning_mod3_results) > 10:
            print(
                f"          ... (+{len(warning_mod3_results) - 10} more)"
            )

    # Validate the complete staged filename/header/length/TSV structure even if
    # biological failures are present. require_all_pass=False lets us preserve
    # a coherent failed run for investigation.
    staged_stats = validate_complete_output(
        staged_clean_dir,
        staged_tsv,
        requested_orthologs,
        included_accessions,
        require_all_pass=False,
    )

    if failed_rows:
        failure_pairs = [
            (
                str(row["accession"]),
                str(row["odb12_id"]),
                str(row["macse_failure_reasons"]),
            )
            for row in failed_rows
        ]
        reason_lines = [
            f"Group {group_size} biological integrity failure.",
            f"Sequences failing final MACSE QC: {len(failed_rows)}",
            "No selected accession or ortholog was silently removed.",
            "The production/test completed output was NOT replaced.",
            "",
            "Failing accession\todb12_id\tmacse_failure_reasons",
            *[
                f"{accession}\t{odb12_id}\t{reasons}"
                for accession, odb12_id, reasons in failure_pairs
            ],
        ]
        failed_dir = archive_failed_run(
            work_root,
            cohort_dir,
            run_id,
            reason_lines,
            test_mode=test_mode,
        )
        raise PipelineError(
            f"Group {group_size} has {len(failed_rows)} sequence(s) failing "
            f"the newly re-run MACSE QC. Failed run archived at: {failed_dir}"
        )

    # Now require that the staged cohort is fully publishable.
    validate_complete_output(
        staged_clean_dir,
        staged_tsv,
        requested_orthologs,
        included_accessions,
        require_all_pass=True,
    )

    published_stats = publish_transactionally(
        staged_clean_dir,
        staged_tsv,
        final_dir,
        final_tsv,
        requested_orthologs,
        included_accessions,
        run_id,
    )

    # Publication moved the staged directory/TSV out of work_root; remove only
    # the remaining temporary pre-alignments/MACSE files.
    if work_root.exists():
        shutil.rmtree(work_root)

    expected_rows = len(requested_orthologs) * len(included_accessions)
    passing_rows = expected_rows  # all-pass is required before publication

    print()
    print(f"{group_size} genomes")
    print(f"  expected accessions per ortholog: {len(included_accessions)}")
    print(f"  requested orthologs:              {len(requested_orthologs)}")
    print(f"  MACSE alignments completed:       {len(results)}")
    print(f"  final FASTAs written:             {published_stats['fasta_count']}")
    print(f"  expected validation rows:         {expected_rows}")
    print(f"  actual validation rows:           {published_stats['validation_rows']}")
    print(f"  sequences passing final MACSE QC: {passing_rows}")
    print(f"  sequences failing final MACSE QC: 0")
    print(
        "  alignments with non-divisible-by-3 length: "
        f"{published_stats['warning_mod3']}"
    )
    print(f"  status: {'PASS (TEST)' if test_mode else 'PASS'}")

    return {
        "group": group_size,
        "status": "pass",
        "requested_orthologs": len(requested_orthologs),
        "validation_rows": published_stats["validation_rows"],
        "warning_mod3": published_stats["warning_mod3"],
        "output_root": str(output_root),
        "staged_warning_mod3": staged_stats["warning_mod3"],
    }


# =============================================================================
# CLI
# =============================================================================

def resolve_macse_executable() -> str:
    executable = shutil.which("macse")
    if not executable:
        raise FileNotFoundError(
            "MACSE executable 'macse' was not found on PATH. Activate the "
            "environment containing MACSE before running this script."
        )
    return str(Path(executable).resolve())


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Re-run MACSE independently for the authoritative 80/70/60/50 "
            "genome cohorts and publish only fully validated nhPhyML-ready "
            "nucleotide alignments."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "genomes_dir",
        help="Project genomes root, e.g. /Users/rossoaa/projects/genomes",
    )
    parser.add_argument(
        "--groups",
        nargs="+",
        type=int,
        choices=sorted(GROUP_DIRS),
        default=list(DEFAULT_GROUPS),
        metavar="N",
        help="One or more cohorts to process. Allowed: 80 70 60 50. Default: all.",
    )
    parser.add_argument(
        "--source-prealignment-dir",
        help=(
            "Override the original unaligned ortholog FASTA directory. Default: "
            "<genomes_dir>/records/compleasm/alignments/01_pre_alignments"
        ),
    )
    parser.add_argument(
        "-j",
        "--threads",
        type=positive_int,
        default=1,
        help="Maximum number of ortholog worker processes. Default: 1.",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help=(
            f"Process the first {TEST_ORTHOLOG_COUNT} ortholog IDs alphabetically "
            "per requested group and write only to the test-specific namespace."
        ),
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)

    if len(set(args.groups)) != len(args.groups):
        raise PipelineError(
            f"Duplicate group value(s) are not allowed: {args.groups}"
        )

    genomes_dir = Path(args.genomes_dir).expanduser().resolve()
    if not genomes_dir.is_dir():
        raise NotADirectoryError(f"Genomes root not found: {genomes_dir}")

    alignments_root = (
        genomes_dir / "records" / "compleasm" / "alignments"
    )
    if not alignments_root.is_dir():
        raise NotADirectoryError(f"Alignments root not found: {alignments_root}")

    source_prealignment_dir = (
        Path(args.source_prealignment_dir).expanduser().resolve()
        if args.source_prealignment_dir
        else alignments_root / "01_pre_alignments"
    )
    if not source_prealignment_dir.is_dir():
        raise NotADirectoryError(
            f"Source unaligned ortholog directory not found: "
            f"{source_prealignment_dir}"
        )

    macse_executable = resolve_macse_executable()

    requested_groups = list(args.groups)

    print("=" * 78)
    print("05_make_clean_alignments.py")
    print("=" * 78)
    print(f"Genomes root:                       {genomes_dir}")
    print(f"Alignments root:                    {alignments_root}")
    print(
        f"Source unaligned ortholog directory: {source_prealignment_dir}"
    )
    print(
        "Requested groups:                   "
        + " ".join(str(group) for group in requested_groups)
    )
    print(f"Threads:                            {args.threads}")
    print(f"MACSE executable:                   {macse_executable}")
    print(f"QC version:                         {QC_VERSION}")
    print(f"Test mode:                          {args.test}")
    if args.test:
        print(
            f"TEST MODE ACTIVE: first {TEST_ORTHOLOG_COUNT} ortholog IDs "
            "alphabetically per group; production outputs cannot be replaced."
        )

    failures: List[Tuple[int, str]] = []
    successes: List[Dict[str, object]] = []

    for group_size in requested_groups:
        try:
            summary = run_group(
                group_size=group_size,
                alignments_root=alignments_root,
                source_prealignment_dir=source_prealignment_dir,
                threads=args.threads,
                macse_executable=macse_executable,
                test_mode=args.test,
            )
            successes.append(summary)
        except Exception as exc:
            message = f"{type(exc).__name__}: {exc}"
            failures.append((group_size, message))
            print()
            print(f"[ERROR] Group {group_size} FAILED")
            print(f"        {message}", file=sys.stderr)

    print()
    print("=" * 78)
    print("RUN SUMMARY")
    print("=" * 78)
    for summary in successes:
        print(
            f"Group {summary['group']}: PASS | "
            f"orthologs={summary['requested_orthologs']} | "
            f"rows={summary['validation_rows']} | "
            f"nondiv3_alignments={summary['warning_mod3']} | "
            f"output={summary['output_root']}"
        )

    for group_size, message in failures:
        print(f"Group {group_size}: FAIL | {message}")

    if failures:
        print(
            f"\nOverall status: FAIL ({len(failures)} requested cohort(s) failed)."
        )
        return 1

    print("\nOverall status: PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\n[ERROR] Interrupted by user.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"[ERROR] {type(exc).__name__}: {exc}", file=sys.stderr)
        raise SystemExit(1)
