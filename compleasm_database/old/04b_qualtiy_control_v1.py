#!/usr/bin/env python3
"""
04b_qualtiy_control_v1.py

Ortholog quality control for the Compleasm -> MACSE -> nhPhyML pipeline.

Pipeline position
-----------------
Run after 04_prealign_to_macse_v3.py and after MACSE has produced the
02_macse_alignments/*_NT.fasta files, but before 05_cleaner.py modifies those
alignments.

This script is intentionally observational. It does not mask, edit, clean,
delete, or otherwise modify upstream FASTA files. It independently evaluates:

1. Raw CDS QC from 01_pre_alignments/<odb12_id>_unaligned.fasta
2. MACSE QC from 02_macse_alignments/<odb12_id>_NT.fasta

and then combines those two independent states into final_qc.

Raw sequence conventions
------------------------
* FASTA sequence whitespace is removed before any calculations and is NOT
  included in raw_length.
* Raw CDS input is expected to be ungapped. Gap-like/special characters such
  as '-', '.', and '!' are not accepted raw DNA symbols; if present, they are
  retained in raw_length and counted as invalid characters.
* Canonical DNA symbols: A C G T
* Accepted IUPAC ambiguity symbols: R Y S W K M B D H V N
  (N is therefore counted as ambiguous, as are all other accepted non-ACGT
  IUPAC ambiguity symbols.)
* raw_ambiguous_fraction = accepted ambiguity symbols / raw_length
* raw_invalid_fraction = characters outside the accepted DNA/IUPAC alphabet /
  raw_length
* Both fractions use raw_length as the denominator. For a zero-length record,
  parsing/validation fails before QC is attempted.

Stop-codon conventions
----------------------
NCBI genetic code 1 (the Standard Code) is used, matching MACSE's default
unless MACSE was explicitly run with a different -gc_def value. Complete
frame-0 codons are translated. Unambiguous TAA/TAG/TGA codons are stops under
this code. Ambiguous or invalid codons are represented as X for stop-QC
purposes rather than being guessed as stops.

For raw CDS, a terminal stop is true only when the raw sequence length is
exactly divisible by three AND its final complete in-frame codon is a stop.
If trailing 1-2 nt remain, a stop in the preceding complete codon is treated
as internal because the sequence continues beyond it.

MACSE conventions
-----------------
MACSE nucleotide alignments use '!' to mark frameshift-containing partial
codons. Ordinary '-' alignment gaps are NOT treated as frameshifts. Because a
single frameshifting codon may contain one or two '!' characters, this script
counts MACSE frameshift EVENTS at the codon level: each aligned codon
containing one or more '!' characters contributes one to macse_frameshift.

MACSE stop codons are represented by their nucleotide codons rather than by a
special NT character. Internal/terminal stop status is therefore determined
from aligned codons. Trailing all-gap codons ('---') are ignored when deciding
which biological codon is terminal.

Version-1 QC rules
------------------
Raw CDS fails only when one or more of these conditions is met:
* one or more internal stop codons;
* raw_ambiguous_fraction > MAX_RAW_AMBIGUOUS_FRACTION;
* raw_invalid_fraction > MAX_RAW_INVALID_FRACTION.

raw_length_mod_3 != 0 is recorded but DOES NOT cause failure in v1.
Presence/absence of a terminal stop DOES NOT cause failure in v1.

MACSE fails only when:
* macse_frameshift > 0; or
* macse_internal_stop > 0.

Presence/absence of a MACSE terminal stop DOES NOT cause failure in v1.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple

from Bio.Data import CodonTable


# ---------------------------------------------------------------------------
# Versioned/configurable QC settings
# ---------------------------------------------------------------------------

QC_VERSION = "v1"
GENETIC_CODE_ID = 1  # NCBI Standard Code; also MACSE's default genetic code.

# v1 raw-QC thresholds. A value equal to the threshold passes; only values
# strictly greater than the threshold fail.
MAX_RAW_AMBIGUOUS_FRACTION = 0.05
MAX_RAW_INVALID_FRACTION = 0.0

SQL_NULL = r"\N"

CANONICAL_DNA: Set[str] = set("ACGT")
IUPAC_AMBIGUITY: Set[str] = set("RYSWKMBDHVN")
ACCEPTED_RAW_ALPHABET: Set[str] = CANONICAL_DNA | IUPAC_AMBIGUITY
ACCEPTED_MACSE_ALPHABET: Set[str] = ACCEPTED_RAW_ALPHABET | {"-", "!"}

CODON_TABLE = CodonTable.unambiguous_dna_by_id[GENETIC_CODE_ID]
STOP_CODONS: Set[str] = set(CODON_TABLE.stop_codons)
FORWARD_TABLE: Dict[str, str] = dict(CODON_TABLE.forward_table)

OUTPUT_COLUMNS = [
    "accession",
    "odb12_id",
    "final_qc",
    "macse_qc",
    "raw_cds_qc",
    "raw_length",
    "raw_length_mod_3",
    "raw_internal_stop_count",
    "raw_terminal_stop",
    "raw_ambiguous_fraction",
    "raw_invalid_fraction",
    "raw_failure_reasons",
    "macse_frameshift",
    "macse_internal_stop",
    "macse_terminal_stop",
    "macse_failure_reasons",
    "qc_version",
]

FINAL_QC_MAP = {
    ("pass", "pass"): "passing",
    ("fail", "fail"): "unsuitable",
    ("fail", "pass"): "rescued",
    ("pass", "fail"): "alignment-problematic",
}

FINAL_QC_ORDER = [
    "passing",
    "unsuitable",
    "rescued",
    "alignment-problematic",
]


class QCValidationError(RuntimeError):
    """Raised when structural input problems make QC assignment unsafe."""


def sanitize_label(text: str) -> str:
    """Match 04_prealign_to_macse_v3.py cohort-label sanitization exactly."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_")


def root_acc(accession: str) -> str:
    """Match 04_prealign_to_macse_v3.py accession-root behavior."""
    accession = str(accession).strip()
    return accession.split(".")[0] if accession else accession


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def format_fraction(value: float) -> str:
    """Stable compact numeric representation for TSV output."""
    return f"{value:.12g}"


def read_project_accessions(path: Path) -> List[str]:
    """
    Read the same CSV/TXT accession input accepted by 04_prealign_to_macse_v3.py.

    CSV requires an 'accession' column. Non-CSV input is treated as a text file
    with one accession per nonblank/noncomment line; the first whitespace-
    separated token is used, matching the upstream script.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Project manifest/accessions file not found: {path}")

    if path.suffix.lower() == ".csv":
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None or "accession" not in reader.fieldnames:
                raise QCValidationError(
                    f"Manifest is missing required 'accession' column: {path}"
                )
            values = []
            for row in reader:
                value = str(row.get("accession") or "").strip()
                if value:
                    values.append(value)
    else:
        values = []
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                value = line.strip()
                if not value or value.startswith("#"):
                    continue
                values.append(value.split()[0])

    accessions = sorted(set(values))
    if not accessions:
        raise QCValidationError(f"No accessions were loaded from: {path}")
    return accessions


def resolve_shared_genes_path(
    genomes_dir: Path,
    cohort_label: str,
    explicit_path: Path | None,
) -> Path:
    """Use the same project-specific/fallback naming logic as the upstream step."""
    if explicit_path is not None:
        return explicit_path.resolve()

    records_dir = genomes_dir / "records/compleasm/records"
    project_specific = records_dir / f"shared_single_genes_{cohort_label}.csv"
    if project_specific.exists():
        return project_specific
    return records_dir / "shared_single_genes.csv"


def load_ortholog_ids(path: Path) -> List[str]:
    if not path.is_file():
        raise FileNotFoundError(f"Shared genes file not found: {path}")

    text = path.read_text(encoding="utf-8").strip()
    values = [item.strip() for item in text.split(",") if item.strip()]
    orthologs = sorted(set(values))
    if not orthologs:
        raise QCValidationError(f"No ODB12 ortholog IDs were loaded from: {path}")
    return orthologs


def parse_fasta_strict(path: Path) -> Dict[str, str]:
    """
    Parse FASTA while enforcing unique first-token identifiers.

    Sequence whitespace is removed. Empty identifiers, duplicate identifiers,
    sequence content before the first header, and empty sequences are structural
    errors because they can corrupt accession/ortholog matching.
    """
    records: Dict[str, str] = {}
    current_id: str | None = None
    sequence_parts: List[str] = []

    def finalize_record() -> None:
        nonlocal current_id, sequence_parts
        if current_id is None:
            return
        sequence = "".join(sequence_parts).upper()
        if not sequence:
            raise QCValidationError(
                f"Empty FASTA sequence for identifier '{current_id}' in {path}"
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
                    raise QCValidationError(
                        f"Empty FASTA header in {path} at line {line_number}"
                    )
                identifier = title.split()[0]
                if identifier in records:
                    raise QCValidationError(
                        f"Duplicate FASTA identifier '{identifier}' in {path}"
                    )
                current_id = identifier
                sequence_parts = []
            else:
                if current_id is None:
                    raise QCValidationError(
                        f"Sequence content before first FASTA header in {path} "
                        f"at line {line_number}"
                    )
                # FASTA whitespace is non-biological formatting and is excluded
                # from all sequence-length/fraction calculations.
                sequence_parts.append(re.sub(r"\s+", "", stripped))

    finalize_record()
    return records


def validate_raw_alphabet_and_measure(sequence: str) -> Tuple[int, int, float, float]:
    """Return ambiguous_count, invalid_count, ambiguous_fraction, invalid_fraction."""
    length = len(sequence)
    if length == 0:
        raise QCValidationError("Internal error: zero-length sequence reached raw QC")

    ambiguous_count = sum(base in IUPAC_AMBIGUITY for base in sequence)
    invalid_count = sum(base not in ACCEPTED_RAW_ALPHABET for base in sequence)
    return (
        ambiguous_count,
        invalid_count,
        ambiguous_count / length,
        invalid_count / length,
    )


def translate_complete_codons(sequence: str) -> List[str]:
    """
    Translate complete frame-0 codons under NCBI genetic code 1.

    Canonical codons use the configured codon table. Ambiguous or invalid
    codons become X for stop-counting purposes; they are never guessed to be a
    stop. Incomplete trailing 1-2 nt are not translated.
    """
    amino_acids: List[str] = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i : i + 3]
        if codon in STOP_CODONS:
            amino_acids.append("*")
        elif all(base in CANONICAL_DNA for base in codon):
            amino_acids.append(FORWARD_TABLE.get(codon, "X"))
        else:
            amino_acids.append("X")
    return amino_acids


def raw_stop_metrics(sequence: str) -> Tuple[int, bool]:
    """Count internal stops and detect a legitimate terminal raw-CDS stop."""
    translated = translate_complete_codons(sequence)
    if not translated:
        return 0, False

    terminal_stop = (len(sequence) % 3 == 0) and (translated[-1] == "*")
    if terminal_stop:
        internal_stop_count = translated[:-1].count("*")
    else:
        internal_stop_count = translated.count("*")
    return internal_stop_count, terminal_stop


def evaluate_raw_sequence(sequence: str) -> Dict[str, object]:
    raw_length = len(sequence)
    raw_length_mod_3 = raw_length % 3

    _, _, ambiguous_fraction, invalid_fraction = validate_raw_alphabet_and_measure(
        sequence
    )
    internal_stop_count, terminal_stop = raw_stop_metrics(sequence)

    reasons: List[str] = []
    if internal_stop_count > 0:
        reasons.append("internal_stop")
    if ambiguous_fraction > MAX_RAW_AMBIGUOUS_FRACTION:
        reasons.append("excessive_ambiguity")
    if invalid_fraction > MAX_RAW_INVALID_FRACTION:
        reasons.append("invalid_nucleotide")

    raw_qc = "fail" if reasons else "pass"

    return {
        "raw_cds_qc": raw_qc,
        "raw_length": raw_length,
        "raw_length_mod_3": raw_length_mod_3,
        "raw_internal_stop_count": internal_stop_count,
        "raw_terminal_stop": bool_text(terminal_stop),
        "raw_ambiguous_fraction": format_fraction(ambiguous_fraction),
        "raw_invalid_fraction": format_fraction(invalid_fraction),
        "raw_failure_reasons": ";".join(reasons) if reasons else SQL_NULL,
    }


def split_aligned_codons(sequence: str, path: Path, accession: str) -> List[str]:
    if len(sequence) % 3 != 0:
        raise QCValidationError(
            f"MACSE NT sequence length is not divisible by 3 for accession "
            f"'{accession}' in {path}: length={len(sequence)}"
        )
    return [sequence[i : i + 3] for i in range(0, len(sequence), 3)]


def validate_macse_alphabet(sequence: str, path: Path, accession: str) -> None:
    invalid = sorted(set(sequence) - ACCEPTED_MACSE_ALPHABET)
    if invalid:
        printable = ", ".join(repr(x) for x in invalid)
        raise QCValidationError(
            f"Unrecognized MACSE NT character(s) for accession '{accession}' "
            f"in {path}: {printable}. Accepted alphabet is IUPAC DNA plus '-' and '!'."
        )


def evaluate_macse_sequence(
    sequence: str,
    path: Path,
    accession: str,
) -> Dict[str, object]:
    validate_macse_alphabet(sequence, path, accession)
    codons = split_aligned_codons(sequence, path, accession)

    # MACSE frameshift event count is codon-level, not raw '!' character count.
    # A partial codon containing one or two ! symbols is one frameshifting codon.
    frameshift_count = sum(1 for codon in codons if "!" in codon)

    # Determine terminal biological codon after ignoring trailing codon-sized
    # alignment gaps. Do not strip partial frameshift codons; they carry signal.
    non_gap_indices = [i for i, codon in enumerate(codons) if codon != "---"]
    terminal_stop = False
    internal_stop_count = 0

    if non_gap_indices:
        last_biological_index = non_gap_indices[-1]
        terminal_stop = codons[last_biological_index] in STOP_CODONS
        for i, codon in enumerate(codons):
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


def validate_alignment_lengths(records: Dict[str, str], path: Path) -> None:
    """A MACSE alignment must have one common aligned length."""
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) > 1:
        details = Counter(len(sequence) for sequence in records.values())
        formatted = ", ".join(
            f"{length}:{count}" for length, count in sorted(details.items())
        )
        raise QCValidationError(
            f"MACSE FASTA does not have a single aligned sequence length: {path}. "
            f"length:count = {formatted}"
        )


def write_tsv_atomic(rows: Sequence[Dict[str, object]], output_path: Path) -> None:
    """Write only after all validation has succeeded, then atomically replace."""
    if not rows:
        raise QCValidationError("No QC rows were generated; refusing to write empty TSV.")

    temp_path = output_path.with_suffix(output_path.suffix + ".tmp")
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
                missing = [column for column in OUTPUT_COLUMNS if column not in row]
                extra = [column for column in row if column not in OUTPUT_COLUMNS]
                raise QCValidationError(
                    f"Internal output-column mismatch. Missing={missing}; extra={extra}"
                )
            writer.writerow(row)

    temp_path.replace(output_path)


def print_startup(
    genomes_dir: Path,
    project_file: Path,
    cohort_label: str,
    shared_genes_path: Path,
    pre_align_dir: Path,
    macse_dir: Path,
    output_dir: Path,
    output_path: Path,
) -> None:
    print("=" * 78)
    print("ORTHOLOG QUALITY CONTROL")
    print("=" * 78)
    print(f"QC version:                 {QC_VERSION}")
    print(f"Genomes directory:          {genomes_dir}")
    print(f"Manifest/accessions file:   {project_file}")
    print(f"Project/cohort label:       {cohort_label}")
    print(f"Shared genes source:        {shared_genes_path}")
    print(f"Raw pre-alignment dir:      {pre_align_dir}")
    print(f"MACSE alignment dir:        {macse_dir}")
    print(f"Output directory:           {output_dir}")
    print(f"Final TSV:                  {output_path}")
    print(f"Genetic code:               NCBI {GENETIC_CODE_ID} (Standard Code)")
    print(
        f"Max raw ambiguity fraction: {MAX_RAW_AMBIGUOUS_FRACTION} "
        "(fail only if greater)"
    )
    print(
        f"Max raw invalid fraction:   {MAX_RAW_INVALID_FRACTION} "
        "(fail only if greater)"
    )
    print("=" * 78)


def print_validation_issues(issues: Sequence[str]) -> None:
    print("\nMISSING / UNMATCHED RAW-MACSE RECORDS")
    print("-" * 78)
    if not issues:
        print("None")
        return
    for issue in issues:
        print(f"- {issue}")


def print_summary(rows: Sequence[Dict[str, object]], validation_issues: Sequence[str]) -> None:
    raw_counts = Counter(str(row["raw_cds_qc"]) for row in rows)
    macse_counts = Counter(str(row["macse_qc"]) for row in rows)
    final_counts = Counter(str(row["final_qc"]) for row in rows)

    accessions = {str(row["accession"]) for row in rows}
    orthologs = {str(row["odb12_id"]) for row in rows}

    print("\n" + "=" * 78)
    print("QC COMPLETE")
    print("=" * 78)
    print(f"Total unique accession/ortholog records evaluated: {len(rows)}")
    print(f"Unique accessions:                              {len(accessions)}")
    print(f"Unique ODB12 orthologs:                         {len(orthologs)}")
    print("\nRaw CDS QC")
    print(f"  pass: {raw_counts.get('pass', 0)}")
    print(f"  fail: {raw_counts.get('fail', 0)}")
    print("\nMACSE QC")
    print(f"  pass: {macse_counts.get('pass', 0)}")
    print(f"  fail: {macse_counts.get('fail', 0)}")
    print("\nFinal QC")
    for category in FINAL_QC_ORDER:
        print(f"  {category}: {final_counts.get(category, 0)}")

    rescued = final_counts.get("rescued", 0)
    rescued_pct = (100.0 * rescued / len(rows)) if rows else 0.0
    print(f"\nRescued records: {rescued}/{len(rows)} ({rescued_pct:.2f}%)")
    print_validation_issues(validation_issues)
    print("=" * 78)


def build_row(
    accession: str,
    odb12_id: str,
    raw_metrics: Dict[str, object],
    macse_metrics: Dict[str, object],
) -> Dict[str, object]:
    raw_qc = str(raw_metrics["raw_cds_qc"])
    macse_qc = str(macse_metrics["macse_qc"])
    final_qc = FINAL_QC_MAP[(raw_qc, macse_qc)]

    # Construct explicitly in the exact requested output order.
    return {
        "accession": accession,
        "odb12_id": odb12_id,
        "final_qc": final_qc,
        "macse_qc": macse_qc,
        "raw_cds_qc": raw_qc,
        "raw_length": raw_metrics["raw_length"],
        "raw_length_mod_3": raw_metrics["raw_length_mod_3"],
        "raw_internal_stop_count": raw_metrics["raw_internal_stop_count"],
        "raw_terminal_stop": raw_metrics["raw_terminal_stop"],
        "raw_ambiguous_fraction": raw_metrics["raw_ambiguous_fraction"],
        "raw_invalid_fraction": raw_metrics["raw_invalid_fraction"],
        "raw_failure_reasons": raw_metrics["raw_failure_reasons"],
        "macse_frameshift": macse_metrics["macse_frameshift"],
        "macse_internal_stop": macse_metrics["macse_internal_stop"],
        "macse_terminal_stop": macse_metrics["macse_terminal_stop"],
        "macse_failure_reasons": macse_metrics["macse_failure_reasons"],
        "qc_version": QC_VERSION,
    }


def run_qc(args: argparse.Namespace) -> Path:
    genomes_dir = Path(args.genomes_dir).expanduser().resolve()
    project_file = Path(args.accessions).expanduser().resolve()

    if not genomes_dir.is_dir():
        raise NotADirectoryError(f"Genomes directory not found: {genomes_dir}")
    if not project_file.is_file():
        raise FileNotFoundError(
            f"Project manifest/accessions file not found: {project_file}"
        )

    cohort_label = sanitize_label(project_file.stem)
    if not cohort_label:
        raise QCValidationError(
            f"Manifest/accessions filename produced an empty cohort label: {project_file.name}"
        )

    pre_align_dir = (
        Path(args.pre_alignment_dir).expanduser().resolve()
        if args.pre_alignment_dir
        else genomes_dir / "records/compleasm/alignments/01_pre_alignments"
    )
    macse_dir = (
        Path(args.macse_alignment_dir).expanduser().resolve()
        if args.macse_alignment_dir
        else genomes_dir / "records/compleasm/alignments/02_macse_alignments"
    )

    if not pre_align_dir.is_dir():
        raise NotADirectoryError(f"Raw pre-alignment directory not found: {pre_align_dir}")
    if not macse_dir.is_dir():
        raise NotADirectoryError(f"MACSE alignment directory not found: {macse_dir}")

    shared_genes_path = resolve_shared_genes_path(
        genomes_dir,
        cohort_label,
        Path(args.shared_genes).expanduser() if args.shared_genes else None,
    )
    ortholog_ids = load_ortholog_ids(shared_genes_path)

    project_accessions = read_project_accessions(project_file)
    allowed_roots = {root_acc(accession) for accession in project_accessions}

    output_dir = (
        genomes_dir
        / "records/compleasm/records/ortholog_validity"
        / cohort_label
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "ortholog_validity.tsv"

    print_startup(
        genomes_dir=genomes_dir,
        project_file=project_file,
        cohort_label=cohort_label,
        shared_genes_path=shared_genes_path,
        pre_align_dir=pre_align_dir,
        macse_dir=macse_dir,
        output_dir=output_dir,
        output_path=output_path,
    )
    print(f"Loaded project accessions:  {len(project_accessions)}")
    print(f"Orthologs requested:        {len(ortholog_ids)}")

    rows: List[Dict[str, object]] = []
    validation_issues: List[str] = []
    seen_pairs: Set[Tuple[str, str]] = set()

    for odb12_id in ortholog_ids:
        raw_path = pre_align_dir / f"{odb12_id}_unaligned.fasta"
        macse_path = macse_dir / f"{odb12_id}_NT.fasta"

        raw_exists = raw_path.is_file()
        macse_exists = macse_path.is_file()

        if not raw_exists:
            validation_issues.append(
                f"odb12_id={odb12_id}: missing raw FASTA: {raw_path}"
            )
        if not macse_exists:
            validation_issues.append(
                f"odb12_id={odb12_id}: missing MACSE FASTA: {macse_path}"
            )
        if not raw_exists or not macse_exists:
            continue

        raw_records = parse_fasta_strict(raw_path)
        macse_records = parse_fasta_strict(macse_path)

        if not raw_records:
            validation_issues.append(
                f"odb12_id={odb12_id}: raw FASTA contains no records: {raw_path}"
            )
            continue
        if not macse_records:
            validation_issues.append(
                f"odb12_id={odb12_id}: MACSE FASTA contains no records: {macse_path}"
            )
            continue

        validate_alignment_lengths(macse_records, macse_path)

        # Raw FASTAs were rewritten for the requested cohort by the upstream
        # script. Reject stale/foreign accession roots rather than silently
        # incorporating them into this project-specific QC table.
        for accession in raw_records:
            if root_acc(accession) not in allowed_roots:
                validation_issues.append(
                    f"odb12_id={odb12_id}: raw accession '{accession}' is outside "
                    f"the requested manifest/accessions cohort"
                )

        raw_ids = set(raw_records)
        macse_ids = set(macse_records)
        raw_only = sorted(raw_ids - macse_ids)
        macse_only = sorted(macse_ids - raw_ids)

        for accession in raw_only:
            validation_issues.append(
                f"odb12_id={odb12_id}: accession '{accession}' present in raw "
                "FASTA but missing from MACSE FASTA"
            )
        for accession in macse_only:
            validation_issues.append(
                f"odb12_id={odb12_id}: accession '{accession}' present in MACSE "
                "FASTA but missing from raw FASTA"
            )

        # Evaluate only exact accession matches. Any mismatch prevents final TSV
        # writing later, so no incomplete/fifth final_qc category is invented.
        for accession in sorted(raw_ids & macse_ids):
            pair = (accession, odb12_id)
            if pair in seen_pairs:
                raise QCValidationError(
                    f"Duplicate accession + odb12_id pair generated: {pair}"
                )
            seen_pairs.add(pair)

            raw_metrics = evaluate_raw_sequence(raw_records[accession])
            macse_metrics = evaluate_macse_sequence(
                macse_records[accession], macse_path, accession
            )
            rows.append(
                build_row(
                    accession=accession,
                    odb12_id=odb12_id,
                    raw_metrics=raw_metrics,
                    macse_metrics=macse_metrics,
                )
            )

    # Explicit final uniqueness validation, independent of the incremental guard.
    final_pairs = [(str(row["accession"]), str(row["odb12_id"])) for row in rows]
    if len(final_pairs) != len(set(final_pairs)):
        duplicates = [
            pair for pair, count in Counter(final_pairs).items() if count > 1
        ]
        raise QCValidationError(
            f"Duplicate accession + odb12_id combinations in final records: {duplicates}"
        )

    if validation_issues:
        print_validation_issues(validation_issues)
        raise QCValidationError(
            f"Found {len(validation_issues)} missing/unmatched/cohort validation "
            "problem(s). No ortholog_validity.tsv was written or overwritten."
        )

    rows.sort(key=lambda row: (str(row["odb12_id"]), str(row["accession"])))
    write_tsv_atomic(rows, output_path)
    print_summary(rows, validation_issues)
    print(f"\nWrote: {output_path}")
    return output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate independent raw-CDS and pre-cleaner MACSE QC for a "
            "manifest-defined ortholog cohort, then assign versioned final QC."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "genomes_dir",
        help="Root genomes directory, e.g. /Users/rossoaa/projects/genomes",
    )
    parser.add_argument(
        "accessions",
        help=(
            "Same project manifest CSV (requires accession column) or accession "
            "TXT file used by 04_prealign_to_macse_v3.py"
        ),
    )
    parser.add_argument(
        "--pre-alignment-dir",
        dest="pre_alignment_dir",
        default=None,
        help=(
            "Override raw pre-alignment directory. Default: "
            "<genomes_dir>/records/compleasm/alignments/01_pre_alignments"
        ),
    )
    parser.add_argument(
        "--macse-alignment-dir",
        dest="macse_alignment_dir",
        default=None,
        help=(
            "Override MACSE NT alignment directory. Default: "
            "<genomes_dir>/records/compleasm/alignments/02_macse_alignments"
        ),
    )
    parser.add_argument(
        "--shared-genes",
        default=None,
        help=(
            "Optional explicit shared-single-genes file. If omitted, use the "
            "same project-specific then shared_single_genes.csv fallback logic "
            "as 04_prealign_to_macse_v3.py."
        ),
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        run_qc(args)
    except (QCValidationError, FileNotFoundError, NotADirectoryError, OSError) as exc:
        print(f"\n[ERROR] {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
