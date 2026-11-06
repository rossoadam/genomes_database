#!/usr/bin/env python3
"""03_raw_ortholog_validity_batch_v4.py

Expand and audit per-genome Compleasm raw CDS QC independent of alignment occupancy.

Unlike v3's separate raw_ortholog_validity output, v4 audits and conservatively
updates the EXISTING 17-column ortholog_validity.tsv produced by
04b_qualtiy_control_v2.py. Reuse that script's evaluate_raw_sequence() through
an explicitly supplied --qc-script; do not recode its QC criteria.

Changes to existing records are restricted to filling missing raw QC fields.
Populated raw values, MACSE values, final_qc, and qc_version are never changed.
Every normal run gets independent diagnostics, and any actual TSV update is
backed up and atomically replaced. Dry runs write nothing.

Example:
  python 03_raw_ortholog_validity_batch_v4.py /Users/rossoaa/projects/genomes \\
    --manifest /Users/rossoaa/projects/genomes/project_manifests/body-size-dna-dynamics_manifest.csv \\
    --qc-script /path/to/04b_qualtiy_control_v2.py --dry-run

Requires Python >= 3.9 and the dependencies of the external QC script (BioPython).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import math
import os
import re
import shutil
import stat
import sys
import tempfile
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

SCRIPT_VERSION = "v4"
SQL_NULL = r"\N"
TRIPLOID_ACCESSION = "GWHFILM00000000.1"
OUTPUT_COLUMNS = [
    "accession", "odb12_id", "final_qc", "macse_qc", "raw_cds_qc",
    "raw_length", "raw_length_mod_3", "raw_internal_stop_count",
    "raw_terminal_stop", "raw_ambiguous_fraction", "raw_invalid_fraction",
    "raw_failure_reasons", "macse_frameshift", "macse_internal_stop",
    "macse_terminal_stop", "macse_failure_reasons", "qc_version",
]
RAW_FIELDS = [
    "raw_cds_qc", "raw_length", "raw_length_mod_3",
    "raw_internal_stop_count", "raw_terminal_stop", "raw_ambiguous_fraction",
    "raw_invalid_fraction", "raw_failure_reasons",
]
MACSE_FIELDS = ["final_qc", "macse_qc", "macse_frameshift",
                "macse_internal_stop", "macse_terminal_stop", "macse_failure_reasons"]
INTEGER_FIELDS = {"raw_length", "raw_length_mod_3", "raw_internal_stop_count"}
FRACTION_FIELDS = {"raw_ambiguous_fraction", "raw_invalid_fraction"}
ACCEPTED_STATUSES = {"single", "duplicated", "duplicate", "fragmented", "missing"}
SOURCE_COLUMNS = ["accession", "odb12_id", "compleasm_status", "source_path",
                  "problem_category", "details", "action_taken"]
NEW_COLUMNS = OUTPUT_COLUMNS + ["source_cds_fasta", "compleasm_status"]
FILL_COLUMNS = ["accession", "odb12_id", "field", "existing_value", "new_value", "source_cds_fasta"]
DISCREPANCY_COLUMNS = ["accession", "odb12_id", "field", "existing_value",
                       "recalculated_value", "discrepancy_category", "source_cds_fasta", "qc_version"]
UNMATCHED_COLUMNS = ["accession", "odb12_id", "reason"]
STATUS_COLUMNS = ["accession", "status", "count", "eligibility"]
Key = Tuple[str, str]


class AuditError(RuntimeError):
    """Unsafe configuration, prior data, or output validation."""


def is_null(value: Any) -> bool:
    return value is None or str(value).strip().lower() in {"", r"\n", "na", "n/a", "nan", "none", "null"}


def string_value(value: Any) -> str:
    return SQL_NULL if value is None else str(value)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for data in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(data)
    return h.hexdigest()


def read_manifest(path: Path) -> List[str]:
    if not path.is_file():
        raise AuditError(f"Manifest not found: {path}")
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None or "accession" not in reader.fieldnames:
            raise AuditError(f"Manifest requires an exact 'accession' column: {path}")
        accessions: List[str] = []
        for line_no, row in enumerate(reader, 2):
            if None in row:
                raise AuditError(f"Malformed manifest record near line {line_no}")
            accession = (row.get("accession") or "").strip()
            if not accession:
                raise AuditError(f"Missing accession in manifest line {line_no}")
            accessions.append(accession)
    if not accessions:
        raise AuditError("Manifest is empty")
    duplicates = sorted(k for k, n in Counter(accessions).items() if n > 1)
    if duplicates:
        raise AuditError(f"Duplicate manifest accessions: {duplicates[:10]}")
    return accessions


def read_existing(path: Path) -> List[Dict[str, str]]:
    if not path.is_file():
        raise AuditError(f"Existing ortholog_validity.tsv is required: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != OUTPUT_COLUMNS:
            raise AuditError(f"Unexpected QC TSV schema in {path}.\n"
                             f"Expected: {OUTPUT_COLUMNS}\nFound: {reader.fieldnames}")
        rows = []
        seen: Set[Key] = set()
        duplicates: List[Key] = []
        for line_no, row in enumerate(reader, 2):
            if None in row or any(v is None for v in row.values()):
                raise AuditError(f"Malformed TSV row near line {line_no}")
            key = ((row["accession"] or "").strip(), (row["odb12_id"] or "").strip())
            if not all(key) or is_null(key[0]) or is_null(key[1]):
                raise AuditError(f"Invalid accession/odb12_id in row {line_no}")
            if key in seen:
                duplicates.append(key)
            seen.add(key)
            rows.append(row)
    if duplicates:
        raise AuditError(f"Duplicate keys in existing TSV; no changes permitted: {duplicates[:12]}")
    return rows


def load_qc_module(path: Path):
    if not path.is_file():
        raise AuditError(f"--qc-script not found: {path}")
    # Prevent imports from creating __pycache__ in the user's project,
    # particularly important for strictly read-only --dry-run execution.
    previous_bytecode_setting = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        try:
            spec = importlib.util.spec_from_file_location("ortholog_qc_v2_external", path)
            if spec is None or spec.loader is None:
                raise AuditError(f"Cannot import QC implementation: {path}")
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
        except (ImportError, SyntaxError, OSError, RuntimeError) as exc:
            raise AuditError(f"Cannot load --qc-script {path}: {exc}") from exc
    finally:
        sys.dont_write_bytecode = previous_bytecode_setting
    expected = {
        "QC_VERSION": "v2", "GENETIC_CODE_ID": 1,
        "MAX_RAW_AMBIGUOUS_FRACTION": 0.05,
        "MAX_RAW_INVALID_FRACTION": 0.0, "SQL_NULL": SQL_NULL,
    }
    for field, value in expected.items():
        if getattr(module, field, None) != value:
            raise AuditError(f"Incompatible QC module {field}: "
                             f"expected {value!r}, found {getattr(module, field, None)!r}")
    if getattr(module, "OUTPUT_COLUMNS", None) != OUTPUT_COLUMNS:
        raise AuditError("The QC module's 17-column schema differs from the expected v2 schema")
    if not callable(getattr(module, "evaluate_raw_sequence", None)):
        raise AuditError("QC module must provide evaluate_raw_sequence(sequence)")
    # Smoke-test core v2 rules against supplied implementation; no local QC replacement.
    probes = [("ATGAAATAA", "pass", 0), ("ATGTAGAAA", "fail", 0),
              ("ATGA", "pass", 1), ("ATGNNNNNN", "fail", 0)]
    for seq, passed, mod3 in probes:
        result = module.evaluate_raw_sequence(seq)
        if result.get("raw_cds_qc") != passed or result.get("raw_length_mod_3") != mod3:
            raise AuditError(f"Unexpected QC behavior in supplied module for {seq!r}")
    return module


def resolve_compleasm_dir(root: Path, accession: str) -> Tuple[Optional[Path], str]:
    folder = root / "records" / "compleasm"
    if not folder.is_dir():
        return None, f"Compleasm root directory not found: {folder}"
    candidates = sorted(p for p in folder.glob(f"{glob_escape(accession)}__*") if p.is_dir()
                        and p.name.split("__", 1)[0] == accession)
    if len(candidates) == 1:
        return candidates[0] / "sauropsida_odb12", ""
    if not candidates:
        return None, f"No exact-accession Compleasm directory for {accession} in {folder}"
    return None, f"Ambiguous directories for {accession}: {', '.join(map(str, candidates))}"


def glob_escape(value: str) -> str:
    """Escape glob metacharacters in an accession (without broadening matching)."""
    return re.sub(r"([*?\[\]])", lambda m: "[" + m.group(1) + "]", value)


def status_normalized(value: str) -> str:
    status = str(value).strip().lower()
    return "duplicated" if status == "duplicate" else status


def eligible_status(accession: str, status: str) -> bool:
    return status == "single" or (status == "duplicated" and accession == TRIPLOID_ACCESSION)


def parse_full_table(path: Path, accession: str, diagnostics: List[Tuple[str, ...]]):
    """Return {odb12_id: (status, row_count)} or raise for unreadable/unsupported table.

    Compleasm files can have a commented header (# Gene\tStatus...) or a bare
    header, or be headerless with Gene/Status as the first two fields. Explicit
    header names take precedence; positional fallback requires valid statuses.
    """
    records: Dict[str, List[str]] = defaultdict(list)
    header: Optional[List[str]] = None
    found_data = False
    recognized_status = 0
    gene_i = status_i = None
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for line_no, line in enumerate(handle, 1):
            text = line.rstrip("\r\n")
            if not text.strip():
                continue
            possible = text.lstrip("#").strip() if text.lstrip().startswith("#") else text
            fields = [f.strip() for f in possible.split("\t")]
            lower = [re.sub(r"[\s_#-]+", "", f.lower()) for f in fields]
            is_header = ("status" in lower and any(
                c in lower for c in ("gene", "buscoid", "orthologid", "odb12id", "orthodbid")))
            if is_header:
                if found_data or header is not None:
                    raise AuditError(f"Repeated or late header in {path}:{line_no}")
                header = fields
                gene_i = next(i for i, c in enumerate(lower) if c in
                              ("gene", "buscoid", "orthologid", "odb12id", "orthodbid"))
                status_i = lower.index("status")
                continue
            if text.lstrip().startswith("#"):
                continue
            found_data = True
            if header is None:
                if len(fields) < 2:
                    diagnostics.append((accession, fields[0] if fields else "", "", str(path),
                                        "malformed_full_table_row", f"Line {line_no}: fewer than 2 fields", "skip_row"))
                    continue
                gene_i, status_i = 0, 1
            if gene_i is None or status_i is None or max(gene_i, status_i) >= len(fields):
                diagnostics.append((accession, fields[0] if fields else "", "", str(path),
                                    "malformed_full_table_row", f"Line {line_no}: missing Gene/Status field", "skip_row"))
                continue
            gene = fields[gene_i]
            status = status_normalized(fields[status_i])
            if not gene or not status:
                diagnostics.append((accession, gene, status, str(path),
                                    "malformed_full_table_row", f"Line {line_no}: empty gene or status", "skip_row"))
                continue
            records[gene].append(status)
            if status in ACCEPTED_STATUSES:
                recognized_status += 1
    if not records or recognized_status == 0:
        raise AuditError(f"Cannot interpret Compleasm Gene/Status records: {path}")
    result: Dict[str, Tuple[str, int]] = {}
    for gene, statuses in records.items():
        unique = sorted(set(statuses))
        if len(unique) > 1:
            diagnostics.append((accession, gene, ";".join(unique), str(path),
                                "conflicting_full_table_status", "Same ortholog has differing statuses",
                                "skip_ambiguous_ortholog"))
            continue
        status = unique[0]
        if status not in ACCEPTED_STATUSES:
            diagnostics.append((accession, gene, status, str(path), "unknown_compleasm_status",
                                "Not one of Single/Duplicated/Fragmented/Missing", "exclude"))
            continue
        if len(statuses) > 1:
            diagnostics.append((accession, gene, status, str(path),
                                "repeated_full_table_entries",
                                f"{len(statuses)} rows for same ODB12 ID; single FASTA record required",
                                "assess_once_if_eligible"))
        result[gene] = (status, len(statuses))
    return result


def select_fasta_file(path: Path) -> Tuple[Optional[Path], str]:
    candidates = sorted(p for pattern in ("*_cds_compleasm.fasta", "*_cds_compleasm.fa")
                        for p in path.glob(pattern) if p.is_file())
    if len(candidates) == 1:
        return candidates[0], ""
    if not candidates:
        return None, f"No *_cds_compleasm.fasta file found in {path}"
    return None, f"Multiple matching CDS FASTA files: {', '.join(map(str, candidates))}"


def read_target_fasta(path: Path, target_ids: Set[str], accession: str,
                      diagnostics: List[Tuple[str, ...]]) -> Dict[str, str]:
    """Read target records without losing valid records to isolated FASTA errors.

    Duplicated first-token IDs are always rejected; do not select first/last.
    Only target sequences are retained in memory (the full file is scanned).
    """
    sequences: Dict[str, str] = {}
    seen: Set[str] = set()
    duplicate: Set[str] = set()
    current: Optional[str] = None
    content: List[str] = []

    def finish():
        nonlocal current, content
        if current is None or current not in target_ids or current in duplicate:
            return
        seq = "".join(content).upper()
        if not seq:
            diagnostics.append((accession, current, "", str(path), "empty_cds_sequence",
                                "FASTA record has no sequence", "skip_ortholog"))
        else:
            sequences[current] = seq

    with path.open("r", encoding="utf-8-sig") as handle:
        for lineno, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                finish()
                tokens = stripped[1:].strip().split()
                if not tokens:
                    diagnostics.append((accession, "", "", str(path), "empty_fasta_header",
                                        f"Line {lineno}", "skip_record"))
                    current = None
                else:
                    current = tokens[0]
                    if current in seen:
                        duplicate.add(current)
                        sequences.pop(current, None)
                        diagnostics.append((accession, current, "", str(path), "duplicate_fasta_identifier",
                                            f"Line {lineno}: repeated FASTA identifier", "skip_ortholog"))
                    seen.add(current)
                content = []
            elif current is None:
                diagnostics.append((accession, "", "", str(path), "fasta_sequence_without_header",
                                    f"Line {lineno}: sequence before a valid header", "skip_line"))
            elif current in target_ids and current not in duplicate:
                content.append(re.sub(r"\s+", "", stripped))
    finish()
    for key in duplicate:
        sequences.pop(key, None)
    return sequences


def normalize_for_compare(field: str, value: str):
    v = str(value).strip()
    if field in INTEGER_FIELDS:
        return int(v)
    if field in FRACTION_FIELDS:
        result = float(v)
        if not math.isfinite(result):
            raise ValueError(v)
        return result
    if field == "raw_terminal_stop":
        lowered = v.lower()
        if lowered in {"true", "t", "1", "yes"}:
            return True
        if lowered in {"false", "f", "0", "no"}:
            return False
        raise ValueError(v)
    if field == "raw_failure_reasons":
        return tuple(sorted(s.strip().lower() for s in v.split(";") if s.strip()))
    return v.lower()


def values_agree(field: str, old: str, recalculated: str) -> bool:
    try:
        a = normalize_for_compare(field, old)
        b = normalize_for_compare(field, recalculated)
    except (ValueError, OverflowError):
        return old == recalculated
    if field in FRACTION_FIELDS:
        return math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-12)
    return a == b


def write_tsv(path: Path, headers: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=headers, delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({h: row[h] for h in headers})
        stream.flush()
        os.fsync(stream.fileno())


def validate_update(original: List[Dict[str, str]], updated: List[Dict[str, str]],
                    new_records: List[Dict[str, str]], fills: List[Dict[str, str]],
                    allowed_accessions: Set[str]) -> None:
    original_keys = [(r["accession"], r["odb12_id"]) for r in original]
    new_keys = [(r["accession"], r["odb12_id"]) for r in new_records]
    updated_keys = [(r["accession"], r["odb12_id"]) for r in updated]
    if len(updated_keys) != len(set(updated_keys)):
        raise AuditError("Updated table contains duplicate accession+odb12_id keys")
    if len(updated) != len(original) + len(new_records) or updated_keys[:len(original)] != original_keys:
        raise AuditError("Historical record order/coverage changed unexpectedly")
    if set(new_keys) & set(original_keys) or len(new_keys) != len(set(new_keys)):
        raise AuditError("New records collide with existing records or one another")
    fills_by_key = {(r["accession"], r["odb12_id"], r["field"]): r for r in fills}
    for old, new in zip(original, updated):
        if set(new) != set(OUTPUT_COLUMNS):
            raise AuditError("Updated existing record does not have 17 columns")
        key = (old["accession"], old["odb12_id"])
        for field in OUTPUT_COLUMNS:
            if old[field] == new[field]:
                continue
            if field not in RAW_FIELDS or not is_null(old[field]):
                raise AuditError(f"Forbidden historical-field edit: {key}: {field}")
            if (key[0], key[1], field) not in fills_by_key:
                raise AuditError(f"Historical-field edit without provenance: {key}: {field}")
    for new in new_records:
        if list(new) != OUTPUT_COLUMNS:
            raise AuditError("New record schema does not match 17-column schema")
        if new["accession"] not in allowed_accessions:
            raise AuditError("New record accession not in manifest")
        if new["raw_cds_qc"] not in {"pass", "fail"}:
            raise AuditError("Invalid raw QC state in new record")
        if new["qc_version"] != "v2":
            raise AuditError("Unexpected QC version in new record")
        if any(new[field] != SQL_NULL for field in MACSE_FIELDS):
            raise AuditError("MACSE fields must be SQL NULL in new records")


def print_summary(stats: Mapping[str, Any], additions: Sequence[Mapping[str, str]],
                  fills: Sequence[Mapping[str, str]], discrepancies: Sequence[Mapping[str, str]],
                  dry_run: bool, backup: Optional[Path], diag_dir: Optional[Path]):
    print("\n" + "=" * 76)
    outcome = "DRY RUN" if dry_run else (
        "FINISHED — PARTIAL SOURCE COVERAGE"
        if int(stats.get("unassessable_eligible_records", 0)) > 0
        or int(stats.get("genomes_with_readable_compleasm_and_fasta", 0)) < int(stats.get("manifest_genomes", 0))
        else "COMPLETE"
    )
    print("RAW COMPLEASM QC EXPANSION — " + outcome)
    print("=" * 76)
    for k, v in stats.items():
        print(f"{k:42} {v}")
    if dry_run:
        for name, rows in (("Proposed new records", additions),
                           ("Proposed missing-field completions", fills),
                           ("Raw QC discrepancies", discrepancies)):
            print(f"\n{name} (first 8):")
            for row in rows[:8]:
                print("  " + ", ".join(f"{k}={str(v)[:55]}" for k, v in row.items()
                                      if k in {"accession", "odb12_id", "field", "existing_value", "new_value", "recalculated_value", "raw_cds_qc"}))
    if backup:
        print(f"Backup: {backup}")
    if diag_dir:
        print(f"Diagnostics: {diag_dir}")
    if dry_run:
        print("No files or directories were modified.")
    print("=" * 76)


def unique_run_dir(parent: Path) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    for n in range(100):
        candidate = parent / (stamp if n == 0 else f"{stamp}_{n}")
        try:
            candidate.mkdir(parents=False, exist_ok=False)
            return candidate
        except FileExistsError:
            continue
    raise AuditError("Could not create unique execution diagnostics directory")


def update_file_atomic(path: Path, updated: List[Dict[str, str]], original_digest: str,
                       run_dir: Path) -> Path:
    if digest(path) != original_digest:
        raise AuditError("ortholog_validity.tsv changed during analysis; refusing replacement")
    backup = run_dir / f"ortholog_validity_before_v4_{run_dir.name}.tsv"
    shutil.copy2(path, backup)
    if digest(backup) != original_digest:
        raise AuditError("Backup did not match original TSV; refusing replacement")
    temp_name: Optional[str] = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", prefix=".ortholog_validity_v4_",
                                         suffix=".tmp", dir=path.parent, delete=False,
                                         encoding="utf-8", newline="") as f:
            temp_name = f.name
            writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS, delimiter="\t",
                                    lineterminator="\n", extrasaction="raise")
            writer.writeheader()
            for row in updated:
                writer.writerow(row)
            f.flush()
            os.fsync(f.fileno())
        staged = Path(temp_name)
        # Atomic replacement should not unexpectedly change TSV permissions.
        os.chmod(staged, stat.S_IMODE(path.stat().st_mode))
        test_rows = read_existing(staged)
        if test_rows != updated:
            raise AuditError("Staged TSV roundtrip verification failed")
        if digest(path) != original_digest:
            raise AuditError("Original TSV changed during staging; refusing replacement")
        os.replace(staged, path)
        temp_name = None
    finally:
        if temp_name and os.path.exists(temp_name):
            os.unlink(temp_name)
    return backup


def run(args: argparse.Namespace) -> None:
    genomes_dir = Path(args.genomes_dir).expanduser().resolve()
    manifest = Path(args.manifest).expanduser().resolve()
    qc_script = Path(args.qc_script).expanduser().resolve()
    if not genomes_dir.is_dir():
        raise AuditError(f"Genomes directory not found: {genomes_dir}")
    accessions = read_manifest(manifest)
    cohort = re.sub(r"[^A-Za-z0-9._-]+", "_", manifest.stem).strip("_")
    output_path = (Path(args.output_tsv).expanduser().resolve() if args.output_tsv else
                   genomes_dir / "records/compleasm/records/ortholog_validity" / cohort /
                   "ortholog_validity.tsv")
    if not output_path.is_file():
        raise AuditError(f"Expected existing TSV was not found: {output_path}")
    qc = load_qc_module(qc_script)
    existing = read_existing(output_path)
    original_hash = digest(output_path)
    updated = [r.copy() for r in existing]
    existing_lookup = {(r["accession"], r["odb12_id"]): i for i, r in enumerate(existing)}
    additions: List[Dict[str, str]] = []
    addition_diagnostics: List[Dict[str, str]] = []
    fills: List[Dict[str, str]] = []
    discrepancies: List[Dict[str, str]] = []
    source_diagnostics: List[Tuple[str, ...]] = []
    genome_status: Counter = Counter()
    exclusion_count = Counter()
    verified: Set[Key] = set()
    discovered_eligible: Set[Key] = set()
    unassessable: Set[Key] = set()
    unverified_reason: Dict[Key, str] = {}
    genome_errors: Dict[str, str] = {}
    genome_dirs = 0
    scanned_genomes = 0
    qc_pass = 0
    qc_fail = 0

    print(f"Manifest: {manifest} ({len(accessions)} accessions)")
    print(f"Existing QC TSV: {output_path} ({len(existing)} records)")
    print(f"Raw CDS QC implementation: {qc_script} (v2; SHA256 {digest(qc_script)[:16]}…)")
    print("Scanning original Compleasm output without an occupancy threshold...")

    for genome_no, accession in enumerate(accessions, 1):
        folder, reason = resolve_compleasm_dir(genomes_dir, accession)
        if folder is None:
            genome_errors[accession] = reason
            source_diagnostics.append((accession, "", "", "", "missing_or_ambiguous_genome_directory",
                                       reason, "skip_genome"))
            continue
        genome_dirs += 1
        full_path = folder / "full_table.tsv"
        if not full_path.is_file():
            reason = f"Missing full_table.tsv: {full_path}"
            genome_errors[accession] = reason
            source_diagnostics.append((accession, "", "", str(full_path), "missing_full_table",
                                       reason, "skip_genome"))
            continue
        try:
            full = parse_full_table(full_path, accession, source_diagnostics)
        except (AuditError, UnicodeError, OSError) as exc:
            reason = f"Unusable full_table.tsv: {exc}"
            genome_errors[accession] = reason
            source_diagnostics.append((accession, "", "", str(full_path), "invalid_full_table",
                                       reason, "skip_genome"))
            continue
        eligible: Dict[str, str] = {}
        for ortholog, (status, nrows) in full.items():
            genome_status[(accession, status)] += 1
            if eligible_status(accession, status):
                eligible[ortholog] = status
                discovered_eligible.add((accession, ortholog))
            else:
                exclusion_count[status] += 1
                source_diagnostics.append((accession, ortholog, status, str(full_path),
                                           "excluded_compleasm_status", "Not eligible under Single/Indotyphlops policy",
                                           "exclude"))
                unverified_reason[(accession, ortholog)] = f"Compleasm status excluded: {status}"
        fasta_path, fasta_problem = select_fasta_file(folder)
        if fasta_path is None:
            genome_errors[accession] = fasta_problem
            source_diagnostics.append((accession, "", "", str(folder), "missing_or_ambiguous_cds_fasta",
                                       fasta_problem, "skip_eligible_records"))
            for oid, status in eligible.items():
                key = (accession, oid)
                unassessable.add(key)
                unverified_reason[key] = fasta_problem
            continue
        scanned_genomes += 1
        try:
            sequences = read_target_fasta(fasta_path, set(eligible), accession, source_diagnostics)
        except (UnicodeError, OSError) as exc:
            reason = f"Unreadable CDS FASTA: {exc}"
            genome_errors[accession] = reason
            source_diagnostics.append((accession, "", "", str(fasta_path), "unreadable_cds_fasta",
                                       reason, "skip_eligible_records"))
            sequences = {}
        for oid, status in eligible.items():
            key = (accession, oid)
            seq = sequences.get(oid)
            if seq is None:
                unassessable.add(key)
                unverified_reason[key] = f"Eligible CDS sequence missing/unassessable in {fasta_path}"
                source_diagnostics.append((accession, oid, status, str(fasta_path),
                                           "unassessable_eligible_ortholog", unverified_reason[key], "skip_ortholog"))
                continue
            try:
                values = {field: string_value(v) for field, v in qc.evaluate_raw_sequence(seq).items()}
                if set(values) != set(RAW_FIELDS) or values["raw_cds_qc"] not in {"pass", "fail"}:
                    raise AuditError("External QC returned unexpected raw columns or classification")
            except Exception as exc:
                unassessable.add(key)
                unverified_reason[key] = f"QC function failed: {exc}"
                source_diagnostics.append((accession, oid, status, str(fasta_path),
                                           "raw_qc_evaluation_error", str(exc), "skip_ortholog"))
                continue
            verified.add(key)
            if values["raw_cds_qc"] == "pass":
                qc_pass += 1
            else:
                qc_fail += 1
            if key not in existing_lookup:
                new = {f: SQL_NULL for f in OUTPUT_COLUMNS}
                new.update({"accession": accession, "odb12_id": oid, "qc_version": qc.QC_VERSION})
                new.update(values)
                additions.append(new)
                addition_diagnostics.append({**new, "source_cds_fasta": str(fasta_path),
                                             "compleasm_status": status})
            else:
                row = updated[existing_lookup[key]]
                for field in RAW_FIELDS:
                    old = row[field]
                    fresh = values[field]
                    if is_null(old):
                        # A NULL failure reason is a REAL, populated result for a
                        # historically passing CDS. If fresh QC fails, that is an
                        # audit disagreement, not a field to backfill.
                        if field == "raw_failure_reasons":
                            historical_raw = existing[existing_lookup[key]]["raw_cds_qc"].strip().lower()
                            if historical_raw == "pass":
                                if fresh != SQL_NULL:
                                    discrepancies.append({
                                        "accession": accession, "odb12_id": oid,
                                        "field": field, "existing_value": old,
                                        "recalculated_value": fresh,
                                        "discrepancy_category": "historical_passing_null_reason_conflict",
                                        "source_cds_fasta": str(fasta_path),
                                        "qc_version": str(row["qc_version"]),
                                    })
                                continue
                        # A computed NULL does not add information to an empty field.
                        if fresh == SQL_NULL:
                            continue
                        row[field] = fresh
                        fills.append({"accession": accession, "odb12_id": oid,
                                      "field": field, "existing_value": old,
                                      "new_value": fresh, "source_cds_fasta": str(fasta_path)})
                    elif not values_agree(field, old, fresh):
                        discrepancies.append({"accession": accession, "odb12_id": oid,
                                              "field": field, "existing_value": old,
                                              "recalculated_value": fresh,
                                              "discrepancy_category": "populated_value_conflict",
                                              "source_cds_fasta": str(fasta_path),
                                              "qc_version": str(row["qc_version"])})
                # Historical MACSE/final classifications must not be rewritten,
                # even if they disagree with newly filled raw information.
                if row["final_qc"] in {"passing", "unsuitable", "rescued", "alignment-problematic"}:
                    raw = row["raw_cds_qc"].strip().lower()
                    macse = row["macse_qc"].strip().lower()
                    map_final = {("pass", "pass"): "passing", ("fail", "fail"): "unsuitable",
                                 ("fail", "pass"): "rescued", ("pass", "fail"): "alignment-problematic"}
                    expected_final = map_final.get((raw, macse))
                    if expected_final and expected_final != row["final_qc"]:
                        discrepancies.append({"accession": accession, "odb12_id": oid,
                                              "field": "final_qc", "existing_value": row["final_qc"],
                                              "recalculated_value": expected_final,
                                              "discrepancy_category": "historical_final_qc_inconsistency_preserved",
                                              "source_cds_fasta": str(fasta_path),
                                              "qc_version": str(row["qc_version"])})
        if genome_no % 10 == 0:
            print(f"  {genome_no}/{len(accessions)} accessions visited; {len(verified):,} CDS records assessed")

    if verified & unassessable or (discovered_eligible != (verified | unassessable)):
        raise AuditError("Internal completeness reconciliation failed for eligible orthologs")
    additions.sort(key=lambda row: (row["accession"], row["odb12_id"]))
    updated.extend(additions)
    validate_update(existing, updated, additions, fills, set(accessions))
    unmatched: List[Dict[str, str]] = []
    for row in existing:
        key = (row["accession"], row["odb12_id"])
        if key in verified:
            continue
        reason = unverified_reason.get(key) or genome_errors.get(row["accession"]) or (
            "Accession is outside current manifest" if row["accession"] not in set(accessions)
            else "ODB12 ID absent, excluded, or ambiguous in current Compleasm scan")
        unmatched.append({"accession": key[0], "odb12_id": key[1], "reason": reason})
    new_keys = {(r["accession"], r["odb12_id"]) for r in additions}
    if new_keys - verified:
        raise AuditError("New records were not verified against Compleasm")
    modified = bool(additions or fills)
    final_counts = Counter(r["final_qc"] for r in updated)
    raw_counts = Counter(r["raw_cds_qc"] for r in updated)
    stats: Dict[str, Any] = {
        "script_version": SCRIPT_VERSION,
        "execution_timestamp": datetime.now().isoformat(timespec="seconds"),
        "genomes_dir": str(genomes_dir), "manifest": str(manifest),
        "manifest_genomes": len(accessions), "genome_directories_found": genome_dirs,
        "genomes_with_readable_compleasm_and_fasta": scanned_genomes,
        "unique_eligible_accession_ortholog_records": len(discovered_eligible),
        "unique_eligible_ortholog_ids": len({o for _, o in discovered_eligible}),
        "excluded_records_by_compleasm_status": sum(exclusion_count.values()),
        "successfully_assessed_eligible_records": len(verified),
        "unassessable_eligible_records": len(unassessable),
        "recalculated_raw_qc_pass": qc_pass,
        "recalculated_raw_qc_fail": qc_fail,
        "original_qc_records": len(existing), "matched_existing_records": len(verified & set(existing_lookup)),
        "newly_added_records": len(additions),
        "existing_records_with_fields_completed": len({(r["accession"], r["odb12_id"]) for r in fills}),
        "individual_raw_fields_completed": len(fills),
        "existing_records_preserved_without_change": len(existing) - len({(r["accession"], r["odb12_id"]) for r in fills}),
        "discrepancies": len(discrepancies), "unmatched_existing_records": len(unmatched),
        "total_updated_qc_records": len(updated),
        "total_stored_raw_pass": raw_counts["pass"], "total_stored_raw_fail": raw_counts["fail"],
        "records_without_macse_assessment": sum(is_null(r["macse_qc"]) for r in updated),
        "final_qc_passing": final_counts["passing"], "final_qc_rescued": final_counts["rescued"],
        "final_qc_unsuitable": final_counts["unsuitable"],
        "final_qc_alignment_problematic": final_counts["alignment-problematic"],
        "qc_script": str(qc_script), "qc_script_sha256": digest(qc_script),
        "qc_version": qc.QC_VERSION, "genetic_code": qc.GENETIC_CODE_ID,
        "max_raw_ambiguous_fraction": qc.MAX_RAW_AMBIGUOUS_FRACTION,
        "max_raw_invalid_fraction": qc.MAX_RAW_INVALID_FRACTION,
        "input_tsv_sha256": original_hash,
        "tsv_changed": str(modified).lower(),
    }
    for status, count in sorted(exclusion_count.items()):
        stats[f"excluded_{status}"] = count

    backup: Optional[Path] = None
    diag_dir: Optional[Path] = None
    if not args.dry_run:
        diag_parent = output_path.parent / "qc_expansion_runs"
        diag_parent.mkdir(parents=True, exist_ok=True)
        diag_dir = unique_run_dir(diag_parent)
        write_tsv(diag_dir / "newly_added_records.tsv", NEW_COLUMNS, addition_diagnostics)
        write_tsv(diag_dir / "completed_raw_qc_fields.tsv", FILL_COLUMNS, fills)
        write_tsv(diag_dir / "raw_qc_discrepancies.tsv", DISCREPANCY_COLUMNS, discrepancies)
        write_tsv(diag_dir / "unmatched_existing_records.tsv", UNMATCHED_COLUMNS, unmatched)
        write_tsv(diag_dir / "source_data_diagnostics.tsv", SOURCE_COLUMNS,
                  (dict(zip(SOURCE_COLUMNS, row)) for row in source_diagnostics))
        write_tsv(diag_dir / "genome_status_summary.tsv", STATUS_COLUMNS,
                  ({"accession": acc, "status": status, "count": count,
                    "eligibility": "eligible" if eligible_status(acc, status) else "excluded"}
                   for (acc, status), count in sorted(genome_status.items())))
        if modified:
            backup = update_file_atomic(output_path, updated, original_hash, diag_dir)
        stats["output_tsv_sha256"] = digest(output_path)
        write_tsv(diag_dir / "run_summary.tsv", ["metric", "value"],
                  ({"metric": k, "value": str(v)} for k, v in stats.items()))
    else:
        if digest(output_path) != original_hash:
            raise AuditError("ortholog_validity.tsv changed during dry run")
    print_summary(stats, addition_diagnostics, fills, discrepancies, args.dry_run, backup, diag_dir)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("genomes_dir", help="Root genomes project directory")
    parser.add_argument("--manifest", required=True, help="Project manifest CSV/TSV with accession column")
    parser.add_argument("--qc-script", required=True, help="Path to 04b_qualtiy_control_v2.py (imported)")
    parser.add_argument("--output-tsv", help="Override existing ortholog_validity.tsv location")
    parser.add_argument("--dry-run", action="store_true", help="Validate and preview without writing files")
    return parser


def main() -> int:
    try:
        run(build_parser().parse_args())
    except (AuditError, OSError, UnicodeError, ValueError) as exc:
        print(f"\n[ERROR] {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
