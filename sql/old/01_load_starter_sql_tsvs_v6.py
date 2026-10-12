#!/usr/bin/env python3
"""
Load TSV outputs from 00b_build_starter_sql_tsvs_v15.py and
00c_build_tsvs_from_gff.py into a MySQL database using the
gc3_dynamics_v8_reduced schema.

This is intentionally standalone for troubleshooting. It supports both legacy
single-table TSVs and the chunked table-specific .tsv.gz layout retained by 00b
v15, loads tables in dependency order, and can optionally truncate the
starter/window tables before loading. If chunk directories are present for
high-volume tables, stale legacy combined TSVs for those same tables are ignored.

Version v6 updates:
  - Add a mandatory fast preflight stage before expensive full-data validation.
  - Report the active Python executable/environment and installed PyMySQL version.
  - Fail immediately if PyMySQL is unavailable.
  - Verify the TSV directory, required table inputs, readability, and TSV headers.
  - Test the MySQL connection before scanning all TSV data rows.
  - Verify all loader-supported MySQL tables and required columns are present.
  - Verify the client requests LOCAL INFILE and the MySQL server has local_infile enabled.
  - Add --skip-data-validation for repeat loads of an already-validated dataset.
  - Add --preflight-only to run only the fast checks and make no database changes.
  - Keep full TSV relationship/duplicate validation as the default and as the
    required behavior of --check-only.
  - Preserve v5 loading, truncation, NULL/timestamp handling, chunk loading,
    post-load orphan checks, row reporting, and rollback behavior.

Typical use (preflight + full validation + load):
  python 01_load_starter_sql_tsvs_v6.py \
    --tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs \
    --mysql-db gc3_dynamics \
    --mysql-user root

Validate everything without loading:
  python 01_load_starter_sql_tsvs_v6.py \
    --tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs \
    --mysql-db gc3_dynamics \
    --mysql-user root \
    --check-only

Fast preflight only:
  python 01_load_starter_sql_tsvs_v6.py \
    --tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs \
    --mysql-db gc3_dynamics \
    --mysql-user root \
    --preflight-only

Repeat a previously validated load without rescanning every TSV row:
  python 01_load_starter_sql_tsvs_v6.py \
    --tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs \
    --mysql-db gc3_dynamics \
    --mysql-user root \
    --skip-data-validation
"""

from __future__ import annotations

import argparse
import csv
import getpass
import gzip
import importlib.metadata
import os
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

try:
    import pymysql
except ImportError:  # handled in main
    pymysql = None


TABLE_COLUMNS: Dict[str, List[str]] = {
    "natural_history": [
        "species_pk",
        "species_normalized",
        "mass_grams_oskyrko",
        "mass_grams_title",
        "mass_grams_ji",
        "max_female_length_svl_mm_oskyrko",
        "svl_mm_title",
        "largest_clutch_size",
        "mean_number_of_eggs_per_clutch",
        "number_clutches_per_year",
        "annual_number_of_eggs",
        "max_longevity_years",
        "range_size",
        "genome_size",
        "ct_min",
        "ct_max",
    ],
    "species_name_audit": [
        "species_name_audit_pk", "species_pk", "source_dataset", "source_species_name",
        "species_normalized", "match_status",
    ],
    "genomes": [
        "genome_pk", "species_ncbi", "species_pk", "accession_id", "is_current",
    ],
    "sequences": [
        "sequence_pk", "genome_pk", "sequence_id", "sequence_length", "sequence_type", "gc",
    ],
    "genes": [
        "gene_pk", "sequence_pk", "gene_id", "g_strand", "g_seq_start", "g_seq_end",
        "length", "gc",
    ],
    "igs": [
        "igs_pk", "sequence_pk", "igs_seq_start", "igs_seq_end", "length", "gc",
    ],
    "gene_summary": [
        "gene_summary_pk", "genome_pk", "species_pk", "n_genes", "callable_bp_total",
        "mean_gc", "weighted_mean_gc", "sd_gc", "var_gc", "median_gc", "mad_gc",
        "iqr_gc", "q05_gc", "q25_gc", "q75_gc", "q95_gc", "mean_gene_length",
        "median_gene_length", "recovery_status", "recovery_notes", "created_at",
    ],
    "igs_summary": [
        "igs_summary_pk", "genome_pk", "species_pk", "n_igs", "callable_bp_total",
        "mean_gc", "weighted_mean_gc", "sd_gc", "var_gc", "median_gc", "mad_gc",
        "iqr_gc", "q05_gc", "q25_gc", "q75_gc", "q95_gc", "mean_igs_length",
        "median_igs_length", "recovery_status", "recovery_notes", "created_at",
    ],
    "analysis_run": [
        "run_pk", "analysis_name", "software", "software_version", "mask_mode",
        "gap_break_bp", "min_callable_frac", "created_at", "notes",
    ],
    "window_set": [
        "window_set_pk", "run_pk", "genome_pk", "tiling_type", "standard_window_size_bp",
        "step_size_bp", "start_offset_bp", "seq_scope", "notes",
    ],
    "genomic_windows": [
        "window_pk", "window_set_pk", "sequence_pk", "window_rank", "start_bp", "end_bp",
        "mid_bp", "standard_width_bp", "width_actual_bp", "callable_frac", "keep_flag",
    ],
    "gc_window_stats": [
        "window_pk", "a_count", "c_count", "g_count", "t_count", "n_count", "other_count",
        "callable_bp", "gc_bp", "gc_prop", "callable_frac", "masked_bp", "gap_bp",
    ],
    "sequence_summary": [
        "sequence_summary_pk", "run_pk", "sequence_pk", "genome_pk", "standard_window_size_bp",
        "step_size_bp", "tiling_type", "mask_mode", "n_windows_total", "n_windows_kept",
        "n_windows_excluded_missing", "n_windows_excluded_short", "mean_gc", "weighted_mean_gc",
        "sd_gc", "var_gc", "median_gc", "mad_gc", "iqr_gc", "q05_gc", "q25_gc",
        "q75_gc", "q95_gc", "mean_callable_fraction", "median_callable_fraction",
        "callable_bp_total", "gap_bp_total",
    ],
    "genome_summary": [
        "genome_summary_pk", "run_pk", "genome_pk", "species_pk", "standard_window_size_bp",
        "step_size_bp", "tiling_type", "mask_mode", "n_windows_total", "n_windows_kept",
        "mean_gc", "weighted_mean_gc", "sd_gc", "var_gc", "median_gc", "mad_gc", "iqr_gc",
        "q05_gc", "q25_gc", "q75_gc", "q95_gc", "mean_callable_fraction", "callable_bp_total",
        "gap_bp_total", "seq_count_used", "largest_seq_fraction",
    ],
}

# Load parents before children. analysis_run remains optional at the file level
# because a compatible run_pk may already exist in SQL, but non-NULL run_pk values
# are validated against analysis_run when the TSV is present.
LOAD_ORDER = [
    "natural_history",
    "genomes",
    "sequences",
    "genes",
    "igs",
    "gene_summary",
    "igs_summary",
    "species_name_audit",
    "analysis_run",
    "window_set",
    "genomic_windows",
    "gc_window_stats",
    "sequence_summary",
    "genome_summary",
]

# Reverse dependency order for safe truncation.
TRUNCATE_ORDER = list(reversed(LOAD_ORDER))

REQUIRED_FILES = [
    "natural_history",
    "genomes",
    "sequences",
    "genes",
    "igs",
    "gene_summary",
    "igs_summary",
    "species_name_audit",
    "window_set",
    "genomic_windows",
    "gc_window_stats",
    "sequence_summary",
    "genome_summary",
]
OPTIONAL_FILES = ["analysis_run"]

# 00b v15 retains the v12 chunk layout for these large/window-derived tables:
#   tsv_dir/genomic_windows/*.tsv.gz
#   tsv_dir/gc_window_stats/*.tsv.gz
#   tsv_dir/sequence_summary/*.tsv.gz
#   tsv_dir/genome_summary/*.tsv.gz
# If stale legacy files like genomic_windows.tsv are still present from v10/v11,
# they must not be read together with chunks because that produces duplicate PKs.
CHUNKABLE_TABLES = {"genomic_windows", "gc_window_stats", "sequence_summary", "genome_summary"}

# Files produced for auditing/bookkeeping rather than represented by tables in
# gc3_dynamics_v8_reduced. They are intentionally never validated or loaded.
IGNORED_NON_SCHEMA_FILES = {
    "window_chunk_manifest.tsv",
    "sequence_type_audit.tsv",
}

NULL_VALUES = {"", r"\N", "NULL", "null", "None", "none", "NaN", "nan", "NA", "na"}


def is_null(value: object) -> bool:
    return value is None or str(value).strip() in NULL_VALUES


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def python_environment_label() -> str:
    """Return a concise description of the active Python environment."""
    conda_env = os.environ.get("CONDA_DEFAULT_ENV")
    conda_prefix = os.environ.get("CONDA_PREFIX")
    virtual_env = os.environ.get("VIRTUAL_ENV")
    if conda_env:
        return f"conda:{conda_env}" + (f" ({conda_prefix})" if conda_prefix else "")
    if virtual_env:
        return f"virtualenv:{virtual_env}"
    return f"sys.prefix:{sys.prefix}"


def pymysql_version() -> str:
    """Return the installed PyMySQL distribution version when available."""
    if pymysql is None:
        return "unavailable"
    try:
        return importlib.metadata.version("PyMySQL")
    except importlib.metadata.PackageNotFoundError:
        return str(getattr(pymysql, "__version__", "unknown"))


def check_python_dependency() -> None:
    """Fail early if PyMySQL is unavailable and report the active interpreter."""
    print(f"[OK] Python executable: {sys.executable}")
    print(f"[OK] Python environment: {python_environment_label()}")
    if pymysql is None:
        raise ImportError(
            "PyMySQL is required but could not be imported by the active Python executable. "
            f"Interpreter: {sys.executable}. Install PyMySQL into this same environment "
            "(for example: python -m pip install pymysql, or install it with conda/mamba)."
        )
    print(f"[OK] Python dependency: pymysql {pymysql_version()}")


def open_text(path: Path):
    """Open plain .tsv or compressed .tsv.gz for text reading."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return path.open("r", newline="")


def table_input_files(tsv_dir: Path, table: str) -> List[Path]:
    """Return input files for one table.

    Supports both:
      1. legacy single-file layout: tsv_dir/table.tsv
      2. v12 chunked layout:      tsv_dir/table/*.tsv or *.tsv.gz

    For high-volume v12 chunkable tables, chunk directories take precedence.
    This avoids reading stale legacy combined files from an older run together
    with the new chunks, which would create duplicate primary keys.
    """
    legacy = tsv_dir / f"{table}.tsv"
    table_dir = tsv_dir / table

    chunk_files: List[Path] = []
    if table_dir.exists() and table_dir.is_dir():
        chunk_files.extend(sorted(table_dir.glob("*.tsv")))
        chunk_files.extend(sorted(table_dir.glob("*.tsv.gz")))

    if table in CHUNKABLE_TABLES and chunk_files:
        return chunk_files

    files: List[Path] = []
    if legacy.exists():
        files.append(legacy)
    files.extend(chunk_files)

    # Remove accidental duplicates while preserving order.
    seen = set()
    unique: List[Path] = []
    for path in files:
        resolved = path.resolve()
        if resolved not in seen:
            unique.append(path)
            seen.add(resolved)
    return unique


def read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"No header found in {path}")
        return [dict(row) for row in reader]


def iter_tsv_rows(path: Path):
    if not path.exists():
        raise FileNotFoundError(path)
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"No header found in {path}")
        for row in reader:
            yield dict(row)


def existing_tables(tsv_dir: Path) -> List[str]:
    """Return only tables represented in TABLE_COLUMNS/LOAD_ORDER.

    Non-schema helper files such as window_chunk_manifest.tsv and
    sequence_type_audit.tsv are intentionally ignored even when present.
    """
    tables = []
    for table in LOAD_ORDER:
        if table_input_files(tsv_dir, table):
            tables.append(table)
    return tables


def validate_required_inputs(tsv_dir: Path) -> List[str]:
    """Fast presence/readability checks without scanning TSV data rows."""
    if not tsv_dir.exists():
        raise FileNotFoundError(f"TSV directory does not exist: {tsv_dir}")
    if not tsv_dir.is_dir():
        raise NotADirectoryError(f"--tsv-dir is not a directory: {tsv_dir}")
    if not os.access(tsv_dir, os.R_OK | os.X_OK):
        raise PermissionError(f"TSV directory is not readable/searchable: {tsv_dir}")

    errors: List[str] = []

    # Non-chunkable required tables are expected as <table>.tsv.
    for table in REQUIRED_FILES:
        legacy = tsv_dir / f"{table}.tsv"
        table_dir = tsv_dir / table

        if table not in CHUNKABLE_TABLES:
            if not legacy.exists():
                errors.append(f"Missing required TSV: {legacy}")
            elif not legacy.is_file():
                errors.append(f"Required TSV path is not a regular file: {legacy}")
            continue

        # Chunkable tables may use either a legacy combined TSV or a chunk
        # directory. When the directory exists, it takes precedence and must
        # contain usable chunks; we do not silently fall back to a stale legacy file.
        if table_dir.exists():
            if not table_dir.is_dir():
                errors.append(f"Expected chunk path to be a directory: {table_dir}")
                continue
            chunk_files = sorted(table_dir.glob("*.tsv")) + sorted(table_dir.glob("*.tsv.gz"))
            if not chunk_files:
                errors.append(
                    f"Chunk directory exists but contains no .tsv or .tsv.gz files: {table_dir}"
                )
        elif not legacy.exists():
            errors.append(
                f"Missing required input for {table}: expected {legacy} or chunk directory {table_dir}"
            )
        elif not legacy.is_file():
            errors.append(f"Required TSV path is not a regular file: {legacy}")

    # analysis_run remains optional at the file level, but if supplied it should
    # use the normal single-TSV layout.
    optional_analysis = tsv_dir / "analysis_run.tsv"
    optional_analysis_dir = tsv_dir / "analysis_run"
    if optional_analysis_dir.exists():
        errors.append(
            f"analysis_run is not a chunked table; use {optional_analysis} instead of {optional_analysis_dir}"
        )
    if optional_analysis.exists() and not optional_analysis.is_file():
        errors.append(f"Optional TSV path is not a regular file: {optional_analysis}")

    tables = existing_tables(tsv_dir)
    tables = [table for table in LOAD_ORDER if table in tables]
    for table in tables:
        for path in table_input_files(tsv_dir, table):
            if not path.is_file():
                errors.append(f"Input is not a regular file: {path}")
            elif not os.access(path, os.R_OK):
                errors.append(f"Input file is not readable: {path}")

    if errors:
        raise FileNotFoundError("Fast input preflight failed:\n  " + "\n  ".join(errors))
    return tables


def report_ignored_inputs(tsv_dir: Path) -> None:
    """Report helper files and stale legacy combined TSVs that will be ignored."""
    ignored_files = [
        tsv_dir / filename
        for filename in sorted(IGNORED_NON_SCHEMA_FILES)
        if (tsv_dir / filename).exists()
    ]
    if ignored_files:
        print("[INFO] Ignoring non-schema helper file(s):", file=sys.stderr)
        for path in ignored_files:
            print(f"  {path}", file=sys.stderr)

    ignored_legacy = []
    for table in CHUNKABLE_TABLES:
        legacy = tsv_dir / f"{table}.tsv"
        table_dir = tsv_dir / table
        if legacy.exists() and table_dir.exists() and list(table_dir.glob("*.tsv*")):
            ignored_legacy.append(legacy)
    if ignored_legacy:
        print(
            "[WARN] Ignoring stale legacy combined TSV(s) because chunk directories are present:",
            file=sys.stderr,
        )
        for path in sorted(ignored_legacy):
            print(f"  {path}", file=sys.stderr)


def validate_headers(tsv_dir: Path, tables: Sequence[str]) -> None:
    errors: List[str] = []
    for table in tables:
        files = table_input_files(tsv_dir, table)
        if not files:
            if table in OPTIONAL_FILES:
                continue
            errors.append(f"Missing required file(s) for table: {table}")
            continue
        expected = TABLE_COLUMNS[table]
        for path in files:
            with open_text(path) as handle:
                reader = csv.reader(handle, delimiter="\t")
                try:
                    header = next(reader)
                except StopIteration:
                    errors.append(f"Empty TSV file: {path}")
                    continue
            if header != expected:
                errors.append(
                    f"Header mismatch for {path}\n"
                    f"  expected: {expected}\n"
                    f"  observed: {header}"
                )
    if errors:
        raise ValueError("\n".join(errors))


def read_table_files(tsv_dir: Path, table: str) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for path in table_input_files(tsv_dir, table):
        rows.extend(read_tsv(path))
    return rows

def ids(rows: Iterable[Dict[str, str]], column: str) -> Set[str]:
    return {clean(row.get(column)) for row in rows if not is_null(row.get(column))}


def find_missing_refs(
    child_rows: Sequence[Dict[str, str]],
    child_table: str,
    child_col: str,
    parent_ids: Set[str],
    parent_table: str,
    allow_null: bool = False,
    max_examples: int = 10,
) -> List[str]:
    missing: List[str] = []
    for row in child_rows:
        val = clean(row.get(child_col))
        if is_null(val):
            if not allow_null:
                missing.append(r"\N")
            continue
        if val not in parent_ids:
            missing.append(val)
    if not missing:
        return []
    unique = []
    seen = set()
    for val in missing:
        if val not in seen:
            unique.append(val)
            seen.add(val)
    examples = ", ".join(unique[:max_examples])
    suffix = " ..." if len(unique) > max_examples else ""
    return [
        f"{child_table}.{child_col} has {len(missing)} row(s) not found in {parent_table}; examples: {examples}{suffix}"
    ]


def validate_relationships(tsv_dir: Path, tables: Sequence[str]) -> Dict[str, int]:
    rows: Dict[str, List[Dict[str, str]]] = {}
    row_counts: Dict[str, int] = {}

    for table in tables:
        files = table_input_files(tsv_dir, table)
        if files:
            table_rows = read_table_files(tsv_dir, table)
            rows[table] = table_rows
            row_counts[table] = len(table_rows)

    errors: List[str] = []

    species_pks = ids(rows.get("natural_history", []), "species_pk")
    genome_pks = ids(rows.get("genomes", []), "genome_pk")
    sequence_pks = ids(rows.get("sequences", []), "sequence_pk")
    run_pks = ids(rows.get("analysis_run", []), "run_pk")
    window_set_pks = ids(rows.get("window_set", []), "window_set_pk")
    window_pks = ids(rows.get("genomic_windows", []), "window_pk")

    if "genomes" in rows and species_pks:
        errors += find_missing_refs(rows["genomes"], "genomes", "species_pk", species_pks, "natural_history", allow_null=True)
    if "species_name_audit" in rows and species_pks:
        errors += find_missing_refs(rows["species_name_audit"], "species_name_audit", "species_pk", species_pks, "natural_history", allow_null=True)
    if "sequences" in rows:
        errors += find_missing_refs(rows["sequences"], "sequences", "genome_pk", genome_pks, "genomes")
    if "genes" in rows:
        errors += find_missing_refs(rows["genes"], "genes", "sequence_pk", sequence_pks, "sequences")
    if "igs" in rows:
        errors += find_missing_refs(rows["igs"], "igs", "sequence_pk", sequence_pks, "sequences")
    if "gene_summary" in rows:
        errors += find_missing_refs(rows["gene_summary"], "gene_summary", "genome_pk", genome_pks, "genomes")
        if species_pks:
            errors += find_missing_refs(rows["gene_summary"], "gene_summary", "species_pk", species_pks, "natural_history", allow_null=True)
    if "igs_summary" in rows:
        errors += find_missing_refs(rows["igs_summary"], "igs_summary", "genome_pk", genome_pks, "genomes")
        if species_pks:
            errors += find_missing_refs(rows["igs_summary"], "igs_summary", "species_pk", species_pks, "natural_history", allow_null=True)

    if "window_set" in rows:
        errors += find_missing_refs(rows["window_set"], "window_set", "genome_pk", genome_pks, "genomes")
        if run_pks:
            errors += find_missing_refs(rows["window_set"], "window_set", "run_pk", run_pks, "analysis_run")
        else:
            nonnull_run_refs = [r.get("run_pk") for r in rows["window_set"] if not is_null(r.get("run_pk"))]
            if nonnull_run_refs:
                errors.append("window_set.tsv contains run_pk values, but analysis_run.tsv was not found. Load/create analysis_run first.")

    if "genomic_windows" in rows:
        errors += find_missing_refs(rows["genomic_windows"], "genomic_windows", "window_set_pk", window_set_pks, "window_set")
        errors += find_missing_refs(rows["genomic_windows"], "genomic_windows", "sequence_pk", sequence_pks, "sequences")

    if "gc_window_stats" in rows:
        errors += find_missing_refs(rows["gc_window_stats"], "gc_window_stats", "window_pk", window_pks, "genomic_windows")
        stats_pks = ids(rows["gc_window_stats"], "window_pk")
        missing_stats = window_pks - stats_pks
        if missing_stats:
            examples = ", ".join(sorted(missing_stats)[:10])
            errors.append(f"{len(missing_stats)} genomic_windows row(s) have no gc_window_stats row; examples: {examples}")

    if "sequence_summary" in rows:
        errors += find_missing_refs(rows["sequence_summary"], "sequence_summary", "sequence_pk", sequence_pks, "sequences")
        errors += find_missing_refs(rows["sequence_summary"], "sequence_summary", "genome_pk", genome_pks, "genomes")
        if run_pks:
            errors += find_missing_refs(rows["sequence_summary"], "sequence_summary", "run_pk", run_pks, "analysis_run")

    if "genome_summary" in rows:
        errors += find_missing_refs(rows["genome_summary"], "genome_summary", "genome_pk", genome_pks, "genomes")
        if species_pks:
            errors += find_missing_refs(rows["genome_summary"], "genome_summary", "species_pk", species_pks, "natural_history", allow_null=True)
        if run_pks:
            errors += find_missing_refs(rows["genome_summary"], "genome_summary", "run_pk", run_pks, "analysis_run")

    # Primary-key duplicate checks across all files/chunks for each table.
    pk_cols = {
        "natural_history": "species_pk",
        "species_name_audit": "species_name_audit_pk",
        "genomes": "genome_pk",
        "sequences": "sequence_pk",
        "genes": "gene_pk",
        "igs": "igs_pk",
        "gene_summary": "gene_summary_pk",
        "igs_summary": "igs_summary_pk",
        "analysis_run": "run_pk",
        "window_set": "window_set_pk",
        "genomic_windows": "window_pk",
        "gc_window_stats": "window_pk",
        "sequence_summary": "sequence_summary_pk",
        "genome_summary": "genome_summary_pk",
    }
    for table, pk_col in pk_cols.items():
        if table not in rows:
            continue
        vals = [clean(r.get(pk_col)) for r in rows[table] if not is_null(r.get(pk_col))]
        if len(vals) != len(set(vals)):
            errors.append(f"Duplicate primary-key values detected in {table}.{pk_col}")

    # Schema-defined unique keys for the 00c outputs.
    unique_keys = {
        "natural_history": ("species_normalized",),
        "genomes": ("accession_id",),
        "sequences": ("genome_pk", "sequence_id"),
        "genes": ("sequence_pk", "gene_id"),
        "igs": ("sequence_pk", "igs_seq_start", "igs_seq_end"),
        "gene_summary": ("genome_pk",),
        "igs_summary": ("genome_pk",),
    }
    for table, columns in unique_keys.items():
        if table not in rows:
            continue
        seen_keys = set()
        duplicate_examples = []
        for row in rows[table]:
            key = tuple(clean(row.get(column)) for column in columns)
            if any(is_null(value) for value in key):
                continue
            if key in seen_keys and len(duplicate_examples) < 10:
                duplicate_examples.append(key)
            seen_keys.add(key)
        if duplicate_examples:
            formatted = ", ".join(str(key) for key in duplicate_examples)
            errors.append(
                f"Duplicate schema unique-key values detected in {table}{columns}; examples: {formatted}"
            )

    if errors:
        raise ValueError("Relationship validation failed:\n  " + "\n  ".join(errors))

    return row_counts


def count_rows_in_file(path: Path) -> int:
    # Count data rows only, excluding header.
    with open_text(path) as handle:
        return max(0, sum(1 for _ in handle) - 1)

def resolve_mysql_password(args: argparse.Namespace) -> str:
    """Prompt at most once, then retain the password in this process for reconnects."""
    if args.mysql_password is None:
        args.mysql_password = getpass.getpass(
            f"MySQL password for {args.mysql_user}@{args.mysql_host}: "
        )
    return args.mysql_password


def connect(args: argparse.Namespace):
    if pymysql is None:
        raise ImportError(
            "PyMySQL is required but is not importable by " + sys.executable
        )
    password = resolve_mysql_password(args)
    return pymysql.connect(
        host=args.mysql_host,
        port=args.mysql_port,
        user=args.mysql_user,
        password=password,
        db=args.mysql_db,
        charset="utf8mb4",
        cursorclass=pymysql.cursors.DictCursor,
        local_infile=True,
        autocommit=False,
    )


def mysql_server_tables(cursor, database: str) -> Set[str]:
    cursor.execute(
        "SELECT TABLE_NAME FROM INFORMATION_SCHEMA.TABLES "
        "WHERE TABLE_SCHEMA=%s AND TABLE_TYPE='BASE TABLE'",
        (database,),
    )
    return {str(row["TABLE_NAME"]) for row in cursor.fetchall()}


def mysql_table_columns(cursor, database: str, table: str) -> List[str]:
    cursor.execute(
        "SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS "
        "WHERE TABLE_SCHEMA=%s AND TABLE_NAME=%s ORDER BY ORDINAL_POSITION",
        (database, table),
    )
    return [str(row["COLUMN_NAME"]) for row in cursor.fetchall()]


def validate_mysql_schema(cursor, database: str) -> None:
    """Verify required SQL tables and loader-required columns exist."""
    required_tables = set(LOAD_ORDER)
    present_tables = mysql_server_tables(cursor, database)
    missing_tables = sorted(required_tables - present_tables)
    if missing_tables:
        raise RuntimeError(
            "Target database is missing required table(s): " + ", ".join(missing_tables)
        )

    column_errors: List[str] = []
    extra_columns: Dict[str, List[str]] = {}
    for table in LOAD_ORDER:
        observed = mysql_table_columns(cursor, database, table)
        expected = TABLE_COLUMNS[table]
        missing = [column for column in expected if column not in observed]
        if missing:
            column_errors.append(f"{table}: missing column(s): {', '.join(missing)}")
        extras = [column for column in observed if column not in expected]
        if extras:
            extra_columns[table] = extras

    if column_errors:
        raise RuntimeError(
            "MySQL table-column preflight failed:\n  " + "\n  ".join(column_errors)
        )

    if extra_columns:
        print(
            "[WARN] MySQL contains additional column(s) not used by this loader:",
            file=sys.stderr,
        )
        for table, columns in extra_columns.items():
            print(f"  {table}: {', '.join(columns)}", file=sys.stderr)


def validate_local_infile(conn, cursor) -> None:
    """Verify LOCAL INFILE is enabled in the client connection and MySQL server."""
    client_enabled: Optional[bool] = None
    try:
        local_files_flag = pymysql.constants.CLIENT.LOCAL_FILES
        client_enabled = bool(getattr(conn, "client_flag", 0) & local_files_flag)
    except (AttributeError, TypeError):
        # Compatibility fallback for alternate/older PyMySQL connection objects.
        private_setting = getattr(conn, "_local_infile", None)
        if private_setting is not None:
            client_enabled = bool(private_setting)

    if client_enabled is not True:
        raise RuntimeError(
            "PyMySQL client connection does not report the LOCAL_FILES capability as enabled. "
            "The loader connects with local_infile=True, so this indicates an incompatible or "
            "unexpected client configuration."
        )

    cursor.execute("SHOW VARIABLES LIKE 'local_infile';")
    row = cursor.fetchone()
    if not row:
        raise RuntimeError("Could not read the MySQL local_infile server variable.")

    # DictCursor key names are normally Variable_name/Value, but use a value scan
    # as a fallback to tolerate case differences across server versions.
    server_value = row.get("Value")
    if server_value is None:
        values = [str(v) for k, v in row.items() if str(k).lower() != "variable_name"]
        server_value = values[0] if values else None
    if str(server_value).strip().lower() not in {"on", "1", "true", "yes"}:
        raise RuntimeError(
            f"MySQL server local_infile is not enabled (reported value: {server_value!r})."
        )


def mysql_preflight(args: argparse.Namespace) -> None:
    """Connect and verify database/schema/LOCAL INFILE without modifying tables."""
    conn = connect(args)
    try:
        with conn.cursor() as cursor:
            cursor.execute("SELECT DATABASE() AS db;")
            selected_db = cursor.fetchone()["db"]
            if selected_db != args.mysql_db:
                raise RuntimeError(
                    f"Connected to unexpected database {selected_db!r}; expected {args.mysql_db!r}."
                )
            print(
                f"[OK] MySQL connection successful: "
                f"{args.mysql_user}@{args.mysql_host}:{args.mysql_port}/{args.mysql_db}"
            )

            validate_mysql_schema(cursor, args.mysql_db)
            print("[OK] Required MySQL tables found")
            print("[OK] MySQL table columns match loader requirements")

            validate_local_infile(conn, cursor)
            print("[OK] LOCAL INFILE enabled for client/server")
    finally:
        conn.close()


def run_fast_preflight(args: argparse.Namespace, tsv_dir: Path) -> List[str]:
    """Run all inexpensive checks before any full TSV row scan."""
    print("[INFO] Running fast preflight checks...")
    check_python_dependency()

    tables = validate_required_inputs(tsv_dir)
    report_ignored_inputs(tsv_dir)
    print("[OK] Required table inputs found and readable")

    validate_headers(tsv_dir, tables)
    print("[OK] TSV headers match expected schema")

    mysql_preflight(args)
    print("[OK] Fast preflight passed")
    return tables


def load_plain_tsv(cursor, table: str, path: Path) -> int:
    """Load a plain TSV while converting blank/\\N fields to SQL NULL.

    00c intentionally leaves unavailable values blank in partial-recovery rows.
    Loading through user variables lets the loader convert those blanks to NULL
    instead of asking MySQL to coerce empty strings into numeric/datetime columns.
    created_at accepts either a standard MySQL datetime or the ISO-8601 UTC
    representation written by 00c (YYYY-MM-DDTHH:MM:SSZ).
    """
    columns = TABLE_COLUMNS[table]
    variables = [f"@v_{i}" for i in range(len(columns))]
    input_sql = ", ".join(variables)

    assignments = []
    for variable, column in zip(variables, columns):
        cleaned = f"NULLIF(NULLIF({variable}, ''), CONCAT(CHAR(92), 'N'))"
        if column == "created_at":
            assignment = (
                f"`{column}` = CASE "
                f"WHEN {cleaned} IS NULL THEN NULL "
                f"WHEN {cleaned} LIKE '%T%Z' THEN "
                f"STR_TO_DATE(REPLACE(REPLACE({cleaned}, 'T', ' '), 'Z', ''), '%Y-%m-%d %H:%i:%s') "
                f"ELSE {cleaned} END"
            )
        else:
            assignment = f"`{column}` = {cleaned}"
        assignments.append(assignment)

    set_sql = ", ".join(assignments)
    sql = (
        f"LOAD DATA LOCAL INFILE %s INTO TABLE `{table}` "
        "FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' "
        f"IGNORE 1 LINES ({input_sql}) SET {set_sql};"
    )
    cursor.execute(sql, (str(path.resolve()),))
    return cursor.rowcount


def load_table_file(cursor, table: str, path: Path) -> int:
    """Load one plain TSV or gzip-compressed TSV chunk into MySQL.

    MySQL LOAD DATA LOCAL INFILE cannot read gzip directly, so compressed
    chunks are decompressed to a temporary TSV and loaded one chunk at a time.
    """
    if str(path).endswith(".gz"):
        with tempfile.NamedTemporaryFile("wb", suffix=f".{table}.tsv", delete=False) as tmp:
            tmp_path = Path(tmp.name)
            with gzip.open(path, "rb") as source:
                shutil.copyfileobj(source, tmp)
        try:
            return load_plain_tsv(cursor, table, tmp_path)
        finally:
            tmp_path.unlink(missing_ok=True)
    return load_plain_tsv(cursor, table, path)

def table_count(cursor, table: str) -> int:
    cursor.execute(f"SELECT COUNT(*) AS n FROM `{table}`;")
    return int(cursor.fetchone()["n"])


def database_fk_checks(cursor) -> List[str]:
    """Run lightweight orphan checks after load."""
    checks = [
        ("genomes.species_pk", "SELECT COUNT(*) AS n FROM genomes g LEFT JOIN natural_history nh ON g.species_pk = nh.species_pk WHERE g.species_pk IS NOT NULL AND nh.species_pk IS NULL"),
        ("species_name_audit.species_pk", "SELECT COUNT(*) AS n FROM species_name_audit sna LEFT JOIN natural_history nh ON sna.species_pk = nh.species_pk WHERE sna.species_pk IS NOT NULL AND nh.species_pk IS NULL"),
        ("sequences.genome_pk", "SELECT COUNT(*) AS n FROM sequences s LEFT JOIN genomes g ON s.genome_pk = g.genome_pk WHERE g.genome_pk IS NULL"),
        ("genes.sequence_pk", "SELECT COUNT(*) AS n FROM genes ge LEFT JOIN sequences s ON ge.sequence_pk = s.sequence_pk WHERE s.sequence_pk IS NULL"),
        ("igs.sequence_pk", "SELECT COUNT(*) AS n FROM igs i LEFT JOIN sequences s ON i.sequence_pk = s.sequence_pk WHERE s.sequence_pk IS NULL"),
        ("gene_summary.genome_pk", "SELECT COUNT(*) AS n FROM gene_summary gs LEFT JOIN genomes g ON gs.genome_pk = g.genome_pk WHERE g.genome_pk IS NULL"),
        ("gene_summary.species_pk", "SELECT COUNT(*) AS n FROM gene_summary gs LEFT JOIN natural_history nh ON gs.species_pk = nh.species_pk WHERE gs.species_pk IS NOT NULL AND nh.species_pk IS NULL"),
        ("igs_summary.genome_pk", "SELECT COUNT(*) AS n FROM igs_summary isy LEFT JOIN genomes g ON isy.genome_pk = g.genome_pk WHERE g.genome_pk IS NULL"),
        ("igs_summary.species_pk", "SELECT COUNT(*) AS n FROM igs_summary isy LEFT JOIN natural_history nh ON isy.species_pk = nh.species_pk WHERE isy.species_pk IS NOT NULL AND nh.species_pk IS NULL"),
        ("window_set.run_pk", "SELECT COUNT(*) AS n FROM window_set ws LEFT JOIN analysis_run ar ON ws.run_pk = ar.run_pk WHERE ws.run_pk IS NOT NULL AND ar.run_pk IS NULL"),
        ("window_set.genome_pk", "SELECT COUNT(*) AS n FROM window_set ws LEFT JOIN genomes g ON ws.genome_pk = g.genome_pk WHERE g.genome_pk IS NULL"),
        ("genomic_windows.window_set_pk", "SELECT COUNT(*) AS n FROM genomic_windows gw LEFT JOIN window_set ws ON gw.window_set_pk = ws.window_set_pk WHERE ws.window_set_pk IS NULL"),
        ("genomic_windows.sequence_pk", "SELECT COUNT(*) AS n FROM genomic_windows gw LEFT JOIN sequences s ON gw.sequence_pk = s.sequence_pk WHERE s.sequence_pk IS NULL"),
        ("gc_window_stats.window_pk", "SELECT COUNT(*) AS n FROM gc_window_stats gcs LEFT JOIN genomic_windows gw ON gcs.window_pk = gw.window_pk WHERE gw.window_pk IS NULL"),
        ("sequence_summary.run_pk", "SELECT COUNT(*) AS n FROM sequence_summary ss LEFT JOIN analysis_run ar ON ss.run_pk = ar.run_pk WHERE ss.run_pk IS NOT NULL AND ar.run_pk IS NULL"),
        ("sequence_summary.sequence_pk", "SELECT COUNT(*) AS n FROM sequence_summary ss LEFT JOIN sequences s ON ss.sequence_pk = s.sequence_pk WHERE s.sequence_pk IS NULL"),
        ("sequence_summary.genome_pk", "SELECT COUNT(*) AS n FROM sequence_summary ss LEFT JOIN genomes g ON ss.genome_pk = g.genome_pk WHERE g.genome_pk IS NULL"),
        ("genome_summary.run_pk", "SELECT COUNT(*) AS n FROM genome_summary gs LEFT JOIN analysis_run ar ON gs.run_pk = ar.run_pk WHERE gs.run_pk IS NOT NULL AND ar.run_pk IS NULL"),
        ("genome_summary.genome_pk", "SELECT COUNT(*) AS n FROM genome_summary gs LEFT JOIN genomes g ON gs.genome_pk = g.genome_pk WHERE g.genome_pk IS NULL"),
        ("genome_summary.species_pk", "SELECT COUNT(*) AS n FROM genome_summary gs LEFT JOIN natural_history nh ON gs.species_pk = nh.species_pk WHERE gs.species_pk IS NOT NULL AND nh.species_pk IS NULL"),
    ]
    problems = []
    for label, sql in checks:
        cursor.execute(sql)
        n = int(cursor.fetchone()["n"])
        if n:
            problems.append(f"{label}: {n} orphan row(s)")
    return problems


def truncate_tables(cursor, tables: Sequence[str]) -> None:
    cursor.execute("SET FOREIGN_KEY_CHECKS=0;")
    try:
        for table in TRUNCATE_ORDER:
            if table in tables:
                cursor.execute(f"TRUNCATE TABLE `{table}`;")
    finally:
        cursor.execute("SET FOREIGN_KEY_CHECKS=1;")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Load 00b/00c SQL TSVs into MySQL using gc3_dynamics_v8_reduced. "
            "Runs a fast environment/file/MySQL preflight before optional full-row validation."
        )
    )
    parser.add_argument(
        "--tsv-dir",
        required=True,
        help="Directory containing TSVs from 00b_build_starter_sql_tsvs_v15.py and 00c_build_tsvs_from_gff.py",
    )
    parser.add_argument("--mysql-db", default="gc3_dynamics", help="MySQL database name. Default: gc3_dynamics")
    parser.add_argument("--mysql-user", default="root", help="MySQL user. Default: root")
    parser.add_argument("--mysql-password", default=None, help="MySQL password. Omit to be prompted securely.")
    parser.add_argument("--mysql-host", default="localhost", help="MySQL host. Default: localhost")
    parser.add_argument("--mysql-port", type=int, default=3306, help="MySQL port. Default: 3306")
    parser.add_argument(
        "--check-only",
        action="store_true",
        help=(
            "Run fast preflight plus complete TSV relationship/duplicate validation, "
            "then exit without changing MySQL."
        ),
    )
    parser.add_argument(
        "--preflight-only",
        action="store_true",
        help=(
            "Run only fast dependency/file/header/MySQL schema/LOCAL INFILE checks; "
            "do not scan all TSV rows or change MySQL."
        ),
    )
    parser.add_argument(
        "--skip-data-validation",
        action="store_true",
        help=(
            "Skip expensive full-row TSV relationship and duplicate-key validation after "
            "the mandatory fast preflight. Intended only for a dataset that already passed full validation."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Run preflight and the selected validation mode, print load order, but do not load.",
    )
    parser.add_argument(
        "--truncate",
        action="store_true",
        help="Truncate target tables before loading, in safe reverse dependency order.",
    )
    parser.add_argument(
        "--no-post-check",
        action="store_true",
        help="Skip post-load database orphan checks. This remains independent of --skip-data-validation.",
    )
    args = parser.parse_args()

    if args.check_only and args.skip_data_validation:
        parser.error("--check-only requires full data validation and cannot be combined with --skip-data-validation")
    if args.preflight_only and args.check_only:
        parser.error("--preflight-only and --check-only are mutually exclusive")
    if args.preflight_only and args.dry_run:
        parser.error("--preflight-only already exits without loading; do not combine it with --dry-run")
    if args.preflight_only and args.skip_data_validation:
        parser.error("--preflight-only does not run full data validation; --skip-data-validation is unnecessary")
    if args.preflight_only and args.truncate:
        parser.error("--preflight-only cannot be combined with --truncate")
    if args.check_only and args.truncate:
        parser.error("--check-only cannot be combined with --truncate")
    if args.dry_run and args.truncate:
        parser.error("--dry-run cannot be combined with --truncate because no SQL changes are made")

    return args


def print_load_order(tsv_dir: Path, tables: Sequence[str]) -> None:
    print("\n[INFO] Load order:")
    for table in tables:
        files = table_input_files(tsv_dir, table)
        print(f"  {table} ({len(files)} file(s))")


def main() -> None:
    args = parse_args()
    tsv_dir = Path(args.tsv_dir).expanduser().resolve()

    # Mandatory fast checks happen before any expensive full-row TSV scan.
    tables = run_fast_preflight(args, tsv_dir)

    if "analysis_run" not in tables:
        print(
            "[WARN] analysis_run.tsv was not found. This is only safe if referenced run_pk "
            "values are nullable or already exist in SQL.",
            file=sys.stderr,
        )

    if args.preflight_only:
        print("\n[OK] Preflight-only checks passed. No full TSV scan or SQL changes made.")
        return

    row_counts: Optional[Dict[str, int]] = None
    if args.skip_data_validation and args.no_post_check:
        print(
            "[WARN] Both full pre-load data validation and post-load orphan checks were explicitly disabled.",
            file=sys.stderr,
        )
    if args.skip_data_validation:
        print("[INFO] Full data validation skipped (--skip-data-validation).")
    else:
        print("[INFO] Running full TSV relationship/duplicate validation...")
        row_counts = validate_relationships(tsv_dir, tables)
        print("[OK] Full TSV data validation passed.")
        for table in tables:
            files = table_input_files(tsv_dir, table)
            file_label = f"{len(files)} file(s)" if len(files) != 1 else files[0].name
            print(f"  {table}\t{row_counts.get(table, 0)} rows\t{file_label}")

    if args.check_only or args.dry_run:
        print_load_order(tsv_dir, tables)
        print("\n[OK] No SQL changes made.")
        return

    # Reconnect after validation. The full scan can be long enough for an idle
    # preflight connection to time out, so the preflight connection is intentionally
    # closed and a fresh connection is used for loading.
    conn = connect(args)
    try:
        with conn.cursor() as cursor:
            # Recheck LOCAL INFILE on the actual load connection. This is inexpensive
            # and protects against a changed server/session configuration.
            validate_local_infile(conn, cursor)

            if args.truncate:
                print("[INFO] Truncating tables in reverse dependency order...")
                truncate_tables(cursor, tables)

            print("[INFO] Loading TSVs...")
            for table in tables:
                table_total = 0
                files = table_input_files(tsv_dir, table)
                for path in files:
                    loaded = load_table_file(cursor, table, path)
                    table_total += loaded
                    print(f"  loaded {path.relative_to(tsv_dir)} into {table}: {loaded} row(s)")
                print(f"  [TABLE TOTAL] {table}: {table_total} row(s)")

            if not args.no_post_check:
                print("[INFO] Running post-load orphan checks...")
                problems = database_fk_checks(cursor)
                if problems:
                    raise RuntimeError(
                        "Post-load relationship checks failed:\n  " + "\n  ".join(problems)
                    )
                print("[OK] Post-load orphan checks passed.")
            else:
                print("[WARN] Post-load orphan checks skipped (--no-post-check).", file=sys.stderr)

            print("[INFO] Database row counts after load:")
            for table in tables:
                print(f"  {table}\t{table_count(cursor, table)} rows")

        conn.commit()
        print("\n[OK] Load committed successfully.")
    except Exception:
        conn.rollback()
        print("\n[ERROR] Load failed. Transaction rolled back where supported by MySQL.", file=sys.stderr)
        raise
    finally:
        conn.close()


if __name__ == "__main__":
    main()
