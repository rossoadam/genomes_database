#!/usr/bin/env python3
"""
Load TSV outputs from 00_build_starter_sql_tsvs_v12.py into a MySQL database.

This is intentionally standalone for troubleshooting. It validates TSV relationships
before loading, supports both legacy single-table TSVs and v12 table-specific
chunk directories containing .tsv.gz files, loads tables in dependency order,
and can optionally truncate the starter/window tables before loading. If v12 chunk directories are present for high-volume tables, stale legacy combined TSVs for those same tables are ignored by default.

Typical use:
  python 01_load_starter_sql_tsvs.py \
    --tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs \
    --mysql-db gc3_dynamics \
    --mysql-user root

Safer first pass:
  python 01_load_starter_sql_tsvs.py \
    --tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs \
    --mysql-db gc3_dynamics \
    --mysql-user root \
    --check-only

Reload from scratch:
  python 01_load_starter_sql_tsvs.py \
    --tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs \
    --mysql-db gc3_dynamics \
    --mysql-user root \
    --truncate
"""

from __future__ import annotations

import argparse
import csv
import getpass
import gzip
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
        "species_pk", "species_normalized", "mass_meiri", "mass_title", "mass_ji",
        "genome_size", "ct_min", "ct_max",
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

# Load parents before children. analysis_run is optional because the user did not
# list it, but v10 writes it and window_set/summary rows refer to run_pk.
LOAD_ORDER = [
    "natural_history",
    "genomes",
    "sequences",
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
    "species_name_audit",
    "window_set",
    "genomic_windows",
    "gc_window_stats",
    "sequence_summary",
    "genome_summary",
]
OPTIONAL_FILES = ["analysis_run"]

# v12 writes these large/window-derived tables as chunk directories:
#   tsv_dir/genomic_windows/*.tsv.gz
#   tsv_dir/gc_window_stats/*.tsv.gz
#   tsv_dir/sequence_summary/*.tsv.gz
#   tsv_dir/genome_summary/*.tsv.gz
# If stale legacy files like genomic_windows.tsv are still present from v10/v11,
# they must not be read together with chunks because that produces duplicate PKs.
CHUNKABLE_TABLES = {"genomic_windows", "gc_window_stats", "sequence_summary", "genome_summary"}

NULL_VALUES = {"", r"\N", "NULL", "null", "None", "none", "NaN", "nan", "NA", "na"}


def is_null(value: object) -> bool:
    return value is None or str(value).strip() in NULL_VALUES


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


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


def existing_tables(tsv_dir: Path, include_audit: bool = False) -> List[str]:
    tables = []
    for table in LOAD_ORDER:
        if table_input_files(tsv_dir, table):
            tables.append(table)
    if include_audit and (tsv_dir / "sequence_type_audit.tsv").exists():
        tables.append("sequence_type_audit")
    return tables


def validate_headers(tsv_dir: Path, tables: Sequence[str]) -> None:
    errors: List[str] = []
    for table in tables:
        if table == "sequence_type_audit":
            continue
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
        if table == "sequence_type_audit":
            continue
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

    if errors:
        raise ValueError("Relationship validation failed:\n  " + "\n  ".join(errors))

    return row_counts


def count_rows_in_file(path: Path) -> int:
    # Count data rows only, excluding header.
    with open_text(path) as handle:
        return max(0, sum(1 for _ in handle) - 1)

def connect(args: argparse.Namespace):
    if pymysql is None:
        raise ImportError("pymysql is required. Install with: pip install pymysql")
    password = args.mysql_password
    if password is None:
        password = getpass.getpass(f"MySQL password for {args.mysql_user}@{args.mysql_host}: ")
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


def load_plain_tsv(cursor, table: str, path: Path) -> int:
    columns = TABLE_COLUMNS[table]
    column_sql = ", ".join(f"`{col}`" for col in columns)
    sql = (
        f"LOAD DATA LOCAL INFILE %s INTO TABLE `{table}` "
        "FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n' "
        f"IGNORE 1 LINES ({column_sql});"
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
        ("sequences.genome_pk", "SELECT COUNT(*) AS n FROM sequences s LEFT JOIN genomes g ON s.genome_pk = g.genome_pk WHERE g.genome_pk IS NULL"),
        ("window_set.genome_pk", "SELECT COUNT(*) AS n FROM window_set ws LEFT JOIN genomes g ON ws.genome_pk = g.genome_pk WHERE g.genome_pk IS NULL"),
        ("genomic_windows.window_set_pk", "SELECT COUNT(*) AS n FROM genomic_windows gw LEFT JOIN window_set ws ON gw.window_set_pk = ws.window_set_pk WHERE ws.window_set_pk IS NULL"),
        ("genomic_windows.sequence_pk", "SELECT COUNT(*) AS n FROM genomic_windows gw LEFT JOIN sequences s ON gw.sequence_pk = s.sequence_pk WHERE s.sequence_pk IS NULL"),
        ("gc_window_stats.window_pk", "SELECT COUNT(*) AS n FROM gc_window_stats gcs LEFT JOIN genomic_windows gw ON gcs.window_pk = gw.window_pk WHERE gw.window_pk IS NULL"),
        ("sequence_summary.sequence_pk", "SELECT COUNT(*) AS n FROM sequence_summary ss LEFT JOIN sequences s ON ss.sequence_pk = s.sequence_pk WHERE s.sequence_pk IS NULL"),
        ("sequence_summary.genome_pk", "SELECT COUNT(*) AS n FROM sequence_summary ss LEFT JOIN genomes g ON ss.genome_pk = g.genome_pk WHERE g.genome_pk IS NULL"),
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
    parser = argparse.ArgumentParser(description="Load starter/window TSVs into MySQL with relationship checks.")
    parser.add_argument("--tsv-dir", required=True, help="Directory containing TSVs from 00_build_starter_sql_tsvs_v12.py")
    parser.add_argument("--mysql-db", default="gc3_dynamics", help="MySQL database name. Default: gc3_dynamics")
    parser.add_argument("--mysql-user", default="root", help="MySQL user. Default: root")
    parser.add_argument("--mysql-password", default=None, help="MySQL password. Omit to be prompted securely.")
    parser.add_argument("--mysql-host", default="localhost", help="MySQL host. Default: localhost")
    parser.add_argument("--mysql-port", type=int, default=3306, help="MySQL port. Default: 3306")
    parser.add_argument("--check-only", action="store_true", help="Validate TSVs and relationships, but do not connect/load.")
    parser.add_argument("--dry-run", action="store_true", help="Validate TSVs and print load order, but do not load.")
    parser.add_argument("--truncate", action="store_true", help="Truncate target tables before loading, in safe reverse dependency order.")
    parser.add_argument("--no-post-check", action="store_true", help="Skip post-load database orphan checks.")
    parser.add_argument("--include-sequence-type-audit", action="store_true", help="Reserved for later; not loaded by default because it is an audit table not listed in the request.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    tsv_dir = Path(args.tsv_dir).expanduser().resolve()
    if not tsv_dir.exists():
        raise FileNotFoundError(f"TSV directory does not exist: {tsv_dir}")

    tables = existing_tables(tsv_dir, include_audit=False)
    missing_required = [table for table in REQUIRED_FILES if table not in tables]
    if missing_required:
        raise FileNotFoundError("Missing required TSV files: " + ", ".join(f"{table}.tsv" for table in missing_required))

    # Keep only recognized/loadable tables and preserve dependency order.
    tables = [table for table in LOAD_ORDER if table in tables]

    ignored_legacy = []
    for table in CHUNKABLE_TABLES:
        legacy = tsv_dir / f"{table}.tsv"
        table_dir = tsv_dir / table
        if legacy.exists() and table_dir.exists() and list(table_dir.glob("*.tsv*")):
            ignored_legacy.append(legacy)
    if ignored_legacy:
        print("[WARN] Ignoring stale legacy combined TSV(s) because v12 chunk directories are present:", file=sys.stderr)
        for path in sorted(ignored_legacy):
            print(f"  {path}", file=sys.stderr)

    print("[INFO] Validating TSV headers...")
    validate_headers(tsv_dir, tables)
    print("[INFO] Validating TSV relationships...")
    row_counts = validate_relationships(tsv_dir, tables)

    print("[OK] TSV validation passed.")
    for table in tables:
        files = table_input_files(tsv_dir, table)
        file_label = f"{len(files)} file(s)" if len(files) != 1 else files[0].name
        print(f"  {table}\t{row_counts.get(table, 0)} rows\t{file_label}")

    if "analysis_run" not in tables:
        print("[WARN] analysis_run.tsv was not found. This is only safe if run_pk is nullable and already exists in SQL.", file=sys.stderr)

    if args.check_only or args.dry_run:
        print("\n[INFO] Load order:")
        for table in tables:
            files = table_input_files(tsv_dir, table)
            print(f"  {table} ({len(files)} file(s))")
        print("\n[OK] No SQL changes made.")
        return

    conn = connect(args)
    try:
        with conn.cursor() as cursor:
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
                    raise RuntimeError("Post-load relationship checks failed:\n  " + "\n  ".join(problems))

            print("[INFO] Database row counts after load:")
            for table in tables:
                print(f"  {table}\t{table_count(cursor, table)} rows")

        conn.commit()
        print("\n[OK] Load committed successfully.")
    except Exception:
        conn.rollback()
        print("\n[ERROR] Load failed. Transaction rolled back.", file=sys.stderr)
        raise
    finally:
        conn.close()


if __name__ == "__main__":
    main()
