#!/usr/bin/env python3
r"""
05_to_sql_compleasm_features.py

Load Compleasm-derived feature TSVs produced by 04_build_compleasm_feature_tsvs_v4.py
into the gc3_dynamics_v6 MySQL database.

Expected pipeline order:
    00_build_starter_sql_tsvs_v14.py
    01_load_starter_sql_tsvs_v3.py
    03_raw_ortholog_validity_batch.py
    04_build_compleasm_feature_tsvs_v4.py
    05_to_sql_compleasm_features.py

Expected input files in --input-dir:
    sauropsida_odb12.tsv
    orthologs.tsv
    ortholog_summary.tsv
    intron_compleasm.tsv
    intron_compleasm_summary.tsv
    flanks_compleasm.tsv
    flank_sets_compleasm.tsv
    flank_summary.tsv

Notes:
    - flank_summary.tsv is loaded into the SQL table flank_compleasm_summary.
    - Extra TSV columns not present in the SQL table are ignored with a warning.
      This allows recovery_status/recovery_notes columns to be loaded if you add
      them to the SQL schema, while remaining compatible with the current v6 DBML.
    - Empty strings, NA, NaN, None, and NULL are converted to MySQL NULL using \N.
    - Boolean text values are normalized to 1/0 for bool/tinyint fields.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pymysql


NULL_STRINGS = {"", "na", "nan", "none", "null", "<na>"}
TRUE_STRINGS = {"true", "t", "1", "yes", "y", "pass", "passed"}
FALSE_STRINGS = {"false", "f", "0", "no", "n", "fail", "failed"}


@dataclass(frozen=True)
class LoadSpec:
    file_name: str
    table_name: str


# FK-safe load order for the Compleasm feature tables.
LOAD_SPECS: List[LoadSpec] = [
    LoadSpec("sauropsida_odb12.tsv", "sauropsida_odb12"),
    LoadSpec("flank_sets_compleasm.tsv", "flank_sets_compleasm"),
    LoadSpec("orthologs.tsv", "orthologs"),
    LoadSpec("ortholog_summary.tsv", "ortholog_summary"),
    LoadSpec("intron_compleasm.tsv", "intron_compleasm"),
    LoadSpec("intron_compleasm_summary.tsv", "intron_compleasm_summary"),
    LoadSpec("flanks_compleasm.tsv", "flanks_compleasm"),
    LoadSpec("flank_summary.tsv", "flank_compleasm_summary"),
]

# Reverse FK-safe order for optional truncation.
TRUNCATE_ORDER = list(reversed([spec.table_name for spec in LOAD_SPECS]))


class LoadError(RuntimeError):
    """Raised for user-facing load/validation errors."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Load Compleasm feature TSVs into gc3_dynamics_v6 MySQL tables."
    )

    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing TSVs from 04_build_compleasm_feature_tsvs_v4.py",
    )
    parser.add_argument("--database", required=True, help="MySQL database name")
    parser.add_argument("--host", default="localhost", help="MySQL host [default: localhost]")
    parser.add_argument("--port", type=int, default=3306, help="MySQL port [default: 3306]")
    parser.add_argument("--user", required=True, help="MySQL user")
    parser.add_argument(
        "--password",
        default=None,
        help="MySQL password. If omitted, MYSQL_PWD is used if available, otherwise you will be prompted.",
    )
    parser.add_argument(
        "--mode",
        choices=["insert", "ignore", "replace"],
        default="insert",
        help="LOAD DATA behavior for duplicate keys [default: insert]",
    )
    parser.add_argument(
        "--truncate-compleasm-tables",
        action="store_true",
        help="Truncate only the Compleasm feature tables before loading. Uses FOREIGN_KEY_CHECKS=0.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate files, headers, and table columns, but do not load data.",
    )
    parser.add_argument(
        "--keep-cleaned-tsvs",
        action="store_true",
        help="Keep temporary cleaned TSVs for debugging instead of deleting them.",
    )
    parser.add_argument(
        "--local-infile",
        action="store_true",
        help="Enable LOAD DATA LOCAL INFILE. Usually required for local TSV loading.",
    )

    return parser.parse_args()


def get_password(args: argparse.Namespace) -> str:
    if args.password is not None:
        return args.password
    env_pwd = os.environ.get("MYSQL_PWD")
    if env_pwd is not None:
        return env_pwd
    import getpass

    return getpass.getpass("MySQL password: ")


def connect(args: argparse.Namespace) -> pymysql.connections.Connection:
    return pymysql.connect(
        host=args.host,
        port=args.port,
        user=args.user,
        password=get_password(args),
        database=args.database,
        charset="utf8mb4",
        local_infile=True,
        autocommit=False,
        cursorclass=pymysql.cursors.DictCursor,
    )


def read_header(path: Path) -> List[str]:
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise LoadError(f"TSV is empty: {path}") from exc
    return [h.strip() for h in header]


def count_data_rows(path: Path) -> int:
    with open(path, "r", newline="") as handle:
        # Subtract header line.
        return max(sum(1 for _ in handle) - 1, 0)


def mysql_identifier(name: str) -> str:
    if not name.replace("_", "").isalnum():
        raise LoadError(f"Unsafe MySQL identifier: {name}")
    return f"`{name}`"


def get_table_columns(conn, table_name: str) -> List[str]:
    with conn.cursor() as cur:
        cur.execute(f"DESCRIBE {mysql_identifier(table_name)}")
        rows = cur.fetchall()
    if not rows:
        raise LoadError(f"Table does not exist or has no columns: {table_name}")
    return [row["Field"] for row in rows]


def get_table_column_types(conn, table_name: str) -> Dict[str, str]:
    with conn.cursor() as cur:
        cur.execute(f"DESCRIBE {mysql_identifier(table_name)}")
        rows = cur.fetchall()
    return {row["Field"]: str(row["Type"]).lower() for row in rows}


def is_boolish_column(column_name: str, mysql_type: str) -> bool:
    if column_name in {"passes_raw_cds_qc", "terminal_stop", "keep_flag", "is_current"}:
        return True
    # MySQL BOOL is usually tinyint(1). Be conservative and only auto-normalize
    # known boolean-like names or true bool/boolean text if present.
    return mysql_type in {"bool", "boolean"}


def clean_value(value: object, column_name: str, mysql_type: str) -> str:
    raw = "" if value is None else str(value).strip()
    if raw.lower() in NULL_STRINGS:
        return r"\N"

    if is_boolish_column(column_name, mysql_type):
        low = raw.lower()
        if low in TRUE_STRINGS:
            return "1"
        if low in FALSE_STRINGS:
            return "0"
        raise LoadError(f"Cannot parse boolean value for {column_name}: {raw!r}")

    return raw


def prepare_clean_tsv(
    source_path: Path,
    cleaned_dir: Path,
    table_name: str,
    table_columns: Sequence[str],
    column_types: Dict[str, str],
) -> Tuple[Path, List[str], List[str], List[str]]:
    """
    Create a temporary TSV with only columns present in the SQL table.

    Returns:
        cleaned_path, load_columns, ignored_tsv_columns, missing_required_table_columns
    """
    source_header = read_header(source_path)
    source_set = set(source_header)
    table_set = set(table_columns)

    load_columns = [col for col in table_columns if col in source_set]
    ignored_columns = [col for col in source_header if col not in table_set]

    if not load_columns:
        raise LoadError(
            f"No overlapping columns between {source_path.name} and SQL table {table_name}"
        )

    cleaned_path = cleaned_dir / f"{table_name}.cleaned.tsv"

    with open(source_path, "r", newline="") as src, open(cleaned_path, "w", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        writer = csv.DictWriter(
            dst,
            fieldnames=load_columns,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in reader:
            cleaned = {
                col: clean_value(row.get(col), col, column_types.get(col, ""))
                for col in load_columns
            }
            writer.writerow(cleaned)

    # Informational only. We do not try to infer SQL nullability here because
    # DESCRIBE output can vary by MySQL version and constraints are already
    # enforced by MySQL during LOAD DATA.
    missing_table_columns = [col for col in table_columns if col not in source_set]

    return cleaned_path, load_columns, ignored_columns, missing_table_columns


def validate_input_files(input_dir: Path) -> None:
    missing = [spec.file_name for spec in LOAD_SPECS if not (input_dir / spec.file_name).exists()]
    if missing:
        raise LoadError(
            "Missing required Compleasm feature TSV(s) in --input-dir:\n  "
            + "\n  ".join(missing)
        )


def truncate_tables(conn) -> None:
    with conn.cursor() as cur:
        cur.execute("SET FOREIGN_KEY_CHECKS=0")
        for table in TRUNCATE_ORDER:
            print(f"Truncating {table}")
            cur.execute(f"TRUNCATE TABLE {mysql_identifier(table)}")
        cur.execute("SET FOREIGN_KEY_CHECKS=1")


def load_cleaned_tsv(
    conn,
    cleaned_path: Path,
    table_name: str,
    load_columns: Sequence[str],
    mode: str,
) -> int:
    priority = {
        "insert": "",
        "ignore": "IGNORE",
        "replace": "REPLACE",
    }[mode]

    columns_sql = ", ".join(mysql_identifier(col) for col in load_columns)
    # PyMySQL parameter binding cannot bind identifiers, but can bind the file path.
    sql = f"""
        LOAD DATA LOCAL INFILE %s
        {priority}
        INTO TABLE {mysql_identifier(table_name)}
        CHARACTER SET utf8mb4
        FIELDS TERMINATED BY '\\t'
        LINES TERMINATED BY '\\n'
        IGNORE 1 LINES
        ({columns_sql})
    """

    with conn.cursor() as cur:
        cur.execute(sql, (str(cleaned_path),))
        return cur.rowcount


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir).resolve()

    if not input_dir.exists():
        raise LoadError(f"--input-dir does not exist: {input_dir}")

    validate_input_files(input_dir)

    conn = connect(args)
    cleaned_paths: List[Path] = []

    try:
        if args.truncate_compleasm_tables and not args.dry_run:
            truncate_tables(conn)

        with tempfile.TemporaryDirectory(prefix="compleasm_sql_clean_") as tmp:
            cleaned_dir = Path(tmp)

            print(f"Input directory: {input_dir}")
            print(f"Database: {args.database}")
            print(f"Mode: {args.mode}")
            if args.dry_run:
                print("DRY RUN: validating only; no rows will be loaded")
            print("")

            for spec in LOAD_SPECS:
                source_path = input_dir / spec.file_name
                table_name = spec.table_name
                table_columns = get_table_columns(conn, table_name)
                column_types = get_table_column_types(conn, table_name)

                cleaned_path, load_columns, ignored_columns, missing_table_columns = prepare_clean_tsv(
                    source_path=source_path,
                    cleaned_dir=cleaned_dir,
                    table_name=table_name,
                    table_columns=table_columns,
                    column_types=column_types,
                )
                cleaned_paths.append(cleaned_path)

                n_source_rows = count_data_rows(source_path)
                print(f"{spec.file_name} -> {table_name}")
                print(f"  source rows: {n_source_rows}")
                print(f"  loading columns: {', '.join(load_columns)}")

                if ignored_columns:
                    print(f"  ignoring TSV-only columns: {', '.join(ignored_columns)}")
                if missing_table_columns:
                    print(f"  table columns not provided by TSV: {', '.join(missing_table_columns)}")

                if args.dry_run:
                    print("  validated")
                else:
                    loaded = load_cleaned_tsv(
                        conn=conn,
                        cleaned_path=cleaned_path,
                        table_name=table_name,
                        load_columns=load_columns,
                        mode=args.mode,
                    )
                    print(f"  loaded/affected rows: {loaded}")
                print("")

                if args.keep_cleaned_tsvs:
                    kept_path = input_dir / f".{table_name}.cleaned_for_load.tsv"
                    kept_path.write_text(cleaned_path.read_text())
                    print(f"  kept cleaned TSV: {kept_path}")

        if not args.dry_run:
            conn.commit()
            print("Finished loading Compleasm feature TSVs. Transaction committed.")
        else:
            conn.rollback()
            print("Dry run finished. No changes made.")

    except Exception:
        conn.rollback()
        print("ERROR: load failed. Transaction rolled back.", file=sys.stderr)
        raise
    finally:
        conn.close()


if __name__ == "__main__":
    try:
        main()
    except LoadError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
