#!/usr/bin/env python3
"""Validate and load the eight Compleasm TSVs for gc3_dynamics_v8_reduced.

Requires Python >= 3.9 and PyMySQL (python -m pip install PyMySQL).
The schema/starter records must already be installed; this script never creates
or alters tables, nor loads/truncates genomes, sequences, or natural_history.

Recommended first run:
  python 05_to_sql_compleasm_features_v2.py \\
    --input-dir /path/to/compleasm_features --mysql-db gc3_test \\
    --mysql-user root --check-only
Then repeat without --check-only to load. Passwords can be supplied through
--mysql-password, MYSQL_PWD, or the interactive password prompt.

--preflight-only reads headers/SQL metadata only. --check-only and --dry-run
make no SQL changes. Full validation is the default. --skip-data-validation
skips relationship/key/recovery checks, but still parses and normalizes every
row before loading; malformed values and row widths are never silently loaded.

Full validation uses a temporary SQLite database and normalized load files on
disk rather than retaining entire datasets/SQL parent tables in Python memory.
Text key comparisons use MySQL's own collation weights (WEIGHT_STRING), so case,
accents and trailing-space rules agree with the target indexes. MySQL 8.x is
recommended; older MySQL/MariaDB collations use the legacy PAD SPACE rule.
Use --temp-dir if the default temporary disk has insufficient free space.

The legacy flank_summary.tsv alias may omit recovery_status/recovery_notes;
these are explicitly mapped to NULL, only when SQL permits NULL. Canonical
headers must contain all v8 columns. Column order may differ; names may not.

DML loads require InnoDB and are committed only after successful checks.
LOCAL INFILE can downgrade MySQL errors to warnings; load warnings therefore
abort/rollback except verified duplicate-key warnings in explicit ignore mode.
REPLACE deletes conflicting rows before inserting replacements and is subject
to normal foreign-key restrictions; it is never selected automatically.
TRUNCATE is irreversible DDL. It is allowed only for the eight target tables,
and is refused when other tables contain references into those targets.
Concurrent writers can invalidate earlier validation; keep the target tables
quiescent during validation/loading. Runtime constraints/warnings remain checked.

Reference: https://dev.mysql.com/doc/refman/8.4/en/load-data.html
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import getpass
import importlib.metadata
import math
import os
from pathlib import Path
import re
import sqlite3
import sys
import tempfile
from typing import Any, Iterable

VERSION = "2.0.0"
COLUMNS = {
    "sauropsida_odb12": "odb12_id name orthodb".split(),
    "orthologs": "ortholog_pk sequence_pk odb12_id passes_raw_cds_qc status strand start end gc gc3 gc4".split(),
    "ortholog_summary": "ortholog_summary_pk genome_pk species_pk n_orthologs callable_bp_total mean_gc weighted_mean_gc mean_gc3 weighted_mean_gc3 mean_gc4 weighted_mean_gc4 sd_gc sd_gc3 sd_gc4 var_gc var_gc3 var_gc4 median_gc median_gc3 median_gc4 mad_gc mad_gc3 mad_gc4 iqr_gc iqr_gc3 iqr_gc4 q05_gc3 q25_gc3 q75_gc3 q95_gc3 mean_ortholog_length median_ortholog_length created_at".split(),
    "intron_compleasm": "intron_pk ortholog_pk parent_id intron_id status strand start end length gc".split(),
    "intron_compleasm_summary": "intron_summary_pk genome_pk species_pk n_introns n_orthologs callable_bp_total mean_gc weighted_mean_gc sd_gc var_gc median_gc mad_gc iqr_gc q05_gc q25_gc q75_gc q95_gc mean_intron_length median_intron_length recovery_status recovery_notes created_at".split(),
    "flanks_compleasm": "flank_pk ortholog_pk flank_set_pk status strand up_start up_end down_start down_end gc".split(),
    "flank_sets_compleasm": "flank_set_pk name upstream_bp downstream_bp created_at notes".split(),
    "flank_compleasm_summary": "flank_summary_pk genome_pk species_pk flank_set_pk n_flanks n_orthologs callable_bp_total mean_gc weighted_mean_gc sd_gc var_gc median_gc mad_gc iqr_gc q05_gc q25_gc q75_gc q95_gc mean_flank_length median_flank_length mean_upstream_gc mean_downstream_gc recovery_status recovery_notes created_at".split(),
}
LOAD_ORDER = (
    "sauropsida_odb12", "flank_sets_compleasm", "orthologs", "ortholog_summary",
    "intron_compleasm", "intron_compleasm_summary", "flanks_compleasm",
    "flank_compleasm_summary",
)
PK = {t: COLUMNS[t][0] for t in COLUMNS}
UNIQUE_KEYS = {
    "sauropsida_odb12": [],
    "orthologs": [("sequence_pk", "odb12_id")],
    "ortholog_summary": [("genome_pk",)],
    "intron_compleasm": [("ortholog_pk", "intron_id")],
    "intron_compleasm_summary": [("genome_pk",)],
    "flanks_compleasm": [("ortholog_pk", "flank_set_pk")],
    "flank_sets_compleasm": [("upstream_bp", "downstream_bp")],
    "flank_compleasm_summary": [("genome_pk", "flank_set_pk")],
}
# (source table, source column, referenced table, referenced column)
RELATIONS = [
    ("orthologs", "sequence_pk", "sequences", "sequence_pk"),
    ("orthologs", "odb12_id", "sauropsida_odb12", "odb12_id"),
    ("intron_compleasm", "ortholog_pk", "orthologs", "ortholog_pk"),
    ("flanks_compleasm", "ortholog_pk", "orthologs", "ortholog_pk"),
    ("flanks_compleasm", "flank_set_pk", "flank_sets_compleasm", "flank_set_pk"),
    ("flank_compleasm_summary", "flank_set_pk", "flank_sets_compleasm", "flank_set_pk"),
]
for _table in ("ortholog_summary", "intron_compleasm_summary", "flank_compleasm_summary"):
    RELATIONS.extend([
        (_table, "genome_pk", "genomes", "genome_pk"),
        (_table, "species_pk", "natural_history", "species_pk"),
    ])
BIGINT_COLS = set(PK.values()) - {"odb12_id"}
BIGINT_COLS.update(("sequence_pk", "ortholog_pk", "genome_pk", "species_pk", "flank_set_pk", "callable_bp_total"))
INT_COLS = set("start end length up_start up_end down_start down_end upstream_bp downstream_bp n_orthologs n_introns n_flanks".split())
TEXT_COLS = set("odb12_id name orthodb status strand parent_id intron_id recovery_status recovery_notes notes".split())
NULL_WORDS = {"", "\\n", "na", "nan", "none", "null"}  # lowercased literal \\N
TRUE_WORDS = {"1", "true", "t", "yes", "y"}
FALSE_WORDS = {"0", "false", "f", "no", "n"}
MAX_EXAMPLES = 10


class LoaderError(Exception):
    """A user-facing failure; SQL has not been committed."""


def log(level: str, message: str) -> None:
    print(f"[{level}] {message}", flush=True)


def ident(name: str) -> str:
    return "`" + name.replace("`", "``") + "`"


def safe_name(name: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9_]+", name):
        raise LoaderError(f"Unsupported SQL charset/collation identifier: {name!r}")
    return name


def query(conn: Any, sql: str, params: Iterable[Any] = ()) -> list[dict]:
    with conn.cursor() as cur:
        cur.execute(sql, tuple(params))
        return list(cur.fetchall())


def command(conn: Any, sql: str, params: Iterable[Any] = ()) -> int:
    with conn.cursor() as cur:
        return cur.execute(sql, tuple(params))


@dataclass(frozen=True)
class InputFile:
    path: Path
    header: tuple[str, ...]
    fingerprint: tuple[int, int, int, int]


def fingerprint(path: Path) -> tuple[int, int, int, int]:
    st = path.stat()
    return st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns


def reader(handle: Any) -> Any:
    # Matches csv.DictWriter(delimiter='\t') in the v4/v5 builder.
    return csv.reader(handle, delimiter="\t", quotechar='"', strict=True)


def inspect_inputs(input_dir: Path) -> dict[str, InputFile]:
    if not input_dir.is_dir() or not os.access(input_dir, os.R_OK | os.X_OK):
        raise LoaderError(f"Input directory does not exist/is not readable: {input_dir}")
    inputs = {}
    canonical = input_dir / "flank_compleasm_summary.tsv"
    alias = input_dir / "flank_summary.tsv"
    if canonical.exists() and alias.exists():
        log("INFO", "Ignoring flank_summary.tsv as a compatibility alias; using flank_compleasm_summary.tsv.")
    for table, expected in COLUMNS.items():
        path = input_dir / f"{table}.tsv"
        is_alias = table == "flank_compleasm_summary" and not path.exists() and alias.exists()
        if is_alias:
            path = alias
            log("WARN", "Canonical flank summary absent; validating flank_summary.tsv fallback.")
        if not path.is_file() or not os.access(path, os.R_OK):
            raise LoaderError(f"{table}: required readable input missing: {path}")
        before = fingerprint(path)
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            header = next(reader(handle), [])
        if len(header) != len(set(header)):
            raise LoaderError(f"{table}: duplicate column names in TSV header: {header}")
        missing, extra = set(expected) - set(header), set(header) - set(expected)
        legacy = is_alias and missing == {"recovery_status", "recovery_notes"} and not extra
        if (missing or extra) and not legacy:
            raise LoaderError(f"{table}: TSV/schema header mismatch; missing={sorted(missing)}, unexpected={sorted(extra)}")
        if legacy:
            log("WARN", "Legacy alias lacks recovery_status/recovery_notes; mapping them to SQL NULL without inventing recovery information.")
        if before != fingerprint(path):
            raise LoaderError(f"Input changed while checking its header: {path}")
        inputs[table] = InputFile(path, tuple(header), before)
    log("OK", "Required Compleasm TSV inputs found and readable")
    log("OK", "TSV headers match gc3_dynamics_v8_reduced (with explicit legacy alias mapping, if used)")
    return inputs


def expected_type(column: str) -> str:
    if column in BIGINT_COLS:
        return "bigint"
    if column in INT_COLS:
        return "int"
    if column == "passes_raw_cds_qc":
        return "bool"
    if column == "created_at":
        return "datetime"
    if column in TEXT_COLS:
        return "string"
    return "double"


def load_sql(table: str, path: Path, mode: str) -> tuple[str, tuple[str, ...]]:
    modifier = {"insert": "", "ignore": "IGNORE ", "replace": "REPLACE "}[mode]
    variables = [f"@v{i}" for i in range(len(COLUMNS[table]))]
    # Blank/null spellings were normalized to MySQL's \\N marker already.
    # Direct assignment avoids collation-sensitive NULLIF changing valid text.
    assignments = ", ".join(f"{ident(c)} = {v}" for c, v in zip(COLUMNS[table], variables))
    # This input is our normalized MySQL-escaped file with no header, not the
    # original CSV-style TSV. No %-formatting of this statement is performed.
    sql = (
        f"LOAD DATA LOCAL INFILE %s {modifier}INTO TABLE {ident(table)} "
        "CHARACTER SET utf8mb4 FIELDS TERMINATED BY %s ESCAPED BY %s "
        "LINES TERMINATED BY %s "
        f"({', '.join(variables)}) SET {assignments}"
    )
    # PyMySQL also quotes these correctly under NO_BACKSLASH_ESCAPES SQL mode.
    return sql, (str(path.resolve()), "\t", "\\", "\n")


@dataclass
class Schema:
    columns: dict[str, dict[str, dict]]
    keys: dict[str, list[tuple[str, ...]]]
    collations: dict[str, str]


def inspect_schema(conn: Any, args: argparse.Namespace, inputs: dict) -> Schema:
    required = set(COLUMNS) | {"genomes", "sequences", "natural_history"}
    marks = ",".join(["%s"] * len(required))
    names = sorted(required)
    tables = query(conn, f"SELECT TABLE_NAME, ENGINE FROM information_schema.TABLES WHERE TABLE_SCHEMA=%s AND TABLE_NAME IN ({marks})", [args.mysql_db] + names)
    present = {r["TABLE_NAME"]: r["ENGINE"] for r in tables}
    if required - present.keys():
        raise LoaderError(f"Missing SQL tables: {sorted(required - present.keys())}. Install the v8 schema and starter data first.")
    for table in COLUMNS:
        if str(present[table]).lower() != "innodb":
            raise LoaderError(f"{table}: InnoDB required for rollback-safe DML; engine={present[table]}")
    rows = query(conn,
        "SELECT TABLE_NAME, COLUMN_NAME, DATA_TYPE, COLUMN_TYPE, IS_NULLABLE, "
        "CHARACTER_MAXIMUM_LENGTH, CHARACTER_SET_NAME, COLLATION_NAME, COLUMN_DEFAULT, EXTRA, DATETIME_PRECISION "
        f"FROM information_schema.COLUMNS WHERE TABLE_SCHEMA=%s AND TABLE_NAME IN ({marks})", [args.mysql_db] + names)
    columns: dict[str, dict[str, dict]] = defaultdict(dict)
    for row in rows:
        columns[row["TABLE_NAME"]][row["COLUMN_NAME"]] = row
    for table, expected in COLUMNS.items():
        missing = set(expected) - columns[table].keys()
        if missing:
            raise LoaderError(f"{table}: required SQL columns missing: {sorted(missing)}")
        # Extra columns are permitted only when omission during INSERT is safe.
        for col, meta in columns[table].items():
            if col not in expected and meta["IS_NULLABLE"] == "NO" and meta["COLUMN_DEFAULT"] is None and not meta["EXTRA"]:
                raise LoaderError(f"{table}.{col}: unexpected required SQL column; schema drift")
        for col in expected:
            meta = columns[table][col]
            kind, actual = expected_type(col), meta["DATA_TYPE"].lower()
            allowed = {"bool": {"tinyint", "bool", "boolean"}, "string": {"varchar", "text", "mediumtext", "longtext"}}.get(kind, {kind})
            if actual not in allowed:
                raise LoaderError(f"{table}.{col}: expected {kind}, SQL has {meta['COLUMN_TYPE']}")
            if "generated" in meta["EXTRA"].lower() and "default_generated" not in meta["EXTRA"].lower():
                raise LoaderError(f"{table}.{col}: generated column cannot accept incoming values")
            if col not in inputs[table].header and meta["IS_NULLABLE"] != "YES":
                raise LoaderError(f"{table}.{col}: legacy alias requires NULL, but SQL forbids it")
    for _, _, table, col in RELATIONS:
        if col not in columns[table]:
            raise LoaderError(f"Parent SQL column missing: {table}.{col}")
    if columns["orthologs"]["sequence_pk"]["IS_NULLABLE"] != "YES":
        raise LoaderError("orthologs.sequence_pk must allow NULL, as required by the v8 Compleasm pipeline")
    for source, col, target, ref in RELATIONS:
        left, right = columns[source][col], columns[target][ref]
        left_unsigned = "unsigned" in left["COLUMN_TYPE"].lower()
        right_unsigned = "unsigned" in right["COLUMN_TYPE"].lower()
        if (left["DATA_TYPE"] != right["DATA_TYPE"] or left["COLLATION_NAME"] != right["COLLATION_NAME"]
                or left_unsigned != right_unsigned):
            raise LoaderError(f"Incompatible relationship types/collations: {source}.{col} -> {target}.{ref}")
    index_rows = query(conn,
        "SELECT TABLE_NAME, INDEX_NAME, SEQ_IN_INDEX, COLUMN_NAME, SUB_PART "
        f"FROM information_schema.STATISTICS WHERE TABLE_SCHEMA=%s AND NON_UNIQUE=0 AND TABLE_NAME IN ({','.join(['%s'] * len(COLUMNS))}) "
        "ORDER BY TABLE_NAME, INDEX_NAME, SEQ_IN_INDEX", [args.mysql_db] + list(COLUMNS))
    indexes: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    for row in index_rows:
        if row["SUB_PART"] is not None or row["COLUMN_NAME"] not in COLUMNS[row["TABLE_NAME"]]:
            raise LoaderError(f"{row['TABLE_NAME']}.{row['INDEX_NAME']}: unsupported prefix/expression/extra-column unique index")
        indexes[row["TABLE_NAME"]][row["INDEX_NAME"]].append(row["COLUMN_NAME"])
    keys = {}
    for table in COLUMNS:
        if indexes[table].get("PRIMARY") != [PK[table]]:
            raise LoaderError(f"{table}: expected primary key ({PK[table]})")
        found = {tuple(cols) for cols in indexes[table].values()}
        for key in UNIQUE_KEYS[table]:
            if key not in found:
                raise LoaderError(f"{table}: missing required UNIQUE index {key}")
        for key in found:
            if any(expected_type(c) not in {"bigint", "int", "bool", "string"} for c in key):
                raise LoaderError(f"{table}: unexpected UNIQUE key {key} uses a type outside the v8 key definitions")
        keys[table] = [(PK[table],)] + sorted(found - {(PK[table],)})
    # Read all declared FKs touching targets. Reject unexpected outbound FKs;
    # inbound FKs are considered again before destructive TRUNCATE.
    declared = query(conn,
        "SELECT TABLE_NAME, COLUMN_NAME, REFERENCED_TABLE_SCHEMA, REFERENCED_TABLE_NAME, REFERENCED_COLUMN_NAME "
        f"FROM information_schema.KEY_COLUMN_USAGE WHERE TABLE_SCHEMA=%s AND REFERENCED_TABLE_NAME IS NOT NULL AND TABLE_NAME IN ({','.join(['%s'] * len(COLUMNS))})",
        [args.mysql_db] + list(COLUMNS))
    expected_relations = set(RELATIONS)
    for row in declared:
        relation = (row["TABLE_NAME"], row["COLUMN_NAME"], row["REFERENCED_TABLE_NAME"], row["REFERENCED_COLUMN_NAME"])
        if row["REFERENCED_TABLE_SCHEMA"] != args.mysql_db or relation not in expected_relations:
            raise LoaderError(f"Unexpected SQL foreign key outside v8 Compleasm relationships: {relation}")
    # A trigger could modify tables outside our requested scope.
    triggers = query(conn, f"SELECT TRIGGER_NAME, EVENT_OBJECT_TABLE FROM information_schema.TRIGGERS WHERE TRIGGER_SCHEMA=%s AND EVENT_OBJECT_TABLE IN ({','.join(['%s'] * len(COLUMNS))})", [args.mysql_db] + list(COLUMNS))
    if triggers:
        raise LoaderError(f"Target table triggers require separate review before loading: {triggers[:MAX_EXAMPLES]}")
    collation_names = sorted({m["COLLATION_NAME"] for t in COLUMNS for m in columns[t].values() if m["COLLATION_NAME"]})
    collations = {}
    if collation_names:
        has_pad = query(conn, "SELECT COLUMN_NAME FROM information_schema.COLUMNS WHERE TABLE_SCHEMA='information_schema' AND TABLE_NAME='COLLATIONS' AND COLUMN_NAME='PAD_ATTRIBUTE'")
        pad_select = "PAD_ATTRIBUTE" if has_pad else "'PAD SPACE' AS PAD_ATTRIBUTE"
        rows = query(conn, f"SELECT COLLATION_NAME, {pad_select} FROM information_schema.COLLATIONS WHERE COLLATION_NAME IN ({','.join(['%s'] * len(collation_names))})", collation_names)
        collations = {r["COLLATION_NAME"]: r["PAD_ATTRIBUTE"] for r in rows}
        if set(collation_names) - collations.keys():
            raise LoaderError("Cannot determine target SQL collation rules")
        tested_collations = set()
        for table in COLUMNS:
            for col, meta in columns[table].items():
                if col in COLUMNS[table] and meta["COLLATION_NAME"] and meta["COLLATION_NAME"] not in tested_collations:
                    charset, coll = safe_name(meta["CHARACTER_SET_NAME"]), safe_name(meta["COLLATION_NAME"])
                    query(conn, f"SELECT HEX(WEIGHT_STRING(CONVERT(%s USING {charset}) COLLATE {coll})) AS w", ["preflight"])
                    tested_collations.add(coll)
    log("OK", "Required SQL tables, columns, primary/unique indexes, and parent relationships found")
    return Schema(dict(columns), keys, collations)


def preflight(args: argparse.Namespace, workdir: Path) -> tuple[Any, dict, Schema]:
    log("INFO", "Running fast preflight checks...")
    log("OK", f"Python executable: {sys.executable} ({sys.version.split()[0]})")
    env = os.environ.get("CONDA_DEFAULT_ENV") or os.environ.get("VIRTUAL_ENV") or "system/default"
    log("OK", f"Python environment: {env}")
    try:
        import pymysql
        from pymysql.constants import CLIENT
    except ImportError as exc:
        raise LoaderError("PyMySQL missing in this Python environment. Install with: python -m pip install PyMySQL") from exc
    try:
        dep_version = importlib.metadata.version("PyMySQL")
    except importlib.metadata.PackageNotFoundError:
        dep_version = getattr(pymysql, "__version__", "unknown")
    log("OK", f"Python dependency: pymysql {dep_version}")
    inputs = inspect_inputs(Path(args.input_dir).expanduser().resolve())
    password = args.mysql_password
    if password is None:
        password = getpass.getpass("MySQL password (blank allowed): ") if sys.stdin.isatty() else ""
    try:
        conn = pymysql.connect(
            host=args.mysql_host, port=args.mysql_port, user=args.mysql_user,
            password=password, database=args.mysql_db, charset="utf8mb4",
            local_infile=True, autocommit=False, connect_timeout=args.connect_timeout,
            cursorclass=pymysql.cursors.DictCursor,
        )
    except Exception as exc:
        raise LoaderError(f"MySQL connection failed: {exc}") from exc
    try:
        row = query(conn, "SELECT DATABASE() AS db, @@GLOBAL.local_infile AS local_infile, @@SESSION.foreign_key_checks AS fk_checks")[0]
        if row["db"] != args.mysql_db:
            raise LoaderError(f"Selected database is {row['db']!r}, expected {args.mysql_db!r}")
        log("OK", f"MySQL connection successful: {args.mysql_user}@{args.mysql_host}:{args.mysql_port}/{args.mysql_db}")
        if not conn.client_flag & CLIENT.LOCAL_FILES:
            raise LoaderError("PyMySQL client LOCAL INFILE capability is not enabled")
        if str(row["local_infile"]).lower() not in {"1", "on", "true"}:
            raise LoaderError("MySQL server local_infile is disabled. Enable it on the server, then retry; the client is already enabled.")
        if not row["fk_checks"]:
            raise LoaderError("SQL session foreign_key_checks is disabled; enable constraints before loading")
        log("OK", "LOCAL INFILE enabled for client/server; SQL foreign-key checking enabled")
        schema = inspect_schema(conn, args, inputs)
        with conn.cursor() as cur:
            for table in LOAD_ORDER:
                sql, params = load_sql(table, workdir / f"{table}.load.tsv", args.mode)
                cur.mogrify(sql, params)
        log("OK", "LOAD DATA SQL parameterization check passed")
        log("OK", "Fast preflight passed")
        return conn, inputs, schema
    except BaseException:
        conn.close()
        raise


class Problems:
    def __init__(self) -> None:
        self.items: dict[tuple[str, str], dict] = {}

    def add(self, table: str, field: str, count: int, examples: Iterable[Any]) -> None:
        if not count:
            return
        item = self.items.setdefault((table, field), {"count": 0, "examples": []})
        item["count"] += count
        item["examples"].extend(list(examples)[:max(0, MAX_EXAMPLES - len(item["examples"]))])

    def fail_if_any(self) -> None:
        if self.items:
            for (table, field), item in self.items.items():
                log("ERROR", f"{table}: {field}; problematic rows={item['count']}; examples={item['examples']}")
            raise LoaderError(f"Validation failed in {len(self.items)} category/categories. No loading performed.")


def normalize(value: str, meta: dict, column: str) -> str | None:
    stripped = value.strip()
    if stripped.lower() in NULL_WORDS:
        if meta["IS_NULLABLE"] != "YES":
            raise ValueError("NULL/blank is forbidden by the SQL schema")
        return None
    kind = expected_type(column)
    if kind in {"bigint", "int"}:
        if not re.fullmatch(r"[+-]?[0-9]+", stripped):
            raise ValueError("expected an integer")
        integer = int(stripped)
        bits = 64 if kind == "bigint" else 32
        unsigned = "unsigned" in meta["COLUMN_TYPE"].lower()
        low, high = (0, 2**bits - 1) if unsigned else (-2**(bits-1), 2**(bits-1) - 1)
        if not low <= integer <= high:
            raise ValueError(f"integer outside SQL range [{low}, {high}]")
        return str(integer)
    if kind == "bool":
        if stripped.lower() in TRUE_WORDS:
            return "1"
        if stripped.lower() in FALSE_WORDS:
            return "0"
        raise ValueError("expected boolean true/false or 1/0")
    if kind == "double":
        try:
            if not re.fullmatch(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?", stripped):
                raise ValueError("invalid numeric spelling")
            number = float(stripped)
            if not math.isfinite(number):
                raise ValueError("non-finite number")
            # Preserve valid decimal/scientific notation, rather than round it.
            decimal = Decimal(stripped)
            if number == 0 and decimal != 0:
                raise ValueError("value underflows SQL DOUBLE")
        except (ValueError, InvalidOperation, OverflowError) as exc:
            raise ValueError("expected a finite number") from exc
        return stripped
    if kind == "datetime":
        try:
            text = stripped[:-1] + "+00:00" if stripped.endswith(("Z", "z")) else stripped
            date = datetime.fromisoformat(text)
            if not re.search(r"[T ]\d{2}:\d{2}:\d{2}", stripped):
                raise ValueError("timestamp must include date and HH:MM:SS")
            if date.tzinfo is not None:
                date = date.astimezone(timezone.utc).replace(tzinfo=None)
            if date.year < 1000:
                raise ValueError("MySQL DATETIME requires year >= 1000")
            precision = int(meta.get("DATETIME_PRECISION") or 0)
            if date.microsecond and date.microsecond % (10 ** (6 - precision)):
                raise ValueError(f"fractional seconds exceed SQL DATETIME({precision}) precision")
            return date.strftime("%Y-%m-%d %H:%M:%S") + (f".{date.microsecond:06d}"[:precision+1] if precision else "")
        except (ValueError, OverflowError) as exc:
            raise ValueError(f"invalid datetime: {exc}") from exc
    limit = meta["CHARACTER_MAXIMUM_LENGTH"]
    if limit is not None and len(value) > int(limit):
        raise ValueError(f"text length {len(value)} exceeds SQL limit {limit}")
    return value  # preserve whitespace, quotes and backslashes in text


def mysql_field(value: str | None) -> str:
    if value is None:
        return "\\N"
    return (value.replace("\\", "\\\\").replace("\0", "\\0").replace("\t", "\\t")
            .replace("\n", "\\n").replace("\r", "\\r").replace("\x1a", "\\Z"))


def local_columns(table: str, schema: Schema) -> list[str]:
    needed = {col for key in schema.keys[table] for col in key}
    needed.update(col for source, col, _, _ in RELATIONS if source == table)
    if table == "orthologs":
        needed.add("odb12_id")
    if table == "intron_compleasm":
        needed.add("parent_id")
    if table in {"intron_compleasm_summary", "flank_compleasm_summary"}:
        needed.update({"recovery_status", "n_orthologs", "callable_bp_total", "n_introns" if table.startswith("intron") else "n_flanks"})
    return [c for c in COLUMNS[table] if c in needed]


def stage_files(inputs: dict, schema: Schema, workdir: Path, full: bool) -> tuple[sqlite3.Connection | None, dict[str, int]]:
    log("INFO", "Parsing/normalizing TSV rows" + (" and building disk-backed validation indexes..." if full else " (relationship/key validation skipped)..."))
    db = sqlite3.connect(workdir / "validation.sqlite") if full else None
    if db:
        db.execute("PRAGMA journal_mode=OFF")
        db.execute("PRAGMA synchronous=OFF")
        db.execute("PRAGMA temp_store=FILE")
        db.execute("PRAGMA cache_size=-32768")
    counts = {}
    errors = Problems()
    try:
        for table in LOAD_ORDER:
            source = inputs[table]
            if source.fingerprint != fingerprint(source.path):
                raise LoaderError(f"{table}: input changed since preflight; restart the run")
            kept = local_columns(table, schema)
            if db:
                db.execute(f"CREATE TABLE {ident(table)} (row_no INTEGER PRIMARY KEY, {', '.join(ident(c) + ' TEXT' for c in kept)})")
            count, pending = 0, []
            with source.path.open("r", encoding="utf-8-sig", newline="") as handle, (workdir / f"{table}.load.tsv").open("w", encoding="utf-8", newline="") as output:
                rows = reader(handle)
                if tuple(next(rows, [])) != source.header:
                    raise LoaderError(f"{table}: header changed after preflight")
                for record_no, fields in enumerate(rows, start=2):
                    count += 1
                    if len(fields) != len(source.header):
                        errors.add(table, "TSV row width", 1, [{"row": record_no, "expected": len(source.header), "found": len(fields)}])
                        continue
                    raw = dict(zip(source.header, fields))
                    values, invalid = {}, False
                    for col in COLUMNS[table]:
                        try:
                            values[col] = normalize(raw.get(col, ""), schema.columns[table][col], col)
                        except ValueError as exc:
                            errors.add(table, col, 1, [{"row": record_no, "value": raw.get(col), "error": str(exc)}])
                            invalid = True
                    if invalid:
                        continue
                    output.write("\t".join(mysql_field(values[c]) for c in COLUMNS[table]) + "\n")
                    if db:
                        pending.append((record_no, *(values[c] for c in kept)))
                        if len(pending) >= 1000:
                            db.executemany(f"INSERT INTO {ident(table)} VALUES ({','.join(['?'] * (len(kept)+1))})", pending)
                            pending.clear()
                    if count % 250000 == 0:
                        log("INFO", f"{table}: parsed {count:,} rows")
            if db and pending:
                db.executemany(f"INSERT INTO {ident(table)} VALUES ({','.join(['?'] * (len(kept)+1))})", pending)
            if source.fingerprint != fingerprint(source.path):
                raise LoaderError(f"{table}: input changed while parsing; restart the run")
            counts[table] = count
            log("OK", f"{table}: {count:,} incoming rows")
        errors.fail_if_any()
        if db:
            db.commit()
        return db, counts
    except BaseException:
        if db:
            db.close()
        raise


def comparison_columns(table: str, schema: Schema) -> list[str]:
    result = {col for key in schema.keys[table] for col in key}
    result.update(col for source, col, _, _ in RELATIONS if source == table)
    return [c for c in COLUMNS[table] if c in result]


def add_comparison_tokens(db: sqlite3.Connection, conn: Any, schema: Schema, batch_size: int) -> None:
    """Obtain exact text-key weights from MySQL using read-only SELECTs."""
    db.execute("CREATE TABLE weights (coll TEXT, raw TEXT, token TEXT, PRIMARY KEY(coll,raw)) WITHOUT ROWID")
    metadata = {}
    for table in LOAD_ORDER:
        for col in comparison_columns(table, schema):
            meta = schema.columns[table][col]
            db.execute(f"ALTER TABLE {ident(table)} ADD COLUMN {ident('k_' + col)} TEXT")
            if meta["COLLATION_NAME"]:
                coll = meta["COLLATION_NAME"]
                metadata[coll] = meta
                db.execute(f"INSERT OR IGNORE INTO weights(coll,raw) SELECT ?, {ident(col)} FROM {ident(table)} WHERE {ident(col)} IS NOT NULL", [coll])
    for coll, meta in metadata.items():
        charset = safe_name(meta["CHARACTER_SET_NAME"])
        safe_name(coll)
        cursor = db.execute("SELECT raw FROM weights WHERE coll=? AND token IS NULL", [coll])
        while True:
            raw_values = [r[0] for r in cursor.fetchmany(batch_size)]
            if not raw_values:
                break
            pad_space = schema.collations[coll] == "PAD SPACE"
            value = "RTRIM(CONVERT(%s USING " + charset + "))" if pad_space else "CONVERT(%s USING " + charset + ")"
            select = f"SELECT %s AS raw, HEX(WEIGHT_STRING({value} COLLATE {coll})) AS token"
            params = [v for raw in raw_values for v in (raw, raw)]
            results = query(conn, " UNION ALL ".join([select] * len(raw_values)), params)
            if len(results) != len(raw_values) or any(r["token"] is None for r in results):
                raise LoaderError(f"Cannot compute MySQL key comparison weights for {coll}")
            db.executemany("UPDATE weights SET token=? WHERE coll=? AND raw=?", [(r["token"], coll, r["raw"]) for r in results])
    for table in LOAD_ORDER:
        for col in comparison_columns(table, schema):
            coll = schema.columns[table][col]["COLLATION_NAME"]
            target = ident("k_" + col)
            if coll:
                db.execute(f"UPDATE {ident(table)} SET {target}=(SELECT token FROM weights WHERE coll=? AND raw={ident(table)}.{ident(col)})", [coll])
            else:
                db.execute(f"UPDATE {ident(table)} SET {target}={ident(col)}")
            db.execute(f"CREATE INDEX {ident('idx_' + table + '_' + col)} ON {ident(table)} ({target})")
    db.commit()


def local_error(db: sqlite3.Connection, errors: Problems, table: str, label: str, sql: str, params: Iterable[Any] = ()) -> int:
    params = tuple(params)
    count = db.execute(f"SELECT COUNT(*) FROM ({sql})", params).fetchone()[0]
    if count:
        cursor = db.execute(sql + f" LIMIT {MAX_EXAMPLES}", params)
        names = [d[0] for d in cursor.description]
        errors.add(table, label, count, [dict(zip(names, row)) for row in cursor.fetchall()])
    return count


def validate_local(db: sqlite3.Connection, schema: Schema, errors: Problems) -> None:
    for table, keys in schema.keys.items():
        for key in keys:
            token = [ident("k_" + c) for c in key]
            nonnull = " AND ".join(c + " IS NOT NULL" for c in token)
            # MySQL UNIQUE permits repeated keys containing any NULL component.
            groups = f"SELECT {', '.join(token)} FROM {ident(table)} WHERE {nonnull} GROUP BY {', '.join(token)} HAVING COUNT(*)>1"
            join = " AND ".join(f"a.{c}=d.{c}" for c in token)
            selected = ",".join("a." + ident(c) for c in key)
            sql = f"SELECT a.row_no,{selected} FROM {ident(table)} a JOIN ({groups}) d ON {join}"
            local_error(db, errors, table, f"duplicate {'primary' if key == (PK[table],) else 'unique'} key {key}", sql)
    for source, col, target, ref in RELATIONS:
        if target not in COLUMNS:
            continue
        sql = (f"SELECT s.row_no,s.{ident(col)} FROM {ident(source)} s WHERE s.{ident(col)} IS NOT NULL AND NOT EXISTS "
               f"(SELECT 1 FROM {ident(target)} p WHERE p.{ident('k_' + ref)}=s.{ident('k_' + col)})")
        local_error(db, errors, source, f"missing incoming FK {col} -> {target}.{ref}", sql)
    local_error(db, errors, "intron_compleasm", "parent_id differs from referenced ortholog odb12_id",
        "SELECT i.row_no,i.ortholog_pk,i.parent_id,o.odb12_id FROM intron_compleasm i JOIN orthologs o ON i.k_ortholog_pk=o.k_ortholog_pk "
        "WHERE i.parent_id IS NOT NULL AND o.odb12_id IS NOT NULL AND i.parent_id<>o.odb12_id")
    for table, count_col in [("intron_compleasm_summary", "n_introns"), ("flank_compleasm_summary", "n_flanks")]:
        cols = [count_col, "n_orthologs", "callable_bp_total"]
        condition = " OR ".join(f"{ident(c)} IS NULL OR {ident(c)}<>'0'" for c in cols)
        local_error(db, errors, table, "not_recoverable_from_compleasm requires zero feature/ortholog/callable counts",
            f"SELECT row_no,recovery_status,{','.join(ident(c) for c in cols)} FROM {ident(table)} WHERE recovery_status='not_recoverable_from_compleasm' AND ({condition})")


def existing_matches(conn: Any, schema: Schema, table: str, key: tuple[str, ...], records: list[tuple]) -> set[int]:
    """Read only matched incoming row numbers, using SQL's native comparison.

    UNION batches avoid reading the entire parent/target table. Explicit casts
    prevent BIGINT keys becoming floating-point comparisons to string parameters.
    """
    selects, params = [], []
    for record in records:
        fields = ["%s AS row_no"]
        params.append(record[0])
        for i, col in enumerate(key):
            meta = schema.columns[table][col]
            if meta["COLLATION_NAME"]:
                charset, coll = safe_name(meta["CHARACTER_SET_NAME"]), safe_name(meta["COLLATION_NAME"])
                cast = f"CONVERT(%s USING {charset}) COLLATE {coll}"
            else:
                cast = "CAST(%s AS UNSIGNED)" if "unsigned" in meta["COLUMN_TYPE"].lower() else "CAST(%s AS SIGNED)"
            fields.append(f"{cast} AS {ident('v' + str(i))}")
            params.append(record[i+1])
        selects.append("SELECT " + ",".join(fields))
    condition = " AND ".join(f"p.{ident(col)}=i.{ident('v' + str(i))}" for i, col in enumerate(key))
    sql = f"SELECT DISTINCT i.row_no FROM ({' UNION ALL '.join(selects)}) i JOIN {ident(table)} p ON {condition}"
    return {int(r["row_no"]) for r in query(conn, sql, params)}


def validate_existing(db: sqlite3.Connection, conn: Any, schema: Schema, args: argparse.Namespace, errors: Problems) -> None:
    for source, col, parent, ref in RELATIONS:
        if parent in COLUMNS:
            continue
        # Group by canonical key token; query only distinct referenced values.
        cursor = db.execute(f"SELECT MIN(row_no),{ident(col)},COUNT(*) FROM {ident(source)} WHERE {ident(col)} IS NOT NULL GROUP BY {ident('k_' + col)}")
        missing_count, examples = 0, []
        while True:
            records = cursor.fetchmany(args.batch_size)
            if not records:
                break
            found = existing_matches(conn, schema, parent, (ref,), [(r[0], r[1]) for r in records])
            for row_no, value, count in records:
                if row_no not in found:
                    missing_count += count
                    if len(examples) < MAX_EXAMPLES:
                        examples.append({"row": row_no, col: value, "rows_with_value": count})
        errors.add(source, f"missing SQL FK {col} -> {parent}.{ref}", missing_count, examples)
        if not missing_count:
            log("OK", f"{source}.{col} -> existing {parent}.{ref}")
    if args.truncate_compleasm_tables:
        log("INFO", "Existing Compleasm key conflicts do not apply: requested TRUNCATE clears those eight tables before loading.")
        return
    for table, keys in schema.keys.items():
        for key in keys:
            cols = ",".join(ident(c) for c in key)
            nonnull = " AND ".join(f"{ident(c)} IS NOT NULL" for c in key)
            cursor = db.execute(f"SELECT row_no,{cols} FROM {ident(table)} WHERE {nonnull}")
            conflict_count, examples = 0, []
            while True:
                records = cursor.fetchmany(args.batch_size)
                if not records:
                    break
                found = existing_matches(conn, schema, table, key, records)
                conflict_count += len(found)
                for record in records:
                    if record[0] in found and len(examples) < MAX_EXAMPLES:
                        examples.append({"row": record[0], **dict(zip(key, record[1:]))})
            if conflict_count and args.mode == "insert":
                errors.add(table, f"existing SQL key conflict {key}", conflict_count, examples)
            elif conflict_count:
                action = "skip conflicting incoming rows" if args.mode == "ignore" else "delete conflicting existing rows and insert replacements"
                log("WARN", f"{table} key {key}: {conflict_count:,} incoming keys already exist; mode={args.mode} will {action}; examples={examples}")
            else:
                log("OK", f"{table} key {key}: no existing SQL conflicts")


def full_validation(db: sqlite3.Connection, conn: Any, schema: Schema, args: argparse.Namespace) -> None:
    log("INFO", "Running full Compleasm TSV relationship/duplicate validation...")
    add_comparison_tokens(db, conn, schema, args.batch_size)
    errors = Problems()
    validate_local(db, schema, errors)
    validate_existing(db, conn, schema, args, errors)
    errors.fail_if_any()
    log("OK", "Incoming FK relationships, intron parent IDs, and recovery counts passed")
    log("OK", "Existing-SQL parent relationships passed")
    log("OK", "Primary/unique key validation passed (NULL-containing UNIQUE keys allowed)")
    log("OK", "Existing SQL key-conflict validation passed for selected load mode")
    log("OK", "Full Compleasm TSV validation passed")


def check_truncate_references(conn: Any, args: argparse.Namespace) -> None:
    """Protect all tables outside the eight targets, including nhPhyML data."""
    rows = query(conn,
        "SELECT TABLE_SCHEMA, TABLE_NAME, COLUMN_NAME, REFERENCED_TABLE_NAME FROM information_schema.KEY_COLUMN_USAGE "
        f"WHERE REFERENCED_TABLE_SCHEMA=%s AND REFERENCED_TABLE_NAME IN ({','.join(['%s'] * len(COLUMNS))})", [args.mysql_db] + list(COLUMNS))
    for row in rows:
        if row["TABLE_SCHEMA"] == args.mysql_db and row["TABLE_NAME"] in COLUMNS:
            continue
        child = ident(row["TABLE_SCHEMA"]) + "." + ident(row["TABLE_NAME"])
        count = query(conn, f"SELECT COUNT(*) AS n FROM {child} WHERE {ident(row['COLUMN_NAME'])} IS NOT NULL")[0]["n"]
        if count:
            raise LoaderError(f"Refusing TRUNCATE: {child}.{row['COLUMN_NAME']} has {count:,} rows referencing {row['REFERENCED_TABLE_NAME']}. Other pipeline data would lose its parent records.")


def truncate_targets(conn: Any, args: argparse.Namespace) -> None:
    # Start a fresh read snapshot before inspecting outside references.
    conn.rollback()
    check_truncate_references(conn, args)
    log("WARN", "TRUNCATE is DDL and cannot be rolled back. Clearing ONLY the eight Compleasm tables.")
    # Even an empty InnoDB child still has an FK definition that blocks a parent
    # TRUNCATE. Temporarily disable checks ONLY for this DDL, restoring immediately.
    conn.rollback()
    command(conn, "SET SESSION FOREIGN_KEY_CHECKS=0")
    try:
        for table in reversed(LOAD_ORDER):
            command(conn, f"TRUNCATE TABLE {ident(table)}")
            log("INFO", f"Truncated {table}")
    except BaseException:
        log("WARN", "TRUNCATE failed; any tables already truncated remain cleared and cannot be restored by rollback.")
        raise
    finally:
        command(conn, "SET SESSION FOREIGN_KEY_CHECKS=1")


def post_checks(conn: Any) -> None:
    errors = []
    for source, col, parent, ref in RELATIONS:
        sql = (f"SELECT COUNT(*) AS n FROM {ident(source)} s LEFT JOIN {ident(parent)} p ON s.{ident(col)}=p.{ident(ref)} "
               f"WHERE s.{ident(col)} IS NOT NULL AND p.{ident(ref)} IS NULL")
        count = query(conn, sql)[0]["n"]
        if count:
            examples = query(conn,
                f"SELECT s.{ident(PK[source])},s.{ident(col)} FROM {ident(source)} s LEFT JOIN {ident(parent)} p ON s.{ident(col)}=p.{ident(ref)} "
                f"WHERE s.{ident(col)} IS NOT NULL AND p.{ident(ref)} IS NULL LIMIT {MAX_EXAMPLES}")
            errors.append(f"{source}.{col} -> {parent}.{ref}: orphan rows={count}, examples={examples}")
    count = query(conn,
        "SELECT COUNT(*) AS n FROM intron_compleasm i JOIN orthologs o ON i.ortholog_pk=o.ortholog_pk "
        "WHERE i.parent_id IS NOT NULL AND o.odb12_id IS NOT NULL AND BINARY i.parent_id<>BINARY o.odb12_id")[0]["n"]
    if count:
        errors.append(f"intron_compleasm.parent_id: {count:,} rows disagree with referenced ortholog odb12_id")
    if errors:
        raise LoaderError("Post-load validation failed: " + "; ".join(errors))
    log("OK", "Post-load SQL orphan checks and intron parent-ID check passed")


def load_tables(conn: Any, args: argparse.Namespace, workdir: Path, counts: dict) -> None:
    log("INFO", "Load order: " + " -> ".join(LOAD_ORDER))
    if args.mode == "replace":
        log("WARN", "Explicit REPLACE mode selected: conflicting existing rows are deleted and replaced; SQL FK constraints remain enabled.")
    if args.mode == "ignore":
        log("WARN", "Explicit IGNORE mode selected: incoming rows conflicting with existing keys are skipped.")
    if args.truncate_compleasm_tables:
        truncate_targets(conn, args)
    # The read-only validation snapshot is over. DML gets its own transaction.
    conn.rollback()
    command(conn, "SET SESSION max_error_count=65535")
    conn.begin()
    try:
        for table in LOAD_ORDER:
            sql, params = load_sql(table, workdir / f"{table}.load.tsv", args.mode)
            affected = command(conn, sql, params)
            # Retrieve warnings immediately; another normal query can reset them.
            warnings = query(conn, "SHOW WARNINGS")
            warning_count = int(query(conn, "SHOW COUNT(*) WARNINGS")[0]["@@session.warning_count"])
            permitted = args.mode == "ignore" and warning_count == len(warnings) and all(int(w["Code"]) == 1062 for w in warnings)
            if warning_count and not permitted:
                raise LoaderError(f"{table}: LOAD DATA reported {warning_count:,} warnings; refusing lossy/skipped/coerced data. Examples={warnings[:MAX_EXAMPLES]}")
            if warnings:
                log("WARN", f"{table}: {warning_count:,} duplicate-key warnings accepted under explicit IGNORE mode; examples={warnings[:MAX_EXAMPLES]}")
            if args.mode == "insert" and affected != counts[table]:
                raise LoaderError(f"{table}: expected {counts[table]:,} inserted rows, LOAD DATA affected {affected:,}")
            log("OK", f"Loaded {table}: incoming={counts[table]:,}, affected={affected:,}")
        if args.no_post_check:
            log("WARN", "Post-load SQL checks skipped by --no-post-check")
        else:
            post_checks(conn)
        final_counts = {t: int(query(conn, f"SELECT COUNT(*) AS n FROM {ident(t)}")[0]["n"]) for t in LOAD_ORDER}
        conn.commit()
        log("OK", "Compleasm load committed successfully")
        for table, count in final_counts.items():
            log("OK", f"{table}: final SQL row count={count:,}")
    except BaseException:
        conn.rollback()
        log("ERROR", "DML transaction rolled back; no load success reported.")
        if args.truncate_compleasm_tables:
            log("WARN", "Prior TRUNCATE operations cannot be undone; target tables may remain empty.")
        raise


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version=VERSION)
    p.add_argument("--input-dir", required=True, help="Directory with all eight logical Compleasm TSVs")
    p.add_argument("--mysql-db", "--database", required=True, dest="mysql_db")
    p.add_argument("--mysql-user", "--user", default="root", dest="mysql_user")
    p.add_argument("--mysql-password", "--password", default=os.environ.get("MYSQL_PWD"), dest="mysql_password", help="Otherwise use MYSQL_PWD or interactive prompt")
    p.add_argument("--mysql-host", "--host", default="localhost", dest="mysql_host")
    p.add_argument("--mysql-port", "--port", type=int, default=3306, dest="mysql_port")
    p.add_argument("--connect-timeout", type=int, default=10)
    p.add_argument("--mode", choices=("insert", "ignore", "replace"), default="insert")
    modes = p.add_mutually_exclusive_group()
    modes.add_argument("--preflight-only", action="store_true", help="Fast environment/header/SQL checks; no full row scan or SQL writes")
    modes.add_argument("--check-only", action="store_true", help="Full validation; no SQL writes")
    modes.add_argument("--dry-run", action="store_true", help="Validate and print intended actions; no SQL writes")
    p.add_argument("--skip-data-validation", action="store_true", help="Skip relationship/key/recovery validation; normalization still reads rows")
    p.add_argument("--truncate-compleasm-tables", action="store_true", help="Irreversible: clear ONLY eight Compleasm tables before loading")
    p.add_argument("--no-post-check", action="store_true", help="Skip post-load orphan checks")
    p.add_argument("--batch-size", type=int, default=250, help="Read-only SQL lookup batch size, 1..1000 (default: 250)")
    p.add_argument("--temp-dir", help="Parent directory for temporary normalized TSVs/SQLite indexes")
    return p


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = parser()
    args = p.parse_args(argv)
    if args.truncate_compleasm_tables and (args.preflight_only or args.check_only or args.dry_run):
        p.error("--truncate-compleasm-tables cannot be combined with --preflight-only, --check-only, or --dry-run")
    if args.skip_data_validation and (args.check_only or args.preflight_only):
        p.error("--skip-data-validation cannot be combined with --check-only or --preflight-only")
    if not 1 <= args.batch_size <= 1000:
        p.error("--batch-size must be in 1..1000")
    if not 1 <= args.mysql_port <= 65535:
        p.error("--mysql-port must be in 1..65535")
    if not 1 <= args.connect_timeout <= 31536000:
        p.error("--connect-timeout must be positive")
    return args


def run(args: argparse.Namespace) -> None:
    conn, db = None, None
    with tempfile.TemporaryDirectory(prefix="compleasm_sql_v2_", dir=args.temp_dir) as temp:
        workdir = Path(temp)
        try:
            conn, inputs, schema = preflight(args, workdir)
            if args.preflight_only:
                log("OK", "Preflight passed. No full TSV scan or SQL changes made.")
                return
            if args.skip_data_validation:
                log("WARN", "Full relationship/duplicate/conflict/recovery validation skipped by --skip-data-validation; use only after a successful --check-only on the same dataset/database.")
            db, counts = stage_files(inputs, schema, workdir, not args.skip_data_validation)
            if db is not None:
                full_validation(db, conn, schema, args)
            if args.check_only:
                log("OK", "Validation passed. No SQL changes made.")
                return
            if args.dry_run:
                log("INFO", "Dry-run intended load order: " + " -> ".join(LOAD_ORDER))
                log("INFO", f"Intended mode={args.mode}; post-load checks={'off' if args.no_post_check else 'on'}")
                log("OK", "Dry-run finished. No SQL changes made.")
                return
            load_tables(conn, args, workdir, counts)
        finally:
            if db is not None:
                db.close()
            if conn is not None:
                try:
                    conn.rollback()  # closes read-only snapshots; committed DML stays committed
                finally:
                    conn.close()


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        # Increase CSV field limits for recovery notes without platform overflow.
        csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
        run(args)
        return 0
    except KeyboardInterrupt:
        log("ERROR", "Interrupted; uncommitted DML rolled back (TRUNCATE cannot be undone).")
        return 130
    except Exception as exc:
        log("ERROR", str(exc))
        return 1


if __name__ == "__main__":
    sys.exit(main())
