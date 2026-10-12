#!/usr/bin/env python3
"""Validate and load the eight Compleasm TSVs for the updated v9 schema.

Requires Python >= 3.9. PyMySQL is needed only for MySQL modes.
This version expects natural_history.max_longevity_years to be DOUBLE.
TSV relationship reports distinguish accepted parent rows, parents rejected
during parsing, and IDs absent from the identifiable source parent keys.
If malformed parent rows have no identifiable key, absence is inconclusive
and reported separately. Every parsing or unresolved-reference error still
fails validation; rejected rows never become usable parent records.

NEW: validate against starter TSVs before uploading or running the 06_ builder:
  python 05_to_sql_compleasm_features_v5.py \\
    --input-dir /path/to/sql_tsvs/compleasm_features \\
    --starter-tsv-dir /path/to/sql_tsvs --tsv-check-only \\
    --check-nhphyml-mappings

TSV-only mode uses Python's standard library and temporary SQLite indexes.
It never connects to MySQL, prompts for a password, or changes source files.
Core offline text keys use exact, case-sensitive comparisons with whitespace
preserved; only NULL spellings and numeric IDs are normalized. Optional nhPhyML
joins trim ODB/accession IDs and trim/lowercase species names like the 06_ script.
Offline success does not establish live SQL row, length, or collation agreement.
Optional nhPhyML readiness checks concern reference mappings, not GC run inputs
or future result files. They are stricter than the 06_ script about unresolved
ortholog genomes and missing focal-genome species mappings, to expose data loss.
A valid NULL sequence_pk remains allowed. Multiple genome+ODB copies are not
silently removed to satisfy a join, including triploid duplicated orthologs.
Starter headers require all v9 columns. The supported legacy species audit
contains the canonical fields plus genome_pk; its presence selects the 06_
direct-genome branch for every audit row, without per-row species fallback.

Python API:
  result = validate_tsvs(input_dir, starter_tsv_dir,
                         check_nhphyml_mappings=True, verbose=False)
  if not result.passed: raise ValueError(result.as_dict())
The API returns counts, source paths, failures, examples and check statistics.
Source file record numbers include the header as record 1; quoted multiline
fields count as one logical record. No pandas or nhPhyML script import is used.

For MySQL modes, install PyMySQL with: python -m pip install PyMySQL.
The schema/starter records must already be installed; this script never creates
or alters tables, nor loads/truncates genomes, sequences, or natural_history.

Recommended first run:
  python 05_to_sql_compleasm_features_v5.py \\
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
headers must contain all v9 columns. Column order may differ; names may not.

V3 added orthologs.genome_pk (the singular name in the supplied schema). It
checks non-NULL genome IDs against genomes.genome_pk and verifies that mapped
sequences belong to the stated genome whenever both genome IDs are known.
NULL sequence_pk remains allowed, and genome_pk nullability follows SQL.
Older ortholog TSVs without genome_pk must be regenerated with the updated
04_ builder. No genome IDs are guessed from coordinates or ortholog IDs.
The DBML includes the genome_pk column but omits its explicit Ref declaration;
the loader still checks that logical relationship before and after loading.

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
from contextlib import contextmanager, nullcontext, redirect_stdout
import csv
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import getpass
import importlib.metadata
import io
import math
import os
from pathlib import Path
import re
import sqlite3
import sys
import tempfile
from typing import Any, Iterable

VERSION = "5.0.0"
COLUMNS = {
    "sauropsida_odb12": "odb12_id name orthodb".split(),
    "orthologs": "ortholog_pk sequence_pk genome_pk odb12_id passes_raw_cds_qc status strand start end gc gc3 gc4".split(),
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
    ("orthologs", "genome_pk", "genomes", "genome_pk"),
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

# All ordinary fields in these v9 DBML tables are nullable; primary keys are not.
# Keep starter definitions separate: these tables are NEVER SQL load targets.
STARTER_COLUMNS = {
    "genomes": "genome_pk species_ncbi species_pk accession_id is_current".split(),
    "sequences": "sequence_pk genome_pk sequence_id sequence_length sequence_type gc".split(),
    "natural_history": "species_pk species_normalized mass_grams_oskyrko mass_grams_title mass_grams_ji max_female_length_svl_mm_oskyrko svl_mm_title largest_clutch_size mean_number_of_eggs_per_clutch number_clutches_per_year annual_number_of_eggs max_longevity_years range_size genome_size ct_min ct_max".split(),
    "species_name_audit": "species_name_audit_pk species_pk source_dataset source_species_name species_normalized match_status".split(),
}
STARTER_KEYS = {
    "genomes": [("genome_pk",), ("accession_id",)],
    "sequences": [("sequence_pk",), ("genome_pk", "sequence_id")],
    "natural_history": [("species_pk",), ("species_normalized",)],
    "species_name_audit": [("species_name_audit_pk",)],
}
STARTER_PK = {t: cols[0] for t, cols in STARTER_COLUMNS.items()}
PK.update(STARTER_PK)
BIGINT_COLS.add("species_name_audit_pk")
INT_COLS.add("sequence_length")
TEXT_COLS.update("species_ncbi accession_id sequence_id sequence_type species_normalized source_dataset source_species_name match_status".split())
STARTER_RELATIONS = [
    ("sequences", "genome_pk", "genomes", "genome_pk"),
    ("genomes", "species_pk", "natural_history", "species_pk"),
]


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


@contextmanager
def expanded_csv_fields():
    """Give the Python API the CLI's field capacity, restoring caller settings."""
    previous = csv.field_size_limit()
    csv.field_size_limit(max(previous, min(sys.maxsize, 2**31 - 1)))
    try:
        yield
    finally:
        csv.field_size_limit(previous)


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
    log("OK", "TSV headers match gc3_dynamics_9_reduced (with explicit legacy alias mapping, if used)")
    return inputs


def expected_type(column: str) -> str:
    if column in BIGINT_COLS:
        return "bigint"
    if column in INT_COLS:
        return "int"
    if column in {"passes_raw_cds_qc", "is_current"}:
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
        raise LoaderError(f"Missing SQL tables: {sorted(required - present.keys())}. Install the v9 schema and starter data first.")
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
    if "genome_pk" not in columns["sequences"]:
        raise LoaderError("Parent SQL column missing: sequences.genome_pk (needed to check ortholog genome assignment)")
    if columns["orthologs"]["sequence_pk"]["IS_NULLABLE"] != "YES":
        raise LoaderError("orthologs.sequence_pk must allow NULL, as required by the v9 Compleasm pipeline")
    type_relations = RELATIONS + [("orthologs", "genome_pk", "sequences", "genome_pk")]
    for source, col, target, ref in type_relations:
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
                raise LoaderError(f"{table}: unexpected UNIQUE key {key} uses a type outside the v9 key definitions")
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
            raise LoaderError(f"Unexpected SQL foreign key outside v9 Compleasm relationships: {relation}")
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
    def __init__(self, inputs: dict[str, InputFile] | None = None) -> None:
        self.items: dict[tuple[str, str], dict] = {}
        self.inputs = inputs or {}

    def add(self, table: str, field: str, count: int, examples: Iterable[Any]) -> None:
        if not count:
            return
        item = self.items.setdefault((table, field), {"count": 0, "examples": []})
        item["count"] += count
        for example in list(examples)[:max(0, MAX_EXAMPLES - len(item["examples"]))]:
            if isinstance(example, dict):
                example = dict(example)
                if table in self.inputs:
                    example.setdefault("file", str(self.inputs[table].path.resolve()))
                if "row_no" in example or "row" in example:
                    example.setdefault("record_number", example.get("row_no", example.get("row")))
            item["examples"].append(example)

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


def columns_for(table: str, schema: Schema) -> list[str]:
    return COLUMNS[table] if table in COLUMNS else list(schema.columns[table])


def local_columns(table: str, schema: Schema) -> list[str]:
    needed = {col for key in schema.keys[table] for col in key}
    needed.update(col for source, col, _, _ in RELATIONS if source == table)
    if table == "orthologs":
        needed.add("odb12_id")
    if table == "intron_compleasm":
        needed.add("parent_id")
    if table in {"intron_compleasm_summary", "flank_compleasm_summary"}:
        needed.update({"recovery_status", "n_orthologs", "callable_bp_total", "n_introns" if table.startswith("intron") else "n_flanks"})
    if table in STARTER_COLUMNS:
        needed.update(set(schema.columns[table]) - {"gc", "sequence_length"})
    return [c for c in columns_for(table, schema) if c in needed]


def stage_files(inputs: dict, schema: Schema, workdir: Path, full: bool, *, write_load_files: bool = True, errors: Problems | None = None) -> tuple[sqlite3.Connection | None, dict[str, int]]:
    log("INFO", "Parsing/normalizing TSV rows" + (" and building disk-backed validation indexes..." if full else " (relationship/key validation skipped)..."))
    db = sqlite3.connect(workdir / "validation.sqlite") if full else None
    if db:
        db.execute("PRAGMA journal_mode=OFF")
        db.execute("PRAGMA synchronous=OFF")
        db.execute("PRAGMA temp_store=FILE")
        db.execute("PRAGMA cache_size=-32768")
    track_rejections = db is not None and not write_load_files
    if track_rejections:
        db.execute("CREATE TABLE rejected_rows (table_name TEXT, row_no INTEGER, parent_pk TEXT, rejected_fields TEXT, PRIMARY KEY(table_name,row_no))")
        db.execute("CREATE INDEX rejected_parent_keys ON rejected_rows(table_name,parent_pk)")
    counts = {}
    own_errors = errors is None
    if errors is None:
        errors = Problems(inputs)
    try:
        order = [*LOAD_ORDER, *(t for t in inputs if t not in COLUMNS)]
        for table in order:
            source = inputs[table]
            if source.fingerprint != fingerprint(source.path):
                raise LoaderError(f"{table}: input changed since preflight; restart the run")
            kept = local_columns(table, schema)
            if db:
                db.execute(f"CREATE TABLE {ident(table)} (row_no INTEGER PRIMARY KEY, {', '.join(ident(c) + ' TEXT' for c in kept)})")
            count, rejected_count, pending, rejected_pending = 0, 0, [], []
            output_context = (workdir / f"{table}.load.tsv").open("w", encoding="utf-8", newline="") if write_load_files else nullcontext(None)
            with source.path.open("r", encoding="utf-8-sig", newline="") as handle, output_context as output:
                rows = reader(handle)
                if tuple(next(rows, [])) != source.header:
                    raise LoaderError(f"{table}: header changed after preflight")
                for record_no, fields in enumerate(rows, start=2):
                    count += 1
                    if len(fields) != len(source.header):
                        errors.add(table, "TSV row width", 1, [{"row": record_no, "expected": len(source.header), "found": len(fields)}])
                        rejected_count += 1
                        if track_rejections:
                            # A malformed width does not establish field alignment.
                            # Do not guess a parent key from shifted/truncated data.
                            rejected_pending.append((table, record_no, None, "TSV row width"))
                            if len(rejected_pending) >= 1000:
                                db.executemany("INSERT INTO rejected_rows VALUES (?,?,?,?)", rejected_pending)
                                rejected_pending.clear()
                        continue
                    raw = dict(zip(source.header, fields))
                    values, invalid_fields = {}, []
                    for col in columns_for(table, schema):
                        try:
                            values[col] = normalize(raw.get(col, ""), schema.columns[table][col], col)
                        except ValueError as exc:
                            errors.add(table, col, 1, [{"row": record_no, "value": raw.get(col), "error": str(exc)}])
                            invalid_fields.append(col)
                    if invalid_fields:
                        rejected_count += 1
                        if track_rejections:
                            rejected_pending.append((table, record_no, values.get(PK[table]), ", ".join(invalid_fields)))
                            if len(rejected_pending) >= 1000:
                                db.executemany("INSERT INTO rejected_rows VALUES (?,?,?,?)", rejected_pending)
                                rejected_pending.clear()
                        continue
                    if output is not None:
                        output.write("\t".join(mysql_field(values[c]) for c in columns_for(table, schema)) + "\n")
                    if db:
                        pending.append((record_no, *(values[c] for c in kept)))
                        if len(pending) >= 1000:
                            db.executemany(f"INSERT INTO {ident(table)} VALUES ({','.join(['?'] * (len(kept)+1))})", pending)
                            pending.clear()
                    if count % 250000 == 0:
                        log("INFO", f"{table}: parsed {count:,} rows")
            if db and pending:
                db.executemany(f"INSERT INTO {ident(table)} VALUES ({','.join(['?'] * (len(kept)+1))})", pending)
            if track_rejections and rejected_pending:
                db.executemany("INSERT INTO rejected_rows VALUES (?,?,?,?)", rejected_pending)
            if source.fingerprint != fingerprint(source.path):
                raise LoaderError(f"{table}: input changed while parsing; restart the run")
            counts[table] = count
            log("ERROR" if rejected_count else "OK", f"{table}: {count:,} incoming rows; accepted={count-rejected_count:,}, rejected during parsing={rejected_count:,}")
        if own_errors:
            errors.fail_if_any()
        if db:
            db.commit()
        return db, counts
    except BaseException:
        if db:
            db.close()
        raise


def comparison_columns(table: str, schema: Schema, relations: Iterable[tuple] = RELATIONS) -> list[str]:
    result = {col for key in schema.keys[table] for col in key}
    result.update(col for source, col, _, _ in relations if source == table)
    return [c for c in columns_for(table, schema) if c in result]


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


def validate_keys(db: sqlite3.Connection, schema: Schema, errors: Problems) -> dict[str, dict]:
    results = {}
    for table, keys in schema.keys.items():
        for key in keys:
            token = [ident("k_" + c) for c in key]
            nonnull = " AND ".join(c + " IS NOT NULL" for c in token)
            # MySQL UNIQUE permits repeated keys containing any NULL component.
            groups = f"SELECT {', '.join(token)} FROM {ident(table)} WHERE {nonnull} GROUP BY {', '.join(token)} HAVING COUNT(*)>1"
            join = " AND ".join(f"a.{c}=d.{c}" for c in token)
            selected = ",".join("a." + ident(c) for c in key)
            sql = f"SELECT a.row_no,{selected} FROM {ident(table)} a JOIN ({groups}) d ON {join}"
            label = f"duplicate {'primary' if key == (PK[table],) else 'unique'} key {key}"
            count = local_error(db, errors, table, label, sql)
            results[f"{table}:{','.join(key)}"] = {"duplicate_rows": count, "columns": list(key)}
    return results


def validate_recovery(db: sqlite3.Connection, errors: Problems) -> None:
    for table, count_col in [("intron_compleasm_summary", "n_introns"), ("flank_compleasm_summary", "n_flanks")]:
        cols = [count_col, "n_orthologs", "callable_bp_total"]
        condition = " OR ".join(f"{ident(c)} IS NULL OR {ident(c)}<>'0'" for c in cols)
        local_error(db, errors, table, "not_recoverable_from_compleasm requires zero feature/ortholog/callable counts",
            f"SELECT row_no,recovery_status,{','.join(ident(c) for c in cols)} FROM {ident(table)} WHERE recovery_status='not_recoverable_from_compleasm' AND ({condition})")


def validate_local(db: sqlite3.Connection, schema: Schema, errors: Problems) -> None:
    validate_keys(db, schema, errors)
    report = ValidationResult()
    for source, col, target, ref in RELATIONS:
        if target not in COLUMNS:
            continue
        relation_report(db, report, errors, (source, col, target, ref), parent_kind="incoming")
    semantic_report(db, report, errors, "intron_compleasm", "parent_id differs from referenced ortholog odb12_id",
        "FROM intron_compleasm i JOIN orthologs o ON i.k_ortholog_pk=o.k_ortholog_pk",
        "i.parent_id IS NOT NULL AND o.odb12_id IS NOT NULL", "i.parent_id<>o.odb12_id",
        "i.row_no,i.ortholog_pk,i.parent_id,o.odb12_id")
    validate_recovery(db, errors)


def existing_rows(conn: Any, schema: Schema, table: str, key: tuple[str, ...], records: list[tuple], extra_columns: tuple[str, ...] = ()) -> list[dict]:
    """Read matched incoming row numbers and requested parent fields.

    UNION batches avoid reading the entire parent/target table. Explicit casts
    prevent BIGINT keys becoming floating-point comparisons to string parameters.
    """
    if not records:
        return []
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
    selected = "DISTINCT i.row_no" + "".join(f",p.{ident(c)} AS {ident(c)}" for c in extra_columns)
    sql = f"SELECT {selected} FROM ({' UNION ALL '.join(selects)}) i JOIN {ident(table)} p ON {condition}"
    return query(conn, sql, params)


def existing_matches(conn: Any, schema: Schema, table: str, key: tuple[str, ...], records: list[tuple]) -> set[int]:
    return {int(r["row_no"]) for r in existing_rows(conn, schema, table, key, records)}


def validate_sequence_genomes(db: sqlite3.Connection, conn: Any, schema: Schema, args: argparse.Namespace, errors: Problems) -> None:
    """Check stated ortholog genomes against known SQL sequence genomes.

    Query distinct incoming pairs in batches, rather than fetching sequences.
    A missing sequence is already reported by the ordinary FK check; NULL
    genome values are permitted when the corresponding SQL column allows NULL.
    """
    cursor = db.execute(
        "SELECT MIN(row_no),sequence_pk,genome_pk,COUNT(*) FROM orthologs "
        "WHERE sequence_pk IS NOT NULL AND genome_pk IS NOT NULL "
        "GROUP BY k_sequence_pk,k_genome_pk")
    total = db.execute("SELECT COUNT(*) FROM orthologs").fetchone()[0]
    eligible, distinct_pairs, matched, comparable, mismatch_count, examples = 0, 0, 0, 0, 0, []
    while True:
        records = cursor.fetchmany(args.batch_size)
        if not records:
            break
        eligible += sum(r[3] for r in records)
        distinct_pairs += len(records)
        incoming = {int(r[0]): r for r in records}
        parents = existing_rows(conn, schema, "sequences", ("sequence_pk",),
                                [(r[0], r[1]) for r in records], ("genome_pk",))
        for parent in parents:
            row_no, sequence_pk, genome_pk, count = incoming[int(parent["row_no"])]
            sequence_genome = parent["genome_pk"]
            matched += count
            if sequence_genome is not None:
                comparable += count
            if sequence_genome is not None and str(sequence_genome) != genome_pk:
                mismatch_count += count
                if len(examples) < MAX_EXAMPLES:
                    examples.append({"row": row_no, "sequence_pk": sequence_pk,
                                     "genome_pk": genome_pk, "sequence_genome_pk": sequence_genome,
                                     "rows_with_pair": count})
    errors.add("orthologs", "genome_pk differs from the referenced sequence's genome_pk", mismatch_count, examples)
    detail = (f"rows={total:,}, non-NULL pairs={eligible:,}, distinct pairs={distinct_pairs:,}, "
              f"matched SQL sequences={matched:,}, comparable={comparable:,}, mismatches={mismatch_count:,}, "
              f"uncomparable={total-comparable:,}")
    if not comparable:
        log("INFO", "Ortholog genome assignments: no comparable rows; incoming values are NULL, SQL sequences are missing, or sequence genomes are NULL; " + detail)
    else:
        log("ERROR" if mismatch_count else "OK", "Ortholog genome assignments against SQL sequences: " + detail)


def validate_existing(db: sqlite3.Connection, conn: Any, schema: Schema, args: argparse.Namespace, errors: Problems) -> None:
    for source, col, parent, ref in RELATIONS:
        if parent in COLUMNS:
            continue
        # Group by canonical key token; query only distinct referenced values.
        cursor = db.execute(f"SELECT MIN(row_no),{ident(col)},COUNT(*) FROM {ident(source)} WHERE {ident(col)} IS NOT NULL GROUP BY {ident('k_' + col)}")
        total = db.execute(f"SELECT COUNT(*) FROM {ident(source)}").fetchone()[0]
        nonnull, distinct, missing_count, examples = 0, 0, 0, []
        while True:
            records = cursor.fetchmany(args.batch_size)
            if not records:
                break
            nonnull += sum(r[2] for r in records)
            distinct += len(records)
            found = existing_matches(conn, schema, parent, (ref,), [(r[0], r[1]) for r in records])
            for row_no, value, count in records:
                if row_no not in found:
                    missing_count += count
                    if len(examples) < MAX_EXAMPLES:
                        examples.append({"row": row_no, col: value, "rows_with_value": count})
        errors.add(source, f"missing SQL FK {col} -> {parent}.{ref}", missing_count, examples)
        log("ERROR" if missing_count else "OK" if nonnull else "INFO",
            f"{source}.{col} -> existing {parent}.{ref}: rows={total:,}, non-NULL={nonnull:,}, "
            f"distinct IDs={distinct:,}, matched={nonnull-missing_count:,}, missing={missing_count:,}, NULL skipped={total-nonnull:,}")
    validate_sequence_genomes(db, conn, schema, args, errors)
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


def full_validation(db: sqlite3.Connection, conn: Any, schema: Schema, args: argparse.Namespace, inputs: dict | None = None) -> None:
    log("INFO", "Running full Compleasm TSV relationship/duplicate validation...")
    add_comparison_tokens(db, conn, schema, args.batch_size)
    errors = Problems(inputs)
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
    genome_condition = (
        "FROM orthologs o JOIN sequences s ON o.sequence_pk=s.sequence_pk "
        "WHERE o.genome_pk IS NOT NULL AND s.genome_pk IS NOT NULL AND o.genome_pk<>s.genome_pk")
    count = query(conn, "SELECT COUNT(*) AS n " + genome_condition)[0]["n"]
    if count:
        examples = query(conn,
            "SELECT o.ortholog_pk,o.sequence_pk,o.genome_pk,s.genome_pk AS sequence_genome_pk "
            + genome_condition + f" LIMIT {MAX_EXAMPLES}")
        errors.append(f"orthologs.genome_pk: {count:,} rows disagree with referenced sequence genome_pk; examples={examples}")
    if errors:
        raise LoaderError("Post-load validation failed: " + "; ".join(errors))
    log("OK", "Post-load SQL orphan checks, intron parent-ID check, and ortholog genome-assignment check passed")


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


@dataclass
class ValidationResult:
    """Serializable, read-only validation report returned by validate_tsvs()."""
    passed: bool = False
    counts: dict[str, int] = field(default_factory=dict)
    valid_rows: dict[str, int] = field(default_factory=dict)
    rejected_rows: dict[str, int] = field(default_factory=dict)
    files: dict[str, str] = field(default_factory=dict)
    issues: list[dict] = field(default_factory=list)
    keys: dict[str, dict] = field(default_factory=dict)
    relationships: dict[str, dict] = field(default_factory=dict)
    semantics: dict[str, dict] = field(default_factory=dict)
    mappings: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return asdict(self)


def inspect_starters(root: Path, check_nhphyml: bool) -> dict[str, InputFile]:
    if not root.is_dir() or not os.access(root, os.R_OK | os.X_OK):
        raise LoaderError(f"Starter TSV directory does not exist/is not readable: {root}")
    tables = ["genomes", "sequences", "natural_history"]
    if check_nhphyml:
        tables.append("species_name_audit")
    inputs = {}
    for table in tables:
        path = root / f"{table}.tsv"
        if not path.is_file() or not os.access(path, os.R_OK):
            raise LoaderError(f"{table}: required readable starter TSV missing: {path}")
        before = fingerprint(path)
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            header = next(reader(handle), [])
        if len(header) != len(set(header)):
            raise LoaderError(f"{path}: duplicate TSV header names")
        missing = set(STARTER_COLUMNS[table]) - set(header)
        extra = set(header) - set(STARTER_COLUMNS[table])
        # Documented legacy extension: all canonical audit fields plus genome_pk.
        if table == "species_name_audit":
            extra.discard("genome_pk")
        if missing or extra:
            raise LoaderError(f"{path}: TSV/schema header mismatch; missing={sorted(missing)}, unexpected={sorted(extra)}")
        if before != fingerprint(path):
            raise LoaderError(f"Input changed while checking header: {path}")
        inputs[table] = InputFile(path, tuple(header), before)
    return inputs


def offline_schema(inputs: dict[str, InputFile]) -> Schema:
    """Expected v9 types/nullability, without inventing live SQL string limits."""
    columns, keys = {}, {}
    for table, source in inputs.items():
        names = COLUMNS[table] if table in COLUMNS else STARTER_COLUMNS[table]
        names = list(names)
        if table == "species_name_audit" and "genome_pk" in source.header:
            names.append("genome_pk")
        columns[table] = {}
        for col in names:
            kind = expected_type(col)
            data_type = {"string": "varchar", "bool": "tinyint"}.get(kind, kind)
            columns[table][col] = {
                "DATA_TYPE": data_type, "COLUMN_TYPE": data_type,
                "IS_NULLABLE": "NO" if col == PK[table] else "YES",
                "CHARACTER_MAXIMUM_LENGTH": None, "COLLATION_NAME": None,
                "DATETIME_PRECISION": 0,
            }
        keys[table] = [(PK[table],), *UNIQUE_KEYS[table]] if table in COLUMNS else list(STARTER_KEYS[table])
    return Schema(columns, keys, {})


def offline_relations(inputs: dict) -> list[tuple]:
    relations = [*RELATIONS, *STARTER_RELATIONS]
    if "species_name_audit" in inputs:
        relations.append(("species_name_audit", "species_pk", "natural_history", "species_pk"))
        if "genome_pk" in inputs["species_name_audit"].header:
            relations.append(("species_name_audit", "genome_pk", "genomes", "genome_pk"))
    return relations


def add_offline_tokens(db: sqlite3.Connection, schema: Schema, relations: list[tuple]) -> None:
    # SQLite BINARY comparison preserves exact text; integer IDs are already
    # normalized decimal strings, including BIGINTs above floating precision.
    for table in schema.keys:
        for col in comparison_columns(table, schema, relations):
            token = ident("k_" + col)
            db.execute(f"ALTER TABLE {ident(table)} ADD COLUMN {token} TEXT")
            db.execute(f"UPDATE {ident(table)} SET {token}={ident(col)}")
            db.execute(f"CREATE INDEX {ident('idx_' + table + '_' + col)} ON {ident(table)} ({token})")
        for i, key in enumerate(schema.keys[table]):
            if len(key) > 1:
                db.execute(f"CREATE INDEX {ident('key_' + table + '_' + str(i))} ON {ident(table)} ({','.join(ident('k_'+c) for c in key)})")
        log("INFO", f"Indexed {table} for offline validation")
    db.commit()


def relation_report(db: sqlite3.Connection, result: ValidationResult, errors: Problems, relation: tuple, *, parent_kind: str = "TSV") -> None:
    """Count child references by parent availability; rejection never satisfies FK.

    missing_rows is the total unavailable-reference count for API compatibility.
    rejected_parent_rows, absent_parent_rows and unresolved_parent_rows partition
    those child rows. parent_rows_rejected counts actual rejected parent records.
    SQL-mode staging aborts on parsing errors before relationship checks, so its
    incoming-parent reports retain the existing lookup path.
    """
    source, col, parent, ref = relation
    name = f"{source}.{col} -> {parent}.{ref}"
    total, nulls, distinct = db.execute(
        f"SELECT COUNT(*),COALESCE(SUM({ident(col)} IS NULL),0),COUNT(DISTINCT {ident('k_'+col)}) FROM {ident(source)}").fetchone()
    selected = f"c.row_no,c.{ident(PK[source])},c.{ident(col)}"
    unavailable = (f"c.{ident(col)} IS NOT NULL AND NOT EXISTS "
                   f"(SELECT 1 FROM {ident(parent)} p WHERE p.{ident('k_'+ref)}=c.{ident('k_'+col)})")
    sql = f"SELECT {selected} FROM {ident(source)} c WHERE {unavailable}"
    tracks_rejections = db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='rejected_rows'").fetchone() is not None
    rejected, absent, unresolved, parent_rejected, unidentified = 0, 0, 0, 0, 0
    if tracks_rejections:
        parent_rejected, unidentified = db.execute(
            "SELECT COUNT(*),COALESCE(SUM(parent_pk IS NULL),0) FROM rejected_rows WHERE table_name=?", (parent,)).fetchone()
        rejected_sql = (
            f"SELECT {selected},r.row_no AS parent_record_number,r.parent_pk,"
            f"r.rejected_fields AS parent_rejected_fields,? AS parent_file FROM {ident(source)} c "
            f"JOIN rejected_rows r ON r.table_name=? AND r.parent_pk=c.{ident('k_'+col)} "
            "AND r.row_no=(SELECT MIN(r2.row_no) FROM rejected_rows r2 "
            "WHERE r2.table_name=r.table_name AND r2.parent_pk=r.parent_pk) "
            f"WHERE {unavailable}")
        rejected = local_error(db, errors, source,
            f"parent exists in source TSV but was rejected during parsing: {col} -> {parent}.{ref}",
            rejected_sql, (result.files.get(parent), parent))
        absent_sql = sql + (" AND NOT EXISTS (SELECT 1 FROM rejected_rows r WHERE r.table_name=? "
                            f"AND r.parent_pk=c.{ident('k_'+col)})")
        if unidentified:
            # With unknown keys in malformed rows, do not claim source absence.
            unresolved = local_error(db, errors, source,
                f"unresolved TSV FK {col} -> {parent}.{ref}: absent from accepted/identified rows; "
                f"{unidentified} rejected parent rows have no identifiable key, so source absence is inconclusive",
                absent_sql, (parent,))
        else:
            absent = local_error(db, errors, source,
                f"missing TSV FK {col} -> {parent}.{ref}: ID genuinely absent from source TSV",
                absent_sql, (parent,))
        missing = rejected + absent + unresolved
    else:
        missing = local_error(db, errors, source, f"missing {parent_kind} FK {col} -> {parent}.{ref}", sql)
        absent = missing
    nonnull = total - nulls
    stats = {"rows_examined": total, "non_null_rows": nonnull, "distinct_ids": distinct,
             "matched_rows": nonnull - missing, "missing_rows": missing, "unavailable_rows": missing, "null_rows_skipped": nulls,
             "rejected_parent_rows": rejected, "absent_parent_rows": absent, "unresolved_parent_rows": unresolved,
             "parent_rows_rejected": parent_rejected, "unidentified_rejected_parent_rows": unidentified,
             "source_rows_read": result.counts.get(source, total),
             "source_rows_rejected": result.rejected_rows.get(source, 0)}
    result.relationships[name] = stats
    level = "ERROR" if missing else "OK" if nonnull else "INFO"
    classification = (f", references to rejected parents={rejected:,}, IDs absent from source={absent:,}, "
                      f"absence inconclusive={unresolved:,}" if tracks_rejections else "")
    unavailable_label = "unavailable" if tracks_rejections else "missing"
    log(level, f"{name}: rows={total:,}, non-NULL={nonnull:,}, distinct={distinct:,}, matched={nonnull-missing:,}, {unavailable_label}={missing:,}{classification} NULL skipped={nulls:,}"
        + ("; no non-NULL references to compare" if not nonnull else ""))


def semantic_report(db: sqlite3.Connection, result: ValidationResult, errors: Problems,
                    table: str, label: str, join: str, known: str, unequal: str,
                    selected: str) -> None:
    total = db.execute(f"SELECT COUNT(*) FROM {ident(table)}").fetchone()[0]
    comparable = db.execute(f"SELECT COUNT(DISTINCT {selected.split(',')[0]}) {join} WHERE {known}").fetchone()[0]
    count = local_error(db, errors, table, label,
                        f"SELECT {selected} {join} WHERE ({known}) AND ({unequal})")
    result.semantics[f"{table}:{label}"] = {
        "rows_examined": total, "comparable_rows": comparable,
        "mismatched_rows": count, "uncomparable_rows": max(0, total-comparable)}
    if comparable:
        log("ERROR" if count else "OK", f"{table}: {label}; comparable={comparable:,}, mismatches={count:,}, uncomparable={max(0,total-comparable):,}")
    else:
        log("INFO", f"{table}: {label}; no comparable rows (NULL values or unmatched parents); no matching-success claim")


def offline_semantics(db: sqlite3.Connection, result: ValidationResult, errors: Problems) -> None:
    semantic_report(db, result, errors, "orthologs", "genome_pk differs from the referenced sequence's genome_pk",
        "FROM orthologs o LEFT JOIN sequences s ON o.k_sequence_pk=s.k_sequence_pk",
        "o.genome_pk IS NOT NULL AND s.genome_pk IS NOT NULL", "o.genome_pk<>s.genome_pk",
        "o.row_no,o.ortholog_pk,o.sequence_pk,o.genome_pk,s.genome_pk AS sequence_genome_pk")
    semantic_report(db, result, errors, "intron_compleasm", "parent_id differs from referenced ortholog odb12_id",
        "FROM intron_compleasm i LEFT JOIN orthologs o ON i.k_ortholog_pk=o.k_ortholog_pk",
        "i.parent_id IS NOT NULL AND o.odb12_id IS NOT NULL", "i.parent_id<>o.odb12_id",
        "i.row_no,i.intron_pk,i.ortholog_pk,i.parent_id,o.odb12_id")
    for table in ("ortholog_summary", "intron_compleasm_summary", "flank_compleasm_summary"):
        semantic_report(db, result, errors, table, "species_pk differs from the referenced genome's species_pk",
            f"FROM {ident(table)} c LEFT JOIN genomes g ON c.k_genome_pk=g.k_genome_pk",
            "c.species_pk IS NOT NULL AND g.species_pk IS NOT NULL", "c.species_pk<>g.species_pk",
            f"c.row_no,c.{ident(PK[table])},c.genome_pk,c.species_pk,g.species_pk AS genome_species_pk")


def nh_error(db: sqlite3.Connection, errors: Problems, table: str, label: str, sql: str) -> int:
    return local_error(db, errors, table, "nhPhyML readiness: " + label, sql)


def duplicate_mapping(db: sqlite3.Connection, errors: Problems, table: str,
                      name: str, mapping: str, key: tuple[str, ...]) -> int:
    columns = ",".join(ident(c) for c in key)
    join = " AND ".join(f"m.{ident(c)} IS d.{ident(c)}" for c in key)
    return nh_error(db, errors, table, name,
        f"SELECT m.* FROM {ident(mapping)} m JOIN "
        f"(SELECT {columns} FROM {ident(mapping)} GROUP BY {columns} HAVING COUNT(*)>1) d ON {join}")


def check_nhphyml_reference_mappings(db: sqlite3.Connection, inputs: dict, result: ValidationResult, errors: Problems) -> None:
    """Mirror the reference routes in 06_ v3, with explicit data-loss checks."""
    log("INFO", "Checking nhPhyML reference mappings (not run inputs or result files)...")
    db.create_function("nh_strip", 1, lambda v: None if v is None else str(v).strip(), deterministic=True)
    db.create_function("nh_species", 1, lambda v: None if v is None else str(v).strip().lower(), deterministic=True)
    nh_error(db, errors, "sequences", "06_ requires non-NULL sequence genome_pk",
             "SELECT row_no,sequence_pk,genome_pk FROM sequences WHERE genome_pk IS NULL")
    db.execute(
        "CREATE TABLE nh_ortholog_map AS SELECT o.row_no,o.ortholog_pk,o.sequence_pk,"
        "COALESCE(o.genome_pk,s.genome_pk) AS genome_pk,nh_strip(o.odb12_id) AS odb12_id "
        "FROM orthologs o LEFT JOIN sequences s ON o.k_sequence_pk=s.k_sequence_pk")
    db.execute("CREATE INDEX nh_ortholog_keys ON nh_ortholog_map(genome_pk,odb12_id)")
    unresolved = nh_error(db, errors, "orthologs", "no resolvable genome (06_ would drop these rows)",
        "SELECT row_no,ortholog_pk,sequence_pk,odb12_id FROM nh_ortholog_map WHERE genome_pk IS NULL")
    nh_error(db, errors, "orthologs", "missing ODB12 ID for genome+ODB12 mapping",
        "SELECT row_no,ortholog_pk,genome_pk,odb12_id FROM nh_ortholog_map WHERE odb12_id IS NULL OR odb12_id=''")
    db.execute("CREATE VIEW nh_usable_orthologs AS SELECT * FROM nh_ortholog_map WHERE genome_pk IS NOT NULL AND odb12_id IS NOT NULL AND odb12_id<>''")
    duplicate_count = duplicate_mapping(db, errors, "orthologs",
        "ambiguous genome_pk+odb12_id -> ortholog_pk; duplicated copies need an explicit analysis selection strategy",
        "nh_usable_orthologs", ("genome_pk", "odb12_id"))
    if duplicate_count:
        log("ERROR", "nhPhyML genome+ODB12 mapping has multiple copies. Biologically valid duplicated orthologs (including triploid I. braminus) are retained; no copy was selected or removed.")

    # 06_ checks accession uniqueness again after stripping whitespace, including
    # duplicate NULLs (pandas duplicated differs from SQL UNIQUE in this respect).
    db.execute("CREATE TABLE nh_accessions AS SELECT row_no,genome_pk,nh_strip(accession_id) AS accession_id FROM genomes")
    duplicate_mapping(db, errors, "genomes", "ambiguous stripped accession_id -> genome_pk", "nh_accessions", ("accession_id",))
    db.execute("CREATE TABLE nh_focal_genomes AS SELECT DISTINCT genome_pk FROM nh_ortholog_map WHERE genome_pk IS NOT NULL")
    db.execute("CREATE UNIQUE INDEX nh_focal_idx ON nh_focal_genomes(genome_pk)")
    legacy_direct = "genome_pk" in inputs["species_name_audit"].header
    target = "genome_pk" if legacy_direct else "species_pk"
    branch = "legacy direct audit.genome_pk" if legacy_direct else "audit.species_pk -> genomes.species_pk"
    log("INFO", "nhPhyML species mapping route: " + branch)
    db.execute(f"CREATE TABLE nh_audit AS SELECT row_no,species_name_audit_pk,nh_species(species_normalized) AS species_name,{ident(target)} AS target_pk FROM species_name_audit")
    usable, excluded = db.execute("SELECT COALESCE(SUM(species_name IS NOT NULL AND target_pk IS NOT NULL),0),COALESCE(SUM(species_name IS NULL OR target_pk IS NULL),0) FROM nh_audit").fetchone()
    log("INFO", f"Audit rows usable for the selected route={usable:,}, NULL name/ID excluded from that route={excluded:,}")
    # Identical observations are collapsed only in this temporary lookup, never
    # in source files. Conflicting mappings are reported with original row IDs.
    nh_error(db, errors, "species_name_audit", "normalized species name has conflicting target IDs",
        "SELECT a.* FROM nh_audit a JOIN "
        "(SELECT species_name FROM nh_audit WHERE species_name IS NOT NULL AND target_pk IS NOT NULL GROUP BY species_name HAVING COUNT(DISTINCT target_pk)>1) d "
        "ON a.species_name=d.species_name WHERE a.target_pk IS NOT NULL")
    db.execute("CREATE TABLE nh_audit_names AS SELECT species_name,target_pk,MIN(row_no) AS row_no FROM nh_audit WHERE species_name IS NOT NULL AND target_pk IS NOT NULL GROUP BY species_name,target_pk")
    if legacy_direct:
        db.execute("CREATE TABLE nh_species_map AS SELECT a.row_no,a.species_name,g.genome_pk FROM nh_audit_names a JOIN genomes g ON a.target_pk=g.k_genome_pk")
    else:
        db.execute(
            "CREATE TABLE nh_selected_genomes AS SELECT g.* FROM genomes g JOIN "
            "(SELECT species_pk,SUM(CASE WHEN is_current='1' THEN 1 ELSE 0 END) AS n_current FROM genomes WHERE species_pk IS NOT NULL GROUP BY species_pk) s "
            "ON g.species_pk=s.species_pk WHERE s.n_current<>1 OR g.is_current='1'")
        db.execute("CREATE INDEX nh_selected_species ON nh_selected_genomes(species_pk)")
        nh_error(db, errors, "genomes", "multiple genome candidates remain after the 06_ current-genome rule",
            "SELECT g.row_no,g.genome_pk,g.species_pk,g.accession_id,g.is_current FROM nh_selected_genomes g JOIN "
            "(SELECT species_pk FROM nh_selected_genomes GROUP BY species_pk HAVING COUNT(*)>1) d ON g.species_pk=d.species_pk")
        db.execute("CREATE TABLE nh_species_map AS SELECT a.row_no,a.species_name,g.genome_pk FROM nh_audit_names a JOIN nh_selected_genomes g ON a.target_pk=g.species_pk")
    db.execute("CREATE INDEX nh_species_genomes ON nh_species_map(genome_pk)")
    duplicate_mapping(db, errors, "species_name_audit", "resolved species name does not identify one genome", "nh_species_map", ("species_name",))
    omitted = nh_error(db, errors, "genomes", "focal genome is absent from the resolved audit species mapping (missing audit or current-genome selection)",
        "SELECT g.row_no,g.genome_pk,g.species_pk,g.accession_id FROM genomes g JOIN nh_focal_genomes f ON f.genome_pk=g.k_genome_pk "
        "WHERE NOT EXISTS (SELECT 1 FROM nh_species_map m WHERE m.genome_pk=g.genome_pk)")
    unmatched_names = db.execute("SELECT COUNT(*) FROM nh_audit_names a WHERE NOT EXISTS (SELECT 1 FROM nh_species_map m WHERE m.species_name=a.species_name)").fetchone()[0]
    focal = db.execute("SELECT COUNT(*) FROM nh_focal_genomes").fetchone()[0]
    mapped_names = db.execute("SELECT COUNT(*) FROM nh_species_map").fetchone()[0]
    mapped_genomes = db.execute("SELECT COUNT(DISTINCT genome_pk) FROM nh_species_map").fetchone()[0]
    resolved_rows = db.execute("SELECT COUNT(*) FROM nh_ortholog_map WHERE genome_pk IS NOT NULL").fetchone()[0]
    result.mappings = {
        "requested": True, "species_route": branch, "unresolved_ortholog_rows": unresolved,
        "ambiguous_ortholog_rows": duplicate_count, "resolved_ortholog_rows": resolved_rows,
        "audit_rows_usable": usable, "audit_null_rows_excluded": excluded,
        "unmatched_audit_names_outside_mapping": unmatched_names,
        "resolved_species_names": mapped_names, "mapped_genomes": mapped_genomes,
        "focal_genomes": focal, "omitted_focal_genomes": omitted,
    }
    log("INFO", f"nhPhyML references: resolved ortholog rows={resolved_rows:,}, focal genomes={focal:,}, resolved species names={mapped_names:,}, audit names without a genome={unmatched_names:,}")
    if unmatched_names:
        result.warnings.append(f"{unmatched_names} audit species names have no genome mapping; broad non-focal trait species are not automatically errors.")
    if not focal:
        result.warnings.append("No focal ortholog genomes exist; nhPhyML checks had no focal-genome coverage to establish.")
        result.skipped.append("nhPhyML focal-genome coverage: no focal genomes")


def collect_issues(errors: Problems, files: dict[str, str]) -> list[dict]:
    issues = []
    for (table, label), item in errors.items.items():
        examples = []
        for example in item["examples"]:
            if isinstance(example, dict):
                example = dict(example)
                if table in files:
                    example.setdefault("file", files[table])
                if "row_no" in example:
                    example.setdefault("record_number", example["row_no"])
                elif "row" in example:
                    example.setdefault("record_number", example["row"])
            examples.append(example)
        issues.append({"table": table, "file": files.get(table), "check": label,
                       "problematic_rows": item["count"], "examples": examples})
    return issues


def validate_tsvs(input_dir: str | Path, starter_tsv_dir: str | Path, *,
                  check_nhphyml_mappings: bool = False, temp_dir: str | Path | None = None,
                  verbose: bool = True) -> ValidationResult:
    """Read-only offline API; return failures rather than loading or raising.

    Every canonical starter column is required. The supported legacy audit is
    the canonical audit header plus genome_pk; arbitrary incomplete headers are
    rejected. No live SQL collation/column-length check is possible offline.
    Errors, including malformed inputs, are returned in result.issues. An
    interrupt still raises KeyboardInterrupt. Temporary indexes are cleaned up.
    """
    result = ValidationResult()
    errors = Problems()
    db = None
    context = nullcontext() if verbose else redirect_stdout(io.StringIO())
    with context, expanded_csv_fields():
        try:
            log("INFO", "Running TSV-only validation; MySQL is not used")
            log("OK", f"Python executable: {sys.executable} ({sys.version.split()[0]})")
            inputs = inspect_inputs(Path(input_dir).expanduser().resolve())
            inputs.update(inspect_starters(Path(starter_tsv_dir).expanduser().resolve(), check_nhphyml_mappings))
            result.files = {table: str(source.path) for table, source in inputs.items()}
            for table, path in result.files.items():
                log("INFO", f"{table}: {path}")
            result.warnings.append("Offline text comparisons are exact; live SQL collations, string lengths, existing SQL records and upload safety must still be checked in MySQL mode.")
            if "recovery_status" not in inputs["flank_compleasm_summary"].header:
                result.warnings.append("Legacy flank_summary.tsv has no recovery fields; its recovery status cannot be validated or inferred.")
                result.skipped.append("Legacy flank summary recovery status: fields absent")
            if not check_nhphyml_mappings:
                result.mappings = {"requested": False}
                result.skipped.append("nhPhyML reference readiness: not requested")
            schema = offline_schema(inputs)
            with tempfile.TemporaryDirectory(prefix="compleasm_tsv_v5_", dir=temp_dir) as temp:
                try:
                    db, result.counts = stage_files(inputs, schema, Path(temp), True,
                                                   write_load_files=False, errors=errors)
                    if db is None:
                        raise LoaderError("Internal error: missing offline validation database")
                    result.valid_rows = {t: db.execute(f"SELECT COUNT(*) FROM {ident(t)}").fetchone()[0] for t in inputs}
                    result.rejected_rows = {t: result.counts[t]-result.valid_rows[t] for t in inputs}
                    for table, count in result.counts.items():
                        if count == 0:
                            result.warnings.append(f"{table}: input contains zero data rows; references cannot match an empty parent.")
                    relations = offline_relations(inputs)
                    add_offline_tokens(db, schema, relations)
                    result.keys = validate_keys(db, schema, errors)
                    for name, stats in result.keys.items():
                        log("ERROR" if stats["duplicate_rows"] else "OK", f"{name}: duplicate rows={stats['duplicate_rows']:,}")
                    for relation in relations:
                        relation_report(db, result, errors, relation)
                    offline_semantics(db, result, errors)
                    validate_recovery(db, errors)
                    if check_nhphyml_mappings:
                        check_nhphyml_reference_mappings(db, inputs, result, errors)
                    # Detect parent/child edits occurring after their initial scan.
                    for table, source in inputs.items():
                        if source.fingerprint != fingerprint(source.path):
                            errors.add(table, "input changed during validation; restart with stable files", 1, [{"file": str(source.path)}])
                finally:
                    if db is not None:
                        db.close()
                        db = None
        except (LoaderError, OSError, csv.Error, sqlite3.Error, ValueError) as exc:
            errors.add("inputs", "offline validation could not finish", 1, [{"error": str(exc)}])
        result.issues = collect_issues(errors, result.files)
        result.passed = not result.issues
        for warning in result.warnings:
            log("WARN", warning)
        for issue in result.issues:
            log("ERROR", f"{issue['table']}: {issue['check']}; problematic rows={issue['problematic_rows']:,}; examples={issue['examples']}")
        if result.passed:
            log("OK", "TSV relationship validation passed. No source TSV or SQL data changed.")
            if check_nhphyml_mappings:
                if result.mappings.get("focal_genomes"):
                    log("OK", "Reference mappings passed the requested nhPhyML readiness checks (reference scope only).")
                else:
                    log("INFO", "No nhPhyML focal-genome comparisons were possible; non-focal reference checks found no errors.")
        else:
            log("ERROR", "TSV validation failed. No source TSV or SQL data changed.")
    return result


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version=VERSION)
    p.add_argument("--input-dir", required=True, help="Directory with all eight logical Compleasm TSVs")
    p.add_argument("--mysql-db", "--database", dest="mysql_db", help="Required in MySQL modes; unused in TSV-only mode")
    p.add_argument("--starter-tsv-dir", help="Required TSV-only parent directory containing genomes, sequences, and natural_history")
    p.add_argument("--check-nhphyml-mappings", action="store_true", help="Additional TSV-only reference mapping checks for 06_; requires species_name_audit.tsv")
    p.add_argument("--mysql-user", "--user", default="root", dest="mysql_user")
    p.add_argument("--mysql-password", "--password", default=os.environ.get("MYSQL_PWD"), dest="mysql_password", help="Otherwise use MYSQL_PWD or interactive prompt")
    p.add_argument("--mysql-host", "--host", default="localhost", dest="mysql_host")
    p.add_argument("--mysql-port", "--port", type=int, default=3306, dest="mysql_port")
    p.add_argument("--connect-timeout", type=int, default=10)
    p.add_argument("--mode", choices=("insert", "ignore", "replace"), default="insert")
    modes = p.add_mutually_exclusive_group()
    modes.add_argument("--tsv-check-only", action="store_true", help="Full offline TSV validation against starter TSVs; no MySQL imports, connections, or writes")
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
    if args.tsv_check_only:
        if not args.starter_tsv_dir:
            p.error("--starter-tsv-dir is required with --tsv-check-only")
        if args.skip_data_validation or args.truncate_compleasm_tables or args.no_post_check or args.mode != "insert":
            p.error("--tsv-check-only cannot use --skip-data-validation, --truncate-compleasm-tables, --no-post-check, or --mode ignore/replace")
    else:
        if not args.mysql_db:
            p.error("--mysql-db is required for MySQL modes")
        if args.check_nhphyml_mappings or args.starter_tsv_dir:
            p.error("--check-nhphyml-mappings and --starter-tsv-dir require --tsv-check-only")
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
    if args.tsv_check_only:
        result = validate_tsvs(args.input_dir, args.starter_tsv_dir,
                               check_nhphyml_mappings=args.check_nhphyml_mappings,
                               temp_dir=args.temp_dir, verbose=True)
        if not result.passed:
            raise LoaderError(f"Offline TSV validation failed in {len(result.issues)} category/categories. No loading performed.")
        return
    conn, db = None, None
    with tempfile.TemporaryDirectory(prefix="compleasm_sql_v5_", dir=args.temp_dir) as temp:
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
                full_validation(db, conn, schema, args, inputs)
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
