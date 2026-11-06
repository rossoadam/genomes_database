#!/usr/bin/env python3
"""
06_nhphyml_to_tsvs_v6.py

Build two SQL-ready TSV files from nhPhyML/GC4 outputs:

    nhphyml_ortholog_results.tsv
    nhphyml_species_summary.tsv

The script uses pandas and follows the nhPhyML-related tables in
`gc3_dynamics_v13_reduced.txt`.

Version v6 updates
------------------
- Updated the authoritative compatibility target to `gc3_dynamics_v13_reduced.txt`.
- `orthologs.genome_pk` and `orthologs.odb12_id` are now required exactly as in
  v13; NULL/blank values are fatal schema violations.
- `orthologs.sequence_pk` remains nullable. When populated, it is used only to
  verify consistency with `orthologs.genome_pk`; it is never used to reconstruct
  or rescue a missing genome identity.
- Populated `sequence_pk` values must exist in sequences.tsv and must resolve to
  the same genome_pk carried directly by orthologs.tsv.
- The v13 `(genome_pk, odb12_id)` unique key is enforced before nhPhyML run
  processing. The nullable `(sequence_pk, odb12_id)` unique key is also checked
  using SQL-like NULL semantics.
- Added lightweight v13 ortholog-contract self-checks while preserving the SQL
  `\\N` self-check, v5 alignment-root inference, run mappings, GC3/GC4 joins,
  output schemas, and deterministic primary-key assignment.

Version v5 updates
------------------
- Alignment roots are no longer hard-coded to unrelated absolute paths.
- When --ortholog-alignments-root and/or --summary-alignments-root are omitted,
  both are derived from --sql-tsv-dir as:

      <sql-tsv-dir parent>/compleasm/alignments

  Thus, for:

      --sql-tsv-dir /Users/rossoaa/projects/genomes/records/sql_tsvs

  the inferred alignment root is:

      /Users/rossoaa/projects/genomes/records/compleasm/alignments

- Explicit --ortholog-alignments-root and --summary-alignments-root overrides
  remain supported.

Version v4 updates
------------------
- Historical v4 change: compatibility target was updated to the then-current gc3_dynamics_v10_reduced SQL schema.
- Canonical nhphyml_runs.tsv metadata now includes the schema's nullable
  `created_at` column. Legacy six-column nhphyml_runs.tsv files remain accepted
  with a warning so existing manually curated metadata is not broken.
- nhphyml_runs.tsv `n_species` and `n_genes` are validated as nullable integers,
  matching the SQL schema.
- Missing values in generated SQL-ready output TSVs are written as the SQL NULL
  marker `\\N` rather than blank fields.
- Existing v3 nhPhyML mapping behavior remains aligned with
  05_to_sql_compleasm_features_v5.py's `--check-nhphyml-mappings` readiness
  checks.

Version v3 updates
------------------
- SQL `\\N` markers are now parsed as NULL during common TSV ingestion.
- Nullable numeric and foreign-key fields accept `\\N` as a legitimate missing
  value.
- Required numeric/integer fields still reject NULL values, including `\\N`.
- Invalid non-numeric strings still fail validation with a clear error.
- SQL-null handling is centralized in `read_tsv()` and reinforced by shared
  numeric normalization/validation helpers.

Version v2 historical compatibility notes
------------------------------------------
Earlier schema versions permitted fallback from `orthologs.sequence_pk` to
`sequences.genome_pk` when `orthologs.genome_pk` was absent. That behavior is
historical only and is intentionally disabled in v6: v13 requires
`orthologs.genome_pk` directly.

The legacy species_name_audit compatibility route from v1 is retained when a
non-schema `genome_pk` column is present, but the canonical v13 route is
`species_name_audit.species_pk -> genomes.species_pk -> genomes.genome_pk`.

No ancestral/internal nhPhyML node is retained unless it resolves to a genome.
"""

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Iterable

import pandas as pd


# -----------------------------------------------------------------------------
# Manual metadata files
# -----------------------------------------------------------------------------

TREE_FILES_COLUMNS = [
    "tree_pk",
    "nhphyml_run_pk",
    "tree_type",
    "tree_text",
    "branch_length_type",
    "notes",
]

NHPHYML_RUNS_COLUMNS = [
    "nhphyml_run_pk",
    "run_name",
    "alignment_name",
    "n_species",
    "n_genes",
    "notes",
    "created_at",
]

NHPHYML_RUNS_LEGACY_COLUMNS = [
    "nhphyml_run_pk",
    "run_name",
    "alignment_name",
    "n_species",
    "n_genes",
    "notes",
]


# -----------------------------------------------------------------------------
# Expected output columns
# -----------------------------------------------------------------------------

ORTHOLOG_OUTPUT_COLUMNS = [
    "nhphyml_ortholog_result_pk",
    "nhphyml_run_pk",
    "ortholog_pk",
    "genome_pk",
    "gc3_nhphyml",
    "gc4",
    "status",
    "created_at",
]

SPECIES_SUMMARY_OUTPUT_COLUMNS = [
    "nhphyml_species_summary_pk",
    "nhphyml_run_pk",
    "genome_pk",
    "mean_gc3",
    "mean_gc4",
    "sd_gc3",
    "sd_gc4",
    "n_orthologs",
    "created_at",
]


# -----------------------------------------------------------------------------
# Run mapping
# -----------------------------------------------------------------------------

RUN_DIRS = {
    1: "t1.0_e1_o5_with_i_braminus",
    2: "t1.0_e1_o15_with_i_braminus",
    3: "t1.0_e1_o25_with_i_braminus",
    4: "t1.0_e1_o35_with_i_braminus",
}


# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------


# SQL-oriented TSVs use the MySQL/LOAD DATA convention `\N` for NULL.
# Match the whole field only; strings that merely contain the letter N are not
# treated as missing.
SQL_NULL_MARKERS = {"\\N"}


def print_manual_file_warning() -> None:
    """Print the required manual-file warning before doing any other work."""
    line = "=" * 79
    print(line, file=sys.stderr)
    print("WARNING: MANUAL nhPhyML METADATA FILES ARE REQUIRED", file=sys.stderr)
    print(line, file=sys.stderr)
    print(
        "You must manually create tree_files.tsv and nhphyml_runs.tsv.\n"
        "The script does NOT generate these two metadata files.\n",
        file=sys.stderr,
    )
    print("tree_files.tsv fields:", file=sys.stderr)
    print("  " + "\t".join(TREE_FILES_COLUMNS), file=sys.stderr)
    print("nhphyml_runs.tsv fields (canonical v13 schema):", file=sys.stderr)
    print("  " + "\t".join(NHPHYML_RUNS_COLUMNS), file=sys.stderr)
    print(r"  created_at is nullable; use \N if you do not want to assign it manually.", file=sys.stderr)
    print("  Legacy six-column nhphyml_runs.tsv files from v1-v3 remain accepted with a warning.", file=sys.stderr)
    print(line, file=sys.stderr)
    print(file=sys.stderr)


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")


def require_columns(df: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(
            f"{label} is missing required column(s): {', '.join(missing)}\n"
            f"Observed columns: {', '.join(map(str, df.columns))}"
        )


def read_tsv(path: Path, label: str) -> pd.DataFrame:
    """
    Read a SQL-oriented TSV.

    The MySQL/LOAD DATA NULL marker `\\N` is treated as a missing value for every
    TSV read through this common ingestion function. Normal pandas NA handling
    remains enabled.
    """
    require_file(path, label)
    try:
        return pd.read_csv(
            path,
            sep="\t",
            dtype=str,
            keep_default_na=True,
            na_values=list(SQL_NULL_MARKERS),
        )
    except Exception as exc:
        raise RuntimeError(f"Could not read {label} as a TSV: {path}\n{exc}") from exc


def format_mtime(path: Path) -> str:
    """Return file modification time in SQL datetime format, using local time."""
    return datetime.fromtimestamp(os.path.getmtime(path)).strftime("%Y-%m-%d %H:%M:%S")


def normalize_sql_nulls(series: pd.Series) -> pd.Series:
    """
    Normalize exact SQL NULL markers to pandas missing values.

    This is intentionally conservative: only whole-field markers in
    SQL_NULL_MARKERS are converted. Ordinary strings containing the letter N are
    left unchanged. Surrounding whitespace is ignored when recognizing a marker.
    """
    normalized = series.astype("string")
    stripped = normalized.str.strip()
    return normalized.mask(stripped.isin(SQL_NULL_MARKERS), pd.NA)


def normalize_text(series: pd.Series) -> pd.Series:
    """Trim and lowercase identifiers used for case-insensitive joins."""
    return series.astype("string").str.strip().str.lower()


def to_numeric_checked(series: pd.Series, label: str) -> pd.Series:
    """
    Convert a column to numeric and fail on non-empty, non-numeric values.

    Legitimate missing values include pandas NA/NaN, blank fields, and SQL `\\N`.
    Invalid non-empty strings such as "abc", "species_17", or "17x" remain
    validation errors.
    """
    normalized = normalize_sql_nulls(series)
    stripped = normalized.str.strip()
    numeric = pd.to_numeric(normalized, errors="coerce")

    bad = numeric.isna() & normalized.notna() & stripped.ne("")
    if bad.any():
        examples = normalized.loc[bad].drop_duplicates().head(10).tolist()
        raise ValueError(
            f"Non-numeric value(s) found in {label}: {examples}"
            + (" ..." if bad.sum() > 10 else "")
        )
    return numeric


def to_integer_checked(series: pd.Series, label: str, nullable: bool = False) -> pd.Series:
    """
    Convert a column to pandas Int64, requiring integer-valued input.

    When nullable=True, pandas missing values, blank fields, and SQL `\\N` are
    permitted. When nullable=False, any missing value remains an error.
    """
    numeric = to_numeric_checked(series, label)
    non_integer = numeric.notna() & ((numeric % 1) != 0)
    if non_integer.any():
        examples = series.loc[non_integer].drop_duplicates().head(10).tolist()
        raise ValueError(f"Non-integer value(s) found in {label}: {examples}")
    if not nullable and numeric.isna().any():
        raise ValueError(f"Null value(s) found in required integer column {label}")
    return numeric.astype("Int64")


def require_nonblank_text(series: pd.Series, label: str) -> pd.Series:
    """Return stripped text after enforcing a v13 required-text constraint."""
    normalized = normalize_sql_nulls(series)
    stripped = normalized.astype("string").str.strip()
    bad = stripped.isna() | stripped.eq("")
    if bad.any():
        raise ValueError(
            f"Null/blank value(s) found in required text column {label}: "
            f"{int(bad.sum()):,} row(s)"
        )
    return stripped


def assert_unique_when_nonnull(
    df: pd.DataFrame,
    columns: list[str],
    nullable_columns: list[str],
    label: str,
) -> None:
    """Enforce a SQL UNIQUE key only where all nullable key components are non-NULL."""
    eligible = df.copy()
    for column in nullable_columns:
        eligible = eligible.loc[eligible[column].notna()]
    dup = eligible.duplicated(columns, keep=False)
    if dup.any():
        examples = eligible.loc[dup, columns].drop_duplicates().head(10)
        raise ValueError(
            f"{label} contains duplicate key(s) for {columns} where nullable key "
            f"components are populated.\nExamples:\n{examples.to_string(index=False)}"
        )


def validate_sql_null_behavior() -> None:
    """
    Lightweight runtime self-check for SQL NULL and integer validation behavior.

    Demonstrates that:
      - ["17", "\\N", "22"] passes when nullable=True;
      - ["17", "abc", "22"] fails numeric validation;
      - ["17", "\\N", "22"] fails when nullable=False.
    """
    nullable_example = pd.Series(["17", "\\N", "22"], dtype="string")
    converted = to_integer_checked(
        nullable_example,
        "internal SQL NULL self-check (nullable)",
        nullable=True,
    )
    if converted.iloc[0] != 17 or not pd.isna(converted.iloc[1]) or converted.iloc[2] != 22:
        raise RuntimeError("Internal SQL NULL self-check failed for nullable integers.")

    try:
        to_integer_checked(
            pd.Series(["17", "abc", "22"], dtype="string"),
            "internal SQL NULL self-check (invalid text)",
            nullable=True,
        )
    except ValueError as exc:
        if "Non-numeric value(s)" not in str(exc):
            raise RuntimeError(
                "Internal SQL NULL self-check raised the wrong error for invalid text."
            ) from exc
    else:
        raise RuntimeError("Internal SQL NULL self-check failed to reject invalid text.")

    try:
        to_integer_checked(
            nullable_example,
            "internal SQL NULL self-check (required)",
            nullable=False,
        )
    except ValueError as exc:
        if "Null value(s) found in required integer column" not in str(exc):
            raise RuntimeError(
                "Internal SQL NULL self-check raised the wrong error for required NULL."
            ) from exc
    else:
        raise RuntimeError("Internal SQL NULL self-check failed to reject required NULL.")


def assert_unique(df: pd.DataFrame, columns: list[str], label: str) -> None:
    dup = df.duplicated(columns, keep=False)
    if dup.any():
        examples = df.loc[dup, columns].drop_duplicates().head(10)
        raise ValueError(
            f"{label} contains duplicate key(s) for {columns}.\n"
            f"Examples:\n{examples.to_string(index=False)}"
        )


def collapse_identical_duplicates(
    df: pd.DataFrame,
    key_columns: list[str],
    value_columns: list[str],
    label: str,
) -> pd.DataFrame:
    """
    Collapse duplicate keys only if all relevant values are identical.

    Conflicting duplicates are unsafe for relational joins and trigger an error.
    """
    dup = df.duplicated(key_columns, keep=False)
    if not dup.any():
        return df

    duplicated = df.loc[dup, key_columns + value_columns].copy()
    nunique = duplicated.groupby(key_columns, dropna=False)[value_columns].nunique(dropna=False)
    conflicts = nunique.gt(1).any(axis=1)
    if conflicts.any():
        bad_keys = conflicts[conflicts].reset_index()[key_columns].head(10)
        raise ValueError(
            f"{label} has conflicting duplicate rows for key(s) {key_columns}.\n"
            f"Examples:\n{bad_keys.to_string(index=False)}"
        )

    return df.drop_duplicates(key_columns, keep="first").copy()


# -----------------------------------------------------------------------------
# Reference mappings
# -----------------------------------------------------------------------------


def validate_ortholog_table_v13(
    orthologs: pd.DataFrame,
    sequences: pd.DataFrame,
    genomes: pd.DataFrame,
    *,
    label: str = "orthologs.tsv",
) -> None:
    """Validate the v13 ortholog identity and foreign-key contract."""
    valid_genome_pks = set(genomes["genome_pk"].dropna().astype(int))
    invalid_genome = ~orthologs["genome_pk"].astype(int).isin(valid_genome_pks)
    if invalid_genome.any():
        examples = orthologs.loc[
            invalid_genome,
            ["ortholog_pk", "genome_pk", "sequence_pk", "odb12_id"],
        ].head(20)
        raise ValueError(
            f"{label} contains genome_pk value(s) that are absent from genomes.tsv.\n"
            f"Examples:\n{examples.to_string(index=False)}"
        )

    # v13 requires one ortholog row for every genome_pk + odb12_id combination.
    dup_genome_odb = orthologs.duplicated(["genome_pk", "odb12_id"], keep=False)
    if dup_genome_odb.any():
        examples = orthologs.loc[
            dup_genome_odb,
            ["ortholog_pk", "genome_pk", "sequence_pk", "odb12_id"],
        ].sort_values(["genome_pk", "odb12_id", "ortholog_pk"]).head(20)
        raise ValueError(
            f"{label} violates the v13 UNIQUE(genome_pk, odb12_id) constraint. "
            "The script will not select, collapse, or drop duplicate ortholog copies.\n"
            f"Examples:\n{examples.to_string(index=False)}"
        )

    # sequence_pk is nullable in v13, so only populated sequence keys participate
    # in the UNIQUE(sequence_pk, odb12_id) check, matching SQL NULL semantics.
    assert_unique_when_nonnull(
        orthologs,
        ["sequence_pk", "odb12_id"],
        ["sequence_pk"],
        f"{label} v13 UNIQUE(sequence_pk, odb12_id)",
    )

    sequence_map = sequences[["sequence_pk", "genome_pk"]].copy().rename(
        columns={"genome_pk": "genome_pk_from_sequence"}
    )
    checked = orthologs[["ortholog_pk", "genome_pk", "sequence_pk", "odb12_id"]].merge(
        sequence_map,
        on="sequence_pk",
        how="left",
        validate="many_to_one",
    )

    missing_sequence = checked["sequence_pk"].notna() & checked["genome_pk_from_sequence"].isna()
    if missing_sequence.any():
        examples = checked.loc[
            missing_sequence,
            ["ortholog_pk", "genome_pk", "sequence_pk", "odb12_id"],
        ].head(20)
        raise ValueError(
            f"{label} contains populated sequence_pk value(s) that are absent from "
            "sequences.tsv. Under v13, sequence_pk may be NULL, but when populated "
            "it must resolve to sequences.sequence_pk.\n"
            f"Examples:\n{examples.to_string(index=False)}"
        )

    conflict = (
        checked["sequence_pk"].notna()
        & checked["genome_pk_from_sequence"].notna()
        & (checked["genome_pk"] != checked["genome_pk_from_sequence"])
    )
    if conflict.any():
        examples = checked.loc[
            conflict,
            [
                "ortholog_pk",
                "sequence_pk",
                "genome_pk",
                "genome_pk_from_sequence",
                "odb12_id",
            ],
        ].head(20)
        raise ValueError(
            f"{label} and sequences.tsv disagree about genome_pk for populated "
            "sequence_pk value(s). v13 treats orthologs.genome_pk as authoritative, "
            "and the sequence relationship must agree.\n"
            f"Examples:\n{examples.to_string(index=False)}"
        )


def load_reference_tables(sql_tsv_dir: Path) -> dict[str, pd.DataFrame]:
    paths = {
        "orthologs": sql_tsv_dir / "compleasm_features" / "orthologs.tsv",
        "species_name_audit": sql_tsv_dir / "species_name_audit.tsv",
        "sequences": sql_tsv_dir / "sequences.tsv",
        "genomes": sql_tsv_dir / "genomes.tsv",
    }

    refs = {name: read_tsv(path, path.name) for name, path in paths.items()}

    require_columns(
        refs["orthologs"],
        ["ortholog_pk", "sequence_pk", "genome_pk", "odb12_id"],
        "orthologs.tsv",
    )
    require_columns(refs["species_name_audit"], ["species_normalized"], "species_name_audit.tsv")
    require_columns(refs["sequences"], ["sequence_pk", "genome_pk"], "sequences.tsv")
    require_columns(refs["genomes"], ["genome_pk", "accession_id"], "genomes.tsv")

    refs["orthologs"]["ortholog_pk"] = to_integer_checked(
        refs["orthologs"]["ortholog_pk"], "orthologs.tsv:ortholog_pk"
    )
    # v13 requires genome_pk directly on every ortholog row. It is no longer
    # reconstructable from sequence_pk.
    refs["orthologs"]["genome_pk"] = to_integer_checked(
        refs["orthologs"]["genome_pk"], "orthologs.tsv:genome_pk", nullable=False
    )
    # sequence_pk remains nullable and is used only for a consistency check.
    refs["orthologs"]["sequence_pk"] = to_integer_checked(
        refs["orthologs"]["sequence_pk"], "orthologs.tsv:sequence_pk", nullable=True
    )
    refs["orthologs"]["odb12_id"] = require_nonblank_text(
        refs["orthologs"]["odb12_id"], "orthologs.tsv:odb12_id"
    )

    refs["sequences"]["sequence_pk"] = to_integer_checked(
        refs["sequences"]["sequence_pk"], "sequences.tsv:sequence_pk"
    )
    refs["sequences"]["genome_pk"] = to_integer_checked(
        refs["sequences"]["genome_pk"], "sequences.tsv:genome_pk"
    )

    refs["genomes"]["genome_pk"] = to_integer_checked(
        refs["genomes"]["genome_pk"], "genomes.tsv:genome_pk"
    )
    refs["genomes"]["accession_id"] = refs["genomes"]["accession_id"].astype("string").str.strip()

    assert_unique(refs["orthologs"], ["ortholog_pk"], "orthologs.tsv")
    assert_unique(refs["sequences"], ["sequence_pk"], "sequences.tsv")
    assert_unique(refs["genomes"], ["genome_pk"], "genomes.tsv")
    assert_unique(refs["genomes"], ["accession_id"], "genomes.tsv")

    validate_ortholog_table_v13(
        refs["orthologs"], refs["sequences"], refs["genomes"]
    )

    return refs


def build_species_to_genome_map(
    species_audit: pd.DataFrame,
    genomes: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build a unique normalized-species -> genome_pk mapping.

    Canonical v13 route:
      species_name_audit.species_pk -> genomes.species_pk -> genomes.genome_pk

    Legacy compatibility:
      If species_pk is absent, a non-schema species_name_audit.genome_pk column
      may be used directly.
    """
    audit = species_audit.copy()
    audit["species_normalized_join"] = normalize_text(audit["species_normalized"])

    # Canonical v13 behavior: prefer species_pk whenever the schema column exists,
    # even if an older extra genome_pk column is also present in the TSV.
    if "species_pk" in audit.columns:
        require_columns(genomes, ["species_pk"], "genomes.tsv for v13 species mapping")

        audit["species_pk"] = to_integer_checked(
            audit["species_pk"], "species_name_audit.tsv:species_pk", nullable=True
        )
        audit_map = audit.loc[
            audit["species_normalized_join"].notna() & audit["species_pk"].notna(),
            ["species_normalized_join", "species_pk"],
        ].copy()
        audit_map = collapse_identical_duplicates(
            audit_map,
            ["species_normalized_join"],
            ["species_pk"],
            "species_name_audit.tsv v13 species mapping",
        )

        genomes_map = genomes.copy()
        genomes_map["species_pk"] = to_integer_checked(
            genomes_map["species_pk"], "genomes.tsv:species_pk", nullable=True
        )
        genomes_map = genomes_map.loc[genomes_map["species_pk"].notna()].copy()

        # If multiple genome rows share a species_pk, prefer exactly one current
        # genome when available; otherwise preserve the ambiguity and fail below.
        if "is_current" in genomes_map.columns:
            current_text = genomes_map["is_current"].astype("string").str.strip().str.lower()
            genomes_map["_is_current_bool"] = current_text.isin(
                {"1", "true", "t", "yes", "y"}
            )

            selected_parts = []
            for _, group in genomes_map.groupby("species_pk", dropna=False, sort=False):
                current = group.loc[group["_is_current_bool"]]
                if len(current) == 1:
                    selected_parts.append(current.iloc[[0]])
                else:
                    selected_parts.append(group)
            genomes_map = (
                pd.concat(selected_parts, ignore_index=True)
                if selected_parts
                else genomes_map.iloc[0:0]
            )

        ambiguous_species = genomes_map.duplicated(["species_pk"], keep=False)
        if ambiguous_species.any():
            examples = genomes_map.loc[
                ambiguous_species, ["species_pk", "genome_pk", "accession_id"]
            ].sort_values(["species_pk", "genome_pk"]).head(20)
            raise ValueError(
                "Cannot uniquely resolve species_pk -> genome_pk from genomes.tsv. "
                "Multiple candidate genomes remain for at least one species.\n"
                f"Examples:\n{examples.to_string(index=False)}"
            )

        genomes_map = genomes_map[["species_pk", "genome_pk"]].copy()
        resolved = audit_map.merge(
            genomes_map, on="species_pk", how="left", validate="many_to_one"
        )
        resolved = resolved.loc[
            resolved["genome_pk"].notna(), ["species_normalized_join", "genome_pk"]
        ]
        assert_unique(
            resolved,
            ["species_normalized_join"],
            "resolved v13 species-to-genome mapping",
        )
        return resolved

    # Backward compatibility only: v13 itself does not define genome_pk on
    # species_name_audit. Use this route only when the canonical species_pk column
    # is absent.
    if "genome_pk" in audit.columns:
        print(
            "[WARNING] species_name_audit.tsv lacks canonical v13 species_pk; "
            "using legacy direct genome_pk compatibility mode.",
            file=sys.stderr,
        )
        audit["genome_pk"] = to_integer_checked(
            audit["genome_pk"], "species_name_audit.tsv:genome_pk", nullable=True
        )
        direct = audit.loc[
            audit["species_normalized_join"].notna() & audit["genome_pk"].notna(),
            ["species_normalized_join", "genome_pk"],
        ].copy()
        direct = collapse_identical_duplicates(
            direct,
            ["species_normalized_join"],
            ["genome_pk"],
            "species_name_audit.tsv legacy direct species-to-genome mapping",
        )

        valid_genome_pks = set(genomes["genome_pk"].dropna().astype(int))
        invalid = ~direct["genome_pk"].astype(int).isin(valid_genome_pks)
        if invalid.any():
            examples = direct.loc[
                invalid, ["species_normalized_join", "genome_pk"]
            ].head(10)
            raise ValueError(
                "Legacy species_name_audit.tsv contains genome_pk value(s) absent "
                "from genomes.tsv.\n"
                f"Examples:\n{examples.to_string(index=False)}"
            )
        return direct

    raise ValueError(
        "species_name_audit.tsv must contain species_pk for the canonical "
        "gc3_dynamics_v13_reduced route, or legacy genome_pk for backward compatibility."
    )


def build_ortholog_reference(
    orthologs: pd.DataFrame,
    sequences: pd.DataFrame,
    genomes: pd.DataFrame,
) -> pd.DataFrame:
    """Build the v13-required genome_pk + odb12_id -> ortholog_pk mapping."""
    # The full v13 contract was already validated during reference loading. Run
    # it again here defensively if this function is called independently.
    validate_ortholog_table_v13(orthologs, sequences, genomes)

    mapping = orthologs[["genome_pk", "odb12_id", "ortholog_pk"]].copy()
    assert_unique(
        mapping,
        ["genome_pk", "odb12_id"],
        "v13 ortholog genome_pk + odb12_id mapping",
    )
    return mapping


def validate_v13_ortholog_contract_behavior() -> None:
    """Lightweight runtime checks for the v13 ortholog contract."""
    genomes = pd.DataFrame(
        {
            "genome_pk": pd.Series([10, 11], dtype="Int64"),
            "accession_id": pd.Series(["GCA_TEST10", "GCA_TEST11"], dtype="string"),
        }
    )
    sequences = pd.DataFrame(
        {
            "sequence_pk": pd.Series([100, 110], dtype="Int64"),
            "genome_pk": pd.Series([10, 11], dtype="Int64"),
        }
    )
    valid = pd.DataFrame(
        {
            "ortholog_pk": pd.Series([1, 2], dtype="Int64"),
            "genome_pk": pd.Series([10, 10], dtype="Int64"),
            "sequence_pk": pd.Series([100, pd.NA], dtype="Int64"),
            "odb12_id": pd.Series(["397969at8457", "414176at8457"], dtype="string"),
        }
    )
    validate_ortholog_table_v13(valid, sequences, genomes, label="internal v13 valid case")

    # Required genome_pk rejects SQL NULL / missing values.
    try:
        to_integer_checked(
            pd.Series(["10", "\\N"], dtype="string"),
            "internal v13 ortholog genome_pk check",
            nullable=False,
        )
    except ValueError:
        pass
    else:
        raise RuntimeError("Internal v13 self-check failed to reject NULL ortholog genome_pk.")

    # Required odb12_id rejects SQL NULL and whitespace-only identifiers.
    for values in (["397969at8457", "\\N"], ["397969at8457", "   "]):
        try:
            require_nonblank_text(
                pd.Series(values, dtype="string"),
                "internal v13 ortholog odb12_id check",
            )
        except ValueError:
            pass
        else:
            raise RuntimeError("Internal v13 self-check failed to reject NULL/blank odb12_id.")

    # Nullable sequence_pk remains valid.
    nullable_sequence = to_integer_checked(
        pd.Series(["100", "\\N"], dtype="string"),
        "internal v13 nullable sequence_pk check",
        nullable=True,
    )
    if nullable_sequence.iloc[0] != 100 or not pd.isna(nullable_sequence.iloc[1]):
        raise RuntimeError("Internal v13 self-check failed for nullable sequence_pk.")

    conflict = valid.iloc[[0]].copy()
    conflict["genome_pk"] = pd.Series([11], dtype="Int64")
    try:
        validate_ortholog_table_v13(conflict, sequences, genomes, label="internal conflict case")
    except ValueError as exc:
        if "disagree about genome_pk" not in str(exc):
            raise RuntimeError("Internal v13 conflict self-check raised the wrong error.") from exc
    else:
        raise RuntimeError("Internal v13 self-check failed to reject sequence/genome disagreement.")

    missing_sequence = valid.iloc[[0]].copy()
    missing_sequence["sequence_pk"] = pd.Series([999], dtype="Int64")
    try:
        validate_ortholog_table_v13(
            missing_sequence, sequences, genomes, label="internal missing-sequence case"
        )
    except ValueError as exc:
        if "absent from sequences.tsv" not in str(exc):
            raise RuntimeError("Internal missing-sequence self-check raised the wrong error.") from exc
    else:
        raise RuntimeError("Internal v13 self-check failed to reject missing sequence_pk.")

    duplicate = pd.DataFrame(
        {
            "ortholog_pk": pd.Series([1, 2], dtype="Int64"),
            "genome_pk": pd.Series([10, 10], dtype="Int64"),
            "sequence_pk": pd.Series([pd.NA, pd.NA], dtype="Int64"),
            "odb12_id": pd.Series(["397969at8457", "397969at8457"], dtype="string"),
        }
    )
    try:
        validate_ortholog_table_v13(duplicate, sequences, genomes, label="internal duplicate case")
    except ValueError as exc:
        if "UNIQUE(genome_pk, odb12_id)" not in str(exc):
            raise RuntimeError("Internal duplicate-key self-check raised the wrong error.") from exc
    else:
        raise RuntimeError("Internal v13 self-check failed to reject duplicate genome+ODB mapping.")


# -----------------------------------------------------------------------------
# Manual metadata validation (if the user has already created the files)
# -----------------------------------------------------------------------------


def validate_manual_metadata_files(sql_tsv_dir: Path) -> None:
    tree_path = sql_tsv_dir / "tree_files.tsv"
    runs_path = sql_tsv_dir / "nhphyml_runs.tsv"

    if runs_path.is_file():
        runs = read_tsv(runs_path, "nhphyml_runs.tsv")
        missing_canonical = [c for c in NHPHYML_RUNS_COLUMNS if c not in runs.columns]
        if missing_canonical == ["created_at"] and all(c in runs.columns for c in NHPHYML_RUNS_LEGACY_COLUMNS):
            print(
                "[WARNING] nhphyml_runs.tsv uses the legacy six-column v1-v3 header. "
                "The v13 schema also contains nullable created_at; treating that field as unspecified.",
                file=sys.stderr,
            )
        elif missing_canonical:
            raise ValueError(
                "nhphyml_runs.tsv is missing required canonical column(s): "
                + ", ".join(missing_canonical)
            )

        runs["nhphyml_run_pk"] = to_integer_checked(
            runs["nhphyml_run_pk"], "nhphyml_runs.tsv:nhphyml_run_pk"
        )
        runs["n_species"] = to_integer_checked(
            runs["n_species"], "nhphyml_runs.tsv:n_species", nullable=True
        )
        runs["n_genes"] = to_integer_checked(
            runs["n_genes"], "nhphyml_runs.tsv:n_genes", nullable=True
        )
        assert_unique(runs, ["nhphyml_run_pk"], "nhphyml_runs.tsv")
        assert_unique(runs, ["run_name"], "nhphyml_runs.tsv")
        expected = set(RUN_DIRS)
        observed = set(runs["nhphyml_run_pk"].dropna().astype(int))
        missing = sorted(expected - observed)
        if missing:
            print(
                f"[WARNING] nhphyml_runs.tsv exists but is missing run PK(s): {missing}",
                file=sys.stderr,
            )
        else:
            header_mode = "canonical v13" if "created_at" in runs.columns else "legacy v1-v3"
            print(f"[OK] nhphyml_runs.tsv ({header_mode}) has the required run PKs 1-4.")
    else:
        print(f"[WARNING] Manual file not found yet: {runs_path}", file=sys.stderr)

    if tree_path.is_file():
        trees = read_tsv(tree_path, "tree_files.tsv")
        require_columns(trees, TREE_FILES_COLUMNS, "tree_files.tsv")
        trees["tree_pk"] = to_integer_checked(trees["tree_pk"], "tree_files.tsv:tree_pk")
        trees["nhphyml_run_pk"] = to_integer_checked(
            trees["nhphyml_run_pk"], "tree_files.tsv:nhphyml_run_pk"
        )
        assert_unique(trees, ["tree_pk"], "tree_files.tsv")
        assert_unique(trees, ["nhphyml_run_pk", "tree_type"], "tree_files.tsv")
        print("[OK] tree_files.tsv has the required columns and unique keys.")
    else:
        print(f"[WARNING] Manual file not found yet: {tree_path}", file=sys.stderr)


# -----------------------------------------------------------------------------
# nhPhyML ortholog-level output
# -----------------------------------------------------------------------------


def process_ortholog_run(
    run_pk: int,
    run_dir: str,
    alignments_root: Path,
    species_to_genome: pd.DataFrame,
    ortholog_reference: pd.DataFrame,
) -> pd.DataFrame:
    gc3_path = alignments_root / run_dir / "12_gc_divergence" / "gc3_per_gene.tsv"
    gc4_path = alignments_root / run_dir / "11_gc4" / "detailed.tsv"

    gc3 = read_tsv(gc3_path, f"run {run_pk} gc3_per_gene.tsv")
    gc4 = read_tsv(gc4_path, f"run {run_pk} detailed.tsv")

    require_columns(gc3, ["ortholog", "node_or_species", "gc3"], str(gc3_path))
    require_columns(gc4, ["alignment", "sequence_id", "gc4"], str(gc4_path))

    # Parse ortholog string as requested.
    gc3["ortholog"] = gc3["ortholog"].astype("string").str.strip()
    gc3["odb12_id"] = gc3["ortholog"].str.split("_", n=1).str[0]
    gc3["alignment_match"] = gc3["ortholog"].str.split("__", n=1).str[0]
    gc3["node_or_species_join"] = normalize_text(gc3["node_or_species"])
    gc3["gc3_nhphyml"] = to_numeric_checked(gc3["gc3"], f"{gc3_path}:gc3")

    # Prepare GC4 keys as requested.
    gc4["alignment_match"] = (
        gc4["alignment"].astype("string").str.strip().str.replace(r"\.fasta$", "", regex=True)
    )
    gc4["node_or_species_join"] = normalize_text(gc4["sequence_id"])
    gc4["gc4"] = to_numeric_checked(gc4["gc4"], f"{gc4_path}:gc4")

    gc3 = collapse_identical_duplicates(
        gc3,
        ["alignment_match", "node_or_species_join"],
        ["odb12_id", "gc3_nhphyml"],
        f"run {run_pk} gc3_per_gene.tsv",
    )
    gc4 = collapse_identical_duplicates(
        gc4,
        ["alignment_match", "node_or_species_join"],
        ["gc4"],
        f"run {run_pk} detailed.tsv",
    )

    # GC3 is the row-driving table; GC4 is nullable in the database schema.
    merged = gc3[
        ["odb12_id", "alignment_match", "node_or_species_join", "gc3_nhphyml"]
    ].merge(
        gc4[["alignment_match", "node_or_species_join", "gc4"]],
        on=["alignment_match", "node_or_species_join"],
        how="left",
        validate="one_to_one",
    )

    missing_gc4 = int(merged["gc4"].isna().sum())
    if missing_gc4:
        print(
            f"[WARNING] Run {run_pk}: {missing_gc4:,} GC3 row(s) did not match a GC4 value; "
            "gc4 will be blank for those rows.",
            file=sys.stderr,
        )

    # Resolve extant species to genome_pk. Internal ancestral nodes remain null
    # after the left join and are then dropped exactly as requested.
    before_species_filter = len(merged)
    merged = merged.merge(
        species_to_genome,
        left_on="node_or_species_join",
        right_on="species_normalized_join",
        how="left",
        validate="many_to_one",
    )
    merged = merged.drop(columns=["species_normalized_join"])
    merged = merged.loc[merged["genome_pk"].notna()].copy()
    dropped_internal = before_species_filter - len(merged)
    print(
        f"[INFO] Run {run_pk}: retained {len(merged):,} extant-species rows; "
        f"dropped {dropped_internal:,} unresolved/internal-node rows."
    )
    merged["genome_pk"] = merged["genome_pk"].astype("Int64")

    # Resolve ortholog_pk by genome_pk + odb12_id.
    merged = merged.merge(
        ortholog_reference,
        on=["genome_pk", "odb12_id"],
        how="left",
        validate="many_to_one",
    )

    unresolved_ortholog = merged["ortholog_pk"].isna()
    if unresolved_ortholog.any():
        n_unresolved = int(unresolved_ortholog.sum())
        examples = merged.loc[
            unresolved_ortholog,
            ["node_or_species_join", "genome_pk", "odb12_id"],
        ].drop_duplicates().head(10)
        print(
            f"[WARNING] Run {run_pk}: dropping {n_unresolved:,} row(s) whose ortholog_pk "
            "could not be resolved from genome_pk + odb12_id.\n"
            f"Examples:\n{examples.to_string(index=False)}",
            file=sys.stderr,
        )
        merged = merged.loc[~unresolved_ortholog].copy()

    merged["ortholog_pk"] = merged["ortholog_pk"].astype("Int64")
    merged["nhphyml_run_pk"] = run_pk
    merged["status"] = pd.NA
    merged["created_at"] = format_mtime(gc3_path)

    out = merged[
        [
            "nhphyml_run_pk",
            "ortholog_pk",
            "genome_pk",
            "gc3_nhphyml",
            "gc4",
            "status",
            "created_at",
        ]
    ].copy()

    assert_unique(
        out,
        ["nhphyml_run_pk", "ortholog_pk", "genome_pk"],
        f"run {run_pk} nhphyml_ortholog_results",
    )
    return out


# -----------------------------------------------------------------------------
# nhPhyML species-summary output
# -----------------------------------------------------------------------------


def process_species_summary_run(
    run_pk: int,
    run_dir: str,
    alignments_root: Path,
    genomes: pd.DataFrame,
) -> pd.DataFrame:
    master_path = alignments_root / run_dir / "12_gc_divergence" / "master_output.tsv"
    master = read_tsv(master_path, f"run {run_pk} master_output.tsv")

    required = ["accession", "mean_gc3", "mean_gc4", "sd_gc3", "sd_gc4", "n_orthologs"]
    require_columns(master, required, str(master_path))

    master["accession_join"] = master["accession"].astype("string").str.strip()
    for col in ["mean_gc3", "mean_gc4", "sd_gc3", "sd_gc4"]:
        master[col] = to_numeric_checked(master[col], f"{master_path}:{col}")
    master["n_orthologs"] = to_integer_checked(
        master["n_orthologs"], f"{master_path}:n_orthologs", nullable=True
    )

    genome_map = genomes[["genome_pk", "accession_id"]].copy()
    genome_map["accession_join"] = genome_map["accession_id"].astype("string").str.strip()
    assert_unique(genome_map, ["accession_join"], "genomes.tsv accession mapping")

    before = len(master)
    merged = master.merge(
        genome_map[["accession_join", "genome_pk"]],
        on="accession_join",
        how="left",
        validate="many_to_one",
    )

    unresolved = merged["genome_pk"].isna()
    if unresolved.any():
        examples = merged.loc[unresolved, ["accession"]].drop_duplicates().head(10)
        print(
            f"[WARNING] Run {run_pk}: dropping {int(unresolved.sum()):,} master_output row(s) "
            "whose accession could not be resolved in genomes.tsv.\n"
            f"Examples:\n{examples.to_string(index=False)}",
            file=sys.stderr,
        )
        merged = merged.loc[~unresolved].copy()

    merged["genome_pk"] = merged["genome_pk"].astype("Int64")
    merged["nhphyml_run_pk"] = run_pk
    merged["created_at"] = format_mtime(master_path)

    out = merged[
        [
            "nhphyml_run_pk",
            "genome_pk",
            "mean_gc3",
            "mean_gc4",
            "sd_gc3",
            "sd_gc4",
            "n_orthologs",
            "created_at",
        ]
    ].copy()

    assert_unique(
        out,
        ["nhphyml_run_pk", "genome_pk"],
        f"run {run_pk} nhphyml_species_summary",
    )
    print(
        f"[INFO] Run {run_pk}: retained {len(out):,}/{before:,} species-summary row(s)."
    )
    return out


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert nhPhyML GC3/GC4 outputs into nhphyml_ortholog_results.tsv "
            "and nhphyml_species_summary.tsv."
        )
    )
    parser.add_argument(
        "--sql-tsv-dir",
        type=Path,
        default=Path("/genomes/records/sql_tsvs"),
        help=(
            "Directory containing orthologs.tsv, species_name_audit.tsv, "
            "sequences.tsv, and genomes.tsv. Default: /genomes/records/sql_tsvs"
        ),
    )
    parser.add_argument(
        "--ortholog-alignments-root",
        type=Path,
        default=None,
        help=(
            "Root used for gc3_per_gene.tsv and detailed.tsv run mappings. "
            "Default: derive from --sql-tsv-dir as "
            "<sql-tsv-dir parent>/compleasm/alignments."
        ),
    )
    parser.add_argument(
        "--summary-alignments-root",
        type=Path,
        default=None,
        help=(
            "Root used for master_output.tsv run mappings. "
            "Default: derive from --sql-tsv-dir as "
            "<sql-tsv-dir parent>/compleasm/alignments."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Default: the value of --sql-tsv-dir.",
    )
    return parser.parse_args()


def main() -> int:
    # This warning intentionally happens before argument parsing and before any
    # file I/O so it is the first thing the user sees when the script runs.
    print_manual_file_warning()

    args = parse_args()
    sql_tsv_dir = args.sql_tsv_dir.expanduser().resolve()

    # Derive both nhPhyML input roots from the selected SQL TSV directory unless
    # the user explicitly overrides them. For a standard project layout such as
    # /Users/rossoaa/projects/genomes/records/sql_tsvs, this resolves to
    # /Users/rossoaa/projects/genomes/records/compleasm/alignments.
    inferred_alignments_root = sql_tsv_dir.parent / "compleasm" / "alignments"
    ortholog_alignments_root = (
        args.ortholog_alignments_root.expanduser().resolve()
        if args.ortholog_alignments_root is not None
        else inferred_alignments_root
    )
    summary_alignments_root = (
        args.summary_alignments_root.expanduser().resolve()
        if args.summary_alignments_root is not None
        else inferred_alignments_root
    )

    output_dir = (args.output_dir or sql_tsv_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] SQL TSV directory: {sql_tsv_dir}")
    if args.ortholog_alignments_root is None or args.summary_alignments_root is None:
        print(f"[INFO] Inferred alignment root from --sql-tsv-dir: {inferred_alignments_root}")
    print(f"[INFO] Ortholog run root: {ortholog_alignments_root}")
    print(f"[INFO] Species-summary run root: {summary_alignments_root}")
    print(f"[INFO] Output directory: {output_dir}")
    print("[INFO] Compatibility target: gc3_dynamics_v13_reduced schema.")
    print("[INFO] v13 ortholog contract: genome_pk and odb12_id are required; sequence_pk is nullable.")
    print("[INFO] v13 ortholog identity: genome_pk + odb12_id must uniquely resolve ortholog_pk.")
    print(r"[INFO] SQL TSV null convention: \N is interpreted as NULL.")

    validate_sql_null_behavior()
    print("[OK] SQL NULL numeric-validation self-checks passed.")
    validate_v13_ortholog_contract_behavior()
    print("[OK] v13 ortholog-contract self-checks passed.")

    refs = load_reference_tables(sql_tsv_dir)
    validate_manual_metadata_files(sql_tsv_dir)

    species_to_genome = build_species_to_genome_map(
        refs["species_name_audit"], refs["genomes"]
    )
    ortholog_reference = build_ortholog_reference(
        refs["orthologs"], refs["sequences"], refs["genomes"]
    )

    print(f"[INFO] Species-to-genome mappings available: {len(species_to_genome):,}")
    print(f"[INFO] Genome+ODB12-to-ortholog mappings available: {len(ortholog_reference):,}")

    ortholog_parts = []
    summary_parts = []

    for run_pk, run_dir in RUN_DIRS.items():
        print(f"\n[INFO] Processing nhPhyML run {run_pk}: {run_dir}")
        ortholog_parts.append(
            process_ortholog_run(
                run_pk,
                run_dir,
                ortholog_alignments_root,
                species_to_genome,
                ortholog_reference,
            )
        )
        summary_parts.append(
            process_species_summary_run(
                run_pk,
                run_dir,
                summary_alignments_root,
                refs["genomes"],
            )
        )

    # Final ortholog-level table and deterministic auto-incremented PK.
    ortholog_results = pd.concat(ortholog_parts, ignore_index=True)
    ortholog_results = ortholog_results.sort_values(
        ["nhphyml_run_pk", "genome_pk", "ortholog_pk"], kind="stable"
    ).reset_index(drop=True)
    assert_unique(
        ortholog_results,
        ["nhphyml_run_pk", "ortholog_pk", "genome_pk"],
        "final nhphyml_ortholog_results",
    )
    ortholog_results.insert(
        0,
        "nhphyml_ortholog_result_pk",
        pd.Series(range(1, len(ortholog_results) + 1), dtype="int64"),
    )
    ortholog_results = ortholog_results[ORTHOLOG_OUTPUT_COLUMNS]

    # Final species-summary table and deterministic auto-incremented PK.
    species_summary = pd.concat(summary_parts, ignore_index=True)
    species_summary = species_summary.sort_values(
        ["nhphyml_run_pk", "genome_pk"], kind="stable"
    ).reset_index(drop=True)
    assert_unique(
        species_summary,
        ["nhphyml_run_pk", "genome_pk"],
        "final nhphyml_species_summary",
    )
    species_summary.insert(
        0,
        "nhphyml_species_summary_pk",
        pd.Series(range(1, len(species_summary) + 1), dtype="int64"),
    )
    species_summary = species_summary[SPECIES_SUMMARY_OUTPUT_COLUMNS]

    ortholog_output_path = output_dir / "nhphyml_ortholog_results.tsv"
    summary_output_path = output_dir / "nhphyml_species_summary.tsv"

    # SQL-oriented TSVs use the explicit MySQL NULL marker.
    ortholog_results.to_csv(ortholog_output_path, sep="\t", index=False, na_rep="\\N")
    species_summary.to_csv(summary_output_path, sep="\t", index=False, na_rep="\\N")

    print("\n[OK] Finished.")
    print(
        f"[OK] Wrote {len(ortholog_results):,} rows to {ortholog_output_path}"
    )
    print(
        f"[OK] Wrote {len(species_summary):,} rows to {summary_output_path}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"\n[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
