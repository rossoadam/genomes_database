#!/usr/bin/env python3
"""
06_nhphyml_to_tsvs_v2.py

Build two SQL-ready TSV files from nhPhyML/GC4 outputs:

    nhphyml_ortholog_results.tsv
    nhphyml_species_summary.tsv

The script uses pandas and follows the nhPhyML-related tables in
`gc3_dynamics_9_reduced.txt`.

Version v2 schema compatibility notes
-------------------------------------
The gc3_dynamics_9_reduced schema adds `genome_pk` directly to `orthologs` while
retaining `sequence_pk`. Some ortholog rows may therefore have a null
`sequence_pk` when the sequence cannot be resolved. For ortholog mapping this
script now:

1. Prefers `orthologs.genome_pk` when it is present.
2. Falls back to `sequences.tsv` through `sequence_pk` only when genome_pk is
   absent.
3. Verifies that direct and sequence-derived genome_pk values agree whenever
   both can be resolved.
4. Drops, with a warning, ortholog reference rows for which neither route can
   resolve a genome_pk.

The species_name_audit compatibility behavior from v1 is retained: the script
supports either direct `species_normalized` + `genome_pk` input or the schema
route `species_normalized` + `species_pk` -> `genomes.genome_pk`.

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
    print("nhphyml_runs.tsv fields:", file=sys.stderr)
    print("  " + "\t".join(NHPHYML_RUNS_COLUMNS), file=sys.stderr)
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
    require_file(path, label)
    try:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=True)
    except Exception as exc:
        raise RuntimeError(f"Could not read {label} as a TSV: {path}\n{exc}") from exc


def format_mtime(path: Path) -> str:
    """Return file modification time in SQL datetime format, using local time."""
    return datetime.fromtimestamp(os.path.getmtime(path)).strftime("%Y-%m-%d %H:%M:%S")


def normalize_text(series: pd.Series) -> pd.Series:
    """Trim and lowercase identifiers used for case-insensitive joins."""
    return series.astype("string").str.strip().str.lower()


def to_numeric_checked(series: pd.Series, label: str) -> pd.Series:
    """Convert a column to numeric and fail on non-empty, non-numeric values."""
    original = series.astype("string")
    numeric = pd.to_numeric(series, errors="coerce")
    bad = numeric.isna() & original.notna() & original.str.strip().ne("")
    if bad.any():
        examples = original.loc[bad].drop_duplicates().head(10).tolist()
        raise ValueError(
            f"Non-numeric value(s) found in {label}: {examples}"
            + (" ..." if bad.sum() > 10 else "")
        )
    return numeric


def to_integer_checked(series: pd.Series, label: str, nullable: bool = False) -> pd.Series:
    """Convert a column to pandas Int64, requiring integer-valued input."""
    numeric = to_numeric_checked(series, label)
    non_integer = numeric.notna() & ((numeric % 1) != 0)
    if non_integer.any():
        examples = series.loc[non_integer].drop_duplicates().head(10).tolist()
        raise ValueError(f"Non-integer value(s) found in {label}: {examples}")
    if not nullable and numeric.isna().any():
        raise ValueError(f"Null value(s) found in required integer column {label}")
    return numeric.astype("Int64")


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
    # gc3_dynamics_9_reduced permits ortholog identity to be carried directly
    # by genome_pk even when sequence_pk cannot be resolved. Therefore both
    # foreign-key columns are nullable at this TSV-building stage.
    refs["orthologs"]["sequence_pk"] = to_integer_checked(
        refs["orthologs"]["sequence_pk"], "orthologs.tsv:sequence_pk", nullable=True
    )
    refs["orthologs"]["genome_pk"] = to_integer_checked(
        refs["orthologs"]["genome_pk"], "orthologs.tsv:genome_pk", nullable=True
    )
    refs["orthologs"]["odb12_id"] = refs["orthologs"]["odb12_id"].astype("string").str.strip()

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

    return refs


def build_species_to_genome_map(
    species_audit: pd.DataFrame,
    genomes: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build a unique normalized-species -> genome_pk mapping.

    Supports:
      A) requested TSV layout: species_name_audit has genome_pk
      B) supplied schema layout: species_name_audit has species_pk; genomes has species_pk
    """
    audit = species_audit.copy()
    audit["species_normalized_join"] = normalize_text(audit["species_normalized"])

    if "genome_pk" in audit.columns:
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
            "species_name_audit.tsv direct species-to-genome mapping",
        )

        valid_genome_pks = set(genomes["genome_pk"].dropna().astype(int))
        invalid = ~direct["genome_pk"].astype(int).isin(valid_genome_pks)
        if invalid.any():
            examples = direct.loc[invalid, ["species_normalized_join", "genome_pk"]].head(10)
            print(
                f"[WARNING] Dropping {int(invalid.sum()):,} species_name_audit mapping(s) "
                "whose genome_pk does not exist in genomes.tsv.\n"
                f"Examples:\n{examples.to_string(index=False)}",
                file=sys.stderr,
            )
            direct = direct.loc[~invalid].copy()
        return direct

    if "species_pk" not in audit.columns:
        raise ValueError(
            "species_name_audit.tsv must contain either genome_pk (requested layout) "
            "or species_pk (layout in gc3_dynamics_9_reduced.txt)."
        )

    require_columns(genomes, ["species_pk"], "genomes.tsv for schema-compatible species mapping")

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
        "species_name_audit.tsv species mapping",
    )

    genomes_map = genomes.copy()
    genomes_map["species_pk"] = to_integer_checked(
        genomes_map["species_pk"], "genomes.tsv:species_pk", nullable=True
    )
    genomes_map = genomes_map.loc[genomes_map["species_pk"].notna()].copy()

    # If the database contains multiple genome records for the same species,
    # prefer the single current genome when is_current is available.
    if "is_current" in genomes_map.columns:
        current_text = genomes_map["is_current"].astype("string").str.strip().str.lower()
        genomes_map["_is_current_bool"] = current_text.isin({"1", "true", "t", "yes", "y"})

        selected_parts = []
        for _, group in genomes_map.groupby("species_pk", dropna=False, sort=False):
            current = group.loc[group["_is_current_bool"]]
            if len(current) == 1:
                selected_parts.append(current.iloc[[0]])
            else:
                selected_parts.append(group)
        genomes_map = pd.concat(selected_parts, ignore_index=True) if selected_parts else genomes_map.iloc[0:0]

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
    resolved = audit_map.merge(genomes_map, on="species_pk", how="left", validate="many_to_one")
    resolved = resolved.loc[resolved["genome_pk"].notna(), ["species_normalized_join", "genome_pk"]]
    assert_unique(resolved, ["species_normalized_join"], "resolved species-to-genome mapping")
    return resolved


def build_ortholog_reference(
    orthologs: pd.DataFrame,
    sequences: pd.DataFrame,
    genomes: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build the genome_pk + odb12_id -> ortholog_pk mapping.

    gc3_dynamics_9_reduced stores genome_pk directly in orthologs.tsv because
    some ortholog rows do not have a resolvable sequence_pk. Direct genome_pk
    is therefore authoritative when present; sequence_pk is retained as a
    fallback and as a cross-check when both routes are available.
    """
    mapping = orthologs[["ortholog_pk", "sequence_pk", "genome_pk", "odb12_id"]].copy()
    mapping = mapping.rename(columns={"genome_pk": "genome_pk_direct"})

    sequence_map = sequences[["sequence_pk", "genome_pk"]].copy().rename(
        columns={"genome_pk": "genome_pk_from_sequence"}
    )
    mapping = mapping.merge(
        sequence_map,
        on="sequence_pk",
        how="left",
        validate="many_to_one",
    )

    # If sequence_pk is populated but absent from sequences.tsv, this is no
    # longer automatically fatal because orthologs.genome_pk can still provide
    # the required genome identity. Report it so the unresolved sequence link
    # remains visible to the user.
    sequence_unresolved = (
        mapping["sequence_pk"].notna()
        & mapping["genome_pk_from_sequence"].isna()
    )
    if sequence_unresolved.any():
        examples = mapping.loc[
            sequence_unresolved,
            ["ortholog_pk", "sequence_pk", "genome_pk_direct", "odb12_id"],
        ].head(10)
        print(
            f"[WARNING] {int(sequence_unresolved.sum()):,} orthologs.tsv row(s) have a "
            "sequence_pk that is not present in sequences.tsv. Direct genome_pk will "
            "be used when available.\n"
            f"Examples:\n{examples.to_string(index=False)}",
            file=sys.stderr,
        )

    # When both routes resolve a genome, they must agree. A disagreement means
    # the reference tables are internally inconsistent and should not be
    # silently reconciled by this script.
    conflict = (
        mapping["genome_pk_direct"].notna()
        & mapping["genome_pk_from_sequence"].notna()
        & (mapping["genome_pk_direct"] != mapping["genome_pk_from_sequence"])
    )
    if conflict.any():
        examples = mapping.loc[
            conflict,
            [
                "ortholog_pk",
                "sequence_pk",
                "genome_pk_direct",
                "genome_pk_from_sequence",
                "odb12_id",
            ],
        ].head(20)
        raise ValueError(
            "orthologs.tsv and sequences.tsv disagree about genome_pk for at least "
            "one ortholog row. Direct genome_pk and the genome_pk resolved through "
            "sequence_pk must agree when both are available.\n"
            f"Examples:\n{examples.to_string(index=False)}"
        )

    mapping["genome_pk"] = mapping["genome_pk_direct"].combine_first(
        mapping["genome_pk_from_sequence"]
    )

    unresolved_genome = mapping["genome_pk"].isna()
    if unresolved_genome.any():
        examples = mapping.loc[
            unresolved_genome,
            ["ortholog_pk", "sequence_pk", "odb12_id"],
        ].head(10)
        print(
            f"[WARNING] Dropping {int(unresolved_genome.sum()):,} ortholog reference "
            "row(s) because genome_pk could not be resolved from either "
            "orthologs.genome_pk or sequence_pk -> sequences.genome_pk.\n"
            f"Examples:\n{examples.to_string(index=False)}",
            file=sys.stderr,
        )
        mapping = mapping.loc[~unresolved_genome].copy()

    mapping["genome_pk"] = mapping["genome_pk"].astype("Int64")

    # Every resolved genome_pk must exist in genomes.tsv. This mirrors the
    # nhphyml_species_summary foreign-key target and catches stale identifiers.
    valid_genome_pks = set(genomes["genome_pk"].dropna().astype(int))
    invalid_genome = ~mapping["genome_pk"].astype(int).isin(valid_genome_pks)
    if invalid_genome.any():
        examples = mapping.loc[
            invalid_genome,
            ["ortholog_pk", "genome_pk", "sequence_pk", "odb12_id"],
        ].head(20)
        raise ValueError(
            "orthologs.tsv contains resolved genome_pk value(s) that are absent from "
            "genomes.tsv.\n"
            f"Examples:\n{examples.to_string(index=False)}"
        )

    # The nhPhyML join requires one ortholog_pk for each genome+odb12_id.
    # Do not silently choose among duplicated ortholog copies.
    dup = mapping.duplicated(["genome_pk", "odb12_id"], keep=False)
    if dup.any():
        examples = mapping.loc[
            dup, ["genome_pk", "odb12_id", "ortholog_pk", "sequence_pk"]
        ].sort_values(["genome_pk", "odb12_id", "ortholog_pk"]).head(20)
        raise ValueError(
            "The genome_pk + odb12_id -> ortholog_pk mapping is ambiguous. At least "
            "one genome has multiple ortholog_pk values for the same odb12_id. The "
            "script will not choose one silently.\n"
            f"Examples:\n{examples.to_string(index=False)}"
        )

    return mapping[["genome_pk", "odb12_id", "ortholog_pk"]].copy()


# -----------------------------------------------------------------------------
# Manual metadata validation (if the user has already created the files)
# -----------------------------------------------------------------------------


def validate_manual_metadata_files(sql_tsv_dir: Path) -> None:
    tree_path = sql_tsv_dir / "tree_files.tsv"
    runs_path = sql_tsv_dir / "nhphyml_runs.tsv"

    if runs_path.is_file():
        runs = read_tsv(runs_path, "nhphyml_runs.tsv")
        require_columns(runs, NHPHYML_RUNS_COLUMNS, "nhphyml_runs.tsv")
        runs["nhphyml_run_pk"] = to_integer_checked(
            runs["nhphyml_run_pk"], "nhphyml_runs.tsv:nhphyml_run_pk"
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
            print("[OK] nhphyml_runs.tsv has the required columns and run PKs 1-4.")
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
        default=Path("/genomes/records/compleasm/alignments"),
        help=(
            "Root used for gc3_per_gene.tsv and detailed.tsv run mappings. "
            "Default follows the paths supplied in the request."
        ),
    )
    parser.add_argument(
        "--summary-alignments-root",
        type=Path,
        default=Path("/Users/rossoaa/projects/genomes/records/compleasm/alignments"),
        help=(
            "Root used for master_output.tsv run mappings. "
            "Default follows the paths supplied in the request."
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
    ortholog_alignments_root = args.ortholog_alignments_root.expanduser().resolve()
    summary_alignments_root = args.summary_alignments_root.expanduser().resolve()
    output_dir = (args.output_dir or sql_tsv_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] SQL TSV directory: {sql_tsv_dir}")
    print(f"[INFO] Ortholog run root: {ortholog_alignments_root}")
    print(f"[INFO] Species-summary run root: {summary_alignments_root}")
    print(f"[INFO] Output directory: {output_dir}")

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

    ortholog_results.to_csv(ortholog_output_path, sep="\t", index=False, na_rep="")
    species_summary.to_csv(summary_output_path, sep="\t", index=False, na_rep="")

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
