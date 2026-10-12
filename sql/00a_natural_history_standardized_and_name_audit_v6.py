#!/usr/bin/env python3
"""
00a_natural_history_standardized_and_name_audit_v6.py

Create a consolidated natural-history table and a species-name resolution audit
for the downstream GC3/body-size database workflow.

Canonical species key:
  species_normalized = genus_species in lowercase, e.g. crotalus_viridis

Natural-history sources supported:
  - genome-size/C-value summaries
  - Oskyrko et al. 2024 ReptTraits (--traits-oskyrko)
  - Title et al. macroevolutionary dataset (--mass-title)
  - Ji-style mass data (--mass-ji)
  - Bennett/GlobTherm thermal limits (--thermal-bennett)

Outputs:
  1. <genomes_dir>/records/natural_history_standardized/
       natural_history_standardized.tsv

     One row per focal project species. All natural-history predictors are
     left-joined to the focal species list using species_normalized.

  2. <genomes_dir>/records/sql_tsvs/
       species_name_audit.tsv

     Name-resolution audit rows for each ingested source dataset. Raw trait
     names are labeled:
       - matched: normalized name is a focal project species
       - unmatched_raw_trait: normalized raw name is not a focal project species
       - missing_trait: focal project species has no normalized-name match in
         that source dataset

The optional --allow-missing-natural-history flag permits unavailable source
datasets. When a source is unavailable, its output columns are blank and all
focal species are represented as missing_trait for that source in the audit.

Version v6 updates:
  - 00a now owns creation and maintenance of the canonical biological species_pk
    because natural_history.species_pk is the project schema's parent key for
    genomes.species_pk and species_name_audit.species_pk.
  - natural_history_standardized.tsv is the authoritative upstream handoff for
    species_pk -> species_normalized. It is written with species_pk as its first
    column and is consumed by 00b.
  - Existing species_pk values are reused from the focal-project input, an
    explicit --species-pk-map, or a previous natural_history_standardized.tsv.
    On a first run, unresolved focal species receive deterministic sequential
    species_pk values in species_normalized order; later new species receive IDs
    above the highest previously used species_pk.
  - species_name_audit_pk is generated only after the complete audit table has
    been assembled, using deterministic ordering before assignment of 1..N.
  - matched and missing_trait audit rows reuse the canonical species_pk;
    unmatched_raw_trait rows keep species_pk NULL and are serialized as \\N.
  - PK/FK relationships are validated before either output is written.
  - 00a never requires downstream genomes.tsv or natural_history.tsv, avoiding
    a circular dependency with 00b.
"""

import argparse
import csv
from decimal import Decimal, InvalidOperation
import re
from pathlib import Path
from typing import Optional

import pandas as pd


BAD_SPECIES_TOKENS = {"sp", "sp.", "spp", "spp.", "cf", "cf.", "aff", "aff.", "nr", "nr."}
HYBRID_MARKERS = {"x", "×"}

BASE_AUDIT_COLUMNS = [
    "source_dataset",
    "source_species_name",
    "species_normalized",
    "match_status",
]

AUDIT_COLUMNS = [
    "species_name_audit_pk",
    "species_pk",
    "source_dataset",
    "source_species_name",
    "species_normalized",
    "match_status",
]

SQL_NULL_TEXT = {"", r"\n", "null", "none", "nan", "na", "n/a", "<na>"}
VALID_AUDIT_STATUSES = {"matched", "unmatched_raw_trait", "missing_trait"}

STANDARDIZED_NATURAL_HISTORY_COLUMNS = [
    "species_pk",
    "species_normalized",
    "genus_species",
    "mean_c_value_pg",
    "mean_genome_mb",
    "mean_genome_gb",
    "genome_size",
    "genome_size_source",
    "mass_grams_oskyrko",
    "mass_grams_title",
    "mass_grams_ji",
    "mass_preferred",
    "mass_source_preferred",
    "max_female_length_svl_mm_oskyrko",
    "svl_mm_title",
    "largest_clutch_size",
    "mean_number_of_eggs_per_clutch",
    "number_clutches_per_year",
    "annual_number_of_eggs",
    "max_longevity_years",
    "range_size",
    "ct_min",
    "ctmin_metric",
    "ct_max",
    "ctmax_metric",
    "thermal_source",
    "n_thermal_rows_for_species",
]

OSKYRKO_COLUMNS = {
    "max_female_length_svl_mm_oskyrko": 'Maximum female length ("SVL", mm)/straight carapace length for turtles ("SCL", mm)',
    "largest_clutch_size": "Largest clutch size",
    "mean_number_of_eggs_per_clutch": "Mean number of offspring per litter or number of eggs per clutch",
    "max_longevity_years": "Maximum Longevity (years)",
    "mass_grams_oskyrko": "Maximum body mass (g)",
    "number_clutches_per_year": "Number of litters or clutches produced per year",
}

TITLE_COLUMNS = {
    "mass_grams_title": "mass",
    "svl_mm_title": "completeSVL",
    "range_size": "rangeSize",
}


def normalize_species_name(value: object) -> str:
    """
    Normalize a source species name to a canonical genus_species key.

    Rules:
      - keep only the first two valid taxonomic tokens
      - convert spaces/underscores to a common representation
      - lowercase the result
      - drop qualifiers such as cf., aff., sp., spp.
      - ignore subspecies/strain text after the binomial
    """
    if value is None or pd.isna(value):
        return "unknown_species"

    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "na", "n/a"}:
        return "unknown_species"

    text = text.replace("_", " ")
    text = re.sub(r"\([^)]*\)", " ", text)
    text = re.sub(r"\[[^]]*\]", " ", text)
    text = re.sub(r"[^A-Za-z×x\s.-]", " ", text)
    text = re.sub(r"\s+", " ", text).strip()

    parts: list[str] = []
    for token in text.split():
        clean = token.strip().strip(".").lower()
        if not clean:
            continue
        if clean in BAD_SPECIES_TOKENS:
            continue
        if clean in HYBRID_MARKERS:
            continue
        parts.append(clean)

    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    if len(parts) == 1:
        return parts[0]
    return "unknown_species"


def read_table_auto(path: Path) -> pd.DataFrame:
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith(".tsv") or suffixes.endswith(".txt"):
        return pd.read_csv(path, sep="\t", encoding="utf-8-sig")
    if suffixes.endswith(".csv"):
        return pd.read_csv(path, encoding="utf-8-sig")

    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        return pd.read_csv(handle, sep=dialect.delimiter)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def write_species_name_audit_tsv(df: pd.DataFrame, path: Path) -> None:
    """
    Write the SQL-facing species-name audit.

    species_pk remains a genuine missing value internally for unmatched raw
    trait rows. pandas serializes those missing values as MySQL-compatible \\N.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep=r"\N")


def choose_first_column(df: pd.DataFrame, candidates: list[str]) -> Optional[str]:
    cleaned_to_real = {str(c).strip().lower(): c for c in df.columns}
    for candidate in candidates:
        key = candidate.strip().lower()
        if key in cleaned_to_real:
            return cleaned_to_real[key]
    return None


def require_column(df: pd.DataFrame, candidates: list[str], label: str, path: Path) -> str:
    col = choose_first_column(df, candidates)
    if col is None:
        raise ValueError(
            f"Could not identify {label} column in {path}. "
            f"Columns found: {list(df.columns)}"
        )
    return col


def default_compleasm_metadata(genomes_dir: Path) -> Path:
    return genomes_dir / "records" / "compleasm" / "records" / "metadata.csv"


def natural_history_dir(genomes_dir: Path) -> Path:
    return genomes_dir / "records" / "natural_history"


def find_default_file(
    directory: Path,
    include_terms: list[str],
    required_columns: Optional[list[str]] = None,
) -> Optional[Path]:
    if not directory.exists():
        return None

    candidates = []
    for path in directory.iterdir():
        if not path.is_file():
            continue
        name = path.name.lower()
        if not name.endswith((".csv", ".tsv", ".txt")):
            continue
        if all(term.lower() in name for term in include_terms):
            candidates.append(path)

    if required_columns is None:
        return sorted(candidates)[0] if candidates else None

    for path in sorted(candidates):
        try:
            df = read_table_auto(path)
        except Exception:
            continue
        found = {str(c).strip().lower() for c in df.columns}
        if all(col.strip().lower() in found for col in required_columns):
            return path

    return None


def find_default_genome_size_file(directory: Path) -> Optional[Path]:
    if not directory.exists():
        return None

    for terms in [["genome", "size"], ["c", "value"], ["summary"]]:
        candidate = find_default_file(
            directory,
            terms,
            required_columns=["mean_c_value_pg"],
        )
        if candidate is not None:
            return candidate

    for path in sorted(directory.iterdir()):
        if not path.is_file() or not path.name.lower().endswith((".csv", ".tsv", ".txt")):
            continue
        try:
            df = read_table_auto(path)
        except Exception:
            continue
        if choose_first_column(df, ["mean_c_value_pg"]) is not None:
            return path

    return None



def is_sql_null(value: object) -> bool:
    """Return True for missing/SQL-null-like scalar values."""
    if value is None:
        return True
    try:
        if pd.isna(value):
            return True
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    return text == r"\N" or text.lower() in SQL_NULL_TEXT


def coerce_nullable_integer(series: pd.Series, label: str) -> pd.Series:
    """
    Strictly convert a PK/FK series to pandas nullable Int64.

    Decimal is used so integer-valued text such as "17.0" is accepted without
    silently rounding non-integer values or passing through Python strings such
    as "nan" or "<NA>".
    """
    values: list[object] = []
    invalid: list[tuple[object, object]] = []

    for idx, value in series.items():
        if is_sql_null(value):
            values.append(pd.NA)
            continue

        text = str(value).strip()
        try:
            number = Decimal(text)
        except InvalidOperation:
            invalid.append((idx, value))
            values.append(pd.NA)
            continue

        if not number.is_finite() or number != number.to_integral_value():
            invalid.append((idx, value))
            values.append(pd.NA)
            continue

        integer = int(number)
        if integer < -(2**63) or integer > (2**63 - 1):
            invalid.append((idx, value))
            values.append(pd.NA)
            continue

        values.append(integer)

    if invalid:
        examples = ", ".join(
            f"row={idx!r} value={value!r}" for idx, value in invalid[:10]
        )
        raise ValueError(
            f"{label} contains {len(invalid)} non-integer, non-numeric, or "
            f"out-of-range BIGINT value(s). Examples: {examples}"
        )

    return pd.Series(
        pd.array(values, dtype="Int64"),
        index=series.index,
        name=series.name,
    )


def validate_species_pk_mapping(mapping: pd.DataFrame, label: str) -> pd.DataFrame:
    """
    Validate an authoritative species_normalized -> species_pk mapping.

    Repeated identical rows are allowed (for example, several genomes for one
    biological species), but a normalized species may not map to multiple PKs.
    Because species_pk is the natural_history primary key, one PK also may not
    identify multiple normalized focal species.
    """
    required = {"species_normalized", "species_pk"}
    missing = required - set(mapping.columns)
    if missing:
        raise ValueError(
            f"{label} is missing required mapping columns: {sorted(missing)}"
        )

    work = mapping[["species_normalized", "species_pk"]].copy()
    work["species_normalized"] = work["species_normalized"].map(
        normalize_species_name
    )
    work["species_pk"] = coerce_nullable_integer(
        work["species_pk"], f"{label}: species_pk"
    )
    work = work.loc[
        (work["species_normalized"] != "unknown_species")
        & work["species_pk"].notna()
    ].copy()

    if (work["species_pk"] <= 0).any():
        bad = work.loc[work["species_pk"] <= 0, ["species_normalized", "species_pk"]]
        raise ValueError(
            f"{label} contains non-positive species_pk value(s). Examples:\n"
            + bad.head(10).to_string(index=False)
        )

    if work.empty:
        return work.reset_index(drop=True)

    n_pk_by_species = work.groupby("species_normalized")["species_pk"].nunique()
    conflicts = n_pk_by_species[n_pk_by_species > 1]
    if not conflicts.empty:
        examples = []
        for species in conflicts.index[:10]:
            values = sorted(
                int(v)
                for v in work.loc[
                    work["species_normalized"] == species, "species_pk"
                ].unique()
            )
            examples.append(f"{species} -> {values}")
        raise ValueError(
            f"{label} associates more than one species_pk with the same "
            f"species_normalized. Examples: {'; '.join(examples)}"
        )

    n_species_by_pk = work.groupby("species_pk")["species_normalized"].nunique()
    reverse_conflicts = n_species_by_pk[n_species_by_pk > 1]
    if not reverse_conflicts.empty:
        examples = []
        for species_pk in reverse_conflicts.index[:10]:
            values = sorted(
                work.loc[
                    work["species_pk"] == species_pk, "species_normalized"
                ].astype(str).unique()
            )
            examples.append(f"{int(species_pk)} -> {values}")
        raise ValueError(
            f"{label} reuses a species_pk for more than one species_normalized. "
            f"Examples: {'; '.join(examples)}"
        )

    return (
        work.drop_duplicates(subset=["species_normalized", "species_pk"])
        .sort_values(["species_normalized", "species_pk"], kind="mergesort")
        .reset_index(drop=True)
    )


def load_authoritative_species_pk_mapping(path: Path) -> pd.DataFrame:
    """
    Read an existing project TSV/CSV containing species_pk.

    The preferred normalized key is species_normalized (or an equivalent
    normalized-name column). genomes.tsv from the project schema is also
    supported by normalizing its species_ncbi column.
    """
    df = read_table_auto(path)

    species_pk_col = choose_first_column(df, ["species_pk"])
    if species_pk_col is None:
        raise ValueError(
            f"Authoritative species_pk mapping {path} has no species_pk column. "
            f"Columns found: {list(df.columns)}"
        )

    species_key_col = choose_first_column(
        df,
        [
            "species_normalized",
            "species_key",
            "genus_species",
            "genus_species_key",
            "normalized_species",
            "normalized_name",
        ],
    )
    species_name_col = choose_first_column(
        df,
        [
            "species_ncbi",
            "organism_name",
            "organismName",
            "Organism Name",
            "species",
            "species_name",
            "binomial",
            "binomial_2020",
        ],
    )

    if species_key_col is None and species_name_col is None:
        raise ValueError(
            f"Authoritative species_pk mapping {path} contains species_pk but "
            "no species name/key that can be normalized. "
            f"Columns found: {list(df.columns)}"
        )

    source_col = (
        species_key_col if species_key_col is not None else species_name_col
    )

    mapping = pd.DataFrame(
        {
            "species_normalized": df[source_col].map(normalize_species_name),
            "species_pk": df[species_pk_col],
        }
    )
    return validate_species_pk_mapping(
        mapping, f"authoritative species_pk mapping {path}"
    )


def unique_existing_paths(paths: list[Path]) -> list[Path]:
    """Return existing unique files while preserving candidate priority."""
    seen: set[Path] = set()
    existing: list[Path] = []

    for path in paths:
        resolved = path.expanduser().resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if resolved.exists() and resolved.is_file():
            existing.append(resolved)

    return existing


def resolve_canonical_species_pks(
    project_species: pd.DataFrame,
    project_source_path: Path,
    existing_natural_history_path: Path,
    explicit_species_pk_map: Optional[Path],
) -> tuple[pd.DataFrame, list[Path]]:
    """
    Resolve or create the canonical focal-species species_pk mapping owned by 00a.

    Authority order:
      1. species_pk already present in the focal-project input;
      2. an explicit --species-pk-map, when supplied;
      3. a previous natural_history_standardized.tsv emitted by 00a;
      4. deterministic generation by 00a for still-unresolved focal species.

    Crucially, this function never depends on downstream SQL TSVs such as
    genomes.tsv or natural_history.tsv. natural_history_standardized.tsv is the
    upstream canonical handoff that 00b must consume.
    """
    mapping_frames: list[pd.DataFrame] = []
    mapping_sources: list[Path] = []

    project_mapping = project_species.loc[
        project_species["species_pk"].notna(),
        ["species_normalized", "species_pk"],
    ].copy()

    if not project_mapping.empty:
        project_mapping = validate_species_pk_mapping(
            project_mapping,
            f"focal project source {project_source_path}",
        )
        project_mapping["_mapping_source"] = str(project_source_path)
        mapping_frames.append(project_mapping)
        mapping_sources.append(project_source_path)
        print(
            f"[INFO] Focal project input supplies species_pk for "
            f"{project_mapping['species_normalized'].nunique():,} focal species"
        )

    explicit_resolved: Optional[Path] = None
    if explicit_species_pk_map is not None:
        explicit_resolved = explicit_species_pk_map.expanduser().resolve()
        if not explicit_resolved.exists():
            raise FileNotFoundError(
                f"--species-pk-map does not exist: {explicit_resolved}"
            )
        mapping = load_authoritative_species_pk_mapping(explicit_resolved)
        if mapping.empty:
            raise ValueError(
                f"--species-pk-map contains no usable species_pk mappings: "
                f"{explicit_resolved}"
            )
        mapping["_mapping_source"] = str(explicit_resolved)
        mapping_frames.append(mapping)
        mapping_sources.append(explicit_resolved)
        print(
            f"[INFO] Loaded {mapping['species_normalized'].nunique():,} species_pk "
            f"mapping(s) from --species-pk-map {explicit_resolved}"
        )

    # Reuse the previous upstream handoff when available. A legacy file without
    # species_pk is ignored rather than treated as an error.
    previous_path = existing_natural_history_path.expanduser().resolve()
    project_source_resolved = project_source_path.expanduser().resolve()
    if (
        previous_path.exists()
        and previous_path.is_file()
        and previous_path != project_source_resolved
        and previous_path != explicit_resolved
    ):
        previous_df = read_table_auto(previous_path)
        if choose_first_column(previous_df, ["species_pk"]) is None:
            print(
                f"[INFO] Existing {previous_path} has no species_pk column "
                "(legacy output); canonical IDs will be initialized by 00a"
            )
        else:
            previous_mapping = load_authoritative_species_pk_mapping(previous_path)
            if not previous_mapping.empty:
                previous_mapping["_mapping_source"] = str(previous_path)
                mapping_frames.append(previous_mapping)
                mapping_sources.append(previous_path)
                print(
                    f"[INFO] Reusing previous 00a species_pk mapping for "
                    f"{previous_mapping['species_normalized'].nunique():,} species "
                    f"from {previous_path}"
                )

    if mapping_frames:
        combined = pd.concat(mapping_frames, ignore_index=True)

        n_pk_by_species = combined.groupby("species_normalized")["species_pk"].nunique()
        conflicts = n_pk_by_species[n_pk_by_species > 1]
        if not conflicts.empty:
            examples = []
            for species in conflicts.index[:10]:
                subset = combined.loc[combined["species_normalized"] == species]
                pks = sorted(int(v) for v in subset["species_pk"].dropna().unique())
                sources = sorted(set(subset["_mapping_source"].astype(str)))
                examples.append(f"{species} -> {pks} from {sources}")
            raise ValueError(
                "Conflicting species_pk assignments were found across upstream "
                f"canonical mappings. Examples: {'; '.join(examples)}"
            )

        n_species_by_pk = combined.groupby("species_pk")["species_normalized"].nunique()
        reverse_conflicts = n_species_by_pk[n_species_by_pk > 1]
        if not reverse_conflicts.empty:
            examples = []
            for species_pk in reverse_conflicts.index[:10]:
                subset = combined.loc[combined["species_pk"] == species_pk]
                names = sorted(set(subset["species_normalized"].astype(str)))
                examples.append(f"{int(species_pk)} -> {names}")
            raise ValueError(
                "An upstream species_pk is associated with multiple normalized "
                f"species. Examples: {'; '.join(examples)}"
            )

        canonical_existing = (
            combined[["species_normalized", "species_pk"]]
            .drop_duplicates()
            .sort_values(["species_normalized", "species_pk"], kind="mergesort")
            .drop_duplicates(subset=["species_normalized"], keep="first")
            .reset_index(drop=True)
        )
        used_pks = {int(v) for v in combined["species_pk"].dropna()}
    else:
        canonical_existing = pd.DataFrame(
            columns=["species_normalized", "species_pk"]
        )
        used_pks: set[int] = set()

    canonical_lookup = (
        canonical_existing.set_index("species_normalized")["species_pk"]
        if not canonical_existing.empty
        else pd.Series(dtype="Int64")
    )

    resolved = project_species.copy()
    mapped_species_pk = resolved["species_normalized"].map(canonical_lookup)

    existing = resolved["species_pk"].notna()
    disagreement = (
        existing
        & mapped_species_pk.notna()
        & (
            resolved["species_pk"].astype("Int64")
            != mapped_species_pk.astype("Int64")
        )
    )
    if disagreement.any():
        examples = resolved.loc[
            disagreement, ["species_normalized", "species_pk"]
        ].copy()
        examples["upstream_species_pk"] = mapped_species_pk.loc[
            disagreement
        ].values
        raise ValueError(
            "Focal-project species_pk conflicts with an existing upstream "
            "canonical mapping. Examples:\n"
            + examples.head(10).to_string(index=False)
        )

    fill_mask = resolved["species_pk"].isna() & mapped_species_pk.notna()
    resolved.loc[fill_mask, "species_pk"] = mapped_species_pk.loc[fill_mask]
    resolved["species_pk"] = coerce_nullable_integer(
        resolved["species_pk"], "resolved focal species_pk"
    )

    unresolved_species = sorted(
        set(
            resolved.loc[
                resolved["species_pk"].isna(), "species_normalized"
            ].astype(str)
        )
    )

    if unresolved_species:
        next_pk = max(used_pks, default=0) + 1
        generated_lookup = {
            species: next_pk + offset
            for offset, species in enumerate(unresolved_species)
        }
        unresolved_mask = resolved["species_pk"].isna()
        resolved.loc[unresolved_mask, "species_pk"] = resolved.loc[
            unresolved_mask, "species_normalized"
        ].map(generated_lookup)
        resolved["species_pk"] = coerce_nullable_integer(
            resolved["species_pk"], "generated focal species_pk"
        )
        print(
            f"[INFO] 00a generated canonical species_pk for "
            f"{len(unresolved_species):,} unresolved focal species "
            f"({min(generated_lookup.values())}..{max(generated_lookup.values())}) "
            "in deterministic species_normalized order"
        )

    final_mapping = validate_species_pk_mapping(
        resolved[["species_normalized", "species_pk"]],
        "00a canonical focal species mapping",
    )
    if len(final_mapping) != len(resolved):
        raise ValueError(
            "Resolved focal species mapping is not one-to-one by species_normalized."
        )

    print(
        f"[OK] Canonical species_pk resolved for all {len(resolved):,} focal species"
    )
    print(
        "[INFO] Authority: 00a owns species_pk; "
        "natural_history_standardized.tsv is the canonical downstream handoff"
    )
    if mapping_sources:
        print(
            "[INFO] Reused species_pk source(s): "
            + "; ".join(str(p) for p in dict.fromkeys(mapping_sources))
        )

    return resolved, list(dict.fromkeys(mapping_sources))


def finalize_species_name_audit(
    audit: pd.DataFrame,
    project_species: pd.DataFrame,
) -> pd.DataFrame:
    """
    Populate canonical species_pk and assign deterministic species_name_audit_pk.

    species_name_audit_pk is intentionally assigned only after all audit frames
    have already been concatenated by build_standardized_natural_history().
    """
    missing_columns = [
        column for column in BASE_AUDIT_COLUMNS if column not in audit.columns
    ]
    if missing_columns:
        raise ValueError(
            f"Final audit is missing required columns: {missing_columns}"
        )

    unknown_statuses = sorted(
        set(audit["match_status"].dropna().astype(str))
        - VALID_AUDIT_STATUSES
    )
    if unknown_statuses:
        raise ValueError(
            f"Audit contains unsupported match_status values: {unknown_statuses}"
        )

    canonical_mapping = validate_species_pk_mapping(
        project_species[["species_normalized", "species_pk"]],
        "resolved focal species mapping",
    )
    canonical_lookup = canonical_mapping.set_index("species_normalized")[
        "species_pk"
    ]

    final_audit = audit[BASE_AUDIT_COLUMNS].copy()
    final_audit["species_pk"] = pd.Series(
        pd.array([pd.NA] * len(final_audit), dtype="Int64"),
        index=final_audit.index,
    )

    focal_status = final_audit["match_status"].isin(
        ["matched", "missing_trait"]
    )
    final_audit.loc[focal_status, "species_pk"] = final_audit.loc[
        focal_status, "species_normalized"
    ].map(canonical_lookup)
    final_audit["species_pk"] = coerce_nullable_integer(
        final_audit["species_pk"], "species_name_audit.species_pk"
    )

    # unmatched_raw_trait is intentionally NULL even if its normalized string
    # happens to resemble a key in another source.
    unmatched = final_audit["match_status"] == "unmatched_raw_trait"
    final_audit.loc[unmatched, "species_pk"] = pd.NA
    final_audit["species_pk"] = final_audit["species_pk"].astype("Int64")

    # Sort before assigning the PK so the same audit rows receive stable IDs on
    # reruns. mergesort preserves input order for completely identical rows.
    final_audit = final_audit.sort_values(
        [
            "source_dataset",
            "species_normalized",
            "match_status",
            "source_species_name",
        ],
        kind="mergesort",
        na_position="last",
    ).reset_index(drop=True)

    final_audit.insert(
        0,
        "species_name_audit_pk",
        pd.Series(
            pd.array(range(1, len(final_audit) + 1), dtype="Int64")
        ),
    )

    final_audit = final_audit[AUDIT_COLUMNS]
    validate_species_name_audit(final_audit, canonical_mapping)
    return final_audit


def validate_species_name_audit(
    audit: pd.DataFrame,
    canonical_mapping: pd.DataFrame,
) -> None:
    """Run the v6 PK/FK and match-status validations before file output."""
    if list(audit.columns) != AUDIT_COLUMNS:
        raise ValueError(
            "species_name_audit column order does not match the SQL schema. "
            f"Expected {AUDIT_COLUMNS}; observed {list(audit.columns)}"
        )
    print("[OK] species_name_audit columns match the required SQL order")

    audit_pk = coerce_nullable_integer(
        audit["species_name_audit_pk"], "species_name_audit_pk"
    )
    if audit_pk.isna().any():
        raise ValueError("species_name_audit_pk contains NULL values")
    if audit_pk.duplicated().any():
        raise ValueError("species_name_audit_pk contains duplicate values")

    observed = [int(value) for value in audit_pk.tolist()]
    expected = list(range(1, len(audit) + 1))
    if observed != expected:
        raise ValueError(
            "species_name_audit_pk does not run sequentially from 1..N"
        )
    print(
        f"[OK] species_name_audit_pk is non-null, numeric, unique, "
        f"and sequential from 1..{len(audit):,}"
    )

    species_pk = coerce_nullable_integer(
        audit["species_pk"], "species_name_audit.species_pk"
    )
    print(
        f"[OK] All {int(species_pk.notna().sum()):,} non-null species_pk "
        "values are numeric and integer-valued"
    )

    matched = audit["match_status"] == "matched"
    missing_trait = audit["match_status"] == "missing_trait"
    unmatched = audit["match_status"] == "unmatched_raw_trait"

    if audit.loc[matched, "species_pk"].isna().any():
        bad = audit.loc[
            matched & audit["species_pk"].isna(),
            ["source_dataset", "species_normalized"],
        ].head(10)
        raise ValueError(
            "matched audit rows must have species_pk. Examples:\n"
            + bad.to_string(index=False)
        )

    if audit.loc[missing_trait, "species_pk"].isna().any():
        bad = audit.loc[
            missing_trait & audit["species_pk"].isna(),
            ["source_dataset", "species_normalized"],
        ].head(10)
        raise ValueError(
            "missing_trait audit rows must have species_pk. Examples:\n"
            + bad.to_string(index=False)
        )

    if audit.loc[unmatched, "species_pk"].notna().any():
        bad = audit.loc[
            unmatched & audit["species_pk"].notna(),
            ["source_dataset", "species_normalized", "species_pk"],
        ].head(10)
        raise ValueError(
            "unmatched_raw_trait audit rows must have NULL species_pk. "
            "Examples:\n" + bad.to_string(index=False)
        )

    print(
        "[OK] species_pk status rules pass: matched and missing_trait are "
        "populated; unmatched_raw_trait is NULL"
    )

    canonical_lookup = canonical_mapping.set_index("species_normalized")[
        "species_pk"
    ]
    focal_rows = audit.loc[matched | missing_trait].copy()
    focal_rows["_canonical_species_pk"] = focal_rows[
        "species_normalized"
    ].map(canonical_lookup)

    unresolved = focal_rows["_canonical_species_pk"].isna()
    if unresolved.any():
        bad = focal_rows.loc[
            unresolved,
            ["source_dataset", "species_normalized", "match_status"],
        ].head(10)
        raise ValueError(
            "Audit rows marked matched/missing_trait contain normalized species "
            "that are absent from the canonical focal-species mapping. Examples:\n"
            + bad.to_string(index=False)
        )

    mismatch = (
        focal_rows["species_pk"].astype("Int64")
        != focal_rows["_canonical_species_pk"].astype("Int64")
    )
    if mismatch.fillna(True).any():
        bad = focal_rows.loc[
            mismatch.fillna(True),
            [
                "species_normalized",
                "species_pk",
                "_canonical_species_pk",
                "match_status",
            ],
        ].head(10)
        raise ValueError(
            "species_name_audit species_pk conflicts with the authoritative "
            "focal-species mapping. Examples:\n"
            + bad.to_string(index=False)
        )

    per_species = focal_rows.groupby("species_normalized")[
        "species_pk"
    ].nunique(dropna=True)
    conflicting = per_species[per_species != 1]
    if not conflicting.empty:
        raise ValueError(
            "A focal species_normalized does not map to exactly one species_pk "
            f"in the audit: {conflicting.to_dict()}"
        )

    canonical_species = set(canonical_mapping["species_normalized"].astype(str))
    represented_species = set(focal_rows["species_normalized"].astype(str))
    absent = sorted(canonical_species - represented_species)
    if absent:
        raise ValueError(
            "Some focal species are absent from all matched/missing_trait audit "
            "rows: " + ", ".join(absent[:20])
        )

    print(
        f"[OK] Each of the {len(canonical_species):,} focal normalized species "
        "maps to exactly one canonical species_pk"
    )
    print(
        "[OK] species_name_audit species_pk values agree with the authoritative "
        "focal-species mapping"
    )


def load_project_species(source_path: Path, source_label: str) -> pd.DataFrame:
    df = read_table_auto(source_path)

    accession_col = choose_first_column(
        df,
        [
            "accession",
            "accession_root",
            "assembly_accession",
            "assemblyAccession",
            "Assembly Accession",
            "genome_accession",
            "ncbi_accession",
        ],
    )
    species_key_col = choose_first_column(
        df,
        [
            "species_normalized",
            "species_key",
            "genus_species",
            "genus_species_key",
            "normalized_species",
            "normalized_name",
        ],
    )
    organism_col = choose_first_column(
        df,
        [
            "organism_name",
            "organismName",
            "Organism Name",
            "species",
            "species_name",
            "binomial",
            "binomial_2020",
        ],
    )
    species_pk_col = choose_first_column(df, ["species_pk"])

    if species_key_col is None and organism_col is None:
        raise ValueError(
            f"Could not find species key or organism/species column in {source_path}. "
            f"Columns found: {list(df.columns)}"
        )

    source_col = species_key_col if species_key_col is not None else organism_col

    out = pd.DataFrame(index=df.index)
    out["species_normalized"] = df[source_col].map(normalize_species_name)

    if organism_col is not None:
        out["organism_name"] = df[organism_col].astype("string").str.strip()
    else:
        out["organism_name"] = out["species_normalized"].str.replace(
            "_", " ", regex=False
        )

    if accession_col is not None:
        out["accession"] = df[accession_col].astype("string").str.strip()
    else:
        out["accession"] = pd.NA

    if species_pk_col is not None:
        out["species_pk"] = coerce_nullable_integer(
            df[species_pk_col],
            f"{source_path}: species_pk",
        )
    else:
        out["species_pk"] = pd.Series(
            pd.array([pd.NA] * len(out), dtype="Int64"),
            index=out.index,
        )

    out["project_species_source"] = source_label
    out = out[out["species_normalized"].astype(str).str.len() > 0]
    out = out[out["species_normalized"] != "unknown_species"]

    # Check all duplicate focal-source rows before collapsing them. If the same
    # normalized biological species carries multiple PKs, fail rather than
    # choosing one arbitrarily.
    present_mapping = out.loc[
        out["species_pk"].notna(),
        ["species_normalized", "species_pk"],
    ]
    if not present_mapping.empty:
        validate_species_pk_mapping(
            present_mapping,
            f"focal project source {source_path}",
        )

    # Preserve v4's one-row-per-focal-species behavior, but do not lose a
    # canonical species_pk merely because the first duplicate row is blank.
    rows = []
    for species_normalized, group in out.groupby(
        "species_normalized", sort=False
    ):
        first = group.iloc[0].copy()
        nonnull_pk = group["species_pk"].dropna()
        first["species_pk"] = (
            nonnull_pk.iloc[0] if len(nonnull_pk) else pd.NA
        )
        rows.append(first)

    out = pd.DataFrame(rows).reset_index(drop=True)
    out["species_pk"] = coerce_nullable_integer(
        out["species_pk"], f"{source_path}: collapsed species_pk"
    )
    out["genus_species"] = out["species_normalized"]

    return out[
        [
            "species_normalized",
            "genus_species",
            "organism_name",
            "accession",
            "project_species_source",
            "species_pk",
        ]
    ]


def audit_raw_species(
    raw_species_names: pd.Series,
    source_dataset: str,
    focal_species: set[str],
) -> pd.DataFrame:
    """
    Audit raw names from one source dataset and add missing focal species rows.

    Each raw input row is retained in the audit. Duplicate raw names therefore
    remain visible rather than being collapsed away with the trait data.
    """
    raw_names = raw_species_names.astype("string")
    normalized = raw_names.map(normalize_species_name)

    audit = pd.DataFrame(
        {
            "source_dataset": source_dataset,
            "source_species_name": raw_names.fillna(""),
            "species_normalized": normalized,
        }
    )
    audit["match_status"] = audit["species_normalized"].map(
        lambda x: "matched" if x in focal_species else "unmatched_raw_trait"
    )

    matched_focal = set(
        audit.loc[audit["match_status"] == "matched", "species_normalized"].astype(str)
    )
    missing = sorted(focal_species - matched_focal)

    if missing:
        missing_df = pd.DataFrame(
            {
                "source_dataset": source_dataset,
                "source_species_name": "",
                "species_normalized": missing,
                "match_status": "missing_trait",
            }
        )
        audit = pd.concat([audit, missing_df], ignore_index=True)

    return audit[BASE_AUDIT_COLUMNS]


def audit_missing_dataset(
    source_dataset: str,
    focal_species: set[str],
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source_dataset": source_dataset,
            "source_species_name": "",
            "species_normalized": sorted(focal_species),
            "match_status": "missing_trait",
        },
        columns=BASE_AUDIT_COLUMNS,
    )


def collapse_numeric_by_species(
    df: pd.DataFrame,
    value_cols: list[str],
) -> pd.DataFrame:
    """Collapse duplicate species rows by taking the first non-null value per column."""
    rows = []
    for species_normalized, group in df.groupby("species_normalized", sort=True):
        row = {"species_normalized": species_normalized}
        for col in value_cols:
            nonnull = group[col].dropna()
            row[col] = nonnull.iloc[0] if len(nonnull) else pd.NA
        rows.append(row)
    return pd.DataFrame(rows)


def standardize_genome_size(
    path: Path,
    focal_species: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = read_table_auto(path)

    cvalue_col = choose_first_column(df, ["mean_c_value_pg"])
    species_key_col = choose_first_column(
        df,
        ["species_normalized", "species_key", "genus_species"],
    )
    species_name_col = choose_first_column(
        df,
        ["species_name", "organism_name", "organismName"],
    )
    mean_mb_col = choose_first_column(df, ["mean_genome_mb"])
    mean_gb_col = choose_first_column(df, ["mean_genome_gb"])

    if cvalue_col is None:
        raise ValueError(
            f"Could not find mean_c_value_pg in {path}. "
            f"Columns found: {list(df.columns)}"
        )
    if species_key_col is None and species_name_col is None:
        raise ValueError(
            f"Could not find species key or species_name in {path}. "
            f"Columns found: {list(df.columns)}"
        )

    raw_name_col = species_name_col if species_name_col is not None else species_key_col
    audit = audit_raw_species(df[raw_name_col], "c_value", focal_species)

    source_names = df[species_key_col] if species_key_col is not None else df[species_name_col]

    out = pd.DataFrame()
    out["species_normalized"] = source_names.map(normalize_species_name)
    out["mean_c_value_pg"] = pd.to_numeric(df[cvalue_col], errors="coerce")
    out["mean_genome_mb"] = (
        pd.to_numeric(df[mean_mb_col], errors="coerce")
        if mean_mb_col is not None
        else pd.NA
    )
    out["mean_genome_gb"] = (
        pd.to_numeric(df[mean_gb_col], errors="coerce")
        if mean_gb_col is not None
        else pd.NA
    )
    out["genome_size_source"] = str(path)

    out = out[out["species_normalized"] != "unknown_species"]
    out = collapse_numeric_by_species(
        out,
        ["mean_c_value_pg", "mean_genome_mb", "mean_genome_gb", "genome_size_source"],
    )
    out["genome_size"] = out["mean_c_value_pg"]

    return out[
        [
            "species_normalized",
            "mean_c_value_pg",
            "mean_genome_mb",
            "mean_genome_gb",
            "genome_size",
            "genome_size_source",
        ]
    ], audit


def standardize_oskyrko_traits(
    path: Path,
    focal_species: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = read_table_auto(path)
    species_col = require_column(
        df,
        ["Species", "species", "species_normalized", "genus_species", "species_name", "binomial"],
        "Oskyrko/ReptTraits species",
        path,
    )

    missing = [
        source_col
        for source_col in OSKYRKO_COLUMNS.values()
        if choose_first_column(df, [source_col]) is None
    ]
    if missing:
        raise ValueError(
            f"Missing required Oskyrko/ReptTraits columns in {path}: {missing}. "
            f"Columns found: {list(df.columns)}"
        )

    audit = audit_raw_species(df[species_col], "oskyrko", focal_species)

    out = pd.DataFrame()
    out["species_normalized"] = df[species_col].map(normalize_species_name)

    for new_col, source_col in OSKYRKO_COLUMNS.items():
        real_col = require_column(df, [source_col], source_col, path)
        out[new_col] = pd.to_numeric(df[real_col], errors="coerce")

    out = out[out["species_normalized"] != "unknown_species"]
    out = collapse_numeric_by_species(out, list(OSKYRKO_COLUMNS.keys()))
    return out, audit


def standardize_title_traits(
    path: Path,
    focal_species: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = read_table_auto(path)
    species_col = require_column(
        df,
        ["treename", "species_normalized", "genus_species", "species", "species_name", "binomial"],
        "Title species",
        path,
    )

    audit = audit_raw_species(df[species_col], "title", focal_species)

    out = pd.DataFrame()
    out["species_normalized"] = df[species_col].map(normalize_species_name)

    for new_col, source_col in TITLE_COLUMNS.items():
        real_col = choose_first_column(df, [source_col])
        if real_col is None:
            out[new_col] = pd.NA
        else:
            out[new_col] = pd.to_numeric(df[real_col], errors="coerce")

    out = out[out["species_normalized"] != "unknown_species"]
    out = collapse_numeric_by_species(out, list(TITLE_COLUMNS.keys()))
    return out, audit


def standardize_ji_mass(
    path: Path,
    focal_species: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = read_table_auto(path)

    species_col = choose_first_column(
        df,
        [
            "species_normalized",
            "genus_species",
            "species_key",
            "species",
            "species_name",
            "binomial",
        ],
    )
    mass_col = choose_first_column(df, ["mass", "body_mass_g", "body mass (g)"])

    if species_col is None or mass_col is None:
        raise ValueError(
            f"Could not identify Ji species/mass columns in {path}. "
            f"Columns found: {list(df.columns)}"
        )

    audit = audit_raw_species(df[species_col], "ji", focal_species)

    out = pd.DataFrame()
    out["species_normalized"] = df[species_col].map(normalize_species_name)
    out["mass_grams_ji"] = pd.to_numeric(df[mass_col], errors="coerce")
    out = out[out["species_normalized"] != "unknown_species"]
    out = collapse_numeric_by_species(out, ["mass_grams_ji"])

    return out, audit


def metric_is_critical_thermal_limit(value: object, desired: str) -> bool:
    if pd.isna(value):
        return False

    text = str(value).strip().lower()
    text = text.replace("-", "").replace("_", "").replace(" ", "")

    if desired == "ctmax":
        return text in {"ctmax", "criticalthermalmaximum", "criticalthermalmax"}
    if desired == "ctmin":
        return text in {"ctmin", "criticalthermalminimum", "criticalthermalmin"}

    raise ValueError(f"Unknown desired metric: {desired}")


def standardize_bennett_thermal(
    path: Path,
    focal_species: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = read_table_auto(path)

    genus_col = choose_first_column(df, ["Genus", "genus"])
    species_col = choose_first_column(df, ["Species", "species"])
    tmax_col = choose_first_column(df, ["Tmax", "tmax"])
    max_metric_col = choose_first_column(df, ["max_metric"])
    tmin_col = choose_first_column(df, ["tmin", "Tmin"])
    min_metric_col = choose_first_column(df, ["min_metric"])

    required = {
        "Genus": genus_col,
        "Species": species_col,
        "Tmax": tmax_col,
        "max_metric": max_metric_col,
        "tmin": tmin_col,
        "min_metric": min_metric_col,
    }
    missing = [name for name, col in required.items() if col is None]
    if missing:
        raise ValueError(
            f"Missing required thermal columns in {path}: {missing}. "
            f"Columns found: {list(df.columns)}"
        )

    raw_species_names = (
        df[genus_col].astype("string").fillna("").str.strip()
        + " "
        + df[species_col].astype("string").fillna("").str.strip()
    ).str.strip()

    audit = audit_raw_species(raw_species_names, "bennett_thermal", focal_species)

    rows = []
    for idx, row in df.iterrows():
        source_name = raw_species_names.loc[idx]
        species_normalized = normalize_species_name(source_name)

        ct_max = pd.NA
        ct_min = pd.NA

        if metric_is_critical_thermal_limit(row[max_metric_col], "ctmax"):
            ct_max = pd.to_numeric(row[tmax_col], errors="coerce")

        if metric_is_critical_thermal_limit(row[min_metric_col], "ctmin"):
            ct_min = pd.to_numeric(row[tmin_col], errors="coerce")

        rows.append(
            {
                "species_normalized": species_normalized,
                "ct_max": ct_max,
                "ctmax_metric": row[max_metric_col],
                "ct_min": ct_min,
                "ctmin_metric": row[min_metric_col],
                "thermal_source": str(path),
            }
        )

    raw = pd.DataFrame(rows)
    raw = raw[raw["species_normalized"] != "unknown_species"]

    collapsed_rows = []
    for species_normalized, group in raw.groupby("species_normalized", sort=True):
        ctmax_series = group["ct_max"].dropna()
        ctmin_series = group["ct_min"].dropna()

        collapsed_rows.append(
            {
                "species_normalized": species_normalized,
                "ct_min": ctmin_series.iloc[0] if len(ctmin_series) else pd.NA,
                "ctmin_metric": (
                    group.loc[group["ct_min"].notna(), "ctmin_metric"].iloc[0]
                    if len(ctmin_series)
                    else ""
                ),
                "ct_max": ctmax_series.iloc[0] if len(ctmax_series) else pd.NA,
                "ctmax_metric": (
                    group.loc[group["ct_max"].notna(), "ctmax_metric"].iloc[0]
                    if len(ctmax_series)
                    else ""
                ),
                "thermal_source": str(path),
                "n_thermal_rows_for_species": len(group),
            }
        )

    return pd.DataFrame(
        collapsed_rows,
        columns=[
            "species_normalized",
            "ct_min",
            "ctmin_metric",
            "ct_max",
            "ctmax_metric",
            "thermal_source",
            "n_thermal_rows_for_species",
        ],
    ), audit


def blank_columns(df: pd.DataFrame, columns: list[str], fill_value: object = pd.NA) -> None:
    for col in columns:
        if col not in df.columns:
            df[col] = fill_value


def build_standardized_natural_history(
    project_species: pd.DataFrame,
    genome_size_path: Optional[Path],
    traits_oskyrko_path: Optional[Path],
    title_traits_path: Optional[Path],
    mass_ji_path: Optional[Path],
    thermal_bennett_path: Optional[Path],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Left-join all natural-history sources onto the focal project species list.
    """
    focal_species = set(project_species["species_normalized"].astype(str))

    out = project_species[["species_pk", "species_normalized", "genus_species"]].copy()
    audit_frames: list[pd.DataFrame] = []

    if genome_size_path is not None:
        genome_df, audit_df = standardize_genome_size(genome_size_path, focal_species)
        out = out.merge(genome_df, on="species_normalized", how="left")
        audit_frames.append(audit_df)
    else:
        blank_columns(
            out,
            [
                "mean_c_value_pg",
                "mean_genome_mb",
                "mean_genome_gb",
                "genome_size",
                "genome_size_source",
            ],
        )
        audit_frames.append(audit_missing_dataset("c_value", focal_species))

    if traits_oskyrko_path is not None:
        oskyrko_df, audit_df = standardize_oskyrko_traits(
            traits_oskyrko_path,
            focal_species,
        )
        out = out.merge(oskyrko_df, on="species_normalized", how="left")
        audit_frames.append(audit_df)
    else:
        blank_columns(out, list(OSKYRKO_COLUMNS.keys()))
        audit_frames.append(audit_missing_dataset("oskyrko", focal_species))

    if title_traits_path is not None:
        title_df, audit_df = standardize_title_traits(
            title_traits_path,
            focal_species,
        )
        out = out.merge(title_df, on="species_normalized", how="left")
        audit_frames.append(audit_df)
    else:
        blank_columns(out, list(TITLE_COLUMNS.keys()))
        audit_frames.append(audit_missing_dataset("title", focal_species))

    if mass_ji_path is not None:
        ji_df, audit_df = standardize_ji_mass(mass_ji_path, focal_species)
        out = out.merge(ji_df, on="species_normalized", how="left")
        audit_frames.append(audit_df)
    else:
        blank_columns(out, ["mass_grams_ji"])
        audit_frames.append(audit_missing_dataset("ji", focal_species))

    out["mass_preferred"] = (
        out["mass_grams_oskyrko"]
        .combine_first(out["mass_grams_title"])
        .combine_first(out["mass_grams_ji"])
    )

    def preferred_mass_source(row: pd.Series) -> str:
        if pd.notna(row["mass_grams_oskyrko"]):
            return "oskyrko"
        if pd.notna(row["mass_grams_title"]):
            return "title"
        if pd.notna(row["mass_grams_ji"]):
            return "ji"
        return ""

    out["mass_source_preferred"] = out.apply(preferred_mass_source, axis=1)
    out["annual_number_of_eggs"] = (
        pd.to_numeric(out["mean_number_of_eggs_per_clutch"], errors="coerce")
        * pd.to_numeric(out["number_clutches_per_year"], errors="coerce")
    )

    if thermal_bennett_path is not None:
        thermal_df, audit_df = standardize_bennett_thermal(
            thermal_bennett_path,
            focal_species,
        )
        out = out.merge(thermal_df, on="species_normalized", how="left")
        audit_frames.append(audit_df)
    else:
        blank_columns(
            out,
            [
                "ct_min",
                "ctmin_metric",
                "ct_max",
                "ctmax_metric",
                "thermal_source",
                "n_thermal_rows_for_species",
            ],
        )
        audit_frames.append(audit_missing_dataset("bennett_thermal", focal_species))

    # Guarantee exact schema/order requested for downstream use.
    blank_columns(out, STANDARDIZED_NATURAL_HISTORY_COLUMNS)
    out = out[STANDARDIZED_NATURAL_HISTORY_COLUMNS]

    audit = pd.concat(audit_frames, ignore_index=True)
    audit = audit[BASE_AUDIT_COLUMNS]

    return out, audit


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create natural_history_standardized.tsv and species_name_audit.tsv "
            "for the GC3/body-size database workflow."
        )
    )

    parser.add_argument(
        "genomes_dir",
        type=Path,
        help="Path to the top-level genomes directory.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="Project manifest CSV/TSV. Supersedes Compleasm metadata.",
    )
    parser.add_argument(
        "--genome-size",
        type=Path,
        default=None,
        help="C-value summary CSV/TSV with mean_c_value_pg.",
    )
    parser.add_argument(
        "--traits-oskyrko",
        type=Path,
        default=None,
        help="Oskyrko et al. 2024 ReptTraits CSV/TSV.",
    )
    parser.add_argument(
        "--mass-title",
        type=Path,
        default=None,
        help="Title et al. macroevolutionary dataset CSV/TSV.",
    )
    parser.add_argument(
        "--mass-ji",
        type=Path,
        default=None,
        help="Ji et al. mass CSV/TSV with species and mass columns.",
    )
    parser.add_argument(
        "--thermal-bennett",
        type=Path,
        default=None,
        help="Bennett/GlobTherm-style thermal limits CSV/TSV.",
    )
    parser.add_argument(
        "--outdir-traits",
        type=Path,
        default=None,
        help=(
            "Directory for natural_history_standardized.tsv. "
            "Default: <genomes_dir>/records/natural_history_standardized/"
        ),
    )
    parser.add_argument(
        "--outdir-audit",
        type=Path,
        default=None,
        help=(
            "Directory for species_name_audit.tsv. "
            "Default: <genomes_dir>/records/sql_tsvs/"
        ),
    )
    parser.add_argument(
        "--allow-missing-natural-history",
        action="store_true",
        help=(
            "Write blank trait values and missing_trait audit rows instead of "
            "failing when natural-history source files are unavailable."
        ),
    )
    parser.add_argument(
        "--species-pk-map",
        type=Path,
        default=None,
        help=(
            "Optional upstream CSV/TSV containing species_pk plus a "
            "species-normalized/name column. 00a also reuses species_pk from "
            "the focal-project input or a previous natural_history_standardized.tsv. "
            "It never reads downstream genomes.tsv/natural_history.tsv."
        ),
    )

    return parser.parse_args()


def resolve_optional_input(path: Optional[Path]) -> Optional[Path]:
    if path is None:
        return None
    resolved = path.expanduser().resolve()
    if not resolved.exists():
        raise FileNotFoundError(f"Input file does not exist: {resolved}")
    return resolved


def main() -> None:
    args = parse_args()

    genomes_dir = args.genomes_dir.expanduser().resolve()
    if not genomes_dir.exists():
        raise FileNotFoundError(f"genomes_dir does not exist: {genomes_dir}")

    outdir_traits = (
        args.outdir_traits.expanduser().resolve()
        if args.outdir_traits
        else genomes_dir / "records" / "natural_history_standardized"
    )
    outdir_audit = (
        args.outdir_audit.expanduser().resolve()
        if args.outdir_audit
        else genomes_dir / "records" / "sql_tsvs"
    )

    outdir_traits.mkdir(parents=True, exist_ok=True)
    outdir_audit.mkdir(parents=True, exist_ok=True)

    nh_dir = natural_history_dir(genomes_dir)

    if args.manifest:
        project_source_path = args.manifest.expanduser().resolve()
        project_source_label = "manifest"
    else:
        project_source_path = default_compleasm_metadata(genomes_dir)
        project_source_label = "compleasm_metadata"

    if not project_source_path.exists():
        raise FileNotFoundError(
            f"Project species source does not exist: {project_source_path}"
        )

    project_species = load_project_species(
        project_source_path,
        project_source_label,
    )
    print(
        f"[OK] Loaded {len(project_species)} focal species from "
        f"{project_source_path} ({project_source_label})"
    )

    project_species, species_pk_sources = resolve_canonical_species_pks(
        project_species=project_species,
        project_source_path=project_source_path,
        existing_natural_history_path=(
            outdir_traits / "natural_history_standardized.tsv"
        ),
        explicit_species_pk_map=args.species_pk_map,
    )

    genome_size_path = (
        resolve_optional_input(args.genome_size)
        if args.genome_size
        else find_default_genome_size_file(nh_dir)
    )
    traits_oskyrko_path = (
        resolve_optional_input(args.traits_oskyrko)
        if args.traits_oskyrko
        else find_default_file(nh_dir, ["oskyrko"])
    )
    title_traits_path = (
        resolve_optional_input(args.mass_title)
        if args.mass_title
        else find_default_file(nh_dir, ["title"])
    )
    mass_ji_path = (
        resolve_optional_input(args.mass_ji)
        if args.mass_ji
        else find_default_file(nh_dir, ["ji"])
    )
    thermal_bennett_path = (
        resolve_optional_input(args.thermal_bennett)
        if args.thermal_bennett
        else find_default_file(nh_dir, ["bennett"])
    )

    if genome_size_path is None and not args.allow_missing_natural_history:
        raise FileNotFoundError(
            "No genome-size/C-value summary file was provided or auto-detected."
        )
    if (
        traits_oskyrko_path is None
        and title_traits_path is None
        and mass_ji_path is None
        and not args.allow_missing_natural_history
    ):
        raise FileNotFoundError(
            "No Oskyrko/ReptTraits, Title, or Ji trait/mass file was "
            "provided or auto-detected."
        )
    if thermal_bennett_path is None and not args.allow_missing_natural_history:
        raise FileNotFoundError(
            "No Bennett/GlobTherm thermal file was provided or auto-detected."
        )

    standardized_df, audit_df = build_standardized_natural_history(
        project_species=project_species,
        genome_size_path=genome_size_path,
        traits_oskyrko_path=traits_oskyrko_path,
        title_traits_path=title_traits_path,
        mass_ji_path=mass_ji_path,
        thermal_bennett_path=thermal_bennett_path,
    )

    # v6: 00a owns the canonical biological species_pk mapping carried by
    # natural_history_standardized.tsv. The audit reuses that mapping, and
    # species_name_audit_pk is assigned only after all source-specific audit
    # frames have been concatenated.
    audit_df = finalize_species_name_audit(audit_df, project_species)

    standardized_out = outdir_traits / "natural_history_standardized.tsv"
    audit_out = outdir_audit / "species_name_audit.tsv"

    write_tsv(standardized_df, standardized_out)
    write_species_name_audit_tsv(audit_df, audit_out)

    print(
        f"[OK] Wrote {standardized_out} "
        f"({len(standardized_df)} focal species; "
        f"genome_size={int(standardized_df['genome_size'].notna().sum())}, "
        f"mass_preferred={int(standardized_df['mass_preferred'].notna().sum())}, "
        f"ct_max={int(standardized_df['ct_max'].notna().sum())}, "
        f"ct_min={int(standardized_df['ct_min'].notna().sum())})"
    )

    status_counts = audit_df["match_status"].value_counts()
    print(
        f"[OK] Wrote {audit_out} "
        f"(rows={len(audit_df):,}, "
        f"matched={int(status_counts.get('matched', 0))}, "
        f"unmatched_raw_trait={int(status_counts.get('unmatched_raw_trait', 0))}, "
        f"missing_trait={int(status_counts.get('missing_trait', 0))})"
    )
    print("[INFO] species_name_audit.tsv SQL NULL representation: \\\\N")
    print(
        "[INFO] 00a owns canonical species_pk assignment; "
        "natural_history_standardized.tsv is the authoritative handoff to 00b"
    )
    if species_pk_sources:
        print(
            "[INFO] reused species_pk source(s): "
            + "; ".join(str(path) for path in species_pk_sources)
        )

    print("[DONE] Natural-history standardization and species-name audit complete.")
    print("[INFO] Canonical join key: species_normalized")
    print("[INFO] Preferred mass priority: oskyrko > title > ji")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[ERROR] {exc}")
        raise
