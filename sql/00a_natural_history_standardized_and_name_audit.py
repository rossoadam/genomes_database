#!/usr/bin/env python3
"""
00a_natural_history_standardized_and_name_audit.py

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
"""

import argparse
import csv
import re
from pathlib import Path
from typing import Optional

import pandas as pd


BAD_SPECIES_TOKENS = {"sp", "sp.", "spp", "spp.", "cf", "cf.", "aff", "aff.", "nr", "nr."}
HYBRID_MARKERS = {"x", "×"}

AUDIT_COLUMNS = [
    "source_dataset",
    "source_species_name",
    "species_normalized",
    "match_status",
]

STANDARDIZED_NATURAL_HISTORY_COLUMNS = [
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

    if species_key_col is None and organism_col is None:
        raise ValueError(
            f"Could not find species key or organism/species column in {source_path}. "
            f"Columns found: {list(df.columns)}"
        )

    source_col = species_key_col if species_key_col is not None else organism_col

    out = pd.DataFrame()
    out["species_normalized"] = df[source_col].map(normalize_species_name)

    if organism_col is not None:
        out["organism_name"] = df[organism_col].astype("string").str.strip()
    else:
        out["organism_name"] = out["species_normalized"].str.replace("_", " ", regex=False)

    if accession_col is not None:
        out["accession"] = df[accession_col].astype("string").str.strip()
    else:
        out["accession"] = pd.NA

    out["project_species_source"] = source_label
    out = out[out["species_normalized"].astype(str).str.len() > 0]
    out = out[out["species_normalized"] != "unknown_species"]
    out = out.drop_duplicates(subset=["species_normalized"], keep="first")
    out["genus_species"] = out["species_normalized"]

    return out[
        [
            "species_normalized",
            "genus_species",
            "organism_name",
            "accession",
            "project_species_source",
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

    return audit[AUDIT_COLUMNS]


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
        columns=AUDIT_COLUMNS,
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

    out = project_species[["species_normalized", "genus_species"]].copy()
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
    audit = audit[AUDIT_COLUMNS]

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

    standardized_out = outdir_traits / "natural_history_standardized.tsv"
    audit_out = outdir_audit / "species_name_audit.tsv"

    write_tsv(standardized_df, standardized_out)
    write_tsv(audit_df, audit_out)

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
        f"(matched={int(status_counts.get('matched', 0))}, "
        f"unmatched_raw_trait={int(status_counts.get('unmatched_raw_trait', 0))}, "
        f"missing_trait={int(status_counts.get('missing_trait', 0))})"
    )

    print("[DONE] Natural-history standardization and species-name audit complete.")
    print("[INFO] Canonical join key: species_normalized")
    print("[INFO] Preferred mass priority: oskyrko > title > ji")


if __name__ == "__main__":
    main()
