#!/usr/bin/env python3
"""
03_make_master_input_tsvs.py

Create standardized TSV input files for the downstream GC divergence workflow.

Canonical species key:
  species_normalized = genus_species in lowercase, e.g. crotalus_viridis

This version can stage natural-history predictors from:
  - genome-size/C-value summaries
  - Oskyrko et al. 2024 ReptTraits (--traits-oskyrko)
  - Title et al. macroevolutionary dataset (--mass-title / title_et_al_data.csv)
  - Ji-style mass data (--mass-ji)
  - Bennett/GlobTherm thermal limits

Outputs, by default to <genomes_dir>/records/gc_metrics_inputs/:
  genome_sizes.tsv
  species_mass.tsv
  thermal_limits.tsv

Each output includes species_normalized as the canonical join key.
For temporary backward compatibility, genus_species is also included as an alias.
"""

import argparse
import csv
import re
from pathlib import Path
from typing import Optional

import pandas as pd


BAD_SPECIES_TOKENS = {"sp", "sp.", "spp", "spp.", "cf", "cf.", "aff", "aff.", "nr", "nr."}
HYBRID_MARKERS = {"x", "×"}


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
    Normalize a source species name to a canonical Genus_species key.

    Rules:
      - keep only the first two valid taxonomic tokens
      - convert spaces to underscores
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


def add_genus_species_alias(df: pd.DataFrame) -> pd.DataFrame:
    if "species_normalized" in df.columns and "genus_species" not in df.columns:
        df.insert(1, "genus_species", df["species_normalized"])
    return df


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
        raise ValueError(f"Could not identify {label} column in {path}. Columns found: {list(df.columns)}")
    return col


def default_compleasm_metadata(genomes_dir: Path) -> Path:
    return genomes_dir / "records" / "compleasm" / "records" / "metadata.csv"


def natural_history_dir(genomes_dir: Path) -> Path:
    return genomes_dir / "records" / "natural_history"


def find_default_file(directory: Path, include_terms: list[str], required_columns: Optional[list[str]] = None) -> Optional[Path]:
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
        candidate = find_default_file(directory, terms, required_columns=["mean_c_value_pg"])
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
        ["accession", "accession_root", "assembly_accession", "assemblyAccession",
         "Assembly Accession", "genome_accession", "ncbi_accession"],
    )
    species_key_col = choose_first_column(
        df,
        ["species_normalized", "species_key", "genus_species", "genus_species_key", "normalized_species", "normalized_name"],
    )
    organism_col = choose_first_column(
        df,
        ["organism_name", "organismName", "Organism Name", "species",
         "species_name", "binomial", "binomial_2020"],
    )

    if species_key_col is None and organism_col is None:
        raise ValueError(
            f"Could not find species key or organism/species column in {source_path}. "
            f"Columns found: {list(df.columns)}"
        )

    out = pd.DataFrame()
    source_col = species_key_col if species_key_col is not None else organism_col
    out["species_normalized"] = df[source_col].map(normalize_species_name)

    if organism_col is not None:
        out["organism_name"] = df[organism_col].astype(str).str.strip()
    else:
        out["organism_name"] = out["species_normalized"].str.replace("_", " ", regex=False)

    if accession_col is not None:
        out["accession"] = df[accession_col].astype(str).str.strip()
    else:
        out["accession"] = pd.NA

    out["project_species_source"] = source_label
    out = out[out["species_normalized"].astype(str).str.len() > 0]
    out = out[out["species_normalized"] != "unknown_species"]
    out = out.drop_duplicates(subset=["species_normalized"], keep="first")
    return add_genus_species_alias(out)


def standardize_genome_size(path: Path, project_species: pd.DataFrame) -> pd.DataFrame:
    df = read_table_auto(path)

    cvalue_col = choose_first_column(df, ["mean_c_value_pg"])
    species_key_col = choose_first_column(df, ["species_normalized", "species_key", "genus_species"])
    species_name_col = choose_first_column(df, ["species_name", "organism_name", "organismName"])
    mean_mb_col = choose_first_column(df, ["mean_genome_mb"])
    mean_gb_col = choose_first_column(df, ["mean_genome_gb"])

    if cvalue_col is None:
        raise ValueError(f"Could not find mean_c_value_pg in {path}. Columns found: {list(df.columns)}")
    if species_key_col is None and species_name_col is None:
        raise ValueError(f"Could not find species key or species_name in {path}. Columns found: {list(df.columns)}")

    rows = []
    for _, row in df.iterrows():
        source_name = row[species_key_col] if species_key_col else row[species_name_col]
        species_normalized = normalize_species_name(source_name)
        rows.append(
            {
                "species_normalized": species_normalized,
                "species_name": row[species_name_col] if species_name_col else species_normalized.replace("_", " "),
                "mean_c_value_pg": pd.to_numeric(row[cvalue_col], errors="coerce"),
                "mean_genome_mb": pd.to_numeric(row[mean_mb_col], errors="coerce") if mean_mb_col else pd.NA,
                "mean_genome_gb": pd.to_numeric(row[mean_gb_col], errors="coerce") if mean_gb_col else pd.NA,
                "genome_size_source": str(path),
            }
        )

    size_df = pd.DataFrame(rows)
    size_df = size_df.dropna(subset=["species_normalized", "mean_c_value_pg"])
    size_df = size_df[size_df["species_normalized"] != "unknown_species"]
    size_df = size_df.drop_duplicates(subset=["species_normalized"], keep="first")

    out = project_species[["species_normalized", "organism_name"]].copy()
    out = out.merge(size_df, on="species_normalized", how="left")
    out["species_name"] = out["species_name"].combine_first(out["organism_name"])
    out["genome_size"] = out["mean_c_value_pg"]
    out = add_genus_species_alias(out)
    return out[["species_normalized", "genus_species", "organism_name", "species_name", "mean_c_value_pg",
                "mean_genome_mb", "mean_genome_gb", "genome_size", "genome_size_source"]]


def collapse_numeric_by_species(df: pd.DataFrame, value_cols: list[str]) -> pd.DataFrame:
    """Collapse duplicate species rows by taking the first non-null value per column."""
    rows = []
    for species_normalized, group in df.groupby("species_normalized", sort=True):
        row = {"species_normalized": species_normalized}
        for col in value_cols:
            nonnull = group[col].dropna()
            row[col] = nonnull.iloc[0] if len(nonnull) else pd.NA
        rows.append(row)
    return pd.DataFrame(rows)


def standardize_oskyrko_traits(path: Path) -> pd.DataFrame:
    df = read_table_auto(path)
    species_col = require_column(
        df,
        ["Species", "species", "species_normalized", "genus_species", "species_name", "binomial"],
        "Oskyrko/ReptTraits species",
        path,
    )

    missing = [source_col for source_col in OSKYRKO_COLUMNS.values() if choose_first_column(df, [source_col]) is None]
    if missing:
        raise ValueError(f"Missing required Oskyrko/ReptTraits columns in {path}: {missing}. Columns found: {list(df.columns)}")

    out = pd.DataFrame()
    out["species_normalized"] = df[species_col].map(normalize_species_name)
    for new_col, source_col in OSKYRKO_COLUMNS.items():
        real_col = require_column(df, [source_col], source_col, path)
        out[new_col] = pd.to_numeric(df[real_col], errors="coerce")

    out = out[out["species_normalized"] != "unknown_species"]
    return collapse_numeric_by_species(out, list(OSKYRKO_COLUMNS.keys()))


def standardize_title_traits(path: Path) -> pd.DataFrame:
    df = read_table_auto(path)
    species_col = require_column(
        df,
        ["treename", "species_normalized", "genus_species", "species", "species_name", "binomial"],
        "Title species",
        path,
    )

    out = pd.DataFrame()
    out["species_normalized"] = df[species_col].map(normalize_species_name)
    for new_col, source_col in TITLE_COLUMNS.items():
        real_col = choose_first_column(df, [source_col])
        if real_col is None:
            out[new_col] = pd.NA
        else:
            out[new_col] = pd.to_numeric(df[real_col], errors="coerce")

    out = out[out["species_normalized"] != "unknown_species"]
    return collapse_numeric_by_species(out, list(TITLE_COLUMNS.keys()))


def standardize_ji_mass(path: Path) -> pd.DataFrame:
    df = read_table_auto(path)
    species_col = choose_first_column(df, ["species_normalized", "genus_species", "species_key", "species", "species_name", "binomial"])
    mass_col = choose_first_column(df, ["mass", "body_mass_g", "body mass (g)"])
    if species_col is None or mass_col is None:
        raise ValueError(f"Could not identify Ji species/mass columns in {path}. Columns found: {list(df.columns)}")
    out = pd.DataFrame()
    out["species_normalized"] = df[species_col].map(normalize_species_name)
    out["mass_grams_ji"] = pd.to_numeric(df[mass_col], errors="coerce")
    out = out.dropna(subset=["species_normalized", "mass_grams_ji"])
    out = out[out["species_normalized"] != "unknown_species"]
    return out.drop_duplicates(subset=["species_normalized"], keep="first")


def build_species_mass_tsv(
    project_species: pd.DataFrame,
    traits_oskyrko_path: Optional[Path],
    title_traits_path: Optional[Path],
    mass_ji_path: Optional[Path],
) -> pd.DataFrame:
    out = project_species[["species_normalized"]].copy()

    if traits_oskyrko_path:
        out = out.merge(standardize_oskyrko_traits(traits_oskyrko_path), on="species_normalized", how="left")
    else:
        for col in OSKYRKO_COLUMNS:
            out[col] = pd.NA

    if title_traits_path:
        out = out.merge(standardize_title_traits(title_traits_path), on="species_normalized", how="left")
    else:
        for col in TITLE_COLUMNS:
            out[col] = pd.NA

    if mass_ji_path:
        out = out.merge(standardize_ji_mass(mass_ji_path), on="species_normalized", how="left")
    else:
        out["mass_grams_ji"] = pd.NA

    out["mass_preferred"] = (
        out["mass_grams_oskyrko"]
        .combine_first(out["mass_grams_title"])
        .combine_first(out["mass_grams_ji"])
    )

    def source(row: pd.Series) -> str:
        if pd.notna(row["mass_grams_oskyrko"]):
            return "oskyrko"
        if pd.notna(row["mass_grams_title"]):
            return "title"
        if pd.notna(row["mass_grams_ji"]):
            return "ji"
        return ""

    out["mass_source_preferred"] = out.apply(source, axis=1)
    out["annual_number_of_eggs"] = out["mean_number_of_eggs_per_clutch"] * out["number_clutches_per_year"]

    out = add_genus_species_alias(out)
    return out[
        [
            "species_normalized",
            "genus_species",
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
        ]
    ]


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


def standardize_bennett_thermal(path: Path) -> pd.DataFrame:
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
        raise ValueError(f"Missing required thermal columns in {path}: {missing}. Columns found: {list(df.columns)}")

    rows = []
    for _, row in df.iterrows():
        species_normalized = normalize_species_name(f"{row[genus_col]} {row[species_col]}")

        ctmax = pd.NA
        ctmin = pd.NA

        if metric_is_critical_thermal_limit(row[max_metric_col], "ctmax"):
            ctmax = pd.to_numeric(row[tmax_col], errors="coerce")

        if metric_is_critical_thermal_limit(row[min_metric_col], "ctmin"):
            ctmin = pd.to_numeric(row[tmin_col], errors="coerce")

        rows.append(
            {
                "species_normalized": species_normalized,
                "ctmax": ctmax,
                "ctmax_metric": row[max_metric_col],
                "ctmin": ctmin,
                "ctmin_metric": row[min_metric_col],
                "thermal_source": str(path),
            }
        )

    raw = pd.DataFrame(rows)
    raw = raw[raw["species_normalized"] != "unknown_species"]

    collapsed_rows = []
    for species_normalized, group in raw.groupby("species_normalized", sort=True):
        ctmax_series = group["ctmax"].dropna()
        ctmin_series = group["ctmin"].dropna()

        collapsed_rows.append(
            {
                "species_normalized": species_normalized,
                "ctmax": ctmax_series.iloc[0] if len(ctmax_series) else pd.NA,
                "ctmax_metric": group.loc[group["ctmax"].notna(), "ctmax_metric"].iloc[0] if len(ctmax_series) else "",
                "ctmin": ctmin_series.iloc[0] if len(ctmin_series) else pd.NA,
                "ctmin_metric": group.loc[group["ctmin"].notna(), "ctmin_metric"].iloc[0] if len(ctmin_series) else "",
                "thermal_source": str(path),
                "n_thermal_rows_for_species": len(group),
            }
        )

    return pd.DataFrame(collapsed_rows)


def build_thermal_limits_tsv(project_species: pd.DataFrame, thermal_bennett_path: Optional[Path]) -> pd.DataFrame:
    out = project_species[["species_normalized"]].copy()
    if thermal_bennett_path:
        out = out.merge(standardize_bennett_thermal(thermal_bennett_path), on="species_normalized", how="left")
    else:
        out["ctmax"] = pd.NA
        out["ctmax_metric"] = ""
        out["ctmin"] = pd.NA
        out["ctmin_metric"] = ""
        out["thermal_source"] = ""
        out["n_thermal_rows_for_species"] = pd.NA

    out = add_genus_species_alias(out)
    return out[["species_normalized", "genus_species", "ctmax", "ctmax_metric", "ctmin", "ctmin_metric", "thermal_source", "n_thermal_rows_for_species"]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create standardized natural-history TSV inputs for GC metric master export.")
    parser.add_argument("genomes_dir", type=Path, help="Path to the top-level genomes directory.")
    parser.add_argument("--manifest", type=Path, default=None, help="Project manifest CSV/TSV. Supersedes Compleasm metadata.")
    parser.add_argument("--genome-size", type=Path, default=None, help="C-value summary CSV/TSV with mean_c_value_pg.")
    parser.add_argument("--traits-oskyrko", type=Path, default=None, help="Oskyrko et al. 2024 ReptTraits CSV/TSV.")
    parser.add_argument("--mass-title", type=Path, default=None, help="Title et al. macroevolutionary dataset CSV/TSV, e.g. title_et_al_data.csv.")
    parser.add_argument("--mass-ji", type=Path, default=None, help="Ji et al. mass TSV with species and mass columns.")
    parser.add_argument("--thermal-bennett", type=Path, default=None, help="Bennett/GlobTherm-style thermal limits CSV/TSV.")
    parser.add_argument("--outdir", type=Path, default=None, help="Output directory. Default: <genomes_dir>/records/gc_metrics_inputs")
    parser.add_argument("--allow-missing-natural-history", action="store_true", help="Write blank values instead of failing if natural-history files are missing.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    genomes_dir = args.genomes_dir.expanduser().resolve()
    if not genomes_dir.exists():
        raise FileNotFoundError(f"genomes_dir does not exist: {genomes_dir}")

    outdir = args.outdir.expanduser().resolve() if args.outdir else genomes_dir / "records" / "gc_metrics_inputs"
    outdir.mkdir(parents=True, exist_ok=True)

    nh_dir = natural_history_dir(genomes_dir)

    if args.manifest:
        project_source_path = args.manifest.expanduser().resolve()
        project_source_label = "manifest"
    else:
        project_source_path = default_compleasm_metadata(genomes_dir)
        project_source_label = "compleasm_metadata"

    if not project_source_path.exists():
        raise FileNotFoundError(f"Project species source does not exist: {project_source_path}")

    project_species = load_project_species(project_source_path, project_source_label)
    print(f"[OK] Loaded {len(project_species)} focal species from {project_source_path} ({project_source_label})")

    genome_size_path = args.genome_size.expanduser().resolve() if args.genome_size else find_default_genome_size_file(nh_dir)
    traits_oskyrko_path = args.traits_oskyrko.expanduser().resolve() if args.traits_oskyrko else find_default_file(nh_dir, ["oskyrko"])
    title_traits_path = args.mass_title.expanduser().resolve() if args.mass_title else find_default_file(nh_dir, ["title"])
    mass_ji_path = args.mass_ji.expanduser().resolve() if args.mass_ji else find_default_file(nh_dir, ["ji"])
    thermal_bennett_path = args.thermal_bennett.expanduser().resolve() if args.thermal_bennett else find_default_file(nh_dir, ["bennett"])

    if genome_size_path is None and not args.allow_missing_natural_history:
        raise FileNotFoundError("No genome-size/C-value summary file was provided or auto-detected.")
    if traits_oskyrko_path is None and title_traits_path is None and mass_ji_path is None and not args.allow_missing_natural_history:
        raise FileNotFoundError("No Oskyrko/ReptTraits, Title, or Ji trait/mass file was provided or auto-detected.")
    if thermal_bennett_path is None and not args.allow_missing_natural_history:
        raise FileNotFoundError("No Bennett/GlobTherm thermal file was provided or auto-detected.")

    if genome_size_path:
        genome_size_df = standardize_genome_size(genome_size_path, project_species)
    else:
        genome_size_df = project_species[["species_normalized", "genus_species", "organism_name"]].copy()
        genome_size_df["species_name"] = genome_size_df["organism_name"]
        genome_size_df["mean_c_value_pg"] = pd.NA
        genome_size_df["mean_genome_mb"] = pd.NA
        genome_size_df["mean_genome_gb"] = pd.NA
        genome_size_df["genome_size"] = pd.NA
        genome_size_df["genome_size_source"] = ""

    genome_size_out = outdir / "genome_sizes.tsv"
    write_tsv(genome_size_df, genome_size_out)
    print(f"[OK] Wrote {genome_size_out} ({int(genome_size_df['genome_size'].notna().sum())}/{len(genome_size_df)} species matched)")

    mass_df = build_species_mass_tsv(project_species, traits_oskyrko_path, title_traits_path, mass_ji_path)
    mass_out = outdir / "species_mass.tsv"
    write_tsv(mass_df, mass_out)
    print(
        f"[OK] Wrote {mass_out} "
        f"(oskyrko_mass={int(mass_df['mass_grams_oskyrko'].notna().sum())}, "
        f"title_mass={int(mass_df['mass_grams_title'].notna().sum())}, "
        f"ji_mass={int(mass_df['mass_grams_ji'].notna().sum())}, "
        f"preferred={int(mass_df['mass_preferred'].notna().sum())}/{len(mass_df)}; "
        "preferred priority=oskyrko>title>ji)"
    )
    print(
        f"[INFO] Trait matches: svl_oskyrko={int(mass_df['max_female_length_svl_mm_oskyrko'].notna().sum())}, "
        f"svl_title={int(mass_df['svl_mm_title'].notna().sum())}, "
        f"annual_eggs={int(mass_df['annual_number_of_eggs'].notna().sum())}, "
        f"range_size={int(mass_df['range_size'].notna().sum())}"
    )

    thermal_df = build_thermal_limits_tsv(project_species, thermal_bennett_path)
    thermal_out = outdir / "thermal_limits.tsv"
    write_tsv(thermal_df, thermal_out)
    print(
        f"[OK] Wrote {thermal_out} "
        f"(ctmax={int(thermal_df['ctmax'].notna().sum())}/{len(thermal_df)}, "
        f"ctmin={int(thermal_df['ctmin'].notna().sum())}/{len(thermal_df)})"
    )

    print("[DONE] Natural-history TSV input staging complete.")
    print(f"[INFO] Output directory: {outdir}")
    print("[INFO] Canonical join key: species_normalized")


if __name__ == "__main__":
    main()
