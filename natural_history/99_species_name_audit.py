#!/usr/bin/env python3
"""
99_species_name_audit.py

Create SQL-ready natural_history.tsv and species_name_audit.tsv files from
natural-history source datasets and optional genomes metadata/project manifest.

This script is intended to live in the genomes_database/natural_history workflow
and write final import tables to:

  <genomes_dir>/records/sql_tsvs/natural_history.tsv
  <genomes_dir>/records/sql_tsvs/species_name_audit.tsv

Recommended DBML schema:

Table natural_history {
  species_pk bigint [primary key]
  species_normalized varchar [unique]
  mass_meiri double
  mass_title double
  mass_ji double
  genome_size double
  ct_min double
  ct_max double
}

Table species_name_audit {
  species_name_audit_pk bigint [primary key]
  species_pk bigint [ref: > natural_history.species_pk]
  source_dataset varchar
  source_species_name varchar
  species_normalized varchar
  match_status varchar
}

Typical usage:

  python 99_species_name_audit.py /path/to/genomes \
    --genomes-metadata /path/to/genomes/records/genomes_metadata.csv \
    --genome-size /path/to/genomes/records/natural_history/reptile_c_value_summary.csv \
    --mass-meiri /path/to/genomes/records/natural_history/meiri.csv \
    --mass-title /path/to/genomes/records/natural_history/title.csv \
    --mass-ji /path/to/genomes/records/natural_history/ji.tsv \
    --thermal-bennett /path/to/genomes/records/natural_history/bennett.csv

If file paths are omitted, the script tries to auto-detect common files in:

  <genomes_dir>/records/natural_history/
  <genomes_dir>/records/genomes_metadata.csv
  <genomes_dir>/records/project_manifests/*.csv

Notes:
  - species_pk and species_name_audit_pk are assigned deterministically by
    alphabetical species_normalized order and audit row order.
  - match_status is evaluated relative to the final natural_history species set.
  - source-specific species names are preserved in species_name_audit.tsv.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd


BAD_TOKENS = {"sp", "sp.", "spp", "spp.", "cf", "cf.", "aff", "aff.", "nr", "nr."}
HYBRID_MARKERS = {"x", "×"}


def normalize_species_name(value: object) -> str:
    """
    Convert a raw species name into the canonical SQL key.

    Examples:
      Crotalus viridis viridis -> crotalus_viridis
      Crotalus_viridis -> crotalus_viridis
      Crotalus cf. viridis -> crotalus_viridis
    """
    if value is None or pd.isna(value):
        return "unknown_species"

    text = str(value).strip()
    if not text:
        return "unknown_species"

    text = text.replace("_", " ")
    text = re.sub(r"\([^)]*\)", " ", text)
    text = re.sub(r"[^A-Za-z×x\s.\-]", " ", text)
    text = re.sub(r"\s+", " ", text).strip()

    parts: list[str] = []
    for token in text.split():
        clean = token.strip().strip(".").lower()
        if not clean:
            continue
        if clean in BAD_TOKENS:
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
    if suffixes.endswith((".tsv", ".txt")):
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


def choose_first_column(df: pd.DataFrame, candidates: Iterable[str]) -> Optional[str]:
    cleaned_to_real = {str(c).strip().lower(): c for c in df.columns}
    for candidate in candidates:
        key = candidate.strip().lower()
        if key in cleaned_to_real:
            return cleaned_to_real[key]
    return None


def natural_history_dir(genomes_dir: Path) -> Path:
    return genomes_dir / "records" / "natural_history"


def default_genomes_metadata(genomes_dir: Path) -> Path:
    return genomes_dir / "records" / "genomes_metadata.csv"


def default_sql_tsv_dir(genomes_dir: Path) -> Path:
    return genomes_dir / "records" / "sql_tsvs"


def find_default_file(directory: Path, include_terms: list[str], required_columns: Optional[list[str]] = None) -> Optional[Path]:
    if not directory.exists():
        return None

    candidates: list[Path] = []
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


def source_label(path: Path, explicit_label: Optional[str] = None) -> str:
    if explicit_label:
        return explicit_label
    return path.stem


def make_audit_rows_from_species_column(
    df: pd.DataFrame,
    species_col: str,
    source_dataset: str,
) -> pd.DataFrame:
    rows = []
    for raw_name in df[species_col].dropna().astype(str):
        raw_name = raw_name.strip()
        if not raw_name:
            continue
        rows.append(
            {
                "source_dataset": source_dataset,
                "source_species_name": raw_name,
                "species_normalized": normalize_species_name(raw_name),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=["source_dataset", "source_species_name", "species_normalized"])
    out = out[out["species_normalized"] != "unknown_species"]
    return out.drop_duplicates(subset=["source_dataset", "source_species_name", "species_normalized"], keep="first")


def load_reference_species_from_file(path: Path, label: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read genomes_metadata, manifest, or similar project source and return:
      1. species list with species_normalized
      2. audit rows containing original source names
    """
    df = read_table_auto(path)
    species_key_col = choose_first_column(
        df,
        ["species_normalized", "species_key", "genus_species", "genus_species_key", "normalized_species", "normalized_name"],
    )
    organism_col = choose_first_column(
        df,
        ["organism_name", "organismName", "Organism Name", "species", "species_name", "binomial", "binomial_2020"],
    )

    if species_key_col is None and organism_col is None:
        raise ValueError(
            f"Could not find species key or organism/species column in {path}. Columns found: {list(df.columns)}"
        )

    key_source_col = species_key_col if species_key_col is not None else organism_col
    species = pd.DataFrame()
    species["species_normalized"] = df[key_source_col].map(normalize_species_name)
    species = species[species["species_normalized"] != "unknown_species"]
    species = species.drop_duplicates(subset=["species_normalized"], keep="first")

    audit_frames = []
    if species_key_col is not None:
        audit_frames.append(make_audit_rows_from_species_column(df, species_key_col, f"{label}:species_key"))
    if organism_col is not None and organism_col != species_key_col:
        audit_frames.append(make_audit_rows_from_species_column(df, organism_col, f"{label}:organism_name"))

    audit = pd.concat(audit_frames, ignore_index=True) if audit_frames else pd.DataFrame()
    return species, audit


def standardize_genome_size(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = read_table_auto(path)
    cvalue_col = choose_first_column(df, ["mean_c_value_pg", "genome_size", "c_value_pg"])
    species_key_col = choose_first_column(df, ["species_normalized", "species_key", "genus_species"])
    species_name_col = choose_first_column(df, ["species_name", "organism_name", "organismName", "binomial"])

    if cvalue_col is None:
        raise ValueError(f"Could not find genome-size column in {path}. Columns found: {list(df.columns)}")
    if species_key_col is None and species_name_col is None:
        raise ValueError(f"Could not find species key/name column in {path}. Columns found: {list(df.columns)}")

    name_col = species_key_col if species_key_col is not None else species_name_col
    out = pd.DataFrame()
    out["species_normalized"] = df[name_col].map(normalize_species_name)
    out["genome_size"] = pd.to_numeric(df[cvalue_col], errors="coerce")
    out = out.dropna(subset=["genome_size"])
    out = out[out["species_normalized"] != "unknown_species"]
    out = out.drop_duplicates(subset=["species_normalized"], keep="first")

    audit_frames = []
    if species_key_col is not None:
        audit_frames.append(make_audit_rows_from_species_column(df, species_key_col, f"{path.stem}:species_key"))
    if species_name_col is not None and species_name_col != species_key_col:
        audit_frames.append(make_audit_rows_from_species_column(df, species_name_col, f"{path.stem}:species_name"))
    audit = pd.concat(audit_frames, ignore_index=True) if audit_frames else pd.DataFrame()
    return out[["species_normalized", "genome_size"]], audit


def standardize_meiri_mass(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = read_table_auto(path)
    species_col = choose_first_column(df, ["binomial_2020", "binomial_(original files)", "species_normalized", "genus_species", "species", "species_name"])
    mass_col = choose_first_column(df, ["body mass (g)", "adult_body_mass (g)", "adult_body_mass_g", "mass", "body_mass_g"])
    if species_col is None or mass_col is None:
        raise ValueError(f"Could not identify Meiri species/mass columns in {path}. Columns found: {list(df.columns)}")

    out = pd.DataFrame()
    out["species_normalized"] = df[species_col].map(normalize_species_name)
    out["mass_meiri"] = pd.to_numeric(df[mass_col], errors="coerce")
    out = out.dropna(subset=["mass_meiri"])
    out = out[out["species_normalized"] != "unknown_species"]
    out = out.drop_duplicates(subset=["species_normalized"], keep="first")
    audit = make_audit_rows_from_species_column(df, species_col, path.stem)
    return out[["species_normalized", "mass_meiri"]], audit


def standardize_title_mass(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = read_table_auto(path)
    species_col = choose_first_column(df, ["treename", "species_normalized", "genus_species", "species", "species_name", "binomial"])
    mass_col = choose_first_column(df, ["mass", "body_mass_g", "body mass (g)"])
    if species_col is None or mass_col is None:
        raise ValueError(f"Could not identify Title species/mass columns in {path}. Columns found: {list(df.columns)}")

    out = pd.DataFrame()
    out["species_normalized"] = df[species_col].map(normalize_species_name)
    out["mass_title"] = pd.to_numeric(df[mass_col], errors="coerce")
    out = out.dropna(subset=["mass_title"])
    out = out[out["species_normalized"] != "unknown_species"]
    out = out.drop_duplicates(subset=["species_normalized"], keep="first")
    audit = make_audit_rows_from_species_column(df, species_col, path.stem)
    return out[["species_normalized", "mass_title"]], audit


def standardize_ji_mass(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = read_table_auto(path)
    species_col = choose_first_column(df, ["species_normalized", "genus_species", "species_key", "species", "species_name", "binomial"])
    mass_col = choose_first_column(df, ["mass", "body_mass_g", "body mass (g)"])
    if species_col is None or mass_col is None:
        raise ValueError(f"Could not identify Ji species/mass columns in {path}. Columns found: {list(df.columns)}")

    out = pd.DataFrame()
    out["species_normalized"] = df[species_col].map(normalize_species_name)
    out["mass_ji"] = pd.to_numeric(df[mass_col], errors="coerce")
    out = out.dropna(subset=["mass_ji"])
    out = out[out["species_normalized"] != "unknown_species"]
    out = out.drop_duplicates(subset=["species_normalized"], keep="first")
    audit = make_audit_rows_from_species_column(df, species_col, path.stem)
    return out[["species_normalized", "mass_ji"]], audit


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


def standardize_bennett_thermal(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
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
    audit_names = []
    for _, row in df.iterrows():
        raw_name = f"{row[genus_col]} {row[species_col]}"
        species_normalized = normalize_species_name(raw_name)
        if species_normalized == "unknown_species":
            continue

        ct_max = pd.NA
        ct_min = pd.NA
        if metric_is_critical_thermal_limit(row[max_metric_col], "ctmax"):
            ct_max = pd.to_numeric(row[tmax_col], errors="coerce")
        if metric_is_critical_thermal_limit(row[min_metric_col], "ctmin"):
            ct_min = pd.to_numeric(row[tmin_col], errors="coerce")

        rows.append({"species_normalized": species_normalized, "ct_max": ct_max, "ct_min": ct_min})
        audit_names.append(raw_name)

    raw = pd.DataFrame(rows)
    if raw.empty:
        thermal = pd.DataFrame(columns=["species_normalized", "ct_min", "ct_max"])
    else:
        collapsed_rows = []
        for species_normalized, group in raw.groupby("species_normalized", sort=True):
            ctmax_series = group["ct_max"].dropna()
            ctmin_series = group["ct_min"].dropna()
            collapsed_rows.append(
                {
                    "species_normalized": species_normalized,
                    "ct_min": ctmin_series.iloc[0] if len(ctmin_series) else pd.NA,
                    "ct_max": ctmax_series.iloc[0] if len(ctmax_series) else pd.NA,
                }
            )
        thermal = pd.DataFrame(collapsed_rows)

    audit = pd.DataFrame(
        {
            "source_dataset": path.stem,
            "source_species_name": audit_names,
            "species_normalized": [normalize_species_name(x) for x in audit_names],
        }
    )
    if not audit.empty:
        audit = audit[audit["species_normalized"] != "unknown_species"]
        audit = audit.drop_duplicates(subset=["source_dataset", "source_species_name", "species_normalized"], keep="first")
    return thermal, audit


def merge_metric(base: pd.DataFrame, metric_df: pd.DataFrame, metric_col: str) -> pd.DataFrame:
    if metric_df.empty:
        base[metric_col] = pd.NA
        return base
    return base.merge(metric_df[["species_normalized", metric_col]], on="species_normalized", how="left")


def classify_match(species_normalized: str, natural_history_species: set[str], natural_history_genera: set[str]) -> str:
    if species_normalized in natural_history_species:
        return "exact_match"
    genus = species_normalized.split("_", 1)[0] if species_normalized else ""
    if genus and genus in natural_history_genera:
        return "genus_match_only"
    return "no_match"


def build_natural_history(
    reference_species: pd.DataFrame,
    genome_size_df: pd.DataFrame,
    meiri_df: pd.DataFrame,
    title_df: pd.DataFrame,
    ji_df: pd.DataFrame,
    thermal_df: pd.DataFrame,
    include_all_source_species: bool,
) -> pd.DataFrame:
    species_frames = [reference_species[["species_normalized"]]]
    if include_all_source_species:
        for df in [genome_size_df, meiri_df, title_df, ji_df, thermal_df]:
            if not df.empty and "species_normalized" in df.columns:
                species_frames.append(df[["species_normalized"]])

    base = pd.concat(species_frames, ignore_index=True)
    base = base.dropna(subset=["species_normalized"])
    base = base[base["species_normalized"] != "unknown_species"]
    base = base.drop_duplicates(subset=["species_normalized"], keep="first")
    base = base.sort_values("species_normalized").reset_index(drop=True)

    base = merge_metric(base, meiri_df, "mass_meiri")
    base = merge_metric(base, title_df, "mass_title")
    base = merge_metric(base, ji_df, "mass_ji")
    base = merge_metric(base, genome_size_df, "genome_size")

    if thermal_df.empty:
        base["ct_min"] = pd.NA
        base["ct_max"] = pd.NA
    else:
        base = base.merge(thermal_df[["species_normalized", "ct_min", "ct_max"]], on="species_normalized", how="left")

    base.insert(0, "species_pk", range(1, len(base) + 1))

    return base[["species_pk", "species_normalized", "mass_meiri", "mass_title", "mass_ji", "genome_size", "ct_min", "ct_max"]]


def build_species_name_audit(audit_frames: list[pd.DataFrame], natural_history_df: pd.DataFrame) -> pd.DataFrame:
    if audit_frames:
        audit = pd.concat([x for x in audit_frames if x is not None and not x.empty], ignore_index=True)
    else:
        audit = pd.DataFrame(columns=["source_dataset", "source_species_name", "species_normalized"])

    if audit.empty:
        return pd.DataFrame(
            columns=["species_name_audit_pk", "species_pk", "source_dataset", "source_species_name", "species_normalized", "match_status"]
        )

    audit = audit.dropna(subset=["species_normalized"])
    audit = audit[audit["species_normalized"] != "unknown_species"]
    audit = audit.drop_duplicates(subset=["source_dataset", "source_species_name", "species_normalized"], keep="first")

    species_to_pk = dict(zip(natural_history_df["species_normalized"], natural_history_df["species_pk"]))
    natural_history_species = set(natural_history_df["species_normalized"])
    natural_history_genera = {x.split("_", 1)[0] for x in natural_history_species if isinstance(x, str) and x}

    audit["species_pk"] = audit["species_normalized"].map(species_to_pk)
    audit["match_status"] = audit["species_normalized"].map(
        lambda x: classify_match(str(x), natural_history_species, natural_history_genera)
    )

    audit = audit.sort_values(["source_dataset", "species_normalized", "source_species_name"]).reset_index(drop=True)
    audit.insert(0, "species_name_audit_pk", range(1, len(audit) + 1))

    return audit[["species_name_audit_pk", "species_pk", "source_dataset", "source_species_name", "species_normalized", "match_status"]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create SQL-ready natural_history.tsv and species_name_audit.tsv files."
    )
    parser.add_argument("genomes_dir", type=Path, help="Path to the top-level genomes directory.")
    parser.add_argument("--genomes-metadata", type=Path, default=None, help="genomes_metadata.csv used as a project species reference.")
    parser.add_argument("--manifest", type=Path, default=None, help="Optional project manifest CSV/TSV used as an additional species reference.")
    parser.add_argument("--genome-size", type=Path, default=None, help="C-value/genome-size summary CSV/TSV.")
    parser.add_argument("--mass-meiri", type=Path, default=None, help="Meiri-style mass CSV/TSV.")
    parser.add_argument("--mass-title", type=Path, default=None, help="Title-style mass CSV/TSV.")
    parser.add_argument("--mass-ji", type=Path, default=None, help="Ji/manual mass TSV/CSV.")
    parser.add_argument("--thermal-bennett", type=Path, default=None, help="Bennett/GlobTherm-style thermal limits CSV/TSV.")
    parser.add_argument("--outdir", type=Path, default=None, help="Output directory. Default: <genomes_dir>/records/sql_tsvs")
    parser.add_argument(
        "--include-all-source-species",
        action="store_true",
        help=(
            "Include species found only in natural-history sources in natural_history.tsv. "
            "Default keeps natural_history.tsv limited to genomes_metadata/manifest species."
        ),
    )
    parser.add_argument(
        "--allow-missing-natural-history",
        action="store_true",
        help="Write outputs with blank metric values instead of failing when natural-history files are missing.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    genomes_dir = args.genomes_dir.expanduser().resolve()
    if not genomes_dir.exists():
        raise FileNotFoundError(f"genomes_dir does not exist: {genomes_dir}")

    nh_dir = natural_history_dir(genomes_dir)
    outdir = args.outdir.expanduser().resolve() if args.outdir else default_sql_tsv_dir(genomes_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    reference_frames: list[pd.DataFrame] = []
    audit_frames: list[pd.DataFrame] = []

    genomes_metadata_path = args.genomes_metadata.expanduser().resolve() if args.genomes_metadata else default_genomes_metadata(genomes_dir)
    if genomes_metadata_path.exists():
        ref, audit = load_reference_species_from_file(genomes_metadata_path, "genomes_metadata")
        reference_frames.append(ref)
        audit_frames.append(audit)
        print(f"[OK] Loaded {len(ref)} reference species from {genomes_metadata_path}")
    elif args.genomes_metadata:
        raise FileNotFoundError(f"genomes_metadata file does not exist: {genomes_metadata_path}")
    else:
        print(f"[WARN] Default genomes_metadata.csv not found: {genomes_metadata_path}")

    if args.manifest:
        manifest_path = args.manifest.expanduser().resolve()
        if not manifest_path.exists():
            raise FileNotFoundError(f"manifest file does not exist: {manifest_path}")
        ref, audit = load_reference_species_from_file(manifest_path, "manifest")
        reference_frames.append(ref)
        audit_frames.append(audit)
        print(f"[OK] Loaded {len(ref)} reference species from {manifest_path}")

    if not reference_frames:
        raise ValueError(
            "No reference species were loaded. Provide --genomes-metadata or --manifest, "
            "or place genomes_metadata.csv at <genomes_dir>/records/genomes_metadata.csv."
        )

    reference_species = pd.concat(reference_frames, ignore_index=True)
    reference_species = reference_species.drop_duplicates(subset=["species_normalized"], keep="first")

    genome_size_path = args.genome_size.expanduser().resolve() if args.genome_size else find_default_genome_size_file(nh_dir)
    mass_meiri_path = args.mass_meiri.expanduser().resolve() if args.mass_meiri else find_default_file(nh_dir, ["meiri"])
    mass_title_path = args.mass_title.expanduser().resolve() if args.mass_title else find_default_file(nh_dir, ["title"])
    mass_ji_path = args.mass_ji.expanduser().resolve() if args.mass_ji else find_default_file(nh_dir, ["ji"])
    thermal_bennett_path = args.thermal_bennett.expanduser().resolve() if args.thermal_bennett else find_default_file(nh_dir, ["bennett"])

    missing_messages = []
    if genome_size_path is None:
        missing_messages.append("No genome-size/C-value summary file was provided or auto-detected.")
    if mass_meiri_path is None and mass_title_path is None and mass_ji_path is None:
        missing_messages.append("No Meiri, Title, or Ji mass file was provided or auto-detected.")
    if thermal_bennett_path is None:
        missing_messages.append("No Bennett/GlobTherm thermal file was provided or auto-detected.")
    if missing_messages and not args.allow_missing_natural_history:
        raise FileNotFoundError("\n".join(missing_messages))
    for message in missing_messages:
        print(f"[WARN] {message}")

    genome_size_df = pd.DataFrame(columns=["species_normalized", "genome_size"])
    meiri_df = pd.DataFrame(columns=["species_normalized", "mass_meiri"])
    title_df = pd.DataFrame(columns=["species_normalized", "mass_title"])
    ji_df = pd.DataFrame(columns=["species_normalized", "mass_ji"])
    thermal_df = pd.DataFrame(columns=["species_normalized", "ct_min", "ct_max"])

    if genome_size_path is not None:
        genome_size_df, audit = standardize_genome_size(genome_size_path)
        audit_frames.append(audit)
        print(f"[OK] Loaded genome_size for {len(genome_size_df)} species from {genome_size_path}")
    if mass_meiri_path is not None:
        meiri_df, audit = standardize_meiri_mass(mass_meiri_path)
        audit_frames.append(audit)
        print(f"[OK] Loaded mass_meiri for {len(meiri_df)} species from {mass_meiri_path}")
    if mass_title_path is not None:
        title_df, audit = standardize_title_mass(mass_title_path)
        audit_frames.append(audit)
        print(f"[OK] Loaded mass_title for {len(title_df)} species from {mass_title_path}")
    if mass_ji_path is not None:
        ji_df, audit = standardize_ji_mass(mass_ji_path)
        audit_frames.append(audit)
        print(f"[OK] Loaded mass_ji for {len(ji_df)} species from {mass_ji_path}")
    if thermal_bennett_path is not None:
        thermal_df, audit = standardize_bennett_thermal(thermal_bennett_path)
        audit_frames.append(audit)
        print(f"[OK] Loaded thermal limits for {len(thermal_df)} species from {thermal_bennett_path}")

    natural_history_df = build_natural_history(
        reference_species=reference_species,
        genome_size_df=genome_size_df,
        meiri_df=meiri_df,
        title_df=title_df,
        ji_df=ji_df,
        thermal_df=thermal_df,
        include_all_source_species=args.include_all_source_species,
    )

    species_name_audit_df = build_species_name_audit(audit_frames, natural_history_df)

    natural_history_out = outdir / "natural_history.tsv"
    species_name_audit_out = outdir / "species_name_audit.tsv"
    write_tsv(natural_history_df, natural_history_out)
    write_tsv(species_name_audit_df, species_name_audit_out)

    print(f"[OK] Wrote {natural_history_out} ({len(natural_history_df)} species)")
    print(f"[OK] Wrote {species_name_audit_out} ({len(species_name_audit_df)} audit rows)")
    print(
        "[INFO] natural_history non-null counts: "
        f"mass_meiri={int(natural_history_df['mass_meiri'].notna().sum())}, "
        f"mass_title={int(natural_history_df['mass_title'].notna().sum())}, "
        f"mass_ji={int(natural_history_df['mass_ji'].notna().sum())}, "
        f"genome_size={int(natural_history_df['genome_size'].notna().sum())}, "
        f"ct_min={int(natural_history_df['ct_min'].notna().sum())}, "
        f"ct_max={int(natural_history_df['ct_max'].notna().sum())}"
    )
    print("[DONE] SQL TSV export complete.")


if __name__ == "__main__":
    main()
