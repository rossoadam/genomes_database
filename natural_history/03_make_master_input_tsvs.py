#!/usr/bin/env python3
"""
01_make_master_input_tsvs.py

Create standardized TSV input files needed by the downstream GC metrics/master
export workflow.

Design goals:
  - The first positional argument is always the top-level genomes directory.
  - Data files can live under:
        <genomes_dir>/records/natural_history/
    while this script can live in the GitHub repo, for example:
        genomes_database/natural_history/
  - The project manifest, if supplied with --manifest, defines the focal species
    set and supersedes Compleasm metadata.
  - If --manifest is not supplied, the script uses:
        <genomes_dir>/records/compleasm/records/metadata.csv
    to define the focal species set.
  - No metadata.tsv is written. The manifest can be used downstream for
    accession/species metadata.
  - species_mass.tsv keeps both Title and Meiri mass columns.
  - genome_sizes.tsv uses C-value summary data, especially mean_c_value_pg,
    rather than calculating genome size from FASTA files.

Main outputs:
  <outdir>/species_mass.tsv
      genus_species
      mass_meiri
      mass_title
      mass_preferred
      mass_source_preferred

  <outdir>/genome_sizes.tsv
      genus_species
      species_name
      mean_c_value_pg
      mean_genome_mb
      mean_genome_gb
      genome_size
      genome_size_source

Default output directory:
  <genomes_dir>/records/gc_metrics_inputs/

Typical usage:
  python 01_make_master_input_tsvs.py /path/to/genomes \
    --manifest /path/to/project_manifest.csv \
    --genome-size /path/to/reptile_c_value_summary.csv \
    --mass-meiri /path/to/meiri.csv \
    --mass-title /path/to/title.csv

If your natural history files are in <genomes_dir>/records/natural_history/,
you can omit --mass-meiri and/or --mass-title if their filenames contain
"meiri" or "title" respectively. You can also omit --genome-size if a summary
file with columns such as mean_c_value_pg is present in that same directory.
"""

import argparse
import csv
from pathlib import Path
from typing import Optional

import pandas as pd


def normalize_species_name(organism_name: str) -> str:
    """
    Normalize names using the same simple Genus_species rule planned for SQL.

    Examples:
      "Sphaerodactylus ariasae" -> "Sphaerodactylus_ariasae"
      "Sphaerodactylus_ariasae" -> "Sphaerodactylus_ariasae"
      "Gekko_kuhli_v2" -> "Gekko_kuhli"
    """
    organism_name = str(organism_name).strip()
    if "_" in organism_name:
        parts = [p for p in organism_name.split("_") if p]
    else:
        parts = organism_name.split()

    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    elif len(parts) == 1:
        return parts[0]
    return "unknown_species"


def read_table_auto(path: Path) -> pd.DataFrame:
    """Read CSV or TSV. Handles UTF-8 BOM in files such as some Meiri CSVs."""
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
    """Choose a column robustly, ignoring case and accidental whitespace."""
    cleaned_to_real = {str(c).strip().lower(): c for c in df.columns}
    for candidate in candidates:
        key = candidate.strip().lower()
        if key in cleaned_to_real:
            return cleaned_to_real[key]
    return None


def default_compleasm_metadata(genomes_dir: Path) -> Path:
    return genomes_dir / "records" / "compleasm" / "records" / "metadata.csv"


def natural_history_dir(genomes_dir: Path) -> Path:
    return genomes_dir / "records" / "natural_history"


def find_default_file(
    directory: Path,
    include_terms: list[str],
    required_columns: Optional[list[str]] = None,
) -> Optional[Path]:
    """
    Find a likely default input file by filename terms and, optionally, columns.

    This lets files live under genomes/records/natural_history/ without forcing
    hard-coded filenames.
    """
    if not directory.exists():
        return None

    candidates = []
    for path in directory.iterdir():
        if not path.is_file():
            continue
        lower_name = path.name.lower()
        if not lower_name.endswith((".csv", ".tsv", ".txt")):
            continue
        if all(term.lower() in lower_name for term in include_terms):
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
    """
    Find a C-value summary file under records/natural_history.

    The key required column is mean_c_value_pg.
    """
    if not directory.exists():
        return None

    preferred_terms = [
        ["genome", "size"],
        ["c", "value"],
        ["summary"],
    ]

    for terms in preferred_terms:
        candidate = find_default_file(directory, terms, required_columns=["mean_c_value_pg"])
        if candidate is not None:
            return candidate

    # Last fallback: inspect all CSV/TSV/TXT files for mean_c_value_pg.
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
    """
    Read the manifest or Compleasm metadata and return the focal species set.

    No metadata.tsv is written; this is only used internally to filter/shape
    the natural-history outputs.
    """
    df = read_table_auto(source_path)

    accession_col = choose_first_column(
        df,
        [
            "accession",
            "assembly_accession",
            "assemblyAccession",
            "Assembly Accession",
            "genome_accession",
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
            "genus_species",
            "binomial",
            "binomial_2020",
        ],
    )

    if organism_col is None:
        raise ValueError(
            f"Could not find an organism/species column in {source_path}. "
            f"Columns found: {list(df.columns)}"
        )

    out = pd.DataFrame()
    out["organism_name"] = df[organism_col].astype(str).str.strip()
    out["genus_species"] = out["organism_name"].map(normalize_species_name)

    if accession_col is not None:
        out["accession"] = df[accession_col].astype(str).str.strip()
    else:
        out["accession"] = pd.NA

    out["project_species_source"] = source_label
    out = out[out["genus_species"].astype(str).str.len() > 0]
    out = out[out["genus_species"] != "unknown_species"]
    out = out.drop_duplicates(subset=["genus_species"], keep="first")
    return out


def standardize_genome_size(path: Path, project_species: pd.DataFrame) -> pd.DataFrame:
    """
    Build genome_sizes.tsv from a C-value summary table.

    The downstream genome_size column is populated from mean_c_value_pg.
    """
    df = read_table_auto(path)

    cvalue_col = choose_first_column(df, ["mean_c_value_pg"])
    species_key_col = choose_first_column(df, ["species_key", "genus_species"])
    species_name_col = choose_first_column(df, ["species_name", "organism_name", "organismName"])
    mean_mb_col = choose_first_column(df, ["mean_genome_mb"])
    mean_gb_col = choose_first_column(df, ["mean_genome_gb"])

    if cvalue_col is None:
        raise ValueError(
            f"Could not find mean_c_value_pg in {path}. "
            f"Columns found: {list(df.columns)}"
        )

    if species_key_col is None and species_name_col is None:
        raise ValueError(
            f"Could not find species_key or species_name in {path}. "
            f"Columns found: {list(df.columns)}"
        )

    rows = []
    for _, row in df.iterrows():
        if species_key_col is not None and pd.notna(row[species_key_col]):
            genus_species = normalize_species_name(row[species_key_col])
        else:
            genus_species = normalize_species_name(row[species_name_col])

        species_name = row[species_name_col] if species_name_col is not None else str(genus_species).replace("_", " ")

        rows.append(
            {
                "genus_species": genus_species,
                "species_name": species_name,
                "mean_c_value_pg": pd.to_numeric(row[cvalue_col], errors="coerce"),
                "mean_genome_mb": pd.to_numeric(row[mean_mb_col], errors="coerce") if mean_mb_col else pd.NA,
                "mean_genome_gb": pd.to_numeric(row[mean_gb_col], errors="coerce") if mean_gb_col else pd.NA,
                "genome_size_source": str(path),
            }
        )

    size_df = pd.DataFrame(rows)
    size_df = size_df.dropna(subset=["genus_species", "mean_c_value_pg"])
    size_df = size_df[size_df["genus_species"] != "unknown_species"]
    size_df = size_df.drop_duplicates(subset=["genus_species"], keep="first")

    out = project_species[["genus_species"]].copy()
    out = out.merge(size_df, on="genus_species", how="left")

    # Generic downstream column requested by prior workflow.
    out["genome_size"] = out["mean_c_value_pg"]

    out = out[
        [
            "genus_species",
            "species_name",
            "mean_c_value_pg",
            "mean_genome_mb",
            "mean_genome_gb",
            "genome_size",
            "genome_size_source",
        ]
    ]
    return out


def standardize_meiri_mass(path: Path) -> pd.DataFrame:
    """
    Parse Meiri-style mass data.

    Example columns:
      binomial_2020
      body mass (g)
      adult_body_mass (g)
    """
    df = read_table_auto(path)

    species_col = choose_first_column(
        df,
        ["binomial_2020", "binomial_(original files)", "genus_species", "species", "species_name"],
    )

    # Prefer the explicit body mass column shown by the user.
    mass_col = choose_first_column(
        df,
        ["body mass (g)", "adult_body_mass (g)", "adult_body_mass_g", "mass"],
    )

    if species_col is None or mass_col is None:
        raise ValueError(
            f"Could not identify Meiri species/mass columns in {path}. "
            f"Columns found: {list(df.columns)}"
        )

    out = pd.DataFrame()
    out["genus_species"] = df[species_col].map(normalize_species_name)
    out["mass_meiri"] = pd.to_numeric(df[mass_col], errors="coerce")
    out = out.dropna(subset=["genus_species", "mass_meiri"])
    out = out[out["genus_species"].astype(str).str.len() > 0]
    out = out[out["genus_species"] != "unknown_species"]
    out = out.drop_duplicates(subset=["genus_species"], keep="first")
    return out


def standardize_title_mass(path: Path) -> pd.DataFrame:
    """
    Parse Title-style mass data.

    Example columns:
      treename
      mass
    """
    df = read_table_auto(path)

    species_col = choose_first_column(df, ["treename", "genus_species", "species", "species_name", "binomial"])
    mass_col = choose_first_column(df, ["mass"])

    if species_col is None or mass_col is None:
        raise ValueError(
            f"Could not identify Title species/mass columns in {path}. "
            f"Columns found: {list(df.columns)}"
        )

    out = pd.DataFrame()
    out["genus_species"] = df[species_col].map(normalize_species_name)
    out["mass_title"] = pd.to_numeric(df[mass_col], errors="coerce")
    out = out.dropna(subset=["genus_species", "mass_title"])
    out = out[out["genus_species"].astype(str).str.len() > 0]
    out = out[out["genus_species"] != "unknown_species"]
    out = out.drop_duplicates(subset=["genus_species"], keep="first")
    return out


def build_species_mass_tsv(
    project_species: pd.DataFrame,
    mass_meiri_path: Optional[Path],
    mass_title_path: Optional[Path],
) -> pd.DataFrame:
    """
    Build species_mass.tsv with both Meiri and Title columns.

    Priority is only used for mass_preferred:
      Meiri > Title
    """
    out = project_species[["genus_species"]].copy()

    if mass_meiri_path is not None:
        meiri = standardize_meiri_mass(mass_meiri_path)
        out = out.merge(meiri, on="genus_species", how="left")
    else:
        out["mass_meiri"] = pd.NA

    if mass_title_path is not None:
        title = standardize_title_mass(mass_title_path)
        out = out.merge(title, on="genus_species", how="left")
    else:
        out["mass_title"] = pd.NA

    out["mass_preferred"] = out["mass_meiri"].combine_first(out["mass_title"])

    def source(row: pd.Series) -> str:
        if pd.notna(row["mass_meiri"]):
            return "meiri"
        if pd.notna(row["mass_title"]):
            return "title"
        return ""

    out["mass_source_preferred"] = out.apply(source, axis=1)

    out = out[
        [
            "genus_species",
            "mass_meiri",
            "mass_title",
            "mass_preferred",
            "mass_source_preferred",
        ]
    ]
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create standardized natural-history TSV inputs for GC metric master export."
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
        help=(
            "Optional project manifest CSV/TSV. If provided, this supersedes "
            "<genomes_dir>/records/compleasm/records/metadata.csv for defining species."
        ),
    )
    parser.add_argument(
        "--genome-size",
        type=Path,
        default=None,
        help=(
            "C-value summary CSV/TSV with mean_c_value_pg. "
            "Default: auto-detect under <genomes_dir>/records/natural_history/."
        ),
    )
    parser.add_argument(
        "--mass-meiri",
        type=Path,
        default=None,
        help=(
            "Meiri-style CSV/TSV. Default: auto-detect a file containing 'meiri' "
            "under <genomes_dir>/records/natural_history/."
        ),
    )
    parser.add_argument(
        "--mass-title",
        type=Path,
        default=None,
        help=(
            "Title-style CSV/TSV. Default: auto-detect a file containing 'title' "
            "under <genomes_dir>/records/natural_history/."
        ),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory. Default: <genomes_dir>/records/gc_metrics_inputs",
    )
    parser.add_argument(
        "--allow-missing-natural-history",
        action="store_true",
        help=(
            "Do not fail if genome-size or mass files are missing. "
            "The script will write rows with blank values when possible."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    genomes_dir = args.genomes_dir.expanduser().resolve()
    if not genomes_dir.exists():
        raise FileNotFoundError(f"genomes_dir does not exist: {genomes_dir}")

    outdir = args.outdir.expanduser().resolve() if args.outdir else genomes_dir / "records" / "gc_metrics_inputs"
    outdir.mkdir(parents=True, exist_ok=True)

    nh_dir = natural_history_dir(genomes_dir)

    if args.manifest is not None:
        project_source_path = args.manifest.expanduser().resolve()
        project_source_label = "manifest"
    else:
        project_source_path = default_compleasm_metadata(genomes_dir)
        project_source_label = "compleasm_metadata"

    if not project_source_path.exists():
        raise FileNotFoundError(
            f"Project species source does not exist: {project_source_path}\n"
            "Provide --manifest, or make sure Compleasm metadata exists at "
            "<genomes_dir>/records/compleasm/records/metadata.csv"
        )

    project_species = load_project_species(project_source_path, project_source_label)
    print(
        f"[OK] Loaded {len(project_species)} focal species from "
        f"{project_source_path} ({project_source_label})"
    )

    genome_size_path = args.genome_size.expanduser().resolve() if args.genome_size else find_default_genome_size_file(nh_dir)
    mass_meiri_path = args.mass_meiri.expanduser().resolve() if args.mass_meiri else find_default_file(nh_dir, ["meiri"])
    mass_title_path = args.mass_title.expanduser().resolve() if args.mass_title else find_default_file(nh_dir, ["title"])

    if genome_size_path is None and not args.allow_missing_natural_history:
        raise FileNotFoundError(
            "No genome-size/C-value summary file was provided or auto-detected.\n"
            "Use --genome-size /path/to/summary.csv, or place a summary with "
            "mean_c_value_pg under <genomes_dir>/records/natural_history/."
        )

    if mass_meiri_path is None and mass_title_path is None and not args.allow_missing_natural_history:
        raise FileNotFoundError(
            "No Meiri or Title mass file was provided or auto-detected.\n"
            "Use --mass-meiri and/or --mass-title, or place files with 'meiri' "
            "and/or 'title' in the filename under <genomes_dir>/records/natural_history/."
        )

    if genome_size_path is not None:
        genome_size_df = standardize_genome_size(genome_size_path, project_species)
    else:
        genome_size_df = project_species[["genus_species"]].copy()
        genome_size_df["species_name"] = pd.NA
        genome_size_df["mean_c_value_pg"] = pd.NA
        genome_size_df["mean_genome_mb"] = pd.NA
        genome_size_df["mean_genome_gb"] = pd.NA
        genome_size_df["genome_size"] = pd.NA
        genome_size_df["genome_size_source"] = ""

    genome_size_out = outdir / "genome_sizes.tsv"
    write_tsv(genome_size_df, genome_size_out)
    n_size = int(genome_size_df["genome_size"].notna().sum())
    print(f"[OK] Wrote {genome_size_out} ({n_size}/{len(genome_size_df)} species matched)")

    mass_df = build_species_mass_tsv(project_species, mass_meiri_path, mass_title_path)
    mass_out = outdir / "species_mass.tsv"
    write_tsv(mass_df, mass_out)

    n_meiri = int(mass_df["mass_meiri"].notna().sum())
    n_title = int(mass_df["mass_title"].notna().sum())
    n_preferred = int(mass_df["mass_preferred"].notna().sum())
    print(
        f"[OK] Wrote {mass_out} "
        f"(meiri={n_meiri}, title={n_title}, preferred={n_preferred}/{len(mass_df)}; "
        "preferred priority=meiri>title)"
    )

    print("[DONE] Natural-history TSV input staging complete.")
    print(f"[INFO] Data input directory expected by default: {nh_dir}")
    print(f"[INFO] Output directory: {outdir}")


if __name__ == "__main__":
    main()
