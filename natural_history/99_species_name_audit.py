#!/usr/bin/env python3
"""
99_species_name_audit.py

Create a SQL-ready species-name audit TSV from one or more CSV/TSV files.

Output schema:
  species_name_audit_pk
  species_pk
  source_dataset
  source_species_name
  species_normalized
  match_status

Intended use:
  Keep this script in genomes_database/natural_history and save the output TSV
  beside the natural-history staging files so it can be uploaded to the final SQL
  species_name_audit table.

Examples:
  # Audit several known natural-history files, without species_pk values yet
  python 99_species_name_audit.py \
    --input meiri.csv --input title.csv --input ji.tsv --input bennett.csv \
    --out species_name_audit.tsv

  # Audit files and fill species_pk using an exported natural_history table
  python 99_species_name_audit.py \
    --input meiri.csv --input title.csv --input ji.tsv --input bennett.csv \
    --species-lookup natural_history.tsv \
    --out species_name_audit.tsv

  # Provide explicit source labels
  python 99_species_name_audit.py \
    --input meiri.csv --label meiri \
    --input title.csv --label title \
    --input ji.tsv --label ji \
    --out species_name_audit.tsv

Notes:
  - If --species-lookup is not supplied, species_pk is left blank.
  - If supplied, --species-lookup must contain species_normalized and species_pk.
  - The script auto-detects common species-name columns.
  - Bennett/GlobTherm-style files with separate Genus and Species columns are handled.
"""

import argparse
import csv
import re
from pathlib import Path
from typing import Optional

import pandas as pd


BAD_SPECIES_TOKENS = {"sp", "sp.", "spp", "spp.", "cf", "cf.", "aff", "aff.", "nr", "nr."}
HYBRID_MARKERS = {"x", "×"}

COMMON_SPECIES_COLUMNS = [
    "source_species_name",
    "species_normalized",
    "species_key",
    "genus_species",
    "genus_species_key",
    "normalized_species",
    "normalized_name",
    "organism_name",
    "organismName",
    "Organism Name",
    "species_name",
    "species",
    "binomial",
    "binomial_2020",
    "binomial_(original files)",
    "treename",
]


def normalize_species_name(value: object) -> str:
    """Normalize a source species name to lowercase genus_species format."""
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


def choose_first_column(df: pd.DataFrame, candidates: list[str]) -> Optional[str]:
    cleaned_to_real = {str(c).strip().lower(): c for c in df.columns}
    for candidate in candidates:
        key = candidate.strip().lower()
        if key in cleaned_to_real:
            return cleaned_to_real[key]
    return None


def infer_source_dataset(path: Path) -> str:
    stem = path.stem.lower()
    for key in ["meiri", "title", "ji", "bennett", "globtherm", "genome_size", "c_value", "manifest", "metadata"]:
        if key in stem:
            return key
    return path.stem


def extract_source_species_names(path: Path, species_column: Optional[str] = None) -> pd.Series:
    df = read_table_auto(path)

    if species_column:
        if species_column not in df.columns:
            raise ValueError(f"Requested --species-column '{species_column}' not found in {path}. Columns: {list(df.columns)}")
        return df[species_column].astype(str).str.strip()

    genus_col = choose_first_column(df, ["Genus", "genus"])
    species_col = choose_first_column(df, ["Species", "species"])
    if genus_col is not None and species_col is not None:
        return (df[genus_col].astype(str).str.strip() + " " + df[species_col].astype(str).str.strip()).str.strip()

    detected = choose_first_column(df, COMMON_SPECIES_COLUMNS)
    if detected is None:
        raise ValueError(
            f"Could not identify a species-name column in {path}. "
            f"Use --species-column if needed. Columns found: {list(df.columns)}"
        )
    return df[detected].astype(str).str.strip()


def load_species_lookup(path: Optional[Path]) -> dict[str, object]:
    if path is None:
        return {}

    df = read_table_auto(path)
    species_col = choose_first_column(df, ["species_normalized", "genus_species", "species_key"])
    pk_col = choose_first_column(df, ["species_pk"])
    if species_col is None or pk_col is None:
        raise ValueError(
            f"species lookup must contain species_pk and species_normalized/genus_species. "
            f"Columns found in {path}: {list(df.columns)}"
        )

    lookup = {}
    for _, row in df.iterrows():
        key = normalize_species_name(row[species_col])
        if key != "unknown_species" and pd.notna(row[pk_col]):
            lookup[key] = row[pk_col]
    return lookup


def build_audit_rows(
    input_paths: list[Path],
    labels: list[str],
    species_lookup: dict[str, object],
    species_column: Optional[str],
    start_pk: int,
    keep_duplicate_source_names: bool,
) -> pd.DataFrame:
    rows = []
    next_pk = start_pk

    for index, path in enumerate(input_paths):
        source_dataset = labels[index] if index < len(labels) and labels[index] else infer_source_dataset(path)
        source_names = extract_source_species_names(path, species_column=species_column)

        source_df = pd.DataFrame({"source_species_name": source_names})
        source_df["source_species_name"] = source_df["source_species_name"].astype(str).str.strip()
        source_df = source_df[source_df["source_species_name"].str.len() > 0]
        source_df = source_df[~source_df["source_species_name"].str.lower().isin({"nan", "none", "na", "n/a"})]
        source_df["species_normalized"] = source_df["source_species_name"].map(normalize_species_name)

        if not keep_duplicate_source_names:
            source_df = source_df.drop_duplicates(subset=["source_species_name", "species_normalized"], keep="first")

        duplicated_normalized = source_df["species_normalized"].duplicated(keep=False)

        for row_i, row in source_df.reset_index(drop=True).iterrows():
            species_normalized = row["species_normalized"]
            species_pk = species_lookup.get(species_normalized, pd.NA)

            if species_normalized == "unknown_species":
                match_status = "unknown_species"
            elif species_lookup and pd.notna(species_pk):
                match_status = "matched_species_pk"
            elif species_lookup and pd.isna(species_pk):
                match_status = "normalized_not_in_species_lookup"
            elif duplicated_normalized.iloc[row_i]:
                match_status = "normalized_duplicate_within_source"
            else:
                match_status = "normalized_no_species_lookup_provided"

            rows.append(
                {
                    "species_name_audit_pk": next_pk,
                    "species_pk": species_pk,
                    "source_dataset": source_dataset,
                    "source_species_name": row["source_species_name"],
                    "species_normalized": species_normalized,
                    "match_status": match_status,
                }
            )
            next_pk += 1

    return pd.DataFrame(rows, columns=[
        "species_name_audit_pk",
        "species_pk",
        "source_dataset",
        "source_species_name",
        "species_normalized",
        "match_status",
    ])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create a SQL-ready species_name_audit TSV from natural-history CSV/TSV files.")
    parser.add_argument("--input", dest="inputs", type=Path, action="append", required=True,
                        help="Input CSV/TSV/TXT file. Repeat this argument for multiple files.")
    parser.add_argument("--label", dest="labels", action="append", default=[],
                        help="Optional source_dataset label corresponding to each --input. Repeat as needed.")
    parser.add_argument("--species-column", default=None,
                        help="Optional explicit species-name column to use for all input files.")
    parser.add_argument("--species-lookup", type=Path, default=None,
                        help="Optional natural_history TSV/CSV containing species_pk and species_normalized.")
    parser.add_argument("--out", type=Path, default=Path("species_name_audit.tsv"),
                        help="Output TSV path. Default: species_name_audit.tsv")
    parser.add_argument("--start-pk", type=int, default=1,
                        help="Starting value for species_name_audit_pk. Default: 1")
    parser.add_argument("--keep-duplicate-source-names", action="store_true",
                        help="Keep repeated source_species_name rows instead of collapsing exact duplicates within each source.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    input_paths = [p.expanduser().resolve() for p in args.inputs]
    for path in input_paths:
        if not path.exists():
            raise FileNotFoundError(f"Input file does not exist: {path}")

    if args.labels and len(args.labels) not in {0, len(input_paths)}:
        raise ValueError("If --label is used, provide either no labels or exactly one --label per --input.")

    species_lookup_path = args.species_lookup.expanduser().resolve() if args.species_lookup else None
    species_lookup = load_species_lookup(species_lookup_path)

    audit_df = build_audit_rows(
        input_paths=input_paths,
        labels=args.labels,
        species_lookup=species_lookup,
        species_column=args.species_column,
        start_pk=args.start_pk,
        keep_duplicate_source_names=args.keep_duplicate_source_names,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    audit_df.to_csv(args.out, sep="\t", index=False)

    print(f"[OK] Wrote {args.out} ({len(audit_df)} rows)")
    print(f"[INFO] source datasets: {', '.join(sorted(audit_df['source_dataset'].unique()))}")
    print("[INFO] match_status counts:")
    for status, count in audit_df["match_status"].value_counts().sort_index().items():
        print(f"  {status}: {count}")


if __name__ == "__main__":
    main()
