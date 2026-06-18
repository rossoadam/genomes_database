#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def read_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def require_columns(df: pd.DataFrame, required: list[str], table_name: str) -> None:
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {table_name}: {missing}")


def add_mass_preferred(natural_history: pd.DataFrame) -> pd.DataFrame:
    natural_history = natural_history.copy()

    natural_history["mass_preferred"] = (
        natural_history["mass_meiri"]
        .combine_first(natural_history["mass_title"])
        .combine_first(natural_history["mass_ji"])
    )

    return natural_history


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Join natural_history.tsv, ortholog_summary.tsv, and "
            "intron_compleasm_summary.tsv by species_pk."
        )
    )

    parser.add_argument("--natural-history", required=True)
    parser.add_argument("--ortholog-summary", required=True)
    parser.add_argument("--intron-compleasm-summary", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    natural_history = read_tsv(args.natural_history)
    ortholog_summary = read_tsv(args.ortholog_summary)
    intron_summary = read_tsv(args.intron_compleasm_summary)

    require_columns(
        natural_history,
        [
            "species_pk",
            "species_normalized",
            "ct_max",
            "genome_size",
            "mass_meiri",
            "mass_title",
            "mass_ji",
        ],
        "natural_history",
    )

    require_columns(
        ortholog_summary,
        [
            "species_pk",
            "n_orthologs",
            "mean_gc3",
            "sd_gc3",
            "mean_gc4",
            "sd_gc4",
        ],
        "ortholog_summary",
    )

    require_columns(
        intron_summary,
        [
            "species_pk",
            "mean_gc",
        ],
        "intron_compleasm_summary",
    )

    natural_history = add_mass_preferred(natural_history)

    natural_keep = [
        "species_pk",
        "species_normalized",
        "ct_max",
        "genome_size",
        "mass_preferred",
    ]

    ortholog_keep = [
        "species_pk",
        "n_orthologs",
        "mean_gc3",
        "sd_gc3",
        "mean_gc4",
        "sd_gc4",
    ]

    intron_keep = [
        "species_pk",
        "mean_gc",
    ]

    combined = (
        ortholog_summary[ortholog_keep]
        .merge(
            intron_summary[intron_keep],
            on="species_pk",
            how="left",
            validate="one_to_one",
        )
        .merge(
            natural_history[natural_keep],
            on="species_pk",
            how="left",
            validate="many_to_one",
        )
    )

    combined = combined.rename(
        columns={
            "mean_gc": "mean_intron_gc",
        }
    )

    final_cols = [
        "species_pk",
        "species_normalized",
        "n_orthologs",
        "mean_gc3",
        "sd_gc3",
        "mean_gc4",
        "sd_gc4",
        "mean_intron_gc",
        "ct_max",
        "genome_size",
        "mass_preferred",
    ]

    combined = combined[final_cols]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    combined.to_csv(output_path, sep="\t", index=False)

    print(f"Wrote: {output_path}")
    print(f"Rows: {len(combined)}")
    print(f"Species without intron summary: {combined['mean_intron_gc'].isna().sum()}")
    print(f"Species without natural history match: {combined['species_normalized'].isna().sum()}")


if __name__ == "__main__":
    main()

