#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def first_non_missing(row, columns):
    for col in columns:
        if col in row.index and pd.notna(row[col]):
            return row[col]
    return pd.NA


def main():
    parser = argparse.ArgumentParser(
        description="Join natural_history.tsv and orthologs_summary.tsv for R analysis."
    )
    parser.add_argument(
        "--natural-history",
        required=True,
        help="Path to natural_history.tsv",
    )
    parser.add_argument(
        "--ortholog-summary",
        required=True,
        help="Path to orthologs_summary.tsv",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output path for combined TSV",
    )

    args = parser.parse_args()

    natural_history_path = Path(args.natural_history)
    ortholog_summary_path = Path(args.ortholog_summary)
    output_path = Path(args.output)

    natural_history = pd.read_csv(natural_history_path, sep="\t")
    orthologs_summary = pd.read_csv(ortholog_summary_path, sep="\t")

    required_natural = [
        "species_pk",
        "species_normalized",
        "ct_max",
        "genome_size",
        "mass_meiri",
        "mass_title",
        "mass_ji",
    ]

    required_orthologs = [
        "species_pk",
        "n_orthologs",
        "mean_gc3",
        "sd_gc3",
        "mean_gc4",
        "sd_gc4",
    ]

    missing_natural = [c for c in required_natural if c not in natural_history.columns]
    missing_orthologs = [c for c in required_orthologs if c not in orthologs_summary.columns]

    if missing_natural:
        raise ValueError(f"Missing columns in natural_history.tsv: {missing_natural}")

    if missing_orthologs:
        raise ValueError(f"Missing columns in orthologs_summary.tsv: {missing_orthologs}")

    natural_history["mass_preferred"] = natural_history.apply(
        first_non_missing,
        axis=1,
        columns=["mass_meiri", "mass_title", "mass_ji"],
    )

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

    combined = orthologs_summary[ortholog_keep].merge(
        natural_history[natural_keep],
        on="species_pk",
        how="left",
        validate="many_to_one",
    )

    final_cols = [
        "species_pk",
        "species_normalized",
        "n_orthologs",
        "mean_gc3",
        "sd_gc3",
        "mean_gc4",
        "sd_gc4",
        "ct_max",
        "genome_size",
        "mass_preferred",
    ]

    combined = combined[final_cols]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_path, sep="\t", index=False)

    print(f"Wrote: {output_path}")
    print(f"Rows: {len(combined)}")
    print(f"Species without natural history match: {combined['species_normalized'].isna().sum()}")


if __name__ == "__main__":
    main()
