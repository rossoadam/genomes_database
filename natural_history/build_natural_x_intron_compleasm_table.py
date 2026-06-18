#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def read_tsv(path):
    return pd.read_csv(path, sep="\t")


def pick_min_median_max_species(ortholog_summary):
    df = ortholog_summary[["species_pk", "sd_gc3"]].copy()
    df = df.dropna(subset=["species_pk", "sd_gc3"])
    df = df.sort_values("sd_gc3").reset_index(drop=True)

    if df.empty:
        raise ValueError("ortholog_summary has no usable sd_gc3 values.")

    min_row = df.iloc[0]
    max_row = df.iloc[-1]
    median_row = df.iloc[len(df) // 2]

    selected = pd.DataFrame(
        [
            {
                "species_pk": min_row["species_pk"],
                "sd_gc3_category": "minimum",
                "sd_gc3": min_row["sd_gc3"],
            },
            {
                "species_pk": median_row["species_pk"],
                "sd_gc3_category": "median",
                "sd_gc3": median_row["sd_gc3"],
            },
            {
                "species_pk": max_row["species_pk"],
                "sd_gc3_category": "maximum",
                "sd_gc3": max_row["sd_gc3"],
            },
        ]
    )

    return selected.drop_duplicates(subset=["species_pk", "sd_gc3_category"])


def main():
    parser = argparse.ArgumentParser(
        description="Create intron GC vs intron length table for min, median, and max sd_gc3 species."
    )

    parser.add_argument("--intron-compleasm", required=True)
    parser.add_argument("--orthologs", required=True)
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--genomes", required=True)
    parser.add_argument("--natural-history", required=True)
    parser.add_argument("--ortholog-summary", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    introns = read_tsv(args.intron_compleasm)
    orthologs = read_tsv(args.orthologs)
    sequences = read_tsv(args.sequences)
    genomes = read_tsv(args.genomes)
    natural_history = read_tsv(args.natural_history)
    ortholog_summary = read_tsv(args.ortholog_summary)

    selected_species = pick_min_median_max_species(ortholog_summary)

    merged = (
        introns[["ortholog_pk", "intron_id", "length", "gc"]]
        .merge(
            orthologs[["ortholog_pk", "sequence_pk"]],
            on="ortholog_pk",
            how="left",
            validate="many_to_one",
        )
        .merge(
            sequences[["sequence_pk", "genome_pk"]],
            on="sequence_pk",
            how="left",
            validate="many_to_one",
        )
        .merge(
            genomes[["genome_pk", "species_pk"]],
            on="genome_pk",
            how="left",
            validate="many_to_one",
        )
        .merge(
            selected_species,
            on="species_pk",
            how="inner",
            validate="many_to_one",
        )
        .merge(
            natural_history[["species_pk", "species_normalized"]],
            on="species_pk",
            how="left",
            validate="many_to_one",
        )
    )

    output = merged.rename(columns={"length": "intron_length"})[
        [
            "species_pk",
            "species_normalized",
            "sd_gc3_category",
            "sd_gc3",
            "intron_id",
            "intron_length",
            "gc",
        ]
    ]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, sep="\t", index=False)

    print(f"Wrote: {output_path}")
    print(f"Rows: {len(output)}")
    print("Selected species:")
    print(
        output[["species_pk", "species_normalized", "sd_gc3_category", "sd_gc3"]]
        .drop_duplicates()
        .sort_values("sd_gc3")
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()