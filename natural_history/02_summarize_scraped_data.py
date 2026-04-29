#!/usr/bin/env python3

"""
summarize_reptile_c_values.py

Summarize Animal Genome Size Database reptile C-values from:
reptile_c_values_detected_summary.csv

Input columns expected:
species_id,species_url,species_name_from_link,species_name_inferred,
species_key,c_values_pg_detected,n_c_values_detected

Example c_values_pg_detected:
2.49;1.25
2.66;3.50;2.60
"""

"""
USAGE:
    python3 02_summarize_scraped_data.py -i ./reptiles/reptile_c_values_detected_summary.csv -o ./analysis_results

DESCRIPTION:
    Processes raw output from '01_scrape_genome_size.py' to generate descriptive 
    statistics for reptile genome sizes. 

    The script 'explodes' multiple C-value measurements per species into a long-format 
    DataFrame, converts picograms (pg) to Megabases (Mb) and Gigabases (Gb) using 
    the standard 1pg = 978Mb conversion factor, and exports species-level and 
    population-level summaries.

OUTPUTS:
    - expanded_c_values.csv: Every detected measurement as an individual record.
    - species_summary.csv: Aggregated stats (mean, median, std, etc.) per species.
    - overall_summary.csv: Global dataset metrics for the study.
"""

import argparse
from pathlib import Path
import pandas as pd
import numpy as np


def parse_cvals(val):
    """Convert semicolon-delimited C-values into float list."""
    if pd.isna(val):
        return []

    vals = []
    for x in str(val).split(";"):
        x = x.strip()
        if not x:
            continue
        try:
            vals.append(float(x))
        except ValueError:
            pass
    return vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input CSV file",
    )
    ap.add_argument(
        "-o",
        "--outdir",
        default="reptile_summary",
        help="Output directory",
    )
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input)

    # explode all c-values into long format
    rows = []
    for _, r in df.iterrows():
        vals = parse_cvals(r["c_values_pg_detected"])
        for v in vals:
            rows.append(
                {
                    "species_id": r["species_id"],
                    "species_key": r["species_key"],
                    "species_name": r["species_name_inferred"],
                    "c_value_pg": v,
                    "genome_size_mb": v * 978,
                    "genome_size_gb": (v * 978) / 1000.0,
                }
            )

    long_df = pd.DataFrame(rows)

    # save raw expanded rows
    long_df.to_csv(outdir / "expanded_c_values.csv", index=False)

    # summarize by species
    summary = (
        long_df.groupby(["species_key", "species_name"])
        .agg(
            n_measurements=("c_value_pg", "count"),
            min_c_value_pg=("c_value_pg", "min"),
            max_c_value_pg=("c_value_pg", "max"),
            mean_c_value_pg=("c_value_pg", "mean"),
            median_c_value_pg=("c_value_pg", "median"),
            std_c_value_pg=("c_value_pg", "std"),
            mean_genome_mb=("genome_size_mb", "mean"),
            mean_genome_gb=("genome_size_gb", "mean"),
        )
        .reset_index()
    )

    summary["std_c_value_pg"] = summary["std_c_value_pg"].fillna(0)

    summary.to_csv(outdir / "species_summary.csv", index=False)

    # overall stats
    overall = {
        "n_species": summary.shape[0],
        "n_measurements": long_df.shape[0],
        "overall_min_pg": long_df["c_value_pg"].min(),
        "overall_max_pg": long_df["c_value_pg"].max(),
        "overall_mean_pg": long_df["c_value_pg"].mean(),
        "overall_median_pg": long_df["c_value_pg"].median(),
    }

    pd.DataFrame([overall]).to_csv(outdir / "overall_summary.csv", index=False)

    # print report
    print("\n=== SUMMARY ===")
    print(f"Species:       {overall['n_species']}")
    print(f"Measurements:  {overall['n_measurements']}")
    print(f"Min C-value:   {overall['overall_min_pg']:.3f} pg")
    print(f"Max C-value:   {overall['overall_max_pg']:.3f} pg")
    print(f"Mean C-value:  {overall['overall_mean_pg']:.3f} pg")
    print(f"Median C-value:{overall['overall_median_pg']:.3f} pg")

    print("\nTop 10 largest genomes:")
    top = summary.sort_values("mean_c_value_pg", ascending=False).head(10)
    print(top[["species_name", "mean_c_value_pg", "mean_genome_gb"]].to_string(index=False))

    print(f"\nSaved outputs to: {outdir}")


if __name__ == "__main__":
    main()
