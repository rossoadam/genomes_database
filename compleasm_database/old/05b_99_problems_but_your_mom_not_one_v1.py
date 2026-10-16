#!/usr/bin/env python3

import argparse
from pathlib import Path
from textwrap import dedent
import pandas as pd


def normalize_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    as_str = series.astype(str).str.strip().str.lower()
    return as_str.isin({"true", "1", "yes", "y"})


def score_recommendations(df: pd.DataFrame, rejected_only: bool, aggressive: bool):
    working_df = df.copy()

    if rejected_only:
        if "sequence_filtered" not in working_df.columns:
            raise ValueError(
                "Rejected-only mode requires a 'sequence_filtered' column in the CSV."
            )
        working_df = working_df[normalize_bool(working_df["sequence_filtered"])]

    if working_df.empty:
        raise ValueError("No rows remain after filtering; cannot compute recommendations.")

    total_unique_genes = df["gene_id"].nunique()
    total_unique_accessions = df["accession"].nunique()

    gene_counts = working_df.groupby("gene_id").size().rename("problem_count")
    gene_frameshifts = working_df.groupby("gene_id")["frameshifts"].sum().rename("total_frameshifts")
    gene_mean_fs = working_df.groupby("gene_id")["frameshifts"].mean().rename("mean_frameshifts")
    gene_df = pd.concat([gene_counts, gene_frameshifts, gene_mean_fs], axis=1).fillna(0)
    gene_df["problem_fraction"] = gene_df["problem_count"] / max(total_unique_accessions, 1)
    gene_df["count_rank"] = gene_df["problem_count"].rank(method="dense", ascending=False)
    gene_df["frameshift_rank"] = gene_df["total_frameshifts"].rank(method="dense", ascending=False)
    gene_df["combined_score"] = gene_df["count_rank"] + gene_df["frameshift_rank"]
    gene_df["severity_score"] = (
        gene_df["problem_count"] * 2
        + gene_df["total_frameshifts"]
        + (gene_df["mean_frameshifts"] * 5)
    )

    accession_counts = working_df.groupby("accession").size().rename("problem_count")
    accession_frameshifts = working_df.groupby("accession")["frameshifts"].sum().rename("total_frameshifts")
    accession_mean_fs = working_df.groupby("accession")["frameshifts"].mean().rename("mean_frameshifts")
    accession_df = pd.concat([accession_counts, accession_frameshifts, accession_mean_fs], axis=1).fillna(0)
    accession_df["problem_fraction"] = accession_df["problem_count"] / max(total_unique_genes, 1)
    accession_df["count_rank"] = accession_df["problem_count"].rank(method="dense", ascending=False)
    accession_df["frameshift_rank"] = accession_df["total_frameshifts"].rank(method="dense", ascending=False)
    accession_df["combined_score"] = accession_df["count_rank"] + accession_df["frameshift_rank"]
    accession_df["severity_score"] = (
        accession_df["problem_count"] * 2
        + accession_df["total_frameshifts"]
        + (accession_df["mean_frameshifts"] * 5)
    )

    if aggressive:
        # Under threshold=1.0, any filtered sequence can make a gene unavailable in the full matrix.
        # Prioritize removing accessions that rescue the most genes, then genes that still block occupancy.
        gene_problem_taxa = working_df.groupby("gene_id")["accession"].nunique().rename("problem_taxa")
        accession_problem_genes = working_df.groupby("accession")["gene_id"].nunique().rename("problem_genes")
        gene_df = gene_df.join(gene_problem_taxa, how="left").fillna({"problem_taxa": 0})
        accession_df = accession_df.join(accession_problem_genes, how="left").fillna({"problem_genes": 0})
        gene_df["aggressive_score"] = (
            gene_df["problem_taxa"] * 3 + gene_df["total_frameshifts"] + gene_df["mean_frameshifts"] * 5
        )
        accession_df["aggressive_score"] = (
            accession_df["problem_genes"] * 3 + accession_df["total_frameshifts"] + accession_df["mean_frameshifts"] * 5
        )
        gene_df = gene_df.sort_values(
            ["aggressive_score", "problem_count", "total_frameshifts"],
            ascending=[False, False, False],
        )
        accession_df = accession_df.sort_values(
            ["aggressive_score", "problem_count", "total_frameshifts"],
            ascending=[False, False, False],
        )
    else:
        gene_df = gene_df.sort_values(
            ["problem_count", "total_frameshifts", "mean_frameshifts"],
            ascending=[False, False, False],
        )
        accession_df = accession_df.sort_values(
            ["problem_count", "total_frameshifts", "mean_frameshifts"],
            ascending=[False, False, False],
        )

    return working_df, gene_df, accession_df, total_unique_genes, total_unique_accessions



def recommend_lists(
    gene_df: pd.DataFrame,
    accession_df: pd.DataFrame,
    gene_fraction_cutoff: float,
    accession_fraction_cutoff: float,
    aggressive: bool,
):
    if aggressive:
        rec_genes = gene_df[
            (gene_df["problem_fraction"] >= min(gene_fraction_cutoff, 0.10))
            | (gene_df["total_frameshifts"] >= gene_df["total_frameshifts"].quantile(0.80))
        ].copy()
        rec_accessions = accession_df[
            (accession_df["problem_fraction"] >= min(accession_fraction_cutoff, 0.05))
            | (accession_df["total_frameshifts"] >= accession_df["total_frameshifts"].quantile(0.80))
        ].copy()
    else:
        rec_genes = gene_df[gene_df["problem_fraction"] >= gene_fraction_cutoff].copy()
        rec_accessions = accession_df[accession_df["problem_fraction"] >= accession_fraction_cutoff].copy()

    # Safety net: always emit something useful.
    if rec_genes.empty:
        rec_genes = gene_df.head(min(10, len(gene_df))).copy()
    if rec_accessions.empty:
        rec_accessions = accession_df.head(min(10, len(accession_df))).copy()

    return rec_genes, rec_accessions



def write_list(path: Path, values):
    with open(path, "w") as handle:
        for value in values:
            handle.write(f"{value}\n")



def write_methods_description(
    outdir: Path,
    csv_path: Path,
    rejected_only: bool,
    aggressive: bool,
    gene_fraction_cutoff: float,
    accession_fraction_cutoff: float,
    total_unique_genes: int,
    total_unique_accessions: int,
    n_gene_recommendations: int,
    n_accession_recommendations: int,
):
    methods_text = dedent(
        f"""
        Frameshift recommendation method summary
        =======================================

        Input CSV:
        {csv_path}

        Overview:
        This analysis used a records_frameshift.csv table to identify recurrently problematic genes and accessions.
        Each row corresponds to one accession-by-gene combination and includes the number of counted frameshifts.

        Filtering mode:
        - rejected_only = {rejected_only}
        - aggressive = {aggressive}

        Working assumptions:
        - When rejected_only is enabled, only rows with sequence_filtered == True are summarized.
        - Genes are ranked by how many accessions are problematic for that gene and by the total number of
          frameshifts accumulated across those problematic accessions.
        - Accessions are ranked by how many genes are problematic for that accession and by the total number
          of frameshifts accumulated across those problematic genes.

        Default recommendation logic:
        - Recommended genes are those with problem_fraction >= {gene_fraction_cutoff:.3f}, where problem_fraction is
          the fraction of accessions represented in the working table that were problematic for a given gene.
        - Recommended accessions are those with problem_fraction >= {accession_fraction_cutoff:.3f}, where problem_fraction is
          the fraction of genes represented in the working table that were problematic for a given accession.
        - If no rows meet these cutoffs, the script falls back to the top-ranked offenders so the output is still usable.

        Aggressive mode:
        - Aggressive mode prioritizes removal candidates that are most likely to increase the number of genes shared across
          all retained accessions under a strict occupancy threshold of 1.0.
        - In practice, this upweights accessions that are problematic across many genes and genes that are problematic across
          many accessions, while still considering total frameshift burden.
        - This mode is appropriate when exploring whether fewer taxa with a denser shared-gene matrix produce more stable
          downstream trends.

        Dataset summary:
        - Unique genes in CSV: {total_unique_genes}
        - Unique accessions in CSV: {total_unique_accessions}
        - Recommended genes written: {n_gene_recommendations}
        - Recommended accessions written: {n_accession_recommendations}

        Practical interpretation:
        - Gene drop recommendations identify loci that repeatedly fail across taxa and may depress matrix occupancy.
        - Accession drop recommendations identify taxa that repeatedly contribute frameshift failures across loci.
        - These recommendations are heuristic and should be treated as candidate exclusions for sensitivity analyses,
          not definitive biological claims.
        """
    ).strip() + "\n"

    with open(outdir / "methods_description.txt", "w") as handle:
        handle.write(methods_text)



def summarize_frameshifts(
    csv_path: Path,
    top_n: int = 20,
    rejected_only: bool = True,
    aggressive: bool = False,
    gene_fraction_cutoff: float = 0.20,
    accession_fraction_cutoff: float = 0.05,
    outdir: Path | None = None,
):
    df = pd.read_csv(csv_path)

    required_cols = {"gene_id", "accession", "frameshifts"}
    missing = required_cols.difference(df.columns)
    if missing:
        raise ValueError(
            f"CSV is missing required columns: {sorted(missing)}. "
            f"Found columns: {list(df.columns)}"
        )

    if outdir is None:
        outdir = csv_path.parent
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loaded frameshift log: {csv_path}")
    print(f"Total rows: {len(df)}")

    if "threshold" in df.columns:
        threshold_values = sorted(df["threshold"].dropna().astype(str).unique())
        if threshold_values:
            print(
                "Threshold values present in CSV: "
                + ", ".join(threshold_values)
                + "  (note: threshold is metadata for occupancy filtering, not frameshift counting)"
            )
        else:
            print("Threshold column is present but empty/NA.")

    if "errors" in df.columns:
        error_values = sorted(df["errors"].dropna().astype(str).unique())
        if error_values:
            print("Errors values present in CSV: " + ", ".join(error_values))

    working_df, gene_df, accession_df, total_unique_genes, total_unique_accessions = score_recommendations(
        df=df,
        rejected_only=rejected_only,
        aggressive=aggressive,
    )

    print(f"Rows summarized: {len(working_df)}")
    print(f"Unique genes summarized: {working_df['gene_id'].nunique()}")
    print(f"Unique accessions summarized: {working_df['accession'].nunique()}")

    print("\nTop accessions by count:")
    print(accession_df["problem_count"].head(top_n).to_string())

    print("\nTop genes by count:")
    print(gene_df["problem_count"].head(top_n).to_string())

    print("\nTop genes by total frameshifts:")
    print(gene_df["total_frameshifts"].sort_values(ascending=False).head(top_n).to_string())

    print("\nTop accessions by total frameshifts:")
    print(accession_df["total_frameshifts"].sort_values(ascending=False).head(top_n).to_string())

    rec_genes, rec_accessions = recommend_lists(
        gene_df=gene_df,
        accession_df=accession_df,
        gene_fraction_cutoff=gene_fraction_cutoff,
        accession_fraction_cutoff=accession_fraction_cutoff,
        aggressive=aggressive,
    )

    gene_path = outdir / "drop_genes.txt"
    accession_path = outdir / "drop_accessions.txt"
    write_list(gene_path, rec_genes.index.tolist())
    write_list(accession_path, rec_accessions.index.tolist())

    write_methods_description(
        outdir=outdir,
        csv_path=csv_path,
        rejected_only=rejected_only,
        aggressive=aggressive,
        gene_fraction_cutoff=gene_fraction_cutoff,
        accession_fraction_cutoff=accession_fraction_cutoff,
        total_unique_genes=total_unique_genes,
        total_unique_accessions=total_unique_accessions,
        n_gene_recommendations=len(rec_genes),
        n_accession_recommendations=len(rec_accessions),
    )

    print("\nRecommended drop genes:")
    print(rec_genes[["problem_count", "total_frameshifts", "problem_fraction"]].head(top_n).to_string())

    print("\nRecommended drop accessions:")
    print(rec_accessions[["problem_count", "total_frameshifts", "problem_fraction"]].head(top_n).to_string())

    print(f"\nWrote gene recommendations to: {gene_path}")
    print(f"Wrote accession recommendations to: {accession_path}")
    print(f"Wrote methods description to: {outdir / 'methods_description.txt'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Summarize problematic genes and accessions from a records_frameshift.csv file "
            "and write recommended drop lists."
        )
    )
    parser.add_argument(
        "csv",
        nargs="?",
        default="/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/projects_local/projects/genomes/records/compleasm/records/records_frameshift.csv",
        help="Path to records_frameshift.csv",
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=20,
        help="Number of top genes/accessions to print.",
    )
    parser.add_argument(
        "--all_rows",
        action="store_true",
        help="Summarize all rows instead of only rows where sequence_filtered == True.",
    )
    parser.add_argument(
        "--aggressive",
        action="store_true",
        help=(
            "Use a more aggressive recommendation mode that prioritizes maximizing shared genes across "
            "retained accessions under a strict threshold=1.0 style objective."
        ),
    )
    parser.add_argument(
        "--gene_fraction_cutoff",
        type=float,
        default=0.20,
        help="Recommend dropping genes problematic in at least this fraction of summarized accessions.",
    )
    parser.add_argument(
        "--accession_fraction_cutoff",
        type=float,
        default=0.05,
        help="Recommend dropping accessions problematic in at least this fraction of summarized genes.",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="Directory where drop_genes.txt, drop_accessions.txt, and methods_description.txt will be written.",
    )
    args = parser.parse_args()

    summarize_frameshifts(
        csv_path=Path(args.csv).resolve(),
        top_n=args.top_n,
        rejected_only=not args.all_rows,
        aggressive=args.aggressive,
        gene_fraction_cutoff=args.gene_fraction_cutoff,
        accession_fraction_cutoff=args.accession_fraction_cutoff,
        outdir=Path(args.outdir).resolve() if args.outdir else None,
    )
