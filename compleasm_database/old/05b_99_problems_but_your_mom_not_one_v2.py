#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
from textwrap import dedent
from typing import Iterable

import pandas as pd


PROJECTION_DROP_COUNTS = (5, 15, 25, 35)


def normalize_bool(series: pd.Series) -> pd.Series:
    """Convert common textual/number boolean representations to True/False."""
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)

    as_str = series.astype(str).str.strip().str.lower()
    return as_str.isin({"true", "1", "yes", "y"})


def parse_keep_accessions(values: list[str] | None) -> set[str]:
    """
    Parse protected accessions supplied directly and/or through text files.

    Each existing file may contain one or more accession IDs per line. Blank
    lines and text after '#' are ignored. Multiple IDs may be separated by
    whitespace or commas.
    """
    protected: set[str] = set()

    for raw_value in values or []:
        token = raw_value.strip()
        if not token:
            continue

        candidate_path = Path(token).expanduser()
        if candidate_path.exists():
            if not candidate_path.is_file():
                raise ValueError(
                    f"--keep_accessions path is not a regular file: {candidate_path}"
                )

            with candidate_path.open("r") as handle:
                for line in handle:
                    content = line.split("#", 1)[0].strip()
                    if not content:
                        continue
                    protected.update(
                        item for item in re.split(r"[\s,]+", content) if item
                    )
            continue

        # A path-looking token is more likely to be a mistyped file than an
        # accession. Fail clearly instead of silently protecting the filename.
        if "/" in token or "\\" in token or token.lower().endswith(
            (".txt", ".tsv", ".csv")
        ):
            raise FileNotFoundError(
                f"Protected-accession file was not found: {candidate_path}"
            )

        protected.add(token)

    return protected


def prepare_pair_table(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """
    Collapse the input to one conservative record per gene-accession pair.

    This prevents an appended records_frameshift.csv from artificially
    inflating counts. For duplicate pairs, the maximum observed frameshift
    count is retained and the pair is considered problematic if any duplicate
    row was marked sequence_filtered=True.
    """
    required_cols = {"gene_id", "accession", "frameshifts"}
    missing = required_cols.difference(df.columns)
    if missing:
        raise ValueError(
            f"CSV is missing required columns: {sorted(missing)}. "
            f"Found columns: {list(df.columns)}"
        )

    clean = df.copy()
    clean["gene_id"] = clean["gene_id"].astype(str).str.strip()
    clean["accession"] = clean["accession"].astype(str).str.strip()
    clean["frameshifts"] = pd.to_numeric(clean["frameshifts"], errors="coerce")

    invalid = (
        clean["gene_id"].eq("")
        | clean["accession"].eq("")
        | clean["frameshifts"].isna()
    )
    if invalid.any():
        print(
            f"Warning: dropping {int(invalid.sum())} rows with an empty gene/accession "
            "or a non-numeric frameshift count."
        )
        clean = clean.loc[~invalid].copy()

    if clean.empty:
        raise ValueError("No valid rows remain after input validation.")

    clean["frameshifts"] = clean["frameshifts"].astype(int)
    if (clean["frameshifts"] < 0).any():
        raise ValueError("Frameshift counts must be non-negative.")

    if "sequence_filtered" in clean.columns:
        clean["problematic"] = normalize_bool(clean["sequence_filtered"])
        problem_definition = "sequence_filtered == True"
    else:
        clean["problematic"] = clean["frameshifts"] > 0
        problem_definition = "frameshifts > 0 (inferred; sequence_filtered absent)"
        print(
            "Warning: sequence_filtered is absent. For this run, a pair is treated "
            "as problematic when frameshifts > 0."
        )

    raw_valid_rows = len(clean)
    pair_table = (
        clean.groupby(["gene_id", "accession"], as_index=False)
        .agg(
            frameshifts=("frameshifts", "max"),
            problematic=("problematic", "max"),
            source_rows=("problematic", "size"),
        )
    )
    duplicates_collapsed = raw_valid_rows - len(pair_table)
    pair_table.attrs["problem_definition"] = problem_definition

    return pair_table, duplicates_collapsed


def score_recommendations(
    pair_table: pd.DataFrame,
    rejected_only: bool,
    aggressive: bool,
):
    """Build gene-level and accession-level problem summaries and rankings."""
    if rejected_only:
        working_df = pair_table[pair_table["problematic"]].copy()
    else:
        working_df = pair_table.copy()

    if working_df.empty:
        raise ValueError("No rows remain after filtering; cannot compute recommendations.")

    total_unique_genes = pair_table["gene_id"].nunique()
    total_unique_accessions = pair_table["accession"].nunique()

    # pair_table is already unique by gene-accession, so size() is a count of
    # unique pairs rather than a count inflated by appended duplicate runs.
    gene_counts = working_df.groupby("gene_id").size().rename("problem_count")
    gene_frameshifts = (
        working_df.groupby("gene_id")["frameshifts"]
        .sum()
        .rename("total_frameshifts")
    )
    gene_mean_fs = (
        working_df.groupby("gene_id")["frameshifts"]
        .mean()
        .rename("mean_frameshifts")
    )
    gene_df = pd.concat(
        [gene_counts, gene_frameshifts, gene_mean_fs], axis=1
    ).fillna(0)
    gene_df["problem_fraction"] = (
        gene_df["problem_count"] / max(total_unique_accessions, 1)
    )
    gene_df["count_rank"] = gene_df["problem_count"].rank(
        method="dense", ascending=False
    )
    gene_df["frameshift_rank"] = gene_df["total_frameshifts"].rank(
        method="dense", ascending=False
    )
    gene_df["combined_score"] = (
        gene_df["count_rank"] + gene_df["frameshift_rank"]
    )
    gene_df["severity_score"] = (
        gene_df["problem_count"] * 2
        + gene_df["total_frameshifts"]
        + gene_df["mean_frameshifts"] * 5
    )

    accession_counts = (
        working_df.groupby("accession").size().rename("problem_count")
    )
    accession_frameshifts = (
        working_df.groupby("accession")["frameshifts"]
        .sum()
        .rename("total_frameshifts")
    )
    accession_mean_fs = (
        working_df.groupby("accession")["frameshifts"]
        .mean()
        .rename("mean_frameshifts")
    )
    accession_df = pd.concat(
        [accession_counts, accession_frameshifts, accession_mean_fs], axis=1
    ).fillna(0)
    accession_df["problem_fraction"] = (
        accession_df["problem_count"] / max(total_unique_genes, 1)
    )
    accession_df["count_rank"] = accession_df["problem_count"].rank(
        method="dense", ascending=False
    )
    accession_df["frameshift_rank"] = accession_df["total_frameshifts"].rank(
        method="dense", ascending=False
    )
    accession_df["combined_score"] = (
        accession_df["count_rank"] + accession_df["frameshift_rank"]
    )
    accession_df["severity_score"] = (
        accession_df["problem_count"] * 2
        + accession_df["total_frameshifts"]
        + accession_df["mean_frameshifts"] * 5
    )

    if aggressive:
        gene_problem_taxa = (
            working_df.groupby("gene_id")["accession"]
            .nunique()
            .rename("problem_taxa")
        )
        accession_problem_genes = (
            working_df.groupby("accession")["gene_id"]
            .nunique()
            .rename("problem_genes")
        )
        gene_df = gene_df.join(gene_problem_taxa, how="left").fillna(
            {"problem_taxa": 0}
        )
        accession_df = accession_df.join(
            accession_problem_genes, how="left"
        ).fillna({"problem_genes": 0})

        gene_df["aggressive_score"] = (
            gene_df["problem_taxa"] * 3
            + gene_df["total_frameshifts"]
            + gene_df["mean_frameshifts"] * 5
        )
        accession_df["aggressive_score"] = (
            accession_df["problem_genes"] * 3
            + accession_df["total_frameshifts"]
            + accession_df["mean_frameshifts"] * 5
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

    return (
        working_df,
        gene_df,
        accession_df,
        total_unique_genes,
        total_unique_accessions,
    )


def recommend_lists(
    gene_df: pd.DataFrame,
    accession_df: pd.DataFrame,
    gene_fraction_cutoff: float,
    accession_fraction_cutoff: float,
    aggressive: bool,
    keep_accessions: set[str],
):
    """Apply recommendation thresholds while excluding protected accessions."""
    eligible_accessions = accession_df.loc[
        ~accession_df.index.isin(keep_accessions)
    ].copy()

    if aggressive:
        rec_genes = gene_df[
            (gene_df["problem_fraction"] >= min(gene_fraction_cutoff, 0.10))
            | (
                gene_df["total_frameshifts"]
                >= gene_df["total_frameshifts"].quantile(0.80)
            )
        ].copy()

        if eligible_accessions.empty:
            rec_accessions = eligible_accessions.copy()
        else:
            rec_accessions = eligible_accessions[
                (
                    eligible_accessions["problem_fraction"]
                    >= min(accession_fraction_cutoff, 0.05)
                )
                | (
                    eligible_accessions["total_frameshifts"]
                    >= eligible_accessions["total_frameshifts"].quantile(0.80)
                )
            ].copy()
    else:
        rec_genes = gene_df[
            gene_df["problem_fraction"] >= gene_fraction_cutoff
        ].copy()
        rec_accessions = eligible_accessions[
            eligible_accessions["problem_fraction"] >= accession_fraction_cutoff
        ].copy()

    # Safety net: always emit useful gene candidates. For accessions, do not
    # violate protection merely to force a non-empty output.
    if rec_genes.empty:
        rec_genes = gene_df.head(min(10, len(gene_df))).copy()
    if rec_accessions.empty and not eligible_accessions.empty:
        rec_accessions = eligible_accessions.head(
            min(10, len(eligible_accessions))
        ).copy()

    return rec_genes, rec_accessions


def count_shared_orthologs(
    usable_matrix: pd.DataFrame,
    retained_accessions: list[str],
) -> int:
    """Count genes present and usable in every retained accession."""
    if not retained_accessions:
        return 0

    retained_matrix = usable_matrix.reindex(
        columns=retained_accessions, fill_value=False
    )
    return int(retained_matrix.all(axis=1).sum())


def project_matrix_occupancy(
    pair_table: pd.DataFrame,
    accession_ranking: Iterable[str],
    keep_accessions: set[str],
    drop_counts: Iterable[int] = PROJECTION_DROP_COUNTS,
) -> pd.DataFrame:
    """
    Project strict shared-ortholog yield after dropping ranked accessions.

    A gene is shared only when every retained accession has a recorded pair and
    that pair is not problematic. Protected accessions remain in every scenario.
    """
    all_accessions = sorted(pair_table["accession"].unique())
    protected_in_data = keep_accessions.intersection(all_accessions)

    usable = pair_table.assign(usable=~pair_table["problematic"]).pivot(
        index="gene_id",
        columns="accession",
        values="usable",
    )
    usable = usable.fillna(False).astype(bool)

    ranked_candidates = [
        accession
        for accession in accession_ranking
        if accession in all_accessions and accession not in protected_in_data
    ]

    baseline_shared = count_shared_orthologs(usable, all_accessions)
    rows = []

    for requested_n in drop_counts:
        dropped = ranked_candidates[:requested_n]
        retained = [a for a in all_accessions if a not in set(dropped)]
        shared = count_shared_orthologs(usable, retained)

        rows.append(
            {
                "requested_accessions_dropped": int(requested_n),
                "actual_accessions_dropped": len(dropped),
                "retained_accessions": len(retained),
                "protected_accessions_retained": len(protected_in_data),
                "shared_orthologs": shared,
                "change_from_baseline": shared - baseline_shared,
                "baseline_shared_orthologs": baseline_shared,
                "dropped_accessions": ",".join(dropped),
            }
        )

    projection = pd.DataFrame(rows)
    projection.attrs["baseline_shared_orthologs"] = baseline_shared
    projection.attrs["total_accessions"] = len(all_accessions)
    projection.attrs["total_genes"] = pair_table["gene_id"].nunique()
    return projection


def print_table_header(title: str, formula: str, explanation: str) -> None:
    print(f"\n{'=' * 80}")
    print(title)
    print(f"Formula: {formula}")
    print(f"Meaning: {explanation}")
    print("-" * 80)


def write_list(path: Path, values: Iterable[str]) -> None:
    with path.open("w") as handle:
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
    keep_accessions: set[str],
    protected_found: set[str],
    protected_missing: set[str],
    duplicates_collapsed: int,
    problem_definition: str,
    projection: pd.DataFrame,
) -> None:
    protected_requested_text = ", ".join(sorted(keep_accessions)) or "None"
    protected_found_text = ", ".join(sorted(protected_found)) or "None"
    protected_missing_text = ", ".join(sorted(protected_missing)) or "None"

    projection_lines = []
    for row in projection.itertuples(index=False):
        projection_lines.append(
            "- Drop top {requested} unprotected accessions "
            "(actual={actual}): retained_accessions={retained}; "
            "shared_orthologs={shared}; change_from_baseline={change}".format(
                requested=row.requested_accessions_dropped,
                actual=row.actual_accessions_dropped,
                retained=row.retained_accessions,
                shared=row.shared_orthologs,
                change=row.change_from_baseline,
            )
        )

    lines = [
        "Frameshift recommendation and matrix-projection method summary",
        "=============================================================",
        "",
        "Input CSV:",
        str(csv_path),
        "",
        "Overview:",
        "This analysis used a records_frameshift.csv table to identify recurrently",
        "problematic genes and accessions. The input was collapsed to one record per",
        "gene-accession combination before calculating statistics. This is important",
        "when the CSV contains appended results from repeated runs.",
        "",
        "Duplicate handling:",
        f"- Duplicate gene-accession rows collapsed: {duplicates_collapsed}",
        "- For duplicate pairs, the maximum frameshift count was retained.",
        "- A duplicate pair was considered problematic if any source row was problematic.",
        f"- Problematic-pair definition: {problem_definition}",
        "",
        "Filtering mode:",
        f"- rejected_only = {rejected_only}",
        f"- aggressive = {aggressive}",
        "",
        "Printed statistics and equations:",
        "- Gene Problem Count(g) = number of unique problematic gene-accession pairs",
        "  associated with gene g.",
        "- Accession Problem Count(a) = number of unique problematic gene-accession",
        "  pairs associated with accession a.",
        "- Gene Total Frameshifts(g) = sum of frameshift counts across the unique",
        "  problematic pairs for gene g.",
        "- Accession Total Frameshifts(a) = sum of frameshift counts across the unique",
        "  problematic pairs for accession a.",
        "- In --all_rows mode, the same calculations use all unique pairs rather than",
        "  only problematic pairs.",
        "",
        "Recommendation logic:",
        f"- Recommended genes have problem_fraction >= {gene_fraction_cutoff:.3f}.",
        "  Gene problem_fraction = Gene Problem Count / total unique accessions.",
        f"- Recommended accessions have problem_fraction >= {accession_fraction_cutoff:.3f}.",
        "  Accession problem_fraction = Accession Problem Count / total unique genes.",
        "- Aggressive mode also admits candidates in the upper 20% of total frameshift",
        "  burden and uses a ranking that upweights problems spread across many taxa or genes.",
        "- If no gene meets the cutoff, the ten highest-ranked genes are returned.",
        "- Protected accessions are excluded from accession recommendations and from all",
        "  projected drop sets, regardless of their problem counts or fraction thresholds.",
        "",
        "Protected accessions:",
        f"- Requested: {protected_requested_text}",
        f"- Matched exactly in the CSV: {protected_found_text}",
        f"- Requested but not found: {protected_missing_text}",
        "",
        "Matrix occupancy projection:",
        "- Usable(g,a) = 1 only when the gene-accession pair exists in the CSV and is",
        "  not problematic; otherwise Usable(g,a) = 0.",
        "- SharedOrthologs(k) = number of genes g for which Usable(g,a) = 1 for every",
        "  accession a retained after dropping the top k ranked, unprotected accessions.",
        "- This is a strict 100% occupancy projection. It evaluates accession removal",
        "  only; recommended gene removals are not applied to the projection.",
        "- Because the projection is derived from records_frameshift.csv, an absent pair",
        "  is treated as unavailable. The result assumes the CSV covers the alignment",
        "  set and taxa intended for analysis.",
        "- Baseline shared orthologs before dropping accessions: "
        f"{projection.attrs.get('baseline_shared_orthologs', 'NA')}",
        "",
        "Projection results:",
        *projection_lines,
        "",
        "Dataset summary:",
        f"- Unique genes in deduplicated CSV: {total_unique_genes}",
        f"- Unique accessions in deduplicated CSV: {total_unique_accessions}",
        f"- Recommended genes written: {n_gene_recommendations}",
        f"- Recommended unprotected accessions written: {n_accession_recommendations}",
        "",
        "Practical interpretation:",
        "Gene recommendations identify loci that fail across many taxa or accumulate a",
        "high frameshift burden. Accession recommendations identify taxa that repeatedly",
        "contribute unusable sequences across loci. These are heuristic candidates for",
        "sensitivity analyses, not definitive evidence that a gene or lineage is biologically",
        "defective.",
        "",
    ]

    with (outdir / "methods_description.txt").open("w") as handle:
        handle.write("\n".join(lines))


def summarize_frameshifts(
    csv_path: Path,
    top_n: int = 20,
    rejected_only: bool = True,
    aggressive: bool = False,
    gene_fraction_cutoff: float = 0.20,
    accession_fraction_cutoff: float = 0.05,
    outdir: Path | None = None,
    keep_accessions: set[str] | None = None,
) -> None:
    if not 0 <= gene_fraction_cutoff <= 1:
        raise ValueError("--gene_fraction_cutoff must be between 0 and 1.")
    if not 0 <= accession_fraction_cutoff <= 1:
        raise ValueError("--accession_fraction_cutoff must be between 0 and 1.")
    if top_n <= 0:
        raise ValueError("--top_n must be greater than zero.")

    keep_accessions = set(keep_accessions or [])
    df = pd.read_csv(csv_path)
    pair_table, duplicates_collapsed = prepare_pair_table(df)

    if outdir is None:
        outdir = csv_path.parent
    outdir.mkdir(parents=True, exist_ok=True)

    all_accessions = set(pair_table["accession"].unique())
    protected_found = keep_accessions.intersection(all_accessions)
    protected_missing = keep_accessions.difference(all_accessions)

    print(f"Loaded frameshift log: {csv_path}")
    print(f"Raw CSV rows: {len(df)}")
    print(f"Unique gene-accession pairs: {len(pair_table)}")
    print(f"Duplicate rows collapsed: {duplicates_collapsed}")

    if "threshold" in df.columns:
        threshold_values = sorted(df["threshold"].dropna().astype(str).unique())
        if threshold_values:
            print(
                "Threshold values present in CSV: "
                + ", ".join(threshold_values)
                + "  (metadata for occupancy filtering; not used to count frameshifts)"
            )
        else:
            print("Threshold column is present but empty/NA.")

    if "errors" in df.columns:
        error_values = sorted(df["errors"].dropna().astype(str).unique())
        if error_values:
            print("Errors values present in CSV: " + ", ".join(error_values))

    if keep_accessions:
        print(
            "Protected accessions matched in CSV: "
            + (", ".join(sorted(protected_found)) or "None")
        )
        if protected_missing:
            print(
                "Warning: protected accessions not found by exact match: "
                + ", ".join(sorted(protected_missing))
            )

    (
        working_df,
        gene_df,
        accession_df,
        total_unique_genes,
        total_unique_accessions,
    ) = score_recommendations(
        pair_table=pair_table,
        rejected_only=rejected_only,
        aggressive=aggressive,
    )

    mode_word = "problematic" if rejected_only else "all"
    print(f"Pairs summarized ({mode_word} unique pairs): {len(working_df)}")
    print(f"Unique genes summarized: {working_df['gene_id'].nunique()}")
    print(f"Unique accessions summarized: {working_df['accession'].nunique()}")

    count_label = "problematic" if rejected_only else "included"

    print_table_header(
        "TOP ACCESSIONS BY PROBLEM COUNT",
        f"ProblemCount(a) = |{{g : (g,a) is a unique {count_label} pair}}|",
        f"Number of unique {count_label} genes associated with each accession.",
    )
    print(accession_df["problem_count"].head(top_n).to_string())

    print_table_header(
        "TOP GENES BY PROBLEM COUNT",
        f"ProblemCount(g) = |{{a : (g,a) is a unique {count_label} pair}}|",
        f"Number of unique {count_label} accessions associated with each gene.",
    )
    print(gene_df["problem_count"].head(top_n).to_string())

    print_table_header(
        "TOP GENES BY TOTAL FRAMESHIFTS",
        f"TotalFrameshifts(g) = sum_a frameshifts(g,a), over unique {count_label} pairs",
        f"Cumulative frameshift events for each gene across {count_label} accessions.",
    )
    print(
        gene_df["total_frameshifts"]
        .sort_values(ascending=False)
        .head(top_n)
        .to_string()
    )

    print_table_header(
        "TOP ACCESSIONS BY TOTAL FRAMESHIFTS",
        f"TotalFrameshifts(a) = sum_g frameshifts(g,a), over unique {count_label} pairs",
        f"Cumulative frameshift events for each accession across {count_label} genes.",
    )
    print(
        accession_df["total_frameshifts"]
        .sort_values(ascending=False)
        .head(top_n)
        .to_string()
    )

    protected_metrics = accession_df.loc[
        accession_df.index.intersection(sorted(protected_found))
    ]
    if not protected_metrics.empty:
        print_table_header(
            "PROTECTED ACCESSIONS (NEVER RECOMMENDED FOR DROP)",
            "Protected(a) = True => a is excluded from recommendation and projection drop sets",
            "Metrics are shown for transparency, but these accessions remain retained.",
        )
        print(
            protected_metrics[
                ["problem_count", "total_frameshifts", "problem_fraction"]
            ].to_string()
        )

    rec_genes, rec_accessions = recommend_lists(
        gene_df=gene_df,
        accession_df=accession_df,
        gene_fraction_cutoff=gene_fraction_cutoff,
        accession_fraction_cutoff=accession_fraction_cutoff,
        aggressive=aggressive,
        keep_accessions=keep_accessions,
    )

    projection = project_matrix_occupancy(
        pair_table=pair_table,
        accession_ranking=accession_df.index.tolist(),
        keep_accessions=keep_accessions,
    )

    print_table_header(
        "STRICT MATRIX OCCUPANCY PROJECTION",
        "SharedOrthologs(k) = |{g : Usable(g,a)=1 for every retained accession a after dropping top k}|",
        "Usable means the pair is present and not problematic; protected accessions are never dropped.",
    )
    projection_display = projection[
        [
            "requested_accessions_dropped",
            "actual_accessions_dropped",
            "retained_accessions",
            "shared_orthologs",
            "change_from_baseline",
        ]
    ]
    print(
        f"Baseline shared orthologs with all {projection.attrs['total_accessions']} "
        f"accessions: {projection.attrs['baseline_shared_orthologs']}"
    )
    print(projection_display.to_string(index=False))

    gene_path = outdir / "drop_genes.txt"
    accession_path = outdir / "drop_accessions.txt"
    projection_path = outdir / "matrix_occupancy_projection.tsv"

    write_list(gene_path, rec_genes.index.tolist())
    write_list(accession_path, rec_accessions.index.tolist())
    projection.to_csv(projection_path, sep="\t", index=False)

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
        keep_accessions=keep_accessions,
        protected_found=protected_found,
        protected_missing=protected_missing,
        duplicates_collapsed=duplicates_collapsed,
        problem_definition=pair_table.attrs.get("problem_definition", "Unknown"),
        projection=projection,
    )

    print_table_header(
        "RECOMMENDED DROP GENES",
        "Recommend(g) = problem_fraction(g) >= gene cutoff (plus aggressive-mode rules when enabled)",
        "Candidate loci for sensitivity analyses.",
    )
    print(
        rec_genes[
            ["problem_count", "total_frameshifts", "problem_fraction"]
        ].head(top_n).to_string()
    )

    print_table_header(
        "RECOMMENDED DROP ACCESSIONS",
        "Recommend(a) = eligible(a) AND problem_fraction(a) >= accession cutoff",
        "eligible(a) requires that accession a is not listed in --keep_accessions.",
    )
    if rec_accessions.empty:
        print("No unprotected accessions are available for recommendation.")
    else:
        print(
            rec_accessions[
                ["problem_count", "total_frameshifts", "problem_fraction"]
            ].head(top_n).to_string()
        )

    print(f"\nWrote gene recommendations to: {gene_path}")
    print(f"Wrote accession recommendations to: {accession_path}")
    print(f"Wrote occupancy projection to: {projection_path}")
    print(f"Wrote methods description to: {outdir / 'methods_description.txt'}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize problematic genes and accessions from a records_frameshift.csv "
            "file, protect specified accessions, project strict shared-ortholog yield, "
            "and write recommended drop lists."
        )
    )
    parser.add_argument(
        "csv",
        nargs="?",
        default=(
            "/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/"
            "home/leptodactylus/projects_local/projects/genomes/records/compleasm/"
            "records/records_frameshift.csv"
        ),
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
        help=(
            "Summarize all unique pairs instead of only pairs where "
            "sequence_filtered == True."
        ),
    )
    parser.add_argument(
        "--aggressive",
        action="store_true",
        help=(
            "Use a more aggressive recommendation mode that prioritizes maximizing "
            "shared genes across retained accessions under a strict occupancy objective."
        ),
    )
    parser.add_argument(
        "--gene_fraction_cutoff",
        type=float,
        default=0.20,
        help=(
            "Recommend dropping genes problematic in at least this fraction of "
            "all unique accessions."
        ),
    )
    parser.add_argument(
        "--accession_fraction_cutoff",
        type=float,
        default=0.05,
        help=(
            "Recommend dropping accessions problematic in at least this fraction "
            "of all unique genes."
        ),
    )
    parser.add_argument(
        "--keep_accessions",
        "--keep-accessions",
        nargs="+",
        default=[],
        metavar="ACCESSION_OR_FILE",
        help=(
            "Space-separated accession IDs and/or path(s) to text files containing "
            "accession IDs. Protected accessions are never recommended or projected "
            "for removal. Matching is exact."
        ),
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help=(
            "Directory where drop_genes.txt, drop_accessions.txt, "
            "matrix_occupancy_projection.tsv, and methods_description.txt are written."
        ),
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    keep_accessions = parse_keep_accessions(args.keep_accessions)

    summarize_frameshifts(
        csv_path=Path(args.csv).expanduser().resolve(),
        top_n=args.top_n,
        rejected_only=not args.all_rows,
        aggressive=args.aggressive,
        gene_fraction_cutoff=args.gene_fraction_cutoff,
        accession_fraction_cutoff=args.accession_fraction_cutoff,
        outdir=Path(args.outdir).expanduser().resolve() if args.outdir else None,
        keep_accessions=keep_accessions,
    )


if __name__ == "__main__":
    main()
