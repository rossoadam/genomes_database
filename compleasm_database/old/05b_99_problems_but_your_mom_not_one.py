#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
from typing import Iterable

import pandas as pd


PROJECTION_DROP_COUNTS = (0, 5, 15, 25, 35)


def normalize_bool(series: pd.Series) -> pd.Series:
    """Convert common textual or numeric boolean representations to bool."""
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)

    as_str = series.astype(str).str.strip().str.lower()
    return as_str.isin({"true", "1", "yes", "y"})


def parse_keep_accessions(values: list[str] | None) -> set[str]:
    """
    Parse protected accessions supplied directly and/or through text files.

    Existing files may contain one or more accession IDs per line. Blank lines
    and text following '#' are ignored. IDs may be separated by whitespace or
    commas.
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

        # A path-looking value is more likely to be a mistyped file than an
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

    This prevents an appended records_frameshift.csv from inflating counts. For
    duplicate pairs, the maximum observed frameshift count is retained, and the
    pair is considered problematic if any source row was marked filtered.
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
            "Warning: sequence_filtered is absent. A pair is being treated as "
            "problematic whenever frameshifts > 0."
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


def build_accession_completeness(pair_table: pd.DataFrame) -> pd.DataFrame:
    """Calculate recorded, missing, problematic, and usable genes per accession."""
    total_genes = pair_table["gene_id"].nunique()

    recorded = pair_table.groupby("accession")["gene_id"].nunique().rename(
        "recorded_genes"
    )
    problematic = pair_table.groupby("accession")["problematic"].sum().rename(
        "problematic_genes"
    )
    all_frameshifts = pair_table.groupby("accession")["frameshifts"].sum().rename(
        "total_frameshifts_all"
    )
    problem_frameshifts = (
        pair_table.loc[pair_table["problematic"]]
        .groupby("accession")["frameshifts"]
        .sum()
        .rename("total_frameshifts_problematic")
    )

    result = pd.concat(
        [recorded, problematic, all_frameshifts, problem_frameshifts], axis=1
    ).fillna(0)

    integer_cols = [
        "recorded_genes",
        "problematic_genes",
        "total_frameshifts_all",
        "total_frameshifts_problematic",
    ]
    result[integer_cols] = result[integer_cols].astype(int)

    result["usable_genes"] = result["recorded_genes"] - result["problematic_genes"]
    result["missing_genes"] = total_genes - result["recorded_genes"]
    result["recorded_fraction"] = result["recorded_genes"] / max(total_genes, 1)
    result["usable_fraction"] = result["usable_genes"] / max(total_genes, 1)
    result["missing_fraction"] = result["missing_genes"] / max(total_genes, 1)
    result["problem_fraction"] = result["problematic_genes"] / max(total_genes, 1)

    # Primary accession-removal ranking: least usable ortholog coverage first.
    # Missingness breaks ties before frameshift burden because the stated goal
    # is to maximize the number of genes represented across retained taxa.
    result = result.sort_values(
        [
            "usable_fraction",
            "recorded_fraction",
            "missing_genes",
            "problematic_genes",
            "total_frameshifts_problematic",
        ],
        ascending=[True, True, False, False, False],
    )
    result.attrs["total_genes"] = total_genes
    return result


def build_gene_completeness(pair_table: pd.DataFrame) -> pd.DataFrame:
    """Calculate recorded, missing, problematic, and usable taxa per gene."""
    total_accessions = pair_table["accession"].nunique()

    recorded = pair_table.groupby("gene_id")["accession"].nunique().rename(
        "recorded_accessions"
    )
    problematic = pair_table.groupby("gene_id")["problematic"].sum().rename(
        "problematic_accessions"
    )
    total_frameshifts = pair_table.groupby("gene_id")["frameshifts"].sum().rename(
        "total_frameshifts_all"
    )

    result = pd.concat([recorded, problematic, total_frameshifts], axis=1).fillna(0)
    result = result.astype(
        {
            "recorded_accessions": int,
            "problematic_accessions": int,
            "total_frameshifts_all": int,
        }
    )
    result["usable_accessions"] = (
        result["recorded_accessions"] - result["problematic_accessions"]
    )
    result["missing_accessions"] = total_accessions - result["recorded_accessions"]
    result["recorded_fraction"] = (
        result["recorded_accessions"] / max(total_accessions, 1)
    )
    result["usable_fraction"] = result["usable_accessions"] / max(total_accessions, 1)

    result = result.sort_values(
        ["usable_fraction", "recorded_fraction", "problematic_accessions"],
        ascending=[True, True, False],
    )
    result.attrs["total_accessions"] = total_accessions
    return result


def build_problem_summaries(
    pair_table: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build problem-only gene and accession summaries."""
    problem_pairs = pair_table.loc[pair_table["problematic"]].copy()
    total_genes = pair_table["gene_id"].nunique()
    total_accessions = pair_table["accession"].nunique()

    if problem_pairs.empty:
        gene_problem = pd.DataFrame(
            columns=[
                "problem_count",
                "total_frameshifts",
                "mean_frameshifts",
                "problem_fraction",
            ]
        )
        accession_problem = gene_problem.copy()
        return problem_pairs, gene_problem, accession_problem

    gene_problem = problem_pairs.groupby("gene_id").agg(
        problem_count=("accession", "nunique"),
        total_frameshifts=("frameshifts", "sum"),
        mean_frameshifts=("frameshifts", "mean"),
    )
    gene_problem["problem_fraction"] = (
        gene_problem["problem_count"] / max(total_accessions, 1)
    )
    gene_problem = gene_problem.sort_values(
        ["problem_count", "total_frameshifts", "mean_frameshifts"],
        ascending=[False, False, False],
    )

    accession_problem = problem_pairs.groupby("accession").agg(
        problem_count=("gene_id", "nunique"),
        total_frameshifts=("frameshifts", "sum"),
        mean_frameshifts=("frameshifts", "mean"),
    )
    accession_problem["problem_fraction"] = (
        accession_problem["problem_count"] / max(total_genes, 1)
    )
    accession_problem = accession_problem.sort_values(
        ["problem_count", "total_frameshifts", "mean_frameshifts"],
        ascending=[False, False, False],
    )

    return problem_pairs, gene_problem, accession_problem


def recommend_genes(
    gene_problem: pd.DataFrame,
    gene_fraction_cutoff: float,
    aggressive: bool,
) -> pd.DataFrame:
    """Recommend genes from problem prevalence and frameshift burden."""
    if gene_problem.empty:
        return gene_problem.copy()

    if aggressive:
        quantile_80 = gene_problem["total_frameshifts"].quantile(0.80)
        recommended = gene_problem.loc[
            (gene_problem["problem_fraction"] >= min(gene_fraction_cutoff, 0.10))
            | (gene_problem["total_frameshifts"] >= quantile_80)
        ].copy()
    else:
        recommended = gene_problem.loc[
            gene_problem["problem_fraction"] >= gene_fraction_cutoff
        ].copy()

    if recommended.empty:
        recommended = gene_problem.head(min(10, len(gene_problem))).copy()

    return recommended


def recommend_accessions(
    accession_completeness: pd.DataFrame,
    keep_accessions: set[str],
    completeness_cutoff: float,
    aggressive: bool,
) -> pd.DataFrame:
    """
    Recommend unprotected accessions with low usable-ortholog completeness.

    Accessions are returned least complete first. In aggressive mode, the
    effective cutoff is at least 0.95 and the bottom 20% by usable fraction are
    also included.
    """
    eligible = accession_completeness.loc[
        ~accession_completeness.index.isin(keep_accessions)
    ].copy()

    if eligible.empty:
        return eligible

    effective_cutoff = max(completeness_cutoff, 0.95) if aggressive else completeness_cutoff
    recommended = eligible.loc[eligible["usable_fraction"] < effective_cutoff].copy()

    if aggressive:
        bottom_20_cutoff = eligible["usable_fraction"].quantile(0.20)
        recommended = eligible.loc[
            (eligible["usable_fraction"] < effective_cutoff)
            | (eligible["usable_fraction"] <= bottom_20_cutoff)
        ].copy()

    if recommended.empty:
        recommended = eligible.head(min(10, len(eligible))).copy()

    return recommended


def build_usable_matrix(pair_table: pd.DataFrame) -> pd.DataFrame:
    """Create a boolean gene-by-accession matrix of usable sequences."""
    usable = pair_table.assign(usable=~pair_table["problematic"]).pivot(
        index="gene_id",
        columns="accession",
        values="usable",
    )
    # Use pandas' nullable BooleanDtype before filling missing combinations to
    # avoid object-dtype downcasting warnings in current and future pandas.
    return usable.astype("boolean").fillna(False).astype(bool)


def count_shared_orthologs(
    usable_matrix: pd.DataFrame,
    retained_accessions: list[str],
) -> int:
    """Count genes present and usable in every retained accession."""
    if not retained_accessions:
        return 0

    retained_matrix = usable_matrix.reindex(
        columns=retained_accessions,
        fill_value=False,
    )
    return int(retained_matrix.all(axis=1).sum())


def project_matrix_occupancy(
    pair_table: pd.DataFrame,
    accession_completeness: pd.DataFrame,
    keep_accessions: set[str],
    drop_counts: Iterable[int] = PROJECTION_DROP_COUNTS,
) -> pd.DataFrame:
    """
    Project strict shared-ortholog yield after dropping least-complete taxa.

    Accessions are ranked by usable_gene_fraction ascending. Missing pairs are
    unavailable, and protected accessions remain in every scenario.
    """
    all_accessions = sorted(pair_table["accession"].unique())
    protected_in_data = keep_accessions.intersection(all_accessions)
    usable_matrix = build_usable_matrix(pair_table)

    ranked_candidates = [
        accession
        for accession in accession_completeness.index
        if accession in all_accessions and accession not in protected_in_data
    ]

    baseline_shared = count_shared_orthologs(usable_matrix, all_accessions)
    rows: list[dict[str, object]] = []

    for requested_n in drop_counts:
        dropped = ranked_candidates[:requested_n]
        dropped_set = set(dropped)
        retained = [a for a in all_accessions if a not in dropped_set]
        shared = count_shared_orthologs(usable_matrix, retained)

        if retained:
            retained_completeness = accession_completeness.loc[retained]
            bottleneck_accession = retained_completeness["usable_fraction"].idxmin()
            bottleneck_usable_genes = int(
                retained_completeness.loc[bottleneck_accession, "usable_genes"]
            )
            bottleneck_usable_fraction = float(
                retained_completeness.loc[bottleneck_accession, "usable_fraction"]
            )
        else:
            bottleneck_accession = ""
            bottleneck_usable_genes = 0
            bottleneck_usable_fraction = 0.0

        rows.append(
            {
                "requested_accessions_dropped": int(requested_n),
                "actual_accessions_dropped": len(dropped),
                "retained_accessions": len(retained),
                "protected_accessions_retained": len(protected_in_data),
                "shared_orthologs": shared,
                "change_from_baseline": shared - baseline_shared,
                "baseline_shared_orthologs": baseline_shared,
                "least_complete_retained_accession": bottleneck_accession,
                "least_complete_retained_usable_genes": bottleneck_usable_genes,
                "least_complete_retained_usable_fraction": bottleneck_usable_fraction,
                "dropped_accessions": ",".join(dropped),
            }
        )

    projection = pd.DataFrame(rows)
    projection.attrs["baseline_shared_orthologs"] = baseline_shared
    projection.attrs["total_accessions"] = len(all_accessions)
    projection.attrs["total_genes"] = pair_table["gene_id"].nunique()
    projection.attrs["ranking_method"] = (
        "ascending usable_gene_fraction, then recorded_fraction; protected taxa excluded"
    )
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
    all_rows_display: bool,
    aggressive: bool,
    gene_fraction_cutoff: float,
    accession_completeness_cutoff: float,
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
            "- Drop {requested} least-complete unprotected accessions "
            "(actual={actual}): retained={retained}; shared_orthologs={shared}; "
            "change_from_baseline={change}; least_complete_retained={bottleneck} "
            "({fraction:.4f} usable).".format(
                requested=row.requested_accessions_dropped,
                actual=row.actual_accessions_dropped,
                retained=row.retained_accessions,
                shared=row.shared_orthologs,
                change=row.change_from_baseline,
                bottleneck=row.least_complete_retained_accession or "None",
                fraction=row.least_complete_retained_usable_fraction,
            )
        )

    lines = [
        "Frameshift, completeness, and matrix-projection method summary",
        "==============================================================",
        "",
        "Input CSV:",
        str(csv_path),
        "",
        "Overview:",
        "The input records_frameshift.csv was collapsed to one conservative record",
        "per gene-accession combination. Gene recommendations are based on recurrent",
        "frameshift failures. Accession recommendations and occupancy projections are",
        "based on ortholog completeness, because the goal is to remove the least-complete",
        "taxa and increase the number of orthologs shared across retained taxa.",
        "",
        "Duplicate handling:",
        f"- Duplicate gene-accession rows collapsed: {duplicates_collapsed}",
        "- For duplicates, the maximum frameshift count was retained.",
        "- A pair was problematic if any source row was problematic.",
        f"- Problematic-pair definition: {problem_definition}",
        "- The script does not independently reapply the --errors tolerance. It trusts",
        "  sequence_filtered from the cleaner when that column is available.",
        "",
        "Display mode:",
        f"- all_rows_display = {all_rows_display}",
        "- Without --all_rows, stdout emphasizes problematic pairs.",
        "- With --all_rows, stdout uses coverage-oriented labels such as recorded genes,",
        "  usable genes, missing genes, and usable fraction. Recommendation logic is",
        "  unchanged by this display option.",
        "",
        "Accession completeness equations:",
        "- RecordedGenes(a) = number of unique genes with a CSV row for accession a.",
        "- ProblematicGenes(a) = number of recorded genes marked problematic.",
        "- UsableGenes(a) = RecordedGenes(a) - ProblematicGenes(a).",
        "- MissingGenes(a) = TotalGenes - RecordedGenes(a).",
        "- UsableFraction(a) = UsableGenes(a) / TotalGenes.",
        "- Accessions are ranked least complete first by UsableFraction, then by",
        "  RecordedFraction, missing-gene count, and problematic-gene count.",
        "",
        "Gene recommendation logic:",
        f"- Recommended genes have problem_fraction >= {gene_fraction_cutoff:.3f}.",
        "- Gene problem_fraction = problematic accessions / total accessions.",
        "- Aggressive mode also includes genes in the upper 20% of total problematic",
        "  frameshift burden and lowers the prevalence cutoff to at most 0.10.",
        "",
        "Accession recommendation logic:",
        f"- Recommended accessions have usable_fraction < {accession_completeness_cutoff:.3f}.",
        "- Recommendations are sorted from least complete to most complete.",
        "- Aggressive mode uses an effective cutoff of at least 0.95 and includes the",
        "  bottom 20% of unprotected accessions by usable fraction.",
        "- Protected accessions are never recommended or projected for removal.",
        "",
        "Protected accessions:",
        f"- Requested: {protected_requested_text}",
        f"- Matched exactly in the CSV: {protected_found_text}",
        f"- Requested but not found: {protected_missing_text}",
        "",
        "Strict matrix occupancy projection:",
        "- Usable(g,a) = 1 only when the pair exists and is not problematic.",
        "- Missing pairs are treated as unavailable.",
        "- SharedOrthologs(k) = number of genes usable in every retained accession",
        "  after dropping the k least-complete unprotected accessions.",
        "- The projection uses strict 100% occupancy and does not apply the gene-drop list.",
        f"- Baseline shared orthologs: {projection.attrs.get('baseline_shared_orthologs', 'NA')}",
        "",
        "Projection results:",
        *projection_lines,
        "",
        "Dataset summary:",
        f"- Unique genes: {total_unique_genes}",
        f"- Unique accessions: {total_unique_accessions}",
        f"- Recommended genes written: {n_gene_recommendations}",
        f"- Recommended unprotected accessions written: {n_accession_recommendations}",
        "",
        "Interpretation:",
        "The accession list is designed for completeness-oriented sensitivity analyses.",
        "A low usable fraction can result from missing orthologs, filtered sequences, or",
        "both. The recommendations are not evidence that a lineage is biologically",
        "defective; they identify taxa that most constrain the shared alignment matrix.",
        "",
    ]

    with (outdir / "methods_description.txt").open("w") as handle:
        handle.write("\n".join(lines))


def summarize_frameshifts(
    csv_path: Path,
    top_n: int = 20,
    all_rows_display: bool = False,
    aggressive: bool = False,
    gene_fraction_cutoff: float = 0.20,
    accession_completeness_cutoff: float = 0.90,
    outdir: Path | None = None,
    keep_accessions: set[str] | None = None,
) -> None:
    if not 0 <= gene_fraction_cutoff <= 1:
        raise ValueError("--gene_fraction_cutoff must be between 0 and 1.")
    if not 0 <= accession_completeness_cutoff <= 1:
        raise ValueError("--accession_completeness_cutoff must be between 0 and 1.")
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

    accession_completeness = build_accession_completeness(pair_table)
    gene_completeness = build_gene_completeness(pair_table)
    problem_pairs, gene_problem, accession_problem = build_problem_summaries(pair_table)

    total_unique_genes = pair_table["gene_id"].nunique()
    total_unique_accessions = pair_table["accession"].nunique()

    print(f"Loaded frameshift log: {csv_path}")
    print(f"Raw CSV rows: {len(df)}")
    print(f"Unique gene-accession pairs: {len(pair_table)}")
    print(f"Duplicate rows collapsed: {duplicates_collapsed}")
    print(f"Unique genes: {total_unique_genes}")
    print(f"Unique accessions: {total_unique_accessions}")
    print(f"Problematic unique pairs: {len(problem_pairs)}")

    if "threshold" in df.columns:
        threshold_values = sorted(df["threshold"].dropna().astype(str).unique())
        if threshold_values:
            print(
                "Threshold values present in CSV: "
                + ", ".join(threshold_values)
                + "  (metadata for alignment occupancy; not used to count frameshifts)"
            )
        else:
            print("Threshold column is present but empty/NA.")

    if "errors" in df.columns:
        error_values = sorted(df["errors"].dropna().astype(str).unique())
        if error_values:
            print(
                "Errors values present in CSV: "
                + ", ".join(error_values)
                + "  (sequence_filtered already reflects this cleaner tolerance)"
            )

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

    if all_rows_display:
        print("Display mode: --all_rows (coverage summaries; not problem-count labels).")

        print_table_header(
            "LEAST COMPLETE ACCESSIONS",
            "UsableFraction(a) = [RecordedGenes(a) - ProblematicGenes(a)] / TotalGenes",
            "Accessions are sorted from lowest to highest usable ortholog coverage.",
        )
        print(
            accession_completeness[
                [
                    "recorded_genes",
                    "usable_genes",
                    "missing_genes",
                    "problematic_genes",
                    "usable_fraction",
                ]
            ]
            .head(top_n)
            .to_string()
        )

        print_table_header(
            "GENES WITH LOWEST USABLE TAXON COVERAGE",
            "UsableFraction(g) = [RecordedAccessions(g) - ProblematicAccessions(g)] / TotalAccessions",
            "Genes are sorted from lowest to highest usable accession coverage.",
        )
        print(
            gene_completeness[
                [
                    "recorded_accessions",
                    "usable_accessions",
                    "missing_accessions",
                    "problematic_accessions",
                    "usable_fraction",
                ]
            ]
            .head(top_n)
            .to_string()
        )

        print_table_header(
            "TOP GENES BY TOTAL OBSERVED FRAMESHIFTS (ALL RECORDED PAIRS)",
            "TotalFrameshiftsAll(g) = sum_a frameshifts(g,a)",
            "Cumulative frameshift symbols across all recorded accession pairs for each gene.",
        )
        print(
            pair_table.groupby("gene_id")["frameshifts"]
            .sum()
            .sort_values(ascending=False)
            .head(top_n)
            .to_string()
        )

        print_table_header(
            "TOP ACCESSIONS BY TOTAL OBSERVED FRAMESHIFTS (ALL RECORDED PAIRS)",
            "TotalFrameshiftsAll(a) = sum_g frameshifts(g,a)",
            "Cumulative frameshift symbols across all recorded genes for each accession.",
        )
        print(
            pair_table.groupby("accession")["frameshifts"]
            .sum()
            .sort_values(ascending=False)
            .head(top_n)
            .to_string()
        )
    else:
        print("Display mode: problematic pairs only.")

        print_table_header(
            "TOP ACCESSIONS BY PROBLEMATIC GENE COUNT",
            "ProblemCount(a) = |{g : (g,a) is a unique problematic pair}|",
            "Number of distinct genes filtered for each accession.",
        )
        if accession_problem.empty:
            print("No problematic accession-gene pairs were found.")
        else:
            print(accession_problem["problem_count"].head(top_n).to_string())

        print_table_header(
            "TOP GENES BY PROBLEMATIC ACCESSION COUNT",
            "ProblemCount(g) = |{a : (g,a) is a unique problematic pair}|",
            "Number of distinct accessions filtered for each gene.",
        )
        if gene_problem.empty:
            print("No problematic gene-accession pairs were found.")
        else:
            print(gene_problem["problem_count"].head(top_n).to_string())

        print_table_header(
            "TOP GENES BY TOTAL FRAMESHIFTS IN PROBLEMATIC PAIRS",
            "TotalFrameshifts(g) = sum_a frameshifts(g,a), for problematic pairs only",
            "Cumulative frameshift events across accessions whose sequences were filtered.",
        )
        if gene_problem.empty:
            print("No problematic gene-accession pairs were found.")
        else:
            print(
                gene_problem["total_frameshifts"]
                .sort_values(ascending=False)
                .head(top_n)
                .to_string()
            )

        print_table_header(
            "TOP ACCESSIONS BY TOTAL FRAMESHIFTS IN PROBLEMATIC PAIRS",
            "TotalFrameshifts(a) = sum_g frameshifts(g,a), for problematic pairs only",
            "Cumulative frameshift events across genes whose sequences were filtered.",
        )
        if accession_problem.empty:
            print("No problematic accession-gene pairs were found.")
        else:
            print(
                accession_problem["total_frameshifts"]
                .sort_values(ascending=False)
                .head(top_n)
                .to_string()
            )

    print_table_header(
        "LEAST COMPLETE ACCESSIONS USED FOR DROP RANKING",
        "UsableGenes(a) = RecordedGenes(a) - ProblematicGenes(a); rank ascending",
        "This ranking drives drop_accessions.txt and the occupancy projection in every display mode.",
    )
    print(
        accession_completeness[
            [
                "recorded_genes",
                "usable_genes",
                "missing_genes",
                "problematic_genes",
                "usable_fraction",
            ]
        ]
        .head(top_n)
        .to_string()
    )

    protected_metrics = accession_completeness.loc[
        accession_completeness.index.intersection(sorted(protected_found))
    ]
    if not protected_metrics.empty:
        print_table_header(
            "PROTECTED ACCESSIONS (NEVER RECOMMENDED OR PROJECTED FOR DROP)",
            "Protected(a) = True => a is excluded from all accession drop sets",
            "Completeness is shown for transparency, but these accessions remain retained.",
        )
        print(
            protected_metrics[
                [
                    "recorded_genes",
                    "usable_genes",
                    "missing_genes",
                    "problematic_genes",
                    "usable_fraction",
                ]
            ].to_string()
        )

    rec_genes = recommend_genes(
        gene_problem=gene_problem,
        gene_fraction_cutoff=gene_fraction_cutoff,
        aggressive=aggressive,
    )
    rec_accessions = recommend_accessions(
        accession_completeness=accession_completeness,
        keep_accessions=keep_accessions,
        completeness_cutoff=accession_completeness_cutoff,
        aggressive=aggressive,
    )

    projection = project_matrix_occupancy(
        pair_table=pair_table,
        accession_completeness=accession_completeness,
        keep_accessions=keep_accessions,
    )

    print_table_header(
        "STRICT MATRIX OCCUPANCY PROJECTION: DROP LEAST-COMPLETE ACCESSIONS",
        "SharedOrthologs(k) = |{g : Usable(g,a)=1 for every retained accession after dropping k least-complete taxa}|",
        "Accessions are removed by ascending usable-gene fraction; protected accessions are retained.",
    )
    projection_display = projection[
        [
            "requested_accessions_dropped",
            "actual_accessions_dropped",
            "retained_accessions",
            "shared_orthologs",
            "change_from_baseline",
            "least_complete_retained_accession",
            "least_complete_retained_usable_genes",
        ]
    ]
    print(
        f"Baseline shared orthologs with all {projection.attrs['total_accessions']} "
        f"accessions: {projection.attrs['baseline_shared_orthologs']}"
    )
    print(projection_display.to_string(index=False))

    gene_path = outdir / "drop_genes.txt"
    accession_path = outdir / "drop_accessions.txt"
    accession_completeness_path = outdir / "accession_completeness.tsv"
    projection_path = outdir / "matrix_occupancy_projection.tsv"

    write_list(gene_path, rec_genes.index.tolist())
    write_list(accession_path, rec_accessions.index.tolist())
    accession_completeness.to_csv(accession_completeness_path, sep="\t", index=True)
    projection.to_csv(projection_path, sep="\t", index=False)

    write_methods_description(
        outdir=outdir,
        csv_path=csv_path,
        all_rows_display=all_rows_display,
        aggressive=aggressive,
        gene_fraction_cutoff=gene_fraction_cutoff,
        accession_completeness_cutoff=accession_completeness_cutoff,
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
        "Recommend(g) = ProblemFraction(g) >= gene cutoff",
        "Genes remain problem-based because recurrent filtered sequences identify troublesome loci.",
    )
    if rec_genes.empty:
        print("No genes were recommended for dropping.")
    else:
        print(
            rec_genes[
                ["problem_count", "total_frameshifts", "problem_fraction"]
            ]
            .head(top_n)
            .to_string()
        )

    print_table_header(
        "RECOMMENDED DROP ACCESSIONS: LEAST COMPLETE FIRST",
        "Recommend(a) = Unprotected(a) AND UsableFraction(a) < completeness cutoff",
        "Accessions are candidates because missing and/or filtered orthologs reduce usable genome coverage.",
    )
    if rec_accessions.empty:
        print("No unprotected accessions were available for recommendation.")
    else:
        print(
            rec_accessions[
                [
                    "recorded_genes",
                    "usable_genes",
                    "missing_genes",
                    "problematic_genes",
                    "usable_fraction",
                ]
            ]
            .head(top_n)
            .to_string()
        )

    print(f"\nWrote gene recommendations to: {gene_path}")
    print(f"Wrote accession recommendations to: {accession_path}")
    print(f"Wrote accession completeness table to: {accession_completeness_path}")
    print(f"Wrote occupancy projection to: {projection_path}")
    print(f"Wrote methods description to: {outdir / 'methods_description.txt'}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize frameshift problems, rank accessions by usable ortholog "
            "completeness, protect specified taxa, project strict shared-ortholog "
            "yield, and write recommended gene/accession drop lists."
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
        help="Number of rows to print in each terminal summary.",
    )
    parser.add_argument(
        "--all_rows",
        action="store_true",
        help=(
            "Print coverage-oriented summaries for all unique gene-accession pairs. "
            "This changes stdout labels only; recommendations remain based on gene "
            "problems and accession completeness."
        ),
    )
    parser.add_argument(
        "--aggressive",
        action="store_true",
        help=(
            "Use broader gene recommendations and a stricter accession completeness "
            "cutoff (at least 0.95 plus the bottom 20 percent by usable fraction)."
        ),
    )
    parser.add_argument(
        "--gene_fraction_cutoff",
        type=float,
        default=0.20,
        help=(
            "Recommend genes problematic in at least this fraction of all accessions."
        ),
    )
    parser.add_argument(
        "--accession_completeness_cutoff",
        "--accession-completeness-cutoff",
        "--accession_fraction_cutoff",
        dest="accession_completeness_cutoff",
        type=float,
        default=0.90,
        help=(
            "Recommend unprotected accessions whose usable-gene fraction is below "
            "this value. The legacy --accession_fraction_cutoff spelling is accepted "
            "as an alias in v3."
        ),
    )
    parser.add_argument(
        "--keep_accessions",
        "--keep-accessions",
        nargs="+",
        default=[],
        metavar="ACCESSION_OR_FILE",
        help=(
            "Space-separated accession IDs and/or text files containing IDs. "
            "Protected accessions are never recommended or projected for removal. "
            "Matching is exact."
        ),
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help=(
            "Directory for drop_genes.txt, drop_accessions.txt, "
            "accession_completeness.tsv, matrix_occupancy_projection.tsv, and "
            "methods_description.txt."
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
        all_rows_display=args.all_rows,
        aggressive=args.aggressive,
        gene_fraction_cutoff=args.gene_fraction_cutoff,
        accession_completeness_cutoff=args.accession_completeness_cutoff,
        outdir=Path(args.outdir).expanduser().resolve() if args.outdir else None,
        keep_accessions=keep_accessions,
    )


if __name__ == "__main__":
    main()
