#!/usr/bin/env python3
"""
04c_99_problems_but_your_mom_not_one.py

Purpose
-------
Use QC results from 04b_qualtiy_control_v2.py to identify nested genome
subsets that maximize the number of shared QC-usable orthologs.

A record is QC-usable when:
    final_qc == "passing" OR final_qc == "rescued"

The optimization uses iterative greedy genome removal:
    85 -> 84 -> ... -> 50 genomes

At each removal step, every currently retained genome is evaluated as a
candidate for removal. The genome producing the largest increase in shared
QC-usable orthologs is removed, then all candidate benefits are recalculated.

Deterministic tie-breaking
--------------------------
If multiple candidate genomes produce the same shared-ortholog gain, remove
the genome with, in order:

    1. lower total number of "passing" orthologs;
    2. higher total number of "rescued" orthologs;
    3. higher total number of "alignment-problematic" orthologs;
    4. higher total number of "unsuitable" orthologs;
    5. accession alphabetically.

The per-genome QC counts used for tie-breaking are calculated once from the
full input table and do not change as genomes are removed.

Input integrity policy
----------------------
By default, this script requires the input to form a complete rectangular
accession x odb12_id grid: every starting accession must have exactly one
record for every odb12_id represented in the input. This is intentionally
strict because missing records would otherwise make shared-ortholog counts
ambiguous.

The expected starting cohort is 85 genomes. To intentionally run a different
starting cohort, supply --expected-starting-genomes N.

Outputs
-------
Figures:
    figure1_per_genome_qc_composition.png/.pdf
    figure2_passing_vs_rescued.png/.pdf
    figure3_genome_removal_benefit.png/.pdf
    figure4_genomes_vs_shared_usable_orthologs.png/.pdf

TSVs:
    per_genome_qc_summary.tsv
    genome_removal_trajectory.tsv
    target_subset_summary.tsv
    initial_genome_removal_benefit.tsv

Subset files:
    <alignments dir>/80_2026oct14/included_accessions.txt
    <alignments dir>/80_2026oct14/dropped_accessions.txt

    <alignments dir>/70__2026oct14/included_accessions.txt
    <alignments dir>/70__2026oct14/dropped_accessions.txt

    <alignments dir>/60_2026oct14/included_accessions.txt
    <alignments dir>/60_2026oct14/dropped_accessions.txt

    <alignments dir>/50_2026oct14/included_accessions.txt
    <alignments dir>/50_2026oct14/dropped_accessions.txt
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd


# =============================================================================
# USER-EDITABLE SETTINGS
# =============================================================================

DEFAULT_INPUT = Path(
    "/Users/rossoaa/projects/genomes/records/compleasm/records/"
    "ortholog_validity/body-size-dna-dynamics_manifest/ortholog_validity.tsv"
)

DEFAULT_OUTPUT_DIR = Path(
    "/Users/rossoaa/projects/genomes/records/compleasm/records/"
    "ortholog_validity"
)

DEFAULT_ALIGNMENTS_DIR = Path(
    "/Users/rossoaa/projects/genomes/records/compleasm/alignments"
)

EXPECTED_STARTING_GENOMES = 85
TARGET_GENOME_COUNTS = (80, 70, 60, 50)

SUBSET_DIRECTORY_NAMES = {
    80: "80_2026oct14",
    70: "70__2026oct14",   # Intentionally TWO underscores, per requested path.
    60: "60_2026oct14",
    50: "50_2026oct14",
}

EXPECTED_FINAL_QC = (
    "passing",
    "rescued",
    "alignment-problematic",
    "unsuitable",
)

# Publication/output settings
PNG_DPI = 400
PDF_DPI = 400
SAVE_PDF = True

# General plot settings
AXES_LINEWIDTH = 1.8
TICK_WIDTH = 1.5
TICK_LABEL_SIZE = 10
AXIS_LABEL_SIZE = 12
TITLE_SIZE = 13
LEGEND_FONT_SIZE = 10
GRID_ALPHA = 0.20

# Figure-specific sizes in inches
FIGURE1_SIZE = (13.0, 22.0)
FIGURE2_SIZE = (10.0, 8.0)
FIGURE3_SIZE = (13.0, 22.0)
FIGURE4_SIZE = (10.0, 8.0)

# Figure 1
FIGURE1_ACCESSION_FONT_SIZE = 7
FIGURE1_BAR_HEIGHT = 0.82

# Figure 2
FIGURE2_POINT_SIZE = 60
FIGURE2_RETAINED50_POINT_SIZE = 80
FIGURE2_ANNOTATE_REMOVED = True
FIGURE2_ANNOTATION_FONT_SIZE = 5.8
FIGURE2_ANNOTATION_OFFSET = (3, 3)

# Figure 3
FIGURE3_ACCESSION_FONT_SIZE = 7
FIGURE3_BAR_HEIGHT = 0.82

# Figure 4
FIGURE4_POINT_SIZE = 28
FIGURE4_TARGET_POINT_SIZE = 90
FIGURE4_LINE_WIDTH = 2.2
FIGURE4_COMPONENT_LINE_WIDTH = 1.5

# Colors/markers are centralized here for easy editing.
QC_COLORS = {
    "passing": "#3B82F6",
    "rescued": "#10B981",
    "alignment-problematic": "#F59E0B",
    "unsuitable": "#EF4444",
}

REMOVAL_STAGE_STYLES = {
    "removed_85_to_80": dict(marker="X", color="#DC2626"),
    "removed_80_to_70": dict(marker="s", color="#EA580C"),
    "removed_70_to_60": dict(marker="D", color="#CA8A04"),
    "removed_60_to_50": dict(marker="^", color="#7C3AED"),
    "retained_50": dict(marker="o", color="#111827"),
}


REQUIRED_COLUMNS = {
    "accession",
    "odb12_id",
    "final_qc",
    "macse_qc",
    "raw_cds_qc",
    "raw_length",
    "raw_length_mod_3",
    "raw_internal_stop_count",
    "raw_terminal_stop",
    "raw_ambiguous_fraction",
    "raw_invalid_fraction",
    "raw_failure_reasons",
    "macse_frameshift",
    "macse_internal_stop",
    "macse_terminal_stop",
    "macse_failure_reasons",
    "qc_version",
}


# =============================================================================
# UTILITIES
# =============================================================================

def die(message: str) -> None:
    raise RuntimeError(message)


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def write_lines(path: Path, values: Iterable[str]) -> None:
    ensure_parent(path)
    with path.open("w", encoding="utf-8") as handle:
        for value in values:
            handle.write(f"{value}\n")


def style_axes(ax: plt.Axes) -> None:
    for spine in ax.spines.values():
        spine.set_linewidth(AXES_LINEWIDTH)
    ax.tick_params(
        axis="both",
        which="major",
        width=TICK_WIDTH,
        labelsize=TICK_LABEL_SIZE,
    )
    ax.tick_params(axis="both", which="minor", bottom=False, left=False)
    ax.grid(axis="both", alpha=GRID_ALPHA, linewidth=0.7)


def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> List[Path]:
    paths = []
    png_path = output_dir / f"{stem}.png"
    fig.savefig(png_path, dpi=PNG_DPI, bbox_inches="tight")
    paths.append(png_path)

    if SAVE_PDF:
        pdf_path = output_dir / f"{stem}.pdf"
        fig.savefig(pdf_path, dpi=PDF_DPI, bbox_inches="tight")
        paths.append(pdf_path)

    plt.close(fig)
    return paths


# =============================================================================
# INPUT + VALIDATION
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Greedily remove genomes to maximize shared QC-usable orthologs "
            "while producing nested 80/70/60/50-genome subsets."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Input ortholog_validity.tsv [default: {DEFAULT_INPUT}]",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Figure/TSV output directory [default: {DEFAULT_OUTPUT_DIR}]",
    )
    parser.add_argument(
        "--alignments-dir",
        type=Path,
        default=DEFAULT_ALIGNMENTS_DIR,
        help=f"Parent directory for subset folders [default: {DEFAULT_ALIGNMENTS_DIR}]",
    )
    parser.add_argument(
        "--expected-starting-genomes",
        type=int,
        default=EXPECTED_STARTING_GENOMES,
        help=(
            "Expected number of unique starting accessions. Default is 85. "
            "Supplying another value is the explicit override for non-85 inputs."
        ),
    )
    return parser.parse_args()


def read_and_validate_input(
    input_path: Path,
    expected_starting_genomes: int,
) -> pd.DataFrame:
    if not input_path.exists():
        die(f"Input file does not exist: {input_path}")

    if input_path.resolve() == DEFAULT_OUTPUT_DIR.resolve():
        die("Input path resolves to the output directory, not a TSV file.")

    df = pd.read_csv(
        input_path,
        sep="\t",
        dtype={"accession": "string", "odb12_id": "string", "final_qc": "string"},
        keep_default_na=False,
        na_values=[],
    )

    missing_columns = sorted(REQUIRED_COLUMNS - set(df.columns))
    if missing_columns:
        die(
            "Input is missing required column(s): "
            + ", ".join(missing_columns)
        )

    if df.empty:
        die("Input TSV contains zero rows.")

    # Strip key categorical/string identifiers.
    for col in ["accession", "odb12_id", "final_qc"]:
        df[col] = df[col].astype("string").str.strip()

    if (df["accession"] == "").any():
        die("One or more rows have a blank accession.")
    if (df["odb12_id"] == "").any():
        die("One or more rows have a blank odb12_id.")
    if (df["final_qc"] == "").any():
        die("One or more rows have a blank final_qc value.")

    # Required uniqueness: exactly one record for accession + ortholog ID.
    duplicate_mask = df.duplicated(["accession", "odb12_id"], keep=False)
    if duplicate_mask.any():
        examples = (
            df.loc[duplicate_mask, ["accession", "odb12_id"]]
            .drop_duplicates()
            .head(20)
        )
        die(
            "Duplicate accession + odb12_id combinations detected. "
            "Examples:\n"
            + examples.to_string(index=False)
        )

    observed_qc = set(df["final_qc"].unique())
    unexpected_qc = sorted(observed_qc - set(EXPECTED_FINAL_QC))
    missing_qc_categories = sorted(set(EXPECTED_FINAL_QC) - observed_qc)
    if unexpected_qc:
        die(
            "Unexpected final_qc categor"
            + ("y" if len(unexpected_qc) == 1 else "ies")
            + ": "
            + ", ".join(unexpected_qc)
            + ". Expected only: "
            + ", ".join(EXPECTED_FINAL_QC)
        )
    if missing_qc_categories:
        print(
            "[WARNING] These expected final_qc categories are absent from the "
            "input: "
            + ", ".join(missing_qc_categories),
            file=sys.stderr,
        )

    n_genomes = df["accession"].nunique()
    if n_genomes != expected_starting_genomes:
        die(
            f"Expected exactly {expected_starting_genomes} unique starting "
            f"accessions, but found {n_genomes}. "
            "Use --expected-starting-genomes N only if this is intentional."
        )

    if n_genomes < max(TARGET_GENOME_COUNTS):
        die(
            f"Starting cohort has {n_genomes} genomes, but the largest requested "
            f"target subset contains {max(TARGET_GENOME_COUNTS)} genomes."
        )
    if min(TARGET_GENOME_COUNTS) >= n_genomes:
        die(
            "The smallest target subset must contain fewer genomes than the "
            "starting cohort."
        )

    # raw_length_mod_3 is required for a requested summary metric; fail if it
    # cannot be interpreted numerically.
    raw_mod3_numeric = pd.to_numeric(df["raw_length_mod_3"], errors="coerce")
    bad_mod3 = raw_mod3_numeric.isna()
    if bad_mod3.any():
        examples = df.loc[
            bad_mod3,
            ["accession", "odb12_id", "raw_length_mod_3"],
        ].head(20)
        die(
            "raw_length_mod_3 contains missing or non-numeric values, so "
            "n_raw_not_divisible_by_3 cannot be calculated reliably. Examples:\n"
            + examples.to_string(index=False)
        )
    df["raw_length_mod_3"] = raw_mod3_numeric.astype(int)

    # Strict rectangular-grid validation.
    accessions = sorted(df["accession"].unique())
    orthologs = sorted(df["odb12_id"].unique())
    expected_rows = len(accessions) * len(orthologs)

    if len(df) != expected_rows:
        counts = (
            df.groupby("accession")["odb12_id"]
            .nunique()
            .sort_values()
        )
        die(
            "Input does not form a complete accession x odb12_id grid. "
            f"Observed {len(df):,} rows, but {len(accessions)} accessions x "
            f"{len(orthologs)} ortholog IDs requires {expected_rows:,} rows.\n"
            "Per-accession ortholog counts:\n"
            + counts.to_string()
        )

    reference_set = set(orthologs)
    inconsistent = []
    for accession, group in df.groupby("accession", sort=False):
        observed = set(group["odb12_id"])
        if observed != reference_set:
            inconsistent.append(
                (
                    accession,
                    len(reference_set - observed),
                    len(observed - reference_set),
                )
            )

    if inconsistent:
        preview = "\n".join(
            f"  {acc}: missing={n_missing}, unexpected={n_extra}"
            for acc, n_missing, n_extra in inconsistent[:20]
        )
        die(
            "Not every genome contains the same odb12_id set. "
            "Shared-ortholog calculations would be ambiguous.\n"
            + preview
        )

    return df


# =============================================================================
# SUMMARIES + MATRICES
# =============================================================================

def build_per_genome_summary(df: pd.DataFrame) -> pd.DataFrame:
    counts = (
        df.groupby(["accession", "final_qc"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=EXPECTED_FINAL_QC, fill_value=0)
    )

    summary = pd.DataFrame(index=counts.index)
    summary["n_passing"] = counts["passing"].astype(int)
    summary["n_rescued"] = counts["rescued"].astype(int)
    summary["n_usable"] = summary["n_passing"] + summary["n_rescued"]
    summary["n_alignment_problematic"] = counts["alignment-problematic"].astype(int)
    summary["n_unsuitable"] = counts["unsuitable"].astype(int)
    summary["n_total"] = counts.sum(axis=1).astype(int)

    for source_col, pct_col in [
        ("n_passing", "pct_passing"),
        ("n_rescued", "pct_rescued"),
        ("n_usable", "pct_usable"),
        ("n_alignment_problematic", "pct_alignment_problematic"),
        ("n_unsuitable", "pct_unsuitable"),
    ]:
        summary[pct_col] = 100.0 * summary[source_col] / summary["n_total"]

    raw_bad = (
        df.assign(raw_not_divisible_by_3=df["raw_length_mod_3"].ne(0))
        .groupby("accession")["raw_not_divisible_by_3"]
        .sum()
        .astype(int)
    )
    summary["n_raw_not_divisible_by_3"] = raw_bad.reindex(summary.index).astype(int)
    summary["pct_raw_not_divisible_by_3"] = (
        100.0
        * summary["n_raw_not_divisible_by_3"]
        / summary["n_total"]
    )

    return summary.reset_index()


def build_status_matrices(
    df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    usable = df["final_qc"].isin(["passing", "rescued"])
    passing = df["final_qc"].eq("passing")
    rescued = df["final_qc"].eq("rescued")

    usable_matrix = (
        df.assign(_usable=usable)
        .pivot(index="odb12_id", columns="accession", values="_usable")
        .astype(bool)
        .sort_index(axis=0)
        .sort_index(axis=1)
    )

    passing_matrix = (
        df.assign(_passing=passing)
        .pivot(index="odb12_id", columns="accession", values="_passing")
        .astype(bool)
        .reindex(index=usable_matrix.index, columns=usable_matrix.columns)
    )

    rescued_matrix = (
        df.assign(_rescued=rescued)
        .pivot(index="odb12_id", columns="accession", values="_rescued")
        .astype(bool)
        .reindex(index=usable_matrix.index, columns=usable_matrix.columns)
    )

    if usable_matrix.isna().any().any():
        die("Internal error: usable matrix contains missing values.")
    if passing_matrix.isna().any().any():
        die("Internal error: passing matrix contains missing values.")
    if rescued_matrix.isna().any().any():
        die("Internal error: rescued matrix contains missing values.")

    return usable_matrix, passing_matrix, rescued_matrix


def shared_usable_count(
    usable_matrix: pd.DataFrame,
    retained_accessions: Sequence[str],
) -> int:
    return int(usable_matrix.loc[:, list(retained_accessions)].all(axis=1).sum())


def shared_component_counts(
    usable_matrix: pd.DataFrame,
    rescued_matrix: pd.DataFrame,
    retained_accessions: Sequence[str],
) -> Tuple[int, int, int]:
    """
    Return:
        total shared QC-usable orthologs,
        shared orthologs for which every retained record is passing,
        shared orthologs that are usable only with >=1 rescued record.
    """
    retained = list(retained_accessions)
    shared_mask = usable_matrix.loc[:, retained].all(axis=1)
    total = int(shared_mask.sum())

    rescue_present = rescued_matrix.loc[:, retained].any(axis=1)
    rescue_dependent = int((shared_mask & rescue_present).sum())
    all_passing = total - rescue_dependent

    return total, all_passing, rescue_dependent


def candidate_after_removal_counts(
    usable_matrix: pd.DataFrame,
    retained_accessions: Sequence[str],
) -> Tuple[int, Dict[str, int]]:
    """
    Efficiently calculate the shared usable count before removal and the count
    after individually removing each currently retained accession.

    For each ortholog, failure_count is the number of retained genomes in which
    that ortholog is not QC-usable. Removing candidate C subtracts one from the
    failure count wherever C is unusable.
    """
    retained = list(retained_accessions)
    current = usable_matrix.loc[:, retained]

    fail_matrix = ~current
    failure_count = fail_matrix.sum(axis=1).to_numpy(dtype=np.int64)
    shared_before = int(np.count_nonzero(failure_count == 0))

    after_counts: Dict[str, int] = {}
    for accession in retained:
        candidate_fail = fail_matrix[accession].to_numpy(dtype=np.int64)
        after_failure_count = failure_count - candidate_fail
        after_counts[accession] = int(
            np.count_nonzero(after_failure_count == 0)
        )

    return shared_before, after_counts


# =============================================================================
# GREEDY REMOVAL
# =============================================================================

def tie_break_key(
    accession: str,
    gain: int,
    per_genome_lookup: pd.DataFrame,
) -> Tuple:
    row = per_genome_lookup.loc[accession]
    return (
        -gain,                                # larger gain first
        int(row["n_passing"]),                # lower passing first
        -int(row["n_rescued"]),               # higher rescued first
        -int(row["n_alignment_problematic"]), # higher alignment problems first
        -int(row["n_unsuitable"]),            # higher unsuitable first
        accession,                            # alphabetical accession
    )


def run_greedy_removal(
    usable_matrix: pd.DataFrame,
    rescued_matrix: pd.DataFrame,
    per_genome_summary: pd.DataFrame,
    final_target: int,
) -> Tuple[pd.DataFrame, Dict[int, List[str]], List[str]]:
    per_genome_lookup = per_genome_summary.set_index("accession")

    retained = sorted(usable_matrix.columns.tolist())
    removal_order: List[str] = []
    trajectory_rows = []
    subsets: Dict[int, List[str]] = {len(retained): retained.copy()}

    removal_step = 0

    while len(retained) > final_target:
        removal_step += 1
        genomes_before = len(retained)

        shared_before, after_counts = candidate_after_removal_counts(
            usable_matrix=usable_matrix,
            retained_accessions=retained,
        )

        candidate_records = []
        for accession in retained:
            shared_after = after_counts[accession]
            gain = shared_after - shared_before
            candidate_records.append((accession, gain, shared_after))

        candidate_records.sort(
            key=lambda x: tie_break_key(
                accession=x[0],
                gain=x[1],
                per_genome_lookup=per_genome_lookup,
            )
        )

        accession_removed, gain, shared_after = candidate_records[0]
        qc = per_genome_lookup.loc[accession_removed]

        retained.remove(accession_removed)
        removal_order.append(accession_removed)
        genomes_after = len(retained)

        # Component counts after the chosen removal.
        total_after, all_passing_after, rescue_dependent_after = (
            shared_component_counts(
                usable_matrix,
                rescued_matrix,
                retained,
            )
        )
        if total_after != shared_after:
            die(
                "Internal consistency error: candidate shared-after count "
                f"({shared_after}) != recomputed count ({total_after}) at "
                f"removal step {removal_step}."
            )

        trajectory_rows.append(
            {
                "removal_step": removal_step,
                "genomes_before_removal": genomes_before,
                "genomes_after_removal": genomes_after,
                "accession_removed": accession_removed,
                "passing_count_of_removed_genome": int(qc["n_passing"]),
                "rescued_count_of_removed_genome": int(qc["n_rescued"]),
                "usable_count_of_removed_genome": int(qc["n_usable"]),
                "alignment_problematic_count_of_removed_genome": int(
                    qc["n_alignment_problematic"]
                ),
                "unsuitable_count_of_removed_genome": int(qc["n_unsuitable"]),
                "shared_usable_orthologs_before": shared_before,
                "shared_usable_orthologs_after": shared_after,
                "shared_all_passing_orthologs_after": all_passing_after,
                "shared_rescue_dependent_orthologs_after": rescue_dependent_after,
                "shared_ortholog_gain": gain,
            }
        )

        subsets[genomes_after] = retained.copy()

    trajectory = pd.DataFrame(trajectory_rows)
    return trajectory, subsets, removal_order


def build_initial_removal_benefit(
    usable_matrix: pd.DataFrame,
    per_genome_summary: pd.DataFrame,
) -> pd.DataFrame:
    retained = sorted(usable_matrix.columns.tolist())
    shared_before, after_counts = candidate_after_removal_counts(
        usable_matrix=usable_matrix,
        retained_accessions=retained,
    )

    summary_lookup = per_genome_summary.set_index("accession")
    rows = []
    for accession in retained:
        shared_after = after_counts[accession]
        qc = summary_lookup.loc[accession]
        rows.append(
            {
                "accession": accession,
                "shared_usable_orthologs_all_starting_genomes": shared_before,
                "shared_usable_orthologs_after_single_removal": shared_after,
                "shared_ortholog_gain": shared_after - shared_before,
                "n_passing": int(qc["n_passing"]),
                "n_rescued": int(qc["n_rescued"]),
                "n_usable": int(qc["n_usable"]),
                "n_alignment_problematic": int(qc["n_alignment_problematic"]),
                "n_unsuitable": int(qc["n_unsuitable"]),
            }
        )

    result = pd.DataFrame(rows)
    result = result.sort_values(
        by=[
            "shared_ortholog_gain",
            "n_passing",
            "n_rescued",
            "n_alignment_problematic",
            "n_unsuitable",
            "accession",
        ],
        ascending=[False, True, False, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)

    result.insert(0, "initial_benefit_rank", np.arange(1, len(result) + 1))
    return result


# =============================================================================
# SUBSET SUMMARIES + VALIDATION
# =============================================================================

def validate_nested_subsets(
    original_accessions: Sequence[str],
    subsets: Dict[int, List[str]],
    removal_order: Sequence[str],
    target_counts: Sequence[int],
) -> None:
    original_set = set(original_accessions)
    start_n = len(original_accessions)

    for target in target_counts:
        if target not in subsets:
            die(f"Internal error: target subset {target} was not captured.")

        included = subsets[target]
        included_set = set(included)
        dropped = list(removal_order[: start_n - target])
        dropped_set = set(dropped)

        if len(included) != target:
            die(
                f"{target}-genome subset contains {len(included)} rows/accessions, "
                f"expected exactly {target}."
            )
        if len(included_set) != target:
            die(
                f"{target}-genome subset does not contain {target} unique accessions."
            )
        if included_set & dropped_set:
            die(
                f"{target}-genome included and dropped accession lists overlap."
            )
        if included_set | dropped_set != original_set:
            die(
                f"{target}-genome included + dropped sets do not exactly reconstruct "
                "the original starting accession set."
            )
        if len(dropped) != start_n - target:
            die(
                f"{target}-genome dropped list has {len(dropped)} entries; "
                f"expected {start_n - target}."
            )
        if len(dropped_set) != len(dropped):
            die(
                f"{target}-genome dropped list contains duplicate accessions."
            )

    # Explicit requested nesting validations.
    if not set(subsets[80]).issubset(original_set):
        die("80-genome subset is not contained in the starting cohort.")
    if not set(subsets[70]).issubset(set(subsets[80])):
        die("70-genome subset is not contained in the 80-genome subset.")
    if not set(subsets[60]).issubset(set(subsets[70])):
        die("60-genome subset is not contained in the 70-genome subset.")
    if not set(subsets[50]).issubset(set(subsets[60])):
        die("50-genome subset is not contained in the 60-genome subset.")


def write_subset_files(
    alignments_dir: Path,
    original_accessions: Sequence[str],
    subsets: Dict[int, List[str]],
    removal_order: Sequence[str],
) -> Dict[int, Path]:
    start_n = len(original_accessions)
    paths = {}

    for target in TARGET_GENOME_COUNTS:
        directory = alignments_dir / SUBSET_DIRECTORY_NAMES[target]
        directory.mkdir(parents=True, exist_ok=True)

        included = sorted(subsets[target])
        dropped = list(removal_order[: start_n - target])

        write_lines(directory / "included_accessions.txt", included)
        write_lines(directory / "dropped_accessions.txt", dropped)
        paths[target] = directory

    return paths


def build_target_subset_summary(
    usable_matrix: pd.DataFrame,
    rescued_matrix: pd.DataFrame,
    subsets: Dict[int, List[str]],
    starting_n: int,
) -> pd.DataFrame:
    targets_in_output = [starting_n, 80, 70, 60, 50]
    rows = []
    previous_n = starting_n

    for n_genomes in targets_in_output:
        retained = subsets[n_genomes]
        total, all_passing, rescue_dependent = shared_component_counts(
            usable_matrix,
            rescued_matrix,
            retained,
        )

        removed_from_previous_target = previous_n - n_genomes
        cumulative_removed = starting_n - n_genomes

        rows.append(
            {
                "n_genomes": n_genomes,
                "n_shared_qc_usable_orthologs": total,
                "n_shared_all_passing_orthologs": all_passing,
                "n_shared_rescue_dependent_orthologs": rescue_dependent,
                "n_genomes_removed": removed_from_previous_target,
                "cumulative_number_genomes_removed": cumulative_removed,
            }
        )
        previous_n = n_genomes

    return pd.DataFrame(rows)


def build_full_cohort_trajectory(
    usable_matrix: pd.DataFrame,
    rescued_matrix: pd.DataFrame,
    subsets: Dict[int, List[str]],
    starting_n: int,
    final_n: int,
) -> pd.DataFrame:
    rows = []
    for n_genomes in range(starting_n, final_n - 1, -1):
        retained = subsets[n_genomes]
        total, all_passing, rescue_dependent = shared_component_counts(
            usable_matrix,
            rescued_matrix,
            retained,
        )
        rows.append(
            {
                "n_genomes": n_genomes,
                "n_shared_qc_usable_orthologs": total,
                "n_shared_all_passing_orthologs": all_passing,
                "n_shared_rescue_dependent_orthologs": rescue_dependent,
            }
        )
    return pd.DataFrame(rows)


# =============================================================================
# FIGURES
# =============================================================================

def plot_figure1(
    per_genome_summary: pd.DataFrame,
    output_dir: Path,
) -> List[Path]:
    plot_df = per_genome_summary.sort_values(
        by=[
            "n_passing",
            "n_rescued",
            "n_alignment_problematic",
            "n_unsuitable",
            "accession",
        ],
        ascending=[False, True, True, True, True],
        kind="mergesort",
    ).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=FIGURE1_SIZE)

    y = np.arange(len(plot_df))
    left = np.zeros(len(plot_df), dtype=float)

    stacks = [
        ("n_passing", "Passing", "passing"),
        ("n_rescued", "Rescued", "rescued"),
        (
            "n_alignment_problematic",
            "Alignment-problematic",
            "alignment-problematic",
        ),
        ("n_unsuitable", "Unsuitable", "unsuitable"),
    ]

    for column, label, color_key in stacks:
        values = plot_df[column].to_numpy()
        ax.barh(
            y,
            values,
            left=left,
            height=FIGURE1_BAR_HEIGHT,
            label=label,
            color=QC_COLORS[color_key],
            edgecolor="none",
        )
        left += values

    ax.set_yticks(y)
    ax.set_yticklabels(
        plot_df["accession"],
        fontsize=FIGURE1_ACCESSION_FONT_SIZE,
    )
    ax.invert_yaxis()
    ax.set_xlabel("Number of evaluated orthologs", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel("Genome accession", fontsize=AXIS_LABEL_SIZE)
    ax.set_title(
        "Per-genome QC composition (ordered by passing ortholog count)",
        fontsize=TITLE_SIZE,
    )
    ax.legend(loc="best", fontsize=LEGEND_FONT_SIZE, frameon=False)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    style_axes(ax)
    fig.tight_layout()

    return save_figure(
        fig,
        output_dir,
        "figure1_per_genome_qc_composition",
    )


def assign_removal_stage(
    accession: str,
    removal_order: Sequence[str],
    starting_n: int,
) -> str:
    if accession not in removal_order:
        return "retained_50"

    step = removal_order.index(accession) + 1
    genomes_after = starting_n - step

    if genomes_after >= 80:
        return "removed_85_to_80"
    if genomes_after >= 70:
        return "removed_80_to_70"
    if genomes_after >= 60:
        return "removed_70_to_60"
    return "removed_60_to_50"


def plot_figure2(
    per_genome_summary: pd.DataFrame,
    removal_order: Sequence[str],
    output_dir: Path,
) -> List[Path]:
    plot_df = per_genome_summary.copy()
    starting_n = len(plot_df)

    plot_df["removal_stage"] = plot_df["accession"].map(
        lambda x: assign_removal_stage(x, removal_order, starting_n)
    )

    labels = {
        "removed_85_to_80": "Removed to reach 80",
        "removed_80_to_70": "Removed 80→70",
        "removed_70_to_60": "Removed 70→60",
        "removed_60_to_50": "Removed 60→50",
        "retained_50": "Retained in final 50",
    }

    stage_order = [
        "retained_50",
        "removed_60_to_50",
        "removed_70_to_60",
        "removed_80_to_70",
        "removed_85_to_80",
    ]

    fig, ax = plt.subplots(figsize=FIGURE2_SIZE)

    for stage in stage_order:
        group = plot_df.loc[plot_df["removal_stage"] == stage]
        style = REMOVAL_STAGE_STYLES[stage]
        point_size = (
            FIGURE2_RETAINED50_POINT_SIZE
            if stage == "retained_50"
            else FIGURE2_POINT_SIZE
        )
        ax.scatter(
            group["n_passing"],
            group["n_rescued"],
            s=point_size,
            marker=style["marker"],
            color=style["color"],
            alpha=0.85,
            label=labels[stage],
            linewidths=0.8,
        )

    if FIGURE2_ANNOTATE_REMOVED:
        removed_set = set(removal_order)
        for row in plot_df.itertuples(index=False):
            if row.accession in removed_set:
                ax.annotate(
                    row.accession,
                    xy=(row.n_passing, row.n_rescued),
                    xytext=FIGURE2_ANNOTATION_OFFSET,
                    textcoords="offset points",
                    fontsize=FIGURE2_ANNOTATION_FONT_SIZE,
                    alpha=0.85,
                )

    ax.set_xlabel(
        'Number of orthologs with final_qc == "passing"',
        fontsize=AXIS_LABEL_SIZE,
    )
    ax.set_ylabel(
        'Number of orthologs with final_qc == "rescued"',
        fontsize=AXIS_LABEL_SIZE,
    )
    ax.set_title(
        "Passing versus rescued orthologs per genome",
        fontsize=TITLE_SIZE,
    )
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend(
        loc="best",
        fontsize=LEGEND_FONT_SIZE,
        frameon=False,
    )
    style_axes(ax)
    fig.tight_layout()

    return save_figure(
        fig,
        output_dir,
        "figure2_passing_vs_rescued",
    )


def plot_figure3(
    initial_benefit: pd.DataFrame,
    output_dir: Path,
) -> List[Path]:
    plot_df = initial_benefit.sort_values(
        by=[
            "shared_ortholog_gain",
            "n_passing",
            "n_rescued",
            "n_alignment_problematic",
            "n_unsuitable",
            "accession",
        ],
        ascending=[False, True, False, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=FIGURE3_SIZE)

    y = np.arange(len(plot_df))
    ax.barh(
        y,
        plot_df["shared_ortholog_gain"],
        height=FIGURE3_BAR_HEIGHT,
        color="#4B5563",
        edgecolor="none",
    )

    ax.set_yticks(y)
    ax.set_yticklabels(
        plot_df["accession"],
        fontsize=FIGURE3_ACCESSION_FONT_SIZE,
    )
    ax.invert_yaxis()

    ax.set_xlabel(
        "Gain in shared QC-usable orthologs after removing one genome",
        fontsize=AXIS_LABEL_SIZE,
    )
    ax.set_ylabel("Genome accession", fontsize=AXIS_LABEL_SIZE)
    ax.set_title(
        "Initial one-genome removal benefit from the original starting cohort\n"
        "(subset optimization itself recalculates benefit after every removal)",
        fontsize=TITLE_SIZE,
    )
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    style_axes(ax)
    fig.tight_layout()

    return save_figure(
        fig,
        output_dir,
        "figure3_genome_removal_benefit",
    )


def plot_figure4(
    full_trajectory: pd.DataFrame,
    target_subset_summary: pd.DataFrame,
    output_dir: Path,
) -> List[Path]:
    fig, ax = plt.subplots(figsize=FIGURE4_SIZE)

    # Total shared QC-usable orthologs.
    ax.plot(
        full_trajectory["n_genomes"],
        full_trajectory["n_shared_qc_usable_orthologs"],
        marker="o",
        markersize=np.sqrt(FIGURE4_POINT_SIZE),
        linewidth=FIGURE4_LINE_WIDTH,
        color="#111827",
        label="Shared QC-usable total",
    )

    # Keep passing vs rescued contributions distinguishable.
    ax.plot(
        full_trajectory["n_genomes"],
        full_trajectory["n_shared_all_passing_orthologs"],
        linewidth=FIGURE4_COMPONENT_LINE_WIDTH,
        color=QC_COLORS["passing"],
        linestyle="--",
        label="Shared orthologs: all records passing",
    )
    ax.plot(
        full_trajectory["n_genomes"],
        full_trajectory["n_shared_rescue_dependent_orthologs"],
        linewidth=FIGURE4_COMPONENT_LINE_WIDTH,
        color=QC_COLORS["rescued"],
        linestyle=":",
        label="Shared orthologs requiring ≥1 rescued record",
    )

    target_only = target_subset_summary[
        target_subset_summary["n_genomes"].isin(TARGET_GENOME_COUNTS)
    ]

    ax.scatter(
        target_only["n_genomes"],
        target_only["n_shared_qc_usable_orthologs"],
        s=FIGURE4_TARGET_POINT_SIZE,
        color="#DC2626",
        marker="D",
        zorder=5,
        label="Target subsets",
    )

    for row in target_only.itertuples(index=False):
        ax.annotate(
            f"{row.n_genomes}: {row.n_shared_qc_usable_orthologs}",
            xy=(row.n_genomes, row.n_shared_qc_usable_orthologs),
            xytext=(5, 7),
            textcoords="offset points",
            fontsize=8,
        )

    ax.set_xlabel("Number of genomes retained", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel(
        "Number of orthologs shared across retained genomes",
        fontsize=AXIS_LABEL_SIZE,
    )
    ax.set_title(
        "Greedy elimination trajectory: genomes retained vs shared usable orthologs",
        fontsize=TITLE_SIZE,
    )

    # Display the trajectory in the requested direction: 85 -> 50 left-to-right.
    ax.invert_xaxis()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend(loc="best", fontsize=LEGEND_FONT_SIZE, frameon=False)
    style_axes(ax)
    fig.tight_layout()

    return save_figure(
        fig,
        output_dir,
        "figure4_genomes_vs_shared_usable_orthologs",
    )


# =============================================================================
# MAIN
# =============================================================================

def main() -> int:
    args = parse_args()

    input_path = args.input.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    alignments_dir = args.alignments_dir.expanduser().resolve()

    output_dir.mkdir(parents=True, exist_ok=True)
    alignments_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Reading: {input_path}")
    df = read_and_validate_input(
        input_path=input_path,
        expected_starting_genomes=args.expected_starting_genomes,
    )

    original_accessions = sorted(df["accession"].unique())
    starting_n = len(original_accessions)
    final_target = min(TARGET_GENOME_COUNTS)

    print(
        f"[INFO] Input validated: {starting_n} genomes, "
        f"{df['odb12_id'].nunique():,} ortholog IDs, "
        f"{len(df):,} accession/ortholog records."
    )

    qc_versions = sorted(
        {
            str(x).strip()
            for x in df["qc_version"].tolist()
            if str(x).strip()
        }
    )
    if qc_versions:
        print("[INFO] QC version value(s): " + ", ".join(qc_versions))

    per_genome_summary = build_per_genome_summary(df)
    usable_matrix, passing_matrix, rescued_matrix = build_status_matrices(df)

    # passing_matrix is intentionally constructed/validated even though the
    # shared component calculation can infer all-passing as total minus
    # rescue-dependent. This ensures the categorical matrices remain internally
    # auditable.
    if not (
        (passing_matrix | rescued_matrix).equals(usable_matrix)
    ):
        die(
            "Internal consistency error: passing OR rescued matrix does not "
            "equal the QC-usable matrix."
        )

    starting_shared, starting_all_passing, starting_rescue_dependent = (
        shared_component_counts(
            usable_matrix,
            rescued_matrix,
            original_accessions,
        )
    )

    print(
        f"[INFO] Starting shared QC-usable orthologs: {starting_shared:,} "
        f"(all-passing={starting_all_passing:,}; "
        f"rescue-dependent={starting_rescue_dependent:,})"
    )

    initial_benefit = build_initial_removal_benefit(
        usable_matrix=usable_matrix,
        per_genome_summary=per_genome_summary,
    )

    trajectory, subsets, removal_order = run_greedy_removal(
        usable_matrix=usable_matrix,
        rescued_matrix=rescued_matrix,
        per_genome_summary=per_genome_summary,
        final_target=final_target,
    )

    # Explicit validation of nested target subsets.
    validate_nested_subsets(
        original_accessions=original_accessions,
        subsets=subsets,
        removal_order=removal_order,
        target_counts=TARGET_GENOME_COUNTS,
    )

    target_subset_summary = build_target_subset_summary(
        usable_matrix=usable_matrix,
        rescued_matrix=rescued_matrix,
        subsets=subsets,
        starting_n=starting_n,
    )

    full_trajectory = build_full_cohort_trajectory(
        usable_matrix=usable_matrix,
        rescued_matrix=rescued_matrix,
        subsets=subsets,
        starting_n=starting_n,
        final_n=final_target,
    )

    # Additional trajectory integrity checks.
    expected_sizes = list(range(starting_n, final_target - 1, -1))
    observed_sizes = full_trajectory["n_genomes"].tolist()
    if observed_sizes != expected_sizes:
        die(
            "Full trajectory does not contain exactly one cohort for every "
            f"size from {starting_n} down to {final_target}."
        )

    # Write TSVs.
    per_genome_path = output_dir / "per_genome_qc_summary.tsv"
    trajectory_path = output_dir / "genome_removal_trajectory.tsv"
    target_summary_path = output_dir / "target_subset_summary.tsv"
    initial_benefit_path = output_dir / "initial_genome_removal_benefit.tsv"

    per_genome_summary.to_csv(per_genome_path, sep="\t", index=False)
    trajectory.to_csv(trajectory_path, sep="\t", index=False)
    target_subset_summary.to_csv(target_summary_path, sep="\t", index=False)
    initial_benefit.to_csv(initial_benefit_path, sep="\t", index=False)

    subset_paths = write_subset_files(
        alignments_dir=alignments_dir,
        original_accessions=original_accessions,
        subsets=subsets,
        removal_order=removal_order,
    )

    # Post-write validation of subset file counts and contents.
    for target, directory in subset_paths.items():
        included_path = directory / "included_accessions.txt"
        dropped_path = directory / "dropped_accessions.txt"

        included = [
            x.strip()
            for x in included_path.read_text(encoding="utf-8").splitlines()
            if x.strip()
        ]
        dropped = [
            x.strip()
            for x in dropped_path.read_text(encoding="utf-8").splitlines()
            if x.strip()
        ]

        expected_dropped = starting_n - target
        if len(included) != target:
            die(
                f"Post-write validation failed: {included_path} contains "
                f"{len(included)} accessions, expected {target}."
            )
        if len(dropped) != expected_dropped:
            die(
                f"Post-write validation failed: {dropped_path} contains "
                f"{len(dropped)} accessions, expected {expected_dropped}."
            )
        if set(included) & set(dropped):
            die(
                f"Post-write validation failed: overlap detected in {directory}."
            )
        if set(included) | set(dropped) != set(original_accessions):
            die(
                f"Post-write validation failed: included+dropped does not "
                f"reconstruct starting accessions in {directory}."
            )

    # Generate figures.
    figure_paths: List[Path] = []
    figure_paths += plot_figure1(per_genome_summary, output_dir)
    figure_paths += plot_figure2(
        per_genome_summary,
        removal_order,
        output_dir,
    )
    figure_paths += plot_figure3(initial_benefit, output_dir)
    figure_paths += plot_figure4(
        full_trajectory,
        target_subset_summary,
        output_dir,
    )

    # Final terminal summary.
    print("\n" + "=" * 78)
    print("GREEDY GENOME-REMOVAL SUMMARY")
    print("=" * 78)
    print(f"Starting genomes: {starting_n}")
    print(f"Starting shared QC-usable orthologs: {starting_shared:,}")
    print("\nRemoval order:")

    for row in trajectory.itertuples(index=False):
        print(
            f"  step {row.removal_step:2d}: "
            f"{row.accession_removed} | "
            f"{row.genomes_before_removal}->{row.genomes_after_removal} genomes | "
            f"gain={row.shared_ortholog_gain:+d} | "
            f"shared={row.shared_usable_orthologs_after:,}"
        )

    print("\nTarget cohort shared QC-usable ortholog counts:")
    for row in target_subset_summary.itertuples(index=False):
        print(
            f"  {row.n_genomes:2d} genomes: "
            f"{row.n_shared_qc_usable_orthologs:,} total "
            f"({row.n_shared_all_passing_orthologs:,} all-passing; "
            f"{row.n_shared_rescue_dependent_orthologs:,} rescue-dependent)"
        )

    print("\nGenerated TSVs:")
    for path in [
        per_genome_path,
        trajectory_path,
        target_summary_path,
        initial_benefit_path,
    ]:
        print(f"  {path}")

    print("\nGenerated figures:")
    for path in figure_paths:
        print(f"  {path}")

    print("\nSubset directories:")
    for target in TARGET_GENOME_COUNTS:
        print(f"  {target}: {subset_paths[target]}")

    print("\n[DONE] All validation checks passed.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
