#!/usr/bin/env python3
"""
04c_99_problems_but_your_mom_not_one_v4.py

Purpose
-------
Use QC results from 04b_qualtiy_control_v2.py to identify nested genome
subsets that maximize the number of shared QC-usable orthologs.

The first positional command-line argument is the genomes project root, for
example:

    /Users/rossoaa/projects/genomes

From that root, the script derives:

    <genomes_dir>/records/compleasm
    <genomes_dir>/records/compleasm/records/ortholog_validity
    <genomes_dir>/records/compleasm/alignments

and, by default, reads:

    <genomes_dir>/records/compleasm/records/ortholog_validity/
        body-size-dna-dynamics_manifest/ortholog_validity.tsv

The --input, --output-dir, and --alignments-dir options remain available as
explicit overrides, but normally only the genomes root needs to be supplied.

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

Input occupancy and integrity policy
------------------------------------
The upstream candidate ortholog pool was generated at a 50% occupancy
threshold across the starting genomes so that genome subsetting can recover
larger fully shared ortholog sets. Therefore the input is intentionally sparse:
not every accession is expected to contain every odb12_id.

An absent accession x odb12_id combination is treated as NOT QC-usable for
shared-ortholog optimization, but absence remains distinct from QC failure in
descriptive summaries. Missing pairs are never relabeled as "unsuitable" or
"alignment-problematic".

Protected genomes
-----------------
Specific accessions may be protected from greedy removal with:

    --protected-accessions ACCESSION1 ACCESSION2 ...

Protected genomes remain in every retained cohort for biological or
phylogenetic reasons. They are excluded from the removal candidate set, but
they still participate fully in shared-ortholog calculations. Therefore an
ortholog must still be present and QC-usable in every protected genome to
count as shared within a retained subset.

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
    <alignments dir>/80_2026oct14/shared_orthologs.txt

    <alignments dir>/70_2026oct14/included_accessions.txt
    <alignments dir>/70_2026oct14/dropped_accessions.txt
    <alignments dir>/70_2026oct14/shared_orthologs.txt

    <alignments dir>/60_2026oct14/included_accessions.txt
    <alignments dir>/60_2026oct14/dropped_accessions.txt
    <alignments dir>/60_2026oct14/shared_orthologs.txt

    <alignments dir>/50_2026oct14/included_accessions.txt
    <alignments dir>/50_2026oct14/dropped_accessions.txt
    <alignments dir>/50_2026oct14/shared_orthologs.txt

For downstream alignment filtering, included_accessions.txt and
shared_orthologs.txt are the authoritative cohort-selection documents intended
for a later 05_cleaner_v6.py operating on 02_macse_alignments_with_i_braminus.
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
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


# =============================================================================
# USER-EDITABLE SETTINGS
# =============================================================================

# Paths are derived from the required genomes_dir positional argument.
# Example genomes_dir:
#     /Users/rossoaa/projects/genomes
#
# Derived compleasm root:
#     <genomes_dir>/records/compleasm
DEFAULT_INPUT_RELATIVE = Path(
    "records/ortholog_validity/"
    "body-size-dna-dynamics_manifest/ortholog_validity.tsv"
)

DEFAULT_OUTPUT_RELATIVE = Path("records/ortholog_validity")
DEFAULT_ALIGNMENTS_RELATIVE = Path("alignments")

EXPECTED_STARTING_GENOMES = 85
TARGET_GENOME_COUNTS = (80, 70, 60, 50)

SUBSET_DIRECTORY_NAMES = {
    80: "80_2026oct14",
    70: "70_2026oct14", 
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
FIGURE2_PROTECTED_POINT_SIZE = 180
FIGURE2_PROTECTED_EDGE_WIDTH = 2.2
FIGURE2_PROTECTED_MARKER = "o"
FIGURE2_ANNOTATE_PROTECTED = True
FIGURE2_PROTECTED_ANNOTATION_FONT_SIZE = 7
FIGURE2_PROTECTED_ANNOTATION_OFFSET = (5, -8)

# Figure 3
FIGURE3_ACCESSION_FONT_SIZE = 7
FIGURE3_BAR_HEIGHT = 0.82
FIGURE3_DEFAULT_BAR_COLOR = "#4B5563"
FIGURE3_PROTECTED_BAR_COLOR = "#111827"
FIGURE3_PROTECTED_HATCH = "///"
FIGURE3_PROTECTED_EDGE_WIDTH = 1.0

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
    "absent": "#9CA3AF",
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
            "while producing nested 80/70/60/50-genome subsets. The required "
            "genomes_dir positional argument establishes the project root."
        )
    )
    parser.add_argument(
        "genomes_dir",
        type=Path,
        help=(
            "Root genomes directory, e.g. /Users/rossoaa/projects/genomes. "
            "The script derives <genomes_dir>/records/compleasm and associated "
            "input/output directories from this root."
        ),
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help=(
            "Optional input TSV override. Default: "
            "<genomes_dir>/records/compleasm/"
            + str(DEFAULT_INPUT_RELATIVE)
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Optional figure/TSV output directory override. Default: "
            "<genomes_dir>/records/compleasm/"
            + str(DEFAULT_OUTPUT_RELATIVE)
        ),
    )
    parser.add_argument(
        "--alignments-dir",
        type=Path,
        default=None,
        help=(
            "Optional parent directory for subset folders. Default: "
            "<genomes_dir>/records/compleasm/"
            + str(DEFAULT_ALIGNMENTS_RELATIVE)
        ),
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
    parser.add_argument(
        "--protected-accessions",
        nargs="+",
        default=[],
        metavar="ACCESSION",
        help=(
            "Optional one-or-more accession IDs that must remain in every "
            "retained subset. Protected genomes participate normally in shared-"
            "ortholog calculations but are never eligible for greedy removal."
        ),
    )
    return parser.parse_args()


def validate_protected_accessions(
    requested: Sequence[str],
    starting_accessions: Sequence[str],
    final_target: int,
) -> set[str]:
    """
    Normalize and validate protected accessions.

    Protected genomes are excluded from the greedy removal candidate set but
    remain in all shared-ortholog calculations.
    """
    normalized_list: List[str] = []
    for raw_value in requested:
        value = str(raw_value).strip()
        if not value:
            die(
                "Protected accession values must be nonblank. "
                "Remove empty values from --protected-accessions."
            )
        normalized_list.append(value)

    protected = set(normalized_list)
    starting_set = set(starting_accessions)

    missing = sorted(protected - starting_set)
    if missing:
        die(
            "Requested protected accession(s) are absent from the starting "
            "cohort: "
            + ", ".join(missing)
        )

    if len(protected) > final_target:
        die(
            f"{len(protected)} accessions are protected, but the final target "
            f"contains only {final_target} genomes. The requested optimization "
            "is impossible."
        )

    return protected


def read_and_validate_input(
    input_path: Path,
    expected_starting_genomes: int,
) -> pd.DataFrame:
    if not input_path.exists():
        die(f"Input file does not exist: {input_path}")

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

    # Sparse occupancy is expected. The upstream candidate ortholog pool was
    # generated using a 50% occupancy threshold so not every accession should
    # be expected to contain every odb12_id.
    accessions = sorted(df["accession"].unique())
    orthologs = sorted(df["odb12_id"].unique())

    if not orthologs:
        die("Input contains zero candidate ortholog IDs.")

    possible_pairs = len(accessions) * len(orthologs)
    observed_pairs = len(df)
    absent_pairs = possible_pairs - observed_pairs

    if absent_pairs < 0:
        die(
            "Internal integrity error: observed accession/ortholog rows exceed "
            "the maximum possible unique accession x odb12_id combinations."
        )

    per_accession_present = df.groupby("accession")["odb12_id"].nunique()
    if (per_accession_present <= 0).any():
        die("One or more starting accessions contain zero ortholog records.")
    if (per_accession_present > len(orthologs)).any():
        die(
            "One or more accessions contain more unique ortholog records than "
            "the candidate ortholog pool size."
        )

    occupancy_pct = (
        100.0 * observed_pairs / possible_pairs if possible_pairs else 0.0
    )
    print(
        "[INFO] Sparse accession x ortholog matrix accepted: "
        f"{observed_pairs:,}/{possible_pairs:,} observed pairs "
        f"({occupancy_pct:.2f}% occupancy); "
        f"{absent_pairs:,} pairs are absent and will be treated as not usable."
    )

    return df


# =============================================================================
# SUMMARIES + MATRICES
# =============================================================================

def build_per_genome_summary(
    df: pd.DataFrame,
    candidate_ortholog_count: int,
    protected_accessions: set[str],
) -> pd.DataFrame:
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

    # n_total remains the number of PRESENT/EVALUATED records so all QC
    # percentages retain their original interpretation.
    summary["n_total"] = counts.sum(axis=1).astype(int)

    summary["n_orthologs_in_candidate_pool"] = int(candidate_ortholog_count)
    summary["n_ortholog_records_present"] = summary["n_total"].astype(int)
    summary["n_ortholog_records_absent"] = (
        candidate_ortholog_count - summary["n_ortholog_records_present"]
    ).astype(int)

    if (summary["n_ortholog_records_absent"] < 0).any():
        die(
            "Internal integrity error: a genome contains more present ortholog "
            "records than the candidate ortholog pool size."
        )

    summary["pct_ortholog_occupancy"] = (
        100.0
        * summary["n_ortholog_records_present"]
        / candidate_ortholog_count
    )

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

    summary["protected"] = summary.index.to_series().isin(protected_accessions)

    return summary.reset_index()


def build_status_matrices(
    df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Build complete candidate-ortholog x accession boolean matrices.

    Sparse input cells are EXPECTED because the upstream candidate pool was
    defined at 50% occupancy. Missing accession x odb12_id combinations are
    explicitly filled with False in all matrices.

    usable=True only for observed passing or rescued records.
    """
    accessions = sorted(df["accession"].unique())
    orthologs = sorted(df["odb12_id"].unique())

    usable = df["final_qc"].isin(["passing", "rescued"])
    passing = df["final_qc"].eq("passing")
    rescued = df["final_qc"].eq("rescued")

    def make_matrix(mask: pd.Series, value_name: str) -> pd.DataFrame:
        matrix = (
            df.assign(**{value_name: mask})
            .pivot(index="odb12_id", columns="accession", values=value_name)
            .reindex(index=orthologs, columns=accessions)
            .astype("boolean")
            .fillna(False)
            .astype(bool)
        )
        return matrix

    usable_matrix = make_matrix(usable, "_usable")
    passing_matrix = make_matrix(passing, "_passing")
    rescued_matrix = make_matrix(rescued, "_rescued")

    if usable_matrix.isna().any().any():
        die("Internal error: usable matrix contains missing values.")
    if passing_matrix.isna().any().any():
        die("Internal error: passing matrix contains missing values.")
    if rescued_matrix.isna().any().any():
        die("Internal error: rescued matrix contains missing values.")

    # Observed records must map exactly according to final_qc.
    observed_lookup = df.set_index(["odb12_id", "accession"])["final_qc"]
    for (odb12_id, accession), final_qc in observed_lookup.items():
        observed_usable = bool(usable_matrix.at[odb12_id, accession])
        expected_usable = final_qc in {"passing", "rescued"}
        if observed_usable != expected_usable:
            die(
                "Internal matrix validation failed for "
                f"{accession} + {odb12_id}: final_qc={final_qc!r}, "
                f"usable_matrix={observed_usable}."
            )

    # Every absent pair must be False in every status matrix.
    observed_index = set(zip(df["odb12_id"], df["accession"]))
    for odb12_id in orthologs:
        for accession in accessions:
            if (odb12_id, accession) not in observed_index:
                if (
                    bool(usable_matrix.at[odb12_id, accession])
                    or bool(passing_matrix.at[odb12_id, accession])
                    or bool(rescued_matrix.at[odb12_id, accession])
                ):
                    die(
                        "Internal matrix validation failed: absent pair "
                        f"{accession} + {odb12_id} was not encoded as False."
                    )

    return usable_matrix, passing_matrix, rescued_matrix


def shared_usable_count(
    usable_matrix: pd.DataFrame,
    retained_accessions: Sequence[str],
) -> int:
    return int(usable_matrix.loc[:, list(retained_accessions)].all(axis=1).sum())


def shared_usable_ortholog_ids(
    usable_matrix: pd.DataFrame,
    retained_accessions: Sequence[str],
) -> List[str]:
    """
    Return deterministically sorted ortholog IDs that are QC-usable in every
    retained accession.

    The usable_matrix already encodes:
        passing/rescued -> True
        alignment-problematic/unsuitable/absent -> False
    """
    retained = list(retained_accessions)
    if not retained:
        die("Cannot calculate shared orthologs for an empty retained cohort.")

    shared_mask = usable_matrix.loc[:, retained].all(axis=1)
    orthologs = usable_matrix.index[shared_mask].astype(str).tolist()
    return sorted(orthologs)


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
    protected_accessions: set[str],
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

        # Protected genomes remain in `retained`, so they continue to constrain
        # which orthologs are shared. They are simply excluded from the removal
        # candidate set.
        removable_candidates = [
            accession
            for accession in retained
            if accession not in protected_accessions
        ]

        if not removable_candidates:
            die(
                "Greedy removal cannot continue because no removable genomes "
                f"remain at cohort size {genomes_before}, but the final target "
                f"is {final_target}."
            )

        candidate_records = []
        for accession in removable_candidates:
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
                "absent_count_of_removed_genome": int(
                    qc["n_ortholog_records_absent"]
                ),
                "occupancy_pct_of_removed_genome": float(
                    qc["pct_ortholog_occupancy"]
                ),
                "protected_accession_removed": bool(
                    accession_removed in protected_accessions
                ),
                "shared_usable_orthologs_before": shared_before,
                "shared_usable_orthologs_after": shared_after,
                "shared_all_passing_orthologs_after": all_passing_after,
                "shared_rescue_dependent_orthologs_after": rescue_dependent_after,
                "shared_ortholog_gain": gain,
            }
        )

        subsets[genomes_after] = retained.copy()

    trajectory = pd.DataFrame(trajectory_rows)

    if set(removal_order) & protected_accessions:
        bad = sorted(set(removal_order) & protected_accessions)
        die(
            "Internal protection error: protected accession(s) were removed: "
            + ", ".join(bad)
        )

    if (
        not trajectory.empty
        and trajectory["protected_accession_removed"].astype(bool).any()
    ):
        die(
            "Internal protection error: trajectory indicates that a protected "
            "accession was removed."
        )

    return trajectory, subsets, removal_order


def build_initial_removal_benefit(
    usable_matrix: pd.DataFrame,
    per_genome_summary: pd.DataFrame,
    protected_accessions: set[str],
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
                "n_ortholog_records_absent": int(
                    qc["n_ortholog_records_absent"]
                ),
                "pct_ortholog_occupancy": float(qc["pct_ortholog_occupancy"]),
                "protected": bool(accession in protected_accessions),
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
    protected_accessions: set[str],
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

    removed_protected = set(removal_order) & protected_accessions
    if removed_protected:
        die(
            "Protected accession(s) appear in removal_order: "
            + ", ".join(sorted(removed_protected))
        )

    for target in target_counts:
        included_set = set(subsets[target])
        missing_protected = protected_accessions - included_set
        if missing_protected:
            die(
                f"Protected accession(s) missing from the {target}-genome subset: "
                + ", ".join(sorted(missing_protected))
            )

        dropped_set = set(removal_order[: start_n - target])
        dropped_protected = protected_accessions & dropped_set
        if dropped_protected:
            die(
                f"Protected accession(s) appear in the {target}-genome dropped set: "
                + ", ".join(sorted(dropped_protected))
            )

    missing_from_final = protected_accessions - set(subsets[min(target_counts)])
    if missing_from_final:
        die(
            "Protected accession(s) are absent from the final target subset: "
            + ", ".join(sorted(missing_from_final))
        )


def write_subset_files(
    alignments_dir: Path,
    original_accessions: Sequence[str],
    subsets: Dict[int, List[str]],
    removal_order: Sequence[str],
    usable_matrix: pd.DataFrame,
) -> Tuple[Dict[int, Path], Dict[int, List[str]]]:
    """
    Write the authoritative downstream cohort-selection files for each target.

    Each target directory contains:
        included_accessions.txt
        dropped_accessions.txt
        shared_orthologs.txt

    shared_orthologs.txt contains exactly the orthologs that are QC-usable in
    100% of the included accessions for that target subset.
    """
    start_n = len(original_accessions)
    paths: Dict[int, Path] = {}
    shared_by_target: Dict[int, List[str]] = {}

    for target in TARGET_GENOME_COUNTS:
        directory = alignments_dir / SUBSET_DIRECTORY_NAMES[target]
        directory.mkdir(parents=True, exist_ok=True)

        included = sorted(subsets[target])
        dropped = list(removal_order[: start_n - target])
        shared_orthologs = shared_usable_ortholog_ids(
            usable_matrix=usable_matrix,
            retained_accessions=included,
        )

        write_lines(directory / "included_accessions.txt", included)
        write_lines(directory / "dropped_accessions.txt", dropped)
        write_lines(directory / "shared_orthologs.txt", shared_orthologs)

        paths[target] = directory
        shared_by_target[target] = shared_orthologs

    return paths, shared_by_target


def validate_shared_ortholog_files(
    usable_matrix: pd.DataFrame,
    subsets: Dict[int, List[str]],
    subset_paths: Dict[int, Path],
    shared_by_target: Dict[int, List[str]],
    target_subset_summary: pd.DataFrame,
) -> None:
    """
    Validate target-specific shared_orthologs.txt files against the exact
    100%-usable matrix mask, target summary counts, uniqueness, and expected
    nesting across progressively smaller genome subsets.
    """
    summary_lookup = target_subset_summary.set_index("n_genomes")
    validated_sets: Dict[int, set[str]] = {}

    for target in TARGET_GENOME_COUNTS:
        directory = subset_paths[target]
        shared_path = directory / "shared_orthologs.txt"
        included = sorted(subsets[target])

        if not shared_path.exists():
            die(f"Expected shared ortholog file was not written: {shared_path}")

        written = [
            x.strip()
            for x in shared_path.read_text(encoding="utf-8").splitlines()
            if x.strip()
        ]

        if len(written) != len(set(written)):
            die(
                f"Duplicate ortholog IDs detected in {shared_path}: "
                f"{len(written)} lines but {len(set(written))} unique IDs."
            )

        expected = shared_usable_ortholog_ids(
            usable_matrix=usable_matrix,
            retained_accessions=included,
        )

        # File content must exactly match the matrix-derived set and order.
        if written != expected:
            written_set = set(written)
            expected_set = set(expected)
            missing = sorted(expected_set - written_set)[:20]
            extra = sorted(written_set - expected_set)[:20]
            die(
                f"{shared_path} does not exactly match the 100%-QC-usable "
                f"ortholog set for the {target}-genome cohort. "
                f"Missing examples={missing}; extra examples={extra}."
            )

        # Also verify against the list returned during writing.
        if written != shared_by_target[target]:
            die(
                f"Internal consistency error: written ortholog list for {target} "
                "genomes differs from the list retained in memory."
            )

        # Every written ortholog must be usable in all included accessions.
        if written:
            matrix_check = usable_matrix.loc[written, included].all(axis=1)
            if not bool(matrix_check.all()):
                bad = matrix_check.loc[~matrix_check].index.astype(str).tolist()[:20]
                die(
                    f"Matrix validation failed for {shared_path}. These orthologs "
                    f"are not QC-usable in every included accession: {bad}"
                )

        if target not in summary_lookup.index:
            die(
                f"target_subset_summary.tsv logic does not contain target {target}."
            )

        expected_count = int(
            summary_lookup.loc[target, "n_shared_qc_usable_orthologs"]
        )
        if len(written) != expected_count:
            die(
                f"Shared ortholog count mismatch for {target} genomes: "
                f"{shared_path} contains {len(written):,} IDs, but "
                "target_subset_summary reports "
                f"{expected_count:,}."
            )

        validated_sets[target] = set(written)

    # Genome subsets are nested 80 > 70 > 60 > 50, so shared ortholog sets must
    # expand monotonically or stay the same as genomes are removed.
    nesting_pairs = [(80, 70), (70, 60), (60, 50)]
    for larger_genome_set, smaller_genome_set in nesting_pairs:
        if not validated_sets[larger_genome_set].issubset(
            validated_sets[smaller_genome_set]
        ):
            lost = sorted(
                validated_sets[larger_genome_set]
                - validated_sets[smaller_genome_set]
            )[:20]
            die(
                "Shared-ortholog nesting validation failed: "
                f"{larger_genome_set}-genome orthologs are not a subset of "
                f"{smaller_genome_set}-genome orthologs. "
                f"Unexpectedly lost ortholog examples: {lost}"
            )


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
            "n_ortholog_records_absent",
            "accession",
        ],
        ascending=[False, True, True, True, True, True],
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
        ("n_ortholog_records_absent", "Absent", "absent"),
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
    ax.set_xlabel("Number of candidate orthologs", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel("Genome accession", fontsize=AXIS_LABEL_SIZE)
    ax.set_title(
        "Per-genome QC composition and occupancy (ordered by passing ortholog count)",
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

    protected_df = plot_df.loc[plot_df["protected"].astype(bool)]
    if not protected_df.empty:
        ax.scatter(
            protected_df["n_passing"],
            protected_df["n_rescued"],
            s=FIGURE2_PROTECTED_POINT_SIZE,
            marker=FIGURE2_PROTECTED_MARKER,
            facecolors="none",
            edgecolors="black",
            linewidths=FIGURE2_PROTECTED_EDGE_WIDTH,
            zorder=8,
            label="Protected genome",
        )

        if FIGURE2_ANNOTATE_PROTECTED:
            for row in protected_df.itertuples(index=False):
                ax.annotate(
                    row.accession,
                    xy=(row.n_passing, row.n_rescued),
                    xytext=FIGURE2_PROTECTED_ANNOTATION_OFFSET,
                    textcoords="offset points",
                    fontsize=FIGURE2_PROTECTED_ANNOTATION_FONT_SIZE,
                    fontweight="bold",
                    alpha=0.95,
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
    bars = ax.barh(
        y,
        plot_df["shared_ortholog_gain"],
        height=FIGURE3_BAR_HEIGHT,
        color=FIGURE3_DEFAULT_BAR_COLOR,
        edgecolor="none",
    )

    for bar, is_protected in zip(
        bars,
        plot_df["protected"].astype(bool).tolist(),
    ):
        if is_protected:
            bar.set_facecolor(FIGURE3_PROTECTED_BAR_COLOR)
            bar.set_edgecolor("black")
            bar.set_linewidth(FIGURE3_PROTECTED_EDGE_WIDTH)
            bar.set_hatch(FIGURE3_PROTECTED_HATCH)

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
        "Initial hypothetical one-genome removal benefit from the original "
        "starting cohort\n"
        "(protected genomes are shown for comparison but are never eligible "
        "for greedy removal)",
        fontsize=TITLE_SIZE,
    )

    if plot_df["protected"].astype(bool).any():
        ax.legend(
            handles=[
                Patch(
                    facecolor=FIGURE3_PROTECTED_BAR_COLOR,
                    edgecolor="black",
                    hatch=FIGURE3_PROTECTED_HATCH,
                    label="Protected genome (hypothetical benefit only)",
                )
            ],
            loc="best",
            fontsize=LEGEND_FONT_SIZE,
            frameon=False,
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

    genomes_dir = args.genomes_dir.expanduser().resolve()
    compleasm_dir = genomes_dir / "records" / "compleasm"

    input_path = (
        args.input.expanduser().resolve()
        if args.input is not None
        else (compleasm_dir / DEFAULT_INPUT_RELATIVE).resolve()
    )
    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else (compleasm_dir / DEFAULT_OUTPUT_RELATIVE).resolve()
    )
    alignments_dir = (
        args.alignments_dir.expanduser().resolve()
        if args.alignments_dir is not None
        else (compleasm_dir / DEFAULT_ALIGNMENTS_RELATIVE).resolve()
    )

    if not genomes_dir.exists():
        die(f"Genomes root directory does not exist: {genomes_dir}")
    if not genomes_dir.is_dir():
        die(f"Genomes root path is not a directory: {genomes_dir}")

    if input_path == output_dir:
        die("Input path resolves to the output directory, not a TSV file.")

    output_dir.mkdir(parents=True, exist_ok=True)
    alignments_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Genomes root: {genomes_dir}")
    print(f"[INFO] Compleasm root: {compleasm_dir}")
    print(f"[INFO] Input TSV: {input_path}")
    print(f"[INFO] Output directory: {output_dir}")
    print(f"[INFO] Alignments directory: {alignments_dir}")
    df = read_and_validate_input(
        input_path=input_path,
        expected_starting_genomes=args.expected_starting_genomes,
    )

    original_accessions = sorted(df["accession"].unique())
    starting_n = len(original_accessions)
    final_target = min(TARGET_GENOME_COUNTS)

    protected_accessions = validate_protected_accessions(
        requested=args.protected_accessions,
        starting_accessions=original_accessions,
        final_target=final_target,
    )

    print(f"[INFO] Protected accessions: {len(protected_accessions)}")
    for accession in sorted(protected_accessions):
        print(f"[INFO]   {accession}")

    print(
        f"[INFO] Input validated: {starting_n} genomes, "
        f"{df['odb12_id'].nunique():,} candidate ortholog IDs, "
        f"{len(df):,} observed accession/ortholog records."
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

    candidate_ortholog_count = df["odb12_id"].nunique()
    possible_pairs = starting_n * candidate_ortholog_count
    observed_pairs = len(df)
    absent_pairs = possible_pairs - observed_pairs
    overall_occupancy_pct = (
        100.0 * observed_pairs / possible_pairs if possible_pairs else 0.0
    )

    print(
        f"[INFO] Candidate ortholog IDs: {candidate_ortholog_count:,}"
    )
    print(
        f"[INFO] Observed accession/ortholog records: {observed_pairs:,}"
    )
    print(
        f"[INFO] Maximum possible accession/ortholog combinations: "
        f"{possible_pairs:,}"
    )
    print(
        f"[INFO] Absent accession/ortholog combinations: {absent_pairs:,}"
    )
    print(
        f"[INFO] Overall accession/ortholog occupancy: "
        f"{overall_occupancy_pct:.2f}%"
    )

    per_genome_summary = build_per_genome_summary(
        df,
        candidate_ortholog_count=candidate_ortholog_count,
        protected_accessions=protected_accessions,
    )
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
        protected_accessions=protected_accessions,
    )

    trajectory, subsets, removal_order = run_greedy_removal(
        usable_matrix=usable_matrix,
        rescued_matrix=rescued_matrix,
        per_genome_summary=per_genome_summary,
        final_target=final_target,
        protected_accessions=protected_accessions,
    )

    # Explicit validation of nested target subsets.
    validate_nested_subsets(
        original_accessions=original_accessions,
        subsets=subsets,
        removal_order=removal_order,
        target_counts=TARGET_GENOME_COUNTS,
        protected_accessions=protected_accessions,
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

    subset_paths, shared_by_target = write_subset_files(
        alignments_dir=alignments_dir,
        original_accessions=original_accessions,
        subsets=subsets,
        removal_order=removal_order,
        usable_matrix=usable_matrix,
    )

    # Post-write validation of subset file counts and contents.
    for target, directory in subset_paths.items():
        included_path = directory / "included_accessions.txt"
        dropped_path = directory / "dropped_accessions.txt"
        shared_path = directory / "shared_orthologs.txt"

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
        shared_orthologs = [
            x.strip()
            for x in shared_path.read_text(encoding="utf-8").splitlines()
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

        missing_protected = protected_accessions - set(included)
        if missing_protected:
            die(
                f"Post-write validation failed: protected accession(s) missing "
                f"from {included_path}: "
                + ", ".join(sorted(missing_protected))
            )

        dropped_protected = protected_accessions & set(dropped)
        if dropped_protected:
            die(
                f"Post-write validation failed: protected accession(s) found "
                f"in {dropped_path}: "
                + ", ".join(sorted(dropped_protected))
            )

        if len(shared_orthologs) != len(set(shared_orthologs)):
            die(
                f"Post-write validation failed: duplicate ortholog IDs in "
                f"{shared_path}."
            )

    # Validate shared-ortholog files against matrix-derived truth and summary
    # counts, including monotonic nesting across 80 -> 70 -> 60 -> 50.
    validate_shared_ortholog_files(
        usable_matrix=usable_matrix,
        subsets=subsets,
        subset_paths=subset_paths,
        shared_by_target=shared_by_target,
        target_subset_summary=target_subset_summary,
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
    print(f"Genomes root: {genomes_dir}")
    print(f"Compleasm root: {compleasm_dir}")
    print(f"Starting genomes: {starting_n}")
    print(f"Candidate ortholog IDs: {candidate_ortholog_count:,}")
    print(f"Observed accession/ortholog records: {observed_pairs:,}")
    print(f"Maximum possible accession/ortholog combinations: {possible_pairs:,}")
    print(f"Absent accession/ortholog combinations: {absent_pairs:,}")
    print(f"Overall accession/ortholog occupancy: {overall_occupancy_pct:.2f}%")
    print(f"Starting shared QC-usable orthologs: {starting_shared:,}")
    print(f"Protected genomes: {len(protected_accessions)}")
    if protected_accessions:
        print("Protected accessions:")
        for accession in sorted(protected_accessions):
            print(f"  {accession}")

    final_protected_ok = protected_accessions.issubset(set(subsets[final_target]))
    print(
        "All protected genomes retained in final "
        f"{final_target}-genome subset: {final_protected_ok}"
    )
    if not final_protected_ok:
        die(
            "Final protection validation failed: one or more protected genomes "
            "are missing from the final target subset."
        )

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

    print("\nSubset directories and authoritative selection files:")
    for target in TARGET_GENOME_COUNTS:
        directory = subset_paths[target]
        included_count = len(subsets[target])
        dropped_count = starting_n - target
        shared_count = len(shared_by_target[target])
        print(f"  {target} genomes: {directory}")
        print(f"    included accessions: {included_count}")
        print(f"    dropped accessions: {dropped_count}")
        print(f"    shared QC-usable orthologs: {shared_count:,}")
        print(f"    included_accessions.txt: {directory / 'included_accessions.txt'}")
        print(f"    shared_orthologs.txt: {directory / 'shared_orthologs.txt'}")

    print("\n[DONE] All validation checks passed.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
