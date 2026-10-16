#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Find shared single-copy Compleasm orthologs for a manifest-defined cohort.

The accession matching is repository-neutral. NCBI accessions (GCA_/GCF_),
CNCB Genome Warehouse accessions (GWH...), and other accession formats are
matched first by their complete accession and then, when needed, by the
unversioned accession root.

Assembly-level filters are normalized so labels such as ``chromosome``,
``Chromosome``, ``Complete Chromosome``, and
``Draft genome in chromosome level`` are treated equivalently.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter
from pathlib import Path
from typing import Iterable

import pandas as pd


# ---------------------------------------------------------------------------
# Example usage
# ---------------------------------------------------------------------------
#
# NCBI and GWH accessions may be mixed in the same manifest:
#
#   python 03_single_orthologs_csv_gwh.py \
#       ../../../genomes \
#       /path/to/project_manifest.csv \
#       chromosome scaffold
#
# Retain every assembly level:
#
#   python 03_single_orthologs_csv_gwh.py \
#       ../../../genomes \
#       /path/to/project_manifest.csv
#
# Require 100% occupancy:
#
#   python 03_single_orthologs_csv_gwh.py \
#       ../../../genomes \
#       /path/to/project_manifest.csv \
#       chromosome \
#       --threshold 1.0
#
# The manifest must contain an ``accession`` column. For example:
#
#   accession
#   GCA_012345678.1
#   GCF_012345678.1
#   GWHFILM00000000.1
#


def normalize_accession(accession: object) -> str:
    """Return a stable accession string without assuming a repository.

    Known NCBI and GWH prefixes are normalized to uppercase. Other identifiers
    retain their original case so future repository-specific identifiers are
    not accidentally changed.
    """

    text = str(accession).strip()
    if not text or text.lower() == "nan":
        return ""

    upper = text.upper()
    if upper.startswith(("GCA_", "GCF_", "GWH")):
        return upper
    return text


def root_acc(accession: object) -> str:
    """Remove the final dot-delimited numeric version from an accession.

    Examples
    --------
    ``GCA_012345678.1`` -> ``GCA_012345678``
    ``GWHFILM00000000.1`` -> ``GWHFILM00000000``

    If the suffix after the final dot is not numeric, the identifier is left
    unchanged. This avoids truncating identifiers whose dots are not versions.
    """

    normalized = normalize_accession(accession)
    if not normalized or "." not in normalized:
        return normalized

    root, suffix = normalized.rsplit(".", 1)
    return root if suffix.isdigit() else normalized


def accession_version(accession: object) -> int:
    """Return the final numeric accession version, or -1 when absent."""

    normalized = normalize_accession(accession)
    if "." not in normalized:
        return -1
    suffix = normalized.rsplit(".", 1)[-1]
    return int(suffix) if suffix.isdigit() else -1


def accession_repository(accession: object) -> str:
    """Return a human-readable repository label for diagnostics only."""

    normalized = normalize_accession(accession)
    if normalized.startswith(("GCA_", "GCF_")):
        return "NCBI"
    if normalized.startswith("GWH"):
        return "CNCB-GWH"
    return "other/unknown"


def normalize_assembly_level(value: object) -> str:
    """Normalize assembly-level labels from NCBI, GWH, and legacy tables."""

    text = str(value).strip().lower()
    if not text or text == "nan":
        return ""

    # Remove punctuation and collapse whitespace before matching aliases.
    simplified = re.sub(r"[^a-z0-9]+", " ", text).strip()

    # Check chromosome before complete because legacy labels may be
    # ``Complete Chromosome``. In that case the useful filtering category is
    # chromosome, not complete.
    if "chromosome" in simplified:
        return "chromosome"
    if "complete" in simplified:
        return "complete"
    if "scaffold" in simplified:
        return "scaffold"
    if "contig" in simplified:
        return "contig"

    return simplified.replace(" ", "_")


def normalize_target_levels(values: Iterable[str]) -> set[str]:
    """Normalize CLI assembly-level filters and discard blank values."""

    return {
        normalized
        for normalized in (normalize_assembly_level(value) for value in values)
        if normalized
    }


class SingleOrthologFinder:
    """Find BUSCO/Compleasm single-copy genes meeting a cohort occupancy cutoff."""

    def __init__(
        self,
        genomes_dir: str,
        accessions_file: str,
        targets: list[str],
        threshold_pct: float = 0.80,
    ) -> None:
        self.genomes_dir = Path(genomes_dir).expanduser().resolve()
        self.accessions_file = Path(accessions_file).expanduser().resolve()
        self.targets_original = targets or []
        self.targets = normalize_target_levels(self.targets_original)
        self.threshold_pct = threshold_pct

        if not 0.0 < self.threshold_pct <= 1.0:
            raise ValueError(
                f"--threshold must be greater than 0 and at most 1; "
                f"received {self.threshold_pct}"
            )

        self.gene_id_dir = (
            self.genomes_dir
            / "records/compleasm/mb_downloads/sauropsida_odb12/hmms"
        )
        self.genomes_meta = self.genomes_dir / "records/genomes_metadata.csv"
        self.compleasm_meta = (
            self.genomes_dir / "records/compleasm/records/metadata.csv"
        )
        self.records_dir = self.genomes_dir / "records/compleasm/records"

        if self.gene_id_dir.exists():
            self.gene_id_list = sorted(
                file.stem for file in self.gene_id_dir.glob("*.hmm")
            )
        else:
            print(f"Error: HMM directory not found at {self.gene_id_dir}")
            self.gene_id_list = []

        self.allowed_accessions, self.allowed_roots = self._load_allowed_accessions()
        self.cohort_label = self._sanitize_label(self.accessions_file.stem)
        self.output_path = (
            self.records_dir / f"shared_single_genes_{self.cohort_label}.csv"
        )
        self.single_gene_ids: list[str] = []

    @staticmethod
    def _sanitize_label(text: str) -> str:
        return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_")

    @staticmethod
    def _validate_table(path: Path, required_columns: set[str], label: str) -> None:
        if not path.is_file():
            raise FileNotFoundError(f"{label} not found: {path}")
        columns = set(pd.read_csv(path, nrows=0).columns)
        missing = required_columns - columns
        if missing:
            raise ValueError(
                f"{label} is missing required column(s) "
                f"{', '.join(sorted(missing))}: {path}"
            )

    def _read_accessions_txt(self, txt_path: Path) -> list[str]:
        accessions: list[str] = []
        with txt_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                value = line.strip()
                if not value or value.startswith("#"):
                    continue
                accession = normalize_accession(value.split()[0])
                if accession:
                    accessions.append(accession)
        return sorted(set(accessions))

    def _read_manifest_accessions(self, manifest_path: Path) -> list[str]:
        dataframe = pd.read_csv(manifest_path)
        if "accession" not in dataframe.columns:
            raise ValueError(
                f"Manifest is missing required 'accession' column: {manifest_path}"
            )

        accessions = [
            normalize_accession(value)
            for value in dataframe["accession"].dropna().tolist()
        ]
        return sorted({accession for accession in accessions if accession})

    def _load_allowed_accessions(self) -> tuple[set[str], set[str]]:
        if not self.accessions_file.is_file():
            raise FileNotFoundError(
                f"Manifest/accession file not found: {self.accessions_file}"
            )

        suffix = self.accessions_file.suffix.lower()
        if suffix == ".csv":
            accessions = self._read_manifest_accessions(self.accessions_file)
            print(
                f"Loaded {len(accessions)} accessions from manifest: "
                f"{self.accessions_file}"
            )
        else:
            accessions = self._read_accessions_txt(self.accessions_file)
            print(
                f"Loaded {len(accessions)} accessions from txt file: "
                f"{self.accessions_file}"
            )

        if not accessions:
            raise ValueError(
                f"No usable accessions were found in {self.accessions_file}"
            )

        repository_counts = Counter(accession_repository(acc) for acc in accessions)
        repository_summary = ", ".join(
            f"{repository}={count}"
            for repository, count in sorted(repository_counts.items())
        )
        print(f"Accession repositories detected: {repository_summary}")

        roots = {root_acc(accession) for accession in accessions}
        return set(accessions), roots

    def _repath(self, path_to_fix: object) -> Path:
        """Rebase saved paths onto the active genomes directory when possible."""

        path_str = str(path_to_fix).strip()
        if not path_str or path_str.lower() == "nan":
            return Path("")

        expanded = Path(path_str).expanduser()
        if expanded.exists():
            return expanded.resolve()

        # Preserve the original behavior for metadata created on another
        # computer: keep everything after the nearest ``genomes`` component and
        # attach it to the active genomes directory.
        parts = expanded.parts
        matching_indexes = [
            index for index, part in enumerate(parts) if part == self.genomes_dir.name
        ]
        if matching_indexes:
            suffix_parts = parts[matching_indexes[-1] + 1 :]
            return self.genomes_dir.joinpath(*suffix_parts)

        return expanded

    def _prepare_accession_columns(
        self,
        dataframe: pd.DataFrame,
        accession_col: str = "accession",
    ) -> pd.DataFrame:
        prepared = dataframe.copy()
        if accession_col not in prepared.columns:
            raise ValueError(
                f"Metadata table is missing required '{accession_col}' column"
            )
        prepared[accession_col] = prepared[accession_col].map(normalize_accession)
        prepared = prepared[prepared[accession_col] != ""].copy()
        prepared["accession_root"] = prepared[accession_col].map(root_acc)
        prepared["accession_version"] = prepared[accession_col].map(
            accession_version
        )
        return prepared

    def _resolve_cohort_rows(
        self,
        dataframe: pd.DataFrame,
        accession_col: str = "accession",
    ) -> tuple[pd.DataFrame, list[str]]:
        """Resolve cohort rows by exact accession, then unversioned root.

        This procedure is generic for GCA_, GCF_, GWH, and future accession
        prefixes. Exact versioned matches are always preferred. A root fallback
        is used only when the requested exact accession is absent.
        """

        prepared = self._prepare_accession_columns(dataframe, accession_col)

        exact = prepared[
            prepared[accession_col].isin(self.allowed_accessions)
        ].copy()
        matched_roots = set(exact["accession_root"])
        unresolved_roots = self.allowed_roots - matched_roots

        fallback = prepared[
            (~prepared[accession_col].isin(self.allowed_accessions))
            & (prepared["accession_root"].isin(unresolved_roots))
        ].copy()

        if not fallback.empty:
            fallback = (
                fallback.sort_values(
                    ["accession_root", "accession_version", accession_col],
                    ascending=[True, False, True],
                )
                .drop_duplicates(subset=["accession_root"], keep="first")
            )

        resolved = pd.concat([exact, fallback], ignore_index=True)
        if resolved.empty:
            return resolved, sorted(self.allowed_roots)

        # If both exact and root-level metadata rows are present, preserve the
        # exact requested version. The boolean sort key explicitly ranks exact
        # matches ahead of fallback rows.
        resolved["_is_exact_requested"] = resolved[accession_col].isin(
            self.allowed_accessions
        )
        resolved = (
            resolved.sort_values(
                [
                    "accession_root",
                    "_is_exact_requested",
                    "accession_version",
                    accession_col,
                ],
                ascending=[True, False, False, True],
            )
            .drop_duplicates(subset=["accession_root"], keep="first")
            .drop(columns=["_is_exact_requested"])
            .reset_index(drop=True)
        )

        missing_roots = sorted(
            self.allowed_roots - set(resolved["accession_root"])
        )
        return resolved, missing_roots

    def mask(self) -> list[str]:
        if not self.gene_id_list:
            print("No BUSCO IDs found in HMM directory, aborting.")
            return []

        self._validate_table(
            self.genomes_meta,
            {"accession", "assembly_level"},
            "Genome metadata",
        )
        self._validate_table(
            self.compleasm_meta,
            {"accession", "full_table"},
            "Compleasm metadata",
        )

        genomes_all = pd.read_csv(self.genomes_meta)
        compleasm_all = pd.read_csv(self.compleasm_meta)

        genomes_df, missing_genome_roots = self._resolve_cohort_rows(
            genomes_all,
            accession_col="accession",
        )
        compleasm_df, missing_compleasm_roots = self._resolve_cohort_rows(
            compleasm_all,
            accession_col="accession",
        )

        if not genomes_df.empty:
            genomes_df["assembly_level_normalized"] = genomes_df[
                "assembly_level"
            ].map(normalize_assembly_level)

        level_by_full = dict(
            zip(
                genomes_df.get("accession", pd.Series(dtype=str)),
                genomes_df.get(
                    "assembly_level_normalized", pd.Series(dtype=str)
                ),
            )
        )
        level_by_root = dict(
            zip(
                genomes_df.get("accession_root", pd.Series(dtype=str)),
                genomes_df.get(
                    "assembly_level_normalized", pd.Series(dtype=str)
                ),
            )
        )

        gene_occurrence_counter: Counter[str] = Counter()
        included_count = 0
        skipped_for_target = 0
        missing_level = 0
        missing_tsv = 0
        processing_errors = 0

        if compleasm_df.empty:
            print("No Compleasm metadata rows matched the manifest/txt cohort.")
            print(f"Accessions requested: {len(self.allowed_accessions)}")
            print(
                "Missing accession roots in Compleasm metadata: "
                f"{len(missing_compleasm_roots)}"
            )
            if missing_compleasm_roots:
                print("  " + "\n  ".join(missing_compleasm_roots))
            return []

        for _, row in compleasm_df.iterrows():
            accession = row["accession"]
            accession_root = row["accession_root"]
            level = level_by_full.get(
                accession,
                level_by_root.get(accession_root, ""),
            )

            if not level:
                missing_level += 1
                print(
                    f"Warning: no assembly_level metadata was found for "
                    f"{accession}; this genome will "
                    + (
                        "be excluded because assembly-level filters were supplied."
                        if self.targets
                        else "remain eligible because no assembly-level filter was supplied."
                    )
                )

            if self.targets and level not in self.targets:
                skipped_for_target += 1
                print(
                    f"Skipped {accession}: normalized assembly level "
                    f"{level or '<missing>'!r} not in {sorted(self.targets)}"
                )
                continue

            tsv_path = self._repath(row["full_table"])
            if not str(tsv_path) or not tsv_path.is_file():
                print(f"Warning: Compleasm full_table missing for {accession}: {tsv_path}")
                missing_tsv += 1
                continue

            try:
                tsv_dataframe = pd.read_csv(
                    tsv_path,
                    sep="\t",
                    comment="#",
                    header=0,
                    usecols=[0, 1],
                    names=["Gene", "Status"],
                )
                current_single = tsv_dataframe.loc[
                    tsv_dataframe["Status"] == "Single",
                    "Gene",
                ].astype(str).tolist()
                gene_occurrence_counter.update(current_single)
                print(
                    f"Processed {accession} "
                    f"({accession_repository(accession)}): "
                    f"assembly_level={level or '<missing>'}"
                )
                included_count += 1
            except Exception as exc:
                processing_errors += 1
                print(f"Error processing {accession}: {type(exc).__name__}: {exc}")

        if included_count == 0:
            print("No genomes passed the filters. Nothing to write.")
            print(f"Accessions requested: {len(self.allowed_accessions)}")
            print(f"Compleasm cohort rows resolved: {len(compleasm_df)}")
            print(f"Rows skipped because of assembly level: {skipped_for_target}")
            print(f"Rows with missing assembly level: {missing_level}")
            return []

        required_count = max(1, math.ceil(included_count * self.threshold_pct))
        gene_id_set = set(self.gene_id_list)
        self.single_gene_ids = sorted(
            gene
            for gene, count in gene_occurrence_counter.items()
            if count >= required_count and gene in gene_id_set
        )

        self.records_dir.mkdir(parents=True, exist_ok=True)
        with self.output_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(self.single_gene_ids)

        print("\nFinished.")
        print(f"Cohort filter file: {self.accessions_file}")
        print(f"Accessions requested: {len(self.allowed_accessions)}")
        print(f"Genome metadata rows resolved: {len(genomes_df)}")
        print(f"Compleasm cohort rows resolved: {len(compleasm_df)}")
        print(
            "Assembly levels retained: "
            f"{sorted(self.targets) if self.targets else 'all'}"
        )
        print(f"Rows skipped because of assembly level: {skipped_for_target}")
        print(f"Rows with missing assembly level: {missing_level}")
        print(f"Rows skipped because full_table was missing: {missing_tsv}")
        print(f"Rows with processing errors: {processing_errors}")
        print(
            "Missing accession roots in genomes metadata: "
            f"{len(missing_genome_roots)}"
        )
        print(
            "Missing accession roots in Compleasm metadata: "
            f"{len(missing_compleasm_roots)}"
        )
        print(
            f"Occupancy Threshold: {self.threshold_pct * 100:.1f}% "
            f"({required_count}/{included_count} genomes)"
        )
        print(f"Genes Meeting Threshold: {len(self.single_gene_ids)}")
        print(f"Shared IDs written to: {self.output_path}")
        return self.single_gene_ids


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Find BUSCO/Compleasm single-copy orthologs across a manifest-defined "
            "cohort containing NCBI, CNCB-GWH, or mixed repository accessions."
        )
    )
    parser.add_argument("genomes_dir", help="Path to genomes directory")
    parser.add_argument(
        "accessions",
        help=(
            "CSV manifest with an accession column, or TXT file with one "
            "accession per line"
        ),
    )
    parser.add_argument(
        "targets",
        nargs="*",
        help=(
            "Optional assembly levels to retain. Values are normalized, so "
            "'Complete Chromosome', 'chromosome', 'Scaffold', and GWH-style "
            "descriptions are supported. Omit to retain all levels."
        ),
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.80,
        help="Occupancy threshold greater than 0 and at most 1 (default: 0.80)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    finder = SingleOrthologFinder(
        args.genomes_dir,
        args.accessions,
        args.targets,
        threshold_pct=args.threshold,
    )
    finder.mask()


if __name__ == "__main__":
    main()
