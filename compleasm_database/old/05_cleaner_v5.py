#!/usr/bin/env python3
"""Clean MACSE alignments using unique manifest-accession occupancy.

With --threshold 1, every non-omitted manifest accession must pass the
frameshift filter. Duplicate expected accessions cause the gene to be dropped
because choosing one record would be ambiguous. This version retains v4's
frameshift filtering rules.
"""
import argparse
from pathlib import Path

import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser


def load_id_file(path):
    """
    Load one ID per line from a plain-text file.
    Blank lines and lines starting with '#' are ignored.
    """
    file_path = Path(path).resolve()
    if not file_path.exists():
        raise FileNotFoundError(f"List file not found: {file_path}")

    items = []
    with open(file_path, "r") as handle:
        for line in handle:
            value = line.strip()
            if not value or value.startswith("#"):
                continue
            items.append(value)
    return items


class AlignmentCleaner:
    def __init__(
        self,
        genomes_dir,
        manifest,
        threshold,
        errors,
        omit=None,
        frameshift_only=False,
        macse_dir=None,
        drop_genes=None,
    ):
        self.genomes_dir = Path(genomes_dir).resolve()
        self.manifest = Path(manifest).resolve()
        self.threshold = threshold
        self.errors = errors
        self.omit = set(omit or [])
        self.frameshift_only = frameshift_only
        self.drop_genes = set(drop_genes or [])

        if not 0 <= self.threshold <= 1:
            raise ValueError("threshold must be between 0 and 1, inclusive.")
        if self.errors < 0:
            raise ValueError("errors must be nonnegative.")

        if macse_dir:
            self.macse_dir = Path(macse_dir).resolve()
        else:
            self.macse_dir = self.genomes_dir / "records/compleasm/alignments/02_macse_alignments"

        self.fs_csv = self.genomes_dir / "records/compleasm/records/records_frameshift.csv"
        self.sum_csv = self.genomes_dir / "records/compleasm/records/records_filter_summary.csv"

        df_manifest = pd.read_csv(self.manifest, dtype=str, keep_default_na=False)

        if "accession" not in df_manifest.columns:
            raise ValueError(
                f"Manifest must contain an 'accession' column. "
                f"Found columns: {list(df_manifest.columns)}"
            )

        self.allowed_taxa = set(df_manifest["accession"].astype(str).str.strip())
        if "" in self.allowed_taxa:
            raise ValueError("Manifest contains a blank accession.")

        self.valid_omit = self.allowed_taxa.intersection(self.omit)
        self.invalid_omit = self.omit.difference(self.allowed_taxa)

        if self.invalid_omit:
            print(
                "Warning: these omitted taxa were not found in the manifest and will be ignored:\n"
                + ", ".join(sorted(self.invalid_omit))
            )

        self.expected_taxa = self.allowed_taxa.difference(self.valid_omit)
        self.total_species = len(self.expected_taxa)

        self.outdir = self.genomes_dir / (
            f"records/compleasm/alignments/03_clean_t{self.threshold}_e{self.errors}_o{len(self.valid_omit)}"
        )
        self.outdir.mkdir(parents=True, exist_ok=True)

        if self.total_species <= 0:
            raise ValueError(
                f"total_species became {self.total_species}. "
                f"Check the manifest and omitted accessions."
            )

        self.occupancy_csv = self.outdir / "records_unique_occupancy.csv"

    @staticmethod
    def normalize_header(header):
        return header.split()[0]

    @staticmethod
    def prompt_csv_mode(csv_path: Path) -> str:
        """
        Ask whether to append to or rewrite an existing frameshift CSV.
        Returns either 'a' or 'w'.
        """
        csv_path.parent.mkdir(parents=True, exist_ok=True)

        if not csv_path.exists():
            print(f"Frameshift CSV does not exist yet. A new file will be created at: {csv_path}")
            return "w"

        while True:
            response = input(
                f"Frameshift CSV already exists at:\n{csv_path}\n"
                "Create a new frameshift CSV by rewriting the old one, or append to the existing file?\n"
                "Enter 'rewrite' or 'append': "
            ).strip().lower()

            if response in {"rewrite", "r"}:
                return "w"
            if response in {"append", "a"}:
                return "a"

            print("Please enter 'rewrite' or 'append'.")

    def iter_nt_files(self):
        nt_files = sorted(self.macse_dir.glob("*_NT.fasta"))
        available_gene_ids = {p.name.replace("_NT.fasta", "") for p in nt_files}
        dropped_files = 0

        for nt_file in nt_files:
            gene_id = nt_file.name.replace("_NT.fasta", "")
            if gene_id in self.drop_genes:
                dropped_files += 1
                continue
            yield nt_file

        if self.drop_genes:
            print(f"--- Dropping {len(self.drop_genes)} requested genes from all processing ---")
            matched = len(available_gene_ids.intersection(self.drop_genes))
            unmatched = sorted(self.drop_genes.difference(available_gene_ids))
            print(f"--- Matched {matched} genes in the alignment directory; skipped {dropped_files} files ---")
            if unmatched:
                print(
                    "Warning: these dropped gene IDs were not found in the alignment directory and will be ignored:\n"
                    + ", ".join(unmatched)
                )

    def collect_frameshift_data(self):
        fs_data = []
        total_processed = 0

        print(f"--- Scanning alignments for frameshifts ---")
        print(f"--- Manifest taxa: {len(self.allowed_taxa)} ---")

        if self.valid_omit:
            print(f"--- Omitting {len(self.valid_omit)} taxa: {', '.join(sorted(self.valid_omit))} ---")

        print(f"--- Total species considered: {self.total_species} ---")

        for nt_file in self.iter_nt_files():
            gene_id = nt_file.name.replace("_NT.fasta", "")
            total_processed += 1

            with open(nt_file, "r") as handle:
                for header, seq in SimpleFastaParser(handle):
                    seq_id = self.normalize_header(header)

                    if seq_id not in self.allowed_taxa:
                        continue

                    if seq_id in self.valid_omit:
                        continue

                    fs_count = seq.count("!")
                    is_filtered = fs_count > int(self.errors)

                    fs_data.append({
                        "threshold": self.threshold,
                        "gene_id": gene_id,
                        "accession": seq_id,
                        "frameshifts": fs_count,
                        "sequence_filtered": is_filtered,
                    })

        return fs_data, total_processed

    def write_frameshift_csv(self, fs_data):
        mode = self.prompt_csv_mode(self.fs_csv)
        write_header = (mode == "w") or (not self.fs_csv.exists())

        pd.DataFrame(fs_data).to_csv(
            self.fs_csv,
            mode=mode,
            index=False,
            header=write_header,
        )

        action = "rewritten" if mode == "w" else "appended to"
        print(f"Frameshift CSV {action}: {self.fs_csv}")

    def clean(self):
        print(f"--- Cleaning with Threshold: {self.threshold} ---")
        print(f"--- Cleaning with Errors: {self.errors} or less ---")
        if self.frameshift_only:
            print("--- Running in --frameshift mode: only the frameshift CSV will be written ---")

        fs_data, total_processed = self.collect_frameshift_data()
        self.write_frameshift_csv(fs_data)

        if self.frameshift_only:
            print(f"Done. Frameshift-only mode processed {total_processed} genes.")
            return

        dropped_ids = []
        occupancy_rows = []

        # Remove outputs for explicitly excluded genes when reusing an output directory.
        for gene_id in self.drop_genes:
            self.output_fasta_path(gene_id).unlink(missing_ok=True)

        for nt_file in self.iter_nt_files():
            gene_id = nt_file.name.replace("_NT.fasta", "")
            clean_seqs = []
            seen_taxa = set()
            duplicate_taxa = set()
            passing_taxa = set()

            with open(nt_file, "r") as handle:
                for header, seq in SimpleFastaParser(handle):
                    seq_id = self.normalize_header(header)

                    if seq_id not in self.expected_taxa:
                        continue

                    if seq_id in seen_taxa:
                        duplicate_taxa.add(seq_id)
                    seen_taxa.add(seq_id)

                    # Evaluate each record directly; duplicate IDs must not overwrite QC.
                    is_filtered = seq.count("!") > self.errors
                    if not is_filtered:
                        passing_taxa.add(seq_id)
                        clean_seqs.append(f">{header}\n{seq.replace('!', '-')}\n")

            occupancy = len(passing_taxa) / self.total_species
            missing_taxa = self.expected_taxa.difference(passing_taxa)
            retained = not duplicate_taxa and occupancy >= self.threshold
            # An explicit set check guarantees complete coverage at threshold 1.
            if self.threshold == 1:
                retained = retained and passing_taxa == self.expected_taxa

            reasons = []
            if duplicate_taxa:
                reasons.append("duplicate_accessions")
                print(
                    f"Warning: dropping {gene_id}; duplicate expected accessions: "
                    + ", ".join(sorted(duplicate_taxa))
                )
            if occupancy < self.threshold:
                reasons.append("below_occupancy_threshold")
            occupancy_rows.append({
                "gene_id": gene_id,
                "expected_accession_count": self.total_species,
                "passing_unique_accession_count": len(passing_taxa),
                "occupancy": occupancy,
                "missing_or_filtered_accessions": ",".join(sorted(missing_taxa)),
                "duplicate_accessions": ",".join(sorted(duplicate_taxa)),
                "retained": retained,
                "drop_reason": ";".join(reasons),
            })
            out_fasta = self.output_fasta_path(gene_id)

            if retained:
                with open(out_fasta, "w") as out_handle:
                    out_handle.writelines(clean_seqs)
            else:
                dropped_ids.append(gene_id)
                # A previous run must not leave a now-rejected alignment in the output.
                out_fasta.unlink(missing_ok=True)

        pd.DataFrame(occupancy_rows, columns=[
            "gene_id", "expected_accession_count", "passing_unique_accession_count",
            "occupancy", "missing_or_filtered_accessions", "duplicate_accessions",
            "retained", "drop_reason",
        ]).to_csv(self.occupancy_csv, index=False)
        print(f"Unique-accession occupancy audit: {self.occupancy_csv}")

        summary = pd.DataFrame([
            {
                "threshold": self.threshold,
                "errors": self.errors,
                "manifest": str(self.manifest),
                "manifest_taxa_count": len(self.allowed_taxa),
                "total_species_considered": self.total_species,
                "omitted_taxa_count": len(self.valid_omit),
                "omitted_taxa": ",".join(sorted(self.valid_omit)),
                "drop_genes_count": len(self.drop_genes),
                "drop_genes": ",".join(sorted(self.drop_genes)),
                "total_genes_analyzed": total_processed,
                "genes_dropped": len(dropped_ids),
                "genes_retained": total_processed - len(dropped_ids),
                "dropped_gene_ids": ",".join(dropped_ids),
            }
        ])

        summary.to_csv(
            self.sum_csv,
            mode="a",
            index=False,
            header=not self.sum_csv.exists(),
        )

        print(f"Done. Retained {total_processed - len(dropped_ids)} genes.")

    def output_fasta_path(self, gene_id):
        return self.outdir / (
            f"{gene_id}_clean_t{self.threshold}_e{self.errors}_o{len(self.valid_omit)}.fasta"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("genomes_dir")
    parser.add_argument(
        "--macse_dir",
        required=False,
        help="Path to directory containing macse alignments.",
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="Path to project manifest CSV containing the allowed taxa.",
    )
    parser.add_argument(
        "--threshold", type=float, default=0.9,
        help="Minimum fraction (0 to 1) of unique expected accessions passing QC. "
             "At 1, every expected accession is required; duplicate accessions drop the gene.",
    )
    parser.add_argument("--errors", type=int, default=1)
    parser.add_argument(
        "--omit",
        nargs="+",
        default=[],
        help="One or more taxa/accessions from the manifest to omit.",
    )
    parser.add_argument(
        "--omit_file",
        required=False,
        help="Path to a text file containing one accession per line to omit.",
    )
    parser.add_argument(
        "--drop_genes",
        nargs="+",
        default=[],
        help="One or more gene IDs to exclude from all processing and output alignments.",
    )
    parser.add_argument(
        "--drop_genes_file",
        required=False,
        help="Path to a text file containing one gene ID per line to exclude.",
    )
    parser.add_argument(
        "--frameshift",
        action="store_true",
        help="Only write the frameshift CSV, without producing cleaned FASTA files or filter summary output.",
    )
    args = parser.parse_args()

    omit_values = list(args.omit)
    if args.omit_file:
        omit_values.extend(load_id_file(args.omit_file))

    drop_gene_values = list(args.drop_genes)
    if args.drop_genes_file:
        drop_gene_values.extend(load_id_file(args.drop_genes_file))

    AlignmentCleaner(
        genomes_dir=args.genomes_dir,
        manifest=args.manifest,
        threshold=args.threshold,
        errors=args.errors,
        omit=sorted(set(omit_values)),
        macse_dir=args.macse_dir,
        frameshift_only=args.frameshift,
        drop_genes=sorted(set(drop_gene_values)),
    ).clean()
