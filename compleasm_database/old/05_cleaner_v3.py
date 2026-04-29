#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser


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

        if macse_dir:
            self.macse_dir = Path(macse_dir).resolve()
        else:
            self.macse_dir = self.genomes_dir / "records/compleasm/alignments/02_macse_alignments"

        self.fs_csv = self.genomes_dir / "records/compleasm/records/records_frameshift.csv"
        self.sum_csv = self.genomes_dir / "records/compleasm/records/records_filter_summary.csv"

        df_manifest = pd.read_csv(self.manifest)

        if "accession" not in df_manifest.columns:
            raise ValueError(
                f"Manifest must contain an 'accession' column. "
                f"Found columns: {list(df_manifest.columns)}"
            )

        self.allowed_taxa = set(df_manifest["accession"].astype(str).str.strip())

        self.valid_omit = self.allowed_taxa.intersection(self.omit)
        self.invalid_omit = self.omit.difference(self.allowed_taxa)

        if self.invalid_omit:
            print(
                "Warning: these --omit taxa were not found in the manifest and will be ignored:\n"
                + ", ".join(sorted(self.invalid_omit))
            )

        self.total_species = len(self.allowed_taxa) - len(self.valid_omit)

        self.outdir = self.genomes_dir / (
            f"records/compleasm/alignments/03_clean_t{self.threshold}_e{self.errors}_o{len(self.valid_omit)}"
        )
        self.outdir.mkdir(parents=True, exist_ok=True)

        if self.total_species <= 0:
            raise ValueError(
                f"total_species became {self.total_species}. "
                f"Check the manifest and --omit values."
            )

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
        dropped_files = 0

        for nt_file in nt_files:
            gene_id = nt_file.name.replace("_NT.fasta", "")
            if gene_id in self.drop_genes:
                dropped_files += 1
                continue
            yield nt_file

        if self.drop_genes:
            print(f"--- Dropping {len(self.drop_genes)} requested genes from all processing ---")
            matched = len({p.name.replace('_NT.fasta', '') for p in nt_files}.intersection(self.drop_genes))
            unmatched = sorted(self.drop_genes.difference({p.name.replace('_NT.fasta', '') for p in nt_files}))
            print(f"--- Matched {matched} genes in the alignment directory; skipped {dropped_files} files ---")
            if unmatched:
                print(
                    "Warning: these --drop_genes IDs were not found in the alignment directory and will be ignored:\n"
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
        fs_lookup = {}
        for row in fs_data:
            fs_lookup.setdefault(row["gene_id"], []).append(row)

        for nt_file in self.iter_nt_files():
            gene_id = nt_file.name.replace("_NT.fasta", "")
            clean_seqs, clean_count = [], 0

            gene_rows = {
                row["accession"]: row["sequence_filtered"]
                for row in fs_lookup.get(gene_id, [])
            }

            with open(nt_file, "r") as handle:
                for header, seq in SimpleFastaParser(handle):
                    seq_id = self.normalize_header(header)

                    if seq_id not in self.allowed_taxa:
                        continue

                    if seq_id in self.valid_omit:
                        continue

                    is_filtered = gene_rows.get(seq_id)
                    if is_filtered is None:
                        continue

                    if not is_filtered:
                        clean_count += 1
                        clean_seqs.append(f">{header}\n{seq.replace('!', '-')}\n")

            occupancy = clean_count / self.total_species

            if occupancy >= self.threshold:
                out_fasta = self.outdir / (
                    f"{gene_id}_clean_t{self.threshold}_e{self.errors}_o{len(self.valid_omit)}.fasta"
                )
                with open(out_fasta, "w") as out_handle:
                    out_handle.writelines(clean_seqs)
            else:
                dropped_ids.append(gene_id)

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("genomes_dir")
    parser.add_argument(
        "--macse_dir",
        required=False,
        help="Path to directory containing macse alignments.")
    parser.add_argument(
        "--manifest",
        required=True,
        help="Path to project manifest CSV containing the allowed taxa.",
    )
    parser.add_argument("--threshold", type=float, default=0.9)
    parser.add_argument("--errors", type=int, default=1)
    parser.add_argument(
        "--omit",
        nargs="+",
        default=[],
        help="One or more taxa/accessions from the manifest to omit.",
    )
    parser.add_argument(
        "--drop_genes",
        nargs="+",
        default=[],
        help="One or more gene IDs to exclude from all processing and output alignments.",
    )
    parser.add_argument(
        "--frameshift",
        action="store_true",
        help="Only write the frameshift CSV, without producing cleaned FASTA files or filter summary output.",
    )
    args = parser.parse_args()

    AlignmentCleaner(
        genomes_dir=args.genomes_dir,
        manifest=args.manifest,
        threshold=args.threshold,
        errors=args.errors,
        omit=args.omit,
        macse_dir=args.macse_dir,
        frameshift_only=args.frameshift,
        drop_genes=args.drop_genes,
    ).clean()
