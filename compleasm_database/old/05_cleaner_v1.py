#!/usr/bin/env python3
import pandas as pd
import argparse
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser


class AlignmentCleaner:
    def __init__(self, genomes_dir, manifest, threshold, errors, omit=None):
        self.genomes_dir = Path(genomes_dir).resolve()
        self.manifest = Path(manifest).resolve()
        self.threshold = threshold
        self.errors = errors
        self.omit = set(omit or [])

        self.macse_dir = self.genomes_dir / "records/compleasm/alignments/02_macse_alignments"
        self.fs_csv = self.genomes_dir / "records/compleasm/records/records_frameshift.csv"
        self.sum_csv = self.genomes_dir / "records/compleasm/records/records_filter_summary.csv"

        # Read project manifest instead of metadata.csv
        df_manifest = pd.read_csv(self.manifest)
        
        self.allowed_taxa = set(df_manifest["accession"].astype(str).str.strip())

        # Only omit taxa that are actually in the manifest
        self.valid_omit = self.allowed_taxa.intersection(self.omit)
        self.invalid_omit = self.omit.difference(self.allowed_taxa)

        if self.invalid_omit:
            print(
                "Warning: these --omit taxa were not found in the manifest and will be ignored:\n"
                + ", ".join(sorted(self.invalid_omit))
            )

        self.total_species = len(self.allowed_taxa) - len(self.valid_omit)

        self.outdir = self.genomes_dir / f"records/compleasm/alignments/03_clean_t{self.threshold}_e{self.errors}_o{len(self.valid_omit)}"
        self.outdir.mkdir(parents=True, exist_ok=True)



        # Assumption: the manifest has an accession column.
        # Change "accession" below if your manifest uses a different column name.
        if "accession" not in df_manifest.columns:
            raise ValueError(
                f"Manifest must contain an 'accession' column. "
                f"Found columns: {list(df_manifest.columns)}"
            )


        if self.total_species <= 0:
            raise ValueError(
                f"total_species became {self.total_species}. "
                f"Check the manifest and --omit values."
            )

    @staticmethod
    def normalize_header(header):
        return header.split()[0]

    def clean(self):
        print(f"--- Cleaning with Threshold: {self.threshold} ---")
        print(f"--- Cleaning with Errors: {self.errors} or less ---")
        print(f"--- Manifest taxa: {len(self.allowed_taxa)} ---")

        if self.valid_omit:
            print(f"--- Omitting {len(self.valid_omit)} taxa: {', '.join(sorted(self.valid_omit))} ---")

        print(f"--- Total species considered: {self.total_species} ---")

        fs_data, dropped_ids = [], []
        total_processed = 0

        for nt_file in self.macse_dir.glob("*_NT.fasta"):
            gene_id = nt_file.name.replace("_NT.fasta", "")
            total_processed += 1
            clean_seqs, clean_count = [], 0

            with open(nt_file, "r") as h:
                for header, seq in SimpleFastaParser(h):
                    seq_id = self.normalize_header(header)

                    # Skip taxa not in the manifest
                    if seq_id not in self.allowed_taxa:
                        continue

                    # Skip explicitly omitted taxa
                    if seq_id in self.valid_omit:
                        continue

                    fs_count = seq.count("!")
                    is_filtered = fs_count > int(self.errors)

                    fs_data.append({
                        "threshold": self.threshold,
                        "gene_id": gene_id,
                        "accession": seq_id,
                        "frameshifts": fs_count,
                        "sequence_filtered": is_filtered
                    })

                    if not is_filtered:
                        clean_count += 1
                        clean_seqs.append(f">{header}\n{seq.replace('!', '-')}\n")

            occupancy = clean_count / self.total_species

            if occupancy >= self.threshold:
                out_fasta = self.outdir / f"{gene_id}_clean_t{self.threshold}_e{self.errors}_o{len(self.valid_omit)}.fasta"
                with open(out_fasta, "w") as f:
                    f.writelines(clean_seqs)
            else:
                dropped_ids.append(gene_id)

        pd.DataFrame(fs_data).to_csv(
            self.fs_csv,
            mode="a",
            index=False,
            header=not self.fs_csv.exists()
        )

        summary = pd.DataFrame([{
            "threshold": self.threshold,
            "errors": self.errors,
            "manifest": str(self.manifest),
            "manifest_taxa_count": len(self.allowed_taxa),
            "total_species_considered": self.total_species,
            "omitted_taxa_count": len(self.valid_omit),
            "omitted_taxa": ",".join(sorted(self.valid_omit)),
            "total_genes_analyzed": total_processed,
            "genes_dropped": len(dropped_ids),
            "genes_retained": total_processed - len(dropped_ids),
            "dropped_gene_ids": ",".join(dropped_ids)
        }])

        summary.to_csv(
            self.sum_csv,
            mode="a",
            index=False,
            header=not self.sum_csv.exists()
        )

        print(f"Done. Retained {total_processed - len(dropped_ids)} genes.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("genomes_dir")
    parser.add_argument(
        "--manifest",
        required=True,
        help="Path to project manifest CSV containing the allowed taxa."
    )
    parser.add_argument("--threshold", type=float, default=0.9)
    parser.add_argument("--errors", type=int, default=1)
    parser.add_argument(
        "--omit",
        nargs="+",
        default=[],
        help="One or more taxa/accessions from the manifest to omit."
    )
    args = parser.parse_args()

    AlignmentCleaner(
        args.genomes_dir,
        args.manifest,
        args.threshold,
        args.errors,
        args.omit
    ).clean()
