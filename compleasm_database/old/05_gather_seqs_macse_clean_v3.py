#!/usr/bin/env python3
import os
import subprocess
import pandas as pd
import argparse
import random
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser

class MacsePipeline:
    def __init__(self, genomes_dir, macse_path="macse", test_mode=False, occupancy_threshold=0.95):
        self.genomes_dir = Path(genomes_dir).resolve()
        self.macse_path = macse_path
        self.test_mode = test_mode
        self.threshold = occupancy_threshold

        # Paths
        self.meta_path = self.genomes_dir / 'records/compleasm/metadata.csv'
        self.shared_genes_path = self.genomes_dir / 'records/compleasm/shared_single_genes.csv'
        self.fs_report_path = self.genomes_dir / 'records/compleasm/records_frameshifts.csv'
        self.filter_summary_path = self.genomes_dir / 'records/compleasm/records/records_filter_summary.txt'

        # Output Directories
        self.pre_align_dir = self.genomes_dir / 'records/compleasm/alignments/01_pre_alignments'
        self.macse_dir = self.genomes_dir / 'records/compleasm/alignments/02_macse_alignments'

        self.pre_align_dir.mkdir(parents=True, exist_ok=True)
        self.macse_dir.mkdir(parents=True, exist_ok=True)

        # Load Metadata to get total species count
        self.df_meta = pd.read_csv(self.meta_path)
        self.total_species_count = len(self.df_meta)

        with open(self.shared_genes_path, 'r') as f:
            all_orthologs = f.read().strip().split(',')

        if self.test_mode:
            print(f"!!! TEST MODE: Selecting 10 random genes !!!")
            self.orthologs = set(random.sample(all_orthologs, min(10, len(all_orthologs))))
        else:
            self.orthologs = set(all_orthologs)

    def _repath(self, original_path):
        path_obj = Path(original_path)
        if path_obj.exists(): return path_obj
        parts = path_obj.parts
        if 'records' in parts:
            return self.genomes_dir / "/".join(parts[parts.index('records'):])
        return path_obj

    def gather_sequences(self):
        print(f"--- Step 1: Gathering sequences for {len(self.orthologs)} genes ---")
        for gene_id in self.orthologs:
            with open(self.pre_align_dir / f"{gene_id}_unaligned.fasta", 'w') as f: pass

        for _, row in self.df_meta.iterrows():
            acc = row['accession']
            cds_path = self._repath(row['cds_fasta'])
            if not cds_path.exists(): continue
            
            with open(cds_path, 'r') as handle:
                for header, seq in SimpleFastaParser(handle):
                    gene_id = header.split()[0]
                    if gene_id in self.orthologs:
                        with open(self.pre_align_dir / f"{gene_id}_unaligned.fasta", 'a') as f:
                            f.write(f">{acc}\n{seq}\n")

    def run_macse(self):
        print("\n--- Step 2: Running MACSE Alignment ---")
        for fasta_file in self.pre_align_dir.glob("*_unaligned.fasta"):
            gene_id = fasta_file.name.replace("_unaligned.fasta", "")
            out_nt = self.macse_dir / f"{gene_id}_NT.fasta"
            if out_nt.exists(): continue
            
            cmd = [self.macse_path, "-prog", "alignSequences", "-seq", str(fasta_file),
                   "-out_NT", str(out_nt), "-out_AA", str(self.macse_dir / f"{gene_id}_AA.fasta")]
            subprocess.run(cmd, check=True)

    def clean_and_filter(self):
        print(f"\n--- Step 3: Filtering (Threshold: {self.threshold*100}%) ---")
        fs_report_data = []
        dropped_genes = []
        total_processed = 0

        for nt_file in self.macse_dir.glob("*_NT.fasta"):
            gene_id = nt_file.name.replace("_NT.fasta", "")
            if gene_id not in self.orthologs: continue
            total_processed += 1

            clean_species_count = 0 
            clean_sequences = []
            
            with open(nt_file, 'r') as handle:
                for acc, seq in SimpleFastaParser(handle):
                    fs_count = seq.count('!')
                    is_filtered = fs_count >= 2
                    
                    fs_report_data.append({
                        'gene_id': gene_id,
                        'accession': acc,
                        'frameshifts': fs_count,
                        'sequence_filtered': is_filtered
                    })
                    
                    if not is_filtered:
                        clean_species_count += 1
                        clean_sequences.append(f">{acc}\n{seq.replace('!', '-')}\n")

            # Occupancy based on CLEAN sequences
            occupancy = clean_species_count / self.total_species_count
            clean_file = self.macse_dir / f"{gene_id}_clean.fasta"

            if occupancy >= self.threshold:
                with open(clean_file, 'w') as f:
                    f.writelines(clean_sequences)
            else:
                dropped_genes.append(f"{gene_id} (Occupancy: {occupancy:.2f})")
                if clean_file.exists(): clean_file.unlink()

        # Save the detailed Frameshift CSV
        pd.DataFrame(fs_report_data).to_csv(self.fs_report_path, index=False)

        # Save the Summary TXT Report
        with open(self.filter_summary_path, 'w') as f:
            f.write(f"MACSE FILTER SUMMARY\n")
            f.write(f"====================\n")
            f.write(f"Total Genes Analyzed: {total_processed}\n")
            f.write(f"Occupancy Threshold: {self.threshold*100}%\n")
            f.write(f"Genes Dropped: {len(dropped_genes)}\n")
            f.write(f"Genes Retained: {total_processed - len(dropped_genes)}\n\n")
            f.write(f"LIST OF DROPPED GENE IDS:\n")
            for item in dropped_genes:
                f.write(f"- {item}\n")

        print(f"Detailed frameshift report: {self.fs_report_path}")
        print(f"Filter summary report: {self.filter_summary_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("genomes_dir")
    parser.add_argument("--test", action="store_true")
    parser.add_argument("--threshold", type=float, default=0.95)
    args = parser.parse_args()

    pipeline = MacsePipeline(args.genomes_dir, test_mode=args.test, occupancy_threshold=args.threshold)
    pipeline.gather_sequences()
    pipeline.run_macse()
    pipeline.clean_and_filter()