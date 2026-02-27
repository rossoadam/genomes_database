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
        self.report_path = self.genomes_dir / 'records/compleasm/records_frameshifts.csv'

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
            self.orthologs = set(random.sample(all_orthologs, min(30, len(all_orthologs))))
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
        report_data = []

        for nt_file in self.macse_dir.glob("*_NT.fasta"):
            gene_id = nt_file.name.replace("_NT.fasta", "")
            if gene_id not in self.orthologs: continue

            species_in_gene = 0
            clean_sequences = []
            
            with open(nt_file, 'r') as handle:
                for acc, seq in SimpleFastaParser(handle):
                    species_in_gene += 1
                    fs_count = seq.count('!')
                    
                    # Track if THIS specific sequence is "bad"
                    is_filtered = fs_count >= 2
                    
                    report_data.append({
                        'gene_id': gene_id,
                        'accession': acc,
                        'frameshifts': fs_count,
                        'sequence_filtered': is_filtered
                    })
                    
                    if not is_filtered:
                        clean_sequences.append(f">{acc}\n{seq.replace('!', '-')}\n")

            # Occupancy Check
            occupancy = species_in_gene / self.total_species_count
            clean_file = self.macse_dir / f"{gene_id}_clean.fasta"

            if occupancy >= self.threshold:
                print(f"Keeping {gene_id}: Occupancy {occupancy:.2f}")
                with open(clean_file, 'w') as f:
                    f.writelines(clean_sequences)
            else:
                print(f"Dropping {gene_id}: Occupancy {occupancy:.2f} too low.")
                if clean_file.exists(): clean_file.unlink()

        pd.DataFrame(report_data).to_csv(self.report_path, index=False)
        print(f"Report saved to {self.report_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("genomes_dir")
    parser.add_argument("--test", action="store_true")
    parser.add_argument("--threshold", type=float, default=0.95, help="Occupancy threshold (0.0 - 1.0)")
    args = parser.parse_args()

    pipeline = MacsePipeline(args.genomes_dir, test_mode=args.test, occupancy_threshold=args.threshold)
    pipeline.gather_sequences()
    pipeline.run_macse()
    pipeline.clean_and_filter()