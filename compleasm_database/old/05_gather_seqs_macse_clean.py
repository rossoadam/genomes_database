#!/usr/bin/env python3
import os
import subprocess
import pandas as pd
import argparse
import random
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser

class MacsePipeline:
    def __init__(self, genomes_dir, macse_path="macse", test_mode=False):
        self.genomes_dir = Path(genomes_dir).resolve()
        self.macse_path = macse_path
        self.test_mode = test_mode

        # Consistent Pathing
        self.meta_path = self.genomes_dir / 'records/compleasm/metadata.csv'
        self.shared_genes_path = self.genomes_dir / 'records/compleasm/shared_single_genes.csv'

        # Output Directories
        self.pre_align_dir = self.genomes_dir / 'records/compleasm/alignments/01_pre_alignments'
        self.macse_dir = self.genomes_dir / 'records/compleasm/alignments/02_macse_alignments'

        self.pre_align_dir.mkdir(parents=True, exist_ok=True)
        self.macse_dir.mkdir(parents=True, exist_ok=True)

        # Load shared orthologs
        if not self.shared_genes_path.exists():
            raise FileNotFoundError(f"Could not find shared genes file at {self.shared_genes_path}")

        with open(self.shared_genes_path, 'r') as f:
            all_orthologs = f.read().strip().split(',')

        # TEST MODE LOGIC: Reduce gene list to 10 if flag is set
        if self.test_mode:
            print(f"!!! TEST MODE ACTIVE: Selecting 10 random genes from {len(all_orthologs)} total !!!")
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
        df_comp = pd.read_csv(self.meta_path)

        for gene_id in self.orthologs:
            with open(self.pre_align_dir / f"{gene_id}_unaligned.fasta", 'w') as f:
                pass

        for _, row in df_comp.iterrows():
            acc = row['accession']
            cds_path = self._repath(row['cds_fasta'])

            if not cds_path.exists():
                continue

            print(f"Extracting from {acc}...")
            with open(cds_path, 'r') as handle:
                for header, seq in SimpleFastaParser(handle):
                    gene_id = header.split()[0]
                    if gene_id in self.orthologs:
                        out_file = self.pre_align_dir / f"{gene_id}_unaligned.fasta"
                        with open(out_file, 'a') as f:
                            f.write(f">{acc}\n{seq}\n")

    def run_macse(self):
    		print("\n--- Step 2: Running MACSE Alignment ---")
    		for fasta_file in self.pre_align_dir.glob("*_unaligned.fasta"):
    			gene_id = fasta_file.name.replace("_unaligned.fasta", "")
    			out_nt = self.macse_dir / f"{gene_id}_NT.fasta"
    			
    			if out_nt.exists(): continue
    
    			# Simplified command using the Conda wrapper
    			cmd = [
    				"macse", 
    				"-prog", "alignSequences",
    				"-seq", str(fasta_file),
    				"-out_NT", str(out_nt),
    				"-out_AA", str(self.macse_dir / f"{gene_id}_AA.fasta")
    			]
    			
    			print(f"Aligning {gene_id}...")
    			# Ensure the conda environment is active when running this
    			subprocess.run(cmd, check=True)

    def clean_for_iqtree(self):
        print("\n--- Step 3: Cleaning for IQ-TREE ---")
        for nt_file in self.macse_dir.glob("*_NT.fasta"):
            gene_id = nt_file.name.replace("_NT.fasta", "")
            if gene_id not in self.orthologs:
                continue

            clean_file = self.macse_dir / nt_file.name.replace("_NT.fasta", "_clean.fasta")
            with open(nt_file, 'r') as f_in, open(clean_file, 'w') as f_out:
                f_out.write(f_in.read().replace("!", "-"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Consolidate sequences and align with MACSE")
    parser.add_argument("genomes_dir", help="Path to the root project directory")
    # parser.add_argument("--macse_jar", default="/usr/local/bin/macse_v2.06.jar", help="Path to MACSE jar")
    parser.add_argument("--test", action="store_true", help="Run only 10 random genes for testing")

    args = parser.parse_args()

    pipeline = MacsePipeline(args.genomes_dir, test_mode=args.test)
    pipeline.gather_sequences()
    pipeline.run_macse()
    pipeline.clean_for_iqtree()
