#!/usr/bin/env python3
import subprocess
import pandas as pd
import argparse
import random
import logging
import shutil
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser
from concurrent.futures import ProcessPoolExecutor

class MacseAligner:
    def __init__(self, genomes_dir, threads=1, test_mode=False):
        self.genomes_dir = Path(genomes_dir).resolve()
        self.threads = threads
        self.test_mode = test_mode
        
        # Paths
        self.records_dir = self.genomes_dir / 'records/compleasm/records'
        self.align_base = self.genomes_dir / 'records/compleasm/alignments'
        self.meta_path = self.records_dir / 'metadata.csv'
        self.shared_genes_path = self.records_dir / 'shared_single_genes_0.8.csv'
        self.pre_align_dir = self.align_base / '01_pre_alignments'
        self.macse_dir = self.align_base / '02_macse_alignments'
        self.log_file = self.records_dir / 'records_macse.log'

        # Ensure directories exist
        self.pre_align_dir.mkdir(parents=True, exist_ok=True)
        self.macse_dir.mkdir(parents=True, exist_ok=True)

        # Setup Logging
        logging.basicConfig(
            filename=self.log_file,
            level=logging.ERROR,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )

        with open(self.shared_genes_path, 'r') as f:
            all_orthologs = f.read().strip().split(',')
        
        # Select orthologs
        if test_mode:
            self.orthologs = set(random.sample(all_orthologs, 300))
            print(f"!!! TEST MODE ACTIVE: Only processing {len(self.orthologs)} genes !!!")
        else:
            self.orthologs = set(all_orthologs)

    def _repath(self, original_path):
        path_obj = Path(original_path)
        if path_obj.exists(): return path_obj
        parts = path_obj.parts
        return self.genomes_dir / "/".join(parts[parts.index('records'):]) if 'records' in parts else path_obj

    def _run_single_macse(self, gene_id):
        """Helper function to run MACSE on a specific gene ID."""
        fasta = self.pre_align_dir / f"{gene_id}_unaligned.fasta"
        out_nt = self.macse_dir / f"{gene_id}_NT.fasta"
        out_aa = self.macse_dir / f"{gene_id}_AA.fasta"
        
        if not fasta.exists():
            return

        if not out_nt.exists():
            try:
                subprocess.run([
                    "macse", "-prog", "alignSequences", 
                    "-seq", str(fasta), 
                    "-out_NT", str(out_nt), 
                    "-out_AA", str(out_aa)
                ], check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                error_msg = f"Gene: {gene_id} | Error: {e.stderr.strip()}"
                logging.error(error_msg)

    def run(self):
        # 1. Clean up old unaligned files if in test mode to avoid confusion
        if self.test_mode:
            print("Cleaning up old test files...")
            for f in self.pre_align_dir.glob("*.fasta"):
                f.unlink()

        # 2. Gathering Sequences
        print(f"--- Gathering Sequences for {len(self.orthologs)} genes ---")
        df_meta = pd.read_csv(self.meta_path)
        
        # Initialize files
        for g in self.orthologs:
            with open(self.pre_align_dir / f"{g}_unaligned.fasta", 'w') as f: pass

        for _, row in df_meta.iterrows():
            cds_path = self._repath(row['cds_fasta'])
            if not cds_path.exists(): continue
            with open(cds_path, 'r') as h:
                for head, seq in SimpleFastaParser(h):
                    gene_id = head.split()[0]
                    if gene_id in self.orthologs:
                        with open(self.pre_align_dir / f"{gene_id}_unaligned.fasta", 'a') as f:
                            f.write(f">{row['accession']}\n{seq}\n")

        # 3. Running MACSE
        print(f"--- Running MACSE in parallel (Threads: {self.threads}) ---")
        
        # Map over the specific list of orthologs we just gathered
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            executor.map(self._run_single_macse, list(self.orthologs))
        
        print(f"--- Finished. Check {self.log_file} for any failures. ---")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("genomes_dir")
    parser.add_argument("--test", action="store_true")
    parser.add_argument("-j", "--threads", type=int, default=1)
    args = parser.parse_args()
    
    MacseAligner(args.genomes_dir, threads=args.threads, test_mode=args.test).run()