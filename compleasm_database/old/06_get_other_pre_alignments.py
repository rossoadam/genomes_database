#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse

class other_pre_align_b:
    def __init__(self, genomes_dir, targets):
        self.genomes_dir = Path(genomes_dir)
        self.targets = targets
        
        # Metadata Paths
        self.compleasm_meta_path = self.genomes_dir / 'records/compleasm/metadata.csv'
        self.genomes_meta_path = self.genomes_dir / 'records/genomes_metadata.csv'
        self.shared_genes_path = self.genomes_dir / 'records/compleasm/shared_single_genes.csv'
        
        # Output Directory
        self.pre_dir = self.genomes_dir / 'records/compleasm/alignments/01_pre_alignments'
        
        # Load shared orthologs
        if not self.shared_genes_path.exists():
            raise FileNotFoundError(f"Missing {self.shared_genes_path}. Run setup script first.")
            
        with open(self.shared_genes_path, 'r') as f:
            self.orthologs = set(f.read().strip().split(','))

    def _repath(self, original_path):
        """Fixes paths if drive mount point has changed."""
        path_obj = Path(original_path)
        if path_obj.exists(): return path_obj
        try:
            parts = path_obj.parts
            idx = parts.index('records')
            return self.genomes_dir / Path(*parts[idx:])
        except ValueError:
            return path_obj

    def _auto_detect_reference(self):
        """Peeks at existing ref files to find the accession to exclude."""
        try:
            # Look for the first directory that actually contains a ref_ file
            for gene_dir in self.pre_dir.iterdir():
                if gene_dir.is_dir():
                    ref_files = list(gene_dir.glob("ref_*.fasta"))
                    if ref_files:
                        with open(ref_files[0], 'r') as f:
                            header = f.readline().strip()
                            return header.replace(">", "")
        except Exception as e:
            print(f"Warning: Could not auto-detect reference: {e}")
        return None

    def run(self):
        # 1. Detect Reference
        ref_accession = self._auto_detect_reference()
        if not ref_accession:
            print("CRITICAL: No reference files found in 01_pre_alignments. Run Script 05 first.")
            return

        # 2. Load and Merge Metadata
        # We need assembly level from genomes_meta and cds_fasta from compleasm_meta
        c_df = pd.read_csv(self.compleasm_meta_path)
        g_df = pd.read_csv(self.genomes_meta_path)
        
        # Create a lookup for assembly level: {accession: level}
        level_map = dict(zip(g_df['accession'], g_df['assembly_level']))

        # 3. Filter for 'Others' that match target quality
        # We exclude the reference and genomes that don't meet assembly targets
        print(f"Filtering genomes. Targets: {self.targets}")
        
        included_count = 0
        for _, row in c_df.iterrows():
            acc = row['accession']
            
            # Skip if it's the reference
            if acc == ref_accession:
                continue
                
            # Check assembly level from the genomes_metadata mapping
            level = level_map.get(acc)
            if level not in self.targets:
                continue

            cds_path = self._repath(row['cds_fasta'])
            if not cds_path.exists():
                print(f"Warning: CDS for {acc} missing at {cds_path}")
                continue

            print(f"Processing {acc} ({level})...")
            included_count += 1
            
            # 4. Parse CDS and append to gene folders
            with open(cds_path, 'r') as handle:
                for header, seq in SimpleFastaParser(handle):
                    gene_id = header.split()[0]
                    
                    if gene_id in self.orthologs:
                        gene_dir = self.pre_dir / gene_id
                        
                        if gene_dir.exists():
                            out_file = gene_dir / f"oth_{gene_id}.fasta"
                            # Append so multiple 'others' end up in the same file
                            with open(out_file, 'a') as f:
                                f.write(f">{acc}\n{seq}\n")

        print(f"\nDone. Processed {included_count} genomes into {self.pre_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("genomes_dir", help="Project root directory")
    parser.add_argument("targets", nargs='+', help="Assembly levels to include (e.g. Complete Chromosome)")
    
    args = parser.parse_args()

    worker = other_pre_align_b(args.genomes_dir, args.targets)
    worker.run()