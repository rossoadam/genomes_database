#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Edited on Tuesday Feb 24 09:14:00 2026

@author: adam
"""
####
#
####

import os
import pandas as pd
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser

class ref_pre_align_b:
    def __init__(self, genomes_dir, reference_input):
        self.genomes_dir = Path(genomes_dir)
        self.reference_input = reference_input
        
        # paths
        self.compleasm_meta_path = self.genomes_dir / 'records/compleasm/metadata.csv'
        self.shared_genes_path = self.genomes_dir / 'records/compleasm/shared_single_genes.csv'
        self.alignments_dir = self.genomes_dir / 'records/compleasm/alignments'
        self.pre_align_dir = self.alignments_dir / '01_pre_alignments'
        
        # load the shared orthologs (assuming it's a single row of IDs)
        # making it into a set should make it faster?
        with open(self.shared_genes_path, 'r') as f:
            content = f.read().strip().split(',') # strip() removes white space tabs newlines
            self.orthologs = set(content)

    def _repath(self, original_path):
        # Fixes paths from metadata if the drive was mounted differently when the metadata was originally created.
        path_obj = Path(original_path)
        # If the file exists as is, just return it
        if path_obj.exists():
            return path_obj
        # If not, we find the 'records' folder in the path and 
        # re-attach it to our current genomes_dir
        try:
            # This finds the index of 'records' in the old path parts
            parts = path_obj.parts # input path_obj to .parts and returns a tuple of strings (/, usr, bin)
            records_index = parts.index('records') # index returns the position in the tuple
            new_path = self.genomes_dir / Path(*parts[records_index:]) # * unpacks everything in parts, then sliced index to the end
            return new_path
        except ValueError:
            # 'records' wasn't in the path, return original and hope for the best
            return path_obj

    def check(self):
        # Does the directory structure exist?
        dir_checklist = [
            self.alignments_dir,
            self.pre_align_dir,
            self.alignments_dir / '02_consolidate'
        ]
        
        for folder in dir_checklist:
            if not folder.exists():
                print(f"Creating directory: {folder}")
                folder.mkdir(parents=True, exist_ok=True)
        print("Setup check complete.")

    def get_ref_data(self):
        # Retrieves the reference accession and its CDS path from metadata.
        df = pd.read_csv(self.compleasm_meta_path)
        query = self.reference_input.lower().replace(" ", "_")
        
        # Check accession or genus_species
        match = df[(df['accession'].astype(str).str.lower() == query) | 
                   (df['genus_species'].astype(str).str.lower().str.replace(" ", "_") == query)]
        
        if match.empty:
            raise ValueError(f"Reference '{self.reference_input}' not found in metadata.")
            
        return match.iloc[0]['accession'], Path(match.iloc[0]['cds_fasta'])

    def run(self):
        # Creates subdirectories and writes ref_[ID].fasta files.
        self.check()
        ref_acc, ref_cds_path = self.get_ref_data()
        
        if not ref_cds_path.exists():
            print(f"**WARNING: Reference CDS file not found at {ref_cds_path}, trying to fix..")
            ref_cds_path = self._repath(ref_cds_path)
            if not ref_cds_path.exists():
               print(f"rebuild unsuccessful: {ref_cds_path}")

        print(f"Processing reference: {ref_acc}")
        print(f"Reading from: {ref_cds_path}")

        # Optimization: We parse the reference FASTA once and pick out the IDs we need
        count = 0
        with open(ref_cds_path, 'r') as handle:
            for header, seq in SimpleFastaParser(handle):
                # The first word in the header is typically the BUSCO ID
                gene_id = header.split()[0]
                
                if gene_id in self.orthologs:
                    # 1. Create the gene-specific subdirectory
                    gene_dir = self.pre_align_dir / gene_id
                    gene_dir.mkdir(exist_ok=True)
                    
                    # 2. Define the output file path
                    out_file = gene_dir / f"ref_{gene_id}.fasta"
                    
                    # 3. Write the file with the reformatted header (Accession only)
                    with open(out_file, 'w') as f:
                        f.write(f">{ref_acc}\n{seq}\n")
                    
                    count += 1

        print(f"Finished. Created {count} reference FASTA files in {self.pre_align_dir}")

if __name__ == "__main__":
    # Example usage:
    # python 05_get_ref_pre-alignments.py /path/to/project "Podarcis muralis"
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("genomes_dir")
    parser.add_argument("reference")
    args = parser.parse_args()

    processor = ref_pre_align_b(args.genomes_dir, args.reference)
    processor.run()