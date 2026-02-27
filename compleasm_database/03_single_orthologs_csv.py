#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pandas as pd
import csv
from pathlib import Path
import argparse
from collections import Counter
################
# EXAMPLE USAGE:
# 
# 1. Standard Run (80% occupancy, only Chromosome and Scaffold assemblies):
#    python 04_single_orthologs_csv.py ../../../genomes "Complete Chromosome" "Scaffold"
# 
# 2. Strict Run (100% occupancy for a high-quality phylogenomic matrix):
#    python 04_single_orthologs_csv.py ../../../genomes "Complete Chromosome" --threshold 1.0
# 
# 3. Relaxed Run (60% occupancy to maximize gene count):
#    python 04_single_orthologs_csv.py ../../../genomes "Complete Chromosome" "Scaffold" --threshold 0.6
# 
# ARGUMENTS:
#   genomes_dir    Path to the root project directory containing 'records/'.
#   targets        One or more assembly levels to include. Common values:
#                  "Complete Chromosome", "Scaffold", "Contig".
#                  (Separate multiple targets with spaces).
# 
# OPTIONAL FLAGS:
#   --threshold    The fraction of genomes (0.0 - 1.0) that must have a gene
#                  marked as 'Single' to include it. Default is 0.80 (80%).
################
class get_single_genes_b:
    def __init__(self, genomes_dir, targets, threshold_pct=0.80):
        self.genomes_dir = Path(genomes_dir).resolve()
        self.targets = targets
        self.threshold_pct = threshold_pct
        
        #### PATHS ####
        self.gene_id_dir = self.genomes_dir / 'records/compleasm/mb_downloads/sauropsida_odb12/hmms/'
        self.genomes_meta = self.genomes_dir / 'records/genomes_metadata.csv'
        self.compleasm_meta = self.genomes_dir / 'records/compleasm/records/metadata.csv'
        
        if self.gene_id_dir.exists():
            self.gene_id_list = sorted([f.stem for f in self.gene_id_dir.glob("*.hmm")])
        else:
            print(f"Error: HMM directory not found at {self.gene_id_dir}")
            self.gene_id_list = []

        self.single_gene_ids = []

    def _repath(self, path_to_fix):
        path_str = str(path_to_fix)
        root_name = self.genomes_dir.name
        if root_name in path_str:
            suffix = path_str.split(root_name)[-1].lstrip('/\\')
            return self.genomes_dir / suffix
        return Path(path_str)

    def mask(self):
        if not self.gene_id_list:
            print("no busco IDs found in HMM directory, aborting.")
            return
        
        genomes_df = pd.read_csv(self.genomes_meta)
        compleasm_df = pd.read_csv(self.compleasm_meta)
        acc_to_level_dict = dict(zip(genomes_df['accession'], genomes_df['assembly_level']))
        
        # Use a Counter to track how many genomes have each gene as 'Single'
        gene_occurrence_counter = Counter()
        included_count = 0

        for _, row in compleasm_df.iterrows():
            acc = row['accession']
            level = acc_to_level_dict.get(acc)
            
            if level not in self.targets:
                continue

            tsv_path = self._repath(row['full_table'])
            if not tsv_path.exists():
                print(f"Warning: TSV file missing for {acc}")
                continue

            try:
                df_tsv = pd.read_csv(tsv_path, sep='\t', comment='#', header=0, usecols=[0,1], names=['Gene','Status'])
                
                # Get genes that are 'Single' in this specific genome
                current_single = df_tsv[df_tsv['Status'] == 'Single']['Gene'].tolist()
                
                # Increment counts for these genes
                gene_occurrence_counter.update(current_single)

                print(f"Processed {acc}: {level}")
                included_count += 1
            except Exception as e:
                print(f"Error processing {acc}: {e}")

        # Calculate the required number of genomes based on the threshold
        required_count = int(included_count * self.threshold_pct)
        
        # Filter genes that met the threshold
        self.single_gene_ids = sorted([
            gene for gene, count in gene_occurrence_counter.items() 
            if count >= required_count and gene in self.gene_id_list
        ])

        print(f"\nFinished. Filtered for assembly levels: {self.targets}")
        print(f"Occupancy Threshold: {self.threshold_pct*100}% ({required_count}/{included_count} genomes)")
        print(f"Genes Meeting Threshold: {len(self.single_gene_ids)}")

        output_path = self.genomes_dir / f"records/compleasm/records/shared_single_genes.csv"
        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(self.single_gene_ids)

        print(f"Shared IDs written to: {output_path}")
        return self.single_gene_ids

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("genomes_dir", help="Path to genomes directory")
    parser.add_argument("targets", nargs='+', help="Target assembly levels")
    parser.add_argument("--threshold", type=float, default=0.80, help="Occupancy threshold (e.g. 0.80 for 80 percent)")
    
    args = parser.parse_args()
    
    # Updated class call (removed reference_input as it wasn't used in your class logic)
    get_genes = get_single_genes_b(args.genomes_dir, args.targets, threshold_pct=args.threshold)
    get_genes.mask()