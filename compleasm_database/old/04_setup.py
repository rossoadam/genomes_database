#!/usr/bin/env python
# -*- coding: utf-8 -*-
# this script is intended to organize the busco_run_results for my pipeline

import os
import sys
import pandas as pd
import csv
import sys
from pathlib import Path
import argparse

#######
#   make sure the setup.mask() is uncommented at the end of the file
#   example usage:
#   
#######

class get_single_genes_b:
   def __init__(self, genomes_dir, reference_input, targets):
      #### INPUTS ####
      self.genomes_dir = Path(genomes_dir).resolve()
      # this the input the user specifies they want to be the reference alignment for kc-align
      self.reference_input = reference_input
      # this is the target assembly level, choose one or more of the following Complete Chromosome Scaffold
      self.targets = targets
      #### PATHS ####
      self.busco_id_dir = self.genomes_dir / 'records/compleasm/mb_downloads/sauropsida_odb12/hmms/'
      self.genomes_meta = self.genomes_dir / 'records/genomes_metadata.csv'
      self.compleasm_meta = self.genomes_dir / 'records/compleasm/metadata.csv'
      # load the master list of busco IDs from the HMM directory
      if self.busco_id_dir.exists():
         self.busco_id_list = sorted([f.stem for f in self.busco_id_dir.glob("*.hmm")])
      else:
         print(f"Error: HMM directory not found at {self.busco_id_dir}")
         self.busco_id_list = []

      # this will be the list complete shared orthologs
      self.single_busco_ids=[]
      # this will be a dictionary with the accession as key and level of the assembly as value
      self.assembly = {}
      # this will be a dictionary with the accession as key and cds_fasta as value
      self.cds_fasta = {}

   def _repath(self, path_to_fix):
      # helper method to swap local paths for current working direcoty path
      path_str = str(path_to_fix)
      # find where genomes starts in the path
      root_name = self.genomes_dir.name # this should just be genomes
      if root_name in path_str:
        # path_str.split(root_name)[-1] --- takes everything to the right of genomes in /Users/rossoaa/projects/genomes/records/compleasm/GCF_004329235.1__Podarcis_muralis/sauropsida_odb12/full_table.tsv 
         suffix = path_str.split(root_name)[-1].lstrip('/\\')
         return self.genomes_dir / suffix
      return Path(path_str)

   # this functions does pattern matching to identify reference alignment
   def retrieve_reference(self):
     # input: genus_species OR accession
     # returns: the CDS fasta file path from compleasm_metadata.csv as string
     try:
        df = pd.read_csv(self.compleasm_meta)
        query = self.reference_input.lower().replace(" ", "_")
        # try to match by accession first
        if 'accession' in df.columns:
           match = df[df['accession'].astype(str).str.lower() == query]
           if not match.empty:
              return match['cds_fasta'].iloc[0]
        # then check genus_species
        if 'genus_species' in df.columns:
           df['norm_name'] = df['genus_species'].astype(str).str.lower().str.replace(" ", "_")
           match = df[df['norm_name'] == query]
           if not match.empty:
              return match['cds_fasta'].iloc[0]
     except Exception as e:
        print(f"Error in retrieve_reference: {e}")
     return None

   def mask(self): # this method reads all the tsv files from compleasm and writes a shared_complete.csv
      if not self.busco_id_list:
         print("no busco IDs found in HMM directory, aborting.")
         return
      # load metadata files
      genomes_df = pd.read_csv(self.genomes_meta)
      compleasm_df = pd.read_csv(self.compleasm_meta)
      
      # create a dictionary of accession,assembly
      # this allows me to filter the compleasm metadata by quality
      acc_to_level_dict = dict(zip(genomes_df['accession'], genomes_df['assembly_level']))
      
      # initialize the shared IDs set with all possible IDs from your HMMs
      shared_single_set = set(self.busco_id_list)

      included_count = 0

      # iterate through the compleasm metadata
      for _, row in compleasm_df.iterrows():
         acc = row['accession']
         tsv_path_str = row['full_table']

         # 1. filter by assembly level
         level = acc_to_level_dict.get(acc)
         if level not in self.targets:
            continue
         # 2. check if the tsv path exists
         tsv_path = self._repath(tsv_path_str)
         if not tsv_path.exists():
            print(f"Warning: TSV file missing for {acc} at {tsv_path}")
            continue
         # 3. Process the TSV
         try:
            # Compleasm TSVs usually have headers starting with #
            # we only need the first two columns Busco id and Status
            df_tsv = pd.read_csv(tsv_path, sep='\t', comment='#', header=0, usecols=[0,1], names=['Gene','Status'])
            # get IDs that are 'Single', other options are missing and duplicated. I think single implies complete?
            current_single = set(df_tsv[df_tsv['Status'] == 'Single']['Gene'])
            # intersection: keep only IDs that have been complete in every genome seen so far
            shared_single_set &= current_single

            print(f"Processed {acc}: {level}")
            included_count += 1
         except Exception as e:
            print(f"Error processing TSV for {acc}: {e}")
      self.single_busco_ids = sorted(list(shared_single_set))

      print(f"\nFinished. Filtered for assembly levels: {self.targets}")
      print(f"Genomes analyzed: {included_count}")
      print(f"Shared 'Single' Genes: {len(self.single_busco_ids)}")

      # save to CSV
      output_path = self.genomes_dir / 'records/compleasm/shared_single_genes.csv'
      with open(output_path, 'w', newline='') as f:
         writer = csv.writer(f)
         writer.writerow(self.single_busco_ids)

      print(f"Shared IDs written to: {output_path}")
      return self.single_busco_ids

if __name__ == "__main__":
   parser = argparse.ArgumentParser(description="Generate, shared single genes list. Eventually organize compleasm run results")
   # positional arguments
   parser.add_argument("genomes_dir", help="Path to genomes directory")
   parser.add_argument("reference_input", help="genus_species OR accession # - GC*####### for reference alingment input")
   # here is where I can add 1 or more arguments nargs='+', help="Assembly targets (e.g., Chromsome Scaffold)" 
   parser.add_argument("targets", nargs='+', help="Please specify targeted assembly level e.g. complete, scaffold, chromosome")
   args = parser.parse_args()
   # pass the parsed arguments to the class
   get_genes = get_single_genes_b(args.genomes_dir, args.reference_input, args.targets)
   # print(setup.retrieve_reference())
   # get_genes.mask()
