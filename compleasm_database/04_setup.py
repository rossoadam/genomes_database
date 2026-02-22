#!/usr/bin/env python
# -*- coding: utf-8 -*-
# this script is intended to organize the busco_run_results for my pipeline

import os
import sys
import pandas as pd
import csv
import sys
from pathlib import Path

#######
# When running on lepidodactylus update with locations on local computer:
#     self.busco
#     self.busco_id_csv
#     self.meta
#     make sure the setup.mask() is uncommented at the end of the file
# the check method needs to be fixed, as of now all buscos get moved to the other folder
#######

class setup_b:
   def __init__(self):
      self.genomes_dir = sys.argv[1]
      # this should be the compleasm folder containing the lineage refseq db
      self.busco_id_dir = Path(self.genomes_dir) / 'records/compleasm/mb_downloads/sauropsida_odb12/hmms/'
      if not self.busco_id_dir.exists():
         print(f"Error: could not find hmm directory at {self.busco_id_dir}")
      self.busco_id_list = [f.stem for f in self.busco_id_dir.glob("*.hmm")]
      # this is the metadata for all genomes - actually I need to add data for busco 03 and 04...
      self.genomes_meta = Path(self.genomes_dir) / 'records/genomes_metadata.csv'
      # this is the metadata for compleasm results
      self.compleasm_metadata = Path(self.genomes_dir) / 'records/compleasm/metadata.csv'
      # this is the directory that will contain ouputs from alignments
      self.alignments = Path(self.genomes_dir) / 'records/compleasm/alignments/'
      # this is the busco output folder for p muralis
      ### edit this self.ref_dir = 'GCF_004329235.1_PodMur_1.0_genomic.fna' # I'm not sure if this needs to be the cds fasta or the full genome podarcis_muralis_cds_compleasm.fasta
      # all other squamate busco results
      ### edit this self.oth_dir = 'GC*_*'
      # this will be the list complete shared orthologs
      self.complete_busco_ids=[]
      # more dirs
      self.dir_checklist = [
         self.alignments / 'reference',
         self.alignments / 'other',
         self.alignments / '01_pre_alignments',
         self.alignments / '02_consolidate']
      # this is the level of the assembly
      self.assembly = {}
      # this is the target assembly level, choose one or more of the following Complete Chromosome Scaffold
      self.target = [sys.argv[2],sys.argv[3],sys.argv[4]]

   def check(self):
      contents = list(self.alignments.iterdir())
      setup_is_complete = all(item in contents for item in self.dir_checklist)
      if setup_is_complete:
         print('setup is complete:',setup_is_complete)
      if not setup_is_complete:
         print('setup is NOT complete:', setup_is_complete)
         for item in self.dir_checklist:
            if item not in contents:
               # Create  a directory for the missing item
               print('creating:', item)
               cmd = f'mkdir {item}'
               os.system(cmd)
      contents = list(self.alignments.iterdir())
      # the following loop moves all busco output files to either reference or other
      for item in contents:
         if item == 'reference':
            if os.listdir(os.path.join(self.busco,item)) == []:
               cmd_ref_mv = f'mv {self.ref_dir} reference'
               os.system(cmd_ref_mv)
         if item == 'other':
            if os.listdir(os.path.join(self.busco,item)) == []:
               cmd_oth_mv = f'mv {self.oth_dir} other'
               os.system(cmd_oth_mv)
   def mask(self): # this method reads all the tsv files from busco and writes a shared_complete.csv
      accession = None
      with open(self.meta, 'r') as  meta:
         meta = pd.read_csv(meta)
         meta_acc = list(meta['accession'])
         meta_ass = list(meta['asemblyLevel'])
      for i in enumerate(meta_acc):
         self.assembly[i[1]] = meta_ass[i[0]]
      with open(self.busco_id_csv, 'r') as busco_ids:
         # read csv as a pandas dataframe
         busco_ids = pd.read_csv(busco_ids)
         # update this variable to only include IDs
         busco_ids = busco_ids['# Busco id']

      # Create a boolean mask with all True assuming all IDs are present initially
      busco_ids_mask = [True] * len(busco_ids)
      # Convert to a panda Series. This prevents warnings when updating the mask
      busco_ids_mask = pd.Series( (busco_ids_mask) )
      # os walk down the directory tree and look for full tsv files
      test = []
      for (root, dirs, file) in os.walk(self.busco, topdown=True):
         # these lines grab the accession from the root...
         text = root
         delimiter_1 = 'GC'
         delimiter_2 = '.'
         index_1 = text.find(delimiter_1)
         index_2 = text.find(delimiter_2)
         if index_1 != -1:
            accession= text[index_1:index_2+2]
         if accession == None:
            pass
         elif accession in meta_acc and self.assembly[accession] in self.target: #... here I use the accession and assembly level to filter which genomes are used
            for name in file:
               # if the name endswith tsv and contains full
               if name.endswith('.tsv') and 'full' in name:
                  # get the pathway - here join intelligently, preserves path
                  path_n_name = os.path.join(root,name)
                  # test.append(path_n_name) a troubleshooting line
                  # now I need to open the tsv and get a list of complete BUSCO
                  with open(path_n_name, 'r') as tsv:
                     tsv = pd.read_csv(tsv, sep='\t', header=2)
                     # remove dulpicates
                     tsv = tsv.drop_duplicates(subset='# Busco id')
                     # store the column as a boolean array
                     status = pd.Series(tsv.Status)
                     # this turns every instance to either True or False, based on presence of Complete
                     status = status == 'Complete'
                     # update mask using logical AND
                     # here I needed to include zip so that both values are evaluated at the same time
                     busco_ids_mask = [mask and stat for mask, stat in zip(status, busco_ids_mask)]
                     #print(sum(busco_ids_mask))

      # filter busco_ids based on the final mask
      self.complete_busco_ids = [busco_id for busco_id, is_complete in zip(busco_ids, busco_ids_mask) if is_complete]
      print('there are '+ str(len(self.complete_busco_ids)) + ' complete busco IDs across these species')
      os.chdir(self.busco)
      with open('shared_complete.csv', 'w') as f:
         writer=csv.writer(f)
         writer.writerow(self.complete_busco_ids)
      return (self.complete_busco_ids)

setup = setup_b()
#setup.check()
# the following command writes a csv file ot the busco directory
#setup.mask()
