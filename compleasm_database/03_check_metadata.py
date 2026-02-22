#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
this code is intended to figure out which genomes are not represented in the meta_genomes_2.csv file
  or to find out which genomes may not have been analyzed by busco
'''
import os
import pandas as pd
import csv
import sys

class checker_b:
   def __init__(self):
      self.genomes_dir = sys.argv[1]
      self.genomes_metadata = os.path.join(self.genomes_dir,'records','genomes_metadata.csv')
      self.accession_list = []
      self.compleasm_metadata = os.path.join(self.genomes_dir,'records/compleasm/','metadata.csv')
   def dif(self):
      # open the metadata file and assign to meta
      # print(self.genomes_metadata, self.compleasm_metadata)
      with open(self.genomes_metadata, 'r') as meta, open(self.compleasm_metadata, 'r') as meta_2:
         genomes_meta = pd.read_csv(meta)
         genomes_meta_accessions = list(genomes_meta['accession'])
         # print(genomes_meta_accessions)
         # print(type(genomes_meta_accessions[0]))
         compleasm_meta = pd.read_csv(meta_2)
         compleasm_meta_accessions = list(compleasm_meta['accession'])
         # print(compleasm_meta_accessions)
         # print(type(genomes_meta_accessions[0]))
      for name in genomes_meta_accessions:
         # print(name)
         if 'GC' in name:
            check_this = name[:15]
            # print(type(check_this))
            if check_this not in compleasm_meta_accessions:
               print("*WARNING*", name, "NOT in COMPLEASM metadata.csv")
               pass
checker = checker_b()
checker.dif()
