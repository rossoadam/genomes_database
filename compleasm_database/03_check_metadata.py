#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
this code is intended to figure out which genomes are not represented in the meta_genomes_2.csv file
  or to find out which genomes may not have been analyzed by busco
'''
import os
import pandas as pd
import csv

class checker_b:
   def __init__(self):
      self.genomes_2_metadata = '/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/genomes_2/metadata_genomes_2.csv'
      self.accession_list = []
      self.busco_dir = '/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/genomes_2/busco_run_05_results (copy)'
   def dif(self):
      # open the metadata file and assign to meta
      with open(self.genomes_2_metadata, 'r') as meta:
         meta = pd.read_csv(meta)
         meta_accessions = list(meta['accession'])
      #print(meta_accessions)
      #print(type(meta_accessions[0]))
      busco_dir_list = os.listdir(self.busco_dir)
      for name in busco_dir_list:
         #print(name)
         if 'GC' in name:
            check_this = name[:15]
            #print(type(check_this))
            if check_this not in meta_accessions:
               print(name)
               pass
checker = checker_b()
checker.dif()
