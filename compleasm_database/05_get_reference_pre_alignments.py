#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 13:42:53 2024

@author: adam
"""
####
# When running on lepidodactylus update:
#     self.busco
#     self.pre_dir
#     self.orthologs
#     self.ref_dir
####

import os
import os.path as p
import pandas as pd
import csv

class ref_pre_align_b:
   def __init__(self):
      # this is the directory of the pipeline - at the moment I can't remeber why i need the pipeline directory
      self.busco='/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/genomes_2/busco_run_05_results'
      # this is the path where pre_alignments will be stored?
      self.pre_dir= '/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/genomes_2/busco_run_05_results/01_pre_alignments/'
      # this needs to be the path of the output from 01_setup.py
      with open('/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/genomes_2/busco_run_05_results/shared_complete.csv', 'r') as orthologs:
         self.orthologs = pd.read_csv(orthologs)
      # this should be the busco results folder for the reference squamate
      self.ref_dir='/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/genomes_2/busco_run_05_results/reference/GCF_004329235.1_PodMur_1.0_genomic.fna/run_sauropsida_odb10/busco_sequences/single_copy_busco_sequences'
   def get_sequence(self, filename):
      fastadata = {}
      name = None
      sequence = ''
      for line in filename:
          if not line.startswith(">"):
              # this is continuing a sequence
              line = line.rstrip()
              sequence = sequence + ''.join(line)

          elif line.startswith(">"):
              # start of a new sequence
              if name is None:
              # the first sequence in the file!
                  line = line.rstrip()
                  name = line[1:len(line)]
              else:
                  # report the previous name and sequence, then refresh the name and sequence
                  fastadata[name] = sequence
                  line = line.rstrip()
                  name = line[1:len(line)]
                  sequence = ''
                  # restart the sequence
      # report the last name and sequence
      fastadata[name] = sequence
      return(sequence)

   '''
   I am now going to use the complete busco ids list to create a single file for each busco ID.
   Within each file there will be sequences with fasta headers that will contain the accession
    of each species from which the sequence came from.
   '''
   def dir_for_each_busco(self):
      # change directory to top
      os.chdir(self.ref_dir)
      # get the current directory
      print(os.getcwd(),'\n', 'I will start writing reference fasta files')

      # for each complete ID in busco IDs
      for ID in self.orthologs:
          # make a dictionary to hold sequence information in the order of:
              # busco_id>GCA#>sequence
          dictionary = {}
          # iterate down the directory tree
          for (root, dirs, file) in os.walk(self.ref_dir, topdown=True):
              # for each item name
              for name in file:
                  # if the file name ends with .fna
                  if name.endswith('.fna'):
                      # partition the name into two pieces the beginning and the extension, os.path.splitext() does the partitioning
                      name_parts = os.path.splitext(name)
                      # name_parts[0] is everything before the extension, I think this is where I filter based on complete busco ID mask
                      if ID == name_parts[0]:
                          # get the path so that I can:
                          path_and_name = root + '/' + name
                          with open(path_and_name, 'r') as fna:
                              # I need to make sure this the fna file has only one sequence in it.
                              sequence = ref_pre_align_b.get_sequence(self, fna)
                          text = path_and_name
                          delimiter_1 = 'GC'
                          delimiter_2 = '.'
                          # I can use .find() to identify the index where the delimiter is and assign this to a variable called delimiter.
                          index_1 = text.find(delimiter_1)
                          index_2 = text.find(delimiter_2)
                          # use this index to get the accession
                          if index_1 != -1:
                              accession = text[index_1:index_2+2]
                          # add the sequence of the accession number to the dictionary that contains sequences from other genomes
                          dictionary[accession] = sequence
          # The follwoing ## lines were used to limi the orthologs to a specific length, if uncommented the indentation will need to be fixed
          ##count = 0
          ##for key in dictionary:
          ##    if len(dictionary[key]) < 10001:
          ##        count +=1
          ################################ this will need to be changed to number of total genomes when I run it on them #########################################################
          ##if count == 1:
              # make an object called newfile to contain sequences from all genomes
          newfile = 'ref_' + ID + '.fasta'
              # change directory to where I want the concatenated files to be saved
          os.chdir(self.pre_dir)
          os.mkdir(ID)
          os.chdir(ID)
              # open the newfile object in append mode
          with open(newfile, 'a') as newfile:
                  # this for loop writes the dictionary in fasta format
              for key in dictionary:
                  newfile.write('>'+ key + '\n' + dictionary[key] + '\n')
          os.chdir(self.ref_dir)
          ##else:
          ##    pass

   '''
   using os.path.splitext was important because I could not just check that ID was in the filename
   this is because , for example: ID == 13at8457 was found in multiple IDs e.g. 9313at8457.fna
   '''
   # name_parts[0] is everything before the extension
   # if ID == name_parts[0]:

ref_pre_align = ref_pre_align_b()
ref_pre_align.dir_for_each_busco()
