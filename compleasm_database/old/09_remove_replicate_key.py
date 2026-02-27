#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 12:25:52 2024

@author: adam

this code is supposed to remove the repeated sequence from the renamed fasta file
it does this by making a dictionary from the fasta file, so any header that appears twice gets replaced
"""
import os

class editor_b():
   def __init__(self):
      self.top = '/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/genomes_2/busco_run_05_results/01_pre_alignments'
      self.consolidate = '/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/genomes_2/busco_run_05_results/02_consolidate'
   def getfasta(self, fastahandle):
       fastadata = {}
       name = None
       sequence = ''
       for line in fastahandle:
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
       return(fastadata)
   def remove_duplicates(self):
      os.getcwd()
      os.chdir(self.top)
      # os.walk down the top directory
      for (root, dirs, file) in os.walk(self.top):
          # for name in the list of files
          for name in file:
              # if the name contains renamed
              if 'renamed_' in name:
                  # change the directory to where the file is contained
                  os.chdir(root)
                  # open the file
                  with open(name, 'r') as fasta:
                      # read it using the get fasta script
                      fasta_dictionary = editor_b.getfasta(self, fasta)
                  # create newfile name
                  text = name
                  delimiter = '_'
                  index = text.find(delimiter)
                  if index != -1:
                      newfile = name[index+1:]
                  os.chdir(self.consolidate)
                  with open(newfile, 'w') as newfile:
                      for key in fasta_dictionary:
                        newfile.write('>'+ key + '\n' + fasta_dictionary[key] + '\n')

editor = editor_b()
editor.remove_duplicates()

