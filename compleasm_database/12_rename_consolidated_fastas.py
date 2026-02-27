#!/usr/bin/env python
import csv
import os

class _rename:
   def __init__(self):
      # this is where the consolidate folder is
      self.dir_consolidate = '/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/genomes_2/busco_run_05_results/02_consolidate'
      self.meta = '/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/genomes_2/metadata_genomes_2.csv'
      #######################
      ####################### This directory needs to be made before running this script, it is a subdirectory if self.dir_conolidate
      #######################
      self.dir_names = '/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/genomes_2/busco_run_05_results/02_consolidate/02_with_names'
      self.accession_to_organism = {}
      self.fasta_list = []

   # Load metadata from CSV file
   def load_metadata(self):
      with open(self.meta, 'r') as file:
         reader = csv.DictReader(file)
         for row in reader:
            accession = row['accession']
            organism_name = row['organismName']
            self.accession_to_organism[accession] = organism_name
      return self.accession_to_organism

   # Load list of fastas which need to have their headers fixed
   def load_fasta_list(self):
      os.chdir(self.dir_consolidate)
      file_list = os.listdir(self.dir_consolidate)
#      print(file_list)
      for file in file_list:
         if file.endswith('fasta'):
            self.fasta_list.append(file)
      #print(self.fasta_list)
   # this function can be made better by including the following
   '''
   text =
   delimiter_1 = GC
   delimiter_2 = .
   index_1 = text.find(delimiter_1)
   index_2 = text.find(delimiter_2)
   if not index = -1:
      accession = text[index_1:index_2+2]
   '''
   def getfasta_replace_headers(self, fastahandle):
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
                   accession = line[1:len(line)]
                   for accession_q in self.accession_to_organism:
                      if accession == accession_q:
                         name = self.accession_to_organism[accession_q]
               else:
                   # report the previous name and sequence, then refresh the name and sequence
                   fastadata[name] = sequence
                   line = line.rstrip()
                   accession = line[1:len(line)]
                   for accession_q in self.accession_to_organism:
                      if accession == accession_q:
                         name = self.accession_to_organism[accession_q]
                   sequence = ''
                   # restart the sequence
       # report the last name and sequence
       fastadata[name] = sequence
       return(fastadata)

   # Replace accession numbers with organism names in a single fasta file
   def write_new_fasta(self, filename):
      os.chdir(self.dir_consolidate)
      with open(filename, 'r') as fasta:
         fasta = self.getfasta_replace_headers(fasta)
      os.chdir(self.dir_names)
      with open(filename, 'w') as new_fasta:
         for header, seq in fasta.items():
            new_fasta.write('>' + header + '\n')
            new_fasta.write(seq + '\n')
   # this iterates through the fasta list and replaces the headers for each file and writes the new file to a different directory
   def loop(self):
      for file in self.fasta_list:
         print('going to rename',file)
         self.write_new_fasta(file)


rename = _rename()
rename.load_metadata()
rename.load_fasta_list()
rename.loop()

#       with open(output_file, 'w') as out_file:
           # Write the first line unchanged (header of the phylip file)
#           out_file.write(lines[0])

           # Replace taxa names (accession numbers) with organism names in the sequences
#           for line in lines[1:]:
#               parts = line.split()
#               accession = parts[0]
#               sequence = parts[1]
#               organism_name = accession_to_organism.get(accession, accession)  # Default to accession if no match
#               out_file.write(f"{organism_name:<10} {sequence}\n")

   # Replace accession numbers with organism names in the PHYLIP file
#   def replace_taxa_in_phylip(phylip_file, accession_to_organism, output_file):
#       with open(phylip_file, 'r') as file:
#           lines = file.readlines()
#
#       with open(output_file, 'w') as out_file:
#           # Write the first line unchanged (header of the phylip file)
#           out_file.write(lines[0])
#
#           # Replace taxa names (accession numbers) with organism names in the sequences
#           for line in lines[1:]:
#               parts = line.split()
#               accession = parts[0]
#               sequence = parts[1]
#               organism_name = accession_to_organism.get(accession, accession)  # Default to accession if no match
#               out_file.write(f"{organism_name:<10} {sequence}\n")

#if __name__ == "__main__":
#    csv_file = '/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/genomes_2/metadata_genomes_2.csv'
#    #phylip_file = 'concatenated_names.phylip'
#    #output_file = 'updated_concatenated_names.phylip'
#    # Load metadata from CSV
#    accession_to_organism = load_metadata(csv_file)

    # Replace taxa names in the PHYLIP file
    #replace_taxa_in_phylip(phylip_file, accession_to_organism, output_file)

