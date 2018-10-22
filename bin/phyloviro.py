#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 11:11:29 2018

@authors: Jeffrey Barrick
Purpose: Import a file of Genbank files downloaded for a virus family into a set of files aligned
         by ORF for each RNA
"""
import argparse
import Bio
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import csv
import os

#------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='')
parser.add_argument('-i',  
    action='store', 
    dest='i',
    required=False,
    default=".",
    type=str,
    help="input path to Genbank files, which  must have file extensions: gbff, gbk, or gb")

parser.add_argument('-o',  
    action='store', 
    dest='o',
    required=False,
    default="output",
    type=str,
    help="prefix for all output files")
##------------------------------------------------------------------------------
options = parser.parse_args()

##Some lengths and tolerances for the different RNAs
## These are for CMV, should make this a command-line option
RNA_descriptions = [
  {"avg" : 3361, "tol": 100, "desc": "RNA1"},
  {"avg" : 3047, "tol": 100, "desc": "RNA2"},
  {"avg" : 2175, "tol": 100, "desc": "RNA3"}
]

directory = os.fsencode(options.i)
genbank_file_name_list = []

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith( (".gbff", ".gbk", ".gb") ): 
        genbank_file_name_list.append(filename)
        

        continue
    else:
        continue


#List of all records we are writing per RNA
output_sequence_lists = []
for RNA_desc in RNA_descriptions:
  output_sequence_lists.append([])

for genbank_file_name in genbank_file_name_list:
  genbank_file_path = os.path.join(options.i, genbank_file_name)
  print(genbank_file_path)
  genbank_records = list(SeqIO.parse(genbank_file_path, "genbank"))
  
  num_accepted = 0
  accepted_list = []
  for RNA_desc in RNA_descriptions:
    accepted_list.append(0)
  
  for genbank_record in genbank_records:
    print("  " + genbank_record.id + "|" + genbank_record.name + "|" + genbank_record.description + "|" + str(len(genbank_record)))
    
    
    assigned = -1
    
    #Check all RNA descriptions
    
    if assigned == -1:
      for i in range(0, len(RNA_descriptions)):
        if genbank_record.description.find(RNA_descriptions[i]["desc"]) != -1:
          assigned = i
          break
          
    if assigned == -1:
      for i in range(0, len(RNA_descriptions)):
        if genbank_record.description.find(" " + str(i+1)) != -1:
          assigned = i
          break

    if assigned == -1:
      for i in range(0, len(RNA_descriptions)):
        RNA_desc = RNA_descriptions[i]
        if (len(genbank_record) >= RNA_desc["avg"]-RNA_desc["tol"]) and (len(genbank_record) <= RNA_desc["avg"]+RNA_desc["tol"]):
          assigned = i
          break
    
    if assigned == -1:
      text = input("Could not determine RNA number, please enter selection:")
      assigned = int(text)-1
    
    # Don't pay attention to partial
    #If genbank_record.description.find("partial") != -1:
    #  assigned = -1
      
    if assigned > -1:
      print("Accepted: " + RNA_descriptions[assigned]["desc"])
      output_sequence_lists[assigned].append(genbank_record)
      accepted_list[assigned] += 1
      num_accepted += 1
    else:
      print("Rejected")
  
  if num_accepted != 0 and num_accepted < len(RNA_descriptions):
    print("Only accepted some RNAs: check settings")
    
  for i in range(0, len(RNA_descriptions)):
    if accepted_list[i] != 1:
      print("WARNING: Some RNAs found multiple times or not at all for input file")

for i in range(0, len(RNA_descriptions)):
  output_file_name = options.o + RNA_descriptions[i]["desc"] + ".fa"
  SeqIO.write(output_sequence_lists[i], output_file_name, "fasta")

  
  ## Classify into the different RNA families
  
# for filename in os.listdir(directory):
# 
# options.i
# 
# fasta = SeqIO.read(options.i,"genbank ",generic_dna)#import fasta, can make into argparse input later
# 
# with open("../data/json-energyRef-ACCUCCUUA.txt","r") as energies: 
#     binding_energies = eval(energies.read()) #evaluates object/data structure of file
# energies.close()
# 
# RNA_fasta = fasta.seq.transcribe()
# uppercase_RNA = RNA_fasta.upper()
# start_site = []
# stop_site = []
# sequence = []
# energy = []
# 
# for i in range(len(uppercase_RNA)-9):
#     query = uppercase_RNA[i:i+9]
#     start_site.append(i+1)
#     stop_site.append(i+9)
#     sequence.append(query)
#     energy.append(binding_energies[query])
#     i += 1
#     
# with open(options.o+'_RBS_energies_F.txt',"w") as outfile:
#     writer = csv.writer(outfile, delimiter = '\t')
#     writer.writerows(zip(start_site, stop_site, sequence, energy))
#     outfile.close()
# 
# #- stranded analysis
# reverse_RNA_fasta = uppercase_RNA.reverse_complement()
# r_start_site= []
# r_stop_site = []
# r_sequence = []
# r_energy = []
# 
# for i in range (len(reverse_RNA_fasta)-9):
#     query = reverse_RNA_fasta[i:i+9]
#     r_start_site.append(len(RNA_fasta)-(i+1))
#     r_stop_site.append(len(RNA_fasta)-(i+9))
#     r_sequence.append(query)
#     r_energy.append(binding_energies[query])
#     
# with open(options.o+'_RBS_energes_R.txt','w') as outfile:
#     writer = csv.writer(outfile,delimiter = '\t')
#     writer.writerows(zip(r_start_site,r_stop_site,r_sequence,r_energy))
#     outfile.close()
