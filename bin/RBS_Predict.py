
# coding: utf-8

# In[37]:


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 11:11:29 2018

@authors: Alexandra Lukasiewicz and Matt McGuffie
Purpose: Using A. Hockenberry SD sequence energy prediction, scan insert fasta 
for potential Ribosome Binding Sites (RBS)
"""
import Bio
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import csv

#------------------------------------------------------------------------------

fasta = SeqIO.read("CMV_RNA1.fa","fasta",generic_dna)#import fasta, can make into argparse input later

with open("json-energyRef-ACCUCCUUA.txt","r") as energies: 
    binding_energies = eval(energies.read()) #evaluates object/data structure of file
energies.close()

RNA_fasta = fasta.seq.transcribe()
uppercase_RNA = RNA_fasta.upper()
start_site = []
stop_site = []
sequence = []
energy = []

for i in range(len(uppercase_RNA)-9):
    query = uppercase_RNA[i:i+9]
    binding_energies[query]
    start_site.append(i)
    stop_site.append(i+9)
    sequence.append(query)
    energy.append(binding_energies[query])
    i += 1
    
with open("RBS_energies.txt","w") as outfile:
    writer = csv.writer(outfile, delimiter = '\t')
    writer.writerows(zip(start_site, stop_site, sequence, energy))
    outfile.close()
    

