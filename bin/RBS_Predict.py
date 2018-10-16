#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 11:11:29 2018

@authors: Alexandra Lukasiewicz and Matt McGuffie
Purpose: Using A. Hockenberry SD sequence energy prediction, scan insert fasta 
for potential Ribosome Binding Sites (RBS)
"""
import argparse
import Bio
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import csv

#------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='import fasta file for RBS calculation')
parser.add_argument('-i',  
    action='store', 
    dest='i',
    required=True,
    type=str,
    help="input file in fasta format")

parser.add_argument('-o',  
    action='store', 
    dest='i',
    required=True,
    type=str,
    help="prefix/title of .txt outfile")
##------------------------------------------------------------------------------
options = parser.parse_args()

fasta = SeqIO.read(options.i,"fasta",generic_dna)#import fasta, can make into argparse input later

with open("../data/json-energyRef-ACCUCCUUA.txt","r") as energies: 
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
    start_site.append(i+1)
    stop_site.append(i+9)
    sequence.append(query)
    energy.append(binding_energies[query])
    i += 1
    
with open(options.o+'_RBS_energies_F.txt',"w") as outfile:
    writer = csv.writer(outfile, delimiter = '\t')
    writer.writerows(zip(start_site, stop_site, sequence, energy))
    outfile.close()

#- stranded analysis
reverse_RNA_fasta = uppercase_RNA.reverse_complement()
r_start_site= []
r_stop_site = []
r_sequence = []
r_energy = []

for i in range (len(reverse_RNA_fasta)-9):
    query = reverse_RNA_fasta[i:i+9]
    r_start_site.append(len(RNA_fasta)-(i+1))
    r_stop_site.append(len(RNA_fasta)-(i+9))
    r_sequence.append(query)
    r_energy.append(binding_energies[query])
    
with open(options.o+'_RBS_energes_R.txt','w') as outfile:
    writer = csv.writer(outfile,delimiter = '\t')
    writer.writerows(zip(r_start_site,r_stop_site,r_sequence,r_energy))
    outfile.close()