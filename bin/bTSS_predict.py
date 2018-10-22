#!/usr/bin/env python

import argparse
import Bio
from Bio.SeqIO import FastaIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess
import csv

#------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='import fasta for bTSS slicing')
parser.add_argument('-i',  
    action='store', 
    dest='i',
    required=True,
    type=str,
    help="input .fasta file to slice upstream of bTSSfinder")

parser.add_argument('-o',  
    action='store', 
    dest='o',
    required=True,
    type=str,
    help="output file prefix")

#------------------------------------------------------------------------------
options = parser.parse_args()

## separates names/coords
name_delimiter = "___"

split_fasta_list = []

#prepare fasta sections of 500 bp for bTSSfinder 
for this_seq in SeqIO.parse(options.i, "fasta"):
  for n in range(0, max(1,int((len(this_seq)/250)))):
    start_1 = 0 + n*250
    end_1 = min(start_1 + 500 - 1, len(this_seq))
    
    #print(str(start_1) + "-" + str(end_1))
    
    split_fasta_name = name_delimiter.join([this_seq.id, str(start_1), str(end_1)])
    split_fasta_seq  = this_seq.seq[start_1:end_1]

    split_fasta_record = SeqRecord(seq = split_fasta_seq, id = split_fasta_name, description="")
    split_fasta_list.append(split_fasta_record)

#write fasta for bTSSfinder
SeqIO.write(split_fasta_list, options.o + ".split.fa", "fasta")

#run bTSSfinder 
subprocess.call('bTSSfinder -i '+ options.o + '.split.fa -o '+ options.o +' -h 2', shell = True)
      
## Read in entries from the .bed file
bedfile_summary = open(options.o + ".bed","r")
lines = bedfile_summary.readlines()
entry_list = []

search_coords = []
sigma_factors = []
start_codons = []
sigma_site = []
direction = []
for line in lines:
  data = line.strip().split('\t')
  
  start_offset = int(data[0].split(name_delimiter)[1])
    
  this_minus_35_position = 0
  this_minus_10_position = 0
  this_sigma_factor_sites_list = data[11].split(',')
  if len(this_sigma_factor_sites_list) > 0:
    this_minus_35_position = this_sigma_factor_sites_list[0]
  if len(this_sigma_factor_sites_list) > 1:
    this_minus_10_position = this_sigma_factor_sites_list[1]
  
  new_entry = {
      "sigma_factor" : data[3],
      "strand" : data[5],
      "TSS_position" : start_offset + int(data[6]),
      "minus_10_position" : start_offset + int(this_minus_10_position),
      "minus_35_position" : start_offset + int(this_minus_35_position),
    }
  entry_list.append(new_entry)

bedfile_summary.close()

with open(options.o + '_corrected.txt','w') as corrected_bed:
  writer = csv.DictWriter(
      corrected_bed, 
      delimiter = '\t',
      fieldnames = ["sigma_factor", "strand", "TSS_position", "minus_35_position", "minus_10_position"]
    )
  writer.writeheader()
  writer.writerows(entry_list)
corrected_bed.close()

    
    

