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
from operator import itemgetter
import re


#------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='input fasta')
parser.add_argument(
    '-i',  
    action='store', 
    dest='i',
    required=True,
    type=str,
    help="input fasta file",
    )

parser.add_argument(
    '-o',  
    action='store', 
    dest='o',
    required=True,
    type=str,
    help="output file prefix",
    )
    

#------------------------------------------------------------------------------

#Generic function for splitting input sequence file
def create_split_sequence_file(input_file_name, output_file_name, window_size = 10000, window_offset = 10000, sequence_format="fasta", name_delimiter="___"):

  split_fasta_list = []

  for this_seq in SeqIO.parse(input_file_name, sequence_format):
    for n in range(0, max(1,int((len(this_seq)/window_offset)))):
      start_0_indexed = 0 + n*window_offset
      end_0_indexed = min(start_0_indexed + window_size - 1, len(this_seq))
    
      #print(str(start_1) + "-" + str(end_1))
    
      split_fasta_name = name_delimiter.join([this_seq.id, str(start_0_indexed+1), str(end_0_indexed+1)])
      split_fasta_seq  = this_seq.seq[start_0_indexed:end_0_indexed+1]

      split_fasta_record = SeqRecord(seq = split_fasta_seq, id = split_fasta_name, description="")
      split_fasta_list.append(split_fasta_record)

  #write split sequence file
  SeqIO.write(split_fasta_list, output_file_name, sequence_format)


#------------------------------------------------------------------------------

options = parser.parse_args()

## global setting - separates names/coords
name_delimiter = "___"


# split into smaller FASTA for processing
create_split_sequence_file(
  input_file_name = options.i, 
  output_file_name = options.o + ".split.fa",
  window_size = options.ws,
  window_offset = options.wo,
  sequence_format = "fasta",
  name_delimiter = name_delimiter,
  )

i=0
for this_seq in SeqIO.parse(options.i, "fasta"):
  i += 1
  if (i>1):
    exit()
  SeqIO.write(this_seq, options.o + '.forward.fa', "fasta")
  SeqIO.write(this_seq.reverse_complement(), options.o + '.reverse.fa', "fasta")


#run BPROM twice. Once for each strand. 
subprocess.call('bprom '+ options.o + '.forward.fa -o '+ options.o +'forward.predictions.txt', shell = True)
subprocess.call('bprom '+ options.o + '.reverse.fa -o '+ options.o +'reverse.predictions.txt', shell = True)

## Parse output and create one summary file


  # gff_line_list = []
#   
#   with open(options.o +'forward.predictions.txt') as f:
#   ()
# 
# 
#   existing_predictions = {}
#   
#   for line in gff_line_list:
#   
#     if re.match('\s?#', line):
#       continue
#     
#     split_line = line.strip().split('\t')
#     
#     start_offset = int(split_line[0].split(name_delimiter)[1])-1
#     
#     # Make a dictionary for all the extra fields
#     extra_field_list = split_line[8].split(';')
#     extra_field_dict = {}
#     for extra_field_item in extra_field_list:
#       extra_field_item_list = extra_field_item.split('=')
#       extra_field_dict[extra_field_item_list[0]] = extra_field_item_list[1]
#   
#     new_entry = {
#       "promoter" : extra_field_dict["promoter"],
#       "strand" : split_line[6],
#       "TSSpos" : start_offset + int(split_line[4]),
#       "score" : split_line[5],
#       "box10pos" : start_offset + int(extra_field_dict["box10pos"]),
#       "box35pos" : start_offset + int(extra_field_dict["box35pos"]),
#       "box10seq" : extra_field_dict["box10seq"],
#       "box35seq" : extra_field_dict["box35seq"],
#       }
#     
#   #sort list by coordinate
#   sorted_entry_list = sorted(entry_list, key=itemgetter('TSSpos')) 
# 
#   with open(options.o + '.final_predictions.tsv','w') as final_predictions_file:
#     writer = csv.DictWriter(
#         final_predictions_file, 
#         delimiter = '\t',
#         fieldnames = ["promoter", "score", "strand", "TSSpos", "box35pos", "box35seq", "box10pos", "box10seq",]
#       )
#     writer.writeheader()
#     writer.writerows(sorted_entry_list)
#   final_predictions_file.close()
# 
