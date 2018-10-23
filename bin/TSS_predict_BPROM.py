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


# Create FASTA for two strands
# (Note that this also wraps the lines of the FASTA if needed, which is required for BROM)

i=0
sequence_length = 0
for this_seq in SeqIO.parse(options.i, "fasta"):
  i += 1
  if (i>1):
    exit()
  sequence_length = len(this_seq)
  SeqIO.write(this_seq, options.o + '.forward.fa', "fasta")
  SeqIO.write(this_seq.reverse_complement(), options.o + '.reverse.fa', "fasta")


#run BPROM twice. Once for each strand. 
#subprocess.call('bprom '+ options.o + '.forward.fa -o '+ options.o +'forward.predictions.txt', shell = True)
#subprocess.call('bprom '+ options.o + '.reverse.fa -o '+ options.o +'reverse.predictions.txt', shell = True)


#returns a list of dictionaries for the rows
def process_bprom_output_file(input_file_name, is_reverse_complement, sequence_length):
  entry_list = []
  
  lines = []
  with open(input_file_name) as f:
    lines = f.readlines()
    
  #remove top four lines
  lines = lines[4:len(lines)+1]
  
  i=0
  new_entry = {}
  for line in lines:
    line = line.strip()
    if (line == ""):
      break
    split_line = line.split()
    i = i+1
    if (i==1):
      new_entry["promoter"] = 'sigma70';
      new_entry["strand"] = '-' if is_reverse_complement else '+'
      new_entry["TSSpos"] = int(split_line[2]) if not is_reverse_complement else sequence_length - int(split_line[2]) + 1
      new_entry["score"] = split_line[4]
    elif (i==2):
      new_entry["box10pos"] = int(split_line[4]) if not is_reverse_complement else sequence_length - int(split_line[4]) + 1
      new_entry["box10seq"] = split_line[5]
    elif (i==3):
      new_entry["box35pos"] = int(split_line[4]) if not is_reverse_complement else sequence_length - int(split_line[4]) + 1
      new_entry["box35seq"] = split_line[5]
    
      entry_list.append(new_entry)
      new_entry = {}
      i=0
  return(entry_list)
  
# Parse output and create one summary file
forward_list = process_bprom_output_file(options.o +'.forward.predictions.txt', False, sequence_length)
reverse_list = process_bprom_output_file(options.o +'.reverse.predictions.txt', True, sequence_length)
final_list = forward_list
final_list.extend(reverse_list)

#sort list by coordinate
final_list = sorted(final_list, key=itemgetter('TSSpos')) 

with open(options.o + '.final_predictions.tsv','w') as final_predictions_file:
  writer = csv.DictWriter(
      final_predictions_file, 
      delimiter = '\t',
      fieldnames = ["promoter", "score", "strand", "TSSpos", "box35pos", "box35seq", "box10pos", "box10seq",]
    )
  writer.writeheader()
  writer.writerows(final_list)
final_predictions_file.close()
