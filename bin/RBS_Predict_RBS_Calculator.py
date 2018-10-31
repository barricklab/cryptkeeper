#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

@authors: Jeffrey Barrick

Uses RBS calculator version 1.0.

Setup:

1) Running the RBA calculator requires the ViennaRNA and NuPACK software program.

Installation of NuPACK 3.2.2 on Mac OSX was successful but mfe function hangs
indefinitely, so there appears to be a problem with the executable. Installation
of NuPACK 3.2.2 on Linux systems works properly.

2) You must add a shebang to the top of the Run_RBS_Calculator.py script:

#!/usr/bin/env python

and put it and all RBS Calculator scripts in your $PATH so that it can be invoked 
from the command line via just that name as a command.

"""
import argparse
import Bio
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import csv
import subprocess
from operator import itemgetter


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
    dest='o',
    required=True,
    type=str,
    help="prefix/title of .txt outfile")
##------------------------------------------------------------------------------
options = parser.parse_args()


i=0
sequence_length = 0
forward_seq = ''
reverse_seq = ''

for this_seq in SeqIO.parse(options.i, "fasta"):
  i += 1
  if (i>1):
    exit()
  sequence_length = len(this_seq)
  forward_seq = str(this_seq.seq)
  reverse_seq = str(this_seq.reverse_complement().seq)


#run RBS Calculator twice. Once for each strand.
subprocess.call('Run_RBS_Calculator.py ' + forward_seq + ' > ' + options.o + '.forward.predictions.txt', shell = True)
subprocess.call('Run_RBS_Calculator.py '+ reverse_seq + ' > ' + options.o + '.reverse.predictions.txt', shell = True)

#returns a list of dictionaries for the rows
def process_RBS_calculator_output_file(input_file_name, is_reverse_complement, sequence_length):
  entry_list = []
  
  lines = []
  with open(input_file_name) as f:
    lines = f.readlines()
      
  for line in lines:
    line = line.strip()
    if (line == ""):
      continue
    split_line = line.split()
    if len(split_line) != 3:
      continue
    
    new_entry = {}
    
    # RBS calculator returns 0-indexed positions
    pos_1 = int(split_line[0]) + 1
    new_entry["position"] =  pos_1 if not is_reverse_complement else sequence_length - pos_1 + 1;
    new_entry["start_codon"] =  forward_seq[pos_1-1:pos_1+1] if not is_reverse_complement else reverse_seq[pos_1-1:pos_1+1];
    new_entry["strand"] = '-' if is_reverse_complement else '+'
    new_entry["score"] = split_line[1]
    new_entry["score2"] = split_line[2]
    entry_list.append(new_entry)
      
  return(entry_list)
  
# Parse output and create one summary file
forward_list = process_RBS_calculator_output_file(options.o +'.forward.predictions.txt', False, sequence_length)
reverse_list = process_RBS_calculator_output_file(options.o +'.reverse.predictions.txt', True, sequence_length)
final_list = forward_list
final_list.extend(reverse_list)

# Sort list by coordinate
final_list = sorted(final_list, key=itemgetter('position')) 

with open(options.o,'w') as final_predictions_file:
  writer = csv.DictWriter(
      final_predictions_file, 
      fieldnames = ["position", "strand", "start_codon", "score", "score2"]
    )
  writer.writeheader()
  writer.writerows(final_list)
final_predictions_file.close()
