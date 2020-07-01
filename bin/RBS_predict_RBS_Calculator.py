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

@TODO: Make it so that RBS Calculator outputs something sane.

"""
import argparse
import Bio
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import csv
import subprocess
from operator import itemgetter
import os
import json
from helpers import timer_decorator
import concurrent.futures

#returns a list of dictionaries for the rows
def process_RBS_calculator_output_file(input_file_name, is_reverse_complement,
                                       sequence_length, forward_seq, reverse_seq):
  entry_list = []
  if not os.path.exists(input_file_name):
    return entry_list
  lines = []
  with open(input_file_name) as f:
    lines = f.readlines()

  failed_prediction = False
  line_index = 0
  while line_index < len(lines):
    line = lines[line_index].strip()

    # First line of pair
    #   1 or 'Warning: There is a leaderless start codon, which is being ignored.'
    if line_index % 2 == 0:

      failed_prediction = False
      #append em
      if line == "'Warning: There is a leaderless start codon, which is being ignored.'":
        failed_prediction = True

    # Second line of pair
    #  pos_1, score, score2
    else:

      split_line = []
      if not failed_prediction:
        split_line = line.split()
        if len(split_line) != 3:
          print("Unexpected number of items in line: " + line)
          continue
      else: #failed prediction: put in dummy score entry so we still sync with ORF predictions
        split_line = [0, 0.00000001, 0]

      new_entry = {}

      # RBS calculator returns 0-indexed positions
      pos_1 = int(split_line[0]) + 1
      new_entry["position"] =  pos_1 if not is_reverse_complement else sequence_length - pos_1 + 1;
      new_entry["start_codon"] =  forward_seq[pos_1-1:pos_1+2] if not is_reverse_complement else reverse_seq[pos_1-1:pos_1+2];
      new_entry["strand"] = '-' if is_reverse_complement else '+'
      new_entry["score"] = split_line[1]
      new_entry["score2"] = split_line[2]
      entry_list.append(new_entry)

    line_index = line_index+1

  return(entry_list)

@timer_decorator
def main(options):
  script_path = os.path.dirname(os.path.realpath(__file__))
  start_codons = ['ATG', 'GTG', 'TTG']
  #start_codons = ['ATG', 'GTG']


  i=0
  sequence_length = 0
  forward_seq = ''
  reverse_seq = ''

  for this_seq in SeqIO.parse(options.i, "fasta"):
    i += 1
    if (i>1):
      exit()
    sequence_length = len(this_seq)
    forward_seq = str(this_seq.upper().seq)
    reverse_seq = str(this_seq.reverse_complement().upper().seq)


  #Three modes possible
  if options.s:
    #This one only tests the requested stop codons (quicker)
    if os.path.exists(options.o + '.forward.predictions.txt'):
      os.remove(options.o + '.forward.predictions.txt')
    if os.path.exists(options.o + '.reverse.predictions.txt'):
      os.remove(options.o + '.reverse.predictions.txt')

    orf_predictions = []
    orf_reader = csv.DictReader(open(options.s))

    #'''  Updated code for parallel processing
    def _rbs_requested_predictions(orf):
      orf["start"] = int(orf["start"])
      if (orf["strand"] == '-'):
        orf["start"] = len(forward_seq) - int(orf["end"]) + 1
      # Note: RBS Calculator expects 0-indexed positions
      sequence = (forward_seq if (orf["strand"] == '+') else reverse_seq)
      start_loc = str(orf["start"] - 1)
      outpath = options.o + ('.forward.predictions.txt' if (orf["strand"] == '+') else '.reverse.predictions.txt')
      subprocess.call('Run_RBS_Calculator.py ' + sequence + ' ' + start_loc + ' >> ' + outpath, shell=True)

    with concurrent.futures.ThreadPoolExecutor() as multiprocessor:
      result = multiprocessor.map(_rbs_requested_predictions, orf_reader)

    '''
    for orf in orf_reader:
      #print(orf)
      orf["start"] = int(orf["start"])
      if (orf["strand"]=='-'):
        orf["start"] = len(forward_seq) - int(orf["end"]) + 1
      #Note: RBS Calculator expects 0-indexed positions
      sequence = (forward_seq if (orf["strand"]=='+') else reverse_seq)
      start_loc = str(orf["start"]-1)
      outpath = options.o + ('.forward.predictions.txt' if (orf["strand"]=='+') else '.reverse.predictions.txt')
      subprocess.call('Run_RBS_Calculator.py ' + sequence + ' ' + start_loc + ' >> ' + outpath, shell = True)
    '''

  elif options.e != None:
    #This one pre-screens for start codons with reasonable RBS sequences upstream

    #Load binding energies
    energies={}
    with open(os.path.join(script_path, "../data/RBS-binding-energies-ACCUCCUUA.json")) as json_file:
      energies = json.load(json_file)

    #searches all 9mers in the window from 1-20 bp upstream
    rbs_max_from_start = 25
    rbs_min_from_start = 0
    rbs_length = 9
    orf_list = []

    for start_pos_0 in range(len(this_seq)):
      if forward_seq[start_pos_0:start_pos_0+3] in start_codons:

        best_energy = 100
        first_rbs_pos_0 = max(0, start_pos_0-rbs_max_from_start)
        last_rbs_pos_0 = start_pos_0-rbs_min_from_start-rbs_length
        #print("Testing range: " + str(first_rbs_pos_0) + "-" + str(last_rbs_pos_0))

        for start_rbs_pos_0 in range(first_rbs_pos_0, last_rbs_pos_0):
          this_energy = energies[forward_seq[start_rbs_pos_0:(start_rbs_pos_0+rbs_length)]]
          best_energy = min(best_energy, this_energy)

        if best_energy <= options.e:
          new_orf = {}
          new_orf["start_codon"] = forward_seq[start_pos_0:start_pos_0+3]
          new_orf["strand"] = "+"
          new_orf["start"] = start_pos_0+1
          new_orf["end"] = 0 #dummy value
          new_orf["energy"] = best_energy
          orf_list.append(new_orf)

    for start_pos_0 in range(len(this_seq)):
      if reverse_seq[start_pos_0:start_pos_0+3] in start_codons:

        best_energy = 100
        first_rbs_pos_0 = max(0, start_pos_0-rbs_max_from_start)
        last_rbs_pos_0 = start_pos_0-rbs_min_from_start-rbs_length
        #print("Testing range: " + str(first_rbs_pos_0) + "-" + str(last_rbs_pos_0))

        for start_rbs_pos_0 in range(first_rbs_pos_0, last_rbs_pos_0):
          this_energy = energies[reverse_seq[start_rbs_pos_0:(start_rbs_pos_0+rbs_length)]]
          best_energy = min(best_energy, this_energy)

        if best_energy <= options.e:
          new_orf = {}
          new_orf["start_codon"] = reverse_seq[start_pos_0:start_pos_0+3]
          new_orf["strand"] = "-"
          new_orf["start"] = len(this_seq) - start_pos_0
          new_orf["end"] = 0 #dummy value
          new_orf["energy"] = best_energy
          orf_list.append(new_orf)

    for orf in orf_list:
      print(orf)
      #print(orf)
      orf["start"] = int(orf["start"])
      if (orf["strand"]=='-'):
        orf["start"] = len(forward_seq) - int(orf["end"]) + 1
      #Note: RBS Calculator expects 0-indexed positions
      from Run_RBS_Calculator import Run_RBS_Calculator
      seq = forward_seq if orf["strand"] == '+' else reverse_seq
      start_pos = str(orf["start"]-1)
      #Run_RBS_Calculator(seq, start_pos)
      subprocess.call('Run_RBS_Calculator.py ' + (forward_seq if (orf["strand"]=='+') else reverse_seq) + ' ' + str(orf["start"]-1) + ' >> ' + options.o + ('.forward.predictions.txt' if (orf["strand"]=='+') else '.reverse.predictions.txt'), shell = True)



  else:
    #Run RBS Calculator twice on entire sequences. Once for each strand.
    subprocess.call('Run_RBS_Calculator.py ' + forward_seq + ' > ' + options.o + '.forward.predictions.txt', shell = True)
    subprocess.call('Run_RBS_Calculator.py ' + reverse_seq + ' >> ' + options.o + '.reverse.predictions.txt', shell = True)


  # Parse output and create one summary file
  forward_list = process_RBS_calculator_output_file(options.o +'.forward.predictions.txt', False, sequence_length,
                                                    forward_seq, reverse_seq)
  reverse_list = process_RBS_calculator_output_file(options.o +'.reverse.predictions.txt', True, sequence_length,
                                                    forward_seq, reverse_seq)
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

def RBS_predict_RBS_Calculator(input, output, starts=None, minenergy=None):
  # Pretend to be an argument parser
  class ObjectClass:
    pass
  options = ObjectClass
  options.i = input
  options.o = output
  if starts:
    options.s = starts
  if minenergy:
    options.e = minenergy
  main(options)

if __name__ == "__main__":
  # ------------------------------------------------------------------------------
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

  parser.add_argument('-s',
                      action='store',
                      dest='s',
                      type=str,
                      help="csv file containing starts to test (generated by ORF_finder.py)")

  parser.add_argument('-e',
                      action='store',
                      dest='e',
                      type=float,
                      default=None,
                      help="minimum energy for RBS binding")
  ##------------------------------------------------------------------------------
  options = parser.parse_args()
  main(options)
