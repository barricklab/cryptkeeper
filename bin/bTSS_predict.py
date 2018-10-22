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
parser = argparse.ArgumentParser(description='import fasta for bTSS slicing')
parser.add_argument(
    '-i',  
    action='store', 
    dest='i',
    required=True,
    type=str,
    help="input .fasta file to slice upstream of bTSSfinder",
    )

parser.add_argument(
    '-o',  
    action='store', 
    dest='o',
    required=True,
    type=str,
    help="output file prefix",
    )
    
parser.add_argument(
    '--window-size',  
    action='store', 
    dest='ws',
    type=int,
    help="size of sequence chunks processed",
    default=500,
    )
    
parser.add_argument(
    '--window-offset',  
    action='store', 
    dest='wo',
    type=int,
    help="offset between starts of sequence chunks",
    default=50,
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

#run bTSSfinder 
subprocess.call('bTSSfinder -i '+ options.o + '.split.fa -o '+ options.o +' -h 2 -x 50', shell = True)
      
## Read in entries from one of the output files
## GFF has more info
mode="GFF"

if (mode == "BED"):

## Note: duplicate removal not implemented for BED

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
        "promoter" : data[3].split('_')[0],
        "strand" : data[5],
        "TSSpos" : start_offset + int(data[6]),
        "box10pos" : start_offset + int(this_minus_10_position),
        "box35pos" : start_offset + int(this_minus_35_position),
      }
    entry_list.append(new_entry)

  bedfile_summary.close()

  #sort list by coordinate
  sorted_entry_list = sorted(entry_list, key=itemgetter('TSSpos')) 

  with open(options.o + '_corrected.txt','w') as corrected_bed:
    writer = csv.DictWriter(
        corrected_bed, 
        delimiter = '\t',
        fieldnames = ["promoter", "strand", "TSSpos", "box35pos", "box10pos"]
      )
    writer.writeheader()
    writer.writerows(sorted_entry_list)
  corrected_bed.close()

elif (mode=="GFF"):
  gff_line_list = []
  with open(options.o + ".gff") as f:
    gff_line_list = f.readlines()
  f.close()
  entry_list = []


  existing_predictions = {}
  
  for line in gff_line_list:
  
    if re.match('\s?#', line):
      continue
    
    split_line = line.strip().split('\t')
    
    start_offset = int(split_line[0].split(name_delimiter)[1])-1
    
    # Make a dictionary for all the extra fields
    extra_field_list = split_line[8].split(';')
    extra_field_dict = {}
    for extra_field_item in extra_field_list:
      extra_field_item_list = extra_field_item.split('=')
      extra_field_dict[extra_field_item_list[0]] = extra_field_item_list[1]
  
    new_entry = {
      "promoter" : extra_field_dict["promoter"],
      "strand" : split_line[6],
      "TSSpos" : start_offset + int(split_line[4]),
      "score" : split_line[5],
      "box10pos" : start_offset + int(extra_field_dict["box10pos"]),
      "box35pos" : start_offset + int(extra_field_dict["box35pos"]),
      "box10seq" : extra_field_dict["box10seq"],
      "box35seq" : extra_field_dict["box35seq"],
      }
    
    # @JEB: Note for uniqueness, sometimes TSSpos are off by one nucleotide for unknown reaasons,
    #       but the -10 and -30 positions always seem to be correct
    new_key = new_entry["promoter"] + new_entry["strand"] + str(new_entry["box35pos"]) + str(new_entry["box10pos"])
    if not (new_key in existing_predictions):
      existing_predictions[new_key] = 1
      entry_list.append(new_entry)
    
  #sort list by coordinate
  sorted_entry_list = sorted(entry_list, key=itemgetter('TSSpos')) 

  with open(options.o + '.final_predictions.tsv','w') as final_predictions_file:
    writer = csv.DictWriter(
        final_predictions_file, 
        delimiter = '\t',
        fieldnames = ["promoter", "score", "strand", "TSSpos", "box35pos", "box35seq", "box10pos", "box10seq",]
      )
    writer.writeheader()
    writer.writerows(sorted_entry_list)
  final_predictions_file.close()

