#!/usr/bin/env python

"""
Created on Thu Oct  18 12:26:00 2018

@authors: Alexandra Lukasiewicz
Purpose: Predict and translate ORF regions in query.
Code adapted from Biopython 1.72 tutorial and cookbook, 20.1.13
"""

import argparse
import Bio
from Bio import SeqIO
import Bio.Data.CodonTable
from operator import itemgetter
import csv
from math import floor

def find_orfs(seq, translation_table_id, minimum_orf_aa_length):
  orfs = []
  seq_len = len(seq)

  #Get the codon table so we know the valid start codons
  translation_table = Bio.Data.CodonTable.unambiguous_dna_by_id[translation_table_id]

  #ignore ultra-rare "UUG" "GUG" start codons (which are not recognized by RBS Calculator)
  if (translation_table_id==11):
    translation_table.start_codons = ['ATG', 'GTG', 'TTG']

  #print(translation_table.start_codons)

  for this_strand, this_seq in [('+', seq), ('-', seq.reverse_complement())]:
    for start_pos_0 in range(len(this_seq)- 2):

      #print(this_seq[start_pos_0:start_pos_0+3].seq)

      #Is this a start codon?
      this_start_codon = this_seq[start_pos_0:start_pos_0+3].seq
      if (this_start_codon not in translation_table.start_codons):
        continue

      start_pos_1 = start_pos_0+1

      #print(this_seq[start_pos_0:])

      #remove bases to make a multiple of three to avoid BioPython warnings
      end_pos_1 = start_pos_0+floor(float(len(this_seq)-start_pos_0)/3.0)*3

      aa_sequence = this_seq[start_pos_0:end_pos_1].seq.translate(translation_table, to_stop=True)
      #print(aa_sequence)
      aa_length = len(aa_sequence)

      end_pos_1 = start_pos_1 + aa_length*3 - 1

      #Is it a long enough reading frame?
      if aa_length < minimum_orf_aa_length:
        continue

      #print(str(this_strand) + " " + str(start_pos_0))

      if (this_strand == '-'):
        start_pos_1 = len(this_seq) - start_pos_1 + 1
        end_pos_1 = len(this_seq) - end_pos_1 + 1
        start_pos_1, end_pos_1 = end_pos_1, start_pos_1

      orfs.append(dict(
        start = start_pos_1,
        end = end_pos_1,
        strand = this_strand,
        start_codon = this_start_codon,
        length = aa_length
        ))

  orfs = sorted(orfs, key=itemgetter('start'))
  return orfs


def main(options):
    translation_table_id = options.t #NCBI translation table for Bacterial, Archaeal, and Plant Plastids
    minimum_orf_aa_length = options.l

    i=0
    main_seq = None
    for this_seq in SeqIO.parse(options.i, "fasta"):
      i += 1
      if (i>1):
        exit()
      main_seq = this_seq.upper()


    orfs = find_orfs(main_seq, translation_table_id, minimum_orf_aa_length)

    #write each tuple of list to outfile
    with open(options.o,'w') as final_predictions_file:
      writer = csv.DictWriter(
          final_predictions_file,
          fieldnames = ["start", "end", "strand", "start_codon", "length"]
        )
      writer.writeheader()
      writer.writerows(orfs)
    final_predictions_file.close()

def ORF_predict(input, output, transtable=11, minlength=0):
    # Pretend to be an argument parser
    class ObjectClass:
        pass
    options = ObjectClass
    options.i = input
    options.o = output
    options.t = transtable
    options.l = minlength
    main(options)

if __name__ == "__main__":
    # ------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='import fasta for ORF detection')
    parser.add_argument('-i',
                        action='store',
                        dest='i',
                        required=True,
                        type=str,
                        help="input .fasta file for ORF search")

    parser.add_argument('-t',
                        action='store',
                        dest='t',
                        required=False,
                        default=11,
                        type=str,
                        help="set NCBI translation table, default = 11: Bacterial, Archaeal, and Plant Plastids")

    parser.add_argument('-l',
                        action='store',
                        dest='l',
                        required=False,
                        default=30,
                        type=str,
                        help="set minimum length of protein (in amino acids)")

    parser.add_argument('-o',
                        action='store',
                        dest='o',
                        required=True,
                        type=str,
                        help="outfile prefix")
    # ------------------------------------------------------------------------------
    options = parser.parse_args()

    main(options)
