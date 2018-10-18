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

#------------------------------------------------------------------------------
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
    default=100,
    type=str,
    help="set minimum protein length of translated product")

parser.add_argument('-o', 
    action='store',
    dest= 'o', 
    required=True,
    type=str,
    help="outfile prefix")
#------------------------------------------------------------------------------
options = parser.parse_args()
records = SeqIO.read(options.i,"fasta")
table = options.t #NCBI translation table for Bacterial, Archaeal, and Plant Plastids
min_pro_len = options.l 


#length of sequence = multiple of 3
#if len(records.seq)%3 >0:
#    records.seq = records.seq + ((3-(len(records.seq)%3)) * "N")
    
def find_orfs_with_trans(seq, trans_table, min_protein_length):
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = nuc[frame:].translate(trans_table)
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end-aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame+aa_start*3
                        end = min(seq_len,frame+aa_end*3+3)
                    else:
                        start = seq_len-frame-aa_end*3-3
                        end = seq_len-frame-aa_start*3
                    answer.append((start, end, strand,
                                   trans[aa_start:aa_end]))
                aa_start = aa_end+1
    answer.sort()
    return answer

orf_list = find_orfs_with_trans(records.seq, table, min_pro_len)

#write each tuple of list to outfile
with open(options.o +'.txt','w') as outfile:
    outfile.write('\n'.join("%s...%s - length %i, strand %i, %i:%i" \
    % (pro[:30], pro[-3:], len(pro), strand, start, end) for (start, end, strand, pro) in orf_list))
    outfile.close()
    
    
    
    