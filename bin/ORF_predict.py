#!/usr/bin/env python

#import argparse
import Bio
from Bio import SeqIO

##------------------------------------------------------------------------------
#parser = argparse.ArgumentParser(description='import fasta for slicing')
#parser.add_argument('-i',  
#    action='store', 
#    dest='i',
#    required=True,
#    type=str,
#    help="input .fa file to slice upstream of bTSSfinder")
#
#parser.add_argument('-l',  
#    action='store', 
#    dest='i',
#    required=True,
#    type=str,
#    help="define length of ORF")
##------------------------------------------------------------------------------
#options = parser.parse_args()
records = SeqIO.read("CMV_RNA1.fa","fasta")
table = 11 #NCBI translation table for Bacterial, Archaeal, and Plant Plastids
min_pro_len = 50 #Arbitrary number based on tutorial, can change to argparse option later

#length of sequence = multiple of 3
if len(records.seq)%3 >0:
    records.seq = records.seq + ((3-(len(records.seq)%3)) * "N")
    
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
print(orf_list)
for start, end, strand, pro in orf_list:
    print("%s...%s - length %i, strand %i, %i:%i" \
          % (pro[:30], pro[-3:], len(pro), strand, start, end))

for i in range(0,len(direction)):
    if direction == "+":
        seq = records.seq[sigma_one[i]:len(seq)]
        forward_ORF = find_orfs_with_trans(seq,table,min_pro_len)
    else:
        revseq = records.seq[0:sigma_one[i]]
        revseq.reverse_complement()

orf_list = find_orfs_with_trans(records.seq, table, min_pro_len)
print(orf_list)

for start, end, strand, pro in orf_list:
    with open('ORF_summary.txt','w') as outfile:
        outfile.write("%s...%s - length %i, strand %i, %i:%i" \
          % (pro[:30], pro[-3:], len(pro), strand, start, end))
        outfile.close()
    
    