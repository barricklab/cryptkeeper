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

records = SeqIO.read(options.i , "fasta")

#prepare fasta sections of 500 bp for bTSSfinder 
fastaSplit={} #create dictionary of sequence regions 
if len(records.seq)%250 >0:
    records.seq = records.seq + ((250-(len(records.seq)%250)) * "N")
    for n in range(0,int((len(records.seq)/250))):
        if n == 0:
            query_a = records.seq[n:n+500] 
            #print(len(query_a))
            name =options.o+ str(n) +"_"+ str(n+500)
            seq=str(query_a)
            fastaSplit[name]=seq
            
        else:
            query_b = records.seq[n*250:n*250+500]
            #print(len(query_b))
            name = options.o+ str(n*250)+"_"+ str(n*250+500) 
            seq = str(query_b)
            fastaSplit[name]=seq

fullFastaList = []
for key in fastaSplit:
    fullFastaList.append(SeqRecord(Seq(fastaSplit[key]), id = key,description=""))
    
#write fasta for bTSSfinder
with open(options.o + "_summary.fasta", "w") as handle:
    SeqIO.write(fullFastaList, handle, "fasta")
    handle.close()

#run bTSSfinder 
subprocess.call('bTSSfinder -i '+ options.o + '_summary.fasta -o '+options.o+' -h 2', shell = True)
      
##how to manage full .bed file summary 
bedfile_summary = open(options.o + ".bed","r")
lines = bedfile_summary.readlines()
search_coords = []
sigma_factors = []
start_codons = []
sigma_site = []
direction = []
for line in lines:
    data = line.strip().split('\t')
    search_coords.append(data[0])
    sigma_factors.append(data[3])
    start_codons.append(data[6])
    sigma_site.append(data[11])
    direction.append(data[5])
bedfile_summary.close()
    
start_coords = [i.split('_',3)[1] for i in search_coords] 
end_coords = [i.split('_',3)[2]for i in search_coords]
print(end_coords)
sigma_one = [i.split(',',1)[0] for i in sigma_site]
sigma_two = [i.split(',',1)[1] for i in sigma_site]

#rectify start codon and binding sites 
corrected_start_codon = []
corrected_sigma_one = []
corrected_sigma_two = []

for i in range(0,len(search_coords)):
    corrected_start_codon.append(int(start_codons[i]) + int(start_coords[i]))
    corrected_sigma_one.append(int(sigma_one[i]) + int(start_coords[i]))
    corrected_sigma_two.append(int(sigma_two[i]) + int(start_coords[i]))
print(corrected_start_codon,corrected_sigma_one,corrected_sigma_two)

with open(options.o + '_corrected.txt','w') as corrected_bed:
    writer = csv.writer(corrected_bed, delimiter = '\t')
    writer.writerows(zip(search_coords,start_coords,end_coords,sigma_factors,direction,corrected_start_codon,corrected_sigma_one,corrected_sigma_two))
corrected_bed.close()

    
    

