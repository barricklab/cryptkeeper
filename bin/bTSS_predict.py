#!/usr/bin/env python

#import argparse
import Bio
from Bio.SeqIO import FastaIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess
import csv

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

#prepare fasta sections of 500 bp for bTSSfinder 
fastaSplit={} #create dictionary of sequence regions 
if len(records.seq)%250 >0:
    records.seq = records.seq + ((250-(len(records.seq)%250)) * "N")
    for n in range(0,int((len(records.seq)/250))):
        if n == 0:
            query_a = records.seq[n:n+500] 
            #print(len(query_a))
            name ="CMV_RNA1_"+ str(n) +"_"+ str(n+500)
            seq=str(query_a)
            fastaSplit[name]=seq
            
        else:
            query_b = records.seq[n*250:n*250+500]
            #print(len(query_b))
            name = "CMV_RNA1_"+ str(n*250)+"_"+ str(n*250+500) 
            seq = str(query_b)
            fastaSplit[name]=seq
print(fastaSplit)
print(len(fastaSplit['CMV_RNA1_250_750']))

fullFastaList = []
for key in fastaSplit:
    fullFastaList.append(SeqRecord(Seq(fastaSplit[key]), id = key,description=""))
    
#write fasta for bTSSfinder
with open("RNA_summary.fasta", "w") as handle:
    SeqIO.write(fullFastaList, handle, "fasta")
    handle.close()

#run bTSSfinder 
subprocess.Popen('bTSSfinder -i RNA1_summary.fa -o RNA1_summary -h 2', shell = True,)
      
##how to manage full .bed file summary 
bedfile_summary = open("RNA1_summary.bed","r")
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

print(search_coords,sigma_factors,start_codons,sigma_site)
    
start_coords = [i.split('_',3)[2] for i in search_coords] 
print(start_coords)
end_coords = [i.split('_',3)[2]for i in search_coords]
print(end_coords)
sigma_one = [i.split(',',1)[0] for i in sigma_site]
print(sigma_one)
sigma_two = [i.split(',',1)[1] for i in sigma_site]
print(sigma_two)

#rectify start codon and binding sites 

for i in range(0,len(search_coords)):
    start_codons.append(int(start_codons[i]) + int(start_coords[i]))
    sigma_one.append(int(sigma_one[i]) + int(start_coords[i]))
    sigma_two.append(int(sigma_two[i]) + int(start_coords[i]))
print(start_codons)
print(sigma_one)
print(sigma_two)


with open('RNA1_summary_corrected.bed','w') as corrected_bed:
    writer = csv.write(corrected_bed, delimiter = '\t', newline = '\n')
    writer.writerows(zip(search_coords,start_coords,end_coords,sigma_factors,direction,start_codons, sigma_one,sigma_two))
    corrected_bed.close()

    
    