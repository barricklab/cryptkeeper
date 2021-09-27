#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

@authors: Jeffrey Barrick

Uses ViennaRNA 2.4.13.

Setup:

1) Running the RBA calculator requires the ViennaRNA.


"""
import argparse
from Bio import SeqIO
import csv
from operator import itemgetter
import os
from helpers import timer_decorator


def process_RBS_calculator_output_file(input_file_name, is_reverse_complement,
                                       sequence_length, forward_seq, reverse_seq):
    '''returns a list of dictionaries for the rows'''
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
        #     1 or 'Warning: There is a leaderless start codon, which is being ignored.'
        if line_index % 2 == 0:

            failed_prediction = False
            #append em
            if line == "'Warning: There is a leaderless start codon, which is being ignored.'":
                failed_prediction = True

        # Second line of pair
        #    pos_1, score, score2
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

            # RBS calculator returns 1-indexed positions
            pos_1 = int(split_line[0]) + 1
            new_entry["position"] =    pos_1 if not is_reverse_complement else sequence_length - pos_1 + 1
            new_entry["start_codon"] =    forward_seq[pos_1-2:pos_1+1] if not is_reverse_complement else reverse_seq[pos_1-2:pos_1+1]
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


    #@TODO: Implement testing requested stop codons and minimum binding screening
        # Three modes possible
    if False:
        return
    else:
            #Run OSTIR Calculator twice on entire sequences. Once for each strand.
            from ostir.ostir import run_ostir

            findings_forward = run_ostir(forward_seq)
            findings_reverse = run_ostir(reverse_seq)

            for finding in findings_forward:
                    finding['strand'] = "+"
            for finding in findings_reverse:
                    finding['strand'] = '-'

            result = findings_forward + findings_reverse

            clean_findings = []    # Remove duplicates
            [clean_findings.append(finding) for finding in result if finding not in clean_findings]
            # Just cleans the output files
            with open(f"{options.o}.forward.predictions.txt", 'w') as _:
                    pass
            with open(f"{options.o}.reverse.predictions.txt", 'w') as _:
                    pass
            forward_file = f"{options.o}.forward.predictions.txt"
            reverse_file = f"{options.o}.reverse.predictions.txt"
            for finding in clean_findings:
                    expr = finding['expression']
                    start_pos = finding['start_position']
                    ks = finding['dG_total']
                    strand = finding['strand']
                    if strand == '+':
                            outfile = forward_file
                    else:
                            outfile = reverse_file
                    with open(outfile, 'a') as file_to_write:
                            file_to_write.writelines(['1\n', f"{start_pos} {expr} {ks}\n"])
    # Parse output and create one summary file
    forward_list = process_RBS_calculator_output_file(options.o +'.forward.predictions.txt', False, sequence_length, forward_seq, reverse_seq)
    reverse_list = process_RBS_calculator_output_file(options.o +'.reverse.predictions.txt', True, sequence_length, forward_seq, reverse_seq)
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
