#!/usr/bin/env python

"""

@authors: Jeffrey Barrick

Uses TransTermHP v2.08.

Setup:

TransTerm requires the path to a file expterms.dat to be passed at the command line
this script assumes that file (expterm.dat) is located at $TRANSTERM_EXPDAT_PATH.
I recommend putting it under $HOME/local/share/transterm/exp.dat and installing the
executable at $HOME/local/bin

"""
import argparse
from Bio import SeqIO
import subprocess
import csv
from operator import itemgetter
import os
from logging import debug, info, warning, error, critical



def process_Trans_Term_calculator_output_file(input_file_name):
    entry_list = []

    last_line_was_prediction = False
    new_entry = {}

    lines = []
    with open(input_file_name) as f:
        lines = f.readlines()

    for line in lines:
        line = line.strip()

        # skip all but lines of interest
        if line == "":
            continue
        split_line = line.split()

        # Pick up the sequence line
        if last_line_was_prediction:
            new_entry["seq_upstream"] = split_line[0]
            new_entry["seq_hairpin_open"] = split_line[1]
            new_entry["seq_hairpin_loop"] = split_line[2]
            new_entry["seq_hairpin_close"] = split_line[3]
            new_entry["seq_tail"] = split_line[4]

            entry_list.append(new_entry)

            new_entry = {}
            last_line_was_prediction = False
            continue

        if len(split_line) < 10:
            continue
        if split_line[0] != 'TERM':
            continue

        # New terminator
        last_line_was_prediction = True
        new_entry = {}

        new_entry["strand"] = split_line[5]  # either '+' or '-'
        if new_entry["strand"] == '+':
            new_entry["start"] = int(split_line[2])
            new_entry["end"] = int(split_line[4])
        else:
            new_entry["start"] = int(split_line[4])
            new_entry["end"] = int(split_line[2])

        new_entry["conf"] = split_line[7]
        new_entry["hairpin_score"] = split_line[8]
        new_entry["tail_score"] = split_line[9]

    return entry_list

def main(options):
    i=0

    for this_seq in SeqIO.parse(options.i, "fasta"):
        i += 1
        if (i>1):
            exit()  # TODO: Asess whether or not exit is appropriate here?
        main_seq = this_seq.upper()

    #unique step, create dummy cords files so that entire sequence is downstream of genes on both strands
    with open(options.o + '.dummy.coords','w') as dummy_coords_file:
        dummy_coords_file.write('gene1 1 2 ' + main_seq.id + '\n')
        dummy_coords_file.write('gene2 ' + str(len(main_seq)) + ' ' + str(len(main_seq)-1) + ' ' + main_seq.id + '\n')

    transterm_expdat_path = ''
    transterm_path = subprocess.check_output(['which', 'transterm']).strip()
    if "TRANSTERM_EXPDAT_PATH" in os.environ:
        transterm_expdat_path = os.getenv("TRANSTERM_EXPDAT_PATH")
    elif transterm_path:
        try:
            # which transterm
            transterm_path = subprocess.check_output(['which', 'transterm']).strip()
            transterm_expdat_path = os.path.normpath(transterm_path)
            transterm_expdat_path = transterm_expdat_path.split(os.sep)
            transterm_expdat_path = os.sep.join(transterm_expdat_path[:-2]) + os.sep +  "data" + os.sep + 'expterm.dat'
        except:  # TODO: catch specific exception
            transterm_path = None
            
    if not transterm_expdat_path:
        warning("Could not detect transterm installation. Please verify that transterm is installed via conda or that the environment variable TRANSTERM_EXPDAT_PATH is set correctly.")
        return 2

    if isinstance(transterm_expdat_path, bytes):
        transterm_expdat_path = transterm_expdat_path.decode('utf-8')

    #run once for predictions on both strands
    print('transterm  -p ' + transterm_expdat_path + ' ' + options.i + options.o + '.dummy.coords > ' + options.o + '.predictions.txt')
    subprocess.call('transterm  -p ' + transterm_expdat_path + ' ' + options.i + ' ' + options.o + '.dummy.coords > ' + options.o + '.predictions.txt', shell = True)

    #returns a list of dictionaries for the rows
  
    # Parse output and create one summary file
    final_list = process_Trans_Term_calculator_output_file(options.o +'.predictions.txt')

    # Sort list by coordinate
    final_list = sorted(final_list, key=itemgetter('start'))

    with open(options.o,'w') as final_predictions_file:
        writer = csv.DictWriter(
          final_predictions_file,
          fieldnames = ["start", "end", "strand", "conf", "hairpin_score", "tail_score", "seq_upstream", "seq_hairpin_open", "seq_hairpin_loop", "seq_hairpin_close", "seq_tail"])
        writer.writeheader()
        writer.writerows(final_list)
    final_predictions_file.close()

def RIT_predict_TransTerm(input, output):
    # Pretend to be an argument parser
    class ObjectClass(object):
        pass
    options = ObjectClass
    options.i = input
    options.o = output
    main(options)

if __name__ == "__main__":
    # ------------------------------------------------------------------------------
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

    # ------------------------------------------------------------------------------

    options = parser.parse_args()
    main(options)
