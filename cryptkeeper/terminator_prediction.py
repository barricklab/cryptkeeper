#!/usr/bin/env python

"""

@authors: Jeffrey Barrick, Cameron Roots

Uses TransTermHP and RhoTermPredict
"""

import argparse
import subprocess
import csv
import os
from Bio import SeqIO
from operator import itemgetter
from rhotermpredict import rho_term_predict
from logging import debug, info, warning, error, critical
from collections import namedtuple


def predict_transterm(infile, outpath):
    """Cryptkeeper wrapper for TransTermHP"""

    def process_transterm_output(input_file_name):
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

    i=0

    for this_seq in SeqIO.parse(infile, "fasta"):
        i += 1
        if (i>1):
            exit()  # TODO: Asess whether or not exit is appropriate here?
        main_seq = this_seq.upper()

    #unique step, create dummy cords files so that entire sequence is downstream of genes on both strands
    with open(outpath + '.dummy.coords','w') as dummy_coords_file:
        dummy_coords_file.write('gene1 1 2 ' + main_seq.id + '\n')
        dummy_coords_file.write('gene2 ' + str(len(main_seq)) + ' ' + str(len(main_seq)-1) + ' ' + main_seq.id + '\n')

    transterm_expdat_path = ''
    transterm_path = subprocess.check_output(['which', 'transterm']).strip()
    if "TRANSTERM_EXPDAT_PATH" in os.environ:
        transterm_expdat_path = os.getenv("TRANSTERM_EXPDAT_PATH")
    elif transterm_path:
        # Find the expterm.dat file
        transterm_path = subprocess.check_output(['which', 'transterm']).strip()
        transterm_expdat_path = os.path.normpath(transterm_path).decode('utf-8')
        transterm_expdat_path = transterm_expdat_path.split(os.sep)
        transterm_expdat_path = os.sep.join(transterm_expdat_path[:-2]) + os.sep +  "data" + os.sep + 'expterm.dat'

    if not transterm_expdat_path:
        warning("Could not detect transterm installation. Please verify that transterm is installed via conda or that the environment variable TRANSTERM_EXPDAT_PATH is set correctly.")
        return 2

    if isinstance(transterm_expdat_path, bytes):
        transterm_expdat_path = transterm_expdat_path.decode('utf-8')

    #run once for predictions on both strands
    print('transterm --min-conf=70  -p ' + transterm_expdat_path + ' ' + infile + outpath + '.dummy.coords > ' + outpath + '.predictions.txt')
    subprocess.call('transterm --min-conf=70  -p ' + transterm_expdat_path + ' ' + infile + ' ' + outpath + '.dummy.coords > ' + outpath + '.predictions.txt', shell = True)

    #returns a list of dictionaries for the rows
  
    # Parse output and create one summary file
    final_list = process_transterm_output(outpath +'.predictions.txt')

    # Sort list by coordinate
    final_list = sorted(final_list, key=itemgetter('start'))

    with open(outpath,'w') as final_predictions_file:
        writer = csv.DictWriter(
          final_predictions_file,
          fieldnames = ["start", "end", "strand", "conf", "hairpin_score", "tail_score", "seq_upstream", "seq_hairpin_open", "seq_hairpin_loop", "seq_hairpin_close", "seq_tail"])
        writer.writeheader()
        writer.writerows(final_list)
    final_predictions_file.close()
    
    return final_list

def predict_rhotermpredict(input_file_name, RhoTermPredict_out):
    """Cryptkeeper wrapper for RhoTermPredict"""
    rho_term_predict(input_file_name, RhoTermPredict_out)  # Run RhoTermPredict 
    # Load RhoTermPredict output into named tuple
    RTP_prediction = namedtuple("RTP_prediction", ["start", "end", "strand", "conf", "hairpin_score", "tail_score", "seq_upstream", "seq_hairpin_open", "seq_hairpin_loop", "seq_hairpin_close", "seq_tail"])
    # Do some processing
    final_list = None
    return final_list

def _unify_predictions(transterm_predictions, rhotermpredict_predictions):
    """Unifies predictions from TransTermHP and RhoTermPredict"""
    pass
