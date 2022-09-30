#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""

@authors: Jeffrey Barrick

Predict cryptic bacterial gene expression signals in an input sequence.

"""

import argparse
import csv
import os
import logging
import tempfile

from Bio import SeqIO
from operator import itemgetter
from copy import deepcopy
from collections import namedtuple

from rhotermpredict import rho_term_predict
from .ORF_predict import ORF_predict, find_orfs
from .dependency_wrappers import ostir, transterm, promocalc
from .export import CryptResults, plot

def main():
    """CLI Entry Point for Cryptkeeper"""
    # ------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='input fasta')
    parser.add_argument(
        '-i', '--input',
        action='store',
        dest='i',
        required=True,
        type=str,
        help="input fasta file",
    )

    parser.add_argument(
        '-c', '--circular',
        action='store_true',
        dest='circular',
        help="The input file is circular. (Note: Increases runtime by 3x)",
    )

    parser.add_argument(
        '-o', '--output',
        action='store',
        dest='o',
        required=True,
        type=str,
        help="output file prefix",
    )

    parser.add_argument(
        '-p', '--plot-only',
        action='store_true',
        dest='plot_only',
        help="plot mode, assumes output files all exist",
    )

    parser.add_argument(
        '-n', '--name',
        action='store',
        dest='name',
        help="name of sample (if not provided the filename is used)",
    )

    parser.add_argument(
        '--rbs-score-cutoff',
        action='store',
        dest='rbs_score_cutoff',
        default=2.0,
        type=float,
        help="Minimum score that is graphed and output to final files (all are used in calculating burden)",
    )

    # parser.add_argument(
    #    '--pdf',
    #    action='store_true',
    #    dest='pdf',
    #    help="write output plot as *.pdf",
    #    )

    # ------------------------------------------------------------------------------
    options = parser.parse_args()

    result = cryptkeeper(input_file = options.i,
                        output=options.o,
                        circular=options.circular,
                        name=options.name,
                        rbs_score_cutoff=options.rbs_score_cutoff)

    plot(result, options.o + "graph.html")

def cryptkeeper(input_file, output=None, circular=False, name=None, rbs_score_cutoff=2.0):

    # @TODO: Many of these steps should be moved to a function

    # ------------------------------------------------------------------------------
    # PROCESS INPUT FILES
    # ------------------------------------------------------------------------------

    # Makes relative file arguments absolute for more robust execution

    input_sequence = os.path.abspath(input_file)
    input_gene_name = input_sequence.split("/")[-1]
    input_file_type = input_gene_name.split(".")[-1]
    input_gene_name = ".".join(input_gene_name.split(".")[:-1])
    input_file_name = input_sequence
    
    if output:
        output_path = os.path.abspath(output)
        output_folder = "/".join(output_path.split("/")[:-1])
        # Prevent execution if output is a file instead of a folder
        if os.path.isfile(output_folder):
            raise ValueError(f"Output folder {output_folder} is a file. Please provide a folder path.")
        # Create the folder if it doesnt exist
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
    else: # If not output
        output_path = tempfile.tempdir()
    

    # Set up log
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    stream_formatter = logging.Formatter('[%(asctime)s] Cryptkeeper: %(message)s')
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(stream_formatter)
    logger.addHandler(stream_handler)

    if output:
        file_formatter = logging.Formatter('[%(asctime)s - %(levelname)s] Cryptkeeper: %(message)s')
        file_handler = logging.FileHandler(filename=f"{output_path}.log")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

    # Convert non-FASTA files to FASTA (supports Genbank) - This is necessary to run TransTermHP
    if input_file_type != "fasta" and input_file_type != "fna":
        logger.info("Non-FASTA file detected. Converting.")
        if input_file_type in ['genbank', 'gb', 'gbk']:
            with open(input_file_name, "r", encoding="utf-8" ) as file_in:
                with open(f"{output_path + '.' + input_gene_name}.fna", "w", encoding="utf-8") as file_converted:
                    sequences = SeqIO.parse(file_in, "genbank")
                    sequences = SeqIO.write(sequences, file_converted, "fasta")
            infile = f"{output_path + '.' + input_gene_name}.fna"
        else:
            raise ValueError(f"Unable to convert file ({input_file_type}?). Supported types: FASTA, GENBANK")
    else:
        infile = input_sequence
    sequences = list(SeqIO.parse(infile, "fasta"))
    input_file_name = infile

    # --- CIRCULAR SEQUENCE ---
    # create a 3x version of the FASTA input file

    output_circular_fasta_file_name = output_path + '.circular.3x.fasta'
    if len(sequences) > 1:
        logger.warning("Input file contains multiple sequences. Only the first will be processed.")
    sequence_length = len(sequences[0])
    if (circular):
        with open(output_circular_fasta_file_name, "w", encoding="utf-8" ) as output_handle:
            sequences[0].seq = sequences[0].seq + sequences[0].seq + sequences[0].seq
            SeqIO.write([sequences[0]], output_handle, "fasta")
            input_file_name = output_circular_fasta_file_name


   # Extract annotation information from genbank files
    feature_tuple = namedtuple('feature', ['name', 'strand', 'start', 'end', 'color', 'nest_level'])
    features_list = []
    if input_file_type in ['genbank', 'gb', 'gbk']:

        features = [feature for feature in [rec for rec in SeqIO.parse(input_sequence, "genbank")][0].features]
        for feature in features:
            if 'ApEinfo_revcolor' in feature.qualifiers:
                color = feature.qualifiers['ApEinfo_revcolor'][0]
            else:
                color = "black"

            feature_hit = feature_tuple(feature.qualifiers['label'][0],
                                        feature.location.strand,
                                        feature.location.start,
                                        feature.location.end,
                                        color,
                                        None)
            features_list.append(feature_hit)


        # Process nest level for displaying features
        strandless_features = [feature for feature in features_list if feature.strand not in [1, -1]]
        forward_features = [feature for feature in features_list if feature.strand == 1]
        reverse_features = [feature for feature in features_list if feature.strand == -1]
        end_positions = dict()
        #for strand in [strandless_features, forward_features, reverse_features]:  Dont care about organizing this
        for strand in [features_list]:
            for i, feature in enumerate(strand):
                if not end_positions:
                    end_positions[str(0)] = feature.end
                    strand[i] = feature._replace(nest_level=0)
                    continue
                for level, position in end_positions.items():
                    level = int(level)
                    if position == 'inf':
                        continue
                    elif int(position) > int(feature.start):
                        continue
                    else:
                        end_positions[str(level)] = feature.end
                        strand[i] = feature._replace(nest_level=level)
                        break
                if not strand[i].nest_level and strand[i].nest_level != 0:  # Layer 0 is falcey
                    lowest_level = max([int(n) for n in end_positions.keys()])
                    end_positions[str(lowest_level+1)] = feature.end
                    strand[i] = feature._replace(nest_level=lowest_level+1)
            for position in end_positions.keys():
                end_positions[position] = "inf"

    # Load sequence
    i = 0
    main_seq = None
    for this_seq in SeqIO.parse(infile, "fasta"):
        i += 1
        if i > 1:
            exit()
        main_seq = this_seq.upper()
        seq_name = this_seq.id


    # ------------------------------------------------------------------------------
    # PERFORM TERMINATOR PREDICTIONS
    # ------------------------------------------------------------------------------

    # Predict terminators

    logger.info('Running RhoTermPredict')
    rhotermpredict_predictions = rho_term_predict(inseq=[str(main_seq.seq)])  # Run RhoTermPredict
    # Load RhoTermPredict output into named tuple
    logger.info('Running RIT_predict_TransTerm')
    transterm_predictions = transterm(input_file_name)
    # ------------------------------------------------------------------------------
    # PERFORM PROMOTER PREDICTION
    # ------------------------------------------------------------------------------

    # Predict promoters
    logger.info('Running TSS_predict_promoter_calculator')
    promoter_calc_predictions = promocalc(main_seq)

    # ------------------------------------------------------------------------------
    # PERFORM ORF PREDICTION / GENERATE ANNOTATION FILE
    # ------------------------------------------------------------------------------

    # Predict ORFs

    logger.info('Running ORF prediction')
    orf_predictions = ORF_predict(input_file_name, output = None)
    [hit.update( {'start_codon_position': int(hit['start'] if hit['strand'] == '+' else hit['end'])}) for hit in orf_predictions]

    # Re-sort by the start codon position so we are in the same order as RBS
    orf_predictions = sorted(orf_predictions, key=itemgetter('start_codon_position'))

    # ------------------------------------------------------------------------------
    # PERFORM RBS PREDICTION
    # ------------------------------------------------------------------------------

    # Predict RBS
    rbs_predictions = ostir(input_file_name)

    # ------------------------------------------------------------------------------

    # Identify expressed ORFs
    expressed_ORFs = []
    expressed_orf = namedtuple("orf_result", "start, end, score, burden, array, array_minus, second_half,")
    ORFs = find_orfs(input_file_name, 11, 0)
    start_codons_fwd = {}   # Value is index in rbs_predictions. Havent checked, probably faster than looping dicts.
    start_codons_rev = {}

    total_burden = 0

    for i, prediction in enumerate(rbs_predictions):
        if prediction.strand == '+':
            start_codons_fwd[prediction.position] = i
        else:
            start_codons_rev[prediction.position] = i

    for ORF in ORFs:

        if ORF['strand'] == '+':
            expressed_starts = start_codons_fwd
            array_minus = 0
            ORF_array = int(ORF['end']) - int(ORF['start']) + 1
            offset = 1
        else:
            expressed_starts = start_codons_rev
            array_minus = int(ORF['start']) - int(ORF['end']) + 1
            ORF_array = 0
            offset = -1

        ORF_info = ORF
        orf_length = abs(ORF_info['end'] - ORF_info['start']) + 1
        if ORF['start']+offset in expressed_starts.keys():
            RBS_info = rbs_predictions[expressed_starts[ORF['start']+offset]]
            rbs_score = RBS_info.score
            if str(rbs_score) in ['inf', '-inf']:
                continue

            if RBS_info.start_codon != ORF_info["start_codon"]:
                logger.debug(f'Codon mismatch {RBS_info.position} {RBS_info.start_codon} and {ORF_info["start"]} {ORF_info["start_codon"]}, {ORF_info["strand"]}\n{ORF_info["translation"]}')

            array = max(int(ORF['end']) - int(ORF['start']), int(ORF['start']) - int(ORF['end'])) + 1
            expressed_ORFs.append({
                "start": int(ORF_info['start']),
                "end": int(ORF_info['end']),
                "score": rbs_score,
                "burden": orf_length * rbs_score,
                "array": array,
                "start_codon": RBS_info.start_codon,
                "strand": ORF_info["strand"]
            })

            total_burden += orf_length * rbs_score

    #   When available, assign RBS info to ORF
            ORF['array'] = ORF_array
            ORF['array_minus'] = array_minus
            ORF['rbs_score'] = 10**RBS_info.score

        else:
            rbs_score = 0
            ORF['array'] = ORF_array
            ORF['array_minus'] = array_minus
            ORF['rbs_score'] = rbs_score

        if circular:  # !!! BUG: IF TRUE THIS WILL ERROR @TODO
            if (int(ORF_info['start']) <= sequence_length) or (int(ORF_info['end']) >= sequence_length * 2 + 1):
                add_to_score = False


    # Filter predictions that we show
    expressed_ORFs = list(filter(lambda x:float(x["score"]) > rbs_score_cutoff, expressed_ORFs))

    # ---- CIRCULAR SEQUENCE ----

    if circular:

        # Remove ones with out of bounds starts/ends/positions (whichever is graphed)
        orf_predictions = list(filter(lambda x: (int(x['start']) >= sequence_length + 1) and
                                                (int(x['start']) <= sequence_length * 2), orf_predictions))
        rbs_predictions = list(filter(lambda x: (int(x['position']) >= sequence_length + 1) and
                                                (int(x['position']) <= sequence_length * 2), rbs_predictions))
        promoter_calc_predictions = list(filter(lambda x: (int(x['TSSpos']) >= sequence_length + 1) and
                                                (int(x['TSSpos']) <= sequence_length * 2), promoter_calc_predictions))
        transterm_predictions = list(filter(lambda x: (int(x.end) >= sequence_length + 1) and
                                                (int(x.end) <= sequence_length * 2), transterm_predictions))

        # Correct coordinates of remaining

        for d in rbs_predictions:
            d['position'] = str(int(d['position']) - sequence_length)

        for d in promoter_calc_predictions:
            d['TSSpos'] = str(int(d['TSSpos']) - sequence_length)

        for d in transterm_predictions:
            d.end = str(int(d.end) - sequence_length)

        # Split reading frames and assign
        added_split_orfs = []
        for d in orf_predictions:
            d['start'] = str(int(d['start']) - sequence_length)
            d['end'] = str(int(d['end']) - sequence_length)

            if int(d['end']) > sequence_length:

                new_split_orf = deepcopy(d)
                d['end'] = str(sequence_length)
                new_split_orf['start'] = "1"
                new_split_orf['end'] = str(int(new_split_orf['end']) - sequence_length)

                added_split_orfs.append(new_split_orf)
        orf_predictions.extend(added_split_orfs)


    # Set up the results object

    result = CryptResults(name = name,
                          sequence = main_seq.seq,
                          transcription_sites = expressed_ORFs,
                          row_dep_terminators = rhotermpredict_predictions,
                          row_ind_terminators = transterm_predictions,
                          promoters =  promoter_calc_predictions,
                          annotations = features_list,
                          burden = total_burden)

    return result


if __name__ == "__main__":
    main()
