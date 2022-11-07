#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""

@authors: Jeffrey Barrick

Predict cryptic bacterial gene expression signals in an input sequence.

"""

import argparse
import os
import logging
import tempfile
from operator import itemgetter
from copy import deepcopy
from collections import namedtuple
from Bio import SeqIO
from rhotermpredict import rho_term_predict
from .orf_predict import orf_predict, find_orfs
from .dependency_wrappers import ostir, transterm, promocalc
from .export import CryptResults, plot, to_csv

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
    
    to_csv(result, options.o)
    plot(result, options.o + "graph.html")

def cryptkeeper(input_file, output=None, circular=False, name=None, rbs_score_cutoff=2.0):
    """Predict cryptic bacterial gene expression signals in an input sequence."""

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
    if circular:
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
                    lowest_level = max([int(n) for n in end_positions])
                    end_positions[str(lowest_level+1)] = feature.end
                    strand[i] = feature._replace(nest_level=lowest_level+1)
            for position in end_positions:
                end_positions[position] = "inf"

    # Load sequence
    i = 0
    main_seq = None
    for this_seq in SeqIO.parse(infile, "fasta"):
        i += 1
        if i > 1:
            exit()
        main_seq = this_seq.upper()


    # ------------------------------------------------------------------------------
    # PERFORM TERMINATOR PREDICTIONS
    # ------------------------------------------------------------------------------

    # Predict terminators

    logger.info('Running RhoTermPredict')
    rhotermpredict_predictions = rho_term_predict(inseq=[str(main_seq.seq)])  # Run RhoTermPredict
    # Load RhoTermPredict output into named tuple
    rdtresult = namedtuple('rdtresult', """strand,
                                        c_over_g,
                                        start_rut,
                                        end_rut,
                                        term_seq,
                                        downstream_seq,
                                        palindromes,
                                        pause_concensus,
                                        scores""")
    rhotermpredict_results = []
    for result in rhotermpredict_predictions:
        rhotermpredict_results.append(rdtresult(result.strand,
                                                result.c_over_g,
                                                result.start_rut,
                                                result.end_rut,
                                                result.term_seq,
                                                result.downstream_seq,
                                                result.palindromes,
                                                result.pause_concensus,
                                                result.scores))


    logger.info('Running RIT_predict_TransTerm')
    transterm_predictions = transterm(input_file_name)
    # ------------------------------------------------------------------------------
    # PERFORM PROMOTER PREDICTION
    # ------------------------------------------------------------------------------

    # Predict promoters
    logger.info('Running TSS_predict_promoter_calculator')
    promoter_calc_predictions = promocalc(main_seq)

    # ------------------------------------------------------------------------------
    # PERFORM orf PREDICTION / GENERATE ANNOTATION FILE
    # ------------------------------------------------------------------------------

    # Predict orfs

    logger.info('Running ORF prediction')
    orf_predictions = orf_predict(input_file_name, output = None)
    [hit.update( {'start_codon_position': int(hit['start'] if hit['strand'] == '+' else hit['end'])}) for hit in orf_predictions]

    # Re-sort by the start codon position so we are in the same order as RBS
    orf_predictions = sorted(orf_predictions, key=itemgetter('start_codon_position'))

    # ------------------------------------------------------------------------------
    # PERFORM RBS PREDICTION
    # ------------------------------------------------------------------------------

    # Predict RBS
    rbs_predictions = ostir(input_file_name)

    # ------------------------------------------------------------------------------

    # Identify expressed orfs
    expressed_orfs = []
    expressed_orf = namedtuple("orf_result", "start, end, expression, burden, dG, array, start_codon, strand")
    orfs = find_orfs(input_file_name, 11, 0)
    start_codons_fwd = {}   # Value is index in rbs_predictions. Havent checked, probably faster than looping dicts.
    start_codons_rev = {}

    total_burden = 0

    for i, prediction in enumerate(rbs_predictions):
        if prediction.strand == '+':
            start_codons_fwd[prediction.position] = i
        else:
            start_codons_rev[prediction.position] = i

    for orf in orfs:

        if orf['strand'] == '+':
            expressed_starts = start_codons_fwd
            array_minus = 0
            orf_array = int(orf['end']) - int(orf['start']) + 1
            offset = 1
        else:
            expressed_starts = start_codons_rev
            array_minus = int(orf['start']) - int(orf['end']) + 1
            orf_array = 0
            offset = -1

        orf_info = orf
        orf_length = abs(orf_info['end'] - orf_info['start']) + 1
        if orf['start']+offset in expressed_starts.keys():
            rbs_info = rbs_predictions[expressed_starts[orf['start']+offset]]
            rbs_score = rbs_info.score
            if str(rbs_score) in ['inf', '-inf']:
                continue

            assert rbs_info.start_codon == orf_info['start_codon']
            
            array = max(int(orf['end']) - int(orf['start']), int(orf['start']) - int(orf['end'])) + 1
            processed_orf = expressed_orf(int(orf_info['start']),
                                          int(orf_info['end']),
                                          rbs_score,
                                          orf_length * rbs_score,
                                          rbs_info.score2,
                                          array,
                                          rbs_info.start_codon,
                                          orf_info["strand"]
                                          )
            expressed_orfs.append(processed_orf)

            total_burden += orf_length * rbs_score

    #   When available, assign RBS info to orf
            orf['array'] = orf_array
            orf['array_minus'] = array_minus
            orf['rbs_score'] = 10**rbs_info.score

        else:
            rbs_score = 0
            orf['array'] = orf_array
            orf['array_minus'] = array_minus
            orf['rbs_score'] = rbs_score

        if circular:  # !!! BUG: IF TRUE THIS WILL ERROR @TODO
            if (int(orf_info['start']) <= sequence_length) or (int(orf_info['end']) >= sequence_length * 2 + 1):
                pass


    # Filter predictions that we show
    expressed_orfs = list(filter(lambda x:float(x.expression) > rbs_score_cutoff, expressed_orfs))

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

        for rbs in rbs_predictions:
            rbs['position'] = str(int(rbs['position']) - sequence_length)

        for rbs in promoter_calc_predictions:
            rbs['TSSpos'] = str(int(rbs['TSSpos']) - sequence_length)

        for rbs in transterm_predictions:
            rbs.end = str(int(rbs.end) - sequence_length)

        # Split reading frames and assign
        added_split_orfs = []
        for rbs in orf_predictions:
            rbs['start'] = str(int(rbs['start']) - sequence_length)
            rbs['end'] = str(int(rbs['end']) - sequence_length)

            if int(rbs['end']) > sequence_length:

                new_split_orf = deepcopy(rbs)
                rbs['end'] = str(sequence_length)
                new_split_orf['start'] = "1"
                new_split_orf['end'] = str(int(new_split_orf['end']) - sequence_length)

                added_split_orfs.append(new_split_orf)
        orf_predictions.extend(added_split_orfs)


    # Set up the results object

    result = CryptResults(name = name,
                          sequence = main_seq.seq,
                          translation_sites = expressed_orfs,
                          row_dep_terminators = rhotermpredict_results,
                          row_ind_terminators = transterm_predictions,
                          promoters =  promoter_calc_predictions,
                          annotations = features_list,
                          burden = total_burden)

    return result


if __name__ == "__main__":
    main()
