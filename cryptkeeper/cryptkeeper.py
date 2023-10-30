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
import numpy as np
from collections import namedtuple
from operator import itemgetter
from copy import deepcopy
from collections import namedtuple
from Bio import SeqIO
from rhotermpredict import rho_term_predict
from .orf_predict import orf_predict, find_orfs
from .dependency_wrappers import ostir, transterm, promocalc
from .export import CryptResults, plot, to_csv, to_summary
from .export_bokeh import export_bokeh
from .constants import COLORS
from .helpers import delay_iterator
import webbrowser

def main() -> None:
    """CLI Entry Point for Cryptkeeper"""

    # Create the argument parser
    parser = argparse.ArgumentParser(description='input fasta')

    # Define the command line arguments
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
        '-j', '-threads',
        action='store',
        dest='j',
        required=False,
        type=int,
        help="Number of threads/processes to use",
        default=1
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
    parser.add_argument(
        '--matplotlib',
        action='store_true',
        dest='matplotlib',
        default=False,
        help="Use MatPlotLib visualization (deprecated)",
    )

    # Parse the command line arguments
    options = parser.parse_args()

    # Call the cryptkeeper function with the provided options
    result = cryptkeeper(input_file=options.i,
                        output=options.o,
                        circular=options.circular,
                        name=options.name,
                        threads=options.j,
                        rbs_score_cutoff=options.rbs_score_cutoff)

    # Perform additional operations on the result
    to_csv(result, options.o)
    to_summary(result, options.o + "_summary.txt")
    result.to_json(options.o + "_results.json")

    # Print a message indicating that the analysis is finished
    # Plot the result using the selected visualization method
    if options.matplotlib:
        DeprecationWarning('Matplotlib base plotting will be replaced with bokah soon')
        plot(result, options.o + "_graph.html")
    else:
        plot_filepath = export_bokeh(result, options.o + "_graph.html")
        webbrowser.open(plot_filepath)

    # Print a message indicating that the process is done


def cryptkeeper(input_file, output=None, circular=False, name=None, threads=1, rbs_score_cutoff=2.0):
    """Predict cryptic bacterial gene expression signals in an input sequence."""

    # @TODO: Many of these steps should be moved to a function

    # Sanity Checking
    if threads < 1:
        raise ValueError("Number of threads must be greater than 0")

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
        outdir = tempfile.TemporaryDirectory()
        output_path = outdir.name


    # Set up log
    logger = logging.getLogger('cryptkeeper')
    logger.setLevel(logging.DEBUG)

    stream_formatter = logging.Formatter('[%(asctime)s] Cryptkeeper: %(message)s')
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(stream_formatter)
    if not logger. hasHandlers():
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
        if input_file_type in ['genbank', 'gb', 'gbk', 'gbff']:
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
    single_sequence = deepcopy(sequences[0].seq)  # We need to keep ahold of the original sequence to revert from circularizing
    input_file_name = infile

    # --- CIRCULAR SEQUENCE ---
    # create a 3x version of the FASTA input file


    if len(sequences) > 1:
        logger.warning("Input file contains multiple sequences. Only the first will be processed.")


    sequence_length = len(sequences[0])

    circular_length = 0
    if circular:
        # If the input sequence is less than 200 bases, use the entire sequence
        # if the largest overlapping ORF is greater than 200 bases, use that
        # otherwise, use 200 bases
        fwd_stop_codons = ['TAA', 'TAG', 'TGA']
        fwd_start_codons = ['ATG', 'GTG', 'TTG']

        rev_stop_codons = ['TTA', 'CTA', 'TCA']
        rev_start_codons = ['CAT', 'CAC', 'CAA']

        circular_length = -1
        junction = sequences[0].seq[-3:] + sequences[0].seq[:3]

        # Find forward strand longest ORF
        for start_codons, stop_codons in [(fwd_start_codons, fwd_stop_codons), (rev_stop_codons, rev_start_codons)]: # Invert start and stop codons for reverse strand
            found_starts = [None, None, None]
            found_stops = [None, None, None]
            if junction[1:4] in start_codons:
                found_starts[1] = -1
            if junction[2:5] in start_codons:
                found_starts[2] = -1

            for e, i in enumerate(range(sequence_length, -1, -1)):
                if sequences[0].seq[i-3:i] in stop_codons:
                    frame = e % 3
                    found_stops[frame] = -1
                    found_starts[frame] = True
                elif sequences[0].seq[i-3:i] in start_codons:
                    frame = e % 3
                    if found_starts[frame] is None:
                        found_starts[frame] = True
                if all(found_starts):
                    break
            else: # If no start codon found fill in the gaps
                for frame, start in enumerate(found_starts):
                    if not start:
                        found_stops[i] = -1
                        found_starts[i] = True


            if junction[0:3] in stop_codons:
                found_stops[0] = -1
            if junction[1:4] in stop_codons:
                found_stops[1] = -1
            if junction[2:5] in stop_codons:
                found_stops[2] = -1
            
            for i in range(0, sequence_length):
                if sequences[0].seq[i:i+3] in stop_codons:
                    frame = i % 3
                    if found_stops[frame] is None:
                        found_stops[frame] = i
                if all(found_stops):
                    break
            else: # If no stop codon found fill in the gaps
                for frame, stop in enumerate(found_stops):
                    if not stop:
                        found_stops[frame] = -1

            overlapping_ORF = max(found_stops)
            if overlapping_ORF > circular_length:
                circular_length = overlapping_ORF

        if circular_length < 197:
            circular_length = 1000
        else:
            circular_length += 3

        output_circular_fasta_file_name = output_path + '.extended.fasta' 
        with open(output_circular_fasta_file_name, "w", encoding="utf-8" ) as output_handle:
            sequences[0].seq = sequences[0].seq[-circular_length:] + sequences[0].seq + sequences[0].seq[:circular_length]
            SeqIO.write([sequences[0]], output_handle, "fasta")
            input_file_name = output_circular_fasta_file_name


   # Extract annotation information from genbank files
    feature_tuple = namedtuple('feature', ['name', 'strand', 'start', 'end', 'color', 'nest_level', 'type'])
    features_list = []
    if input_file_type in ['genbank', 'gb', 'gbk', 'gbff']:

        features = [feature for feature in [rec for rec in SeqIO.parse(input_sequence, "genbank")][0].features]
        
        delayable_features = delay_iterator(features)  # Works like appending to a list you're looping through, but lets you know when you're getting delayed objects
        for feature, is_delayed in delayable_features:
            # Make note of the feature type
            feature_type = feature.type.lower()

            # Save ORFs for last
            translation_names = ['cds', 'translation', 'orf']
            if feature_type in translation_names and not is_delayed:
                delayable_features.delay(feature)
                continue

            # If the color is defined in the genbank file, use it
            if 'ApEinfo_revcolor' in feature.qualifiers:
                color = feature.qualifiers['ApEinfo_revcolor'][0]
            elif 'ApEinfo_fwdcolor' in feature.qualifiers:
                color = feature.qualifiers['ApEinfo_fwdcolor'][0]

            # Otherwise, we assign color based on the feature type
            elif feature_type in ["gene"]:
                color = COLORS['cds']
            elif feature_type in ['promoter']:
                color = COLORS['promoter']
            elif feature_type in ['terminator']:
                color = COLORS['terminator']
            elif feature_type in ['ori', 'rep_origin']:
                color = COLORS['ori']
            elif feature_type in ['ncRNA']:
                color = COLORS['ncrna']
            elif feature_type in ['primer_bind']:
                color = COLORS['primer_bind']
            elif feature_type in translation_names:
                color = COLORS['cds']
            else:
                color = COLORS['misc']

            feature_name = feature.qualifiers.get('label', [f'Unnamed {feature_type}'])

            if not feature_name:
                feature_name = f'Unnamed {feature_type}'
            else:
                feature_name = feature_name[0]

            # Filter out CDSs that match the start and end of a gene
            if feature_type in translation_names:
                for gene in features_list:
                    if gene.type == 'gene' and feature.location.start == gene.start and feature.location.end == gene.end:
                        break
                else:
                    feature_hit = feature_tuple(feature_name,
                                        int(feature.location.strand),
                                        feature.location.start,
                                        feature.location.end,
                                        color,
                                        None,
                                        feature_type)
                    features_list.append(feature_hit)
            else:
                # Add the processed feature to the list
                feature_hit = feature_tuple(feature_name,
                                            int(feature.location.strand),
                                            feature.location.start,
                                            feature.location.end,
                                            color,
                                            None,
                                            feature_type)
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
    if not name and sequences[0].name:
        # Get the name from the first contig in the genbank file
        name = sequences[0].name
    # Load sequence
    i = 0
    main_seq = None
    for this_seq in SeqIO.parse(input_file_name, "fasta"):
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
                                        score""")
    rhotermpredict_results = []
    for result in rhotermpredict_predictions:
        rhotermpredict_results.append(rdtresult(result.strand,
                                                result.c_over_g,
                                                result.start_rut-circular_length,
                                                result.end_rut-circular_length,
                                                result.term_seq,
                                                result.downstream_seq,
                                                result.palindromes,
                                                result.pause_concensus,
                                                result.score))


    logger.info('Running RIT_predict_TransTerm')
    transterm_predictions = transterm(input_file_name, circular_length)
    # ------------------------------------------------------------------------------
    # PERFORM PROMOTER PREDICTION
    # ------------------------------------------------------------------------------

    # Predict promoters
    logger.info('Running TSS_predict_promoter_calculator')
    promoter_calc_predictions = promocalc(main_seq, circular_length=circular_length, threads=threads)

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
    logger.info('Running RBS prediction')
    rbs_predictions = ostir(input_file_name, threads=threads)

    # ------------------------------------------------------------------------------

    # Identify expressed orfs
    expressed_orfs = []
    expressed_orf = namedtuple("orf_result", "start, end, expression, burden, dG, array, start_codon, strand")
    orfs = orf_predictions
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
            processed_orf = expressed_orf(int(orf_info['start'])-circular_length,
                                            int(orf_info['end'])-circular_length,
                                            rbs_score,  # Log Version
                                            orf_length * rbs_info.score,  # Burden
                                            rbs_info.score2,
                                            array,
                                            rbs_info.start_codon,
                                            orf_info["strand"]
                                            )
            expressed_orfs.append(processed_orf)

            total_burden += orf_length * rbs_info.score

    #   When available, assign RBS info to orf
            orf['array'] = orf_array
            orf['array_minus'] = array_minus
            orf['rbs_score'] = rbs_score

        else:
            rbs_score = 0
            orf['array'] = orf_array
            orf['array_minus'] = array_minus
            orf['rbs_score'] = rbs_score


    # Filter predictions that we show
    expressed_orfs = list(filter(lambda x:float(x.expression) > rbs_score_cutoff, expressed_orfs))

    # ---- CIRCULAR SEQUENCE ----

    if circular:


        # Remove ones with out of bounds starts/ends/positions (whichever is graphed)

        expressed_orfs = list(filter(lambda x: (0 <= int(x.start) < len(single_sequence)), expressed_orfs))

        promoter_calc_predictions = list(filter(lambda x: 0< int(x.TSSpos) < len(single_sequence), promoter_calc_predictions))

        transterm_predictions = list(filter(lambda x: 0< int(x.start) < len(single_sequence), transterm_predictions))

        rhotermpredict_results = list(filter(lambda x: 0< int(x.start_rut) < len(single_sequence), rhotermpredict_results))

    # Set up the results object
    result = CryptResults(name = name,
                          sequence = str(single_sequence),
                          translation_sites = expressed_orfs,
                          rho_dep_terminators = rhotermpredict_results,
                          rho_ind_terminators = transterm_predictions,
                          promoters =  promoter_calc_predictions,
                          annotations = features_list,
                          burden = total_burden)
    return result


if __name__ == "__main__":
    main()
