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
from collections import namedtuple
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
        outdir = tempfile.TemporaryDirectory()
        output_path = outdir.name


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
    single_sequence = deepcopy(sequences[0].seq)  # We need to keep ahold of the original sequence to revert from circularizing
    input_file_name = infile

    # --- CIRCULAR SEQUENCE ---
    # create a 3x version of the FASTA input file


    if len(sequences) > 1:
        logger.warning("Input file contains multiple sequences. Only the first will be processed.")

    
    sequence_length = len(sequences[0])


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
        for start_codons, stop_codons in [(fwd_start_codons, fwd_stop_codons), (rev_start_codons, rev_stop_codons)]:
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

            overlapping_ORF = max(found_starts)
            if overlapping_ORF > circular_length:
                circular_length = overlapping_ORF

        if circular_length < 197:
            circular_length = 200

        output_circular_fasta_file_name = output_path + '.extended.fasta'
        with open(output_circular_fasta_file_name, "w", encoding="utf-8" ) as output_handle:
            sequences[0].seq = sequences[0].seq[-circular_length:] + sequences[0].seq + sequences[0].seq[:circular_length]
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
            if not feature.qualifiers.get('label'): # Some features have no label (ex. transcripts from Benchling)
                            continue
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


    # Filter predictions that we show
    expressed_orfs = list(filter(lambda x:float(x.expression) > rbs_score_cutoff, expressed_orfs))

    # ---- CIRCULAR SEQUENCE ----

    if circular:

        # Remove ones with out of bounds starts/ends/positions (whichever is graphed)
        orf_predictions = list(filter(lambda x: (int(x['end']) >= circular_length + 1) and
                                                (int(x['start']) <= circular_length + len(single_sequence)), orf_predictions))

        rbs_predictions = list(filter(lambda x: (int(x.position) >= circular_length + 1) and
                                                (int(x.position) <= circular_length + len(single_sequence)), rbs_predictions))

        expressed_orfs = list(filter(lambda x: (int(x.end) >= circular_length + 1) and
                                                (int(x.start) <= circular_length + len(single_sequence)), expressed_orfs))
        
        promoter_calc_predictions = list(filter(lambda x: (int(x.TSSpos) >= circular_length + 1) and
                                                (int(x.TSSpos) <= circular_length + len(single_sequence)), promoter_calc_predictions))
        transterm_predictions = list(filter(lambda x: (int(x.end) >= circular_length + 1) and
                                                (int(x.end) <= circular_length + len(single_sequence)), transterm_predictions))
        rhotermpredict_results = list(filter(lambda x: (int(x.end_rut) >= circular_length + 1) and
                                                (int(x.end_rut) <= circular_length + len(single_sequence)), rhotermpredict_results))
        # Correct coordinates of remaining

        for i, rbs in enumerate(rbs_predictions):
            rbs_hit = namedtuple("rbs_hit", "position, start_codon, strand, score, score2")
            rbs = rbs_hit(str(int(rbs.position) - circular_length),
                            rbs.start_codon,
                            rbs.strand,
                            rbs.score,
                            rbs.score2)
            rbs_predictions[i] = rbs

        for i, tss in enumerate(promoter_calc_predictions):
            promoter_calc_result = namedtuple("promoter_calc_result",
                                                "seq, score, strand, TSSpos, box35pos, box35seq, box10pos, box10seq")
            tss = promoter_calc_result(tss.seq,
                                        tss.score,
                                        tss.strand,
                                        int(tss.TSSpos) - circular_length,
                                        [int(tss.box35pos[0]) - circular_length, int(tss.box35pos[1]) - circular_length],
                                        tss.box35seq,
                                        [int(tss.box10pos[0]) - circular_length, int(tss.box10pos[1]) - circular_length],
                                        tss.box10seq)
            promoter_calc_predictions[i] = tss

        for i, rit in enumerate(transterm_predictions):
            transterm_result = namedtuple("transterm_result", 
                                            "start, end, strand, conf, hairpin_score, tail_score, seq_upstream, seq_hairpin_open, seq_hairpin_close, seq_tail, seq_hairpin_loop")
            rit = transterm_result(int(rit.start) - circular_length,
                                    int(rit.end) - circular_length,
                                    rit.strand,
                                    rit.conf,
                                    rit.hairpin_score,
                                    rit.tail_score,
                                    rit.seq_upstream,
                                    rit.seq_hairpin_open,
                                    rit.seq_hairpin_close,
                                    rit.seq_tail,
                                    rit.seq_hairpin_loop)
            transterm_predictions[i] = rit

        for i, rdt in enumerate(rhotermpredict_results):
            rhotermpredict_result = namedtuple("rhotermpredict_result",
                                                "strand, c_over_g, start_rut, end_rut, term_seq, palindromes, pause_concensus, scores")
            rdt = rhotermpredict_result(rdt.strand,
                                        rdt.c_over_g,
                                        int(rdt.start_rut) - circular_length,
                                        int(rdt.end_rut) - circular_length,
                                        rdt.term_seq,
                                        rdt.palindromes,
                                        rdt.pause_concensus,
                                        rdt.scores)
            rhotermpredict_results[i] = rdt

        for i, orf in enumerate(expressed_orfs):
            orf_object = namedtuple("orf_result", "start, end, expression, burden, dG, array, start_codon, strand")
            orf = orf_object(orf.start - circular_length,
                        orf.end - circular_length,
                        orf.expression,
                        orf.burden,
                        orf.dG,
                        orf.array,
                        orf.start_codon,
                        orf.strand)
            expressed_orfs[i] = orf

        # Split reading frames and assign
        added_split_orfs = []
        for orf in orf_predictions:
            orf['start'] = int(orf['start']) - circular_length
            orf['end'] = int(orf['end']) - circular_length

            if int(orf['end']) > circular_length:

                new_split_orf = deepcopy(orf)
                orf['end'] = circular_length
                new_split_orf['start'] = "1"
                new_split_orf['end'] = int(new_split_orf['end']) - circular_length

                added_split_orfs.append(new_split_orf)
        orf_predictions.extend(added_split_orfs)


    # Set up the results object

    result = CryptResults(name = name,
                          sequence = single_sequence,
                          translation_sites = expressed_orfs,
                          row_dep_terminators = rhotermpredict_results,
                          row_ind_terminators = transterm_predictions,
                          promoters =  promoter_calc_predictions,
                          annotations = features_list,
                          burden = total_burden)

    print(f"Found promoters: {len(promoter_calc_predictions)}")
    print(f"Found row-independent terminators: {len(transterm_predictions)}")
    print(f"Found row-dependent terminators: {len(rhotermpredict_results)}")
    print(f"Found expressed ORFs: {len(expressed_orfs)}")


    return result


if __name__ == "__main__":
    main()
