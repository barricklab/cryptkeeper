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
import plotly
import plotly.graph_objs as go
import numpy

from Bio import SeqIO
from operator import itemgetter
from copy import deepcopy
from collections import namedtuple

from .TTS_predict_promoter_calculator import tts_predict
from rhotermpredict import rho_term_predict
from .RIT_predict_TransTerm import RIT_predict_TransTerm
from .ORF_predict import ORF_predict
from .RBS_calc_via_OSTIR import RBS_predict_RBS_Calculator
from .ORF_predict import find_orfs


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

    parser.add_argument(
        '-w', '--web',
        action='store_true',
        dest='web',
        help="web mode. Writes output as a div without including javascript to file *.div instead of normal *.html",
    )

    # parser.add_argument(
    #    '--pdf',
    #    action='store_true',
    #    dest='pdf',
    #    help="write output plot as *.pdf",
    #    )

    # ------------------------------------------------------------------------------
    options = parser.parse_args()

    cryptkeeper(input_file = options.i,
                output=options.o,
                circular=options.circular,
                plot_only=options.plot_only,
                name=options.name,
                rbs_score_cutoff=options.rbs_score_cutoff,
                web=options.web)


def cryptkeeper(input_file, output=None, circular=False, plot_only=False, name=None, rbs_score_cutoff=2.0, web=False):

    # @TODO: Many of these steps should be moved to a function

    # ------------------------------------------------------------------------------
    # PROCESS INPUT FILES
    # ------------------------------------------------------------------------------

    # Makes relative file arguments absolute for more robust execution

    input_sequence = os.path.abspath(input_file)
    output_path = os.path.abspath(output)
    output_folder = "/".join(output_path.split("/")[:-1])
    input_gene_name = input_sequence.split("/")[-1]
    input_file_type = input_gene_name.split(".")[-1]
    input_gene_name = ".".join(input_gene_name.split(".")[:-1])
    input_file_name = input_sequence

    # Set up log
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    stream_formatter = logging.Formatter('[%(asctime)s] Cryptkeeper: %(message)s')
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(stream_formatter)
    logger.addHandler(stream_handler)

    file_formatter = logging.Formatter('[%(asctime)s - %(levelname)s] Cryptkeeper: %(message)s')
    file_handler = logging.FileHandler(filename=f"{output_path}.log")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)


    # Convert non-FASTA files to FASTA (supports Genbank)
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
    features_dict = {}
    if input_file_type in ['genbank', 'gb', 'gbk']:

        features = [feature for feature in [rec for rec in SeqIO.parse(input_sequence, "genbank")][0].features]
        for i, feature in enumerate(features):
            if 'ApEinfo_revcolor' in feature.qualifiers:
                color = feature.qualifiers['ApEinfo_revcolor'][0]
            else:
                color = "black"

            features_dict[i] = {
                "name": feature.qualifiers['label'][0],
                "strand": feature.location.strand,
                "start": feature.location.start,
                "end": feature.location.end,
                "color": color,
                'nest_level': None
            }


        # Process nest level for displaying features
        strandless_features = [feature for feature in features_dict if features_dict[feature]['strand'] not in [1, -1]]
        forward_features = [feature for feature in features_dict if features_dict[feature]['strand'] == 1]
        reverse_features = [feature for feature in features_dict if features_dict[feature]['strand'] == -1]
        end_positions = dict()
        #for strand in [strandless_features, forward_features, reverse_features]:  Dont care about organizing this
        for strand in [features_dict]:
            for feature in strand:
                if not end_positions:
                    end_positions[str(0)] = features_dict[feature]['end']
                    features_dict[feature]['nest_level'] = 0
                    continue
                for level, position in end_positions.items():
                    level = int(level)
                    if position == 'inf':
                        continue
                    elif int(position) > int(features_dict[feature]['start']):
                        continue
                    else:
                        end_positions[str(level)] = features_dict[feature]['end']
                        features_dict[feature]['nest_level'] = level
                        break
                if not features_dict[feature]['nest_level'] and features_dict[feature]['nest_level'] != 0:  # Layer 0 is falcey
                    lowest_level = max([int(n) for n in end_positions.keys()])
                    end_positions[str(lowest_level+1)] = features_dict[feature]['end']
                    features_dict[feature]['nest_level'] = lowest_level+1
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


    # @TODO: These should directly report data instead of going through an intermediate file

    # ------------------------------------------------------------------------------
    # PERFORM TERMINATOR PREDICTIONS
    # ------------------------------------------------------------------------------

    # Predict terminators
    rit_prediction_file_name = output_path + '.RIT.csv'


    RhoTermPredict_out = output_path+'_rtp.csv'
    if not plot_only:
        logger.info('Running RhoTermPredict')
        rho_term_predict(inseq=[str(main_seq.seq)], csv_out = RhoTermPredict_out)  # Run RhoTermPredict
        # Load RhoTermPredict output into named tuple
    RhoTermPredict_out += f"RhoTermPredict_{seq_name}.csv"

    if not plot_only:
        logger.info('Running RIT_predict_TransTerm')
        RIT_predict_TransTerm(input_file_name, rit_prediction_file_name)

    rit_transterm_predictions = []  # Transterm is the only dependency that requires file IO.
    if os.path.isfile(rit_prediction_file_name):
        rit_reader = csv.DictReader(open(rit_prediction_file_name, encoding="utf-8" ))
        for row in rit_reader:
            rit_transterm_predictions.append(row)


    rit_rhotermpredict_predictions = []
    if os.path.isfile(RhoTermPredict_out):
        rit_reader = csv.DictReader(open(RhoTermPredict_out, encoding="utf-8" ))
        for row in rit_reader:
            rit_rhotermpredict_predictions.append(row)

    # ------------------------------------------------------------------------------
    # PERFORM PROMOTER PREDICTION
    # ------------------------------------------------------------------------------

    # Predict promoters
    tss_prediction_file_name = output_path + '.TSS.csv'

    if not plot_only:
        #logger.info('Running TSS_predict_BPROM')
        #TTS_predict_BPROM(input_file_name, tss_prediction_file_name)
        logger.info('Running TSS_predict_promoter_calculator')
        tss_predictions = tts_predict(main_seq, tss_prediction_file_name)

    # ------------------------------------------------------------------------------
    # PERFORM ORF PREDICTION / GENERATE ANNOTATION FILE
    # ------------------------------------------------------------------------------

    # Predict ORFs

    orf_prediction_file_name = output_path + '.ORF.csv'

    if not plot_only:
        logger.info('Running ORF prediction')
        ORF_predict(input_file_name, orf_prediction_file_name)

    orf_predictions = []
    if os.path.isfile(orf_prediction_file_name):
        orf_reader = csv.DictReader(open(orf_prediction_file_name, encoding="utf-8" ))
        for row in orf_reader:
            row['start_codon_position'] = int(row['start'] if row['strand'] == '+' else row['end'])
            orf_predictions.append(row)

    # Re-sort by the start codon position so we are in the same order as RBS
    orf_predictions = sorted(orf_predictions, key=itemgetter('start_codon_position'))

    # Create a gene annotation file
    gene_annotation_file_name = output_path + '.gene.csv'

    # @todo: add gene annotation file

    # ------------------------------------------------------------------------------
    # PERFORM RBS PREDICTION
    # ------------------------------------------------------------------------------

    # Predict RBS
    rbs_prediction_file_name = output_path + '.RBS.csv'

    if not plot_only:
        logger.info("Running RBS prediction")
        RBS_predict_RBS_Calculator(input_file_name, rbs_prediction_file_name, starts=orf_prediction_file_name)

    # Combine ORF and RBS predictions AND calculate the crypt score
    rbs_predictions = []
    rbs_reader = csv.DictReader(open(rbs_prediction_file_name, encoding="utf-8" ))

    # Find ORFs with predicted RBS and append their scores. @croots: Bad things happen if you throw the whole seq at OSTIR
    for row in rbs_reader:
        row["score"] = numpy.log10(float(row["score"]))
        rbs_predictions.append(row)

    total_burden = 0


    # ------------------------------------------------------------------------------


    summary_file_name = output_path + '.summary.txt'
    summary_file = open(summary_file_name, 'w', encoding="utf-8" )


    # Identify expressed ORFs
    expressed_ORFs = []
    ORFs = find_orfs(input_file_name, 11, 0)
    start_codons_fwd = {}   # Value is index in rbs_predictions. Havent checked, probably faster than looping dicts.
    start_codons_rev = {}

    for i, prediction in enumerate(rbs_predictions):
        if prediction['strand'] == '+':
            start_codons_fwd[prediction['position']] = i
        else:
            start_codons_rev[prediction['position']] = i

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

        if str(ORF['start']+offset) in expressed_starts.keys():
            RBS_info = rbs_predictions[expressed_starts[str(ORF['start']+offset)]]
            rbs_score = RBS_info['score']
            if str(rbs_score) in ['inf', '-inf']:
                continue

            if RBS_info["start_codon"] != ORF_info["start_codon"]:
                logger.debug(f'Codon mismatch {RBS_info["position"]} {RBS_info["start_codon"]} and {ORF_info["start"]} {ORF_info["start_codon"]}, {ORF_info["strand"]}\n{ORF_info["translation"]}')

            expressed_ORFs.append({
                "start": int(ORF_info['start']),
                "end": int(ORF_info['end']),
                "score": rbs_score,
                "burden": orf_length * rbs_score,
                "array": ORF_array,
                "array_minus": array_minus,
                "second_half": "T",
                "start_codon": RBS_info["start_codon"],
                "strand": ORF_info["strand"]
            })
    #   When available, assign RBS info to ORF
            ORF['array'] = ORF_array
            ORF['array_minus'] = array_minus
            ORF['rbs_score'] = 10**RBS_info['score']
            ORF['second_half'] = "T"

        else:
            rbs_score = 0
            ORF['array'] = ORF_array
            ORF['array_minus'] = array_minus
            ORF['rbs_score'] = rbs_score
            ORF['second_half'] = "T"

        add_to_score = True
        if circular:  # !!! BUG: IF TRUE THIS WILL ERROR @TODO
            if (int(ORF_info['start']) <= sequence_length) or (int(ORF_info['end']) >= sequence_length * 2 + 1):
                add_to_score = False

        if add_to_score:
            total_burden += rbs_score * orf_length
            summary_file.write(f"  {str(ORF_info['start'])} Burden: {str(orf_length * rbs_score)} [ RBS Score: " +
                               f"{str(rbs_score)} ORF Length: {str(orf_length)} ]\n")

    summary_file.write("Total Burden: " + str(total_burden))

    # Filter predictions that we show
    expressed_ORFs = list(filter(lambda x:float(x["score"]) > rbs_score_cutoff, expressed_ORFs))

    #reload the gene annotation files
    gene_annotations = []
    if os.path.isfile(gene_annotation_file_name):
        gene_reader = csv.DictReader(open(gene_annotation_file_name, encoding="utf-8" ))
        for row in gene_reader:
            row['position'] = row['end'] if row['strand'] == '+' else row['start']
            row['array'] = 0 if row['strand'] == '+' else int(row['end']) - int(row['start'])
            row['array_minus'] = int(row['end']) - int(row['start']) if row['strand'] == '+' else 0
            gene_annotations.append(row)

    # ---- CIRCULAR SEQUENCE ----

    if circular:

        # Remove ones with out of bounds starts/ends/positions (whichever is graphed)
        orf_predictions = list(filter(lambda x: (int(x['start']) >= sequence_length + 1) and
                                                (int(x['start']) <= sequence_length * 2), orf_predictions))
        rbs_predictions = list(filter(lambda x: (int(x['position']) >= sequence_length + 1) and
                                                (int(x['position']) <= sequence_length * 2), rbs_predictions))
        tss_predictions = list(filter(lambda x: (int(x['TSSpos']) >= sequence_length + 1) and
                                                (int(x['TSSpos']) <= sequence_length * 2), tss_predictions))
        rit_transterm_predictions = list(filter(lambda x: (int(x['end']) >= sequence_length + 1) and
                                                (int(x['end']) <= sequence_length * 2), rit_transterm_predictions))

        # Correct coordinates of remaining

        for d in rbs_predictions:
            d['position'] = str(int(d['position']) - sequence_length)

        for d in tss_predictions:
            d['TSSpos'] = str(int(d['TSSpos']) - sequence_length)

        for d in rit_transterm_predictions:
            d['end'] = str(int(d['end']) - sequence_length)

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

                if d['strand'] != '+':
                    d['second_half'] = 'T'
                else:
                    new_split_orf['second_half'] = 'T'

                added_split_orfs.append(new_split_orf)
        orf_predictions.extend(added_split_orfs)

    # logger.normal out entire changed lists - for graphing outside of this program in R, etc.
    with open(output_path + ".final.ORF.csv",'w', encoding="utf-8" ) as final_orf_predictions_file:
        writer = csv.DictWriter(
            final_orf_predictions_file,
            fieldnames=["start", "end", "strand", "start_codon", 'length', 'rbs_score',
                        'array', 'array_minus', 'start_codon_position', 'second_half', 'translation']
            )
        writer.writeheader()
        writer.writerows(orf_predictions)
        final_orf_predictions_file.close()

    with open(output_path + ".final.RBS.csv", 'w', encoding="utf-8" ) as final_rbs_predictions_file:
        writer = csv.DictWriter(
            final_rbs_predictions_file,
            fieldnames=["position", "strand", "start_codon", "score", "score2", "burden", "array", "array_minus"]
            )
        writer.writeheader()
        writer.writerows(rbs_predictions)
    final_rbs_predictions_file.close()

    # Display predictions

    #settings
    TSS_color = '#2ab717'
    RIT_color = '#ff0000'
    RBS_color = '#0d30e0'
    marker_size = 20
    tss_score_max = 8000


    tss_series = go.Scatter(
      x = [d.TSSpos for d in tss_predictions],
      y = [d.score for d in tss_predictions],
      mode = 'markers',
      name = 'TSS',
      text = [ (f'Score: {d.score}<br>-35: {d.box35seq}<br>-10: {d.box10seq}')  for d in tss_predictions],  # TODO Change this up for Rho dependant?
      marker = dict(
        symbol=[ ('triangle-right'  if d.strand == '+' else 'triangle-left')  for d in tss_predictions],
        size=marker_size,
        color=TSS_color,
      ),

      error_y=dict(
                type='data',
                array=[0],
                arrayminus=[d.score for d in tss_predictions],
                width=0,
                color=TSS_color,
            )
    )

    rit_transterm_series = go.Scatter(
      x = [int(e['end']) for e in rit_transterm_predictions],
      y = [int(e['conf']) for e in rit_transterm_predictions],
      mode = 'markers',
      name = 'RIT',
      yaxis='y2',
      text = [ ('hairpin score:' + d['hairpin_score'] + '; tail score:' + d['tail_score'] + '<br>' +
                d['seq_hairpin_open'] + '-' + d['seq_hairpin_loop'] + '-' + d['seq_hairpin_close'] +
                '-' + d['seq_tail'][0:7] )  for d in rit_transterm_predictions],
      marker = dict(
        symbol=[ ('triangle-right'  if d['strand'] == '+' else 'triangle-left')  for d in rit_transterm_predictions],
        size=marker_size,
        color=RIT_color,
      ),
      error_y=dict(
              type='data',
              array=[0],
              arrayminus=[e['conf'] for e in rit_transterm_predictions],
              width=0,
              color=RIT_color,
          )
    )


    rdt_rhotermpredict_series = go.Scatter(
      x = [int(e['end_rut']) for e in rit_rhotermpredict_predictions],
      y = [int(e['score_sum']) for e in rit_rhotermpredict_predictions],
      mode = 'markers',
      name = 'RDT',
      yaxis='y6',
      text = [ ('placeholder')  for d in rit_rhotermpredict_predictions],
      marker = dict(
        symbol=[ ('triangle-right'  if d['strand'] == '+' else 'triangle-left')  for d in rit_rhotermpredict_predictions],
        size=marker_size,
        color=RIT_color,
      ),
      error_y=dict(
              type='data',
              array=[0],
              arrayminus=[e['score_sum'] for e in rit_rhotermpredict_predictions],
              width=0,
              color=RIT_color,
          )
    )

    negatives_exist = True if "-" in [f['strand'] for f in expressed_ORFs] else False
    rbs_series = go.Scatter(
      x = [f['start'] for f in expressed_ORFs],
      y = [f['score'] for f in expressed_ORFs],
      mode = 'markers',
      text = [ ('start codon:' + d['start_codon'] + "<br>Burden ×10<sup>6</sup>: " +
                "{0:.6f}".format(d['burden']/1000000) + " (" + "{0:.2f}".format(100*d['burden']/total_burden) + "%}")
               for d in expressed_ORFs],
      name = 'RBS',
      yaxis='y3',
      marker = dict(
        symbol=[ ('triangle-right'  if d['strand'] == '+' else 'triangle-left')  for d in expressed_ORFs],
        size=marker_size,
        color=RBS_color,
      ),
      error_y=dict(
              type='data',
              array=[0],
              arrayminus=[f['score'] for f in expressed_ORFs],
              width=0,
              color=RBS_color,
      ),
      error_x=dict(
              type='data',
              array=[f['array'] for f in expressed_ORFs],
              arrayminus=[f['array_minus'] for f in expressed_ORFs],
              width=3,
              color=RBS_color,
      ),
    )

    nucleotide_sequence  = go.Scatter(
      x=list(range(0,len(main_seq)-1)),
      y=[0.97]*len(main_seq),
      mode='text',
      name = 'Sequence',
      text=list(main_seq.seq),
      yaxis='y4',
    )

    # Draw features
    minimum_annotation_length = 10
    base_y_level = 0.9
    feature_annotations = []
    feature_annotations.append(go.Scatter(x=[-1], y=[-1], fill="toself", yaxis='y5',  # Placeholder
                                          text=None,
                                          name="Annotations",
                                          mode='none',
                                          marker=dict(
                                              line=dict(width=0.5)),
                                          hoveron='fills',
                                          legendgroup='1',
                                          fillcolor="#000000"
                                          ))

    for feature, orf in features_dict.items():  # Draw features on Cryptkeeper Results
        start_pos_x = orf['start']-1
        end_pos_x = orf['end']-1
        y_pos = base_y_level-orf['nest_level']*0.025

        if end_pos_x - start_pos_x < minimum_annotation_length:
            continue

        # Draw box and arrow
        arrow_offset = 100
        '''
        if orf['strand'] == 1:
            x_positions = [start_pos_x, end_pos_x - arrow_offset, end_pos_x - arrow_offset, start_pos_x, start_pos_x, None,
                       end_pos_x, end_pos_x - arrow_offset, end_pos_x - arrow_offset, end_pos_x]
            y_positions = [round(y_pos + 0.005, 4), round(y_pos + 0.005, 4), round(y_pos - 0.005, 4),
                           round(y_pos - 0.005, 4), round(y_pos + 0.005, 4), None,
                           y_pos, round(y_pos + 0.01, 4), round(y_pos - 0.01, 4), y_pos]
        elif orf['strand'] == -1:
            x_positions = [start_pos_x + arrow_offset, end_pos_x, end_pos_x, start_pos_x + arrow_offset, start_pos_x + arrow_offset, None,
                       start_pos_x, start_pos_x + arrow_offset, start_pos_x + arrow_offset, start_pos_x]
            y_positions = [round(y_pos + 0.005, 4), round(y_pos + 0.005, 4), round(y_pos - 0.005, 4),
                           round(y_pos - 0.005, 4), round(y_pos + 0.005, 4), None,
                           y_pos, round(y_pos + 0.01, 4), round(y_pos - 0.01, 4), y_pos]
        '''
        x_positions = [start_pos_x + arrow_offset, end_pos_x, end_pos_x, start_pos_x + arrow_offset,
                        start_pos_x + arrow_offset]
        y_positions = [round(y_pos + 0.005, 4), round(y_pos + 0.005, 4), round(y_pos - 0.005, 4),
                        round(y_pos - 0.005, 4), round(y_pos + 0.005, 4)]


        feature_annotations.append(go.Scatter(x=x_positions, y=y_positions, fill="toself", yaxis='y5',
                                              text=orf['name'],
                                              name="Annotations",
                                              mode='none',
                                              marker=dict(
                                                  line=dict(width=0.5)),
                                              hoveron='fills',
                                              legendgroup='1',
                                              showlegend=False,
                                              fillcolor=features_dict[feature]['color']
                                              ))

    ## Use display name at command line if it was given, otherwise file name
    if name:
        display_name = name
    else:
        display_name = input_file_name

    layout = go.Layout(
      title='CryptKeeper Results: ' + display_name + "  Burden ×10<sup>6</sup>: " + "{0:.2f}".format(total_burden/1000000),
      hovermode = 'closest',
      autosize=True,

      xaxis=dict(
            domain=[0.20, 1.0]
      ),
      yaxis=dict(
        title='TSS score',
        range = [0, tss_score_max+0.2],
        titlefont=dict(
            color=TSS_color
        ),
        tickfont=dict(
            color=TSS_color
        )
      ),
      yaxis2=dict(
          title='RIT score',
          range = [0, 200],
          titlefont=dict(
              color=RIT_color
          ),
          tickfont=dict(
              color=RIT_color
          ),
          anchor='free',
          overlaying='y',
          side='left',
          position=0.10,
      ),
      yaxis3=dict(
          title='log10 RBS score',
          range = [2, 6],   # TODO: Dynamicly set this
          titlefont=dict(
              color=RBS_color
          ),
          tickfont=dict(
              color=RBS_color
          ),
          anchor='free',
          overlaying='y',
          side='left',
          position=0.15,
      ),
      yaxis4=dict(
          title='sequence',
          range = [0, 1],
          titlefont=dict(
              color='black'
          ),
          tickfont=dict(
              color='black'
          ),
          anchor='free',
          overlaying='y',
          side='left',
          position=0.05,
          visible=False,
          showticklabels=False
      ),
        yaxis5=dict(
            title='Annotations',
            range=[0, 1],
            titlefont=dict(
                color='black'
            ),
            tickfont=dict(
                color='black'
            ),
            anchor='free',
            side='left',
            position=0.05,
            visible=False,
            showticklabels=False,
            overlaying='y',
        ),
        yaxis6=dict(
            title='RDT Score',
            range = [0, 300],
          titlefont=dict(
              color=RIT_color
          ),
          tickfont=dict(
              color=RIT_color
          ),
          anchor='free',
          overlaying='y',
          side='left',
          position=0.05,
        ),
        
    )
    data = [tss_series, rit_transterm_series, rbs_series, nucleotide_sequence, rdt_rhotermpredict_series]
    data.extend(feature_annotations)

    fig = dict( data=data, layout=layout )

    #if (pdf):
    #  fig.write_image(output_path+'.pdf')
    if (web):
        dev_content = plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')
        f = open(output_path+'.plot.div', 'w')
        f.write(dev_content)
        f.close()
    else:
        plotly.offline.plot(fig, filename=output_path+'.plot.html')


    return expressed_ORFs, tss_predictions, rit_rhotermpredict_predictions, rit_transterm_predictions


if __name__ == "__main__":
    main()
