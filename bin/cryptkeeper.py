#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""

@authors: Jeffrey Barrick

Predict cryptic bacterial gene expression signals in an input sequence.

"""

import argparse
from Bio import SeqIO
import csv
from operator import itemgetter
import os
from copy import deepcopy

import plotly
import plotly.graph_objs as go
import numpy
from helpers import Logger

logger = Logger()

def main(options):
    # Makes relative file arguments absolute for more robust execution
    options.i = os.path.abspath(options.i)
    options.o = os.path.abspath(options.o)
    output_folder = "/".join(options.o.split("/")[:-1])
    input_gene_name = options.i.split("/")[-1]
    input_file_type = input_gene_name.split(".")[-1]
    input_gene_name = ".".join(input_gene_name.split(".")[:-1])
    input_file_name = options.i

    # Convert non-FASTA files to FASTA (supports Genbank)
    if input_file_type != "fasta" and input_file_type != "fna":
        logger.normal("Non-FASTA file detected. Converting.")
        if input_file_name == "genbank" or input_file_name == "gb":
            with open(input_file_name, "r") as file_in:
                with open(f"{output_folder + input_gene_name}.fna", "w") as file_converted:
                    sequences = SeqIO.parse(file_in, "genbank")
                    sequences = SeqIO.write(sequences, file_converted, "fasta")
        else:
            raise ValueError("Unable to convert file. Supported types: FASTA, GENBANK")
    else:
        sequences = list(SeqIO.parse(options.i, "fasta"))

    # --- CIRCULAR SEQUENCE ---
    # create a 3x version of the FASTA input file

    output_circular_fasta_file_name = options.o + '.circular.3x.fasta'
    if len(sequences) > 1:
        logger.normal("Input file contains multiple sequences. Only the first will be processed.")
    sequence_length = len(sequences[0])
    if (options.circular):
        with open(output_circular_fasta_file_name, "w") as output_handle:
            sequences[0].seq = sequences[0].seq + sequences[0].seq + sequences[0].seq
            SeqIO.write([sequences[0]], output_handle, "fasta")
            input_file_name = output_circular_fasta_file_name
    '''
    i = 0
    for this_seq in SeqIO.parse(options.i, "fasta"):
        i += 1
        if i > 1:
            logger.normal("FASTA contains multiple sequences. Only the first will be processed.")
            break
        sequence_length = len(this_seq)
        if (options.circular):
            with open(output_circular_fasta_file_name, "w") as output_handle:
                this_seq.seq = this_seq.seq + this_seq.seq + this_seq.seq
                SeqIO.write([this_seq], output_handle, "fasta")
                input_file_name = output_circular_fasta_file_name
    '''


    # Predict promoters
    tss_prediction_file_name = options.o + '.TSS.csv'

    if not options.plot_only:
        logger.normal('Running TTS_predict_BPROM')
        from TSS_predict_BPROM import TTS_predict_BPROM
        TTS_predict_BPROM(input_file_name, tss_prediction_file_name)
        '''
        tss_command = f'python3 TSS_predict_BPROM.py ' \
                      f'-i {input_file_name} -o {tss_prediction_file_name}'
        logger.normal(tss_command)
        try:
            subprocess.check_call([tss_command], shell=True)
        except subprocess.CalledProcessError as e:
            logger.normal(e)
        except OSError as e:
            logger.normal(e)
        '''


    # Create a gene annotation file
    gene_annotation_file_name = options.o + '.gene.csv'

    # #=> to be implemented to read if input is genbank or GFF

    # Predict terminators
    rit_prediction_file_name = options.o + '.RIT.csv'

    if not options.plot_only:
        logger.normal('Running RIT_predict_TransTerm')
        from RIT_predict_TransTerm import RIT_predict_TransTerm
        RIT_predict_TransTerm(input_file_name, rit_prediction_file_name)
        '''
        rit_command = 'python3 RIT_predict_TransTerm.py -i ' + input_file_name + ' -o ' + rit_prediction_file_name
        logger.normal(rit_command)
        try:
            subprocess.check_call([rit_command], shell=True)
        except subprocess.CalledProcessError as e:
            logger.normal(e)
        '''

    # Predict ORFs

    orf_prediction_file_name = options.o + '.ORF.csv'

    if not options.plot_only:
        logger.normal('Running ORF prediction')
        from ORF_predict import ORF_predict
        ORF_predict(input_file_name, orf_prediction_file_name)
        '''
        orf_command = 'python3 ORF_predict.py -i ' + input_file_name + ' -o ' + orf_prediction_file_name
        logger.normal(orf_command)
        try:
            subprocess.check_call(orf_command, shell=True)
        except subprocess.CalledProcessError as e:
            logger.normal(e)
        except OSError as e:
            logger.normal(e)
        '''

    # Predict RBS
    rbs_prediction_file_name = options.o + '.RBS.csv'

    if not options.plot_only:
        logger.normal("Running RBS prediction")
        from RBS_calc_via_OSTIR import RBS_predict_RBS_Calculator
        RBS_predict_RBS_Calculator(input_file_name, rbs_prediction_file_name, starts=orf_prediction_file_name)


    # Load sequence
    i = 0
    main_seq = None
    for this_seq in SeqIO.parse(options.i, "fasta"):
        i += 1
        if i > 1:
            exit()
        main_seq = this_seq.upper()

    # Load predictions
    tss_predictions = []
    if os.path.isfile(tss_prediction_file_name):
        tss_reader = csv.DictReader(open(tss_prediction_file_name))
        for row in tss_reader:
            tss_predictions.append(row)

    rit_predictions = []
    if os.path.isfile(rit_prediction_file_name):
        rit_reader = csv.DictReader(open(rit_prediction_file_name))
        for row in rit_reader:
            rit_predictions.append(row)

    orf_predictions = []
    if os.path.isfile(orf_prediction_file_name):
        orf_reader = csv.DictReader(open(orf_prediction_file_name))
        for row in orf_reader:
            row['start_codon_position'] = int(row['start'] if row['strand'] == '+' else row['end'])
            orf_predictions.append(row)

    # Re-sort by the start codon position so we are in the same order as RBS
    orf_predictions = sorted(orf_predictions, key=itemgetter('start_codon_position'))
    # Combine ORF and RBS predictions AND calculate the crypt score
    rbs_predictions = []
    rbs_reader = csv.DictReader(open(rbs_prediction_file_name))

    # Find ORFs with predicted RBS and append their scores. @croots: Bad things happen if you throw the whole seq at OSTIR
    for row in rbs_reader:
        row["score"] = numpy.log10(float(row["score"]))
        rbs_predictions.append(row)

    total_burden = 0

    summary_file_name = options.o + '.summary.txt'
    summary_file = open(summary_file_name, 'w')


    # Identify expressed ORFs
    expressed_ORFs = []
    from ORF_predict import find_orfs
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
                print(f'Codon mismatch {RBS_info["position"]} {RBS_info["start_codon"]} and {ORF_info["start"]} {ORF_info["start_codon"]}, {ORF_info["strand"]}\n{ORF_info["translation"]}')

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
        if options.circular:  # !!! BUG: IF TRUE THIS WILL ERROR @TODO
            if (int(ORF_info['start']) <= sequence_length) or (int(ORF_info['end']) >= sequence_length * 2 + 1):
                add_to_score = False

        if add_to_score:
            total_burden += rbs_score * orf_length
            summary_file.write(f"  {str(ORF_info['start'])} Burden: {str(orf_length * rbs_score)} [ RBS Score: " +
                               f"{str(rbs_score)} ORF Length: {str(orf_length)} ]\n")

    summary_file.write("Total Burden: " + str(total_burden))

    # Filter predictions that we show
    expressed_ORFs = list(filter(lambda x:float(x["score"]) > options.rbs_score_cutoff, expressed_ORFs))

    #reload the gene annotation files
    gene_annotations = []
    if os.path.isfile(gene_annotation_file_name):
        gene_reader = csv.DictReader(open(gene_annotation_file_name))
        for row in gene_reader:
            row['position'] = row['end'] if row['strand'] == '+' else row['start']
            row['array'] = 0 if row['strand'] == '+' else int(row['end']) - int(row['start'])
            row['array_minus'] = int(row['end']) - int(row['start']) if row['strand'] == '+' else 0
            gene_annotations.append(row)

    # ---- CIRCULAR SEQUENCE ----

    if options.circular:

        # Remove ones with out of bounds starts/ends/positions (whichever is graphed)
        orf_predictions = list(filter(lambda x: (int(x['start']) >= sequence_length + 1) and
                                                (int(x['start']) <= sequence_length * 2), orf_predictions))
        rbs_predictions = list(filter(lambda x: (int(x['position']) >= sequence_length + 1) and
                                                (int(x['position']) <= sequence_length * 2), rbs_predictions))
        tss_predictions = list(filter(lambda x: (int(x['TSSpos']) >= sequence_length + 1) and
                                                (int(x['TSSpos']) <= sequence_length * 2), tss_predictions))
        rit_predictions = list(filter(lambda x: (int(x['end']) >= sequence_length + 1) and
                                                (int(x['end']) <= sequence_length * 2), rit_predictions))

        # Correct coordinates of remaining

        for d in rbs_predictions:
            d['position'] = str(int(d['position']) - sequence_length)

        for d in tss_predictions:
            d['TSSpos'] = str(int(d['TSSpos']) - sequence_length)

        for d in rit_predictions:
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
    with open(options.o + ".final.ORF.csv",'w') as final_orf_predictions_file:
        writer = csv.DictWriter(
            final_orf_predictions_file,
            fieldnames=["start", "end", "strand", "start_codon", 'length', 'rbs_score',
                        'array', 'array_minus', 'start_codon_position', 'second_half', 'translation']
            )
        writer.writeheader()
        writer.writerows(orf_predictions)
        final_orf_predictions_file.close()

    with open(options.o + ".final.RBS.csv", 'w') as final_rbs_predictions_file:
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
    gene_color = '#000000'
    gene_text_color = '#FFFFFF'
    marker_size = 20
    tss_score_max = 25

    tss_series = go.Scatter(
      x = [d['TSSpos'] for d in tss_predictions],
      y = [d['score'] for d in tss_predictions],
      mode = 'markers',
      name = 'TSS',
      text = [ ('-35: ' + d['box35seq'] + "<br>" + '-10: ' + d['box10seq'])  for d in tss_predictions],

      marker = dict(
        symbol=[ ('triangle-right'  if d['strand'] == '+' else 'triangle-left')  for d in tss_predictions],
        size=marker_size,
        color=TSS_color,
      ),

      error_y=dict(
                type='data',
                array=[0],
                arrayminus=[d['score'] for d in tss_predictions],
                width=0,
                color=TSS_color,
            )
    )

    rit_series = go.Scatter(
      x = [e['end'] for e in rit_predictions],
      y = [e['conf'] for e in rit_predictions],
      mode = 'markers',
      name = 'RIT',
      text = [ ('hairpin score:' + d['hairpin_score'] + '; tail score:' + d['tail_score'] + '<br>' +
                d['seq_hairpin_open'] + '-' + d['seq_hairpin_loop'] + '-' + d['seq_hairpin_close'] +
                '-' + d['seq_tail'][0:7] )  for d in rit_predictions],
      yaxis='y2',
      marker = dict(
        symbol=[ ('triangle-right'  if d['strand'] == '+' else 'triangle-left')  for d in rit_predictions],
        size=marker_size,
        color=RIT_color,
      ),
      error_y=dict(
              type='data',
              array=[0],
              arrayminus=[e['conf'] for e in rit_predictions],
              width=0,
              color=RIT_color,
          )
    )

    negatives_exist = True if "-" in [f['strand'] for f in expressed_ORFs] else False
    print(f"Negatives exist: {negatives_exist}")
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

    gene_series  = go.Scatter(
      x=[f['position'] for f in gene_annotations],
      y=[0.93]*len(gene_annotations),
      mode='markers+text',
      text=[ d['gene'] for d in gene_annotations ],
      name = 'Gene',
      yaxis='y4',
      marker = dict(
      symbol=[ ('triangle-right' if d['strand'] == '+' else 'triangle-left') for d in gene_annotations],
        size=marker_size,
        color=gene_color,
        ),
      error_x=dict(
              type='data',
              array=[f['array'] for f in gene_annotations],
              arrayminus=[f['array_minus'] for f in gene_annotations],
              thickness=10,
              width=0,
              color=gene_color,
      ),
    )


    ## Use display name at command line if it was given, otherwise file name
    if options.name:
      display_name = options.name
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
          range = [0, 110],
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
      ),
    )

    data = [tss_series, rit_series, rbs_series, nucleotide_sequence, gene_series]

    fig = dict( data=data, layout=layout )

    #if (options.pdf):
    #  fig.write_image(options.o+'.pdf')
    if (options.web):
      dev_content = plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')
      f = open(options.o+'.plot.div', 'w')
      f.write(dev_content)
      f.close()
    else:
      plotly.offline.plot(fig, filename=options.o+'.plot.html')


if __name__ == "__main__":
    '''Called if function is run as executable. Parses arguments and then launches main script'''
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
    main(options)
