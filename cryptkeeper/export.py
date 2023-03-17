import dataclasses
from collections import namedtuple
from typing import NamedTuple
from typing import List
import plotly
import plotly.graph_objs as go
import csv
import json

CDS_COLOR = '#4b61d1'
RBS_COLOR = "#d14b4b"
PROMOTER_COLOR = '#2ab717' 
TERMINATOR_COLOR = '#ff0000' 

@dataclasses.dataclass
class CryptResults:
    """Cryptkeeper results class"""
    name: str
    sequence: str
    translation_sites: List[NamedTuple]
    rho_dep_terminators: List[NamedTuple]
    rho_ind_terminators: List[NamedTuple]
    promoters: List[NamedTuple]
    annotations: List[NamedTuple]
    burden: float

    def plot(self, output_path: str):
        plot(self, output_path)

    def to_csv(self, output_path: str):
        '''Export results to csv files'''
        to_csv(self, output_path)

    def to_summary(self, output_path: str):
        '''Export results to summary file'''
        to_summary(self, output_path)

    def to_json(self, output_path: str):
        '''Export results to json file'''
        to_json(self, output_path)


def plot(results: CryptResults, output_path: str):
    """Plot results"""
    
    if not isinstance(results, CryptResults):
        raise TypeError("Results must be of type CryptResults")

    # Unpack results from CryptResults
    name = results.name
    promoter_calc_predictions = results.promoters
    transterm_predictions = results.rho_ind_terminators
    rhotermpredict_predictions = results.rho_dep_terminators
    expressed_ORFs = results.translation_sites
    features_list = results.annotations
    total_burden = results.burden
    sequence = results.sequence

    # Display predictions

    #settings
    TSS_color = PROMOTER_COLOR
    RIT_color = TERMINATOR_COLOR
    RBS_color = CDS_COLOR
    marker_size = 20
    tss_score_max = 8000


    tss_series = go.Scatter(
      x = [d.TSSpos for d in promoter_calc_predictions],
      y = [d.score for d in promoter_calc_predictions],
      mode = 'markers',
      name = 'TSS',
      text = [ (f'Score: {d.score}<br>-35: {d.box35seq}<br>-10: {d.box10seq}')  for d in promoter_calc_predictions],  # TODO Change this up for Rho dependant?
      marker = dict(
        symbol=[ ('triangle-right'  if d.strand == '+' else 'triangle-left')  for d in promoter_calc_predictions],
        size=marker_size,
        color=TSS_color,
      ),

      error_y=dict(
                type='data',
                array=[0],
                arrayminus=[d.score for d in promoter_calc_predictions],
                width=0,
                color=TSS_color,
            )
    )

    rit_transterm_series = go.Scatter(
      x = [int(e.end) for e in transterm_predictions],
      y = [int(e.conf) for e in transterm_predictions],
      mode = 'markers',
      name = 'RIT',
      yaxis='y2',
      text = [ ('hairpin score:' + d.hairpin_score + '; tail score:' + d.tail_score + '<br>' +
                d.seq_hairpin_open + '-' + d.seq_hairpin_loop + '-' + d.seq_hairpin_close +
                '-' + d.seq_tail[0:7] )  for d in transterm_predictions],
      marker = dict(
        symbol=[ ('triangle-right'  if d.strand == '+' else 'triangle-left')  for d in transterm_predictions],
        size=marker_size,
        color=RIT_color,
      ),
      error_y=dict(
              type='data',
              array=[0],
              arrayminus=[e.conf for e in transterm_predictions],
              width=0,
              color=RIT_color,
          )
    )

    rdt_rhotermpredict_series = go.Scatter(
      x = [int(e.end_rut) for e in rhotermpredict_predictions],
      y = [int(sum(e.scores)) for e in rhotermpredict_predictions],
      mode = 'markers',
      name = 'RDT',
      yaxis='y6',
      text = [ ('placeholder')  for d in rhotermpredict_predictions],
      marker = dict(
        symbol=[ ('triangle-right'  if d.strand == '+' else 'triangle-left')  for d in rhotermpredict_predictions],
        size=marker_size,
        color=RIT_color,
      ),
      error_y=dict(
              type='data',
              array=[0],
              arrayminus=[sum(e.scores) for e in rhotermpredict_predictions],
              width=0,
              color=RIT_color,
          )
    )

    negatives_exist = True if "-" in [f.strand for f in expressed_ORFs] else False
    rbs_series = go.Scatter(
      x = [f.start for f in expressed_ORFs],
      y = [f.expression for f in expressed_ORFs],
      mode = 'markers',
      text = [ ('start codon:' + d.start_codon + "<br>Burden ×10<sup>6</sup>: " +
                "{0:.6f}".format(d.burden/1000000) + " (" + "{0:.2f}".format(100*d.burden/total_burden) + "%}")
               for d in expressed_ORFs],
      name = 'RBS',
      yaxis='y3',
      marker = dict(
        symbol=[ ('triangle-right'  if d.strand == '+' else 'triangle-left')  for d in expressed_ORFs],
        size=marker_size,
        color=RBS_color,
      ),
      error_y=dict(
              type='data',
              array=[0],
              arrayminus=[f.expression for f in expressed_ORFs],
              width=0,
              color=RBS_color,
      ),
      error_x=dict(
              type='data',
              array=[f.array if f.strand == '+' else 0 for f in expressed_ORFs],
              arrayminus=[f.array if f.strand == '-' else 0 for f in expressed_ORFs],
              width=3,
              color=RBS_color,
      ),
    )

    nucleotide_sequence  = go.Scatter(
      x=list(range(0,len(sequence)-1)),
      y=[0.97]*len(sequence),
      mode='text',
      name = 'Sequence',
      text=list(sequence),
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

    for feature in features_list:  # Draw features on Cryptkeeper Results
        start_pos_x = feature.start-1
        end_pos_x = feature.end-1
        y_pos = base_y_level-feature.nest_level*0.025

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
                                              text=feature.name,
                                              name="Annotations",
                                              mode='none',
                                              marker=dict(
                                                  line=dict(width=0.5)),
                                              hoveron='fills',
                                              legendgroup='1',
                                              showlegend=False,
                                              fillcolor=feature.color
                                              ))

    ## Use display name at command line if it was given, otherwise file name
    if name:
        display_name = name
    else:
        display_name = output_path

    layout = go.Layout(
      title='CryptKeeper Results: ' + display_name + "  Burden ×10<sup>6</sup>: " + "{0:.2f}".format(total_burden/1000000),
      hovermode = 'closest',
      autosize=True,

      xaxis=dict(
            domain=[0.20, 1.0]
      ),
      yaxis=dict(
        title='TSS score',
        range = [1000, tss_score_max+0.2],
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

    plotly.offline.plot(fig, filename=output_path)

def to_csv(results: CryptResults, output_path: str):
    """Write results to csv files"""
    def _csv_writer():
        if not dataset:
            return
        with open(writepath, 'w', encoding='utf8') as file:
            header = dataset[0]._fields
            writer = csv.writer(file)
            writer.writerow(header)
            for data in dataset:
                writer.writerow(data._asdict().values())

    # Write RBS predictions to CSV
    # Headers: Position, Strand, Start_codon, score1, score2
    # @TODO: Make scores descriptive
    dataset = results.translation_sites
    writepath = output_path + '_rbs.csv'
    _csv_writer()

    # Rho Independent Termination to CSV
    # Headers: Start, end, strand, coef, hairpin score, tail score
    # old system also has hairpin upstream, open, loop, close, and tail
    dataset = results.rho_ind_terminators
    writepath = output_path + '_rit.csv'
    _csv_writer()

    # Promoter Sites:
    # Headers: Score, strand, TTSpos, box35pos, box35seq, box10pos, box10seq
    dataset = results.promoters
    writepath = output_path + '_tss.csv'
    _csv_writer()

    # Also need to output Rho Dependent Termination to CSV, which is new
    dataset = results.rho_dep_terminators
    writepath = output_path + '_rdt.csv'
    _csv_writer()

def to_summary(results: CryptResults, output_path: str):
    """Write summary of results to file"""
    with open(output_path, 'w') as file:
        file.write("CryptKeeper Results Summary\n")
        file.write("Total Burden: " + str(results.burden) + "\n")
        file.write("Predicted Promoters: " + str(len(results.promoters)) + "\n")
        file.write("Predicted Translation: " + str(len(results.translation_sites)) + "\n")
        file.write("Predicted Rho Independent Terminators: " + str(len(results.rho_ind_terminators)) + "\n")
        file.write("Total Rho Dependent Terminators: " + str(len(results.rho_dep_terminators)) + "\n")

def to_json(results: CryptResults, output_path: str):
    """Write results to json file"""
    with open(output_path, 'w') as file:
        json.dump(dataclasses.asdict(results), file, indent=4)

def from_json(input_path: str) -> CryptResults:
    """Read results from json file and load them into a CryptResults object"""
    with open(input_path, 'r') as file:
        # @TODO: Add sanity checking to json file
        data = json.load(file)
        results = CryptResults(**data)

    # Correct the data type from the json
    expressed_orf = namedtuple("orf_result", "start, end, expression, burden, dG, array, start_codon, strand")
    promocalc_hit = namedtuple('promoter_calculator_result', 'seq score strand TSSpos box35pos box35seq box10pos box10seq')
    transterm_hit = namedtuple("transterm_result","start end strand conf hairpin_score tail_score seq_upstream seq_hairpin_open seq_hairpin_loop seq_hairpin_close seq_tail")
    rtp_hit = namedtuple('rdtresult', """strand,
                                        c_over_g,
                                        start_rut,
                                        end_rut,
                                        term_seq,
                                        downstream_seq,
                                        palindromes,
                                        pause_concensus,
                                        scores""")
    annotations = namedtuple('feature', ['name', 'strand', 'start', 'end', 'color', 'nest_level'])

    for i, hit in enumerate(results.translation_sites):
        results.translation_sites[i] = expressed_orf(*hit)
    for i, hit in enumerate(results.promoters):
        results.promoters[i] = promocalc_hit(*hit)
    for i, hit in enumerate(results.rho_ind_terminators):
        results.rho_ind_terminators[i] = transterm_hit(*hit)
    for i, hit in enumerate(results.rho_dep_terminators):
        results.rho_dep_terminators[i] = rtp_hit(*hit)
    for i, hit in enumerate(results.annotations):
        results.annotations[i] = annotations(*hit)

    return results


# Blank EOF line
